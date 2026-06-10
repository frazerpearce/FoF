#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#define OUTPUT_DIR "outdir"
#define TIMINGS_PATH "outdir/c_fof_timings.txt"

#define BOX_SIZE 1.0
#define B 0.2
#define MIN_GROUP_SIZE 10
#define N_TRIALS 3
#define FIRST_POWER 12
#define LAST_POWER 20

typedef struct {
    double x[3];
} Particle;

typedef struct {
    int parent;
    unsigned char rank;
} UFNode;

typedef struct {
    int root;
    int size;
    int *members;
} Group;

typedef struct {
    int64_t *keys;
    int *heads;
    size_t capacity;
} CellHash;

static void die_errno(const char *message)
{
    fprintf(stderr, "Error: %s: %s\n", message, strerror(errno));
    exit(EXIT_FAILURE);
}

static void die_message(const char *message)
{
    fprintf(stderr, "Error: %s\n", message);
    exit(EXIT_FAILURE);
}

static void *xmalloc(size_t nbytes)
{
    void *ptr = malloc(nbytes);
    if (ptr == NULL && nbytes != 0) die_errno("malloc failed");
    return ptr;
}

static void *xrealloc(void *old, size_t nbytes)
{
    void *ptr = realloc(old, nbytes);
    if (ptr == NULL && nbytes != 0) die_errno("realloc failed");
    return ptr;
}

static size_t next_power_of_two(size_t value)
{
    size_t capacity = 1;
    while (capacity < value) capacity <<= 1;
    return capacity;
}

static uint64_t hash_u64(uint64_t x)
{
    x ^= x >> 33;
    x *= UINT64_C(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x *= UINT64_C(0xc4ceb9fe1a85ec53);
    x ^= x >> 33;
    return x;
}

static size_t hash_slot(const CellHash *hash, int64_t key)
{
    size_t mask = hash->capacity - 1;
    size_t slot = (size_t)hash_u64((uint64_t)key) & mask;
    while (hash->keys[slot] != -1 && hash->keys[slot] != key) {
        slot = (slot + 1) & mask;
    }
    return slot;
}

static CellHash cell_hash_create(int n_particles)
{
    CellHash hash;
    hash.capacity = next_power_of_two((size_t)n_particles * 4u + 16u);
    hash.keys = xmalloc(hash.capacity * sizeof(*hash.keys));
    hash.heads = xmalloc(hash.capacity * sizeof(*hash.heads));
    for (size_t i = 0; i < hash.capacity; i++) {
        hash.keys[i] = -1;
        hash.heads[i] = -1;
    }
    return hash;
}

static void cell_hash_free(CellHash *hash)
{
    free(hash->keys);
    free(hash->heads);
    hash->keys = NULL;
    hash->heads = NULL;
    hash->capacity = 0;
}

static int *cell_hash_head_ptr(CellHash *hash, int64_t key)
{
    size_t slot = hash_slot(hash, key);
    if (hash->keys[slot] == -1) {
        hash->keys[slot] = key;
        hash->heads[slot] = -1;
    }
    return &hash->heads[slot];
}

static int cell_hash_get_head(const CellHash *hash, int64_t key)
{
    size_t mask = hash->capacity - 1;
    size_t slot = (size_t)hash_u64((uint64_t)key) & mask;
    while (hash->keys[slot] != -1) {
        if (hash->keys[slot] == key) return hash->heads[slot];
        slot = (slot + 1) & mask;
    }
    return -1;
}

static int uf_find(UFNode *nodes, int item)
{
    int root = item;
    while (nodes[root].parent != root) root = nodes[root].parent;

    while (nodes[item].parent != item) {
        int next = nodes[item].parent;
        nodes[item].parent = root;
        item = next;
    }
    return root;
}

static void uf_union(UFNode *nodes, int left, int right)
{
    int left_root = uf_find(nodes, left);
    int right_root = uf_find(nodes, right);

    if (left_root == right_root) return;

    if (nodes[left_root].rank < nodes[right_root].rank) {
        nodes[left_root].parent = right_root;
    } else if (nodes[left_root].rank > nodes[right_root].rank) {
        nodes[right_root].parent = left_root;
    } else {
        nodes[right_root].parent = left_root;
        nodes[left_root].rank++;
    }
}

static Particle *load_positions(const char *path, int *n_particles)
{
    FILE *file = fopen(path, "r");
    if (file == NULL) die_errno(path);

    int capacity = 1024;
    int count = 0;
    Particle *particles = xmalloc((size_t)capacity * sizeof(*particles));

    char line[4096];
    int line_number = 0;
    while (fgets(line, sizeof(line), file) != NULL) {
        double x, y, z;
        char extra;
        line_number++;
        if (sscanf(line, " %lf %lf %lf %c", &x, &y, &z, &extra) != 3) {
            fprintf(stderr, "Error: %s:%d: expected 3 numeric columns\n", path, line_number);
            exit(EXIT_FAILURE);
        }
        if (count == capacity) {
            capacity *= 2;
            particles = xrealloc(particles, (size_t)capacity * sizeof(*particles));
        }
        particles[count].x[0] = x;
        particles[count].x[1] = y;
        particles[count].x[2] = z;
        count++;
    }

    if (ferror(file)) die_errno(path);
    fclose(file);

    *n_particles = count;
    return particles;
}

static double periodic_delta(double left, double right)
{
    double delta = fabs(left - right);
    if (delta > 0.5 * BOX_SIZE) delta = BOX_SIZE - delta;
    return delta;
}

static int periodic_cell_index(double coordinate, int n_cells)
{
    int index = (int)floor(coordinate * (double)n_cells);
    if (index < 0) index = 0;
    if (index >= n_cells) index = n_cells - 1;
    return index;
}

static int wrap_cell_index(int index, int n_cells)
{
    if (index < 0) return index + n_cells;
    if (index >= n_cells) return index - n_cells;
    return index;
}

static int64_t cell_key(int ix, int iy, int iz, int n_cells)
{
    return (int64_t)ix + (int64_t)n_cells * ((int64_t)iy + (int64_t)n_cells * (int64_t)iz);
}

static void find_groups_union(const Particle *particles, int n_particles, UFNode *nodes)
{
    for (int i = 0; i < n_particles; i++) {
        nodes[i].parent = i;
        nodes[i].rank = 0;
    }

    if (n_particles == 0) return;

    const double linking_length = B * pow(BOX_SIZE * BOX_SIZE * BOX_SIZE / (double)n_particles, 1.0 / 3.0);
    const double linking_length2 = linking_length * linking_length;

    int n_cells = (int)floor(BOX_SIZE / linking_length);
    if (n_cells < 1) n_cells = 1;

    int *cell_x = xmalloc((size_t)n_particles * sizeof(*cell_x));
    int *cell_y = xmalloc((size_t)n_particles * sizeof(*cell_y));
    int *cell_z = xmalloc((size_t)n_particles * sizeof(*cell_z));
    int *next = xmalloc((size_t)n_particles * sizeof(*next));
    CellHash hash = cell_hash_create(n_particles);

    for (int i = 0; i < n_particles; i++) {
        cell_x[i] = periodic_cell_index(particles[i].x[0], n_cells);
        cell_y[i] = periodic_cell_index(particles[i].x[1], n_cells);
        cell_z[i] = periodic_cell_index(particles[i].x[2], n_cells);

        const int64_t key = cell_key(cell_x[i], cell_y[i], cell_z[i], n_cells);
        int *head = cell_hash_head_ptr(&hash, key);
        next[i] = *head;
        *head = i;
    }

    for (int i = 0; i < n_particles; i++) {
        const Particle *pi = particles + i;
        const int ix0 = cell_x[i];
        const int iy0 = cell_y[i];
        const int iz0 = cell_z[i];

        for (int dz_cell = -1; dz_cell <= 1; dz_cell++) {
            const int iz = wrap_cell_index(iz0 + dz_cell, n_cells);
            for (int dy_cell = -1; dy_cell <= 1; dy_cell++) {
                const int iy = wrap_cell_index(iy0 + dy_cell, n_cells);
                for (int dx_cell = -1; dx_cell <= 1; dx_cell++) {
                    const int ix = wrap_cell_index(ix0 + dx_cell, n_cells);
                    const int64_t key = cell_key(ix, iy, iz, n_cells);

                    for (int j = cell_hash_get_head(&hash, key); j >= 0; j = next[j]) {
                        if (j <= i) continue;

                        const Particle *pj = particles + j;
                        const double dx = periodic_delta(pi->x[0], pj->x[0]);
                        if (dx * dx > linking_length2) continue;
                        const double dy = periodic_delta(pi->x[1], pj->x[1]);
                        if (dx * dx + dy * dy > linking_length2) continue;
                        const double dz = periodic_delta(pi->x[2], pj->x[2]);
                        if (dx * dx + dy * dy + dz * dz <= linking_length2) {
                            uf_union(nodes, i, j);
                        }
                    }
                }
            }
        }
    }

    cell_hash_free(&hash);
    free(next);
    free(cell_z);
    free(cell_y);
    free(cell_x);
}

static int compare_groups(const void *left, const void *right)
{
    const Group *a = (const Group *)left;
    const Group *b = (const Group *)right;

    if (a->size != b->size) return b->size - a->size;
    return a->members[0] - b->members[0];
}

static Group *build_kept_groups(UFNode *nodes, int n_particles, int *n_groups)
{
    int *sizes = xmalloc((size_t)n_particles * sizeof(*sizes));
    int *offsets = xmalloc((size_t)n_particles * sizeof(*offsets));
    int *used = xmalloc((size_t)n_particles * sizeof(*used));

    for (int i = 0; i < n_particles; i++) {
        sizes[i] = 0;
        offsets[i] = 0;
        used[i] = 0;
    }

    for (int i = 0; i < n_particles; i++) {
        int root = uf_find(nodes, i);
        sizes[root]++;
    }

    int kept = 0;
    for (int i = 0; i < n_particles; i++) {
        if (sizes[i] >= MIN_GROUP_SIZE) kept++;
    }

    Group *groups = xmalloc((size_t)kept * sizeof(*groups));
    int index = 0;
    for (int i = 0; i < n_particles; i++) {
        if (sizes[i] >= MIN_GROUP_SIZE) {
            groups[index].root = i;
            groups[index].size = sizes[i];
            groups[index].members = xmalloc((size_t)sizes[i] * sizeof(*groups[index].members));
            used[i] = index + 1;
            index++;
        }
    }

    for (int i = 0; i < n_particles; i++) {
        int root = uf_find(nodes, i);
        int group_slot = used[root];
        if (group_slot != 0) {
            Group *group = groups + (group_slot - 1);
            group->members[offsets[root]++] = i + 1;
        }
    }

    qsort(groups, (size_t)kept, sizeof(*groups), compare_groups);

    free(sizes);
    free(offsets);
    free(used);

    *n_groups = kept;
    return groups;
}

static void make_output_dir(void)
{
    if (mkdir(OUTPUT_DIR, 0777) != 0 && errno != EEXIST) die_errno(OUTPUT_DIR);
}

static void free_groups(Group *groups, int n_groups)
{
    if (groups == NULL) return;
    for (int i = 0; i < n_groups; i++) free(groups[i].members);
    free(groups);
}

static int count_raw_groups(UFNode *nodes, int n_particles)
{
    int count = 0;
    for (int i = 0; i < n_particles; i++) {
        if (uf_find(nodes, i) == i) count++;
    }
    return count;
}

static double seconds_now(void)
{
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) die_errno("clock_gettime failed");
    return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

static int file_exists(const char *path)
{
    FILE *file = fopen(path, "r");
    if (file == NULL) return 0;
    fclose(file);
    return 1;
}

static void benchmark_file(FILE *timings_file, const char *input_path)
{
    int n_particles = 0;
    Particle *particles = load_positions(input_path, &n_particles);

    printf("Timing %s (%d particles)\n", input_path, n_particles);

    for (int trial = 1; trial <= N_TRIALS; trial++) {
        UFNode *nodes = xmalloc((size_t)n_particles * sizeof(*nodes));

        const double start = seconds_now();
        find_groups_union(particles, n_particles, nodes);
        const int n_raw_groups = count_raw_groups(nodes, n_particles);
        int n_kept_groups = 0;
        Group *groups = build_kept_groups(nodes, n_particles, &n_kept_groups);
        const double end = seconds_now();

        fprintf(timings_file, "%s %d %d %.9e %d %d\n",
                input_path, n_particles, trial, end - start, n_raw_groups, n_kept_groups);
        fflush(timings_file);

        printf("  trial %d: %.6f s, raw=%d, kept=%d\n",
               trial, end - start, n_raw_groups, n_kept_groups);

        free_groups(groups, n_kept_groups);
        free(nodes);
    }

    free(particles);
}

int main(void)
{
    make_output_dir();

    FILE *timings_file = fopen(TIMINGS_PATH, "w");
    if (timings_file == NULL) die_errno(TIMINGS_PATH);

    fprintf(timings_file, "# input_file n_particles trial analysis_seconds n_raw_groups n_kept_groups\n");

    int n_files = 0;
    for (int power = FIRST_POWER; power <= LAST_POWER; power++) {
        const int n_particles = 1 << power;
        char input_path[256];
        snprintf(input_path, sizeof(input_path), "data/ics_%07d.txt", n_particles);

        if (!file_exists(input_path)) {
            printf("Skipping missing %s\n", input_path);
            continue;
        }

        benchmark_file(timings_file, input_path);
        n_files++;
    }

    fclose(timings_file);

    if (n_files == 0) die_message("no cascade IC files found: expected data/ics_0004096.txt ... data/ics_1048576.txt");

    printf("Wrote timings for %d IC files to %s\n", n_files, TIMINGS_PATH);

    return EXIT_SUCCESS;
}
