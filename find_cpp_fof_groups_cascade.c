#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr const char *DATA_DIR = "data";
constexpr const char *OUTPUT_DIR = "outdir";
constexpr const char *TIMINGS_PATH = "outdir/cpp_fof_timings.txt";
constexpr const char *GROUPS_PATH = "outdir/cpp_groups.txt";
constexpr const char *MEMBERS_PATH = "outdir/cpp_members.txt";

constexpr double BOX_SIZE = 1.0;
constexpr double B = 0.2;
constexpr int MIN_GROUP_SIZE = 10;
constexpr int FIRST_POWER = 12;
constexpr int LAST_POWER = 20;
constexpr int N_TRIALS = 3;

class Particle {
public:
    Particle() = default;
    Particle(double x, double y, double z) : x_{x, y, z} {}

    double operator[](std::size_t axis) const { return x_[axis]; }

private:
    double x_[3]{};
};

class ParticleSet {
public:
    explicit ParticleSet(const std::string &path) { load(path); }

    std::size_t size() const { return particles_.size(); }
    const Particle &operator[](std::size_t index) const { return particles_[index]; }

private:
    std::vector<Particle> particles_;

    void load(const std::string &path) {
        particles_.clear();
        std::ifstream input(path);
        if (!input) throw std::runtime_error("cannot open " + path);

        std::string line;
        int line_number = 0;
        while (std::getline(input, line)) {
            ++line_number;
            std::istringstream parser(line);
            double x = 0.0, y = 0.0, z = 0.0;
            std::string extra;
            if (!(parser >> x >> y >> z) || (parser >> extra)) {
                throw std::runtime_error(path + ":" + std::to_string(line_number) +
                                         ": expected 3 numeric columns");
            }
            particles_.emplace_back(x, y, z);
        }

        if (input.bad()) throw std::runtime_error("error while reading " + path);
    }
};

class UnionFind {
public:
    explicit UnionFind(std::size_t size) : parent_(size), rank_(size, 0) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    int find(int item) {
        int root = item;
        while (parent_[root] != root) root = parent_[root];

        while (parent_[item] != item) {
            const int next = parent_[item];
            parent_[item] = root;
            item = next;
        }
        return root;
    }

    void unite(int left, int right) {
        int left_root = find(left);
        int right_root = find(right);
        if (left_root == right_root) return;

        if (rank_[left_root] < rank_[right_root]) {
            parent_[left_root] = right_root;
        } else if (rank_[left_root] > rank_[right_root]) {
            parent_[right_root] = left_root;
        } else {
            parent_[right_root] = left_root;
            ++rank_[left_root];
        }
    }

private:
    std::vector<int> parent_;
    std::vector<unsigned char> rank_;
};

class Group {
public:
    explicit Group(int root = -1) : root_(root) {}

    void reserve(std::size_t size) { members_.reserve(size); }
    void add_member(int one_based_particle_id) { members_.push_back(one_based_particle_id); }

    int root() const { return root_; }
    std::size_t size() const { return members_.size(); }
    int first_member() const { return members_.front(); }
    const std::vector<int> &members() const { return members_; }

private:
    int root_;
    std::vector<int> members_;
};

struct FoFResult {
    std::vector<Group> kept_groups;
    int raw_group_count = 0;
};

class CellList {
public:
    CellList(const ParticleSet &particles, double linking_length)
        : particles_(particles) {
        const double target_cell_size = 4.0 * linking_length;
        ncell_ = static_cast<int>(std::floor(BOX_SIZE / target_cell_size));
        if (ncell_ < 1) ncell_ = 1;
        cell_size_ = BOX_SIZE / static_cast<double>(ncell_);
        ncell2_ = ncell_ * ncell_;
        const std::size_t ncells = static_cast<std::size_t>(ncell_) *
                                   static_cast<std::size_t>(ncell_) *
                                   static_cast<std::size_t>(ncell_);

        head_.assign(ncells, -1);
        next_.assign(particles_.size(), -1);

        for (std::size_t i = 0; i < particles_.size(); ++i) {
            const int cell = cell_index(particles_[i]);
            next_[i] = head_[static_cast<std::size_t>(cell)];
            head_[static_cast<std::size_t>(cell)] = static_cast<int>(i);
        }
    }

    template <typename Callback>
    void for_candidate_pairs(Callback &&callback) const {
        for (int cx = 0; cx < ncell_; ++cx) {
            for (int cy = 0; cy < ncell_; ++cy) {
                for (int cz = 0; cz < ncell_; ++cz) {
                    const int cell = flat_index(cx, cy, cz);
                    for (int i = head_[static_cast<std::size_t>(cell)]; i >= 0; i = next_[static_cast<std::size_t>(i)]) {
                        for (int dx = -1; dx <= 1; ++dx) {
                            const int nx = wrap_cell(cx + dx);
                            for (int dy = -1; dy <= 1; ++dy) {
                                const int ny = wrap_cell(cy + dy);
                                for (int dz = -1; dz <= 1; ++dz) {
                                    const int nz = wrap_cell(cz + dz);
                                    const int neighbour_cell = flat_index(nx, ny, nz);
                                    for (int j = head_[static_cast<std::size_t>(neighbour_cell)]; j >= 0; j = next_[static_cast<std::size_t>(j)]) {
                                        if (j > i) callback(i, j);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

private:
    const ParticleSet &particles_;
    int ncell_ = 1;
    int ncell2_ = 1;
    double cell_size_ = BOX_SIZE;
    std::vector<int> head_;
    std::vector<int> next_;

    int wrap_cell(int value) const {
        if (value < 0) return value + ncell_;
        if (value >= ncell_) return value - ncell_;
        return value;
    }

    int coordinate_to_cell(double x) const {
        int cell = static_cast<int>(std::floor(x / cell_size_));
        if (cell < 0) cell = 0;
        if (cell >= ncell_) cell = ncell_ - 1;
        return cell;
    }

    int flat_index(int cx, int cy, int cz) const {
        return (cx * ncell_ + cy) * ncell_ + cz;
    }

    int cell_index(const Particle &particle) const {
        return flat_index(
            coordinate_to_cell(particle[0]),
            coordinate_to_cell(particle[1]),
            coordinate_to_cell(particle[2])
        );
    }
};

class FoFFinder {
public:
    explicit FoFFinder(const ParticleSet &particles)
        : particles_(particles), union_find_(particles.size()) {}

    FoFResult find_groups() {
        link_particles();
        return build_groups();
    }

private:
    const ParticleSet &particles_;
    UnionFind union_find_;

    static double periodic_delta(double left, double right) {
        double delta = std::fabs(left - right);
        if (delta > 0.5 * BOX_SIZE) delta = BOX_SIZE - delta;
        return delta;
    }

    void link_particles() {
        const std::size_t n = particles_.size();
        if (n == 0) return;

        const double linking_length = B * std::pow((BOX_SIZE * BOX_SIZE * BOX_SIZE) / static_cast<double>(n), 1.0 / 3.0);
        const double linking_length2 = linking_length * linking_length;
        const CellList cells(particles_, linking_length);

        cells.for_candidate_pairs([&](int i, int j) {
            const Particle &pi = particles_[static_cast<std::size_t>(i)];
            const Particle &pj = particles_[static_cast<std::size_t>(j)];
            const double dx = periodic_delta(pi[0], pj[0]);
            if (dx * dx > linking_length2) return;
            const double dy = periodic_delta(pi[1], pj[1]);
            const double dxy2 = dx * dx + dy * dy;
            if (dxy2 > linking_length2) return;
            const double dz = periodic_delta(pi[2], pj[2]);
            if (dxy2 + dz * dz <= linking_length2) {
                union_find_.unite(i, j);
            }
        });
    }

    FoFResult build_groups() {
        const std::size_t n = particles_.size();
        std::vector<int> sizes(n, 0);
        std::vector<int> group_slot(n, -1);

        for (std::size_t i = 0; i < n; ++i) {
            ++sizes[static_cast<std::size_t>(union_find_.find(static_cast<int>(i)))];
        }

        FoFResult result;
        for (std::size_t root = 0; root < n; ++root) {
            if (sizes[root] > 0) ++result.raw_group_count;
            if (sizes[root] >= MIN_GROUP_SIZE) {
                group_slot[root] = static_cast<int>(result.kept_groups.size());
                result.kept_groups.emplace_back(static_cast<int>(root));
                result.kept_groups.back().reserve(static_cast<std::size_t>(sizes[root]));
            }
        }

        for (std::size_t i = 0; i < n; ++i) {
            const int root = union_find_.find(static_cast<int>(i));
            const int slot = group_slot[static_cast<std::size_t>(root)];
            if (slot >= 0) result.kept_groups[static_cast<std::size_t>(slot)].add_member(static_cast<int>(i) + 1);
        }

        std::sort(result.kept_groups.begin(), result.kept_groups.end(), [](const Group &a, const Group &b) {
            if (a.size() != b.size()) return a.size() > b.size();
            return a.first_member() < b.first_member();
        });

        return result;
    }
};

class OutputWriter {
public:
    OutputWriter(std::string groups_path, std::string members_path)
        : groups_path_(std::move(groups_path)), members_path_(std::move(members_path)) {}

    void write(const std::vector<Group> &groups, const ParticleSet &particles) const {
        std::filesystem::create_directories(OUTPUT_DIR);

        std::ofstream groups_file(groups_path_);
        if (!groups_file) throw std::runtime_error("cannot open " + groups_path_);
        std::ofstream members_file(members_path_);
        if (!members_file) throw std::runtime_error("cannot open " + members_path_);

        groups_file << "# group_id size com_x com_y com_z\n";
        members_file << "# group_id size member_ids...\n";

        groups_file << std::scientific << std::setprecision(16);
        for (std::size_t i = 0; i < groups.size(); ++i) {
            const int group_id = static_cast<int>(i) + 1;
            const Group &group = groups[i];
            const auto center = periodic_center(group, particles);

            groups_file << group_id << ' ' << group.size() << ' '
                        << center[0] << ' ' << center[1] << ' ' << center[2] << '\n';

            members_file << group_id << ' ' << group.size();
            for (const int member : group.members()) members_file << ' ' << member;
            members_file << '\n';
        }
    }

private:
    std::string groups_path_;
    std::string members_path_;

    static double mod_box(double x) {
        double y = std::fmod(x, BOX_SIZE);
        if (y < 0.0) y += BOX_SIZE;
        return y;
    }

    static std::array<double, 3> periodic_center(const Group &group, const ParticleSet &particles) {
        std::array<double, 3> center{0.0, 0.0, 0.0};
        const Particle &reference = particles[static_cast<std::size_t>(group.first_member() - 1)];

        for (std::size_t axis = 0; axis < 3; ++axis) {
            double sum = 0.0;
            for (const int member_id : group.members()) {
                const Particle &particle = particles[static_cast<std::size_t>(member_id - 1)];
                double delta = std::fmod(particle[axis] - reference[axis] + 0.5 * BOX_SIZE, BOX_SIZE);
                if (delta < 0.0) delta += BOX_SIZE;
                delta -= 0.5 * BOX_SIZE;
                sum += reference[axis] + delta;
            }
            center[axis] = mod_box(sum / static_cast<double>(group.size()));
        }
        return center;
    }
};

std::string cascade_path(int power) {
    const int n = 1 << power;
    std::ostringstream path;
    path << DATA_DIR << "/ics_" << std::setw(7) << std::setfill('0') << n << ".txt";
    return path.str();
}

std::vector<std::string> find_cascade_inputs() {
    std::vector<std::string> paths;
    for (int power = FIRST_POWER; power <= LAST_POWER; ++power) {
        std::string path = cascade_path(power);
        if (std::filesystem::exists(path)) paths.push_back(path);
    }
    return paths;
}

}  // namespace

int main() {
    try {
        std::filesystem::create_directories(OUTPUT_DIR);
        const std::vector<std::string> input_paths = find_cascade_inputs();
        if (input_paths.empty()) {
            throw std::runtime_error("no cascade IC files found: expected data/ics_0004096.txt ... data/ics_1048576.txt");
        }

        std::ofstream timings(TIMINGS_PATH);
        if (!timings) throw std::runtime_error(std::string("cannot open ") + TIMINGS_PATH);
        timings << "# input_file n_particles trial analysis_seconds n_raw_groups n_kept_groups\n";
        timings << std::scientific << std::setprecision(8);

        FoFResult last_result;
        ParticleSet last_particles(input_paths.front());

        for (const std::string &input_path : input_paths) {
            const ParticleSet particles(input_path);
            std::cout << "Timing " << input_path << " (" << particles.size() << " particles)\n";

            for (int trial = 1; trial <= N_TRIALS; ++trial) {
                const auto start = std::chrono::steady_clock::now();
                FoFFinder finder(particles);
                FoFResult result = finder.find_groups();
                const auto stop = std::chrono::steady_clock::now();

                const double seconds = std::chrono::duration<double>(stop - start).count();
                timings << input_path << ' '
                        << particles.size() << ' '
                        << trial << ' '
                        << seconds << ' '
                        << result.raw_group_count << ' '
                        << result.kept_groups.size() << '\n';

                if (&input_path == &input_paths.back() && trial == N_TRIALS) {
                    last_result = std::move(result);
                    last_particles = particles;
                }
            }
        }

        const OutputWriter writer(GROUPS_PATH, MEMBERS_PATH);
        writer.write(last_result.kept_groups, last_particles);

        std::cout << "Wrote timings to " << TIMINGS_PATH << '\n';
        std::cout << "Wrote final cascade groups to " << GROUPS_PATH << '\n';
        std::cout << "Wrote final cascade memberships to " << MEMBERS_PATH << '\n';
    } catch (const std::exception &error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
