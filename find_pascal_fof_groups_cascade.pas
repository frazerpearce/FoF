program FoFCascadeTimer;

{$mode objfpc}{$H+}

uses
  SysUtils, Math;

const
  MinKeep = 10;
  BoxSize = 1.0;
  B = 0.2;
  FirstPower = 12;
  LastPower = 20;
  NTrials = 3;

type
  TRealArray = array of Double;
  TIntArray = array of LongInt;
  TInt64Array = array of Int64;

var
  X, Y, Z: TRealArray;
  Parent, Rank, Root, GroupSize: TIntArray;
  GroupFirst, GroupId: TIntArray;
  GroupHead, GroupTail, GroupNext: TIntArray;
  Order, Work, Members: TIntArray;
  NextParticle, CellX, CellY, CellZ: TIntArray;
  UsedSlot: TIntArray;
  HashKey: TInt64Array;
  HashHead: TIntArray;
  TimingFile, InputFile, GroupFile, MemberFile: Text;
  N, NCell, NUsed, NGroups, NKeep: LongInt;
  HashSize: LongInt;
  LinkLength, LinkLength2: Double;
  CurrentInputPath: String;

function NextPowerOfTwo(Value: LongInt): LongInt;
var
  ResultValue: LongInt;
begin
  ResultValue := 1;
  while ResultValue < Value do
    ResultValue := ResultValue * 2;
  NextPowerOfTwo := ResultValue
end;

function PositiveMod64(Value, ModulusValue: Int64): Int64;
begin
  PositiveMod64 := Value mod ModulusValue;
  if PositiveMod64 < 0 then
    PositiveMod64 := PositiveMod64 + ModulusValue
end;

function WrapUnit(Value: Double): Double;
begin
  WrapUnit := Value - Floor(Value / BoxSize) * BoxSize;
  if WrapUnit < 0.0 then
    WrapUnit := WrapUnit + BoxSize
end;

function FindRoot(Item: LongInt): LongInt;
var
  Current, Following: LongInt;
begin
  Current := Item;
  while Parent[Current] <> Current do
  begin
    Following := Parent[Current];
    Parent[Current] := Parent[Following];
    Current := Parent[Current]
  end;
  FindRoot := Current
end;

procedure Join(LeftItem, RightItem: LongInt);
var
  LeftRoot, RightRoot: LongInt;
begin
  LeftRoot := FindRoot(LeftItem);
  RightRoot := FindRoot(RightItem);
  if LeftRoot <> RightRoot then
    if Rank[LeftRoot] < Rank[RightRoot] then
      Parent[LeftRoot] := RightRoot
    else if Rank[LeftRoot] > Rank[RightRoot] then
      Parent[RightRoot] := LeftRoot
    else
    begin
      Parent[RightRoot] := LeftRoot;
      Rank[LeftRoot] := Rank[LeftRoot] + 1
    end
end;

procedure AllocateArrays(ParticleCount: LongInt);
begin
  SetLength(X, ParticleCount + 1);
  SetLength(Y, ParticleCount + 1);
  SetLength(Z, ParticleCount + 1);

  SetLength(Parent, ParticleCount + 1);
  SetLength(Rank, ParticleCount + 1);
  SetLength(Root, ParticleCount + 1);
  SetLength(GroupSize, ParticleCount + 1);
  SetLength(GroupFirst, ParticleCount + 1);
  SetLength(GroupId, ParticleCount + 1);
  SetLength(GroupHead, ParticleCount + 1);
  SetLength(GroupTail, ParticleCount + 1);
  SetLength(GroupNext, ParticleCount + 1);
  SetLength(Order, ParticleCount + 1);
  SetLength(Work, ParticleCount + 1);
  SetLength(Members, ParticleCount + 1);
  SetLength(NextParticle, ParticleCount + 1);
  SetLength(CellX, ParticleCount + 1);
  SetLength(CellY, ParticleCount + 1);
  SetLength(CellZ, ParticleCount + 1);
  SetLength(UsedSlot, ParticleCount + 1);

  HashSize := NextPowerOfTwo(ParticleCount * 4 + 1024);
  SetLength(HashKey, HashSize + 1);
  SetLength(HashHead, HashSize + 1)
end;

procedure ReadParticles(Path: String);
var
  Tx, Ty, Tz: Double;
  Capacity: LongInt;
begin
  CurrentInputPath := Path;
  Assign(InputFile, Path);
  Reset(InputFile);

  Capacity := 1024;
  AllocateArrays(Capacity);
  N := 0;

  while not Eof(InputFile) do
  begin
    ReadLn(InputFile, Tx, Ty, Tz);
    N := N + 1;
    if N > Capacity then
    begin
      Capacity := Capacity * 2;
      SetLength(X, Capacity + 1);
      SetLength(Y, Capacity + 1);
      SetLength(Z, Capacity + 1)
    end;
    X[N] := Tx;
    Y[N] := Ty;
    Z[N] := Tz
  end;
  Close(InputFile);

  AllocateArrays(N)
end;

function HashSlot(Key: Int64; Create: Boolean): LongInt;
var
  Slot: LongInt;
begin
  Slot := LongInt(PositiveMod64(Key, HashSize)) + 1;
  while (HashKey[Slot] <> Key) and (HashKey[Slot] <> -1) do
  begin
    Slot := Slot + 1;
    if Slot > HashSize then
      Slot := 1
  end;

  if Create and (HashKey[Slot] = -1) then
  begin
    HashKey[Slot] := Key;
    NUsed := NUsed + 1;
    UsedSlot[NUsed] := Slot
  end;

  if HashKey[Slot] = Key then
    HashSlot := Slot
  else
    HashSlot := 0
end;

function Better(LeftGroup, RightGroup: LongInt): Boolean;
begin
  if GroupSize[LeftGroup] > GroupSize[RightGroup] then
    Better := True
  else if GroupSize[LeftGroup] < GroupSize[RightGroup] then
    Better := False
  else
    Better := GroupFirst[LeftGroup] <= GroupFirst[RightGroup]
end;

procedure SortGroups;
var
  Width, LeftEdge, Middle, RightEdge: LongInt;
  LeftPos, RightPos, OutPos, I: LongInt;
begin
  Width := 1;
  while Width < NGroups do
  begin
    LeftEdge := 1;
    while LeftEdge <= NGroups do
    begin
      Middle := LeftEdge + Width - 1;
      if Middle > NGroups then
        Middle := NGroups;
      RightEdge := LeftEdge + 2 * Width - 1;
      if RightEdge > NGroups then
        RightEdge := NGroups;
      LeftPos := LeftEdge;
      RightPos := Middle + 1;
      OutPos := LeftEdge;

      while (LeftPos <= Middle) and (RightPos <= RightEdge) do
      begin
        if Better(Order[LeftPos], Order[RightPos]) then
        begin
          Work[OutPos] := Order[LeftPos];
          LeftPos := LeftPos + 1
        end
        else
        begin
          Work[OutPos] := Order[RightPos];
          RightPos := RightPos + 1
        end;
        OutPos := OutPos + 1
      end;

      while LeftPos <= Middle do
      begin
        Work[OutPos] := Order[LeftPos];
        LeftPos := LeftPos + 1;
        OutPos := OutPos + 1
      end;

      while RightPos <= RightEdge do
      begin
        Work[OutPos] := Order[RightPos];
        RightPos := RightPos + 1;
        OutPos := OutPos + 1
      end;

      LeftEdge := RightEdge + 1
    end;

    for I := 1 to NGroups do
      Order[I] := Work[I];
    Width := Width * 2
  end
end;

procedure Initialise;
var
  I: LongInt;
begin
  for I := 1 to N do
  begin
    Parent[I] := I;
    Rank[I] := 0;
    Root[I] := 0;
    GroupSize[I] := 0;
    GroupFirst[I] := 0;
    GroupId[I] := 0;
    GroupHead[I] := 0;
    GroupTail[I] := 0;
    GroupNext[I] := 0;
    NextParticle[I] := 0
  end
end;

procedure BuildCells;
var
  I, Slot: LongInt;
  Key: Int64;
begin
  for I := 1 to HashSize do
  begin
    HashKey[I] := -1;
    HashHead[I] := 0
  end;

  NUsed := 0;
  for I := 1 to N do
  begin
    CellX[I] := Trunc(X[I] * NCell);
    CellY[I] := Trunc(Y[I] * NCell);
    CellZ[I] := Trunc(Z[I] * NCell);
    if CellX[I] = NCell then CellX[I] := NCell - 1;
    if CellY[I] = NCell then CellY[I] := NCell - 1;
    if CellZ[I] = NCell then CellZ[I] := NCell - 1;

    Key := Int64(CellX[I]) + Int64(NCell) *
           (Int64(CellY[I]) + Int64(NCell) * Int64(CellZ[I]));
    Slot := HashSlot(Key, True);
    NextParticle[I] := HashHead[Slot];
    HashHead[Slot] := I
  end
end;

procedure TestPair(I, J: LongInt);
var
  DX, DY, DZ, R2: Double;
begin
  DX := Abs(X[I] - X[J]);
  if DX > 0.5 * BoxSize then DX := BoxSize - DX;
  if DX <= LinkLength then
  begin
    DY := Abs(Y[I] - Y[J]);
    if DY > 0.5 * BoxSize then DY := BoxSize - DY;
    if DY <= LinkLength then
    begin
      DZ := Abs(Z[I] - Z[J]);
      if DZ > 0.5 * BoxSize then DZ := BoxSize - DZ;
      if DZ <= LinkLength then
      begin
        R2 := DX * DX + DY * DY + DZ * DZ;
        if R2 <= LinkLength2 then Join(I, J)
      end
    end
  end
end;

procedure FindLinks;
var
  U, I, J: LongInt;
  IX, IY, IZ, JX, JY, JZ: LongInt;
  KX, KY, KZ, ThisSlot, OtherSlot: LongInt;
  Key, Remainder: Int64;
begin
  for U := 1 to NUsed do
  begin
    ThisSlot := UsedSlot[U];
    Key := HashKey[ThisSlot];
    IZ := LongInt(Key div (Int64(NCell) * Int64(NCell)));
    Remainder := Key - Int64(IZ) * Int64(NCell) * Int64(NCell);
    IY := LongInt(Remainder div NCell);
    IX := LongInt(Remainder - Int64(IY) * NCell);

    for KZ := 0 to 1 do
      for KY := -1 to 1 do
        if (KZ <> 0) or (KY >= 0) then
          for KX := -1 to 1 do
            if (KZ <> 0) or (KY <> 0) or (KX >= 0) then
            begin
              JX := (IX + KX + NCell) mod NCell;
              JY := (IY + KY + NCell) mod NCell;
              JZ := (IZ + KZ + NCell) mod NCell;
              Key := Int64(JX) + Int64(NCell) *
                     (Int64(JY) + Int64(NCell) * Int64(JZ));
              OtherSlot := HashSlot(Key, False);

              if OtherSlot <> 0 then
              begin
                I := HashHead[ThisSlot];
                while I <> 0 do
                begin
                  if (KX = 0) and (KY = 0) and (KZ = 0) then
                    J := NextParticle[I]
                  else
                    J := HashHead[OtherSlot];

                  while J <> 0 do
                  begin
                    TestPair(I, J);
                    J := NextParticle[J]
                  end;
                  I := NextParticle[I]
                end
              end
            end
  end
end;

procedure MakeGroups;
var
  I, G: LongInt;
begin
  NGroups := 0;
  for I := 1 to N do
  begin
    Root[I] := FindRoot(I);
    G := GroupId[Root[I]];
    if G = 0 then
    begin
      NGroups := NGroups + 1;
      G := NGroups;
      GroupId[Root[I]] := G;
      GroupFirst[G] := I;
      GroupHead[G] := I;
      GroupTail[G] := I;
      Order[G] := G
    end
    else
    begin
      GroupNext[GroupTail[G]] := I;
      GroupTail[G] := I
    end;
    GroupSize[G] := GroupSize[G] + 1
  end;
  SortGroups
end;

procedure CountKeptGroups;
var
  I, G: LongInt;
begin
  NKeep := 0;
  for I := 1 to NGroups do
  begin
    G := Order[I];
    if GroupSize[G] >= MinKeep then
      NKeep := NKeep + 1
  end
end;

procedure Analyse;
begin
  LinkLength := B * Exp(Ln(BoxSize * BoxSize * BoxSize / N) / 3.0);
  LinkLength2 := LinkLength * LinkLength;
  NCell := Trunc(BoxSize / LinkLength);
  if NCell < 1 then NCell := 1;

  Initialise;
  BuildCells;
  FindLinks;
  MakeGroups;
  CountKeptGroups
end;

procedure WriteGroups;
var
  I, J, K, G, Count: LongInt;
  RX, RY, RZ, DX, DY, DZ, SX, SY, SZ: Double;
  CX, CY, CZ: Double;
begin
  if not DirectoryExists('outdir') then
    CreateDir('outdir');

  Assign(GroupFile, 'outdir/pascal_groups.txt');
  Rewrite(GroupFile);
  Assign(MemberFile, 'outdir/pascal_members.txt');
  Rewrite(MemberFile);

  WriteLn(GroupFile, '# group_id size com_x com_y com_z');
  WriteLn(MemberFile, '# group_id size member_ids...');

  NKeep := 0;
  for I := 1 to NGroups do
  begin
    G := Order[I];
    if GroupSize[G] >= MinKeep then
    begin
      NKeep := NKeep + 1;
      J := GroupFirst[G];
      RX := X[J]; RY := Y[J]; RZ := Z[J];
      SX := 0.0; SY := 0.0; SZ := 0.0;
      Count := 0;
      J := GroupHead[G];

      while J <> 0 do
      begin
        Count := Count + 1;
        Members[Count] := J;

        DX := WrapUnit(X[J] - RX + 0.5 * BoxSize);
        DY := WrapUnit(Y[J] - RY + 0.5 * BoxSize);
        DZ := WrapUnit(Z[J] - RZ + 0.5 * BoxSize);

        SX := SX + RX + DX - 0.5 * BoxSize;
        SY := SY + RY + DY - 0.5 * BoxSize;
        SZ := SZ + RZ + DZ - 0.5 * BoxSize;
        J := GroupNext[J]
      end;

      CX := WrapUnit(SX / GroupSize[G]);
      CY := WrapUnit(SY / GroupSize[G]);
      CZ := WrapUnit(SZ / GroupSize[G]);

      WriteLn(GroupFile, NKeep, ' ', GroupSize[G], ' ',
              CX:0:16, ' ', CY:0:16, ' ', CZ:0:16);
      Write(MemberFile, NKeep, ' ', GroupSize[G]);
      for K := 1 to Count do
        Write(MemberFile, ' ', Members[K]);
      WriteLn(MemberFile)
    end
  end;

  Close(GroupFile);
  Close(MemberFile)
end;

function ICPath(Power: LongInt): String;
var
  ParticleCount: LongInt;
begin
  ParticleCount := 1 shl Power;
  ICPath := Format('data/ics_%.7d.txt', [ParticleCount])
end;

function SecondsNow: Double;
begin
  SecondsNow := GetTickCount64 * 1e-3
end;

procedure TimeCascade;
var
  Power, Trial: LongInt;
  Path: String;
  StartTime, Elapsed: Double;
begin
  if not DirectoryExists('outdir') then
    CreateDir('outdir');

  Assign(TimingFile, 'outdir/pascal_fof_timings.txt');
  Rewrite(TimingFile);
  WriteLn(TimingFile, '# input_file n_particles trial analysis_seconds n_raw_groups n_kept_groups');

  for Power := FirstPower to LastPower do
  begin
    Path := ICPath(Power);
    if not FileExists(Path) then
    begin
      WriteLn('Skipping missing ', Path);
      Continue
    end;

    ReadParticles(Path);
    if N = 0 then
    begin
      WriteLn('Skipping empty ', Path);
      Continue
    end;

    for Trial := 1 to NTrials do
    begin
      StartTime := SecondsNow;
      Analyse;
      Elapsed := SecondsNow - StartTime;

      WriteLn(TimingFile, Path, ' ', N, ' ', Trial, ' ',
              Elapsed:0:8, ' ', NGroups, ' ', NKeep);
      Flush(TimingFile)
    end;

    WriteGroups;
    WriteLn('Timed ', Path, ' (', N, ' particles; ', NKeep, ' kept groups)')
  end;

  Close(TimingFile)
end;

begin
  TimeCascade
end.
