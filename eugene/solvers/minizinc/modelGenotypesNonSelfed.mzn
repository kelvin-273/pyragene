include "globals.mzn";

int: nLoci;
int: nGenotypes;
set of int: NLoci = 1..nLoci;
set of int: NLociDiff = NLoci diff {1};
set of int: NGenotypes = 1..nGenotypes;
set of int: NGenotypesOpt = {0} union NGenotypes;
int: maxCrossovers;

int: nTreeCells;
set of int: NTreeCells = 1..nTreeCells;

enum NodeType = {Node, Leaf, Null};
array[NTreeCells] of var NodeType: treeType;
array[NTreeCells] of var NTreeCells union {0}: treeLeft;
array[NTreeCells] of var NTreeCells union {0}: treeRight;

% Tree structure
constraint increasing(treeType);
constraint forall (i in NTreeCells where treeType[i] == Node) (
  treeLeft[i] > i /\
  treeRight[i] > i /\
  treeType[treeLeft[i]] != Null /\
  treeType[treeRight[i]] != Null /\
  true
);
constraint forall (i in NTreeCells where treeType[i] != Node) (
  treeLeft[i] = 0 /\ treeRight[i] = 0
);

set of int: B = 0..1;
array[NGenotypes, 1..2, NLoci] of B: genotypes;
array[NTreeCells, 1..2, NLoci] of var B: xs;

array[NTreeCells] of var NGenotypesOpt: index;
constraint forall (i in NTreeCells) (
  (treeType[i] = Leaf <-> index[i] > 0) /\
  (treeType[i] = Null <-> forall (c in 1..2) (
    xs[i, c, NLoci] = [0 | _ in NLoci])
  )
);
constraint all_different_except_0(index);
constraint forall (i in NTreeCells where treeType[i] == Leaf) (
  forall (c in 1..2, j in NLoci) (
    xs[i, c, j] = genotypes[index[i], c, j]
  )
);

% First plant is the ideotype
predicate isIdeotype(NTreeCells: i) = ( 
  xs[i, 1, NLoci] == [1 | _ in NLoci] /\
  xs[i, 2, NLoci] == [1 | _ in NLoci] /\
  treeType[i] != Null
);
constraint isIdeotype(1);

% Each internal node and its child nodes (genetic parents in reality)
% must be related by crossing
constraint forall (i in NTreeCells where treeType[i] = Node) (
  crossedGenotype(
    treeLeft[i],
    treeRight[i],
    i,
  )
);
predicate crossedGenotype(
  var NTreeCells: igx,
  var NTreeCells: igy,
  var NTreeCells: igz,
) = (
  crossedGamete(igx, 1, igx, 2, igz, 1) /\
  crossedGamete(igy, 1, igy, 2, igz, 2)
);
predicate crossedGameteSPC(
  var NTreeCells: igx, var 1..2: cgx,
  var NTreeCells: igy, var 1..2: cgy,
  var NTreeCells: igz, var 1..2: cgz,
) = (
  let {
    var NLoci: k;
  } in (forall (j in NLoci where j <= k) (
    xs[igz, cgz, j] = xs[igx, cgx, j]
    ) /\ forall (j in NLoci where j > k) (
      xs[igz, cgz, j] = xs[igy, cgy, j]
    )) \/ (forall (j in NLoci where j <= k) (
    xs[igz, cgz, j] = xs[igy, cgy, j]
    ) /\ forall (j in NLoci where j > k) (
      xs[igz, cgz, j] = xs[igx, cgx, j]
    ))

);
predicate crossedGamete(
  var NTreeCells: igx, var 1..2: cgx,
  var NTreeCells: igy, var 1..2: cgy,
  var NTreeCells: igz, var 1..2: cgz,
) = (
  let {
    array[NLoci] of var 1..2: source;
    array[NLociDiff] of var 0..1: swap;

    constraint forall (j in NLoci diff {1}) (
      swap[j] = 1 <-> source[j-1] != source[j]
    );
    constraint sum(swap) <= maxCrossovers;
    constraint count(swap, 1) <= maxCrossovers;
    % dominance
  } in forall (j in NLoci) (
    xs[igz, cgz, j] = [xs[igx, cgx, j], xs[igy, cgy, j]][source[j]]
  )
);

constraint symmetry_breaking_constraint(
  forall (i in NTreeCells where treeType[i] == Node) (
    treeLeft[i] <= treeRight[i]
  )
);

% ideotype not in genotypes -> not selfing
constraint (
  forall (i in NGenotypes) (
    exists (c in 1..2, j in NLoci) (
      genotypes[i, c, j] == 0
    )
  )
) -> (treeLeft[1] = 2 xor treeRight[1] = 2);
constraint forall (i in NTreeCells where i > 1 /\ treeType[i] != Null) (
  not isIdeotype(i)
);
