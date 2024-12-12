include "globals.mzn";

int: nLoci;
set of int: NLoci = 1..nLoci;
set of int: NLociDiff = NLoci diff {1};
set of int: B = 0..1;   % allele values
int: maxCrossovers;

int: nParents;
set of int: NParents = 1..nParents;
array[NParents, 1..2, NLoci] of B: parents;

int: nTreeCells;
set of int: NTreeCells = 1..nTreeCells;

enum NodeType = {Node, Leaf, Null};
array[NTreeCells] of var NodeType: treeType;
array[NTreeCells, 1..2, NLoci] of var B: treeData;
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

% Mapping leaves to parents surjectively
array[NTreeCells] of var {0} union NParents: index;
constraint forall (i in NTreeCells) (treeType[i] = Leaf <-> index[i] > 0);
constraint all_different_except_0(index);
constraint forall (i in NTreeCells where treeType[i] == Leaf) (
  forall (c in 1..2, j in NLoci) (
    treeData[i, c, j] = parents[index[i], c, j]
  )
);

constraint forall (i in NTreeCells) (
  (treeType[i] = Null <-> (
      treeData[i, 1, NLoci] = [0 | _ in NLoci] /\
      treeData[i, 2, NLoci] = [0 | _ in NLoci]
    )
  )
);

% First plant is the desired plant
constraint 
  treeData[1, 1, NLoci] == [1 | _ in NLoci] /\
  treeData[1, 2, NLoci] == [1 | _ in NLoci] /\
  treeType[1] != Null;

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
    treeData[igz, cgz, j] = treeData[igx, cgx, j]
    ) /\ forall (j in NLoci where j > k) (
      treeData[igz, cgz, j] = treeData[igy, cgy, j]
    )) \/ (forall (j in NLoci where j <= k) (
    treeData[igz, cgz, j] = treeData[igy, cgy, j]
    ) /\ forall (j in NLoci where j > k) (
      treeData[igz, cgz, j] = treeData[igx, cgx, j]
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
  } in forall (j in NLoci) (
    treeData[igz, cgz, j] = [treeData[igx, cgx, j], treeData[igy, cgy, j]][source[j]]
  )
);

% ideotype not in parents -> selfing
constraint forall (i in NParents) (
  sum(parents[i, 1..2, NLoci]) < 2*nLoci
) -> (treeLeft[1] = 2 /\ treeRight[1] = 2);

% lexicographical ordering between chromosomes
constraint symmetry_breaking_constraint(
  forall (i in NTreeCells where treeType[i] == Node /\ i > 1) (
    lex_less(treeData[i, 1, NLoci], treeData[i, 2, NLoci])
  )
);

% dominance between predecessors
constraint forall(iz in NTreeCells, ix in NTreeCells where treeType[iz] = Node /\ iz < ix) (
  forall(cz in 1..2, cx in 1..2) (
    treeData[iz, cz, NLoci] = treeData[ix, cx, NLoci] \/
    exists(j in NLoci) (
      treeData[iz, cz, j] > treeData[ix, cx, j]
    )
  )
);
% dominance between chromosomes
constraint forall(i in NTreeCells where treeType[i] = Node /\ i > 1) (
  exists(j1 in NLoci, j2 in NLoci where j1 != j2) (
    treeData[i, 1, j1] > treeData[i, 2, j1] /\
    treeData[i, 2, j2] > treeData[i, 1, j2]
  )
)
