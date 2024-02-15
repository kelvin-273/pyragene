include "modelGenotypes.mzn";

var NLoci union {0}: nCrossings;
constraint nCrossings = count(treeType, Node);
%constraint nCrossings >= ceil((nGenotypes + nLoci) / 2);
%constraint nCrossings <= ceil(((nGenotypes+1)*nGenotypes/2 + nLoci) / 2);
constraint nCrossings <= nLoci;

array[NTreeCells] of var NLoci union {0}: gen;
constraint forall (i in NTreeCells) (
  (treeType[i] == Null -> gen[i] = 0) /\
  (treeType[i] == Leaf -> gen[i] = 0) /\
  (treeType[i] == Node -> gen[i] = max(gen[treeLeft[i]], gen[treeRight[i]]) + 1)
);

var int: obj;
constraint obj = nCrossings + 10 * 2^gen[1];
solve minimize obj;
