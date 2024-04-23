include "modelGenotypes.mzn";

var NLoci union {0}: obj;
constraint obj = count(treeType, Node);
constraint obj <= nLoci;
solve minimize obj;
