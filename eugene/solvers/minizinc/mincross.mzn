include "modelGenotypes.mzn";

var NLoci union {0}: obj;
constraint obj = count(treeType, Node);
constraint obj >= ceil((nGenotypes + nLoci) / 2);
constraint obj <= ceil(((nGenotypes+1)*nGenotypes/2 + nLoci) / 2);
constraint obj <= nLoci;
solve minimize obj;
