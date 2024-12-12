include "mincross.mzn";

constraint obj >= ceil((nParents+ nLoci) / 2);
constraint obj <= ceil(((nParents+1)*nParents/2 + nLoci) / 2);
