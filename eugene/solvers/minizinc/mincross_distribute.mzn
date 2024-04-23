include "mincross.mzn";

constraint obj >= ceil((nGenotypes + nLoci) / 2);
constraint obj <= ceil(((nGenotypes+1)*nGenotypes/2 + nLoci) / 2);
