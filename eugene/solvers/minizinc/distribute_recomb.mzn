include "globals.mzn";

int: nLoci;
set of int: NLoci = 1..nLoci;

int: outMax = nLoci-1;
int: nPop;
set of int: NIn = 0..nPop-1;
set of int: NOut = 0..outMax;
array[NLoci] of NIn: instance;

var NIn: gx;
var NIn: gy;
array[NOut] of var 0..nLoci: recLocus::no_output;
array[NLoci] of var NOut: xs;
array[NLoci] of var NOut: ys::no_output;
%array[NLoci] of var NOut: ysFreeMax;

var 1..outMax: nGametes;
var 1..nLoci: nSegments;

array[NOut] of var NIn: strtSource::no_output;
array[NOut] of var NIn: stopSource::no_output;
array[NOut] of var bool: switchesSource::no_output;
array[NOut] of var NLoci: strt::no_output;
array[NOut] of var NLoci: stop::no_output;
array[NOut] of var {0} union NLoci: preRec::no_output;
array[NOut] of var {0} union NLoci: pstRec::no_output;


% assert instance is a distribute instance
constraint instance[1] = 0;
constraint forall (i in 2..nLoci) (
  instance[i] != instance[i-1]
);
constraint forall (i in 2..nLoci) (
  (instance[i] > 0) -> member(instance[1..i-1], instance[i] - 1)
);

% Choose gametes to recombine
constraint gx < gy;
predicate isgxgy(var NLoci: i) = instance[i] in {gx, gy};

constraint nGametes = card({xs[i] | i in NLoci});
constraint nSegments = nLoci - sum(i in 2..nLoci)(xs[i] == xs[i-1]);

any: my_arr = [gx, gy] ++ [j | j in nPop..(nPop + nLoci -2)];

% free alleles owned by gx and gy, and fix everything else
constraint forall (i in NLoci) (
  if isgxgy(i)
  then xs[i] == my_arr[ys[i]]
  else xs[i] == instance[i] /\ ys[i] == 0
  endif
);
constraint seq_precede_chain(ys);

% constrain strt and stop
constraint forall (i in NLoci) (
  forall([
    strt[xs[i]] <= i,
    stop[xs[i]] >= i,
    xs[strt[xs[i]]] == xs[i],
    xs[stop[xs[i]]] == xs[i],
  ])
);

% Recombinations: recLocus and switchesSource

% for any for any output gamete any pair of alleles sourced from different
% input gametes should have a recombination point between their respective loci
constraint forall (i in NLoci, j in NLoci where i < j) (
  (xs[i] == xs[j] /\ instance[i] != instance[j]) -> (
    i <= recLocus[xs[i]] /\ recLocus[xs[i]] < j
  )
);

% symmetry breaking
constraint forall (i in nGametes..outMax) (
  recLocus[i] = 0 /\
  true
);

%%%%%%%%%%%%%%%%%%%%%%%%
%%  Dominance pruning  %
%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1. Any redistribution `ys` whose gametes are all subsets of gametes in
% `instance` can be pruned
constraint xs != instance;
% TODO: the full version of this constraint

% 2. Given redistribution `ys`, if there exists another redistribution `zs`
% such that `ys` gametes are all subsets of gametes in `zsz, then `ys` can be
% pruned

% constrain strtSource and stopSource to be the gametes
%   that xs[i] started stoped being sourced from
constraint forall (i in NLoci) (
  strtSource[xs[i]] = instance[strt[xs[i]]] /\
  stopSource[xs[i]] = instance[stop[xs[i]]]
);
% switchesSource[gx] is whether gx is created by switching source at some point
constraint forall (i in NLoci) (
  switchesSource[xs[i]] = (strtSource[xs[i]] != stopSource[xs[i]])
);

% preRec[gx] and pstRec[gx]
constraint forall (i in NLoci) (
  switchesSource[xs[i]] = (preRec[xs[i]] > 0) /\
  switchesSource[xs[i]] = (pstRec[xs[i]] > 0) /\
  true
);
constraint forall (i in NLoci where forall([
  switchesSource[xs[i]]
])) (
  forall([
    instance[preRec[xs[i]]] = strtSource[xs[i]],
    instance[pstRec[xs[i]]] = stopSource[xs[i]],
    preRec[xs[i]] <= recLocus[xs[i]],
    pstRec[xs[i]] >= recLocus[xs[i]],
    xs[preRec[xs[i]]] = xs[i],
    xs[pstRec[xs[i]]] = xs[i],
    (i <= recLocus[xs[i]] -> preRec[xs[i]] >= i), % should the first <= be < ?
    (i > recLocus[xs[i]] -> pstRec[xs[i]] <= i),
  ])
);
%% testing an assertion that should be true given that the previous constraints
%% are correct
%constraint forall (i in NLoci) (
%  switchesSource[xs[i]] -> (
%    forall([
%      strt[xs[i]] <= preRec[xs[i]],
%      preRec[xs[i]] <= recLocus[xs[i]],
%      recLocus[xs[i]] <= pstRec[xs[i]],
%      pstRec[xs[i]] <= stop[xs[i]]
%    ])
%  )
%);

% There shall be no two gametes that can start from the same source,
% end from the same source, and share a recombination point
constraint forall (i in NLoci, j in NLoci where i < j) (
  forall([
    xs[i] != xs[j],
    switchesSource[xs[i]],
    switchesSource[xs[j]],
  ]) -> (
    not forall([
      strtSource[xs[i]] == strtSource[xs[j]],
      stopSource[xs[i]] == stopSource[xs[j]],
      max(preRec[xs[i]], preRec[xs[j]]) < min(pstRec[xs[i]], pstRec[xs[j]]),
    ])
  )
);

% There shall be no two gametes that don't require recombination and share the
% same source
constraint forall (i in NLoci, j in NLoci where i < j) (
  forall([
    xs[i] != xs[j],
    not switchesSource[xs[i]],
    not switchesSource[xs[j]],
  ]) -> (
    not forall([
      strtSource[xs[i]] == stopSource[xs[i]],
      strtSource[xs[j]] == stopSource[xs[j]],
      strtSource[xs[i]] == stopSource[xs[j]],
    ])
  )
);

% There shall be no two single-source gametes whose (strt, stop) ranges do not
% overlap
constraint forall (i in NLoci, j in NLoci where i < j) (
  forall([
    xs[i] != xs[j],
    not switchesSource[xs[i]],
    not switchesSource[xs[j]],
    instance[i] in {gx, gy},
    instance[j] in {gx, gy},
  ]) -> (
    max(strt[xs[i]], strt[xs[j]]) <= min(stop[xs[i]], stop[xs[j]])
  )
);

% There shall be no two gametes where
% one gamete requires crossover,
% one gamete is single source
% and the single source can fit in the range of that source provided by the crossover
constraint forall (i in NLoci, j in NLoci) (
  forall([
    xs[i] != xs[j],
    switchesSource[xs[i]],
    not switchesSource[xs[j]],
    instance[i] in {gx, gy},
    instance[j] in {gx, gy},
  ]) -> (
    % single-source is on the left
    not (
      stop[xs[j]] < pstRec[xs[i]] /\
      strtSource[xs[i]] == strtSource[xs[j]]
    ) /\
    % single-source is on the right
    not (
      preRec[xs[i]] < strt[xs[j]] /\
      stopSource[xs[i]] == stopSource[xs[j]]
    )
  )
);

%solve minimize nGametes + nSegments;
%constraint nGametes + nSegments < nPop + nLoci;

output [
  %"my_arr: ", show(my_arr),
  %"\nin: ", show(instance),
  "xs: ", show(xs),
  " ", show((gx, gy)),
  ": ", show((nPop, nLoci)),
  " -> ", show((nGametes, nSegments)),
  %"\nys: ", show(ys),
  %"\nym: ", show(ysFreeMax),
  %"\ngn: ", show([j | j in {gx, gy} union nPop..nPop + nLoci - 2]),
  %"\ngn[ys[2]]: ", show([j | j in {gx, gy} union nPop..nPop + nLoci - 2][ys[2]]),
  %"\nys > 0 <-> isgxgy(ys): ", show([ys[i] > 0 <-> isgxgy(i) | i in NLoci])
];

%constraint xs = [0,0,1,2,1]
%        \/ xs = [0,1,0,0,2]
%        \/ xs = [0,1,0,1,0]
%        \/ xs = [0,1,0,2,1]
%        \/ xs = [0,1,0,2,2]
%        \/ xs = [0,1,1,2,0]
%        \/ xs = [0,1,1,2,1]
%        \/ xs = [0,1,2,2,0];
%constraint xs = [0,1,0,2,0]; % should be unsat
%constraint xs = [0,1,0,1,0]
%        /\ ys = [0,1,0,1,0];

%output [ show(xs) ];

%constraint gx = 0;
%constraint gy = 5;
%constraint xs = [0, 1, 0, 2, 3, 6, 6, 1, 4, 1, 0, 5, 0];
%constraint xs = [0, 1, 0, 2, 3, 2, 6, 1, 4, 1, 6, 5, 0];

%constraint gx = 2;
%constraint gy = 3;
