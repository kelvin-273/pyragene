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
array[NLoci] of var NOut: xs;
array[NOut] of var 0..nLoci: recLocus::no_output;

array[NLoci] of var NOut: max_instance::no_output;
%array[NLoci] of var NOut: max_xs;

var 1..outMax: nGametes;
var 1..nLoci: nSegments;

array[NOut] of var NIn: startSource::no_output;
array[NOut] of var NIn: stopSource::no_output;
array[NOut] of var bool: switchesSource::no_output;
array[NOut] of var NLoci: start::no_output;
array[NOut] of var NLoci: stop::no_output;
array[NOut] of var NLoci: preRec::no_output;
array[NOut] of var NLoci: postRec::no_output;


% distribute instance
constraint instance[1] = 0;
constraint forall (i in 2..nLoci) (
  instance[i] != instance[i-1]
);

% each element in instance should be as at most one more than any of the
% previously seen elements
constraint max_instance[1] = instance[1];
constraint forall (i in 2..nLoci) (
  max_instance[i] = max(max_instance[i-1], instance[i]) /\
  instance[i] <= max_instance[i-1] + 1
);

% Choose gametes to recombine
constraint gx < gy;

% fix 0s and free 1s and 2s
constraint forall (i in NLoci) (
  (instance[i] != gx /\ instance[i] != gy) -> (xs[i] == instance[i])
);
constraint forall (i in NLoci) (
  (instance[i] == gx \/ instance[i] == gy) -> (
    xs[i] in {gx, gy} union nPop..nLoci
  )
);
constraint forall (i in NLoci) (
  xs[i] > 0 -> member(xs, xs[i] - 1)
);
constraint forall (i in NLoci, j in NLoci where instance[i] in {gx, gy} /\ instance[j] in {gx, gy}) (
  start[xs[i]] < start[xs[j]] <-> xs[i] < xs[j]
);

predicate distribute_structure(array[$$X] of var int: zs) = (
  zs[1] = 0 /\ forall (i in NLoci) (
    zs[i] > 0 -> member(zs[1..i-1], zs[i] - 1)
  )
);
constraint distribute_structure(xs);

% each element in xs should be as at most one more than any of the
% previously seen elements
constraint xs[1] = 0;

constraint nGametes = max(xs) + 1;
constraint nSegments = nLoci - sum(i in 2..nLoci)(xs[i] == xs[i-1]);
constraint forall (i in nGametes..outMax) (
  recLocus[i] = 0
);

% for any for any output gamete any pair of alleles sourced from different
% input gametes should have a recombination point between their respective loci
constraint forall (i in NLoci, j in NLoci where i < j) (
  (xs[i] == xs[j] /\ instance[i] != instance[j]) -> (
    i <= recLocus[xs[i]] /\ recLocus[xs[i]] < j
  )
);

%%%%%%%%%%%%%%%%%%%%%%%
%  Dominance pruning  %
%%%%%%%%%%%%%%%%%%%%%%%

% 1. Any redistribution `ys` whose gametes are all subsets of gametes in
% `instance` can be pruned
constraint xs != instance;
constraint forall (i in NLoci where instance[i] in {gx, gy}) (
  true
);

% 2. Given redistribution `ys`, if there exists another redistribution `zs`
% such that `ys` gametes are all subsets of gametes in `zsz, then `ys` can be
% pruned

%% constrain start and stop
%constraint forall (i in NLoci) (
%  forall([
%    start[xs[i]] <= i,
%    stop[xs[i]] >= i,
%    xs[start[xs[i]]] == xs[i],
%    xs[stop[xs[i]]] == xs[i],
%  ])
%);
%% constrain startSource and stopSource
%constraint forall (i in NLoci) (
%  startSource[xs[i]] = instance[start[xs[i]]] /\
%  stopSource[xs[i]] = instance[stop[xs[i]]]
%);
%constraint forall (i in NLoci) (
%  switchesSource[xs[i]] = (startSource[xs[i]] != stopSource[xs[i]])
%);
%
%% preRec and postRec
%constraint forall (i in NLoci) (
%  forall([
%    switchesSource[xs[i]]
%  ]) -> (
%    forall([
%      instance[preRec[xs[i]]] = startSource[xs[i]],
%      preRec[xs[i]] <= recLocus[xs[i]],
%      xs[preRec[xs[i]]] = xs[i],
%    ])
%  )
%);
%constraint forall (i in NLoci) (
%  forall([
%    switchesSource[xs[i]]
%  ]) -> (
%    forall([
%      instance[postRec[xs[i]]] = stopSource[xs[i]],
%      postRec[xs[i]] >= recLocus[xs[i]],
%      xs[postRec[xs[i]]] = xs[i],
%    ])
%  )
%);
%
%% testing an assertion that should be true given that the previous constraints
%% are correct
%%constraint forall (i in NLoci) (
%%  switchesSource[xs[i]] -> (
%%    forall([
%%      start[xs[i]] <= preRec[xs[i]],
%%      preRec[xs[i]] <= recLocus[xs[i]],
%%      recLocus[xs[i]] <= postRec[xs[i]],
%%      postRec[xs[i]] <= stop[xs[i]]
%%    ])
%%  )
%%);
%
%% There shall be no two gametes that can start from the same source,
%% end from the same source, and share a recombination point
%constraint forall (i in NLoci, j in NLoci where i < j) (
%  forall([
%    xs[i] != xs[j],
%    switchesSource[xs[i]],
%    switchesSource[xs[j]],
%  ]) -> (
%    not forall([
%      startSource[xs[i]] == startSource[xs[i]],
%      stopSource[xs[i]] == stopSource[xs[i]],
%      startSource[xs[i]] != stopSource[xs[i]],
%      max(preRec[xs[i]], preRec[xs[j]]) < min(postRec[xs[i]], postRec[xs[j]]),
%    ])
%  )
%);
%
%% There shall be no two gametes that don't require recombination and share the
%% same source
%constraint forall (i in NLoci, j in NLoci where i < j) (
%  forall([
%    xs[i] != xs[j],
%    not switchesSource[xs[i]],
%    not switchesSource[xs[j]],
%  ]) -> (
%    not forall([
%      startSource[xs[i]] == stopSource[xs[i]],
%      startSource[xs[j]] == stopSource[xs[j]],
%      startSource[xs[i]] == stopSource[xs[j]],
%    ])
%  )
%);
%
%% There shall be no two single-source gametes whose (start, stop) ranges do not
%% overlap
%constraint forall (i in NLoci, j in NLoci where i < j) (
%  forall([
%    xs[i] != xs[j],
%    not switchesSource[xs[i]],
%    not switchesSource[xs[j]],
%    instance[i] in {gx, gy},
%    instance[j] in {gx, gy},
%  ]) -> (
%    max(start[xs[i]], start[xs[j]]) <= min(stop[xs[i]], stop[xs[j]])
%  )
%);
%
%% There shall be no two gametes where one gamete requires crossover, one gamete is single source and the single source can fit in the range of that source provided by the crossover
%constraint forall (i in NLoci, j in NLoci) (
%  forall([
%    xs[i] != xs[j],
%    switchesSource[xs[i]],
%    not switchesSource[xs[j]],
%    instance[i] in {gx, gy},
%    instance[j] in {gx, gy},
%  ]) -> (
%    % single-source is on the left
%    not (
%      stop[xs[j]] < postRec[xs[i]] /\
%      startSource[xs[i]] == startSource[xs[j]]
%    ) /\
%    % single-source is on the right
%    not (
%      preRec[xs[i]] < start[xs[j]] /\
%      stopSource[xs[i]] == stopSource[xs[j]]
%    )
%  )
%);

%array[NLoci] of var NOut: ys::no_output;
%array[NLoci] of var NLoci: ysIndex::no_output;

%% create an array that simplifies the representation of `xs`.
%% when this is the output, gecode will remove solutions that are equivalent.
%constraint ysIndex[1] = 1;
%constraint forall (j in NLoci where j > nSegments) (
%  ysIndex[j] = max(NLoci) /\ ys[j] = 0
%);
%constraint forall (j in NLoci where j <= nSegments) (
%  ys[j] = xs[ysIndex[j]]
%);
%constraint forall (j in NLoci where j < nSegments) (
%  forall([
%    ysIndex[j] > ysIndex[j+1],
%    ys[j] != ys[j+1]
%  ])
%);

%solve minimize nGametes + nSegments;
%constraint nGametes + nSegments < nPop + nLoci;

%output [
%  show(xs),
%  " ", show((gx, gy)),
%  ": ", show((nPop, nLoci)),
%  " -> ", show((nGametes, nSegments)),
%];

%output [ show(xs) ];

%constraint gx = 0;
%constraint gy = 5;
%constraint xs = [0, 1, 0, 2, 3, 6, 6, 1, 4, 1, 0, 5, 0];
%constraint xs = [0, 1, 0, 2, 3, 2, 6, 1, 4, 1, 6, 5, 0];
