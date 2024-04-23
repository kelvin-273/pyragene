import json
import matplotlib.pyplot as plt
import pandas as pd

import eugene.solvers.base_min_crossings_wedges as ew
import eugene.db as ed

with open("./distribute_data_2.json") as f:
    db = ed.DistributeDB(json.load(f))

n_locis = [len(eval(case)) for case in db.db]
opts = [db.get_base_solution(case).objective for case in db.db]
heus = [ew.breeding_program(len(eval(case)), eval(case)).objective for case in db.db]

df = pd.DataFrame({"n_loci": n_locis, "opt": opts, "heu": heus})

f, (ax0, ax1) = plt.subplots(2, 1, sharex=True)

# plots of the optimal and heuristic value
opts_μ = [df.opt[df.n_loci == n_loci].mean() for n_loci in df.n_loci.unique()]
opts_σ = [df.opt[df.n_loci == n_loci].std() for n_loci in df.n_loci.unique()]

heus_μ = [df.heu[df.n_loci == n_loci].mean() for n_loci in df.n_loci.unique()]
heus_σ = [df.heu[df.n_loci == n_loci].std() for n_loci in df.n_loci.unique()]

ax0.set_title("$cross(T)$ values from optimal and 2-approximation algorithms")
ax0.set_xlabel("$n_l$")
ax0.set_ylabel("$cross(T)$")
ax0.errorbar(df.n_loci.unique(), opts_μ, yerr=opts_σ, fmt='o', alpha=0.7, capsize=3, label="optimal")
ax0.errorbar(df.n_loci.unique(), heus_μ, yerr=heus_σ, fmt='o', alpha=0.7, capsize=3, label="2-approx")
ax0.legend(loc="lower right")
ax0.grid(alpha=0.3)

# scatter of the relative error
errs_μ = [(df.heu / df.opt)[df.n_loci == n_loci].mean() for n_loci in df.n_loci.unique()]
errs_σ = [(df.heu / df.opt)[df.n_loci == n_loci].std() for n_loci in df.n_loci.unique()]
ax1.set_title("Average relative error of 2-approximation for $n_l$ loci")
ax1.set_xlabel("$n_l$")
ax1.set_ylabel("relative error of $cross(T)$")
ax1.errorbar(df.n_loci.unique(), errs_μ, yerr=errs_σ, fmt='o', alpha=0.7, capsize=3)
ax1.grid(alpha=0.3)

# SHOWTIME!
f.show()
