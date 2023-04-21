import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
from os import system


def my_newfile(outfile, infiles):
    return system("cat {} | sort -k 3 -n > {}".format(
        ' '.join(infiles),
        outfile))

my_newfile("/tmp/res.tsv", [
    "results/benchmark_modsim_2_2.tsv"
])

my_newfile("/tmp/res_improved_2.tsv", [
    # "results/benchmark_modsim_improved2_1.tsv",
    # "results/benchmark_modsim_improved2_2.tsv",
    "results/benchmark_modsim_improved3_2.tsv",
    "results/benchmark_modsim_improved3_3.tsv",
    "results/benchmark_modsim_improved3_4.tsv",
])

my_newfile("/tmp/res_plsgod.tsv", [
    "results/benchmark_modsim_plsgod.tsv",
    "results/benchmark_modsim_plsgod2.tsv",
    # "results/benchmark_modsim_plsgod3.tsv",
])

header = [
    'Type',
    'Algo',
    'Nloci',
    'Runtime_mean',
    'Runtime_std',
]

name_dict = {
    "Enum": "Naive1",
    "EnumDom": "Naive2",
    "Greedy": "Segment",
}

df = pd.read_table('/tmp/res.tsv', names=header)
# df2 = pd.read_table('results/benchmark_modsim_5.tsv', names=header)
df3 = pd.read_table('/tmp/res_improved_2.tsv', names=header)

# f = lambda df1: sn.scatterplot(x=df1['Nloci'], y=df1['Runtime_mean'])
# f = lambda df1: plt.plot(df1['Nloci'], df1['Runtime_mean'])
# f(df.loc[df['Algo'] == 'Enum'])
# f(df.loc[df['Algo'] == 'EnumDom'])
# f(df.loc[df['Algo'] == 'Greedy'])

# plt.title("Mean runtime (s) with respect to $n_l$")
# plt.xlabel("Number of loci $n_l$")
# plt.ylabel("Mean runtime (s)")
# plt.show()


def g(df1, label):
    plt.loglog(df1['Nloci'], df1['Runtime_mean'], label=label)


def my_plot(df: pd.DataFrame):
    for ty in df['Algo'].unique():
        g(df.loc[df['Algo'] == ty], name_dict[ty])


my_plot(df)
my_plot(df3)

# x = np.logspace(0, 4, 1001)
# plt.plot(x, x+1)

plt.title("Mean runtime (s) with respect to $n_l$")
plt.xlabel("Number of loci $n_l$")
plt.ylabel("Mean runtime (s)")
plt.legend()
plt.show()
