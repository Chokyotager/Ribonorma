import scipy
from scipy import stats
from scipy import integrate
import numpy as np
import warnings
from ribonorma import stat_transformations as strf

warnings.filterwarnings("ignore")

def twosample_independent (sample1, sample2, genes):

    # Welch's t-test
    significances = list()

    for i in range(len(genes)):

        stat, p = stats.ttest_ind(sample1[i], sample2[i], equal_var=False)
        fold_change = np.mean(sample2[i]) / (np.mean(sample1[i]) + 1e-128)
        significances.append({"statistic": stat, "p-value": p, "gene": genes[i], "fold-change": fold_change})

    return sorted(significances, key=lambda x : x["p-value"])

def twosample_dependent (sample1, sample2, genes):

    # Wilcoxon signed ranks
    significances = list()

    for i in range(len(genes)):
        stat, p = stats.wilcoxon(sample1[i], sample2[i], correction=True)
        fold_change = (np.median(sample2[i]) + 1e-128) / (np.median(sample1[i]) + 1e-128)
        significances.append({"statistic": stat, "p-value": p, "gene": genes[i], "fold-change": fold_change})

    return sorted(significances, key=lambda x : x["p-value"])

def ANOVA (samples, genes):

    # Kruskal-Wallis test
    significances = list()

    for i in range(len(genes)):

        current_set = [x[i] for x in samples]

        means = [np.mean(x) + 1e-128 for x in current_set]
        geometric_mean = strf.geometric_mean(means)

        stat, p = stats.kruskal(*current_set)
        fold_change = np.mean([x/geometric_mean for x in means])

        significances.append({"statistic": stat, "p-value": p, "gene": genes[i], "fold-change": fold_change})

    return sorted(significances, key=lambda x : x["p-value"])

def calculateOVL (samples):

    def get_x_distribution (X, pdf):

        s = stats.tstd(X)
        x_bar = sum(X) / len(X)
        df = len(X) - 1

        if s == 0:
            s = 1e-128

        return lambda x : pdf(x, df=df, loc=x_bar, scale=s)

    functions = list()

    for sample in samples:
        functions.append(get_x_distribution(sample, stats.t.pdf))

    output = integrate.quad(lambda x : min([function(x) for function in functions]), -np.inf, np.inf, full_output=1)

    return output[0]
