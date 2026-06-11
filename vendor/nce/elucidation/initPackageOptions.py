from warnings import filterwarnings
import pandas as pd
import seaborn as sns

pd.options.mode.chained_assignment = None
pd.set_option("display.max_rows", 5000)
pd.set_option("display.max_colwidth", 100)
pd.set_option("display.max_columns", 120)
pd.set_option("display.width", 120)
pd.DataFrame.dropcol = lambda self, col: self.drop(col, axis=1)

sns.set(font_scale=1.5)

filterwarnings("ignore", category=FutureWarning, append=False)


def shut_up(module, line):
    filterwarnings("ignore", category=Warning, module=module, lineno=line, append=False)


shut_up("sklearn.metrics.regression", 543)
shut_up("statsmodels.nonparametric._kernel_base", 516)
shut_up("numpy.core._methods", 140)
shut_up("seaborn.distributions", 423)
shut_up("numpy.core.fromnumeric", 3118)
shut_up("matplotlib.colors", 479)
shut_up("numpy.lib.histograms", 824)
shut_up("numpy.lib.histograms", 825)
shut_up("statsmodels.stats.weightstats", 128)
shut_up("statsmodels.nonparametric.kernels", 128)
shut_up("sklearn.base", 329)
shut_up("sklearn.linear_model.logistic", 947)
shut_up("sklearn.linear_model.logistic", 1437)
# This one happens when there are too few data points to draw contours in a plot.
shut_up("matplotlib.contour", 0)
# This one is a FutureWarning.  I should probably investigate.
shut_up("scipy.stats.stats", 1713)
shut_up("sklearn.linear_model.logistic", 947)  # ... increase number of iterations.
# see https://github.com/numpy/numpy/pull/432
filterwarnings("ignore", "numpy.ufunc size changed")
