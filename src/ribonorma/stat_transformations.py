import statistics
import numpy as np
import math
import random

def trim_percentiles (data, percentiles=[5, 95]):

    # Two-tailed trimming
    lower = np.percentile(data, percentiles[0])
    upper = np.percentile(data, percentiles[1])

    return [x for x in data if x >= lower and x <= upper]

def stochastic_transformation (data, sample_size=3, samples=300):

    # Uses the central limit theorem
    output = list()
    for i in range(samples):

        choiced = random.choices(data, k=sample_size)
        average = sum(choiced) / sample_size

        output.append(average)

    return output

def geometric_mean (data):

    sum_of_logs = sum([math.log(x) for x in data]) / len(data)
    return math.exp(sum_of_logs)

def standard_error_of_mean (data):
    return statistics.variance(data) / math.sqrt(len(data))
