import math
import statistics

from ribonorma import stat_transformations as strf

def tpm (reads, gene_length):

    """
    Transcripts Per Million.
    """

    assert (len(reads) == len(gene_length))

    kappa = list()

    for i in range(len(reads)):
        read = reads[i]
        current_length = gene_length[i]

        a = (read * 1e3)/current_length
        kappa.append(a)

    kappa_sum = sum(kappa)
    normalised = list()

    for i in range(len(kappa)):
        tpm_count = kappa[i] / kappa_sum * 1e6
        normalised.append(tpm_count)

    return normalised

def fpkm (reads, gene_length):

    assert (len(reads) == len(gene_length))

    kappa = sum([reads[i] * gene_length[i] for i in range(len(reads))])
    normalised = list()

    for i in range(len(reads)):

        fpkm_count = (reads[i] * 1e9) / kappa
        normalised.append(fpkm_count)

    return normalised

def tpmr (samples, gene_length=[], experimental_conditions=[], percent_housekeep=10):

    assert (len(samples) == len(experimental_conditions))

    all_samples_tpm = list()

    for sample in samples:

        sample_tpm = tpm(sample, gene_length=gene_length)
        all_samples_tpm.append(sample_tpm)

    transposed_tpm = [list(x) for x in list(zip(*all_samples_tpm))]

    unique_conditions = list(dict.fromkeys(experimental_conditions))
    experimental_averaged_tpm = list()

    for i in range(len(transposed_tpm)):

        condition_average = [0] * len(unique_conditions)

        for j in range(len(experimental_conditions)):

            index = unique_conditions.index(experimental_conditions[j])

            condition_average[index] += transposed_tpm[i][j] / experimental_conditions.count(experimental_conditions[j])

        experimental_averaged_tpm.append(condition_average)

    geomeans = list()
    for k1 in range(len(unique_conditions)):

        condition_1 = unique_conditions[k1]
        condition_factor = 0

        ratio_means = list()

        for k2 in range(len(unique_conditions)):

            condition_2 = unique_conditions[k2]

            # Perform normalisation
            ratios = list()

            for gene_data in experimental_averaged_tpm:
                ratio = gene_data[k2] / (gene_data[k1] + 1e-128)

                ratios.append(ratio)

            trimmed_ratios = strf.trim_percentiles(ratios, [50 - percent_housekeep/2, 50 + percent_housekeep/2])

            # Mean values
            mean = strf.geometric_mean(trimmed_ratios)

            ratio_means.append(mean)

        geomean = strf.geometric_mean(ratio_means)
        geomeans.append(geomean)

    for j in range(len(samples)):

        # Find the sample type and geomean
        index = unique_conditions.index(experimental_conditions[j])
        norm_factor = geomeans[index]

        for i in range(len(samples[j])):
            all_samples_tpm[j][i] *= norm_factor

    return all_samples_tpm

    # Central limit theorems
    #sample_size = math.floor(0.1 * len(trimmed_ratios))
    #transform = strf.stochastic_transformation(trimmed_ratios, sample_size=sample_size, samples=300)

    #return geomeans

def tpmm (samples, gene_length):

    all_samples_tpm = list()

    for sample in samples:

        sample_tpm = tpm(sample, gene_length=gene_length)
        all_samples_tpm.append(sample_tpm)

    transposed_tpm = [list(x) for x in list(zip(*all_samples_tpm))]

    sample_medians = list()
    for sample_tpm in all_samples_tpm:

        # Sum of TPM is always 1e6
        median = statistics.median(sample_tpm)
        sample_medians.append(median)

    # Geometric mean normalisation to medians
    nfs = list()
    for j in range(len(sample_medians)):

        # Calculate geometric mean
        other_substituents = sample_medians[:j] + sample_medians[j+1:]
        norm_factor = strf.geometric_mean(other_substituents)

        nf = (norm_factor/sample_medians[j])
        nfs.append(nf)

    # De-adjust
    for j in range(len(sample_medians)):

        for i in range(len(all_samples_tpm[j])):
            all_samples_tpm[j][i] *= nfs[j]

    return all_samples_tpm

def tpmr_2 (samples, gene_length, experimental_conditions=[], percent_housekeep=10):

    assert (len(samples) == len(experimental_conditions))

    all_samples_tpm = tpmm(samples=samples, gene_length=gene_length)
    transposed_tpm = [list(x) for x in list(zip(*all_samples_tpm))]

    unique_conditions = list(dict.fromkeys(experimental_conditions))
    experimental_averaged_tpm = list()

    for i in range(len(transposed_tpm)):

        condition_average = [0] * len(unique_conditions)

        for j in range(len(experimental_conditions)):

            index = unique_conditions.index(experimental_conditions[j])

            condition_average[index] += transposed_tpm[i][j] / experimental_conditions.count(experimental_conditions[j])

        experimental_averaged_tpm.append(condition_average)

    geomeans = list()
    for k1 in range(len(unique_conditions)):

        condition_1 = unique_conditions[k1]
        condition_factor = 0

        ratio_means = list()

        for k2 in range(len(unique_conditions)):

            condition_2 = unique_conditions[k2]

            # Perform normalisation
            ratios = list()

            for gene_data in experimental_averaged_tpm:
                ratio = gene_data[k2] / (gene_data[k1] + 1e-128)

                ratios.append(ratio)

            trimmed_ratios = strf.trim_percentiles(ratios, [50 - percent_housekeep/2, 50 + percent_housekeep/2])

            # Mean values
            mean = strf.geometric_mean(trimmed_ratios)
            ratio_means.append(mean)

        geomean = strf.geometric_mean(ratio_means)
        geomeans.append(geomean)

    for j in range(len(samples)):

        # Find the sample type and geomean
        index = unique_conditions.index(experimental_conditions[j])
        norm_factor = geomeans[index]

        for i in range(len(samples[j])):
            all_samples_tpm[j][i] *= norm_factor

    return all_samples_tpm
