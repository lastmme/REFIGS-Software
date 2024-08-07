import numpy as np

from numba import njit

from timeit import default_timer as timer
from datetime import timedelta


@njit
def approx(x, y, ppmTol):
    if x == y:
        return 1e-7
    ppmDiff = ((x-y)*1000000)/x
    return (ppmDiff if abs(ppmDiff) < ppmTol else 0)


@njit
def find_matching_peaks(libMzs, libIntensities, libTags, queMzs, queIntensities, queTags, ppmTol):
    lenLib = len(libMzs)
    lenQue = len(queMzs)
    matchLibTags = []
    matchLibIntensities = []
    matchQueTags = []
    matchQueIntensities = []
    ppmMatches = []
    i, j = 0, 0
    while i < lenLib and j < lenQue:
        if not approx(libMzs[i], queMzs[j], ppmTol):
            if libMzs[i] > queMzs[j]:
                j += 1
                continue
            if libMzs[i] < queMzs[j]:
                i += 1
                continue
        p = i + 0
        while (p < lenLib):
            ppm = approx(libMzs[p], queMzs[j], ppmTol)
            if p == lenLib or not ppm:
                break
            matchLibTags.append(libTags[p])
            matchLibIntensities.append(libIntensities[p])
            matchQueTags.append(queTags[j])
            matchQueIntensities.append(queIntensities[j])
            ppmMatches.append(ppm)
            p += 1
        j += 1
    return matchLibTags, matchLibIntensities, matchQueTags, matchQueIntensities, ppmMatches


def find_offset_tolerance(data, stdev, mean=True):
    if len(data) == 0:
        return 0, 10
    hist, bins = np.histogram(data, bins=200)
    center = (bins[:-1] + bins[1:]) / 2

    # offset is calculated as the mean or median value of the provided data, though this is overwritten if no corrected standard deviation is provided
    if mean:
        offset = sum(data)/len(data)
    else:
        offset = data[len(data)//2]

    # If a corrected standard deviation is provided, it is used to set the tolerance. If not, see below conditional statement
    tolerance = np.std(data)*stdev

    # If no corrected standard deviation is provided, one is customized.
    #  Offset is considered the highest point of the histogram,
    #  tolerance the range necessary before the peaks are the same size as the noise on either end.
    if not stdev:
        index_max = max(range(len(hist)), key=hist.__getitem__)
        noise = np.mean(hist[:10] + hist[-10:])
        min_i = 0
        max_i = len(hist)-1
        for i in range(index_max, 0, -1):
            if hist[i] < noise:
                min_i = i
                break
        for i in range(index_max, len(hist)):
            if hist[i] < noise:
                max_i = i
                break

        offset = center[index_max]
        if index_max - min_i >= max_i - index_max:
            tolerance = offset - center[min_i]
        else:
            tolerance = center[max_i] - offset
    return offset, tolerance


def print_milestone(text):
    return f'{text}\n{str(timedelta(seconds=timer()))}\n'


def format_spectra_for_pooling(spectrum, scanNumber, sqrt=True):
    scanNumber = int(scanNumber)
    if sqrt:
        intensity = np.power(spectrum['intensity array'], 0.5)
    else:
        intensity = spectrum['intensity array']
    peakIDs = np.array([scanNumber] * spectrum['peaksCount'])
    mz = np.expand_dims(spectrum['m/z array'], axis=1)
    intensity = np.expand_dims(intensity, axis=1)
    ids = np.expand_dims(peakIDs, axis=1)
    return np.concatenate((mz, intensity, ids), axis=1)


def calculate_heavy_mz(seq, mz, z):
    hK = 8.014199  # mass of heavy lysine
    hR = 10.00827  # mass of heavy arg

    nK = seq.count('K')
    nR = seq.count('R')
    heavyMz = mz + (nK*hK)/z + (nR*hR)/z
    return heavyMz


def approx_list(x, l, ppmTol=10):
    for i in range(len(l)):
        if approx(x, l[i], ppmTol):
            return i
    return -1
