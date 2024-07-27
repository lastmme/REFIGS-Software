from collections import defaultdict
from multiprocessing import Queue

import numpy as np
from pyteomics import mzxml
from bisect import bisect
from timeit import default_timer as timer

from . import IdentificationSpectraMatcher
from . import spectra_matcher_functions as smf


def perform_spectra_pooling_and_analysis(
    querySpectraFile: str,
    outFile: str,
    lib: dict[tuple[str, int], dict],
    text_queue: Queue,
    tolerance=30,
    maxQuerySpectraToPool=np.inf,
    is_corrected: bool = True
):
    text_queue.put(
        smf.print_milestone('Begin Grouping Scans by m/z Windows:')
    )
    queWindowDict, queScanValuesDict = pool_scans_by_mz_windows(
        querySpectraFile)

    text_queue.put(
        'Number of Unpooled MS/MS Query Spectra: ' +
        str(len(queScanValuesDict))
    )

    text_queue.put(
        'Number of Pooled MS/MS Query Spectra/Mz Windows: ' +
        str(len(queWindowDict))
    )

    # To enhance the print experience, status prints will be given at intervals tailored to the number of identified windows.
    #  example: if there are 1-99 pooled query spectra, print statements are made after every pooled query spectra analysis is complete.
    #           if there are 100-999, print after every 10 pooled spectra. And so on.
    printFriendlyCounter = 100
    while printFriendlyCounter < len(queWindowDict):
        printFriendlyCounter *= 10
    printFriendlyCounter /= 100

    allLibKeys, libIdToKeyDict, libIdToDecoyDict = gather_library_metadata(lib)
    allSpectraMatches = IdentificationSpectraMatcher.IdentificationSpectraMatcher()
    numWindowsAnalyzed = 0

    prevtime = timer()
    text_queue.put(
        smf.print_milestone('Begin Pooled Spectra Analysis:')
    )
    with mzxml.read(querySpectraFile, use_index=True) as spectra:

        for precMz_win, scans in queWindowDict.items():
            top_mz = precMz_win[0] + precMz_win[1] / 2
            bottom_mz = precMz_win[0] - precMz_win[1] / 2
            libKeys = identify_lib_spectra_in_window(
                top_mz, bottom_mz, allLibKeys)
            if len(libKeys) == 0:
                continue
            pooledLibSpectra = pool_lib_spectra(lib, libKeys)
            pooledQueSpectra = []

            for i in range(len(scans)):
                scanNumber = scans[i]
                queSpectrum = spectra.get_by_id(scanNumber)
                pooledQueSpectra.append(
                    smf.format_spectra_for_pooling(queSpectrum, scanNumber)
                )

                if (i % maxQuerySpectraToPool == 0 and i != 0) or i == len(scans)-1:
                    AfterpooledQueSpectra = np.concatenate(
                        pooledQueSpectra, axis=0)
                    indices = np.lexsort(
                        (AfterpooledQueSpectra[:, 2], AfterpooledQueSpectra[:, 1], AfterpooledQueSpectra[:, 0]))
                    AfterpooledQueSpectra = AfterpooledQueSpectra[indices]
                    windowSpectraMatches = IdentificationSpectraMatcher.IdentificationSpectraMatcher()
                    windowSpectraMatches.compare_spectra(
                        pooledLibSpectra, AfterpooledQueSpectra, tolerance, libIdToDecoyDict)
                    allSpectraMatches.extend_all_spectra(windowSpectraMatches)
                    pooledQueSpectra.clear()

            numWindowsAnalyzed += 1
            # if numWindowsAnalyzed % printFriendlyCounter == 0:
            #     time = timer()
            #     # print('\nNumber of Pooled Experimental Spectra Analyzed: ' +
            #     #       str(numWindowsAnalyzed))
            #     # print('Number of Spectra in Current Pooled Spectra: ' +
            #     #       str(len(scans)))
            #     # print('Time Since Last Checkpoint: ' +
            #     #       str(round(time-prevtime, 2)) + ' Seconds', flush=True)
            #     prevtime = time

    if is_corrected:
        text_queue.put(
            smf.print_milestone('Begin corrected:')
        )
        maccCutoff = allSpectraMatches.find_score_fdr_cutoff()
        allSpectraMatches.filter_by_corrected_ppm_window(maccCutoff)

    text_queue.put(
        smf.print_milestone('\nBegin Writing to File: ')
    )
    allSpectraMatches.write_output(outFile.replace(
        '.csv', 'NoFilter.csv'), querySpectraFile, -1, queScanValuesDict, libIdToKeyDict, lib
    )
    return outFile.replace('.csv', 'NoFilter.csv')


def pool_scans_by_mz_windows(querySpectraFile):
    queWindowDict = defaultdict(list)
    queScanValuesDict = defaultdict(dict)

    with mzxml.read(querySpectraFile, use_index=True) as spectra:
        for spec in spectra:
            if 'precursorMz' not in spec:
                continue
            scan = spec['num']
            precMz = spec['precursorMz'][0]['precursorMz']
            windowWidth = spec['precursorMz'][0]['windowWideness']
            queWindowDict[precMz, windowWidth].append(scan)

            queScanValuesDict[scan]['precursorMz'] = precMz
            queScanValuesDict[scan]['windowWideness'] = windowWidth
            queScanValuesDict[scan]['peaksCount'] = spec['peaksCount']
            if 'compensationVoltage' in spec:
                CV = spec['compensationVoltage']
            else:
                CV = ''
            queScanValuesDict[scan]['CV'] = CV

    return dict(queWindowDict), dict(queScanValuesDict)


def identify_lib_spectra_in_window(top_mz, bottom_mz, sortedLibKeys):
    temp = sortedLibKeys[:]
    top_key = (top_mz, "z")
    bottom_key = (bottom_mz, "")

    i1 = bisect(temp, bottom_key)
    temp.insert(i1, bottom_key)
    i2 = bisect(temp, top_key)
    if i2-i1 == 1:
        return []
    temp.insert(i2, top_key)

    return temp[i1+1:i2]


def gather_library_metadata(lib):
    allLibKeys = sorted(lib.keys())
    libIdToKeyDict = {}
    libIdToDecoyDict = {}
    for key in allLibKeys:
        libIdToKeyDict[lib[key]['ID']] = key
        libIdToDecoyDict[lib[key]['ID']] = lib[key]['Decoy']
    return allLibKeys, libIdToKeyDict, libIdToDecoyDict


def pool_lib_spectra(lib, libKeys):
    final = np.concatenate([lib[key]['Peaks'] for key in libKeys])
    indices = np.lexsort((final[:, 2], final[:, 1], final[:, 0]))
    return final[indices]
