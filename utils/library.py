import copy
from typing import Sequence
from multiprocessing import Queue

import numpy as np
import numpy.typing as npt
import pandas as pd

from REFIGS.myms import GenerateDecoyLibrary
from CsoDIAq import spectra_matcher_functions as smf


def transform_cosdiq_format(libs: dict):
    temp_dict = {}
    for key, metadata in libs.items():
        sequence = key[0]
        id = metadata['ID']
        spectrum = metadata['Spectrum']
        spectrum[:, 1] = np.power(spectrum[:, 1], 0.5)
        new_key = (metadata['PrecursorMZ'], sequence)
        temp_dict[new_key] = {
            'PrecursorCharge': metadata['PrecursorCharge'],
            'transition_group_id': metadata['transition_group_id'],
            'Peaks': np.concatenate((spectrum, np.array([[id] for _ in spectrum[:, 0]])), axis=1),
            'ProteinName': '',
            'ID': metadata['ID'],
            'Decoy': metadata['Decoy']
        }
    return temp_dict


def entrapment(
    library: dict[tuple[str, int], dict],
    distance: int,
    text_queue: Queue
):
    text_queue.put(
        smf.print_milestone(
            f"{'#' * 15}  generating entrapment library...  {'#' * 15}")
    )
    decoy_library = GenerateDecoyLibrary(
        library, distance, text_queue
    )
    entrap_library = {}
    seq: str
    for (seq, charge), metadata in decoy_library.items():
        entrap_library[(seq.replace('DECOY-', 'TRAP-'), charge)] = metadata
    return entrap_library


def decoy(
    library: dict[tuple[str, int], dict],
    distance: int,
    text_queue: Queue
):
    text_queue.put(
        smf.print_milestone(
            f"{'#' * 15}  generating decoy library...  {'#' * 15}")
    )
    decoy_library = GenerateDecoyLibrary(
        library, distance, text_queue)
    return decoy_library


def entrapment_decoy_strategy(
    library: dict[tuple[str, int], dict],
    entrap_distance: int,
    decoy_distance: int,
    text_queue: Queue
):
    entrap_library = entrapment(
        library, entrap_distance, text_queue)
    library.update(entrap_library)
    decoy_library = decoy(
        library, decoy_distance, text_queue)
    library.update(decoy_library)
    return library


def decoy_strategy(
    library: dict[tuple[str, int], dict],
    decoy_distance: int,
    text_queue: Queue
):
    decoy_library = decoy(
        library, decoy_distance, text_queue)
    library.update(decoy_library)
    return library


def mgf2pkl(path: str):
    with open(path, 'r') as f:
        content = f.read()
        entries = content.split('END IONS\n\n')

        target = {}
        decoy = {}
        for id, entry in enumerate(entries, 1):
            if len(entry) < 5:
                continue
            else:
                title = entry.split('TITLE=')[1].split(' ')[0]
                charge = int(entry.split('CHARGE=')[1].split('\n')[0][0])
                mz = float(entry.split('PEPMASS=')[1].split('\n')[0])
                seq = entry.split('SEQ=')[1].split('\n')[0]
                spectrum_string = entry.split('=')[-1].split('\n')[1:-1]
                spectrum = np.array([peak.split(' ')
                                    for peak in spectrum_string]).astype('float')
                spectrum = spectrum[np.argsort(spectrum[:, 0])]
                if 'DECOY' in title:
                    seq = f'DECOY-{seq}'
                    is_decoy = 1
                else:
                    is_decoy = 0
                key = (seq, charge)
                tmp_dict = {
                    'PrecursorCharge': charge,
                    'transition_group_id': title,
                    'Spectrum': spectrum,
                    'PrecursorMZ': mz,
                    'ID': id,
                    'Decoy': is_decoy
                }
                tmp_dict['PrecursorMZ'] = mz
                tmp_dict['Spectrum'] = spectrum
                if 'DECOY' in title:
                    decoy[key] = tmp_dict
                else:
                    target[key] = tmp_dict
        target.update(decoy)
    return target


def extend_library(
    mgf_path: str,
    entrap_distance: int,
    decoy_distance: int,
    text_queue: Queue,
    is_entrap: bool = True,
    is_top10: bool = True
):
    text_queue.put(
        smf.print_milestone(f"{'#'*15}  Loading library...  {'#'*15}")
    )
    library = mgf2pkl(mgf_path)
    text_queue.put(
        smf.print_milestone(f"{'#'*15}  End  {'#'*15}")
    )
    if is_entrap:
        library = entrapment_decoy_strategy(
            library, entrap_distance, decoy_distance, text_queue)
    else:
        library = decoy_strategy(library, decoy_distance, text_queue)
    text_queue.put(
        smf.print_milestone(f"{'#'*15}  End *** decoy library {'#'*15}")
    )
    if is_top10:
        for metadata in library.values():
            spectrum = metadata['Spectrum']
            indices = np.argsort(spectrum[:, 1])[::-1][:10]
            spectrum = spectrum[indices]
            spectrum = spectrum[np.argsort(spectrum[:, 0])]
            metadata['Spectrum'] = spectrum

    # 转换为 CsoDIAq 库的格式
    for id, ((seq, charge), metadata) in enumerate(library.items()):
        is_decoy = 0
        if 'DECOY-' in seq:
            is_decoy = 1
        metadata.update({
            'PrecursorCharge': charge,
            'transition_group_id': 'DECOY' if is_decoy else 'TARGET',
            'ID': id,
            'ProteinName': '',
            'Decoy': is_decoy
        })

    csodiaq_library = transform_cosdiq_format(copy.deepcopy(library))

    return library, csodiaq_library


def fill_peaks(peaks: npt.NDArray, n_top: int = 6):
    length = len(peaks)
    if length < n_top:
        peaks = np.concatenate((np.zeros((n_top-length, 2)), peaks), axis=0)
    return peaks


def fill_fragments(fragments: npt.NDArray, n_top: int = 6):

    length = len(fragments)
    if length < n_top:
        fragments = np.concatenate((
            np.full((n_top - length, 4), None),
            fragments
        ), axis=0)
    return fragments


def fill_peaks_fragments(peaks: npt.NDArray, fragments: npt.NDArray, n_top: int = 6):
    peaks = fill_peaks(peaks, n_top)
    fragments = fill_fragments(fragments, n_top)
    return peaks, fragments


def selected_intensity_fragment_topN(peaks: npt.NDArray, fragments: npt.NDArray, n_top: int = 6):
    """
        选取峰强度前六的峰以及其离子类型

        可能出现有多个相同的 b/y 离子出现在筛选的峰中
    """
    selected_indexs = np.argsort(-peaks[:, 1])[:n_top]
    peaks, fragments = peaks[selected_indexs], fragments[selected_indexs]
    sort_indexs = np.argsort(peaks[:, 0])
    peaks = peaks[sort_indexs]
    fragments = fragments[sort_indexs]
    # peaks, fragments = fill_peaks_fragments(peaks, fragments, n_top)
    return peaks, fragments


def selected_intensity_topN(peaks: npt.NDArray, n_top: int = 6):
    selected_indexs = np.argsort(-peaks[:, 1])[:n_top]
    peaks = peaks[selected_indexs]
    sort_indexs = np.argsort(peaks[:, 0])
    peaks = peaks[sort_indexs]
    # peaks = fill_peaks(peaks, n_top)
    return peaks


def get_only_value(hashDict: dict, df: pd.DataFrame, columns: Sequence[str]):
    """
        选取第一个数值, 因为所有数都是相同的
    """
    for column in columns:
        hashDict[column] = df[column].values[0]


def process_fragment(hashDict: dict, df: pd.DataFrame, peak_columns: Sequence[str], fragment_columns: Sequence[str]):
    """
        将数据分为两部分

        1. 峰 [mz, intensity]

        2. 峰的类型以及原始形式 [type, number, lossType, charge] 
    """
    peaks = df[peak_columns].values
    fragments = df[fragment_columns].values
    # 如果无自定义的筛选函数, 则默认是选峰强度前六的峰
    # 存入字典中
    hashDict['Spectrum'] = peaks
    hashDict['Fragment'] = fragments


def target_library(raw_library_path: str):

    raw_library = pd.read_csv(raw_library_path, sep='\t', low_memory=False)
    group_columns = ['ModifiedPeptide', 'PrecursorCharge']
    #
    df = raw_library[['FragmentCharge', 'FragmentLossType', 'ExcludeFromAssay', 'ProteinGroups', 'ModifiedPeptide', 'StrippedPeptide',
                      'PrecursorCharge', 'iRT', 'IonMobility', 'PrecursorMz', 'FragmentNumber', 'FragmentType', 'FragmentMz', 'RelativeIntensity']]
    split_group = df.groupby(group_columns)
    # generate target library
    library = {}
    for key in split_group.groups.keys():
        modifiedPeptideProperty = split_group.get_group(key)
        library[key] = {}
        get_only_value(library[key], modifiedPeptideProperty, [
                       'StrippedPeptide', 'PrecursorMz', 'iRT', 'IonMobility', 'ProteinGroups'])
        process_fragment(library[key], modifiedPeptideProperty, ['FragmentMz', 'RelativeIntensity'], [
                         'FragmentType', 'FragmentNumber', 'FragmentLossType', 'FragmentCharge'])
    return library


def extend_library_tsv(
    tsv_path: str,
    entrap_distance: int,
    decoy_distance: int,
    is_entrap: bool = True,
    is_top10: bool = True
):
    smf.print_milestone(f"{'#'*15}  Loading library...  {'#'*15}")
    library = target_library(tsv_path)
    smf.print_milestone(f"{'#'*15}  End  {'#'*15}")
    if is_entrap:
        library = entrapment_decoy_strategy(
            library, entrap_distance, decoy_distance)
    else:
        library = decoy_strategy(library, decoy_distance)

    smf.print_milestone(f"{'#'*15}  End *** decoy library {'#'*15}")

    for id, ((seq, charge), metadata) in enumerate(library.items()):
        entrap = False
        decoy = 0
        if 'DECOY-' in seq:
            decoy = 1
        if 'TRAP-' in seq:
            entrap = True

        if is_top10 and not entrap and not decoy:
            peaks, fragment = selected_intensity_fragment_topN(
                metadata['Spectrum'], metadata['Fragment'], 10)
            metadata['Spectrum'] = peaks
            metadata['Fragment'] = fragment
        elif is_top10:
            peaks = selected_intensity_topN(metadata['Spectrum'], 10)
            metadata['Spectrum'] = peaks
            del metadata['Fragment']

        metadata.update({
            'PrecursorCharge': charge,
            'transition_group_id': 'DECOY' if decoy else 'TARGET',
            'ID': id,
            'ProteinName': metadata['ProteinGroups'],
            'Decoy': decoy,
            'Trap': entrap
        })

    csodiaq_library = transform_cosdiq_format(copy.deepcopy(library))

    return library, csodiaq_library
