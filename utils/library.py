import copy
from multiprocessing import Queue

import numpy as np

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
