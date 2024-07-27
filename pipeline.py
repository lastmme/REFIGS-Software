import os,sys
from typing import Sequence
from multiprocessing import Queue
import traceback

import numpy as np
from utils.library import extend_library
from CsoDIAq.csodiaq_identification_functions import perform_spectra_pooling_and_analysis

from REFIGS.identification import refigs_identification
from REFIGS.quantification import refigs_quantification


def pipeline(
    mzxml_files: Sequence[str],
    result_dir: str,
    library_file: str,
    is_entrap: bool,
    entrap_distance: int,
    decoy_distance: int,
    tolerance: int,
    is_corrected: bool,
    fdr_value: float,
    is_top10: bool,
    start_cycle: int,
    end_cycle: int,
    good_shared_limit: int,
    good_cos_sim_limit: int,
    good_sqrt_cos_sim_limit: int,
    good_count_within_cycle_limit: int,
    # scans_per_cycle: int,
    seed: int,
    text_queue: Queue,
    max_num_spectra_pool: int = np.inf 
):
    """
    TODO: 
    # 由于需要实现多进程，因此后面这个函数需要将 mzxml_files 这个参数单个的 mzxml_file
    -----------------------------------------------------------------------------
    # 因为采用多进程，你并不知道哪个进程先结束，因此 mzxml_file 应该是由一个主进程 (可以是 view_interface) 分发过来的数据
    # 我的想法是 
                             view_interface
                            /      |       \   ...      \
                         Worker1  Worker2  Worker3 ... Workern
    每个 Worker 之中包含了 pipeline 的这些参数 (可以通过 **kwargs 进行传输)

    每个 Worker 结束之后需要通过结束的信号触发 view_interface 中的

    Worker 之中需要实现的东西我也在相应的 TODO 中标注出来了

    大概就是参数 function (pipeline 函数), **kwargs 

    run 函数中
    -----
    1. 调用 QProcess 类，将函数传进去 
    2. QProcess.start()  进程启动
    3. QProcess.finished() ... 进程结束触发什么信号
    """
    try:
        
        # library, csodiaq_library = extend_library(
        #     library_file, entrap_distance,
        #     decoy_distance, text_queue, is_entrap, is_top10
        # )

        # np.save(
        #     os.path.join(result_dir, 'library.npy'),
        #     library
        # )

        # np.save(
        #     os.path.join(result_dir, 'csodiaq_library.npy'),
        #     csodiaq_library
        # )
        
        library=np.load(os.path.join(result_dir,'library.npy'),allow_pickle=True).item()
        csodiaq_library=np.load(os.path.join(result_dir,'csodiaq_library.npy'),allow_pickle=True).item()
        # text_queue.put('over success')
        # sys.exit()
        
        mzxml_dir, _ = os.path.split(mzxml_files[0])

        for mzxml_file in mzxml_files:
            _, file_name = os.path.split(mzxml_file)
            text_queue.put(
                f"{'#' * 15}  Analyzing the file {file_name.split('.')[0]}...  {'#' * 15}"
            )

            csv_file = file_name.replace('.mzXML', '.csv')

            csodiaq_savepath = perform_spectra_pooling_and_analysis(
                mzxml_file,
                os.path.join(result_dir, csv_file),
                csodiaq_library,
                text_queue,
                tolerance,
                max_num_spectra_pool,
                is_corrected
            )
            # csodiaq_savepath=os.path.join(result_dir,csv_file).replace('.csv', 'NoFilter.csv')

            refigs_ids_savepath = refigs_identification(
                fdr_value,
                text_queue,
                mzxml_file,
                csodiaq_savepath,
                library,
                start_cycle=start_cycle,
                end_cycle=end_cycle,
                good_shared_limit=good_shared_limit,
                good_cos_sim_limit=good_cos_sim_limit,
                good_sqrt_cos_sim_limit=good_sqrt_cos_sim_limit,
                good_count_within_cycle_limit=good_count_within_cycle_limit,
                tol=tolerance,
                # scans_per_cycle=scans_per_cycle,
                seed=seed
            )

            refigs_quantification(
                text_queue,
                refigs_ids_savepath,
                csodiaq_savepath,
                mzxml_file,
                library,
                tol=tolerance
            )
            text_queue.put(
                f"{'#' * 15}  End *** analyzing the file {file_name.split('.')[0]}...  {'#' * 15}"
            )
        text_queue.put('over success')        
        # text_queue.put('over')
    except Exception as e:
        # 将错误信息发送到 text_queue
        error_info = traceback.format_exc()
        text_queue.put(f'Process Error: {str(e)}')
        print(error_info)
        text_queue.put('over error')
        sys.exit()

def cal_library(
    result_dir: str,
    library_file: str,
    is_entrap: bool,
    entrap_distance: int,
    decoy_distance: int,
    is_top10: bool,
    text_queue:Queue,
    **kwargs
    ):
    
    # text_queue.put('cal library over')  
    # return
    
    library, csodiaq_library = extend_library(
        library_file, entrap_distance,
        decoy_distance, text_queue, is_entrap, is_top10
    )

    np.save(
        os.path.join(result_dir, 'library.npy'),
        library
    )

    np.save(
        os.path.join(result_dir, 'csodiaq_library.npy'),
        csodiaq_library
    )  
    text_queue.put('cal library over')  