import sys
from multiprocessing import Queue

import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

from .myms import LoadMS2, get_scanNum_per_cycle, Peaks_match, spectrum_normalization
from CsoDIAq import spectra_matcher_functions as smf


def refigs_identification(
    fdr_value: float,
    text_queue: Queue,
    mzxml_file: str,
    ssm_file: str,
    library_raw: dict[str, dict],
    start_cycle: int = 1,
    end_cycle: int = 14,
    good_shared_limit: int = 5,
    good_cos_sim_limit: float = 0.8,
    good_sqrt_cos_sim_limit: float = 0.9,
    good_count_within_cycle_limit: int = 5,
    tol: int = 20,
    seed: int = 13,
):
    # print(fdr_value, mzxml_file, ssm_file, start_cycle, end_cycle, good_shared_limit, good_cos_sim_limit, good_sqrt_cos_sim_limit, good_count_within_cycle_limit, tol, seed)
    tol = tol * 1e-6
    cycle_range = range(start_cycle, end_cycle)
    cycle_number = len(cycle_range)

    # libray
    library = spectrum_normalization(library_raw)
    text_queue.put(
        smf.print_milestone('loading SSM file')
    )

    scans_per_cycle = get_scanNum_per_cycle(mzxml_file)

    # loading ssm result
    df = pd.read_csv(ssm_file)
    if not 'cycle' in df.columns:
        df['cycle'] = df['scan'].apply(lambda x: (x-1)//scans_per_cycle+1)

    # cycle filtering
    df = df[df['cycle'].isin(cycle_range)]
    df = df.reset_index()
    df['Precursor'] = df[['peptide', 'zLIB']].apply(
        lambda x: x['peptide']+'_'+str(x['zLIB']), axis=1)

    text_queue.put(
        smf.print_milestone("loading MS2!")
    )
    MS2 = LoadMS2(mzxml_file)

    MS2_dict = {}
    for ms2 in (MS2):
        MS2_dict[ms2[-1]] = ms2[:-1]

    # Feature extraction(calculation)
    match_peaks = pd.DataFrame([])
    scan_list = df['scan'].values
    peptide_list = df['peptide'].values
    zLIB_list = df['zLIB'].values
    cosine_list = df['cosine'].values
    name_list = df['name'].values
    Peaks_Library_list = df['Peaks(Library)'].values
    shared_list = df['shared'].values
    MaCC_list = df['MaCC_Score'].values
    cycle_list = df['cycle'].values

    cos_sim_list = np.zeros(len(df))
    sqrt_cos_sim_list = np.zeros(len(df))
    norm_mse_list = np.zeros(len(df))
    norm_mae_list = np.zeros(len(df))
    match_it_ratio_list = np.zeros(len(df))
    match_number_list = np.zeros(len(df))

    text_queue.put(
        smf.print_milestone("Calculating features!")
    )
    col_map = {
        col: i
        for i, col in enumerate(df.columns)
    }

    for idx, row in enumerate(df.values):
        if idx % 1000 == 0:
            text_queue.put(
                smf.print_milestone(
                    f'{idx}/{df.shape[0]}]\tcalculating features')
            )
        scan = row[col_map['scan']]
        precursor = row[col_map['Precursor']]
        key = (row[col_map['peptide']], row[col_map['zLIB']])
        lib_spectrum = library[key]
        exp_spectrum = MS2_dict[scan]
        lib_spc = lib_spectrum['Spectrum']
        exp_spc = exp_spectrum[0]

        match_mz, match_lib_it, match_exp_it = Peaks_match(
            lib_spc, exp_spc, tol)
        # if len(match_mz) != 0:
        # smf.print_milestone(f'{key} No matching peaks!')
        match_lib_it = np.array(match_lib_it)
        match_exp_it = np.array(match_exp_it)

        cos_sim = cosine_similarity(match_lib_it.reshape(
            1, len(match_lib_it)), match_exp_it.reshape(1, len(match_exp_it)))[0][0]

        sqrt_match_lib_it = np.sqrt(match_lib_it)
        sqrt_match_exp_it = np.sqrt(match_exp_it)
        sqrt_cos_sim = cosine_similarity(sqrt_match_lib_it.reshape(1, len(
            sqrt_match_lib_it)), sqrt_match_exp_it.reshape(1, len(sqrt_match_exp_it)))[0][0]

        l2_norm_sqrt_match_lib_it = match_lib_it / \
            np.linalg.norm(match_lib_it + 1e-6)
        l2_norm_sqrt_match_exp_it = match_exp_it / \
            np.linalg.norm(match_exp_it + 1e-6)
        norm_mse = mean_squared_error(
            l2_norm_sqrt_match_lib_it, l2_norm_sqrt_match_exp_it)
        norm_mae = np.mean(
            abs(l2_norm_sqrt_match_lib_it-l2_norm_sqrt_match_exp_it))

        match_it_ratio = np.sum(match_lib_it)/np.sum(lib_spc[:, 1])
        # else:
        #     cos_sim = 0
        #     sqrt_cos_sim = 0
        #     norm_mse = 0
        #     norm_mae = 0
        #     match_it_ratio = 0

        cos_sim_list[idx] = cos_sim
        sqrt_cos_sim_list[idx] = sqrt_cos_sim
        norm_mse_list[idx] = norm_mse
        norm_mae_list[idx] = norm_mae
        match_it_ratio_list[idx] = match_it_ratio
        match_number_list[idx] = int(len(match_lib_it))

    match_peaks['scan'] = scan_list
    match_peaks['peptide'] = peptide_list
    match_peaks['zLIB'] = zLIB_list
    match_peaks['cosine'] = cosine_list
    match_peaks['name'] = name_list
    match_peaks['Peaks(Library)'] = Peaks_Library_list
    match_peaks['shared'] = shared_list
    match_peaks['MaCC_Score'] = MaCC_list
    match_peaks['cycle'] = cycle_list
    match_peaks['cos_sim'] = cos_sim_list
    match_peaks['sqrt_cos_sim'] = sqrt_cos_sim_list
    match_peaks['norm_mse'] = norm_mse_list
    match_peaks['norm_mae'] = norm_mae_list
    match_peaks['match_it_ratio'] = match_it_ratio_list
    match_peaks['match_number'] = match_number_list
    feature_filepath = ssm_file.replace(
        '.csv', '_withFeature_'+str(cycle_number)+'cycle.csv')
    match_peaks.to_csv(feature_filepath, index=False)

    # RE-FIGS distinguish decoy by the peptide sequence.
    # Decoy's name starts with "DECOY"

    df = pd.read_csv(feature_filepath)
    df['protein'] = df['name'].apply(
        lambda x: 'DECOY_null' if x.startswith('DECOY') else 'TARGET')
    df['label'] = df['protein'].apply(lambda x: 0 if x == 'DECOY_null' else 1)

    cycle_max = max(df['cycle'].values)

    # calculate the number of peptide within cycle. and keep peptide with best MaCC_Score
    df_dup_rm = pd.DataFrame([])
    for i in (range(1, cycle_max+1)):
        df_cycle = df[df['cycle'] == i]
        count = df_cycle['peptide'].value_counts()
        count_df = count.reset_index()
        count_df['peptide'] = count.index
        count_df['count_within_cycle'] = count.values
        count_df = pd.DataFrame(
            count_df[['peptide', 'count_within_cycle']].values)
        count_df.columns = ['peptide', 'count_within_cycle']
        df_cycle = df_cycle.sort_values('MaCC_Score', ascending=False).groupby(
            'peptide', as_index=False).first()
        df_cycle = pd.merge(df_cycle, count_df)
        df_dup_rm = pd.concat([df_dup_rm, df_cycle], sort=False)

    df_run = df_dup_rm.copy()
    # calculate the number of peptide between cycle. and keep peptide with best MaCC_Score
    count = df_run['peptide'].value_counts()
    count_df = count.reset_index()
    count_df['peptide'] = count.index
    count_df['cycle_count'] = count.values
    count_df = pd.DataFrame(count_df[['peptide', 'cycle_count']].values)
    count_df.columns = ['peptide', 'cycle_count']
    df_run = df_run.sort_values('MaCC_Score', ascending=False).groupby(
        'peptide', as_index=False).first()
    df_run = pd.merge(df_run, count_df)

    # seed for reproduction
    np.random.seed(seed)
    final_df = df_run.copy()
    target = final_df[final_df['protein'] == 'TARGET']
    decoy = final_df[final_df['protein'] != 'TARGET']

    # choose good target peptide with higher requirements.
    # the choosen good target can affect the classification.
    good_target = target[target['shared'] > good_shared_limit]
    good_target = good_target[good_target['cos_sim'] > good_cos_sim_limit]
    good_target = good_target[good_target['sqrt_cos_sim']
                              > good_sqrt_cos_sim_limit]
    good_target = good_target[good_target['count_within_cycle']
                              >= good_count_within_cycle_limit]
    if len(good_target) < 2:
        text_queue.put(
            smf.print_milestone("The program has ended ahead of schedule, the number of good target is "+"%d" %
                                len(good_target)+", less than 2! Please choose another parameter.")
        )
        text_queue.put('Over Fail')
        sys.exit()
    if len(good_target) < 500:
        text_queue.put(
            smf.print_milestone("Not a warning, the number of good target is "+"%d" %
                                len(good_target)+", less than 500!")
        )
    size = len(good_target)
    weights_list = []
    text_queue.put(
        smf.print_milestone("Ensemble learning LDA!")
    )
    # Ensemble learning。Train 10 LDA model and average the weights.
    for j in range(10):
        randidx = np.random.randint(0, len(decoy), size)
        choosen_decoy = decoy.iloc[randidx]
        DataSet_df = pd.concat([good_target, choosen_decoy])
        col = ['shared', 'MaCC_Score', 'cos_sim', 'sqrt_cos_sim', 'norm_mse',
               'norm_mae', 'match_it_ratio', 'count_within_cycle', 'cycle_count', 'label']
        DataSet = np.array(DataSet_df[col])
        X = DataSet[:, :-1]
        y = DataSet[:, -1].astype(int)
        clf = LinearDiscriminantAnalysis()
        clf.fit(X, y)
        weights = clf.coef_[0]
        weights_list.append(weights)

    weights_list = np.array(weights_list)
    avg_weights = np.mean(weights_list, axis=0)
    # calculate the final LDA score by weights the other features.
    total_X = np.array(final_df[col])[:, :-1]
    scores = np.dot(total_X, avg_weights)
    final_df['LDA_Score'] = scores

    # FDR calculate
    final_df = final_df.sort_values(by="LDA_Score", ascending=False)
    FDR = []
    decoy_number = 0
    target_number = 0
    for idx, row in final_df.iterrows():
        protein = row['protein']
        if protein == 'DECOY_null':
            decoy_number += 1
        else:
            target_number += 1
        # for decoy/target in library is 1:1
        if target_number == 0:
            FDR.append(1)
        else:
            FDR.append(decoy_number/target_number)
    final_df['FDR'] = FDR

    # p-value calculate
    final_df = final_df.sort_values(by="LDA_Score", ascending=False)
    pvalue = []
    decoy_number = 0
    total_decoy = len(decoy)
    for idx, row in final_df.iterrows():
        protein = row['protein']
        if protein == 'DECOY_null':
            decoy_number += 1
        else:
            target_number += 1
        # for decoy/target in library is 1:1
        pvalue.append((decoy_number+1)/(total_decoy+1))
    final_df['p-value'] = pvalue
    final_df.to_csv(feature_filepath.replace(
        '.csv', '_pvalue.csv'), index=False)

    # Q-value filter.
    # find the least LDA_score which FDR<= fdf_value
    cutoff = 1e9
    reverse_df_dup_rm = final_df.iloc[::-1]
    for idx, row in reverse_df_dup_rm.iterrows():
        if row['FDR'] <= fdr_value:
            cutoff = row['LDA_Score']
            break
    good = final_df[final_df['LDA_Score'] >= cutoff]
    good.to_csv(feature_filepath.replace('.csv', '_'+str(good_shared_limit) +
                '_'+str(good_count_within_cycle_limit)+'_LDA_ID.csv'), index=False)
    text_queue.put(
        smf.print_milestone('finished')
    )
    return feature_filepath.replace('.csv', '_'+str(good_shared_limit) +
                                    '_'+str(good_count_within_cycle_limit)+'_LDA_ID.csv')


def refigs_identification_for_args(
    fdr_value:float,
    text_queue:Queue,
    mzxml_file:str,
    ssm_file:str,
    library_raw: dict[str,dict],
    start_cycle:int=1,
    end_cycle:int=14,
    tol:int=20,
    seed:int=13,
):
    tol=tol * 1e-6
    
    library=spectrum_normalization(library_raw)
    scans_per_cycle=get_scanNum_per_cycle(mzxml_file)
    df=pd.read_csv(ssm_file)
    
    if not 'cycle' in df.columns:
        df['cycle']=df['scan'].apply(lambda x:(x-1)//scans_per_cycle+1)
    
    
    
    # TODO: cycle range
    cycle_range=range(start_cycle,end_cycle)
    cycle_number=len(cycle_range)
    
     
    max_cycle_count=cycle_number
    
    df=df[df['cycle'].isin(cycle_range)]
    df=df.reset_index()
    df['Precursor'] = df[['peptide', 'zLIB']].apply(
        lambda x: x['peptide']+'_'+str(x['zLIB']), axis=1)
    
    MS2=LoadMS2(mzxml_file)
    
    MS2_dict = {}
    for ms2 in (MS2):
        MS2_dict[ms2[-1]] = ms2[:-1]

    # Feature extraction(calculation)
    match_peaks = pd.DataFrame([])
    scan_list = df['scan'].values
    peptide_list = df['peptide'].values
    zLIB_list = df['zLIB'].values
    cosine_list = df['cosine'].values
    name_list = df['name'].values
    Peaks_Library_list = df['Peaks(Library)'].values
    shared_list = df['shared'].values
    MaCC_list = df['MaCC_Score'].values
    cycle_list = df['cycle'].values

    cos_sim_list = np.zeros(len(df))
    sqrt_cos_sim_list = np.zeros(len(df))
    norm_mse_list = np.zeros(len(df))
    norm_mae_list = np.zeros(len(df))
    match_it_ratio_list = np.zeros(len(df))
    match_number_list = np.zeros(len(df))

    
    col_map = {
        col: i
        for i, col in enumerate(df.columns)
    }

    for idx, row in enumerate(df.values):
        
        scan = row[col_map['scan']]
        precursor = row[col_map['Precursor']]
        key = (row[col_map['peptide']], row[col_map['zLIB']])
        lib_spectrum = library[key]
        exp_spectrum = MS2_dict[scan]
        lib_spc = lib_spectrum['Spectrum']
        exp_spc = exp_spectrum[0]

        match_mz, match_lib_it, match_exp_it = Peaks_match(
            lib_spc, exp_spc, tol)
        # if len(match_mz) != 0:
        # smf.print_milestone(f'{key} No matching peaks!')
        match_lib_it = np.array(match_lib_it)
        match_exp_it = np.array(match_exp_it)

        cos_sim = cosine_similarity(match_lib_it.reshape(
            1, len(match_lib_it)), match_exp_it.reshape(1, len(match_exp_it)))[0][0]

        sqrt_match_lib_it = np.sqrt(match_lib_it)
        sqrt_match_exp_it = np.sqrt(match_exp_it)
        sqrt_cos_sim = cosine_similarity(sqrt_match_lib_it.reshape(1, len(
            sqrt_match_lib_it)), sqrt_match_exp_it.reshape(1, len(sqrt_match_exp_it)))[0][0]

        l2_norm_sqrt_match_lib_it = match_lib_it / \
            np.linalg.norm(match_lib_it + 1e-6)
        l2_norm_sqrt_match_exp_it = match_exp_it / \
            np.linalg.norm(match_exp_it + 1e-6)
        norm_mse = mean_squared_error(
            l2_norm_sqrt_match_lib_it, l2_norm_sqrt_match_exp_it)
        norm_mae = np.mean(
            abs(l2_norm_sqrt_match_lib_it-l2_norm_sqrt_match_exp_it))

        match_it_ratio = np.sum(match_lib_it)/np.sum(lib_spc[:, 1])

        cos_sim_list[idx] = cos_sim
        sqrt_cos_sim_list[idx] = sqrt_cos_sim
        norm_mse_list[idx] = norm_mse
        norm_mae_list[idx] = norm_mae
        match_it_ratio_list[idx] = match_it_ratio
        match_number_list[idx] = int(len(match_lib_it))

    match_peaks['scan'] = scan_list
    match_peaks['peptide'] = peptide_list
    match_peaks['zLIB'] = zLIB_list
    match_peaks['cosine'] = cosine_list
    match_peaks['name'] = name_list
    match_peaks['Peaks(Library)'] = Peaks_Library_list
    match_peaks['shared'] = shared_list
    match_peaks['MaCC_Score'] = MaCC_list
    match_peaks['cycle'] = cycle_list
    match_peaks['cos_sim'] = cos_sim_list
    match_peaks['sqrt_cos_sim'] = sqrt_cos_sim_list
    match_peaks['norm_mse'] = norm_mse_list
    match_peaks['norm_mae'] = norm_mae_list
    match_peaks['match_it_ratio'] = match_it_ratio_list
    match_peaks['match_number'] = match_number_list
    feature_filepath = ssm_file.replace(
        '.csv', '_withFeature_'+str(cycle_number)+'cycle.csv')
    match_peaks.to_csv(feature_filepath, index=False)

    # RE-FIGS distinguish decoy by the peptide sequence.
    # Decoy's name starts with "DECOY"

    df = pd.read_csv(feature_filepath)
    df['protein'] = df['name'].apply(
        lambda x: 'DECOY_null' if x.startswith('DECOY') else 'TARGET')
    df['label'] = df['protein'].apply(lambda x: 0 if x == 'DECOY_null' else 1)

    cycle_max = max(df['cycle'].values)

    # calculate the number of peptide within cycle. and keep peptide with best MaCC_Score
    df_dup_rm = pd.DataFrame([])
    for i in (range(1, cycle_max+1)):
        df_cycle = df[df['cycle'] == i]
        count = df_cycle['peptide'].value_counts()
        count_df = count.reset_index()
        count_df['peptide'] = count.index
        count_df['count_within_cycle'] = count.values
        count_df = pd.DataFrame(
            count_df[['peptide', 'count_within_cycle']].values)
        count_df.columns = ['peptide', 'count_within_cycle']
        df_cycle = df_cycle.sort_values('MaCC_Score', ascending=False).groupby(
            'peptide', as_index=False).first()
        df_cycle = pd.merge(df_cycle, count_df)
        df_dup_rm = pd.concat([df_dup_rm, df_cycle], sort=False)

    df_run = df_dup_rm.copy()
    # calculate the number of peptide between cycle. and keep peptide with best MaCC_Score
    count = df_run['peptide'].value_counts()
    count_df = count.reset_index()
    count_df['peptide'] = count.index
    count_df['cycle_count'] = count.values
    count_df = pd.DataFrame(count_df[['peptide', 'cycle_count']].values)
    count_df.columns = ['peptide', 'cycle_count']
    df_run = df_run.sort_values('MaCC_Score', ascending=False).groupby(
        'peptide', as_index=False).first()
    df_run = pd.merge(df_run, count_df)

    # seed for reproduction
    np.random.seed(seed)
    final_df = df_run.copy()
    target = final_df[final_df['protein'] == 'TARGET']
    decoy = final_df[final_df['protein'] != 'TARGET']
    
    # TODO: args range
    max_lib_and_exp_shareds=6
    good_shared_limits=range(0,max_lib_and_exp_shareds,1)
    good_cos_sim_limits=np.arange(0.5,1,0.1)
    good_sqrt_cos_sim_limits=np.arange(0.5,1,0.1)
    good_count_within_cycle_limits=range(0,max_cycle_count,1)
    
    good_shared_limit=5
    good_cos_sim_limit=0.7
    good_sqrt_cos_sim_limit=0.7
    good_count_within_cycle_limit=3
    
    max_identification_count=0
    
    
    args_columns=['max_identification_count','good_shared_limit','good_cos_sim_limit','good_sqrt_cos_sim_limit','good_count_within_cycle_limit']
    args_list=[]
    
    for arg_shared_limit in good_shared_limits:
        good_target = target[target['shared'] > arg_shared_limit]
        good_target = good_target[good_target['cos_sim'] > good_cos_sim_limit]
        good_target = good_target[good_target['sqrt_cos_sim']
                                > good_sqrt_cos_sim_limit]
        good_target = good_target[good_target['count_within_cycle']
                                >= good_count_within_cycle_limit]
        size = len(good_target)
        if size < 2:
            continue
        identification_count=LDA_learn_and_FDR_filter(final_df,decoy,good_target,size,fdr_value)
        if identification_count>max_identification_count:
            max_identification_count=identification_count
            good_shared_limit=arg_shared_limit
            args_list.append([max_identification_count,good_shared_limit,good_cos_sim_limit,good_sqrt_cos_sim_limit,good_count_within_cycle_limit])
            
    for arg_cos_sim_limit in good_cos_sim_limits:
        good_target = target[target['shared'] > good_shared_limit]
        good_target = good_target[good_target['cos_sim'] > arg_cos_sim_limit]
        good_target = good_target[good_target['sqrt_cos_sim']
                                > good_sqrt_cos_sim_limit]
        good_target = good_target[good_target['count_within_cycle'] >= good_count_within_cycle_limit]
        size = len(good_target)
        if size < 2:
            continue    
        identification_count=LDA_learn_and_FDR_filter(final_df,decoy,good_target,size,fdr_value)    
        if identification_count>max_identification_count:
            max_identification_count=identification_count
            good_cos_sim_limit=arg_cos_sim_limit
            args_list.append([max_identification_count,good_shared_limit,good_cos_sim_limit,good_sqrt_cos_sim_limit,good_count_within_cycle_limit])
    for arg_sqrt_cos_sim_limit in good_sqrt_cos_sim_limits:
        good_target = target[target['shared'] > good_shared_limit]
        good_target = good_target[good_target['cos_sim'] > good_cos_sim_limit]
        good_target = good_target[good_target['sqrt_cos_sim']
                                > arg_sqrt_cos_sim_limit]
        good_target = good_target[good_target['count_within_cycle'] >= good_count_within_cycle_limit]
        size = len(good_target)
        if size < 2:
            continue    
        identification_count=LDA_learn_and_FDR_filter(final_df,decoy,good_target,size,fdr_value)    
        if identification_count>max_identification_count:
            max_identification_count=identification_count
            good_sqrt_cos_sim_limit=arg_sqrt_cos_sim_limit  
            args_list.append([max_identification_count,good_shared_limit,good_cos_sim_limit,good_sqrt_cos_sim_limit,good_count_within_cycle_limit])
    for arg_count_within_cycle_limit in good_count_within_cycle_limits:
        good_target = target[target['shared'] > good_shared_limit]
        good_target = good_target[good_target['cos_sim'] > good_cos_sim_limit]
        good_target = good_target[good_target['sqrt_cos_sim']
                                > good_sqrt_cos_sim_limit]
        good_target = good_target[good_target['count_within_cycle'] >= arg_count_within_cycle_limit]
        size = len(good_target)
        if size < 2:
            continue    
        identification_count=LDA_learn_and_FDR_filter(final_df,decoy,good_target,size,fdr_value)    
        if identification_count>max_identification_count:
            max_identification_count=identification_count
            good_count_within_cycle_limit=arg_count_within_cycle_limit
            args_list.append([max_identification_count,good_shared_limit,good_cos_sim_limit,good_sqrt_cos_sim_limit,good_count_within_cycle_limit])
    
    args_df=pd.DataFrame(args_list,columns=args_columns)
    print(args_df)
    return [max_identification_count,good_shared_limit,good_cos_sim_limit,good_sqrt_cos_sim_limit,good_count_within_cycle_limit]


def LDA_learn_and_FDR_filter(
    final_df:pd.DataFrame,
    decoy:pd.DataFrame,
    good_target:pd.DataFrame,
    size:int,
    fdr_value:float
    
):
    weights_list = []
        
    for j in range(10):
        randidx = np.random.randint(0, len(decoy), size)
        choosen_decoy = decoy.iloc[randidx]
        DataSet_df = pd.concat([good_target, choosen_decoy])
        col = ['shared', 'MaCC_Score', 'cos_sim', 'sqrt_cos_sim', 'norm_mse',
            'norm_mae', 'match_it_ratio', 'count_within_cycle', 'cycle_count', 'label']
        DataSet = np.array(DataSet_df[col])
        X = DataSet[:, :-1]
        y = DataSet[:, -1].astype(int)
        clf = LinearDiscriminantAnalysis()
        clf.fit(X, y)
        weights = clf.coef_[0]
        weights_list.append(weights)

    weights_list = np.array(weights_list)
    avg_weights = np.mean(weights_list, axis=0)
    # calculate the final LDA score by weights the other features.
    total_X = np.array(final_df[col])[:, :-1]
    scores = np.dot(total_X, avg_weights)
    
    test_df=final_df.copy()
    test_df['LDA_Score'] = scores

    # FDR calculate
    test_df = test_df.sort_values(by="LDA_Score", ascending=False)
    FDR = []
    decoy_number = 0
    target_number = 0
    for idx, row in test_df.iterrows():
        protein = row['protein']
        if protein == 'DECOY_null':
            decoy_number += 1
        else:
            target_number += 1
        # for decoy/target in library is 1:1
        if target_number == 0:
            FDR.append(1)
        else:
            FDR.append(decoy_number/target_number)
    test_df['FDR'] = FDR

    # p-value calculate
    test_df = test_df.sort_values(by="LDA_Score", ascending=False)
    decoy_number = 0
    
    cutoff = 1e9
    reverse_df_dup_rm = test_df.iloc[::-1]
    for idx, row in reverse_df_dup_rm.iterrows():
        if row['FDR'] <= fdr_value:
            cutoff = row['LDA_Score']
            break
    good = test_df[test_df['LDA_Score'] >= cutoff]
    
    identification_count=len(good[good['name']=='TARGET'])
    
    return identification_count