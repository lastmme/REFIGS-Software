from multiprocessing import Queue

import pandas as pd
import numpy as np
from scipy.signal import find_peaks
from scipy.stats import moment
from statsmodels.stats.diagnostic import acorr_ljungbox

from .myms import spectrum_normalization, LoadMS2, deconvolute
from CsoDIAq import spectra_matcher_functions as smf


def refigs_quantification(
    text_queue: Queue,
    ids_file: str,
    SSM_file: str,
    mzxml_file: str,
    raw_library: dict[str, dict],
    tol: int = 20
):
    tol = tol * 1e-6

    library = spectrum_normalization(raw_library)
    text_queue.put(
        smf.print_milestone(
            f"{'#' * 15}  loading SSM file...  {'#' * 15}")
    )
    Identify = pd.read_csv(SSM_file, sep=',')
    Identify['Precursor'] = Identify[['peptide', 'zLIB']].apply(
        lambda x: x['peptide']+'_'+str(x['zLIB']), axis=1)
    smf.print_milestone(
        f"{'#' * 15}'  loading mzML file...  '{'#' * 15}")
    MS2 = LoadMS2(mzxml_file)
    MS2_dict = {}
    for ms2 in MS2:
        MS2_dict[ms2[-1]] = ms2

    text_queue.put(
        smf.print_milestone(f"{'#' * 15}  claculating coeffs...  {'#' * 15}")
    )
    output = []
    scans = set(Identify['scan'].values)
    for i, scan in enumerate(scans):
        if i % 1000 == 0:
            text_queue.put(
                smf.print_milestone(
                    f'{i}/{len(scans)}\t FIGS is calculating coeffs')
            )
        tmp_df = Identify[Identify['scan'] == scan]
        keys = tmp_df['Precursor'].values
        tmp_library = {}
        for key in keys:
            tmp_library[(key.split('_')[0], int(key.split('_')[1]))
                        ] = library[(key.split('_')[0], int(key.split('_')[1]))]
        tmp_ms2 = MS2_dict[scan]
        output.extend(deconvolute(tmp_ms2, tmp_library, tol, False))

    Coeffs = pd.DataFrame(output,columns=['Coeff', 'scan', 'peptide', 'charge', 'windowMZ', 'spectrumRT', 'level', 'corr'])
    # print(Coeffs)
    
    ids = pd.read_csv(ids_file, sep=',')
    Coeffs = Coeffs[Coeffs['peptide'].apply(
        lambda x: x in set(ids['peptide']))
    ]

    Coeffs.to_csv(SSM_file.replace('.csv', '_Coeffs.csv'), index=False, header=[
        'Coeff', 'scan', 'peptide', 'charge', 'windowMZ', 'spectrumRT', 'level', 'corr'])
    text_queue.put(
        smf.print_milestone('calculate coeffs over!')
    )
    
    text_queue.put(
        smf.print_milestone(f"{'#' * 15}  quantifying...  {'#' * 15}")
    )
    # 定量计算
    header = [[(x[1][0] + x[1][1]) / 2, x[2], x[3]] for x in MS2]
    header = pd.DataFrame(header,columns=['precursorMZ', 'RetentionTime', 'scan'])
    
    quantsWithDecoy=QuantifyAllFromCoeffs(Coeffs, header,MaxFillNum=3)
    quantsWithDecoy.to_csv(SSM_file.replace('.csv', '_Quants.csv'), index=False)
    text_queue.put(
        smf.print_milestone('quantification over!')
    )
    

def QuantifyAllFromCoeffs(
    coeffs:pd.DataFrame,
    header:pd.DataFrame,
    fill_flag:bool = True,
    MaxFillNum:int = 3
):
    IDs=set(coeffs[['peptide', 'charge','level']].apply(tuple, axis=1))
    results=[]
    for id in IDs:
        peptide_df=coeffs[(coeffs['peptide']==id[0]) & (coeffs['charge']==id[1]) & (coeffs['level']==id[2])]
        if EnoughData(peptide_df):
            MaxScan=peptide_df['scan'].max()
            MinScan=peptide_df['scan'].min()
            peptide_header=header[(header['scan'].between(MinScan, MaxScan)) & (header['precursorMZ'].isin(peptide_df['windowMZ']))]
            result=QuantifyPeptide(id,peptide_df, peptide_header, fill_flag=fill_flag, maxFillNum=MaxFillNum)
            results.append(result)
        else:
            result.append([id[0], id[1],id[2], 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]) 
            continue 
    #TODO: add columns
    columns=['peptide', 'charge','level', 'Quantify', 'LBPval', 'SNR', 'MaxCoeffTime', 'MaxCoeff', 'PeakWidth', 'LBPvalOnGrid', 'Variance', 'Skewness', 'Kurtosis', 'PeakStart', 'PeakEnd', 'Correlation']
    result_df=pd.DataFrame(results, columns=columns)
    return result_df
    

def EnoughData(
    peptide_coeffs: pd.DataFrame,
    RTwindow: bool = False
):
    # if RTwindow:
    #     peptide_coeffs = peptide_coeffs[peptide_coeffs['spectrumRT'] - peptide_coeffs['RefspectrumRT'] < 300]  
    if len(peptide_coeffs) > 0:
        return True
    else:
        return False

def QuantifyPeptide(
    id:tuple,
    peptide_coeffs:pd.DataFrame,
    peptide_header:pd.DataFrame,
    RTwindow:bool = False,
    fill_flag:bool = True,
    maxFillNum:int = 3
):
    # if RTwindow:
    #     peptide_coeffs=peptide_coeffs[peptide_coeffs['spectrumRT']-peptide_coeffs['RefspectrumRT']<300]
        
    peptide_coeffs=peptide_coeffs[~peptide_coeffs['spectrumRT'].duplicated()]
    peptide_coeffs['Coeff']=peptide_coeffs['Coeff'].apply(lambda x: 0 if x<1 else x)
    CoeffList=[0]*len(peptide_header)
    peptide_header=peptide_header.sort_values('scan',ascending=True)
    for i in range(len(peptide_header)):
        scan=peptide_header['scan'].iloc[i]
        if scan in set(peptide_coeffs['scan']):
            CoeffList[i]=peptide_coeffs[peptide_coeffs['scan']==scan]['Coeff'].values[0]   

    
    """
    TODO: 在补全的系数列表中截取第最高峰，计算面积
    最终返回的信息有：
        Peptide: 肽段序列。
        Charge: 电荷。
        level: 图谱等级。
        Quantity: 截取系数的梯形积分。
        LBPval: ljung-box 统计检验的 p-value ，未补全系数进行。
        SNR: 信噪比（Signal-to-Noise Ratio），表示信号相对于噪声的比率。
        MaxCoeffTime: 最大系数出现的保留时间。
        MaxCoeff: 最大系数值。
        PeakWidth: 峰宽，峰的起始到截止的保留时间时长。
        LBPvalOnGrid: 在均匀网格上的 ljung-box p-value，补全系数进行。
        PeakVariance: 峰的方差，表示峰的变化程度。
        PeakSkewness: 峰的偏度，表示峰的对称性。
        PeakKurtosis: 峰的峰度，表示峰的尖锐程度。
        PeakStart: 峰的起始位置，保留时间。
        PeakEnd: 峰的结束位置，保留时间。
        Correlation: 相关性。
    """ 
    BoxTestPval = 1
    BoxTestPvalOnGrid = 1
    MaxCoeff = 0
    TimeAtTopOfPeak = 0
    snr = 0
    PeakCentralMoments = [0] * 3
    area = 0

    columns = ['Peptide', 'Charge', 'Quantify', 'BoxTestPval', 'SNR', 'RT', 'MaxCoeff', 'PeakWidth', 'BoxTestPvalOnGrid', 'Variance', 'Skewness', 'Kurtosis', 'PeakStart', 'PeakEnd', 'Correlation']
    result = dict()
    # result['Peptide'].astype(str)
    result['Peptide'] = id[0]
    result['Charge'] = id[1]
    result['level'] = id[2]
    result['Quantify'] = area
    result['BoxTestPval'] = BoxTestPval
    result['SNR'] = snr
    result['RT'] = TimeAtTopOfPeak
    result['MaxCoeff'] = MaxCoeff
    result['PeakWidth'] = 0
    result['BoxTestPvalOnGrid'] = BoxTestPvalOnGrid
    result['Variance'] = 0
    result['Skewness'] = 0
    result['Kurtosis'] = 0
    result['PeakStart'] = 0
    result['PeakEnd'] = 0
    result['Correlation'] = 0    
    
    if len([i for i in range(len(CoeffList)) if CoeffList[i]>0])>3:
        if fill_flag:
            CoeffList=FillCoeffs(CoeffList,MaxFillNum=maxFillNum)
        SmoothPeak_max=FindPeaksInData(CoeffList)
        if SmoothPeak_max is not None:
            start=SmoothPeak_max[2]
            end=SmoothPeak_max[3]
            peak_idx=SmoothPeak_max[1]
            peakDataOnUniformGrid=peptide_header[start:end+1].copy()
            peakDataOnUniformGrid.loc[:,'Coeff']=CoeffList[start:end+1]
            
            area=np.trapz(peakDataOnUniformGrid['Coeff'], peakDataOnUniformGrid['RetentionTime'])
            TimeAtTopOfPeak=peptide_header['RetentionTime'].iloc[peak_idx]
            MaxCoeff=SmoothPeak_max[0]
            
            scaled_data = peakDataOnUniformGrid['Coeff'] / np.std(peakDataOnUniformGrid['Coeff'])
            PeakCentralMoments[0] = moment(scaled_data, moment=3)
            PeakCentralMoments[1] = moment(scaled_data, moment=4)
            PeakCentralMoments[2] = moment(scaled_data, moment=5)
            
            if peptide_coeffs['Coeff'].values.size>0:
                BoxTest=acorr_ljungbox(peptide_coeffs['Coeff'].values,return_df=True,auto_lag=True)
                BoxTestPval=round(BoxTest['lb_pvalue'].values[0], 5)
            if len(CoeffList)>0:
                BoxTestOnGrid=acorr_ljungbox(CoeffList,return_df=True,auto_lag=True)
                BoxTestPvalOnGrid=round(BoxTestOnGrid['lb_pvalue'].values[0], 5)    
            
            snr=area/np.std(np.concatenate((CoeffList[:start], CoeffList[end+1:])))   
            peakWidth=peptide_header['RetentionTime'].iloc[end]-peptide_header['RetentionTime'].iloc[start]
            Correlation=np.mean(peptide_coeffs['corr'])
            
            result['Quantify'] = area
            result['BoxTestPval'] = BoxTestPval
            result['SNR'] = snr
            result['RT'] = TimeAtTopOfPeak
            result['MaxCoeff'] = MaxCoeff
            result['PeakWidth'] = peakWidth
            result['BoxTestPvalOnGrid'] = BoxTestPvalOnGrid
            result['Variance'] = PeakCentralMoments[0]
            result['Skewness'] = PeakCentralMoments[1]
            result['Kurtosis'] = PeakCentralMoments[2]
            result['PeakStart'] = peptide_header['RetentionTime'].iloc[start]
            result['PeakEnd'] = peptide_header['RetentionTime'].iloc[end]
            result['Correlation'] = Correlation
    
    return [result['Peptide'],result['Charge'],result['level'],result['Quantify'],result['BoxTestPval'],result['SNR'],result['RT'],result['MaxCoeff'],result['PeakWidth'],result['BoxTestPvalOnGrid'],result['Variance'],result['Skewness'],result['Kurtosis'],result['PeakStart'],result['PeakEnd'],result['Correlation']]   


def FillCoeffs(
    coeffsList:list[float],
    MaxFillNum:int = 3
):
    if len(coeffsList)<=2:
        return coeffsList
    for i in range(1, len(coeffsList)):
        if coeffsList[i]==0 and coeffsList[i-1]>0:
            j=i+1
            while j<len(coeffsList) and coeffsList[j]==0 and j-i+1<=MaxFillNum:
                j+=1
            if j<len(coeffsList) and coeffsList[j]>0 and j-i<=MaxFillNum:
                grouth=(coeffsList[j]-coeffsList[i-1])/(j-i+1)
                for k in range(i, j-1):
                    coeffsList[k]=coeffsList[k-1]+grouth
    return coeffsList 

def FindPeaksInData(
    CoeffList:list[float],
    IntensityCutoff:int=0,
    QuantikeCutoff:int=0,
    smooth:str='rollmean',
    FilterWindow:int=3,
    KZiter:int=3,
    fill_flag:bool=True,
    MaxFillNum:int=3
):
    peptide_peaks=[]
    
    if smooth=='rollmean':
        KZData=kz_filter(CoeffList, FilterWindow, KZiter)
        KZData[np.where(KZData<1)]=0
        if fill_flag:
            KZData=FillCoeffs(KZData,MaxFillNum=MaxFillNum)
            
        peaks, _ = find_peaks(KZData)
        for peak in peaks:
            up=peak+1
            down=peak-1
            while up<len(KZData) and KZData[up]<=KZData[up-1] and KZData[up]>1:
                up+=1
            while down>=0 and KZData[down]<=KZData[down+1] and KZData[down]>1:
                down-=1
            if up-peak-1>=2 and peak-down-1>=2:
                peptide_peaks.append([KZData[peak], peak, down+1, up-1])
    if len(peptide_peaks)==0:
        return None
    else:
        peak_max = np.argmax(np.array(peptide_peaks)[:, 0])
        return peptide_peaks[peak_max]
    



def kz_filter(data, window_size, iterations):
    for _ in range(iterations):
        data = np.convolve(data, np.ones(window_size) / window_size, mode='valid')
    return data

if __name__ == '__main__':
    ids_file='../../../data/refigs_test/temp/480_20210929_noRT_FAIMS_30to80_120ms_40min_250nL_A_1NoFilter_withFeature_19cycle_4_3_LDA_ID.csv'
    ssm_file='../../../data/refigs_test/temp/480_20210929_noRT_FAIMS_30to80_120ms_40min_250nL_A_1NoFilter.csv'
    text_queue = Queue()
    mzxml_file = '../../../data/refigs_test/mzxml/480_20210929_noRT_FAIMS_30to80_120ms_40min_250nL_A_1.mzXML'


    MS2=LoadMS2(mzxml_file)
    header = [[(x[1][0] + x[1][1]) / 2, x[2], x[3]] for x in MS2]
    header = pd.DataFrame(header,columns=['precursorMZ', 'RetentionTime', 'scan'])
        
    Coeffs=pd.read_csv('../../../data/refigs_test/temp/480_20210929_noRT_FAIMS_30to80_120ms_40min_250nL_A_1NoFilter_Coeffs.csv')

    quants=QuantifyAllFromCoeffs(Coeffs,header)
    quants.to_csv('../../../data/refigs_test/temp/480_20210929_noRT_FAIMS_30to80_120ms_40min_250nL_A_1NoFilter_Quants.csv',index=False)