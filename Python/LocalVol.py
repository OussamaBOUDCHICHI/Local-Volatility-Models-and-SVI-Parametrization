# -*- coding : utf8 -*-
# author : BOUDCHICHI Oussama
# Local volatility calibration using finite differences


from scipy.interpolate import CubicSpline
import numpy as np
import pandas as pd


def cubic_Spline(data, N):
    '''
    Cubic Splines interpolation
    
    params:
    =======
    data : pandas.DataFrame
           data containing strikes in *rows index* and maturities in *columns index*
    N    : int
           Grid size.
    
    returns:
    ========
    final : pandas.DataFrame
            interpolated values (it preseverves the same structure as the initial dataframe)  
    '''
    
    # Strike Interpolation :

    kx = data.index
    Ks = np.linspace(kx[0], kx[-1], N)
    inter = {}
    for j in range(data.shape[1]):

        yx = data.iloc[:,j]
        cs = CubicSpline(kx, yx)
        inter[data.columns[j]] = pd.Series(cs(Ks), index=Ks)

    temp = pd.concat(inter, axis=1)

    # Maturity interpolation

    tx = data.columns
    Ts = np.linspace(tx[0], tx[-1], N)
    inter = {}
    for j in range(temp.shape[0]):

        yx = temp.iloc[j,:]
        cs = CubicSpline(tx, yx)
        inter[temp.index[j]] = pd.Series(cs(Ts), index=Ts)

    final = pd.concat(inter,axis=1).T 
    return final  




def LocalVol(data, r, method = 'Price'):
    '''
    Local Volatility calibration function
    
    params:
    =======
    data      : pandas.DataFrame
                data containing strikes in *rows index* and maturities in *columns index*
    method    : string
                string indicating which formula to use (Price / ImpliedVol) .
    
    returns:
    ========
    Sig : pandas.DataFrame
          **Squared** (note this) local volatilies (it preseverves the same structure as the initial dataframe)  
    '''
    
    C = data.copy()
    Sig_Dup ={}

    if method == 'Price':
        for j in range(C.shape[1]-1):
            dC_dT = np.array((C.iloc[:,j+1] - C.iloc[:,j]) / (C.columns[j+1] - C.columns[j]))
            
            dC_dK = [(C.iloc[0+1,j] - C.iloc[0,j]) / (C.index[1] - C.index[0])]
            dC_dKsq = []

            for i in range(1,C.shape[0]-1):
                dC_dK.append((C.iloc[i+1,j] - C.iloc[i,j]) / (C.index[i+1]-C.index[i]))
                dC_dKsq.append((C.iloc[i+1,j] - 2 * C.iloc[i,j] + C.iloc[i-1,j])  /((C.index[i+1]-C.index[i])**2))

            
            Denom = (C.index[1:-1]**2)* np.array(dC_dKsq)

            Sig_Dup[str(np.round(C.columns[j],2))] = pd.Series(2 * (dC_dT[1:-1] + r * C.index[1:-1] * np.array(dC_dK[:-1])) / Denom , index=C.index[1:-1])
    
    if method == 'ImpliedVol':
        for j in range(C.shape[1]-1):
            dsig_dT = np.array((C.iloc[:,j+1] - C.iloc[:,j]) / (C.columns[j+1] - C.columns[j]))
            
            dsig_dK = [(C.iloc[0+1,j] - C.iloc[0,j]) / (C.index[1] - C.index[0])]
            dsig_dKsq = []

            for i in range(1,C.shape[0]-1):
                dsig_dK.append((C.iloc[i+1,j] - C.iloc[i,j]) / (C.index[i+1]-C.index[i]))
                dsig_dKsq.append((C.iloc[i+1,j] - 2 * C.iloc[i,j] + C.iloc[i-1,j])  /((C.index[i+1]-C.index[i])**2))

            
            Num  = C.iloc[1:-1,j]**2 + 2 * C.iloc[1:-1,j] * C.columns[j] * (dsig_dT[1:-1] + r * C.index[1:-1] * np.array(dsig_dK[:-1]))
            x = np.log(100.0 / C.index[1:-1]) - r * C.columns[j]
            Denom = (1 - (C.index[1:-1] * x)  / C.iloc[1:-1,j])**2 + C.index[1:-1] * C.iloc[1:-1,j] * C.columns[j] * (dsig_dK[:-1] - 0.25 * C.index[1:-1] * C.columns[j]* C.iloc[1:-1,j] * (np.array(dsig_dK[:-1]))**2)
            Denom += C.index[1:-1] * dsig_dKsq 

            Sig_Dup[str(np.round(C.columns[j],2))] = pd.Series( Num / Denom , index=C.index[1:-1])

    Sig = pd.concat(Sig_Dup, axis=1)

    return Sig