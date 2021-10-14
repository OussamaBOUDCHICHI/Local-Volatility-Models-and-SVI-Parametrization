# -*- coding : utf8 -*-
# author : BOUDCHICHI Oussama
# Raw SVI parametrization and calibration




import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.optimize import Bounds, minimize
import matplotlib.pyplot as plt

def x(K,T, r, S_0):
    
    '''
    log-moneyness
    
    params:
    =======
    K : float
        Strike
    T : float
        Maturity
    r : float
        Risk-free rate
    S_0 : float
        Initial value of price
    
    returns:
    ========
    x : float
        log-moneyness  
    '''
    return np.log(K / S_0) - r * T


def SVI(theta, x):
    '''
    SVI parametrization
    
    params:
    =======
    theta : tuple
            vector of parameters
    x     : float 
            log-moneyness
    
    returns:
    ========
    w : float
        total implied variance  
    '''
    a, b, rho, m, sigma = theta
    return a + b * (rho * (x-m) + np.sqrt((x-m)**2 + sigma**2))

def g(theta, x):
    '''
    g function for convexity test
    
    params:
    =======
    theta : tuple
            vector of parameters
    x     : float 
            log-moneyness
    
    returns:
    ========
    g : float 
    '''
    w_x = SVI(theta, x)

    a, b, rho, m, sigma = theta
    w_prime = b *rho +  b *(x-m) / np.sqrt((x-m)**2 + sigma**2)
    w_second = (b * sigma**2) / (((x-m)**2 + sigma **2) ** (3/2))

    g_x = 0.5 * w_second + (1- (x * w_prime) / (2*w_x)) **2 - (((w_prime)/2)**2) * (1/(w_x) + 0.25)

    return g_x

def prepareData(y, TT, r, S_0):
    '''
    Function for data preparation
    
    params:
    =======
    y : dict of DataFrames
        input data. DF columns ['C', 'K', 'sig'] : C price, K : Strike, sig : implied volatiliyu
    returns:
    ========
    temp : dict of DataFrames
           output data : DF columns ['x', 'w'] : x : log-moneyness, w : total implied variance  
    '''
    temp = {}
    for key, t in zip(y.keys(), TT):
        o = x(y[key]['K'], t, r, S_0)
        w = (y[key]['sig']**2) * t
        temp[key] = pd.DataFrame({'x': o, 'w': w})
    return temp  



def getFitted(params, data):
    
    fitted = {}
    butterfly = {}
    for t, key in zip(params.keys(), data.keys()):
        fitted[t] = pd.DataFrame({'w': SVI(params[t],data[key]['x']), 'x' : data[key]['x']})
        butterfly[t] = pd.DataFrame({'g': g(params[t],data[key]['x']), 'x' : data[key]['x']})
        
    return fitted, butterfly




def check_calendar_plot(fitted, figsize = None):
    fig = plt.figure(figsize=figsize)
    
    for key in fitted.keys():
        plt.plot(fitted[key]['x'], fitted[key]['w'])
    return None    
  
    
    
def check_butterfly_plot(butterfly, figsize = None):
    fig = plt.figure(figsize=figsize)
    
    for key in butterfly.keys():
        plt.plot(butterfly[key]['x'], butterfly[key]['g'])
    return None   


def LSQ(theta, x, w_m):
    return np.linalg.norm(SVI(theta, x) - w_m, 2)


def fitOneSlice(w_market, x_market, epsilon_g, 
                   theta_0 = None, 
                   lower = None, upper = None,  
                   automatic_bounds = True, 
                   automatic_init=True, 
                   const = None):
 
    
    if automatic_init : 
        w_min = w_market.min()
        theta_0 = (0.5 * w_min, 0.1, -0.5, 0.1, 0.1)

    if automatic_bounds :
            x_min = x_market.min()
            x_max = x_market.max()
            w_max = w_market.max()
            w_min = w_market.min()

            lower = [1e-5, 1e-2, - 0.99999, 2 * x_min , 1e-2]
            upper = [w_max, 1., 0.99999, 2 * x_max, 1]
            bounds = Bounds(lower, upper)
    else :
            assert(lower <= upper)
            bounds = Bounds(lower, upper)
        
    if const is None:
        const_ineq = {'type': 'ineq', 'fun': lambda x: g(x, x_market) - epsilon_g}
        const = [const_ineq]

    res = minimize(lambda x : LSQ(x, x_market, w_market) , theta_0, method='SLSQP', bounds=bounds, constraints=const)
    
    return res.x



def fitSliceperSlice(data, TT, epsilon_g, 
                     theta_0 = None, 
                     lower = None, upper = None,  
                     automatic_bounds = True, 
                     automatic_init=True, 
                     const = None):
    
    params = {}

    for key, t in zip(data.keys(), TT):
        w_market = data[key]['w']
        x_market = data[key]['x']
        params[t] = fitOneSlice(w_market, x_market, epsilon_g, 
                                theta_0, lower, upper, const = const)
    
    return params



def fitMultiSlice(data, TT,  epsilon_g, 
                   theta_0 = None, 
                   lower = None, upper = None,  
                   automatic_bounds = True, 
                   automatic_init=True,
                   penalBounds = None,
                   calendar = 'Default'):
    
    # Fit for the first maturity:
    params = {}
    Ts = list(data.keys())
    theta_0  = fitOneSlice(w_market = data[Ts[0]]['w'],
                          x_market = data[Ts[0]]['x'],
                          epsilon_g = epsilon_g)
    
    params[TT[0]] = theta_0
    
    x_min, x_max = penalBounds
    
    PenalGrid = np.linspace(x_min, x_max, num=200)
    
    # Number of maturities 
    N = len(TT)
    
    for date, j in zip(Ts[1:], range(1,N)):
        
        theta = fitOneSlice(w_market = data[date]['w'],
                          x_market = data[date]['x'],
                           epsilon_g = epsilon_g, theta_0 = params[TT[j-1]], automatic_init=False)
        
        try:
            
            w_current = SVI(theta, PenalGrid)
            w_before = SVI(params[TT[j-1]], PenalGrid)
            
        except KeyError:
            
            print('Problem with : theta :', theta, 'before : ', params[TT[j-1]], 'iter : ', j)
            pass
        
        # Check Calendar Spread :
        if (w_current < w_before).any():
            
            # if it exists, add Calendar constraint and refit
            eps = 5e-3
            if calendar == 'Default':
                const_ineq_cal = {'type':'ineq', 'fun': lambda x : SVI(x, PenalGrid) - (w_before + eps)}
            if calendar =='Penalty':
                const_ineq_cal = {'type':'ineq', 'fun': lambda x : np.mean(np.abs(np.clip(SVI(x, PenalGrid) - (w_before + eps), 0.0, None)))}
            
            const_ineq_butter = {'type': 'ineq', 'fun': lambda x: g(x, data[date]['x']) - epsilon_g}
            
            const = [const_ineq_butter, const_ineq_cal]
            
            # Refit :
            theta = fitOneSlice(w_market = data[date]['w'],
                                x_market = data[date]['x'],
                                epsilon_g = epsilon_g, 
                               const = const)
            
        params[TT[j]] = theta      
    return params