# -*- coding : utf8 -*-
# author : BOUDCHICHI Oussama
# Calls and Puts Data



import numpy as np
import pandas as pd
from yahoo_fin import options
from datetime import datetime


def getData(ticker, flag = 'C', compute_TT = True):
    
    '''
    function that retreives Calls and Puts data from : https://finance.yahoo.com/
    
    params:
    =======
    ticker     : string
                 Underlying ticker
    flag       : string
                 Call or Put
    compute_TT : bool
                 compute time to muturity or not.
    
    returns:
    ========
    y           : dict
                  dictionnary of DataFrames.  
    '''

    exp = options.get_expiration_dates(ticker)
    dt = [datetime.strptime(x, '%B %d, %Y').date() for x in exp]

    y = {}

    if flag == 'C':
        for date in dt:
            y[date] = options.get_calls(ticker, date)

    if flag == 'P':
        for date in dt:
            y[date] = options.get_puts(ticker, date)



    for key in y.keys() :
            y[key] = y[key][['Last Price', 'Strike', 'Implied Volatility']]
            y[key].columns = ['C', 'K', 'sig']
            y[key]['sig'] =   y[key]['sig'].apply(lambda x: float(x.replace('%', '').replace(',',''))) / 100
    
            


    
    if compute_TT:
        t_0 = datetime.today().date()
        TT = [(y-t_0).days / 365 for y in dt]
        return y, TT
    else:

        return y    

    
