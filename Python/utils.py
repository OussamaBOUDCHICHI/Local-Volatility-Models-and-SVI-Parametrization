# -*- coding : utf8 -*-
# author : BOUDCHICHI Oussama
# Utils (Black-Scholes dependencies)





from os import stat
import numpy as np
import scipy.stats as stats
import pandas as pd
import matplotlib.pyplot as plt

N = stats.norm.cdf
n = stats.norm.pdf


def PriceBS(S, t, K, r, T, sig, flag='C'):

    d_1  = (np.log(S/K) + (r+(sig**2/2))*np.sqrt(T-t))  / (sig * np.sqrt(T-t))
    d_2 = d_1 - sig * np.sqrt(T-t)
    C = S * N(d_1) - K * np.exp(-r*(T-t)) * N(d_2)
    if flag =='C':
        return C
    if flag=='P':
        return  C - S + K * np.exp(-r*(T-t))
    if flag not in ['C', 'P']:
        raise ValueError('Enter C or P')  


# Greeks :

class Greeks:
    
    def __init__(self, S, t, K, r, T, sig):
        self.S = S
        self.t = t
        self.K = K
        self.r = r
        self.T = T
        self.sig = sig
        return None

    def __d_1(self):
        
        d_1  = (np.log(self.S / self.K) + (self.r+(self.sig**2/2))*np.sqrt(self.T-self.t))  / (self.sig * np.sqrt(self.T-self.t))
        return d_1
    def __d_2(self):
        
        d_2 = Greeks.__d_1(self) - self.sig * np.sqrt(self.T-self.t)
        return d_2
    

    def Delta(self):
        
        return N(Greeks.__d_1(self))

    def Gamma(self):
        
        return n(Greeks.__d_1(self)) / (self.S * self.sig * np.sqrt(self.T-self.t))

    def Theta(self):
        
        return (self.S * self.sig * n(Greeks.__d_1(self)))/ (2 * np.sqrt(self.T - self.t)) + self.r * self.K * np.exp(-self.r * (self.T - self.t)) * N(Greeks.__d_2(self))

    def Vega(self):
        
        return self.S * np.sqrt(self.T - self.t) *    n(Greeks.__d_1(self))
    
    def Rho(self):
        
        return self.K * (self.T - self.t) *  np.exp(-self.r * (self.T - self.t)) * N(Greeks.__d_2(self))

    def dCdK(self):
        
        return - np.exp(-self.r * (self.T - self.t)) * N(Greeks.__d_2(self))    


def implidVolatility(S, t, K, r, T, sig_0, 
                     C_market, maxIter = 50, 
                     tolerance = 1e-5, 
                     method = 'N-R', 
                     flag = 'C', 
                     a = .0001, b = 2.0):
    
    if method == 'N-R':

        sig = sig_0
        C = PriceBS(S, t, K, r, T, sig, flag)

        stopping_criterion = np.abs(C - C_market)
        iter = 0

        while((stopping_criterion > tolerance) & (iter < maxIter)):
            iter += 1
            Vega = Greeks(S, t, K, r, T, sig).Vega()
            
            if Vega == float(0):
                message = 'Vega equals ', 0 , 'at iteration :', iter, '. I Suggest another method. Sigma will be put to 0.'
                sig = 0.
                break
            else :
                message = 'Algorithm Converged in : ', iter, ' iterations' 

            sig = sig - (C - C_market)  / Vega
            
            C = PriceBS(S, t, K, r, T, sig, flag)
            stopping_criterion = np.abs(C - C_market)
        
        print(message)
        return sig    

    if method == 'Dichotomy':
        C_min = PriceBS(S, t, K, r, T, a, flag)
        C_max = PriceBS(S, t,K, r, T, b, flag)

        
        try:
            assert((C_min <= C_market) & (C_market <= C_max))

        except AssertionError:
            eps = 0.1
            a = np.maximum(a - 0.1, 0.001)
            b = np.minimum(b + 0.1, 3.0)
            

        sig_min = a
        sig_max = b

        sig = (sig_min + sig_max)  / 2
        C = PriceBS(S, t, K, r, T, sig, flag)
        stopping_criterion = np.abs(C - C_market)
        iter = 0

        while((stopping_criterion > tolerance) & (iter < maxIter)):
            iter += 1

            if C - C_market > 0 :
                sig_max = sig
                sig = (sig_min + sig_max) / 2
            else :
                sig_min = sig
                sig = (sig_min + sig_max) / 2
            C = PriceBS(S, t, K, r, T, sig, flag)
            stopping_criterion = np.abs(C - C_market)

        print('Algorithm Converged in : ', iter, ' iterations.')
        return sig    

def deltaHedging(S_values, K, r, ts,  T, sig, flag = 'C', plot = True, fig_size=(10,8)):

    C_0 = PriceBS(S_values[0], ts[0], K, r, T, sig, flag)
    Delta_0 = Greeks(S_values[0], ts[0], K, r, T, sig).Delta()
    Bank_avant = [0]
    Bank_apres = [C_0 - Delta_0 * S_values[0]]
    S_avant = [0]
    S_apres = [S_values[0] * Delta_0]
    VL_avant = [C_0]
    VL_apres = [C_0]
    OPT_VALUES = np.maximum(S_values - K, 0)

    Deltas = [Greeks(S, tt, K, r, T, sig).Delta() for S,tt in zip(S_values, ts)]
    
    for i in range(1,len(ts)):

        S_avant.append(Deltas[i-1] * S_values[i])
        S_apres.append(Deltas[i] * S_values[i])
        Bank_avant.append(Bank_apres[i-1] * np.exp(r * (ts[i]- ts[i-1])))
        Bank_apres.append(Bank_avant[i] + (S_avant[i] - S_apres[i]))
        VL_apres.append(S_apres[i] + Bank_apres[i])
        VL_avant.append(S_avant[i] + Bank_avant[i])

    Data = pd.DataFrame({'OPT': OPT_VALUES, 'VL_avant': VL_avant, 'VL_apres':VL_apres}, index=ts)
    if plot == True:
        plt.rcParams['figure.figsize']=fig_size
        fig, ax = plt.subplots(2)
        ax[0].plot(ts, OPT_VALUES, color='darkblue', label='Option Value', alpha=0.8, linewidth=2)
        ax[0].plot(ts, VL_avant, color='darkred', linestyle='dashed', label='Delta-Hedged Portfolio Value', linewidth=2)
        ax[1].plot(ts, OPT_VALUES/ VL_avant, color='darkgreen', label="Hedge-Ratio")
        ax[0].legend(prop={'weight':'bold'})
        ax[1].legend(prop={'weight':'bold'})

    return Data    

def getFreeRate(Cs, Ps, Ks, T):
    
    B = stats.linregress(Ks, Cs - Ps).slope
    return -np.log(-B) / T





