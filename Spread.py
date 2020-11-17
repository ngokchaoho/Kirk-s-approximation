# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 17:08:24 2020

@author: Flipped
"""


import numpy as np
from scipy.stats import norm


# Define Spread Option Class

class Spread(object):
    def __init__(self, S1, S2, K, T, r, sigma1, sigma2, rho, steps, paths):
        self.S1 = float(S1)
        self.S2 = float(S2)
        self.K = float(K)
        self.T = float(T)
        self.r = float(r)
        self.sigma1 = float(sigma1)
        self.sigma2 = float(sigma2)
        self.rho = float(rho)
        self.steps = int(steps)
        self.paths = int(paths)

    @property
    def price_mc(self):
        prices1 = np.zeros(self.paths, dtype=np.float64)
        prices1[:] = self.S1
        prices2 = np.zeros(self.paths, dtype=np.float64)
        prices2[:] = self.S2
        for t in range(0, self.steps-1):
            np.random.seed()
            z1 =np.random.standard_normal(self.paths)
            z2 =np.random.standard_normal(self.paths)
            dt = self.T/self.steps
            dw1 = np.sqrt(dt)*z1
            dw2 = np.sqrt(dt)*z2
            prices1[:] = prices1[:] * np.exp((self.r -0.5*self.sigma1**2)*dt + self.sigma1*dw1)
            prices2[:] = prices2[:] * np.exp((self.r-0.5*self.sigma2**2)*dt + self.sigma2*(self.rho*dw1+np.sqrt(1-self.rho**2)*dw2))
        return prices1, prices2


    @property
    def payoff_mc(self):
        pm1, pm2 = self.price_mc
        payoff = np.maximum(pm1 - pm2 - self.K, 0)
        return payoff

    @property
    def option_value(self):
        payoff = self.payoff_mc
        values = np.zeros(self.paths)
        values = payoff * np.exp(-self.r*self.T)
        value = np.mean(values)
        return value


    @property
    def price_kirk(self):
        z = self.S2 / (self.S2 + self.K*np.exp(-self.r*self.T))
        sigma = np.sqrt(self.sigma1**2 + self.sigma2**2*z**2  - 2*self.rho*self.sigma1*self.sigma2*z)
        d1 = (np.log(self.S1 / (self.S2 + self.K * np.exp(-self.r * self.T) ) )
              / (sigma*np.sqrt(self.T)) + 0.5*sigma*np.sqrt(self.T) )
        d2 = d1 - sigma*np.sqrt(self.T)
        price = self.S1*norm.cdf(d1) - (self.S2+self.K*np.exp(-self.r*self.T) )*norm.cdf(d2)
        return price

spread_call = Spread(100, 90, 10, 1, 0.05, 0.2, 0.3, 0.4, 365, 100000)
print (spread_call.option_value)
print (spread_call.price_kirk)