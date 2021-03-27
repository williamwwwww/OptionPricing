from os import path
from datetime import datetime, date
import numpy as np
import pandas as pd
from math import log, sqrt, pi, exp, e
from scipy.stats import norm
import scipy.stats as si
from pandas import DataFrame
import json
import csv

# C(S, t) = SN(d1) - Ke^(-r(T-t))N(d2)
# P(S, t) = Ke^(-r(T-t))N(-d2) - SN(-d1)

def cal_call(S,K,t1,t2,sigma,r,q=0):
    result = S*exp(-q*(t2-t1))*norm.cdf(d1(S,K,t1,t2,sigma,r))-K*exp(-r*(t2-t1))*norm.cdf(d2(S,K,t1,t2,sigma,r))
    return result
  
def cal_put(S,K,t1,t2,sigma,r,q=0):
    result = K*exp(-r*(t2-t1))*norm.cdf(-1*d2(S,K,t1,t2,sigma,r)) - S*exp(-q*(t2-t1))*norm.cdf(-1*d1(S,K,t1,t2,sigma,r))
    return result

def d1(S,K,t1,t2,sigma,r,q=0):
    result = (log(S/K) + ((r-q)+0.5*sigma**2) * (t2-t1)) / (sigma*sqrt(t2-t1))
    return result

def d2(S,K,t1,t2,sigma,r,q=0):
    result = d1(S,K,t1,t2,sigma,r) - sigma*sqrt(t2-t1)
    return result

#Q2
def normal_random_number_generator():
    x = np.random.normal(size=(2, 200))
    return x

def generate_x_z_data(rho, size=200):
    tmp_data = normal_random_number_generator()
    x_list = tmp_data[0]
    y_list = tmp_data[1]
    z_result = []

    for i in range(size):
        z_result.append(rho*x_list[i] + sqrt(1-rho**2)*y_list[i])

    # print(z_result)
    return x_list, z_result

def cal_covariance(x_list, z_list, size=200):
    # print(len(x_list), len(z_list))
    x_mean = sum(x_list)/size
    z_mean = sum(z_list)/size

    tmp_sum = 0
    for i in range(size):
        tmp_sum += (x_list[i]-x_mean)*(z_list[i]-z_mean)

    return tmp_sum/(size-1)

def cal_sd(tmp_list, size=200):
    data_mean = sum(tmp_list)/size
    
    tmp_sum = 0
    for i in range(size):
        tmp_sum += (tmp_list[i]-data_mean)**2

    return sqrt(tmp_sum/size)

def cal_rho(theoretical_rho = 0.5):
    x_result, z_result = generate_x_z_data(theoretical_rho)
    x_z_covariance = cal_covariance(x_result, z_result)
    x_sd = cal_sd(x_result)
    z_sd = cal_sd(z_result)
    print("x_result", x_result)
    print("z_result", z_result)
    x_z_rho = x_z_covariance/(x_sd*z_sd)

    return x_z_rho

#Q3
def cal_vega(S,K,t1,t2,sigma,r,q):
    d1 = (log(S/K) + (r-q)*(t2-t1))/(sigma*sqrt(t2-t1)) + (1/2)*sigma*sqrt(t2-t1)
    d_d1 = norm.pdf(d1)
    return S*e**(-q*(t2-t1))*sqrt(t2-t1)*d_d1
    
    
def cal_implied_vol(S, K, C_true, option_type, r=0.04, q=0.2, t1=0, t2=((24-16)/365)):
    sigmahat = sqrt(2*abs((log(S/K) + (r-q)*(t2-t1)/(t2-t1))))
    sigma = sigmahat
    sigmadiff = 1
    nmax = 100
    n = 1
    tol = 1e-8
    
    while (sigmadiff >= tol and n < nmax):
        if option_type == 'C':
            C = cal_call(S,K,t1,t2,sigma,r,q)
        else:
            C = cal_put(S,K,t1,t2,sigma,r,q)
        Cvega = cal_vega(S,K,t1,t2,sigma,r,q)
        increment = (C-C_true)/Cvega
        sigma = sigma - increment
        n = n+1
        sigmadiff = abs(increment)
    if n >= nmax and sigmadiff >= tol:
        return 'NaN'
    # print(sigma)
    return sigma

def generate_bid_ask_vol_result():
    with open(path.join(path.dirname(__file__), 'mkt_date_result.json')) as json_file:
        mkt_data = json.load(json_file)
    
    fields = ['Strike', 'BidVolP', 'AskVolP', 'BidVolC', 'AskVolC']
    for t in ['31', '32', '33']:
        result_dict = {}
        for imnt in mkt_data:
            if mkt_data[imnt]['type'] == "Equity":
                pass
            else:
                strike = float(mkt_data[imnt]['strike'])
                spot = (float(mkt_data['510050'][t]['bid']) + float(mkt_data['510050'][t]['ask']))/2
            
                if strike not in result_dict:
                    result_dict[strike] = [0.0,0.0,0.0,0.0]
                if mkt_data[imnt]['option_type'] == 'P':
                    # Put Bid Vol
                    result_dict[strike][0] = cal_implied_vol(spot, strike, float(mkt_data[imnt][t]['bid']), 'P')
                    # Put Ask Vol
                    result_dict[strike][1] = cal_implied_vol(spot, strike, float(mkt_data[imnt][t]['ask']), 'P')
                elif mkt_data[imnt]['option_type'] == 'C':
                    # Call Bid Vol
                    result_dict[strike][2] = cal_implied_vol(spot, strike, float(mkt_data[imnt][t]['bid']), 'C')
                    # Call Ask Vol
                    result_dict[strike][3] = cal_implied_vol(spot, strike, float(mkt_data[imnt][t]['ask']), 'C')
                
        final_list = []
        for r in result_dict:
            result_dict[r].insert(0,r)
            final_list.append(result_dict[r])
        with open(f"{t}.csv", 'w', newline='') as f: 
            write = csv.writer(f) 
            write.writerow(fields) 
            write.writerows(final_list)

if __name__ == "__main__":
    #Q1
    # print(f"S=50, K=50, t=0, T=0.5, sigma=20%, r=1% ==> Call = {cal_call(50,50,0,0.5,0.2,0.01)}, Put = {cal_put(50,50,0,0.5,0.2,0.01)}")
    # print(f"S=50, K=60, t=0, T=0.5, sigma=20%, r=1% ==> Call = {cal_call(50,60,0,0.5,0.2,0.01)}, Put = {cal_put(50,60,0,0.5,0.2,0.01)}")
    # print(f"S=50, K=50, t=0, T=1.0, sigma=20%, r=1% ==> Call = {cal_call(50,50,0,1,0.2,0.01)}, Put = {cal_put(50,50,0,1,0.2,0.01)}")
    # print(f"S=50, K=50, t=0, T=0.5, sigma=30%, r=1% ==> Call = {cal_call(50,50,0,0.5,0.3,0.01)}, Put = {cal_put(50,50,0,0.5,0.3,0.01)}")
    # print(f"S=50, K=50, t=0, T=0.5, sigma=20%, r=2% ==> Call = {cal_call(50,50,0,0.5,0.2,0.02)}, Put = {cal_put(50,50,0,0.5,0.2,0.02)}")

    #Q2.2
    # x_z_rho = cal_rho()
    # print("x_z_rho: ", x_z_rho)

    #Q3
    # generate_bid_ask_vol_result()