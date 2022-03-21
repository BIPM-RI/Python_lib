# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 16:34:53 2021

This module contains the "Procedure A" from Maurice Cox used to evaluate \
a reference value and degrees of equivalence associated with an inter-\
laboratory comparison.

@author: romain.coulon
"""

import unittest, xmlrunner
# from numpy import mean, asarray, average, sqrt, zeros
from scipy.stats import chi2
from numpy import asarray, sqrt

def CoxProcedureA(x,u,*,alpha=0.05,noFilter=True,k=2.5):
    """ This function calculate the reference value 
    following the procedure A in Cox 2002 Metrologia 39 589 
    
    Reference:
    https://doi.org/10.1088/0026-1394/39/6/10
    
    this function uses:
    _ the Chi-squared test from the module scipy.stats (scipy.py) 
    _ asarray and sqrt methods from the module numpy.py
    
    
    Parameters
    ----------
    x : list of floats
        Measurement values.
    u : list of floats
        Standard uncertainties of measurement values.
    alpha : float
        Risk of first species (Defaul = 0.05).
    noFilter: boolean
        Inactivate or activate the filtering of outliers (Default = True).

    Returns
    -------
    y : float
        Reference value.
    u_y : float
        Standard uncertainty of the reference value.
    chi2TestResult : boolean
        Result of the consistency check.
    d : list of floats
        list of degrees of equivalence.
    U_di : list of floats
        Extended uncertainty of the degrees of equivalance.
    """
    
    N=len(x)
    
    numerator_y=0
    denominator=0
    for i, ix in enumerate(x):
        numerator_y+=ix/(u[i]**2)
        denominator+=u[i]**-2
    y=numerator_y/denominator
    u_y=denominator**-0.5
    
    khi2_obs=0
    for i, ix in enumerate(x):
        khi2_obs+=(ix-y)**2/(u[i]**2)
    chi2TestResult=chi2.cdf(khi2_obs, N-1)<1-alpha
    
    if chi2TestResult or noFilter:
        d=asarray(x)-y
        U_di=2*sqrt(asarray(u)**2-asarray(u_y)**2)
    else:
        d=asarray(x)-y
        U_di=2.0*sqrt(asarray(u)**2-asarray(u_y)**2)
        x_filter=[]; u_filter=[]
        for i, di in enumerate(d):
            if abs(di)<k*U_di[i]/2:
                x_filter.append(x[i])
                u_filter.append(u[i])
            else:
                print(i)
        numerator_y=0
        denominator=0
        for i, ix in enumerate(x_filter):
            numerator_y+=ix/(u_filter[i]**2)
            denominator+=u_filter[i]**-2
        y=numerator_y/denominator
        u_y=denominator**-0.5
        d=asarray(x_filter)-y
        U_di=2*sqrt(asarray(u_filter)**2-asarray(u_y)**2)        
    print("p-value: ",1-chi2.cdf(khi2_obs, N-1)) 
    return y, u_y, chi2TestResult, d, U_di


class TestCoxProcedureA(unittest.TestCase): # test the PMM function
    def test_ProcA1(self):
        """
        TESTS FOR THE FUNCTION CoxProcedureA(x,u,*,alpha=0.05,noFilter=True,k=2.5)
        """   
        x=[43772,43910,43947,43910,43900]
        u=[97,100,104.5,160,280]
        result = CoxProcedureA(x,u)
        self.assertEqual((round(result[0]), round(result[1])) ,(43877,53))
                
if __name__ == '__main__':
    #out = io.BytesIO()
    #unittest.main()
    with open("Unit-tests_result_PMM.xml", 'w') as output:
        unittest.main(testRunner=xmlrunner.XMLTestRunner(output=output), failfast=False, buffer=False, catchbreak=False)
    f=open("Unit-tests_result_PMM.xml", 'r')
    d=f.read()
    f.close()

