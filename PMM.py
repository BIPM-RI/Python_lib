# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 10:42:33 2021

This module contains the Power Moderated Mean function used in CCRI(II).K1 key comparisons for radionuclide metrology with validation tests.

@author: Romain Coulon
email: romain.coulon@bipm.org
"""

import unittest, xmlrunner
from numpy import asarray, mean, average, sqrt, zeros

def PMM(x,u,*,autoRej=False,k=2.5,conv=1e-4): # PMM calculation
    """This function calculates the Power-Modered Mean (PMM) on a data sample.
    
    References: 
        [Accred Qual Assur (2008)13:83-89, Metrologia 52(2015)S200]
        https://link.springer.com/article/10.1007/s00769-007-0330-1
        https://iopscience.iop.org/article/10.1088/0026-1394/52/3/S200/pdf
    
    This method is applied to calculate de Key Comparison Reference Value (KCRV)
    related to international inter-laboratory comparisons of the CIPM/CCRI(II). 
    It was validated from PMM v1.xlsm provided by S. Pommé [email of the 14 Jan 2020]
    
    This function uses the method asarray, mean, average, sqrt and zeros from the module numpy.py
    
    :param x: Sample of values
    :type x: array of floats
    :param u: Sample of standard uncertainties related to the values
    :type u: array of floats
    :param r: (Optional) Set if the auto-rejection capability is activated or not (False by default)
    :type tau: bool 
    :param k: (Optional) Threshold (see as quantiles of a Normal density) of the frequentist test (default value = 2.5)
    :type k: float 
    :param conv: (Optional) Convergence paramater (set ny default equal to 1e-4)
    :type conv: float    
    :param x_ref: Estimation of the reference value
    :type x_ref: float
    :param ux_ref: Estimation of uncertainty of the reference value
    :type ux_ref: float
    :param d: Estimation of the degrees of equivalence
    :type d: array of floats
    :param ud: Estimation of the uncertainties related to the degrees of equivalence
    :type ud: array of floats    
    :param dr: Estimation of the relative degrees of equivalence
    :type dr: array of floats
    :param udr: Estimation of the uncertainties related to the relative degrees of equivalence
    :type udr: array of floats  
    :param w: Weights associated to each data point
    :type w: array of floats
    :param alpha: Value of the power paramater alpha 
    :type alpha: float
    :param birge0: Value of the Birge ratio (before adjustment) 
    :type birge0: float
    :param s: Estimation of the inter-laboratory standard deviation
    :type s: float
    :param rej: Index of the rejected points
    :type rej: list of int  

    :return y: x_ref, ux_ref, d, ud, dr, udr, w, alpha, birge0, s, rej
    :rtype y: tuple
    """
   
    def BirgeAdjust(b,x,u,s,n,*,c=1e-4):
        """This function minimize the cost function (Birge-1).
        
        :param b: Initial Birge ratio
        :type b: floats
        :param x: Sample of values
        :type x: array of floats
        :param u: Sample of standard uncertainties related to the values
        :type u: array of floats
        :param s: Initial inter-laboratory standard deviation
        :type s: floats
        :param n: Length of the data set
        :type n: int
        :param c: (Optional) Convergence paramater (set ny default equal to 1e-4)
        :type c: float
       
        :return s: Estimation of the inter-laboratory standard deviation
        :rtype s: float    
        """
        
        while (b-1)>1e-4:       # Make to Birge ratio equal to 1
            x_mp=average(x,weights=1/(u**2+s**2)) # average weighted by the inverse of the variance
            b=(1/(n-1))*sum((x-x_mp)**2/((u**2+s**2)))  # Calculation of the Birge ratio
            if (b-1)>0:
                s=s+ds      # increase the value of the inter-laboratory standart deviation 
            else:
                s=s-ds      # decrease the value of the inter-laboratory standart deviation 
            if s<0:         # detection of correltation between values
                warn("The inter-laboratory uncertainty is negative means that correlation between points could exists")
                break # break the calculation
        return s

    x=asarray(x)        # format input data
    u=asarray(u)        # format input data
    n=len(x)            # length of the data set
    alpha=2-3/n         # initialization of the power paramater alpha    
    s=0                 # initialization of the inter-laboratory standart deviation 
    ds=mean(u)*conv*0.1 # motion step of the inter-laboratory standart deviation 
    x_mp0=average(x,weights=1/(u**2+s**2)) # initial calculation of average weighted by the inverse of the variance
    birge0=sqrt((1/(n-1))*sum((x-x_mp0)**2/((u**2+s**2))))  # Calculation of initial the Birge ratio
    s=BirgeAdjust(10,x,u,s,n,c=conv) # estimation of the inter-laboratory standart deviation 
    barx=sum(x)/n      # calculation of the unweighted empirical mean
    u2_barx=sum((x-barx)**2/(n*(n-1))) # calculation of the unweighted empirical variance (external variance)
    u2_xmp=1/(sum(1/(u**2+s**2)))  # calculation of the weighted variance (internal variance)
    S=sqrt(n*max(u2_barx,u2_xmp))  # characteristic uncertainty per datum
    u2_xref=1/sum(1/((u**2+s**2)**(alpha/2)*S**(2-alpha))) # variance of the reference value
    ux_ref=sqrt(u2_xref) # standard deviation of the reference value
    w=u2_xref/((u**2+s**2)**(alpha/2)*S**(2-alpha))  # weighing factors
    x_ref=sum(w*x) # reference value
    
    rej=[]  # initial the rejection list
    if autoRej==True: # assess and reject outliers
        e=x-x_ref # euclidian distance from the reference value
        u2e_incl=u2_xref*(w**-1-1) # variance of DE including all value
        Xf=zeros(n)
        uf=zeros(n)
        xf=[]   # initial filtered data list
        uf=[]   # initial filtered data related standard deviation
        for i in range(n):
            if abs(e[i])>k*sqrt(u2e_incl[i]): # frequentist test on data
                rej.append(i) # index of the rejected point
                warn("Rejection of the point "+str(i))
            else: # the data point is kept
                xf.append(x[i])
                uf.append(u[i])
        if rej!=[]:
            xf=asarray(xf)  # format the list to array
            uf=asarray(uf)  # format the list to array
            n2=len(xf) # length of the new data set
            s=0
            alpha=2-3/n2         # recalculation of the power paramater alpha
            s=BirgeAdjust(10,xf,uf,s,n,c=conv) # estimation of the inter-laboratory standart deviation 
            barx=sum(Xf)/n2      # calculation of the unweighted empirical mean
            u2_barx=sum((xf-barx)**2/(n2*(n2-1))) # calculation of the unweighted empirical variance (external variance)
            u2_xmp=1/(sum(1/(uf**2+s**2)))  # calculation of the weighted variance (internal variance)
            S=sqrt(n2*max(u2_barx,u2_xmp))  # characteristic uncertainty per datum
            u2_xref=1/sum(1/((uf**2+s**2)**(alpha/2)*S**(2-alpha))) # variance of the reference value
            ux_ref=sqrt(u2_xref) # standard deviation of the reference value
            w=u2_xref/((uf**2+s**2)**(alpha/2)*S**(2-alpha))  # weighing factors
            x_ref=sum(w*xf) # reference value
    if s<0: # the inter-laboratory uncertainty is negative = correlation
        warn("The inter-laboratory uncertainty is negative means that correlation between points could exists")
    
    
    
    def DoE(x,u,x_ref,ux_ref,*,w=[],k=2):
        """This function aims to calculate Degrees of equivalence with Power-Modered Mean (PMM) estimation of reference values.
        
        References: 
            [Accred Qual Assur (2008)13:83-89, Metrologia 52(2015)S200]
            https://link.springer.com/article/10.1007/s00769-007-0330-1
            https://iopscience.iop.org/article/10.1088/0026-1394/52/3/S200/pdf
        
        :param x: Sample of values
        :type x: array of floats
        :param u: Sample of standard uncertainties related to the values
        :type u: array of floats
        :param x_ref: Estimation of the reference value
        :type x_ref: float
        :param ux_ref: Estimation of uncertainty of the reference value
        :type ux_ref: float
        
        :param w: (Optional) Weights associated to each data point.
        :type w: array of floats
        :param k: (Optional) Coverage factor (set by default equal to 2)
        :type k: float    
               
        :return y: d, ud, dr, udr
        :rtype y: tuple
        """
        if w==[]: w=zeros(len(x))
        x=asarray(x) # format input data
        u=asarray(u) # format input data
        w=asarray(w) # format input data
        k=2
        d=x-x_ref  # euclidian distance from the reference value
        u2d=(1-2*w)*u**2+ux_ref**2 # variance associated with DE (the weight factor is available)
        ud=k*u2d**0.5     # enlarged standard deviation associated with DoE
        dr=d/x_ref        # relative DoE
        udr=ud/x_ref      # relative u(DoE)
        return d, ud, dr, udr
    
    d, ud, dr, udr=DoE(x,u,x_ref,ux_ref,w=w)

    return x_ref, ux_ref, d, ud, dr, udr, w, alpha, birge0, s, rej

class TestPMM(unittest.TestCase): # test the PMM function
    def test_PMM1(self): # Ba133 before 2015
        """
        TESTS FOR THE FUNCTION PMM(x,u,*,autoRej=False,k=2.5,conv=1e-4)
        # Validation of PMM with PMM v1.xlsm provided by S. Pommé [email of the 14 Jan 2020]
        """   
        x=[43772,43910,43947,43910,43900,mean([43788,43773]),43910,43990,mean([43784,43769,43748]),mean([43430,43360]),44030,44300,44340,44440,43590,43880]
        u=[97,100,104.5,160,280,118,120,420,mean([59,58,58]),350,100,250,140,370,140,190]
        result = PMM(x,u)
        self.assertEqual((round(result[0]), round(result[1])) ,(43906,55))
        
    def test_PMM2(self): # Cs137 before 2015
        x=[mean([27596,27583]),27470,27380,mean([27613,27601,27615]),27628,27340,27514,mean([27255,27270]),27634,27720,27710,mean([27307,27269]),27908,mean([27914,27930]),27600]
        u=[mean([66,65]),370,220,mean([123,123,123]),99,160,84,mean([170,170]),98,140,90,mean([530,520]),169,mean([76,76]),100]
        result = PMM(x,u)
        self.assertEqual((round(result[0]), round(result[1])) ,(27618,47))
        
    def test_PMM3(self): # Sr85 in 2015
        x=[30060,30130,30000,29850,29970,29920,29890,30150,29782]
        u=[54.0,170.0,79.0,140.0,340.0,200.0,140.0,100.0,71.0]
        result = PMM(x,u)
        self.assertEqual((round(result[0]), round(result[1])) ,(29975,46)) # from past report and PMM.xlsx, accuracy 1 kBq
        
    def test_PMM4(self): # Sr85 in 2015 (strange)
        x=[29970.0,30060.0,29850.0,30000.0,29890.0,30150.0,29920.0,29782.0,30133.0]
        u=[340.0, 54.0, 140.0, 79.0, 140.0, 100.0, 200.0, 71.0, 170.0]
        result = PMM(x,u,autoRej=True)
        self.assertEqual((round(result[0]), round(result[1])) ,(29975,46))
        
if __name__ == '__main__':
    #out = io.BytesIO()
    #unittest.main()
    with open("Unit-tests_result_PMM.xml", 'w') as output:
        unittest.main(testRunner=xmlrunner.XMLTestRunner(output=output), failfast=False, buffer=False, catchbreak=False)
    f=open("Unit-tests_result_PMM.xml", 'r')
    d=f.read()
    f.close()