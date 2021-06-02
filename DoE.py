# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 11:10:27 2021

This module contains the function to calculation Degrees of Equivalence in the framework of CCRI(II).K1 key comparisons for radionuclide metrology with validation tests.

@author: Romain Coulon
email: romain.coulon@bipm.org
"""

import unittest, xmlrunner
from numpy import zeros, asarray


def DoE(x,u,x_ref,ux_ref,*,w=[],k=2):
    """This function aims to calculate Degrees of equivalence.
    
    This calculation applicable in the frame of international inter-laboratory comparisons of the CIPM/CCRI(II).
    This calculation is associated with the use of the Power-Modered Mean (PMM) estimation of reference values.
    
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
    
    :param d: Estimation of the degrees of equivalence
    :type d: array of floats
    :param ud: Estimation of the uncertainties related to the degrees of equivalence
    :type ud: array of floats    
    :param dr: Estimation of the relative degrees of equivalence
    :type dr: array of floats
    :param udr: Estimation of the uncertainties related to the relative degrees of equivalence
    :type udr: array of floats  
    
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


class TestDoE(unittest.TestCase): # test the DoE function
    def test_DoE1(self): # Ba-133
        """
        TESTS FOR THE FUNCTION DoE(x,u,x_ref,ux_ref,*,w=[],k=2)
        """   
        x=[43772,43910,44030,44340,44440,43590,43880,44060,43920,43304,43881,44115,43941]
        u=[97,120,100,140,370,140,190,300,100,255,263,157,695]
        d=[16154, 16292, 16412, 16722, 16822, 15972, 16262, 16442, 16302, 15686, 16263, 16497, 16323]
        ud=[215.5736533067063, 257.75181861628056, 220.9886874932742, 295.357410606201, 745.9463787699489, 295.357410606201, 391.4537009660274, 607.3186972257647, 220.9886874932742, 518.5903971343859, 534.3332293616036, 327.76821078316914, 1393.1747916180511]
        result = DoE(x,u,27618,47)        
        self.assertEqual((result[0].tolist(), result[1].tolist()) ,(d,ud))
        
    def test_DoE2(self): # Sr-85 in 2015        
        x=[30060.0, 30086.0, 30130.0, 30133.0]
        u=[54.0, 91.0, 150.0, 170.0]
        n=len(x)
        w=[0.2133,0,0, 0.06695]
        d=[0.085,0.11,0.16,0.16] # from past report and PMM.xlsx
        ud=[0.12,0.20,0.31,0.33] # from past report and PMM.xlsx
        result = DoE(x,u,29974.862,45.62,w=w)
        for i in range(n):
            self.assertAlmostEqual(round(result[0][i],2),1000*d[i],places=-1)         # accuracy 10 kBq
            self.assertAlmostEqual(round(result[1][i],2),1000*ud[i],places=-1)        # accuracy 10 kBq
            
    def test_DoE3(self): # Sr-85 in 2020        
        x=[30086.0, 30130.0, 30133.0, 30177.0]
        u=[91.0, 150.0, 170.0, 108.0]
        n=len(x)
        w=[0, 0, 0.07881554190977942, 0.13189462767657084]
        d=[0.10, 0.15, 0.15, 0.19] # from past report and PMM.xlsx
        ud=[0.21, 0.32, 0.33, 0.21] # from past report and PMM.xlsx
        result = DoE(x,u,29983.06051896807,52.354677501912256,w=w)
        for i in range(n):
            self.assertAlmostEqual(round(result[0][i],2),1000*d[i],places=-1)         # accuracy 10 kBq
            self.assertAlmostEqual(round(result[1][i],2),1000*ud[i],places=-1)        # accuracy 10 kBq
            
    def test_DoE4(self): # Sr-85 in 2020 with round        
        x=[30086.0, 30130.0, 30133.0, 30180.0]
        u=[91.0, 150.0, 170.0, 110.0]
        n=len(x)
        w=[0, 0, 0.07881554190977942, 0.13189462767657084]
        d=[0.10, 0.15, 0.15, 0.20] # from past report and PMM.xlsx
        ud=[0.21, 0.32, 0.33, 0.22] # from past report and PMM.xlsx
        result = DoE(x,u,29983.06051896807,52.354677501912256,w=w)
        for i in range(n):
            self.assertAlmostEqual(round(result[0][i],2),1000*d[i],places=-1)         # accuracy 10 kBq
            self.assertAlmostEqual(round(result[1][i],2),1000*ud[i],places=-1)        # accuracy 10 kBq
            
    def test_DoE5(self): # Na-22 in 2020 with round        
        x=[7570.0, 7490.0, 7523.0, 7533.0, 7568.0, 7523.0]
        u=[27.0, 31.0, 29.0, 17.0, 14.0, 15.0]
        n=len(x)
        #w=[0.060583292097284885, 0, 0.060583292097284885, 0.13493713184451323, 0.17401885385593027, 0.1594362447191583]
        w=[0.060583292097284885, 0, 0.060583292097284885, 0.13493713184451323, 0.17401885385593027, 0.1594362447191583]
        w=[0, 0, 0, 0, 0, 0]
        d=[36, -44, -11, -1, 34, -11] # from PMM.xlsx calculted by CM
        ud=[56, 64, 60, 37, 32, 33] # from PMM.xlsx calculated by CM
        result = DoE(x,u,7534.047302022265,7.2856962130888485,w=w)
        for i in range(n):
            self.assertAlmostEqual(round(result[0][i]),d[i],places=0)         # accuracy 1 kBq
            self.assertAlmostEqual(round(result[1][i]),ud[i],places=0)        # accuracy 1 kBq
            
if __name__ == '__main__':
    #out = io.BytesIO()
    #unittest.main()
    with open("Unit-tests_result_DoE.xml", 'w') as output:
        unittest.main(testRunner=xmlrunner.XMLTestRunner(output=output), failfast=False, buffer=False, catchbreak=False)
    f=open("Unit-tests_result_DoE.xml", 'r')
    d=f.read()
    f.close()