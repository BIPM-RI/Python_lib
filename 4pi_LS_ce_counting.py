# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 08:51:43 2021

@author: romain.coulon
"""

import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
sys.path.insert(1, 'G:\Python_modules\BIPM_RI_PyModules')
import DataProcessing as dp
import TDCRLabZy as td
from scipy.optimize import curve_fit
import xml.etree.ElementTree as ET 


pathDir="G:/Activity_mesurements/LS_spectrometry/"

listMode=True


### Mass measurement
m=[0.035371,0.047671,0.059204,0.072800] # /g


epsilon = 1

#k_gamma=0.997  # Kossert
k_gamma=0.999320 # Coulon
u_k_gamma=0.000018 # Coulon

p_ce = 0.96337 # Kossert
u_p_ce = 0.00033  # Kossert

#p_ce = (41.8+44.1+9.04+1.413)/100      # DDEP
#u_p_ce=np.sqrt(0.8**2+0.9**2+0.19**2+0.029**2)/100 # DDEP

T05=461.9 # days Monographie BIPM-5 vol. 8
u_T05=0.4
ref_date=dt.date(2021,9,1)
meas_date=dt.date(2021,12,14)
days_lasted=(meas_date-ref_date).days
time_lasted=days_lasted*24*3600
Lambda=np.log(2)/(T05*24*3600)
Q=np.exp(Lambda*time_lasted)
uQ=np.sqrt((-Q*np.log(Q)/T05)**2*u_T05**2)





### TRICARB

folder="Tricarb-Cd-109bis/"

def readTricarbSpectrum(path,fileN,m,e_cal):
    f=open(path+fileN, "r") 
    df=f.read() # read the file
    f.close()
    df=df.split("\n") # split the lines
    TopData=[]
    for i, x in enumerate(df):
        if i>0: TopData.append(x.split(","))
        else: 
            header=x.split(",")
    
#    print(TopData[0][0])
#    TR=float(TopData[0][0][10:])*60 # s
    TR=20*60 # s
    S=[]
    for i in range(4,len(TopData)):
        if int(TopData[i][0])==-1: break
        S.append(int(TopData[i][0]))
        
    E=np.arange(0,len(S),1)
    e_cal=1
    E=E*e_cal   
    
    f=open(path+"Prot.dat", "r") 
    df=f.read() # read the file
    f.close()
    
    t_flag=df.find("DATE")
    meas_date=dt.date(int(df[t_flag+11:t_flag+13])+2000,int(df[t_flag+8:t_flag+10]),int(df[t_flag+5:t_flag+7]))
    ref_date=dt.date(2021,9,1)
    days_lasted=(meas_date-ref_date).days
    
    t_flag=df.find("TIME")
    hour_lasted=int(df[t_flag+5:t_flag+7])-12
    min_mes=int(df[t_flag+8:t_flag+10])
    
    time_lasted=days_lasted*24*3600+hour_lasted*3600+min_mes*60
    
    Lambda=np.log(2)/(461.9*24*3600)
    
    ktm=(Lambda*TR)/(1-np.exp(-Lambda*TR))
    Qtm=np.exp(Lambda*time_lasted)
        
    S=np.asarray(S)*ktm*Qtm/(TR*m)
    return E, S

e_cal=1

E_4, S_4 = readTricarbSpectrum(pathDir+folder,"Cd-109_4.Spectrum",m[3],e_cal)
S_4 = S_4 - readTricarbSpectrum(pathDir+folder,"Blank.Spectrum",m[3],e_cal)[1]

E_3, S_3 = readTricarbSpectrum(pathDir+folder,"Cd-109_3.Spectrum",m[2],e_cal)
S_3 = S_3 - readTricarbSpectrum(pathDir+folder,"Blank.Spectrum",m[2],e_cal)[1]

E_2, S_2 = readTricarbSpectrum(pathDir+folder,"Cd-109_2.Spectrum",m[1],e_cal)
S_2 = S_2 - readTricarbSpectrum(pathDir+folder,"Blank.Spectrum",m[1],e_cal)[1]

E_1, S_1 = readTricarbSpectrum(pathDir+folder,"Cd-109_1.Spectrum",m[0],e_cal)
S_1 = S_1 - readTricarbSpectrum(pathDir+folder,"Blank.Spectrum",m[0],e_cal)[1]


#### DT5730S

#def readDT5730S(folder,channel,m):
#    f=open(folder+"CH"+str(channel)+"@DT5730SB_14158_E"+".txt", "r") 
#    df2=f.read() # read the file
#    f.close()
#    
#    f = open(pathDir+folder[:-4]+"run.info",'r')
#    run_file=f.read() # read the file
#    f.close() # close the file
#    t_flag=run_file.find("time.real=")
#    TR=float(run_file[t_flag+10:t_flag+12])*3600+float(run_file[t_flag+13:t_flag+15])*60+float(run_file[t_flag+16:t_flag+18])
#    flg1=run_file.find(str(channel)+".throughput")
#    throughput=float(run_file[flg1+13:flg1+20])
#    flg1=run_file.find(str(channel)+".icr")
#    icr=float(run_file[flg1+6:flg1+13])
#    flg1=run_file.find(str(channel)+".ocr")
#    ocr=float(run_file[flg1+6:flg1+13])
#
#    
#    TL=TR#*throughput/icr
#    
#    t_mes=run_file.find("time.start=")
#    min_mes=float(run_file[t_mes+25:t_mes+27])
#    sec_mes=float(run_file[t_mes+28:t_mes+30])
#    ref_date=dt.date(2021,9,1)
#    meas_date=dt.date(int(run_file[t_mes+11:t_mes+15]),int(run_file[t_mes+16:t_mes+18]),int(run_file[t_mes+19:t_mes+21]))
#    days_lasted=(meas_date-ref_date).days
#    hour_lasted=int(run_file[t_mes+22:t_mes+24])-1-12
#    
#    time_lasted=days_lasted*24*3600+hour_lasted*3600+min_mes*60+sec_mes
#    
#    Lambda=np.log(2)/(461.9*24*3600)
#    
#    ktm=(Lambda*TR)/(1-np.exp(-Lambda*TR))
#    Qtm=np.exp(Lambda*time_lasted)
#    
#    Slist=df2.split('\n')
#    S=[]
#    for i in Slist:
#        if i=='': break
#        S.append(float(i))
#    
#    E=np.arange(0,len(S),1)
#    e_cal=1
#    E=E*e_cal 
#    S=np.asarray(S)*ktm*Qtm/(TL*m)
#    return E, S

#def readDT5730S(folder,channel,m):
#    f=open(folder+"CH"+str(channel)+"@DT5730SB_14158_E"+".txt", "r") 
#    df2=f.read() # read the file
#    f.close()
#    
#    f = open(pathDir+folder[:-4]+"run.info",'r')
#    run_file=f.read() # read the file
#    f.close() # close the file
#    t_flag=run_file.find("time.real=")
#    TR=float(run_file[t_flag+10:t_flag+12])*3600+float(run_file[t_flag+13:t_flag+15])*60+float(run_file[t_flag+16:t_flag+18])
#    flg1=run_file.find(str(channel)+".throughput")
#    throughput=float(run_file[flg1+13:flg1+20])
#    flg1=run_file.find(str(channel)+".icr")
#    icr=float(run_file[flg1+6:flg1+13])
#    flg1=run_file.find(str(channel)+".ocr")
#    ocr=float(run_file[flg1+6:flg1+13])
#
#    
#    TL=TR#*throughput/icr
#    
#    t_mes=run_file.find("time.start=")
#    min_mes=float(run_file[t_mes+25:t_mes+27])
#    sec_mes=float(run_file[t_mes+28:t_mes+30])
#    ref_date=dt.date(2021,9,1)
#    meas_date=dt.date(int(run_file[t_mes+11:t_mes+15]),int(run_file[t_mes+16:t_mes+18]),int(run_file[t_mes+19:t_mes+21]))
#    days_lasted=(meas_date-ref_date).days
#    hour_lasted=int(run_file[t_mes+22:t_mes+24])-1-12
#    
#    time_lasted=days_lasted*24*3600+hour_lasted*3600+min_mes*60+sec_mes
#    
#    Lambda=np.log(2)/(461.9*24*3600)
#    
#    ktm=(Lambda*TR)/(1-np.exp(-Lambda*TR))
#    Qtm=np.exp(Lambda*time_lasted)
#    
#    Slist=df2.split('\n')
#    S=[]
#    for i in Slist:
#        if i=='': break
#        S.append(float(i))
#    
#    E=np.arange(0,len(S),1)
#    e_cal=1
#    E=E*e_cal 
#    S=np.asarray(S)*ktm*Qtm/(TL*m)
#    return E, S


def readDT5730Sn42(folder,channel,m):
    tree = ET.parse(folder+"CH"+str(channel)+"@DT5730SB_14158_E.n42")
    root = tree.getroot()
    t_mes=root[3][1].text
    TR=root[3][2].text
    TR=float(TR[2])*3600+float(TR[4:6])*60+float(TR[7:-1])
    TL=root[3][3][0].text
    TL=float(TL[2])*3600+float(TL[4:6])*60+float(TL[7:-1])
    print(TL,TR)
    ref_date=dt.date(2021,9,1)
    
    meas_date=dt.date(int(t_mes[:4]),int(t_mes[5:7]),int(t_mes[8:10]))
    days_lasted=(meas_date-ref_date).days
    hour_lasted=int(t_mes[11:13])-12
    min_mes=float(t_mes[14:16])-0
    sec_mes=float(t_mes[17:19])-0
    
    time_lasted=days_lasted*24*3600+hour_lasted*3600+min_mes*60+sec_mes
    
    Lambda=np.log(2)/(461.9*24*3600)
    
    ktm=(Lambda*TR)/(1-np.exp(-Lambda*TR))
    Qtm=np.exp(Lambda*time_lasted)
    
    
    Spectrum=root[3][3][1].text
    S2=Spectrum.replace("\n","")
    Slist=S2.split(' ')
    S=[]
    for i in Slist:
        if i=='': break
        S.append(float(i))
    
    E=np.arange(0,len(S),1)
    e_cal=1
    E=E*e_cal
    S=np.asarray(S)*ktm*Qtm/(TR*m)
    return E, S





E_1_c1, S_1_c1=readDT5730Sn42("DT5730S-Cd-109-1ter/RAW/",0,m[0])
E_0_c1, S_0_c1=readDT5730Sn42("DT5730S-blank-bis/RAW/",0,m[0])
S_1_c1=S_1_c1-S_0_c1
E_1_c2, S_1_c2=readDT5730Sn42("DT5730S-Cd-109-1ter/RAW/",1,m[0])
E_0_c2, S_0_c2=readDT5730Sn42("DT5730S-blank-bis/RAW/",1,m[0])
S_1_c2=S_1_c2-S_0_c2
E_1_c3, S_1_c3=readDT5730Sn42("DT5730S-Cd-109-1ter/RAW/",2,m[0])
E_0_c3, S_0_c3=readDT5730Sn42("DT5730S-blank-bis/RAW/",2,m[0])
S_1_c3=S_1_c3-S_0_c3

E_2_c1, S_2_c1=readDT5730Sn42("DT5730S-Cd-109-2ter/RAW/",0,m[1])
E_0_c1, S_0_c1=readDT5730Sn42("DT5730S-blank-bis/RAW/",0,m[1])
S_2_c1=S_2_c1-S_0_c1
E_2_c2, S_2_c2=readDT5730Sn42("DT5730S-Cd-109-2ter/RAW/",1,m[1])
E_0_c2, S_0_c2=readDT5730Sn42("DT5730S-blank-bis/RAW/",1,m[1])
S_2_c2=S_2_c2-S_0_c2
E_2_c3, S_2_c3=readDT5730Sn42("DT5730S-Cd-109-2ter/RAW/",2,m[1])
E_0_c3, S_0_c3=readDT5730Sn42("DT5730S-blank-bis/RAW/",2,m[1])
S_2_c3=S_2_c3-S_0_c3


E_3_c1, S_3_c1=readDT5730Sn42("DT5730S-Cd-109-3ter/RAW/",0,m[2])
E_0_c1, S_0_c1=readDT5730Sn42("DT5730S-blank-bis/RAW/",0,m[2])
S_3_c1=S_3_c1-S_0_c1
E_3_c2, S_3_c2=readDT5730Sn42("DT5730S-Cd-109-3ter/RAW/",1,m[2])
E_0_c2, S_0_c2=readDT5730Sn42("DT5730S-blank-bis/RAW/",1,m[2])
S_3_c2=S_3_c2-S_0_c2
E_3_c3, S_3_c3=readDT5730Sn42("DT5730S-Cd-109-3ter/RAW/",2,m[2])
E_0_c3, S_0_c3=readDT5730Sn42("DT5730S-blank-bis/RAW/",2,m[2])
S_3_c3=S_3_c3-S_0_c3

E_4_c1, S_4_c1=readDT5730Sn42("DT5730S-Cd-109-4ter/RAW/",0,m[3])
E_0_c1, S_0_c1=readDT5730Sn42("DT5730S-blank-bis/RAW/",0,m[3])
S_4_c1=S_4_c1-S_0_c1
E_4_c2, S_4_c2=readDT5730Sn42("DT5730S-Cd-109-4ter/RAW/",1,m[3])
E_0_c2, S_0_c2=readDT5730Sn42("DT5730S-blank-bis/RAW/",1,m[3])
S_4_c2=S_4_c2-S_0_c2
E_4_c3, S_4_c3=readDT5730Sn42("DT5730S-Cd-109-4ter/RAW/",2,m[3])
E_0_c3, S_0_c3=readDT5730Sn42("DT5730S-blank-bis/RAW/",2,m[3])
S_4_c3=S_4_c3-S_0_c3



# read nanoTDCR spectrum

def readNanoTDCR(path,m):
    
    f=open(path, "r") # open the file
    d=f.read() # read the file
    f.close() # close the file
    t1=d.find("  <data>\n") # find the start flag key word
    t2=d.find("  </data>") # find the end flag key word
    C=(d[t1+9:t2-1]) # Counting data
    C=C.split("\n") # Separate the bins
    n=len(C) # number of bins
    S=np.empty(n) # define an empty vector
    for x in range(n):
        S[x]=int(C[x]) # convert string to integer
    
    
    t3=d.find("<real>")
    t4=d.find("</real>")
    TR=float(d[t3+6:t4])
    t5=d.find("<live>")
    t6=d.find("</live>")
    TL=float(d[t5+6:t6])
    t7=d.find("<date>")
    t8=d.find("</date>")
    t_mes=d[t7+6:t8]
    
    ref_date=dt.date(2021,9,1)
    
    meas_date=dt.date(int(t_mes[6:10]),int(t_mes[3:5]),int(t_mes[:2]))
    days_lasted=(meas_date-ref_date).days
    hour_lasted=int(t_mes[12:14])-12
    min_mes=float(t_mes[15:17])-0
    sec_mes=float(t_mes[18:20])-0
    
    time_lasted=days_lasted*24*3600+hour_lasted*3600+min_mes*60+sec_mes
    
    Lambda=np.log(2)/(461.9*24*3600)
    
    ktm=(Lambda*TR)/(1-np.exp(-Lambda*TR))
    Qtm=np.exp(Lambda*time_lasted)
        
    E=np.arange(0,len(S),1)
    e_cal=1
    E=E*e_cal
    print(sum(S))
    S=np.asarray(S)*ktm*Qtm/(TL*m)
    return E, S


E_4_td, S_4_td=readNanoTDCR("nanoTDCR/SpectrumS4_ABC.lts",m[3])
E_0_td, S_0_td=readNanoTDCR("nanoTDCR/Blank.lts",m[3])
S_4_td=S_4_td-S_0_td
#ur_bck=100*np.sqrt(sum(S_0_td[1000:]))/sum(S_4_td)



E_3_td, S_3_td=readNanoTDCR("nanoTDCR/SpectrumS3_ABC.lts",m[2])
E_0_td, S_0_td=readNanoTDCR("nanoTDCR/Blank.lts",m[2])
S_3_td=S_3_td-S_0_td

E_2_td, S_2_td=readNanoTDCR("nanoTDCR/SpectrumS2_ABC.lts",m[1])
E_0_td, S_0_td=readNanoTDCR("nanoTDCR/Blank.lts",m[1])
S_2_td=S_2_td-S_0_td

E_1_td, S_1_td=readNanoTDCR("nanoTDCR/SpectrumS1_ABC.lts",m[0])
E_0_td, S_0_td=readNanoTDCR("nanoTDCR/Blank.lts",m[0])
S_1_td=S_1_td-S_0_td
ur_bck=100*(sum(S_0_td[1000:3200]))/sum(S_1_td[1000:3200])

## LIST mode analysis

def resetTriplet():
    return [False,False,False]


def listMode(file):
    

    
    SumW=3
    baseTime = 1.0888302125545235e-12    
    trigger_time=0
    comdtW=0
    thres=0
    ExtDT=10*1e-6/baseTime
    ResolTime=50*1e-9/baseTime
    E_list_A_raw=[]
    E_list_B_raw=[]
    E_list_C_raw=[]
    triplet=[False, False, False]
    Acount=0
    Bcount=0
    Ccount=0
    Scount=0
    Dcount=0
    Tcount=0
    ABcount=0
    BCcount=0
    ACcount=0
    
    Arate=[]
    Brate=[]
    Crate=[]
    ABrate=[]
    BCrate=[]
    ACrate=[]
    Srate=[]
    Drate=[]
    Trate=[]
    Elist=[]
    
    liveTime=0
    countrejpulse=0
    busy=0
    w=1

    with open(file) as csv_file:
        csv_k = csv.reader(csv_file, delimiter=';')
        for row in csv_k:
            if row[0] != "BOARD":
    #            print(row)
                # ALL EVENTS
                charge=float(row[SumW])  # read the pulse fast charge
                #t_memo=t
                t=int(row[2])      # arrival time of the event
                realTime=t*baseTime
                          
    #            realTime+=t-t_memo
                if int(row[1])==0: E_list_A_raw.append(charge*0.071) # record charge to build channel A spectrum 
                if int(row[1])==1: E_list_B_raw.append(charge*0.067) # record charge to build channel B spectrum 
                if int(row[1])==2: E_list_C_raw.append(charge*0.088) # record charge to build channel C spectrum 
                
    #            if row[5] == "0xc480": print("saturation",row[3])
                
                # (1) DETECT EVENT WHEN IDLE
                if charge>thres and t>trigger_time+comdtW:
                    # RECORD THE PREVIOUS COUNT                
                    if sum(triplet)>=1:
                        Scount+=1 # record the single events
                        
                        if triplet==[True, False, False]:
                            Acount+=1 # record event in channel A only
    #                        Elist.append(E_list_A_raw[-1])
                        if triplet==[False, True, False]:
                            Bcount+=1 # record event in channel B only
    #                        Elist.append(E_list_B_raw[-1])
                        if triplet==[False, False, True]:
                            Ccount+=1 # record event in channel C only
    #                        Elist.append(E_list_C_raw[-1])
                    if sum(triplet)>=2:
                        Dcount+=1 # record the logic-sum of double coincidences
    #                    Elist.append(E_list_A_raw[-1]+E_list_B_raw[-1]+E_list_C_raw[-1])
                        if triplet==[True, True, False]:
                            ABcount+=1 # record coincident event in channel A and B
                            Elist.append(E_list_A_raw[-1]+E_list_B_raw[-1])
                        if triplet==[False, True, True]:
                            BCcount+=1 # record coincident event in channel B and C
                            Elist.append(E_list_B_raw[-1]+E_list_C_raw[-1])
                        if triplet==[True, False, True]:
                            ACcount+=1 # record coincident event in channel A and C
                            Elist.append(E_list_A_raw[-1]+E_list_C_raw[-1])
                        if sum(triplet)==3:
                            Tcount+=1; ABcount+=1; BCcount+=1; ACcount+=1 # record the triple coincidences
                            Elist.append(E_list_A_raw[-1]+E_list_B_raw[-1]+E_list_C_raw[-1])
                    
                    # liveTime+=(t-(trigger_time+comdtW+ExtDT))*baseTime
                    liveTime+=(t-(trigger_time+comdtW+countrejpulse*busy))*baseTime
                    
                    # print(" *** ")
                    # print("InterPulse dead time = ",round((comdtW+countrejpulse*busy)*baseTime*1e6,3)," µs")
                    # print("InterPulse real time = ",round((t-trigger_time)*baseTime*1e6,3)," µs")
                    # print(countrejpulse)
                                   
                    if realTime>t_mes*w: # end of a run
                        
                        print("\n****")
                        print(int(realTime/60)," min processed.")
                        print("run #",w)
                        print('****')
                        print("A,B,C ",round(Acount/liveTime,1),round(Bcount/liveTime,1),round(Ccount/liveTime,1))
                        print("AB, BC, AC", round(ABcount/liveTime,1),round(BCcount/liveTime,1),round(ACcount/liveTime,1))
                        print("S, D, T", round(Scount/liveTime,1), round(Dcount/liveTime,1), round(Tcount/liveTime,1))
                        print("Check sum = ",ABcount+BCcount+ACcount-2*Tcount==Dcount)
                        print("TDCR = ",round(Tcount/(Dcount+1),3))
                        print("Dead time = ", round(100*(1-liveTime/realTime)), " %")
                        print("LiveTime = ", round(liveTime,2), " s")
                        print("Real Time = ", round(realTime,2), " s")
        
        
        
                        
                        # CALCULATE COUNT RATE
                        Arate.append(Acount/liveTime)
                        Brate.append(Bcount/liveTime)
                        Crate.append(Ccount/liveTime)
                        ABrate.append(ABcount/liveTime)
                        BCrate.append(BCcount/liveTime)
                        ACrate.append(ACcount/liveTime)
                        Srate.append(Scount/liveTime)
                        Drate.append(Dcount/liveTime)
                        Trate.append(Tcount/liveTime)
                                   
    #                    Acount=0
    #                    Bcount=0
    #                    Ccount=0
    #                    ABcount=0
    #                    BCcount=0
    #                    ACcount=0
    #                    Scount=0
    #                    Dcount=0
    #                    Tcount=0
    #                    
    #                    realTime=0
    #                    liveTime=0
                                        
                        w+=1
                    
                    # RESTART FOR A NEW EVENT ANALYSIS
                    triplet=resetTriplet() # reset the coincidence detection vector
                    trigger_time=t # triggering time of the detected event in ns
                    triplet[int(row[1])]=True # record of the event
                    comdtW=ExtDT
                    countrejpulse=1
                # (2) DETECT EVENT WHEN BUSY
                elif charge>thres:
                    #(2-1) DETECT EVENT IN THE COINCIDENCE WINDOW
                    if t<trigger_time+ResolTime:
                        triplet[int(row[1])]=True # record of the event
                        comdtW+=t-trigger_time
                    #(2-1) DETECT EVENT OUT OF THE COINCIDENCE WINDOW
                    else:
                        countrejpulse+=1
                        # print((t-trigger_time)*baseTime*1e6)
                        comdtW=t-trigger_time+ExtDT
    #            else:
    #                print("Warning: underfined status")
            else: print(row)

    return  E_list_A_raw, E_list_B_raw, E_list_C_raw, Elist



def listHisto(E_list_A, E_list_B, E_list_C, E_list, liveTime):
    plt.figure("Spectra raw")
    plt.clf()
    sA_r=plt.hist(np.asarray(E_list_A),bins=4000,label="channel A raw")[0]
    sB_r=plt.hist(np.asarray(E_list_B),bins=4000,label="chanel B raw")[0]
    sC_r=plt.hist(np.asarray(E_list_C),bins=4000,label="channel C raw")[0]
    s_r=plt.hist(np.asarray(E_list),bins=4000,label="all coinc")[0]
    plt.legend()
    plt.show()
    
    def filterHF(x):
        xf=[]
        k=5
        dth=10
        s=10000
        for i, y in enumerate(x):
            if i>dth:
                s=np.std(x[i-dth:i])
            if y-x[i-1]>k*s:
                xf.append(np.random.normal(x[i-1],s))
            else:
                xf.append(y)
        return xf
    
    sA_r=filterHF(sA_r)
    sB_r=filterHF(sB_r)
    sC_r=filterHF(sC_r)
    s_r=filterHF(s_r)
        
    plt.figure("Spectra raw 2")
    plt.clf()
    plt.plot(sA_r,label="channel A raw")
    plt.plot(sB_r,label="channel B raw")
    plt.plot(sC_r,label="channel C raw")
    plt.plot(s_r,label="all coinc")
    plt.legend()
    plt.show()

#E_list_A, E_list_B, E_list_C, E_list=listMode("DT5730S-Cd-109-4/RAW/SDataR_Cd-109-1-0_17.csv")




def findValley(S,l,n):
    dS=S[l:]-S[0:-l]
    dS=dS.tolist()+np.zeros(l).tolist()
    countZcrossing=0
    dSm=1
    h=0
    for i, s in enumerate(dS):
        if s>0:
            if dSm==-1:
                countZcrossing+=1
            dSm=1
        else:
            dSm=-1
        if countZcrossing==n:
            h=i
    return h+1, dS

def extrapolate(S,c,l):
    slope=(S[c+l]-S[c])/(l)
    intercept=S[c]-slope*(c)
    if intercept>0:
        nb=0.5*slope*c**2+intercept*c
    else:
        nb=0.5*slope*c**2+intercept*c+intercept**2/(2*slope)
        print("error")
    return slope, intercept, nb

def cdFct(x, A, x0, sigma0, B, x1, sigma1):
    return A * np.exp(-(x - x0) ** 2 / (2 * sigma0 ** 2)) + B * np.exp(-(x - x1) ** 2 / (2 * sigma1 ** 2)) 

def cdFct2(x, A, x0, sigma0, B, x1, sigma1, C, x2, sigma2):
    return A * np.exp(-(x - x0) ** 2 / (2 * sigma0 ** 2)) + B * np.exp(-(x - x1) ** 2 / (2 * sigma1 ** 2)) + C * np.exp(-(x - x2) ** 2 / (2 * sigma2 ** 2)) 

def gauss(x, A, x0, sigma0):
    return A * np.exp(-(x - x0) ** 2 / (2 * sigma0 ** 2))




"TRICARB CALCULATION"

#c_1,dS1=findValley(S_1,1,1)
#c_2,dS2=findValley(S_2,1,1)
#c_3,dS3=findValley(S_3,1,1)
#c_4,dS4=findValley(S_4,1,1)
#c_1=47
#c_2=45
#c_3=45
#c_4=42
#na_1=sum(S_2[c_1:200])
#na_2=sum(S_2[c_2:200])
#na_3=sum(S_3[c_3:200])
#na_4=sum(S_3[c_4:200])
#
#slope_1,intercept_1,nb_1=extrapolate(S_1,c_1,2)
#slope_2,intercept_2,nb_2=extrapolate(S_2,c_2,2)
#slope_3,intercept_3,nb_3=extrapolate(S_3,c_3,2)
#slope_4,intercept_4,nb_4=extrapolate(S_4,c_4,2)

#A_1=(na_1+nb_1)*k_tm*k_gamma/(epsilon*p_ce)
#A_2=(na_2+nb_2)*k_tm*k_gamma/(epsilon*p_ce)
#A_3=(na_3+nb_3)*k_tm*k_gamma/(epsilon*p_ce)
#A_4=(na_4+nb_4)*k_tm*k_gamma/(epsilon*p_ce)



guess=np.array([sum(S_1[10:50]),23,8,sum(S_1[50:250]),100,40])
fit_1=curve_fit(cdFct,E_1[10:250],S_1[10:250],guess)
fit_1_c=gauss(E_1,fit_1[0][0],fit_1[0][1],fit_1[0][2])
fit_1_b=gauss(E_1,fit_1[0][3],fit_1[0][4],fit_1[0][5])
na_1=sum(fit_1_b)

fit_2=curve_fit(cdFct,E_2[10:250],S_2[10:250],guess)
fit_2_c=gauss(E_2,fit_2[0][0],fit_2[0][1],fit_2[0][2])
fit_2_b=gauss(E_2,fit_2[0][3],fit_2[0][4],fit_2[0][5])
na_2=sum(fit_2_b)

fit_3=curve_fit(cdFct,E_3[10:250],S_3[10:250],guess)
fit_3_c=gauss(E_3,fit_3[0][0],fit_3[0][1],fit_3[0][2])
fit_3_b=gauss(E_3,fit_3[0][3],fit_3[0][4],fit_3[0][5])
na_3=sum(fit_3_b)

fit_4=curve_fit(cdFct,E_4[10:250],S_4[10:250],guess)
fit_4_c=gauss(E_4,fit_4[0][0],fit_4[0][1],fit_4[0][2])
fit_4_b=gauss(E_4,fit_4[0][3],fit_4[0][4],fit_4[0][5])
na_4=sum(fit_4_b)







A_1=(na_1)*k_gamma/(epsilon*p_ce)
A_2=(na_2)*k_gamma/(epsilon*p_ce)
A_3=(na_3)*k_gamma/(epsilon*p_ce)
A_4=(na_4)*k_gamma/(epsilon*p_ce)

print("\n**** TRICARB ***")
print("SOURCE #1")
print("CE emission = ",round(na_1),' s-1 g-1')
print("Activity = ",round(A_1)," Bq g-1")

print("SOURCE #2")
print("CE emission = ",round(na_2),' s-1 g-1')
print("Activity = ",round(A_2)," Bq g-1")

print("SOURCE #3")
print("CE emission = ",round(na_3),' s-1 g-1')
print("Activity = ",round(A_3)," Bq g-1")

print("SOURCE #4")
print("CE emission = ",round(na_4),' s-1 g-1')
print("Activity = ",round(A_4)," Bq g-1")

print("MEAN")
print("Activity = ",round(np.mean([A_1,A_2,A_3,A_4])),"+/-",round(np.std([A_1,A_2,A_3,A_4]))," Bq g-1", round(100*np.std([A_1,A_2,A_3,A_4])/np.mean([A_1,A_2,A_3,A_4]),3)," %")


#factor_y=20

plt.figure("spectrum Tricarb"); plt.clf()

# TRICARB
plt.plot(E_1,S_1,'-r',label=r'Tricarb source #1')
plt.plot(E_1,fit_1_b,"--r",label='fit')
plt.plot(E_2,S_2,'-k',label=r'Tricarb source #2')
plt.plot(E_2,fit_2_b,"--k",label='fit')
plt.plot(E_3,S_3,'-g',label=r'Tricarb source #3')
plt.plot(E_3,fit_3_b,"--g",label='fit')
plt.plot(E_4,S_4,'-b',label=r'Tricarb source #4')
plt.plot(E_4,fit_4_b,"--b",label='fit')

plt.xlabel(r'$E$ /keV', fontsize=20)
plt.ylabel(r'$S$ /(s$^{-1}$ g$^{-1}$)', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()
plt.xlim([0,300])
plt.ylim([0,6000])
plt.show()


"DT5730S CALCULATION"

guess=np.array([sum(S_1_c1[80:350]),180,50,sum(S_1_c1[350:1500]),800,150])
fit_1=curve_fit(cdFct,E_1_c1[80:1500],S_1_c1[80:1500],guess)
fit_1_c=gauss(E_1_c1,fit_1[0][0],fit_1[0][1],fit_1[0][2])
fit_1_b=gauss(E_1_c1,fit_1[0][3],fit_1[0][4],fit_1[0][5])
na_1_1=sum(fit_1_b)

guess=np.array([sum(S_2_c1[80:350]),180,50,sum(S_2_c1[350:1500]),800,150])
fit_2=curve_fit(cdFct,E_2_c1[80:1500],S_2_c1[80:1500],guess)
fit_2_c=gauss(E_2_c1,fit_2[0][0],fit_2[0][1],fit_2[0][2])
fit_2_b=gauss(E_2_c1,fit_2[0][3],fit_2[0][4],fit_2[0][5])
na_2_1=sum(fit_2_b)

guess=np.array([sum(S_3_c1[80:350]),180,50,sum(S_3_c1[350:1500]),800,150])
fit_3=curve_fit(cdFct,E_3_c1[80:1500],S_3_c1[80:1500],guess)
fit_3_c=gauss(E_3_c1,fit_3[0][0],fit_3[0][1],fit_3[0][2])
fit_3_b=gauss(E_3_c1,fit_3[0][3],fit_3[0][4],fit_3[0][5])
na_3_1=sum(fit_3_b)

guess=np.array([sum(S_4_c1[80:350]),180,50,sum(S_4_c1[350:1500]),800,150])
fit_4=curve_fit(cdFct,E_4_c1[80:1500],S_4_c1[80:1500],guess)
fit_4_c=gauss(E_4_c1,fit_4[0][0],fit_4[0][1],fit_4[0][2])
fit_4_b=gauss(E_4_c1,fit_4[0][3],fit_4[0][4],fit_4[0][5])
na_4_1=sum(fit_4_b)

A_1_1=(na_1_1)*k_gamma/(epsilon*p_ce)
A_2_1=(na_2_1)*k_gamma/(epsilon*p_ce)
A_3_1=(na_3_1)*k_gamma/(epsilon*p_ce)
A_4_1=(na_4_1)*k_gamma/(epsilon*p_ce)

guess=np.array([sum(S_1_c2[80:350]),180,50,sum(S_1_c2[350:1500]),800,150])
fit_1=curve_fit(cdFct,E_1_c2[80:1500],S_1_c2[80:1500],guess)
fit_1_c=gauss(E_1_c2,fit_1[0][0],fit_1[0][1],fit_1[0][2])
fit_1_b=gauss(E_1_c2,fit_1[0][3],fit_1[0][4],fit_1[0][5])
na_1_2=sum(fit_1_b)

guess=np.array([sum(S_2_c2[80:350]),180,50,sum(S_2_c2[350:1500]),800,150])
fit_2=curve_fit(cdFct,E_2_c2[80:1500],S_2_c2[80:1500],guess)
fit_2_c=gauss(E_2_c2,fit_2[0][0],fit_2[0][1],fit_2[0][2])
fit_2_b=gauss(E_2_c2,fit_2[0][3],fit_2[0][4],fit_2[0][5])
na_2_2=sum(fit_2_b)

guess=np.array([sum(S_3_c2[80:350]),180,50,sum(S_3_c2[350:1500]),800,150])
fit_3=curve_fit(cdFct,E_3_c2[80:1500],S_3_c2[80:1500],guess)
fit_3_c=gauss(E_3_c2,fit_3[0][0],fit_3[0][1],fit_3[0][2])
fit_3_b=gauss(E_3_c2,fit_3[0][3],fit_3[0][4],fit_3[0][5])
na_3_2=sum(fit_3_b)

guess=np.array([sum(S_4_c2[80:350]),180,50,sum(S_4_c2[350:1500]),800,150])
fit_4=curve_fit(cdFct,E_4_c2[80:1500],S_4_c2[80:1500],guess)
fit_4_c=gauss(E_4_c2,fit_4[0][0],fit_4[0][1],fit_4[0][2])
fit_4_b=gauss(E_4_c2,fit_4[0][3],fit_4[0][4],fit_4[0][5])
na_4_2=sum(fit_4_b)

A_1_2=(na_1_2)*k_gamma/(epsilon*p_ce)
A_2_2=(na_2_2)*k_gamma/(epsilon*p_ce)
A_3_2=(na_3_2)*k_gamma/(epsilon*p_ce)
A_4_2=(na_4_2)*k_gamma/(epsilon*p_ce)

guess=np.array([sum(S_1_c3[80:350]),180,50,sum(S_1_c3[350:1500]),800,150])
fit_1=curve_fit(cdFct,E_1_c3[80:1500],S_1_c3[80:1500],guess)
fit_1_c=gauss(E_1_c3,fit_1[0][0],fit_1[0][1],fit_1[0][2])
fit_1_b=gauss(E_1_c3,fit_1[0][3],fit_1[0][4],fit_1[0][5])
na_1_3=sum(fit_1_b)

guess=np.array([sum(S_2_c3[80:350]),180,50,sum(S_2_c3[350:1500]),800,150])
fit_2=curve_fit(cdFct,E_2_c3[80:1500],S_2_c3[80:1500],guess)
fit_2_c=gauss(E_2_c3,fit_2[0][0],fit_2[0][1],fit_2[0][2])
fit_2_b=gauss(E_2_c3,fit_2[0][3],fit_2[0][4],fit_2[0][5])
na_2_3=sum(fit_2_b)

guess=np.array([sum(S_3_c3[80:350]),180,50,sum(S_3_c3[350:1500]),800,150])
fit_3=curve_fit(cdFct,E_3_c3[80:1500],S_3_c3[80:1500],guess)
fit_3_c=gauss(E_3_c3,fit_3[0][0],fit_3[0][1],fit_3[0][2])
fit_3_b=gauss(E_3_c3,fit_3[0][3],fit_3[0][4],fit_3[0][5])
na_3_3=sum(fit_3_b)

guess=np.array([sum(S_4_c3[80:350]),180,50,sum(S_4_c3[350:1500]),800,150])
fit_4=curve_fit(cdFct,E_4_c3[80:1500],S_4_c3[80:1500],guess)
fit_4_c=gauss(E_4_c3,fit_4[0][0],fit_4[0][1],fit_4[0][2])
fit_4_b=gauss(E_4_c3,fit_4[0][3],fit_4[0][4],fit_4[0][5])
na_4_3=sum(fit_4_b)

A_1_3=(na_1_3)*k_gamma/(epsilon*p_ce)
A_2_3=(na_2_3)*k_gamma/(epsilon*p_ce)
A_3_3=(na_3_3)*k_gamma/(epsilon*p_ce)
A_4_3=(na_4_3)*k_gamma/(epsilon*p_ce)


print("\n**** DT5730S ***")
print("SOURCE #1")
print("CE emission = ",round(na_1_1),round(na_1_2),round(na_1_3),' s-1 g-1')
#print("Hiden CE emission = ",round(nb_1)," s-1 g-1 (", round(100*nb_1/na_1,1),"%)")
print("Activity = ",round(A_1_1),round(A_1_2),round(A_1_3)," Bq g-1")

print("SOURCE #2")
print("CE emission = ",round(na_2_1),round(na_2_2),round(na_2_3),' s-1 g-1')
#print("Hiden CE emission = ",round(nb_2)," s-1 g-1 (", round(100*nb_2/na_2,1),"%)")
print("Activity = ",round(A_2_1),round(A_2_2),round(A_2_3)," Bq g-1")

print("SOURCE #3")
print("CE emission = ",round(na_3_1),round(na_3_2),round(na_3_3),' s-1 g-1')
#print("Hiden CE emission = ",round(nb_3)," s-1 g-1 (", round(100*nb_3/na_3,1),"%)")
print("Activity = ",round(A_3_1),round(A_3_2),round(A_3_3)," Bq g-1")

print("SOURCE #4")
print("CE emission = ",round(na_4_1),round(na_4_2),round(na_4_3),' s-1 g-1')
#print("Hiden CE emission = ",round(nb_4)," s-1 g-1 (", round(100*nb_4/na_4,1),"%)")
print("Activity = ",round(A_4_1),round(A_4_2),round(A_4_3)," Bq g-1")

print("MEAN")
print("Activity = ",round(np.mean([A_1_1,A_2_1,A_3_1,A_4_1])),"+/-",round(np.std([A_1_1,A_2_1,A_3_1,A_4_1]))," Bq g-1", round(100*np.std([A_1_1,A_2_1,A_3_1,A_4_1])/np.mean([A_1_1,A_2_1,A_3_1,A_4_1]),3)," %")




plt.figure("spectrum DT5730S"); plt.clf()
# DT5730S Spectrum
#plt.plot(E_1_c1[20:1500], S_1_c1[20:1500],'-r',label=r'DT5730S source #1 PM1')
#plt.plot(E_1_c1[80:1500],fit_1_b[80:1500],"--r",label='fit')
#plt.plot(E_2_c1[20:1500], S_2_c1[20:1500],'-k',label=r'DT5730S source #2 PM1')
#plt.plot(E_2_c1[80:1500],fit_2_b[80:1500],"--k",label='fit')
#plt.plot(E_3_c1[20:1500], S_3_c1[20:1500],'-g',label=r'DT5730S source #3 PM1')
#plt.plot(E_3_c1[80:1500],fit_3_b[80:1500],"--g",label='fit')
#plt.plot(E_4_c1[20:1500], S_4_c1[20:1500],'-b',label=r'DT5730S source #4 PM1')
#plt.plot(E_4_c1[80:1500],fit_4_b[80:1500],"--b",label='fit')

plt.plot(E_1_c3[20:1500], S_1_c3[20:1500],'-r',label=r'DT5730S source #1 PM1')
plt.plot(E_1_c3[80:1500],fit_1_b[80:1500],"--r",label='fit')
plt.plot(E_2_c3[20:1500], S_2_c3[20:1500],'-k',label=r'DT5730S source #2 PM1')
plt.plot(E_2_c3[80:1500],fit_2_b[80:1500],"--k",label='fit')
plt.plot(E_3_c3[20:1500], S_3_c3[20:1500],'-g',label=r'DT5730S source #3 PM1')
plt.plot(E_3_c3[80:1500],fit_3_b[80:1500],"--g",label='fit')
plt.plot(E_4_c3[20:1500], S_4_c3[20:1500],'-b',label=r'DT5730S source #4 PM1')
plt.plot(E_4_c3[80:1500],fit_4_b[80:1500],"--b",label='fit')


#plt.plot(E_1_c2, S_1_c2*factor_y,'-r',label=r'DT5730S source #1 PM2')
#plt.plot(E_2_c2, S_2_c2*factor_y,'-k',label=r'DT5730S source #2 PM2')
#plt.plot(E_3_c2, S_3_c2*factor_y,'-g',label=r'DT5730S source #3 PM2')
#plt.plot(E_4_c2, S_4_c2*factor_y,'-b',label=r'DT5730S source #4 PM2')
#         
#plt.plot(E_1_c3, S_1_c3*factor_y,'-k',label=r'DT5730S source #1 PM3')
#plt.plot(E_2_c3, S_2_c3*factor_y,'-g',label=r'DT5730S source #2 PM3')
#plt.plot(E_3_c3, S_3_c3*factor_y,'-r',label=r'DT5730S source #3 PM3')
#plt.plot(E_4_c3, S_4_c3*factor_y,'-b',label=r'DT5730S source #4 PM3')

#if listMode: plt.plot(E2_1,S2_1,'-r',label=r'DT5730S list mode')
#plt.plot([26,26],[0,4000],'--k',label=r'threshold')
plt.xlabel(r'$E$ /keV', fontsize=20)
plt.ylabel(r'$S$ /(s$^{-1}$ g$^{-1}$)', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()
#plt.xlim([0,1500])
#plt.ylim([0,6000])
plt.show()


"nanoTDCR CALCULATION"

guess=np.array([sum(S_1_td[200:1000]),500,100,sum(S_1_td[1000:1900]),1650,100,sum(S_1_td[1900:3300]),2250,100])
fit_1=curve_fit(cdFct2,E_1_td[200:3300],S_1_td[200:3300],guess)
fit_1_c=gauss(E_1_td,fit_1[0][0],fit_1[0][1],fit_1[0][2])
fit_1_b=gauss(E_1_td,fit_1[0][3],fit_1[0][4],fit_1[0][5])
fit_1_d=gauss(E_1_td,fit_1[0][6],fit_1[0][7],fit_1[0][8])
na_1=sum(fit_1_b+fit_1_d)
u_na_1=sum(S_1_td[3200:])+sum(fit_1_b[:1100])+sum(fit_1_d[:1100])

guess=np.array([sum(S_2_td[200:1000]),500,100,sum(S_2_td[1000:1900]),1650,100,sum(S_2_td[1900:3300]),2250,100])
fit_2=curve_fit(cdFct2,E_2_td[200:3300],S_2_td[200:3300],guess)
fit_2_c=gauss(E_2_td,fit_2[0][0],fit_2[0][1],fit_2[0][2])
fit_2_b=gauss(E_2_td,fit_2[0][3],fit_2[0][4],fit_2[0][5])
fit_2_d=gauss(E_2_td,fit_2[0][6],fit_2[0][7],fit_2[0][8])
na_2=sum(fit_2_b+fit_2_d)
u_na_2=sum(S_2_td[3200:])+sum(fit_2_b[:1100])+sum(fit_2_d[:1100])


guess=np.array([sum(S_3_td[200:1000]),500,100,sum(S_3_td[1000:1900]),1650,100,sum(S_3_td[1900:3300]),2250,100])
fit_3=curve_fit(cdFct2,E_3_td[200:3300],S_3_td[200:3300],guess)
fit_3_c=gauss(E_3_td,fit_3[0][0],fit_3[0][1],fit_3[0][2])
fit_3_b=gauss(E_3_td,fit_3[0][3],fit_3[0][4],fit_3[0][5])
fit_3_d=gauss(E_3_td,fit_3[0][6],fit_3[0][7],fit_3[0][8])
na_3=sum(fit_3_b+fit_3_d)
u_na_3=sum(S_3_td[3200:])+sum(fit_3_b[:1100])+sum(fit_3_d[:1100])

guess=np.array([sum(S_4_td[200:1000]),500,100,sum(S_4_td[1000:1900]),1650,100,sum(S_4_td[1900:3300]),2250,100])
fit_4=curve_fit(cdFct2,E_4_td[200:3300],S_4_td[200:3300],guess)
fit_4_c=gauss(E_4_td,fit_4[0][0],fit_4[0][1],fit_4[0][2])
fit_4_b=gauss(E_4_td,fit_4[0][3],fit_4[0][4],fit_4[0][5])
fit_4_d=gauss(E_4_td,fit_4[0][6],fit_4[0][7],fit_4[0][8])
na_4=sum(fit_4_b+fit_4_d)
u_na_4=sum(S_4_td[3200:])+sum(fit_4_b[:1100])+sum(fit_4_d[:1100])

A_1=(na_1)*k_gamma/(epsilon*p_ce)
A_2=(na_2)*k_gamma/(epsilon*p_ce)
A_3=(na_3)*k_gamma/(epsilon*p_ce)
A_4=(na_4)*k_gamma/(epsilon*p_ce)

print("\n**** nanoTDCR ***")
print("SOURCE #1")
print("CE emission = ",round(na_1),' s-1 g-1')
#print("Hiden CE emission = ",round(nb_1)," s-1 g-1 (", round(100*nb_1/na_1,1),"%)")
print("Activity = ",round(A_1)," Bq g-1")
#
print("SOURCE #2")
print("CE emission = ",round(na_2),' s-1 g-1')
#print("Hiden CE emission = ",round(nb_2)," s-1 g-1 (", round(100*nb_2/na_2,1),"%)")
print("Activity = ",round(A_2)," Bq g-1")

print("SOURCE #3")
print("CE emission = ",round(na_3),' s-1 g-1')
#print("Hiden CE emission = ",round(nb_3)," s-1 g-1 (", round(100*nb_3/na_3,1),"%)")
print("Activity = ",round(A_3)," Bq g-1")

print("SOURCE #4")
print("CE emission = ",round(na_4),' s-1 g-1')
#print("Hiden CE emission = ",round(nb_4)," s-1 g-1 (", round(100*nb_4/na_4,1),"%)")
print("Activity = ",round(A_4)," Bq g-1")


print("MEAN")
print("Activity = ",round(np.mean([A_1,A_2,A_3,A_4])),"+/-",round(np.std([A_1,A_2,A_3,A_4]))," Bq g-1", round(100*np.std([A_1,A_2,A_3,A_4])/np.mean([A_1,A_2,A_3,A_4]),3)," %")


print('ratio branching',fit_1[0][6]/fit_1[0][3],fit_2[0][6]/fit_2[0][3],fit_3[0][6]/fit_3[0][3],fit_4[0][6]/fit_4[0][3])


plt.figure("spectrum nanoTDCR"); plt.clf()

plt.plot(E_1_td, S_1_td,'-r',label=r'nanoTDCR source #1')
plt.plot(E_1_td[200:3300],fit_1_b[200:3300]+fit_1_d[200:3300],"--r",label='fit')
plt.plot(E_2_td, S_2_td,'-k',label=r'nanoTDCR source #2')
plt.plot(E_2_td[200:3300],fit_2_b[200:3300]+fit_2_d[200:3300],"--k",label='fit')
plt.plot(E_3_td, S_3_td,'-g',label=r'nanoTDCR source #3')
plt.plot(E_3_td[200:3300],fit_3_b[200:3300]+fit_3_d[200:3300],"--g",label='fit')
plt.plot(E_4_td, S_4_td,'-b',label=r'nanoTDCR source #4')
plt.plot(E_3_td[200:3300],fit_4_b[200:3300]+fit_4_d[200:3300],"--b",label='fit')

#if listMode: plt.plot(E2_1,S2_1,'-r',label=r'DT5730S list mode')
#plt.plot([26,26],[0,4000],'--k',label=r'threshold')
plt.xlabel(r'$E$ /keV', fontsize=20)
plt.ylabel(r'$S$ /(s$^{-1}$ g$^{-1}$)', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()
#plt.xlim([0,1500])
#plt.ylim([0,6000])
plt.show()


plt.figure("spectrum nanoTDCR - deconvol"); plt.clf()
e_cal=0.038
plt.plot(e_cal*E_1_td, S_1_td,'-k',label=r'nanoTDCR source #1')
plt.plot(e_cal*E_1_td[200:3300],fit_1_c[200:3300]+fit_1_b[200:3300]+fit_1_d[200:3300],'-r',label=r'fit total')
plt.plot(e_cal*E_1_td[200:3300],fit_1_b[200:3300]+fit_1_d[200:3300],"--g",label='fit ce')
#plt.plot(E_1_td[200:3300],fit_1_c[200:3300],"--b",label='fit ec')


#if listMode: plt.plot(E2_1,S2_1,'-r',label=r'DT5730S list mode')
#plt.plot([26,26],[0,4000],'--k',label=r'threshold')
plt.xlabel(r'$E$ /keV', fontsize=20)
plt.ylabel(r'$S$ /(s$^{-1}$ g$^{-1}$)', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()
#plt.xlim([0,1500])
#plt.ylim([0,6000])
plt.show()








# counting statistic
#ur_stat=100/np.sqrt(np.mean([na_1,na_2,na_3,na_4])*2000) # type B (Poisson hypothesis)
ur_stat=round(100*np.std([A_1,A_2,A_3,A_4])/np.mean([A_1,A_2,A_3,A_4]),3)

u_deconv=100*np.mean([u_na_1/na_1,u_na_2/na_2,u_na_3/na_3,u_na_4/na_4])


print("Uncertainty budget")
CountStat=ur_stat
print("Counting statistics = ",CountStat, " %")
print("background=", ur_bck, " %")
uLT=np.mean([0.6015,0.4926,0.3913,0.2957])
print("")
print("Live time = ",uLT," %")
print("Decay data =")
DecayData=100*u_p_ce/p_ce
print("\tP_ce prob =", DecayData, " %")
DecayCorr=100*uQ/Q
print("\tDecay corr =", DecayCorr, " %")
Weighing=100*5.5*4/35334
print("Weighing = ", Weighing, " %")
MC_corr=100*u_k_gamma/k_gamma
print("\tCorr gamma =", MC_corr, " %")
print("\tDeconv ec and ce event in spectrum",u_deconv," %")



u_combined=np.sqrt(CountStat**2+ur_bck**2+DecayData**2+uLT**2+DecayCorr**2+Weighing**2+MC_corr**2+u_deconv**2)
print("u combined", u_combined)
