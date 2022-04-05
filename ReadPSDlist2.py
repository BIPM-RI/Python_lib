# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 09:05:04 2021

@author: romain.coulon
"""

import matplotlib.pyplot as plt
import numpy as np
import csv
import seaborn as sb
import sys
sys.path.insert(1, 'G:\Python_modules\BIPM_RI_PyModules')
import DataProcessing as dp
import TDCRcalculation as tc

import sys, subprocess, configparser


print("*************************************")
print("****** PROCESS CAEN DT5730SB ********")
print("****** LIST MODE DATA        ********")
print("****** FOR TDCR MEASUREMENT  ********")
print("*************************************")

print("\n Reading of the configuration file...")
config = configparser.ConfigParser()
config.read('ConfigRunDigitizer.ini')

path=config["Inputs"].get("path")
file=config["Inputs"].get("file")
tdc_file=config["Inputs"].get("tdc_file")
lts_fileA=config["Inputs"].get("lts_fileA")
lts_fileB=config["Inputs"].get("lts_fileB")
lts_fileC=config["Inputs"].get("lts_fileC")

Bkg_mode=config["Options"].getboolean("Bkg_mode")
n_run=config["Options"].getint("n_run")

ResolTime=config["Parameters"].getfloat("ResolTime")
ExtDT=config["Parameters"].getfloat("ExtDT")
thres=config["Parameters"].getfloat("thres")
SumW=config["Parameters"].getint("SumW")
baseTime=config["Parameters"].getfloat("baseTime")
typeOfFilter=config["Parameters"].get("typeOfFilter")
scaling=config["Parameters"].getint("scaling")

print("\tInput file :", path+file)
print("\tOnput file :", path+tdc_file)
# BACKGROUND Expected values from previous TDCR measurements
#(Expected R_D=1.70(2) s-1)
#(Expected R_T=1.2(3) s-1)
#(Expected R_AB=1.37(1) s-1)
#(Expected R_BC=1.37(1) s-1)
#(Expected R_AC=1.36(1) s-1)
#(Expected R_A=350(80) s-1)
#(Expected R_B=360(60) s-1)
#(Expected R_C=320(70) s-1)
# PMTA SEP (Expected V=514(15) mV)
# PMTB SEP (Expected V=556(15) mV)
# PMTC SEP (Expected V=445(15) mV)
# PMTA VAL (Expected V=125 mV)
# PMTB VAL (Expected V=125 mV)
# PMTC VAL (Expected V=125 mV)

# HMI Expected values from previous TDCR measurements
#(Expected R_D=2280(1) s-1)
#(Expected TDCR=0.7389(4) s-1)


print("\nReading of the date of the measurement...")
f = open(path+"run.info",'r')
run_file=f.read() # read the file
f.close() # close the file
t_flag=run_file.find("time.start=")
date0=run_file[t_flag+11:t_flag+39]
date=date0[8:10]+"/"+date0[5:7]+"/"+date0[0:4]+"  "+date0[11:13]+":"+date0[14:16]+":"+date0[17:19]+"."+date0[20:23]
print("Date of the measurement: ", date)

print("\nReading of the duration of the measurement...") 
f = open(path+"run.info",'r')
run_file=f.read() # read the file
f.close() # close the file
t_flag=run_file.find("time.real=")
t_mes=float(run_file[t_flag+10:t_flag+12])*3600+float(run_file[t_flag+13:t_flag+15])*60+float(run_file[t_flag+16:t_flag+18])
t_mes/=n_run

t_flag=run_file.find("0.throughput=")
throughput=run_file[t_flag+13:t_flag+20]
t_flag=run_file.find("0.icr=")
icr=run_file[t_flag+6:t_flag+13]
t_flag=run_file.find("0.ocr=")
ocr=run_file[t_flag+6:t_flag+13]

print("\nInitialize the calculation...") 
# Convert resol and extended times
ResolTime=ResolTime*1e-9/baseTime # Coincidence resolution time in ns
ExtDT=ExtDT*1e-6/baseTime # Extented dead time in µs
#busy=(53+70+70%8)*1e-9/baseTime # see note about dead time
busy=0 # see note about dead time

hold=False
coincW=False
triplet=[False,False,False]
trigger_time=0

Acount=0
Bcount=0
Ccount=0
ABcount=0
BCcount=0
ACcount=0
Scount=0
Dcount=0
Tcount=0
comdtW=0
realTime=0
liveTime=0
count_pulse=0
t=0

E_list_A=[]
E_list_B=[]
E_list_C=[]

E_list_A_raw=[]
E_list_B_raw=[]
E_list_C_raw=[]


Arate=[]
Brate=[]
Crate=[]
ABrate=[]
BCrate=[]
ACrate=[]
Srate=[]
Drate=[]
Trate=[]
lt=[]

w=1

def resetTriplet():
    """
    This function aims to reset the triplet buffer recording events within the 3 channels.
    """
    return [False,False,False]


def writeTDCfile(f,date,w,t_mes,liveTime,Acount,Arate,Bcount,Brate,Ccount,Crate,ABcount,ABrate,BCcount,BCrate,ACcount,ACrate,Dcount,Drate,Tcount,Trate):
    """
    This function aims to generate a nanoTDCR like .tdc file as an output of the TDCR list mode processing.
    """
    f.write("Run"+str(w)+": START TIME: UTC,"+date+"\n")  #TBE
    f.write("Run"+str(w)+": Real Time [s],   "+str(round(t_mes,3))+"\n")
    f.write("Run"+str(w)+": Live Time PMT A E1 [s],   "+str(round(liveTime,3))+"\n") 
    f.write("Run"+str(w)+": Live Time PMT A E2 [s],   "+str(round(liveTime,3))+"\n")
    f.write("Run"+str(w)+": Live Time PMT B E1 [s],   "+str(round(liveTime,3))+"\n")
    f.write("Run"+str(w)+": Live Time PMT B E2 [s],   "+str(round(liveTime,3))+"\n")
    f.write("Run"+str(w)+": Live Time PMT C E1 [s],   "+str(round(liveTime,3))+"\n")
    f.write("Run"+str(w)+": Live Time PMT C E2 [s],   "+str(round(liveTime,3))+"\n")
    f.write("\n")
    f.write("Run"+str(w)+": PMT A Raw [counts], "+str(int(Acount))+"\n")
    f.write("Run"+str(w)+": PMT A Raw [cps],  "+str(round(Arate[-1],3))+"\n")
    f.write("Run"+str(w)+": PMT B Raw [counts], "+str(int(Bcount))+"\n")
    f.write("Run"+str(w)+": PMT B Raw [cps],  "+str(round(Brate[-1],3))+"\n")
    f.write("Run"+str(w)+": PMT C Raw [counts], "+str(int(Ccount))+"\n")
    f.write("Run"+str(w)+": PMT C Raw [cps],  "+str(round(Crate[-1],3))+"\n")
    f.write("Run"+str(w)+": PMT A Ext1 [counts], "+str(int(Acount))+"\n")
    f.write("Run"+str(w)+": PMT A Ext1 [cps],  "+str(round(Arate[-1],3))+"\n")
    f.write("Run"+str(w)+": PMT B Ext1 [counts], "+str(int(Bcount))+"\n")
    f.write("Run"+str(w)+": PMT B Ext1 [cps],  "+str(round(Brate[-1],3))+"\n")
    f.write("Run"+str(w)+": PMT C Ext1 [counts], "+str(int(Ccount))+"\n")
    f.write("Run"+str(w)+": PMT C Ext1 [cps],  "+str(round(Crate[-1],3))+"\n")
    f.write("Run"+str(w)+": PMT A Ext2 [counts], "+str(int(Acount))+"\n")
    f.write("Run"+str(w)+": PMT A Ext2 [cps],  "+str(round(Arate[-1],3))+"\n")
    f.write("Run"+str(w)+": PMT B Ext2 [counts], "+str(int(Bcount))+"\n")
    f.write("Run"+str(w)+": PMT B Ext2 [cps],  "+str(round(Brate[-1],3))+"\n")
    f.write("Run"+str(w)+": PMT C Ext2 [counts], "+str(int(Ccount))+"\n")
    f.write("Run"+str(w)+": PMT C Ext2 [cps],  "+str(round(Crate[-1],3))+"\n")
    f.write("\n")
    f.write("Run"+str(w)+": Live Time AB Ext1 [s],   "+str(round(liveTime,3))+"\n")
    f.write("Run"+str(w)+": Live Time BC Ext1 [s],   "+str(round(liveTime,3))+"\n")
    f.write("Run"+str(w)+": Live Time AC Ext1 [s],   "+str(round(liveTime,3))+"\n")
    f.write("Run"+str(w)+": Live Time T Ext1 [s],   "+str(round(liveTime,3))+"\n")
    f.write("\n")
    f.write("Run"+str(w)+": Coincidence AB_N1 [counts], "+str(int(ABcount))+"\n")
    f.write("Run"+str(w)+": Coincidence AB_N1 [cps],  "+str(round(ABrate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence BC_N1 [counts], "+str(int(BCcount))+"\n")
    f.write("Run"+str(w)+": Coincidence BC_N1 [cps],  "+str(round(BCrate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence AC_N1 [counts], "+str(int(ACcount))+"\n")
    f.write("Run"+str(w)+": Coincidence AC_N1 [cps],  "+str(round(ACrate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence D_N1 [counts], "+str(int(Dcount))+"\n")
    f.write("Run"+str(w)+": Coincidence D_N1 [cps],  "+str(round(Drate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence T_N1 [counts], "+str(int(Tcount))+"\n")
    f.write("Run"+str(w)+": Coincidence T_N1 [cps],   "+str(round(Trate[-1],3))+"\n")
    f.write("\n")
    f.write("Run"+str(w)+": Coincidence AB_M1 [counts], "+str(int(ABcount))+"\n")
    f.write("Run"+str(w)+": Coincidence AB_M1 [cps],  "+str(round(ABrate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence BC_M1 [counts], "+str(int(BCcount))+"\n")
    f.write("Run"+str(w)+": Coincidence BC_M1 [cps],  "+str(round(BCrate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence AC_M1 [counts], "+str(int(ACcount))+"\n")
    f.write("Run"+str(w)+": Coincidence AC_M1 [cps],  "+str(round(ACrate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence D_M1 [counts], "+str(int(Dcount))+"\n")
    f.write("Run"+str(w)+": Coincidence D_M1 [cps],  "+str(round(Drate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence T_M1 [counts], "+str(int(Tcount))+"\n")
    f.write("Run"+str(w)+": Coincidence T_M1 [cps],   "+str(round(Trate[-1],3))+"\n")
    f.write("\n")
    f.write("Run"+str(w)+": Live Time AB Ext2 [s],   "+str(round(liveTime,3))+"\n")
    f.write("Run"+str(w)+": Live Time BC Ext2 [s],   "+str(round(liveTime,3))+"\n")
    f.write("Run"+str(w)+": Live Time AC Ext2 [s],   "+str(round(liveTime,3))+"\n")
    f.write("Run"+str(w)+": Live Time T Ext2 [s],   "+str(round(liveTime,3))+"\n")
    f.write("\n")
    f.write("Run"+str(w)+": Coincidence AB_N2 [counts], "+str(int(ABcount))+"\n")
    f.write("Run"+str(w)+": Coincidence AB_N2 [cps],  "+str(round(ABrate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence BC_N2 [counts], "+str(int(BCcount))+"\n")
    f.write("Run"+str(w)+": Coincidence BC_N2 [cps],  "+str(round(BCrate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence AC_N2 [counts], "+str(int(ACcount))+"\n")
    f.write("Run"+str(w)+": Coincidence AC_N2 [cps],  "+str(round(ACrate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence D_N2 [counts], "+str(int(Dcount))+"\n")
    f.write("Run"+str(w)+": Coincidence D_N2 [cps],  "+str(round(Drate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence T_N2 [counts], "+str(int(Tcount))+"\n")
    f.write("Run"+str(w)+": Coincidence T_N2 [cps],   "+str(round(Trate[-1],3))+"\n")
    f.write("\n")
    f.write("Run"+str(w)+": Coincidence AB_M2 [counts], "+str(int(ABcount))+"\n")
    f.write("Run"+str(w)+": Coincidence AB_M2 [cps],  "+str(round(ABrate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence BC_M2 [counts], "+str(int(BCcount))+"\n")
    f.write("Run"+str(w)+": Coincidence BC_M2 [cps],  "+str(round(BCrate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence AC_M2 [counts], "+str(int(ACcount))+"\n")
    f.write("Run"+str(w)+": Coincidence AC_M2 [cps],  "+str(round(ACrate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence D_M2 [counts], "+str(int(Dcount))+"\n")
    f.write("Run"+str(w)+": Coincidence D_M2 [cps],  "+str(round(Drate[-1],3))+"\n")
    f.write("Run"+str(w)+": Coincidence T_M2 [counts], "+str(int(Tcount))+"\n")
    f.write("Run"+str(w)+": Coincidence T_M2 [cps],   "+str(round(Trate[-1],3))+"\n")
    f.write("\n")
    
    
    
    





print("\nOpen the .tdc file for writting...") 
f = open(path+typeOfFilter+tdc_file,'w')

f.write("Hardware, CAEN DT5730S\n")
f.write("Serial Number, 2-11-14158\n")  # TBE
f.write("Firmware Version, 4.23 build 4B06\n") # TBE
f.write("FPGA Version, Arria V GX\n") # TBE
f.write("Measurement, Settings\n") # TBE
f.write("Measurement Tag,TDCR Measurement\n") # TBE
f.write("Operator,\n") # TBE
f.write("Nuclides,\n") # TBE
f.write("LS Coctail,\n") # TBE
f.write("High Voltage,\n") # TBE
f.write("Threshold A [mV],-"+str(round(thres,2))+"\n")
f.write("Threshold B [mV],-"+str(round(thres,2))+"\n")
f.write("Threshold C [mV],-"+str(round(thres,2))+"\n")
f.write("Dead Time Extension 1 [us], "+str(round(ExtDT*baseTime*1e6,3))+"\n")    
f.write("Dead Time Extension 2 [us], "+str(round(ExtDT*baseTime*1e6,3))+"\n")
f.write("Coincidence Window N [ns],"+str(int(ResolTime*baseTime*1e9))+"\n")   
f.write("Coincidence Window M [ns],"+str(int(ResolTime*baseTime*1e9))+"\n")
f.write("Preset Sequential Runs,"+str(n_run)+"\n")
f.write("Preset Time per Run [s], "+str(round(t_mes,3))+"\n")   
f.write("Time Gap between Runs [s],  0.000\n") 
f.write("Stop Timer,REAL\n")
f.write(" \n")
f.write("Measurement, Data\n")
f.write("Runs Completed, "+str(n_run)+" \n")
f.write(" \n")

print("\n Process the TDCR signal processing...") 
with open(path+typeOfFilter+file) as csv_file:
    csv = csv.reader(csv_file, delimiter=';')
    for row in csv:
        if row[0] != "BOARD": # skip the header
            t=int(row[2])            # arrival time of the pulse event
            realTime=t*baseTime      # arrival time of the pulse event in second
            charge=float(row[SumW])  # charge of the pulse event

            if Bkg_mode: # record every event to build energy spectra
                if int(row[1])==0: E_list_A_raw.append(charge) # record charge to build channel A spectrum 
                if int(row[1])==1: E_list_B_raw.append(charge) # record charge to build channel B spectrum 
                if int(row[1])==2: E_list_C_raw.append(charge) # record charge to build channel C spectrum
            
            # (1) DETECT EVENT WHEN IDLE
            # ie. when the charge is above the energy threshold and when the
            # is set available after the appropriate extended time period (comdtW)
            # since the previous triggering.
            if charge>=thres and t>trigger_time+comdtW:
                # (1.1) RECORD THE PREVIOUS EVENT IN COUNTERS AND UPDATE THE LIVETIME
                # Read the triplet buffer
                if sum(triplet)>=1: # at least one event recorded among the channels 
                    Scount+=1 # record the single events in the counter S
                    if triplet==[True, False, False]: Acount+=1 # record event in channel A only
                    if triplet==[False, True, False]: Bcount+=1 # record event in channel B only
                    if triplet==[False, False, True]: Ccount+=1 # record event in channel C only
                    if sum(triplet)>=2: # at least two events recorded among the channels
                        Dcount+=1 # record the double coincidence event in the counter D
                        if triplet==[True, True, False]: ABcount+=1 # record coincident event in channel A and B
                        if triplet==[False, True, True]: BCcount+=1 # record coincident event in channel B and C
                        if triplet==[True, False, True]: ACcount+=1 # record coincident event in channel A and C
                        if sum(triplet)==3: # at least three events recorded in each od the channels
                            Tcount+=1; ABcount+=1; BCcount+=1; ACcount+=1 # record the triple coincidence event in the counter T
                else:
                    print("Warning: no count recorded when triggering in idle")
                
                # (1.1.(opt))PRODUCE AN INTERMEDIATE MEASUREMENT RESULT                                            
                if realTime>t_mes*w: # end of a run
                    print("\n****")
                    print(int(realTime*w/60)," min processed.")
                    print("run #",w)
                    print('****')
#                    print("A,B,C ",round(Acount/liveTime,1),round(Bcount/liveTime,1),round(Ccount/liveTime,1))
#                    print("AB, BC, AC", round(ABcount/liveTime,1),round(BCcount/liveTime,1),round(ACcount/liveTime,1))
#                    print("S, D, T", round(Scount/liveTime,1), round(Dcount/liveTime,1), round(Tcount/liveTime,1))
#                    print("Check sum = ",ABcount+BCcount+ACcount-2*Tcount==Dcount)
#                    print("TDCR = ",round(Tcount/(Dcount+1),3))
#                    print("Dead time = ", round(100*(1-liveTime/realTime)), " %")
#                    print("LiveTime = ", round(liveTime,2), " s")
#                    print("Real Time = ", round(realTime,2), " s")
                    
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
                    lt.append(liveTime)
                    # WRITE VALUES IN THE REPORT
                    writeTDCfile(f,date,w,t_mes,liveTime,Acount,Arate,Bcount,Brate,Ccount,Crate,ABcount,ABrate,BCcount,BCrate,ACcount,ACrate,Dcount,Drate,Tcount,Trate)
                    
                    # REINITIALIZE THE COUNTERS
                    Acount=0; Bcount=0; Ccount=0; ABcount=0; BCcount=0; ACcount=0; Scount=0; Dcount=0; Tcount=0
                    liveTime=0
                    w+=1 # new measurement
                
                # (1.2) UPDATE THE LIVETIME CALCULATION
                # add the real time interval (t-trigger_time)
                # and substract the extended dead time (comdtW)
                # and substract the intrinsic dead time (count_pulse*busy)
                liveTime+=(t-(trigger_time+comdtW+count_pulse*busy))*baseTime # livetime in s
                
                
                # (1.3) REINITIALIZE BUFFERS FOR A NEW EVENT ANALYSIS
                triplet=resetTriplet()      # reset the coincidence event buffer
                trigger_time=t              # set the new triggering time
                triplet[int(row[1])]=True   # record of the event in the corresponding channel in the triplet buffer
                comdtW=ExtDT                # impose the extended dead time
                count_pulse=1               # initialize the pulse counter used to evaluate the intrinsic dead time
            
            # (2) DETECT EVENT WHEN BUSY
            # Events arriving during the coincidence resolving time are recorded in the triplet event buffer
            # Cumulate the parralytime time by adding an extended dead time period starting from the pulse arrival
            elif charge>=thres: # only pulse event with charge above the defined threshold
                count_pulse+=1 # increment the pulse counter
                # (2-1) DETECT EVENT IN THE COINCIDENCE WINDOW
                if t<trigger_time+ResolTime: # until the arrival is within the resolving time since the primary trigger
                    triplet[int(row[1])]=True # record of the event in the corresponding channel in the triplet buffer
                    comdtW+=t-trigger_time # reconduct the parallyzing time by the period from the primary trigger
                # (2-1) DETECT EVENT OUT OF THE COINCIDENCE WINDOW
                else:
                    # print((t-trigger_time)*baseTime*1e6)
                    comdtW=t-trigger_time+ExtDT # reconduct the parrallyzing time by adding an extended dead time period from the hidden event arrival time  
            else:
#                comdtW=t-trigger_time+ExtDT
                count_pulse+=1 # increment the pulse counter
#                print("Warnng: detect event below the thresold", charge," < ",thres) 
        else:
            print("\tHeader: ",row) # print the header of the list mode file 
         
            
                
        
# Write LabZy-like output file
f.close()
print("\nEnd of the TDCR signal processing and recording of the .tdc file...") 


tdcr=np.mean(np.asarray(Trate)/np.asarray(Drate))
ab=np.mean(ABrate)
bc=np.mean(BCrate)
ac=np.mean(ACrate)
t=np.mean(Trate)
d=np.mean(Drate)

print("\n****")
print("RESULTS")
print("****")
print("thres",thres)
print("A,B,C ",round(np.mean(Arate),1),round(np.mean(Brate),1),round(np.mean(Crate),1))
print("u(A),u(B),u(C) ",round(np.std(Arate),1),round(np.std(Brate),1),round(np.std(Crate),1))
print("AB, BC, AC", round(ab,1),round(bc,1),round(ac,1))
print("u(AB), u(BC), u(AC)", round(np.std(ABrate),1),round(np.std(BCrate),1),round(np.std(ACrate),1))
print("S, D, T", round(np.mean(Srate),1),round(d,1),round(t,1))
print("u(S), u(D), u(T)", round(np.std(Srate),1),round(np.std(Drate),1),round(np.std(Trate),1))
print("Check sum = ",ABcount+BCcount+ACcount-2*Tcount==Dcount)
print("TDCR = ",round(tdcr,4))
print("u(TDCR) = ",round(np.std(np.asarray(Trate)/np.asarray(Drate)),4))
print("Dead time = ", round(100-100*sum(lt)/realTime,1), " %")
print("LiveTime = ", round(sum(lt),2), " s")
print("Real Time = ", round(realTime,2), " s")


#out_i2=tc.I2calc(tdcr,t/ab,t/bc,t/ac,"H-3",0.0126461904139022) #HMI♣
#out_i2=tc.I2calc(tdcr,t/ab,t/bc,t/ac,"C-14",0.007546825) #CMF
out_i2=tc.I2calc(tdcr,t/ab,t/bc,t/ac,"C-14",0.006496422) #CCT

print("I2 = ",d/out_i2[2],d/out_i2[3])


plt.figure("Single count rates")
plt.clf()
plt.plot(Arate,'-r',label="channel A")
plt.plot(Brate,'-k',label="channel B")
plt.plot(Crate,'-b',label="channel C")
plt.xlabel("run #")
plt.ylabel(r"count rate /s$^{-1}$")
plt.legend()
plt.show()

plt.figure("Double coincidence count rates")
plt.clf()
plt.plot(Drate,'-r',label="channel D")
plt.plot(ABrate,'-g',label="channel AB")
plt.plot(BCrate,'-k',label="channel BC")
plt.plot(ACrate,'-b',label="channel AC")
plt.xlabel("run #")
plt.ylabel(r"count rate /s$^{-1}$")
plt.legend()
plt.show()

plt.figure("Triplet coincidence count rates")
plt.clf()
plt.plot(Trate,'-g',label="channel T")
plt.xlabel("run #")
plt.ylabel(r"count rate /s$^{-1}$")
plt.legend()
plt.show()



if Bkg_mode:
    print("process spectrometric analysis...")
    # Display Spectrum
    # sb.set_style("whitegrid")  # Setting style(Optional)
    # plt.figure(figsize = (10,5)) #Specify the size of figure we want(Optional)
    # sb.distplot(x = E_list_A  ,  bins = 4000 , kde = False , color = 'teal'\
    #              , kde_kws=dict(linewidth = 4 , color = 'black'))
    # [h.get_height() for h in sb.distplot(E_list_A).patches]
    # plt.show()

    
    plt.figure("Spectra raw")
    plt.clf()
    sA_r=plt.hist(np.asarray(E_list_A_raw),bins=scaling,label="channel A raw")[0]
    sB_r=plt.hist(np.asarray(E_list_B_raw),bins=scaling,label="chanel B raw")[0]
    sC_r=plt.hist(np.asarray(E_list_C_raw),bins=scaling,label="channel C raw")[0]
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
    
    sA_r[-1]=0
    sA_r=sA_r/sum(sA_r)
    sA_r=dp.MAfilter(sA_r,7)
    sB_r[-1]=0
    sB_r=sB_r/sum(sB_r)
    sB_r=dp.MAfilter(sB_r,7)
    sC_r[-1]=0
    sC_r=sC_r/sum(sC_r)
    sC_r=dp.MAfilter(sC_r,7)
                
    MaxSA_r=np.argmax(sA_r[250:])+250
    MaxSB_r=np.argmax(sB_r[250:])+250
    MaxSC_r=np.argmax(sC_r[250:])+250
    
    ValSA_r=np.argmin(sA_r[50:int(MaxSA_r)])+50
    ValSB_r=np.argmin(sB_r[50:int(MaxSB_r)])+50
    ValSC_r=np.argmin(sC_r[50:int(MaxSC_r)])+50
    
    
    
    plt.figure("Spectra raw 2")
    plt.clf()
    plt.plot(sA_r,label="channel A raw")
    plt.plot(sB_r,label="channel B raw")
    plt.plot(sC_r,label="channel C raw")
    plt.legend()
    plt.show()
    
        
    
    plt.figure("Spectra coinc")
    plt.clf()
    sA=plt.hist(np.asarray(E_list_A),bins=scaling,label="channel A")[0]
    sB=plt.hist(np.asarray(E_list_B),bins=scaling,label="chanel B")[0]
    sC=plt.hist(np.asarray(E_list_C),bins=scaling,label="channel C")[0]
    plt.legend()
    plt.show()
    
    
    sA=filterHF(sA)
    sB=filterHF(sB)
    sC=filterHF(sC)
    
    sA[-1]=0
    sA=sA/sum(sA)
    sA=dp.MAfilter(sA,15)
    sB[-1]=0
    sB=sB/sum(sB)
    sB=dp.MAfilter(sB,15)
    sC[-1]=0
    sC=sC/sum(sC)
    sC=dp.MAfilter(sC,15)
                
    MaxSA=np.argmax(sA[250:])+250
    MaxSB=np.argmax(sB[250:])+250
    MaxSC=np.argmax(sC[250:])+250
    
    ValSA=np.argmin(sA[50:int(MaxSA)])+50
    ValSB=np.argmin(sB[50:int(MaxSB)])+50
    ValSC=np.argmin(sC[50:int(MaxSC)])+50
    
    plt.figure("Spectra coinc 2")
    plt.clf()
    plt.plot(sA,label="channel A raw")
    plt.plot(sB,label="channel B raw")
    plt.plot(sC,label="channel C raw")
    plt.legend()
    plt.show()
    
    print("Write .lts files...")
    fA = open(path+typeOfFilter+lts_fileA,'w')
    fA.write("<?xml version=\"1.0\"?>\n")
    fA.write("<caen DT5730>\n")
    fA.write(" <serialnumber>10030</serialnumber>\n")
    fA.write("<spectrum>\n")
    fA.write("  <tag>labZY-TDCR spectrum</tag>\n")
    fA.write("  <hardsize>4096</hardsize>\n")
    fA.write("  <softsize>4096</softsize>\n")
    fA.write("  <offset>4096</offset>\n")
    fA.write("  <data>\n")
    for i in sA_r:    fA.write(str(i)+"\n")
    fA.write("  </data>\n")
    fA.write("</registers>\n")
    fA.write("<calibration>\n")
    fA.write("  <enabled>NO</enabled>\n")
    fA.write("  <units>0</units>\n")
    fA.write("  <useall>NO</useall>\n")
    fA.write("</calibration>\n")
    fA.write("<volatile>\n")
    fA.write("  <firmware>50.20</firmware>\n")
    fA.write("  <intemp>20</intemp>\n")
    fA.write("  <slowadc> 0.82</slowadc>\n")
    fA.write("</volatile>\n")
    fA.write("</caen DT5730>\n")
    fA.close()
    
    fB = open(path+typeOfFilter+lts_fileB,'w')
    fB.write("<?xml version=\"1.0\"?>\n")
    fB.write("<caen DT5730>\n")
    fB.write(" <serialnumber>10030</serialnumber>\n")
    fB.write("<spectrum>\n")
    fB.write("  <tag>labZY-TDCR spectrum</tag>\n")
    fB.write("  <hardsize>4096</hardsize>\n")
    fB.write("  <softsize>4096</softsize>\n")
    fB.write("  <offset>4096</offset>\n")
    fB.write("  <data>\n")
    for i in sB_r:    fB.write(str(i)+"\n")
    fB.write("  </data>\n")
    fB.write("</registers>\n")
    fB.write("<calibration>\n")
    fB.write("  <enabled>NO</enabled>\n")
    fB.write("  <units>0</units>\n")
    fB.write("  <useall>NO</useall>\n")
    fB.write("</calibration>\n")
    fB.write("<volatile>\n")
    fB.write("  <firmware>50.20</firmware>\n")
    fB.write("  <intemp>20</intemp>\n")
    fB.write("  <slowadc> 0.82</slowadc>\n")
    fB.write("</volatile>\n")
    fB.write("</caen DT5730>\n")
    fB.close()
    
    fC = open(path+typeOfFilter+lts_fileC,'w')
    fC.write("<?xml version=\"1.0\"?>\n")
    fC.write("<caen DT5730>\n")
    fC.write(" <serialnumber>10030</serialnumber>\n")
    fC.write("<spectrum>\n")
    fC.write("  <tag>labZY-TDCR spectrum</tag>\n")
    fC.write("  <hardsize>4096</hardsize>\n")
    fC.write("  <softsize>4096</softsize>\n")
    fC.write("  <offset>4096</offset>\n")
    fC.write("  <data>\n")
    for i in sC_r:    fC.write(str(i)+"\n")
    fC.write("  </data>\n")
    fC.write("</registers>\n")
    fC.write("<calibration>\n")
    fC.write("  <enabled>NO</enabled>\n")
    fC.write("  <units>0</units>\n")
    fC.write("  <useall>NO</useall>\n")
    fC.write("</calibration>\n")
    fC.write("<volatile>\n")
    fC.write("  <firmware>50.20</firmware>\n")
    fC.write("  <intemp>20</intemp>\n")
    fC.write("  <slowadc> 0.82</slowadc>\n")
    fC.write("</volatile>\n")
    fC.write("</caen DT5730>\n")
    fC.close()
    
    print("****")
    print("Spectra analysis")
    print("****")
    print("All pulses")
    print("SEP centroid A,B,C = ", MaxSA_r, MaxSB_r, MaxSC_r)
    print("SEP valley A,B,C = ", ValSA_r, ValSB_r, ValSC_r)
    print("coincident events")
    print("SEP centroid A,B,C = ", MaxSA, MaxSB, MaxSC)
    print("SEP valley A,B,C = ", ValSA, ValSB, ValSC)

print("\nEND.")