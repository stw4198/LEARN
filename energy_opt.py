from cmath import tan
import subprocess
import os
import pandas as pd
from scipy import integrate
import numpy as np
import sys
from sig_choice import sig_choice,sig_dict
from tqdm import tqdm

path = sys.argv[1]
tank = int(sys.argv[2])
sig = sys.argv[3]
file = path + "/results_ML.csv"

signal_components,background_components = sig_choice(sig)
components = signal_components+background_components
muon_eff = 0.95
t_veto = 1
t_half_n17 = 4.173
t_half_li9 = 0.1783
def integrand_n17(t):
    return np.exp(-np.log(2)*t/t_half_n17)
def integrand_li9(t):
    return np.exp(-np.log(2)*t/t_half_li9)
R_n17_cor = 1-muon_eff+muon_eff*integrate.quad(integrand_n17,t_veto,np.infty)[0]/integrate.quad(integrand_n17,0,np.infty)[0]
R_li9_cor = 1-muon_eff+muon_eff*integrate.quad(integrand_li9,t_veto,np.infty)[0]/integrate.quad(integrand_li9,0,np.infty)[0]

E_lower = np.arange(21,22,1)
E_upper = np.arange(66,67,1)
dwell_times = []
s_total = []
b_total = []
b_err_total = []
b_str_total = []

progress = 0
df_list = []

for x in tqdm(E_lower, desc='Total'):
    df_list_1 = []

    dwell_times_upper = []
    s_upper = []
    b_upper = []
    b_err_upper = []
    b_str_upper = []
    for y in tqdm(E_upper,desc='Inner Loop',leave=False):
        df = pd.read_csv(file)
        df = df.set_index('Component')
        sig_output_removed = []
        sig_output_total = []

        for i in components:
            output = subprocess.run([os.path.dirname(__file__)+"/energy",path+"/"+"ML_classified.root",str(x),str(y),str(sig_dict[i])],\
                stdout=subprocess.PIPE,universal_newlines = True).stdout
            sig_output_removed.append(int(output.split()[1]))
            sig_output_total.append(int(output.split()[4]))
            df.loc[i,str(tank)+" kept"]=int(output.split()[4])-int(output.split()[1])

        rates = []
        fn_cl = 3.69
        for i in components:
            if i == 'fn':
                if df.loc[i,"%i kept"%tank] == 0:
                    df.loc[i,"%i kept"%tank]=fn_cl
                rates.append(df.loc[i,"%i kept"%tank]*df.loc[i,'%i rate'%tank]*86400/df.loc[i,'%i MC'%tank])
            elif i == 'li9':
                df.loc[i,"%i kept"%tank]=df.loc[i,"%i kept"%tank]*R_li9_cor
                rates.append(df.loc[i,"%i kept"%tank]*df.loc[i,'%i rate'%tank]*86400/df.loc[i,'%i MC'%tank])
            elif i == 'n17':
                df.loc[i,"%i kept"%tank]=df.loc[i,"%i kept"%tank]*R_n17_cor
                rates.append(df.loc[i,"%i kept"%tank]*df.loc[i,'%i rate'%tank]*86400/df.loc[i,'%i MC'%tank])
            else:
                if df.loc[i,"%i kept"%tank]==0 or df.loc[i,'%i rate'%tank]==0 or df.loc[i,'%i MC'%tank]==0:
                    rates.append(0)
                else:
                    rates.append(df.loc[i,"%i kept"%tank]*df.loc[i,'%i rate'%tank]*86400/df.loc[i,'%i MC'%tank])
        df_list_1.append(df)
        d = {'components':components,'rates':rates}
        df_analysis = pd.DataFrame(data=d)
        df_analysis = df_analysis.set_index('components')
        s = 0
        b = 0
        f = 0
        ferr = 0.27
        li9 = 0
        li9err = 0.002
        n17 = 0
        n17err = 0.002
        world = 0
        hinkley_c = 0
        sizewell_b = 0
        gravelines = 0
        world_err = 0.06
        geo = 0
        geoerr = 0.25
        b_str=""
        for i in signal_components:
            s+=df_analysis.loc[i,'rates']
            b_str+="%s: %s per day\n"%(i,df_analysis.loc[i,'rates'])
        for k in background_components:
            if k=="world":
                b+=df_analysis.loc[k,'rates']
                world+=df_analysis.loc[k,'rates']
            elif k=="hinkley_c":
                b+=df_analysis.loc[k,'rates']
                hinkley_c+=df_analysis.loc[k,'rates']
            elif k=="gravelines":
                b+=df_analysis.loc[k,'rates']
                gravelines+=df_analysis.loc[k,'rates']
            elif k=="sizewell_b":
                b+=df_analysis.loc[k,'rates']
                sizewell_b+=df_analysis.loc[k,'rates']
            elif k=="geo":
                b+=df_analysis.loc[k,'rates']
                geo+=df_analysis.loc[k,'rates']
            elif k=="fn":
                b+=df_analysis.loc[k,'rates']#/2
                f+=df_analysis.loc[k,'rates']#/2
                #b+=0.04645677948717948*0.254
                #f+=0.04645677948717948*0.254
            elif k=="li9":
                b+=df_analysis.loc[k,'rates']
                li9+=df_analysis.loc[k,'rates']
            elif k=="n17":
                b+=df_analysis.loc[k,'rates']
                n17+=df_analysis.loc[k,'rates']
            else:
                b+=df_analysis.loc[k,'rates']
            b_str+="%s: %s per day\n"%(k,df_analysis.loc[k,'rates'])
        b_str_upper.append(b_str)
        b_err = np.sqrt((li9err*li9)**2 + (n17err*n17)**2 + (world_err*world)**2 + (geoerr*geo)**2 + (ferr*f)**2 + (world_err*hinkley_c)**2 + (world_err*gravelines)**2 + (world_err*sizewell_b)**2)
        #if (sig=='hartlepool') or (sig=='hartlepool_1' and tank==22):
        if (sig=='hartlepool' or sig=='hartlepool_1'):
            t = np.arange(0.01,30,0.01)
            sigma = np.sqrt(2*((s*t + b*t)*np.log((s*t + b*t)*(b*t + b_err**2*t**2)/(b**2*t**2 + (s*t + b*t)*b_err**2*t**2))\
              - (b**2*t**2)*np.log(1 + (b_err**2*t**2*s*t)/(b*t*(b*t + b_err**2*t**2)))/(b_err**2*t**2)))
            if max(sigma)>3:
                time_index = [ n for n,i in enumerate(sigma) if i>3][0]
                t3sigma = t[time_index]
            else:
                t3sigma = 100000
        else:
            #t3sigma = 9*(b+s)/(s**2 - 9*(world_err*s)**2 - 9*(b_err)**2)
            t3sigma = 9*b/(s**2 - 9*(b_err)**2)   
        s_upper.append(s)
        b_upper.append(b)
        b_err_upper.append(b_err)
        if t3sigma < 0:
            dwell_times_upper.append(100000)
        else:
            dwell_times_upper.append(t3sigma)
    dwell_times.append(dwell_times_upper)
    s_total.append(s_upper)
    b_total.append(b_upper)
    b_err_total.append(b_err_upper)
    progress += 1
    df_list.append(df_list_1)
    b_str_total.append(b_str_upper)

dwell_times = np.reshape(dwell_times,(len(E_lower),len(E_upper)))
s_total = np.reshape(s_total,(len(E_lower),len(E_upper)))
b_total = np.reshape(b_total,(len(E_lower),len(E_upper)))
b_err_total = np.reshape(b_err_total,(len(E_lower),len(E_upper)))
b_str_total = np.reshape(b_str_total,(len(E_lower),len(E_upper)))
min_time_idx = np.unravel_index(dwell_times.argmin(), dwell_times.shape)
print(df_list[min_time_idx[0]][min_time_idx[1]])
df_list[min_time_idx[0]][min_time_idx[1]].to_csv("results_energy_"+sig+".csv",sep=',',float_format='%.9e')
print("Dwell time = %f days"%dwell_times[min_time_idx])
print("Signal rate =",s_total[min_time_idx],"per day")
print("Background =",b_total[min_time_idx],"+/-",b_err_total[min_time_idx],"per day")
print("n100 = %f, n100 = %f"%(E_lower[min_time_idx[0]],E_upper[min_time_idx[1]]))
resultsfile = "results_%s.txt"%sig
results = "Signal = %s\nDwell time = %.3f days\nSignal rate = %.5f per day\nBackground rate = %.5f +/- %.5f per day\nn100 = %.3f\nn100 = %.3f\n%s"%(sig,dwell_times[min_time_idx],s_total[min_time_idx],b_total[min_time_idx],b_err_total[min_time_idx],E_lower[min_time_idx[0]],E_upper[min_time_idx[1]],b_str_total[min_time_idx])
with open(resultsfile,'a') as resfile:
    resfile.write(results+"\n")

for i in components:
    output = subprocess.run([os.path.dirname(__file__)+"/energy_write",path+"/"+"ML_classified.root",path+"/"+i+"_classified_energy.root",str(E_lower[min_time_idx[0]]),str(E_upper[min_time_idx[1]]),str(sig_dict[i])],\
        stdout=subprocess.PIPE,universal_newlines = True).stdout
