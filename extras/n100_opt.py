import subprocess
import os
import pandas as pd
from scipy import integrate
import numpy as np

path = "../learn_test/subid_test/analysis_heysham_2"
tank = 22
file = path + "/results_learn.csv"

components = ["heysham_2",\
"torness_full",\
"world",\
"geo",\
"li9",\
"n17",\
"fn",\
]
other = ["heysham_2",\
"torness_full",\
"world",\
"geo",\
"fn",\
]
signal_components = ["heysham_2"]
background_components = ["torness_full",\
"world",\
"geo",\
"li9",\
"n17",\
"fn",\
]

df = pd.read_csv(file)
df = df.set_index('Component')

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

E_lower = np.arange(20,26,1)
E_upper = np.arange(65,68,1)
dwell_times = []
s_total = []
b_total = []
b_err_total =[]

progress = 0

for x in E_lower:

    if progress/len(E_lower) % 0.05:
        print("%f%%"%(100*progress/len(E_lower)))

    dwell_times_upper = []
    s_upper = []
    b_upper = []
    b_err_upper = []
    for y in E_upper:
        print(x,y)
        sig_output_removed = []
        sig_output_total = []

        for i in range(len(signal_components)):
            output = subprocess.run([os.getcwd()+"/n100",path+"/"+signal_components[i]+"_classified.root",str(x),str(y)],\
                stdout=subprocess.PIPE,universal_newlines = True).stdout
            sig_output_removed.append(int(output.split()[1]))
            sig_output_total.append(int(output.split()[4]))
            df.loc[signal_components[i],str(tank)+" kept"]=int(output.split()[4])-int(output.split()[1])
        for i in range(len(background_components)):
            output = subprocess.run([os.getcwd()+"/n100",path+"/"+background_components[i]+"_classified.root",str(x),str(y)],\
                stdout=subprocess.PIPE,universal_newlines = True).stdout
            sig_output_removed.append(int(output.split()[1]))
            sig_output_total.append(int(output.split()[4]))
            df.loc[background_components[i],str(tank)+" kept"]=int(output.split()[4])-int(output.split()[1])

        rates = []
        fn_cl = 3.69
        for i in components:
            if i == 'fn':
                if df.loc[i,"%i kept"%tank] == 0:
                    df.loc[i,"%i kept"%tank]=fn_cl
                rates.append(df.loc[i,"%i kept"%tank]*df.loc[i,'%i rate'%tank]*86400/df.loc[i,'%i MC'%tank])
            else:
                if df.loc[i,"%i kept"%tank]==0 or df.loc[i,'%i rate'%tank]==0 or df.loc[i,'%i MC'%tank]==0:
                    rates.append(0)
                else:
                    rates.append(df.loc[i,"%i kept"%tank]*df.loc[i,'%i rate'%tank]*86400/df.loc[i,'%i MC'%tank])
        d = {'components':components,'rates':rates}
        df_analysis = pd.DataFrame(data=d)
        df_analysis = df_analysis.set_index('components')
        df_analysis.loc["li9"] = df_analysis.loc["li9"]*R_li9_cor
        df_analysis.loc["n17"] = df_analysis.loc["n17"]*R_n17_cor
        print(df_analysis)
        s = 0
        b = 0
        f = 0
        ferr = 0.27
        li9 = 0
        li9err = 0.002
        n17 = 0
        n17err = 0.002
        world = 0
        world_err = 0.06
        geo = 0
        geoerr = 0.25

        for i in signal_components:
            s+=df_analysis.loc[i,'rates']
        for k in background_components:
            if k!="li9" and k!="n17" and k!="world" and k!="geo" and k!="fn":
                b+=df_analysis.loc[k,'rates']
            elif k=="world":
                b+=df_analysis.loc[k,'rates']
                world+=df_analysis.loc[k,'rates']
            elif k=="geo":
                b+=df_analysis.loc[k,'rates']
                geo+=df_analysis.loc[k,'rates']
            elif k=="fn":
                b+=df_analysis.loc[k,'rates']
                f+=df_analysis.loc[k,'rates']
            elif k=="li9":
                b+=df_analysis.loc[k,'rates']
                li9+=df_analysis.loc[k,'rates']
            elif k=="n17":
                b+=df_analysis.loc[k,'rates']
                n17+=df_analysis.loc[k,'rates']

        #b_err_frac = np.sqrt((li9err*li9)**2 + (n17err*n17)**2 + (world_err*world)**2 + (geoerr*geo)**2 + (ferr*f)**2)
        #b_err = b*b_err_frac
        b_err = np.sqrt((li9err*li9)**2 + (n17err*n17)**2 + (world_err*world)**2 + (geoerr*geo)**2 + (ferr*f)**2)
        #t3sigma = 9*b/(s**2 - 9*(b_err))
        t3sigma = 9*b/(s**2 - 9*(b_err)**2)
        # print(df_analysis)
        # print("Signal rate =",s,"per day")
        # print("Background =",b,"+/-",b_err,"per day")
        # if t3sigma<0:
        #     print("Dwell time does not converge (Gaussian) for veto of %f s"%t_veto)
        # else:
        #     print("Dwell time (Gaussian) =",t3sigma,"days")
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

dwell_times = np.reshape(dwell_times,(len(E_lower),len(E_upper)))
s_total = np.reshape(s_total,(len(E_lower),len(E_upper)))
b_total = np.reshape(b_total,(len(E_lower),len(E_upper)))
b_err_total = np.reshape(b_err_total,(len(E_lower),len(E_upper)))
min_time_idx = np.unravel_index(dwell_times.argmin(), dwell_times.shape)
print("Dwell time = %f days"%dwell_times[min_time_idx])
print("Signal rate =",s_total[min_time_idx],"per day")
print("Background =",b_total[min_time_idx],"+/-",b_err_total[min_time_idx],"per day")
print("n100_min = %f, n100_max = %f"%(E_lower[min_time_idx[0]],E_upper[min_time_idx[1]]))