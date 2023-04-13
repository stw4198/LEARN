import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from scipy import integrate

file = sys.argv[1]
tank = int(sys.argv[2])
sig = sys.argv[3]

muon_rate = 4.09e-8
if tank == 16:
    pmt = 5.7
elif tank == 22:
    pmt = 9.0
elif tank == 28:
    pmt = 12.0
    
#Read in rates as csv
df = pd.read_csv(file)
df = df.set_index('Component')

if sig == "heysham_2":
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
elif sig == "heysham_full":
    components = ["heysham_full",\
    "torness_full",\
    "world",\
    "geo",\
    "li9",\
    "n17",\
    "fn",\
    ]
    other = ["heysham_full",\
    "heysham_2",\
    "torness_full",\
    "world",\
    "geo",\
    "fn",\
    ]
    signal_components = ["heysham_full"]
    background_components = ["torness_full",\
    "world",\
    "geo",\
    "li9",\
    "n17",\
    "fn",\
    ]
elif sig == "torness_full":
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
    signal_components = ["torness_full"]
    background_components = ["heysham_2",\
    "world",\
    "geo",\
    "li9",\
    "n17",\
    "fn",\
    ]
elif sig == "torness_only":
    components = ["torness_full",\
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
    signal_components = ["torness_full"]
    background_components = ["world",\
    "geo",\
    "li9",\
    "n17",\
    "fn",\
    ]
elif sig == "heytor":
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
    signal_components = ["heysham_2","torness_full"]
    background_components = ["world",\
    "geo",\
    "li9",\
    "n17",\
    "fn",\
    ]
elif sig == "hartlepool":
    components = ["hartlepool_1",\
    "hartlepool_2",\
    "heysham_full",\
    "torness_full",\
    "world",\
    "geo",\
    "li9",\
    "n17",\
    "fn",\
    ]
    other = ["hartlepool_1",\
    "hartlepool_2",\
    "heysham_full",\
    "torness_full",\
    "world",\
    "geo",\
    "fn",\
    ]
    signal_components = ["hartlepool_1","hartlepool_2"]
    background_components = ["heysham_full",\
    "torness_full",\
    "world",\
    "geo",\
    "li9",\
    "n17",\
    "fn",\
    ]
elif sig == "hartlepool_1":
    components = ["hartlepool_1",\
    "hartlepool_2",\
    "heysham_full",\
    "torness_full",\
    "world",\
    "geo",\
    "li9",\
    "n17",\
    "fn",\
    ]
    other = ["hartlepool_1",\
    "hartlepool_2",\
    "heysham_full",\
    "torness_full",\
    "world",\
    "geo",\
    "fn",\
    ]
    signal_components = ["hartlepool_1"]
    background_components = ["hartlepool_2",\
    "heysham_full",\
    "torness_full",\
    "world",\
    "geo",\
    "li9",\
    "n17",\
    "fn",\
    ]

radio = ["li9",\
"n17",\
]

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

signal = 0
background = 0

d = {'components':components,'rates':rates}
df_muon = pd.DataFrame(data=d)
df_muon = df_muon.set_index('components')

muon_eff = 0.95
t_half_n17 = 4.173
t_half_li9 = 0.1783
def integrand_n17(t):
    return np.exp(-np.log(2)*t/t_half_n17)
def integrand_li9(t):
    return np.exp(-np.log(2)*t/t_half_li9)
t_dead = np.arange(0,15.5,1)

t3sigma = []
s_ = []
b_ = []
b_err_ = []
for i in range(len(t_dead)):
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
    R_n17_cor = 1-muon_eff+muon_eff*integrate.quad(integrand_n17,t_dead[i],np.infty)[0]/integrate.quad(integrand_n17,0,np.infty)[0]
    R_li9_cor = 1-muon_eff+muon_eff*integrate.quad(integrand_li9,t_dead[i],np.infty)[0]/integrate.quad(integrand_li9,0,np.infty)[0]

    for i in signal_components:
        s+=df_muon.loc[i,'rates']
    for k in background_components:
        if k!="li9" and k!="n17" and k!="world" and k!="geo" and k!="fn":
            b+=df_muon.loc[k,'rates']
        elif k=="world":
            b+=df_muon.loc[k,'rates']
            world+=df_muon.loc[k,'rates']
        elif k=="geo":
            b+=df_muon.loc[k,'rates']
            geo+=df_muon.loc[k,'rates']
        elif k=="fn":
            b+=df_muon.loc[k,'rates']
            f+=df_muon.loc[k,'rates']
        elif k=="li9":
            b+=df_muon.loc[k,'rates']*R_li9_cor
            li9+=df_muon.loc[k,'rates']*R_li9_cor
        elif k=="n17":
            b+=df_muon.loc[k,'rates']*R_n17_cor
            n17+=df_muon.loc[k,'rates']*R_n17_cor

    b_err_1 = np.sqrt((li9err*li9)**2 + (n17err*n17)**2 + (world_err*world)**2 + (geoerr*geo)**2 + (ferr*f)**2)
    b_err = b*b_err_1
    t_3 = 9*b/(s**2 - 9*(b_err)**2)
    t3sigma.append(t_3)
    s_.append(s)
    b_.append(b)
    b_err_.append(b_err)

plt.plot(t_dead,t3sigma)
plt.ylim(0,max(t3sigma))
plt.xlim(np.argmax(t3sigma),max(t_dead))
plt.grid()
plt.ylabel("Time to 3\u03C3 [days]")
plt.xlabel("Deadtime [s]")
#plt.title("%s"%sig)
plt.clf()

print("Signal rate =",s_[1],"per day")
print("Background =",b_[1],"+/-",b_err,"per day")

if t3sigma[1]<0:
    print("Dwell time does not converge (Gaussian) for veto of %f s"%t_dead[1])
else:
    print("Dwell time (Gaussian) =",t3sigma[1],"days")

if sig=='hartlepool' or sig=='hartlepool_1':
    t = np.arange(0.01,30,0.01)
    s = s_[1]
    b = b_[1]
    b_err = b_err_[1]
    sigma = np.sqrt(2*((s*t + b*t)*np.log((s*t + b*t)*(b*t + b_err**2*t**2)/(b**2*t**2 + (s*t + b*t)*b_err**2*t**2))\
        - (b**2*t**2)*np.log(1 + (b_err**2*t**2*s*t)/(b*t*(b*t + b_err**2*t**2)))/(b_err**2*t**2)))
    if max(sigma)<3:
        print("Dwell time does not converge (Poisson) for veto of %f s"%t_dead[1])
    else:
        time_index = [ n for n,i in enumerate(sigma) if i>3][0]
        print("Dwell time (Poisson) =",t[time_index],"days")