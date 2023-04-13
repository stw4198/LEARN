import pandas as pd
import numpy as np
import sys
try:
    import ROOT
except:
    print("Please source thisroot.sh")
    exit()
from array import array
from tqdm import tqdm
import os
sys.path.insert(1, os.path.dirname(os.path.realpath(__file__))+"/../")
from sig_choice import sig_choice,sig_dict

file = sys.argv[1]
sig = sys.argv[2]
tank = sys.argv[3]

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)

def mlparser(file,learn_file,sig,tank):
    df = pd.read_csv(file)
    f = ROOT.TFile("ML_classified.root", "RECREATE")
    tree = ROOT.TTree("data", "ML Classified Data")
    mc_energy = array('f', [0.])
    subid = array('f', [0.])
    n100 = array('f', [0.])
    n100_prev = array('f', [0.])
    mcid = array('f', [0.])
    innerPE = array('f', [0.])
    code = array('i', [0])
    closestPMT = array('f', [0.])
    tree.Branch("mc_energy", mc_energy, 'mc_energy/F')
    tree.Branch("subid", subid, 'subid/F')
    tree.Branch("n100", n100, 'n100/F')
    tree.Branch("n100_prev", n100_prev, 'n100_prev/F')
    tree.Branch("mcid", mcid, 'mcid/F')
    tree.Branch("innerPE", innerPE, 'innerPE/F')
    tree.Branch("code", code, 'code/I')
    tree.Branch("closestPMT", closestPMT, 'closestPMT/F')

    fn_tot = fn_keep = 0
    heysham_2_tot = heysham_2_keep = 0
    torness_tot = torness_keep = 0
    sizewell_b_tot = sizewell_b_keep = 0
    hinkley_c_tot = hinkley_c_keep = 0
    gravelines_tot = gravelines_keep = 0
    world_tot = world_keep = 0
    geo_tot = geo_keep = 0
    li9_tot = li9_keep = 0
    n17_tot = n17_keep = 0
    heysham_full_tot = heysham_full_keep = 0
    hartlepool_1_tot = hartlepool_1_keep = 0
    hartlepool_2_tot = hartlepool_2_keep = 0
    energy_tot = energy_keep = 0

    for i in tqdm(range(len(df)),desc='Parsing %s'%file):
        if df['source'][i]==sig_dict['fn']:
            fn_tot+=1
            if df['classifier'][i]==1:
                continue
            elif i>0 and df['classifier'][i]==0 and df['classifier'][i-1]==1 and df['mcid'][i]==df['mcid'][i-1]:
                continue
            elif i<len(df)-1 and df['classifier'][i]==0 and df['classifier'][i+1]==1 and df['mcid'][i]==df['mcid'][i+1]:
                continue
            else:
                fn_keep+=1
                mc_energy[0] = df['mc_energy'][i]
                n100[0] = df['n100'][i]
                n100_prev[0] = df['n100_prev'][i]
                subid[0] = df['subid'][i]
                mcid[0] = df['mcid'][i]
                innerPE[0] = df['innerPE'][i]
                code[0] = df['source'][i]
                closestPMT[0] = df['closestPMT'][i]
                tree.Fill()
        elif df['source'][i]==sig_dict['geo']:
            geo_tot+=1
            if df['classifier'][i]==1:
                continue
            elif i>0 and df['classifier'][i]==0 and df['classifier'][i-1]==1 and df['mcid'][i]==df['mcid'][i-1]:
                continue
            elif i<len(df)-1 and df['classifier'][i]==0 and df['classifier'][i+1]==1 and df['mcid'][i]==df['mcid'][i+1]:
                continue
            else:
                geo_keep+=1
                mc_energy[0] = df['mc_energy'][i]
                n100[0] = df['n100'][i]
                n100_prev[0] = df['n100_prev'][i]
                subid[0] = df['subid'][i]
                mcid[0] = df['mcid'][i]
                innerPE[0] = df['innerPE'][i]
                code[0] = df['source'][i]
                closestPMT[0] = df['closestPMT'][i]
                tree.Fill()
        elif df['source'][i]==sig_dict['heysham_2']:
            heysham_2_tot+=1
            if df['classifier'][i]==1:
                continue
            elif i>0 and df['classifier'][i]==0 and df['classifier'][i-1]==1 and df['mcid'][i]==df['mcid'][i-1]:
                continue
            elif i<len(df)-1 and df['classifier'][i]==0 and df['classifier'][i+1]==1 and df['mcid'][i]==df['mcid'][i+1]:
                continue
            else:
                heysham_2_keep+=1
                mc_energy[0] = df['mc_energy'][i]
                n100[0] = df['n100'][i]
                n100_prev[0] = df['n100_prev'][i]
                subid[0] = df['subid'][i]
                mcid[0] = df['mcid'][i]
                innerPE[0] = df['innerPE'][i]
                code[0] = df['source'][i]
                closestPMT[0] = df['closestPMT'][i]
                tree.Fill()
        elif df['source'][i]==sig_dict['hinkley_c']:
            hinkley_c_tot+=1
            if df['classifier'][i]==1:
                continue
            elif i>0 and df['classifier'][i]==0 and df['classifier'][i-1]==1 and df['mcid'][i]==df['mcid'][i-1]:
                continue
            elif i<len(df)-1 and df['classifier'][i]==0 and df['classifier'][i+1]==1 and df['mcid'][i]==df['mcid'][i+1]:
                continue
            else:
                hinkley_c_keep+=1
                mc_energy[0] = df['mc_energy'][i]
                n100[0] = df['n100'][i]
                subid[0] = df['subid'][i]
                mcid[0] = df['mcid'][i]
                innerPE[0] = df['innerPE'][i]
                code[0] = df['source'][i]
                closestPMT[0] = df['closestPMT'][i]
                tree.Fill()
        elif df['source'][i]==sig_dict['sizewell_b']:
            sizewell_b_tot+=1
            if df['classifier'][i]==1:
                continue
            elif i>0 and df['classifier'][i]==0 and df['classifier'][i-1]==1 and df['mcid'][i]==df['mcid'][i-1]:
                continue
            elif i<len(df)-1 and df['classifier'][i]==0 and df['classifier'][i+1]==1 and df['mcid'][i]==df['mcid'][i+1]:
                continue
            else:
                sizewell_b_keep+=1
                mc_energy[0] = df['mc_energy'][i]
                n100[0] = df['n100'][i]
                n100_prev[0] = df['n100_prev'][i]
                subid[0] = df['subid'][i]
                mcid[0] = df['mcid'][i]
                innerPE[0] = df['innerPE'][i]
                code[0] = df['source'][i]
                closestPMT[0] = df['closestPMT'][i]
                tree.Fill()
        elif df['source'][i]==sig_dict['gravelines']:
            gravelines_tot+=1
            if df['classifier'][i]==1:
                continue
            elif i>0 and df['classifier'][i]==0 and df['classifier'][i-1]==1 and df['mcid'][i]==df['mcid'][i-1]:
                continue
            elif i<len(df)-1 and df['classifier'][i]==0 and df['classifier'][i+1]==1 and df['mcid'][i]==df['mcid'][i+1]:
                continue
            else:
                gravelines_keep+=1
                mc_energy[0] = df['mc_energy'][i]
                n100[0] = df['n100'][i]
                n100_prev[0] = df['n100_prev'][i]
                subid[0] = df['subid'][i]
                mcid[0] = df['mcid'][i]
                innerPE[0] = df['innerPE'][i]
                code[0] = df['source'][i]
                closestPMT[0] = df['closestPMT'][i]
                tree.Fill()
        elif df['source'][i]==sig_dict['li9']:
            li9_tot+=1
            if df['classifier'][i]==1:
                continue
            elif i>0 and df['classifier'][i]==0 and df['classifier'][i-1]==1 and df['mcid'][i]==df['mcid'][i-1]:
                continue
            elif i<len(df)-1 and df['classifier'][i]==0 and df['classifier'][i+1]==1 and df['mcid'][i]==df['mcid'][i+1]:
                continue
            else:
                li9_keep+=1
                mc_energy[0] = df['mc_energy'][i]
                n100[0] = df['n100'][i]
                n100_prev[0] = df['n100_prev'][i]
                subid[0] = df['subid'][i]
                mcid[0] = df['mcid'][i]
                innerPE[0] = df['innerPE'][i]
                code[0] = df['source'][i]
                closestPMT[0] = df['closestPMT'][i]
                tree.Fill()
        elif df['source'][i]==sig_dict['n17']:
            n17_tot+=1
            if df['classifier'][i]==1:
                continue
            elif i>0 and df['classifier'][i]==0 and df['classifier'][i-1]==1 and df['mcid'][i]==df['mcid'][i-1]:
                continue
            elif i<len(df)-1 and df['classifier'][i]==0 and df['classifier'][i+1]==1 and df['mcid'][i]==df['mcid'][i+1]:
                continue
            else:
                n17_keep+=1
                mc_energy[0] = df['mc_energy'][i]
                n100[0] = df['n100'][i]
                n100_prev[0] = df['n100_prev'][i]
                subid[0] = df['subid'][i]
                mcid[0] = df['mcid'][i]
                innerPE[0] = df['innerPE'][i]
                code[0] = df['source'][i]
                closestPMT[0] = df['closestPMT'][i]
                tree.Fill()
        elif df['source'][i]==sig_dict['torness']:
            torness_tot+=1
            if df['classifier'][i]==1:
                continue
            elif i>0 and df['classifier'][i]==0 and df['classifier'][i-1]==1 and df['mcid'][i]==df['mcid'][i-1]:
                continue
            elif i<len(df)-1 and df['classifier'][i]==0 and df['classifier'][i+1]==1 and df['mcid'][i]==df['mcid'][i+1]:
                continue
            else:
                torness_keep+=1
                mc_energy[0] = df['mc_energy'][i]
                n100[0] = df['n100'][i]
                n100_prev[0] = df['n100_prev'][i]
                subid[0] = df['subid'][i]
                mcid[0] = df['mcid'][i]
                innerPE[0] = df['innerPE'][i]
                code[0] = df['source'][i]
                closestPMT[0] = df['closestPMT'][i]
                tree.Fill()
        elif df['source'][i]==sig_dict['world']:
            world_tot+=1
            if df['classifier'][i]==1:
                continue
            elif i>0 and df['classifier'][i]==0 and df['classifier'][i-1]==1 and df['mcid'][i]==df['mcid'][i-1]:
                continue
            elif i<len(df)-1 and df['classifier'][i]==0 and df['classifier'][i+1]==1 and df['mcid'][i]==df['mcid'][i+1]:
                continue
            else:
                world_keep+=1
                mc_energy[0] = df['mc_energy'][i]
                n100[0] = df['n100'][i]
                n100_prev[0] = df['n100_prev'][i]
                subid[0] = df['subid'][i]
                mcid[0] = df['mcid'][i]
                innerPE[0] = df['innerPE'][i]
                code[0] = df['source'][i]
                closestPMT[0] = df['closestPMT'][i]
                tree.Fill()
        elif df['source'][i]==sig_dict['hartlepool_1']:
            hartlepool_1_tot+=1
            if df['classifier'][i]==1:
                continue
            elif i>0 and df['classifier'][i]==0 and df['classifier'][i-1]==1 and df['mcid'][i]==df['mcid'][i-1]:
                continue
            elif i<len(df)-1 and df['classifier'][i]==0 and df['classifier'][i+1]==1 and df['mcid'][i]==df['mcid'][i+1]:
                continue
            else:
                hartlepool_1_keep+=1
                mc_energy[0] = df['mc_energy'][i]
                n100[0] = df['n100'][i]
                n100_prev[0] = df['n100_prev'][i]
                subid[0] = df['subid'][i]
                mcid[0] = df['mcid'][i]
                innerPE[0] = df['innerPE'][i]
                code[0] = df['source'][i]
                closestPMT[0] = df['closestPMT'][i]
                tree.Fill()
        elif df['source'][i]==sig_dict['hartlepool_2']:
            hartlepool_2_tot+=1
            if df['classifier'][i]==1:
                continue
            elif i>0 and df['classifier'][i]==0 and df['classifier'][i-1]==1 and df['mcid'][i]==df['mcid'][i-1]:
                continue
            elif i<len(df)-1 and df['classifier'][i]==0 and df['classifier'][i+1]==1 and df['mcid'][i]==df['mcid'][i+1]:
                continue
            else:
                hartlepool_2_keep+=1
                mc_energy[0] = df['mc_energy'][i]
                n100[0] = df['n100'][i]
                n100_prev[0] = df['n100_prev'][i]
                subid[0] = df['subid'][i]
                mcid[0] = df['mcid'][i]
                innerPE[0] = df['innerPE'][i]
                code[0] = df['source'][i]
                closestPMT[0] = df['closestPMT'][i]
                tree.Fill()
        elif df['source'][i]==sig_dict['heysham_full']:
            heysham_full_tot+=1
            if df['classifier'][i]==1:
                continue
            elif i>0 and df['classifier'][i]==0 and df['classifier'][i-1]==1 and df['mcid'][i]==df['mcid'][i-1]:
                continue
            elif i<len(df)-1 and df['classifier'][i]==0 and df['classifier'][i+1]==1 and df['mcid'][i]==df['mcid'][i+1]:
                continue
            else:
                heysham_full_keep+=1
                mc_energy[0] = df['mc_energy'][i]
                n100[0] = df['n100'][i]
                n100_prev[0] = df['n100_prev'][i]
                subid[0] = df['subid'][i]
                mcid[0] = df['mcid'][i]
                innerPE[0] = df['innerPE'][i]
                code[0] = df['source'][i]
                closestPMT[0] = df['closestPMT'][i]
                tree.Fill()
        # elif df['source'][i]==sig_dict['energy']:
        #     energy_tot+=1
        #     if df['classifier'][i]==1:
        #         continue
        #     elif i>0 and df['classifier'][i]==0 and df['classifier'][i-1]==1 and df['mcid'][i]==df['mcid'][i-1]:
        #         continue
        #     elif i<len(df)-1 and df['classifier'][i]==0 and df['classifier'][i+1]==1 and df['mcid'][i]==df['mcid'][i+1]:
        #         continue
        #     else:
        #         energy_keep+=1
        #         mc_energy[0] = df['mc_energy'][i]
        #         n100[0] = df['n100'][i]
        #         n100_prev[0] = df['n100_prev'][i]
        #         subid[0] = df['subid'][i]
        #         mcid[0] = df['mcid'][i]
        #         innerPE[0] = df['innerPE'][i]
        #         code[0] = df['source'][i]
        #         tree.Fill()
    tree.Write("data", ROOT.TObject.kOverwrite)
    f.Close()

    #sig='heysham_2'
    signal_components,background_components = sig_choice(sig)
    components = signal_components+background_components

    #tank = 22
    #learn_file = 'results_learn.csv'
    df = pd.read_csv(learn_file)
    df = df.set_index('Component')

    keep_dict = {"fn":fn_keep,\
                "geo":geo_keep,\
                "heysham_2":heysham_2_keep,\
                "torness":torness_keep,\
                "sizewell_b":sizewell_b_keep,\
                "hinkley_c":hinkley_c_keep,\
                "gravelines":gravelines_keep,\
                "world":world_keep,\
                "li9":li9_keep,\
                "n17":n17_keep,\
                "energy":energy_keep,\
                "hartlepool_1":hartlepool_1_keep,\
                "hartlepool_2":hartlepool_2_keep,\
                "heysham_full":heysham_full_keep}
    for i in components:
        keep = keep_dict[i]
        df.loc[i,"%i kept"%tank] = keep
    
    print("\nWriting ML data to results_ML.csv\n")
    
    print(df)
    df.to_csv("results_ML.csv",sep=',',float_format='%.9e')
    return

if __name__ == "__main__":
    file = sys.argv[1]
    learn_file = 'results_learn.csv'
    sig = sys.argv[2]
    tank = int(sys.argv[3])
    mlparser(file,learn_file,sig,tank)