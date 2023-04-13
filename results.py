import numpy as np
import pandas as pd
import os
import sys

file = sys.argv[1]
component = sys.argv[2]
tank = sys.argv[3]
kept = sys.argv[4]
MC = sys.argv[5]
rate = sys.argv[6]

def results_csv(file, component,tank,kept,MC,rate):

    path = os.getcwd() + "/" + file

    key_2_MC = tank + " MC"
    key_2_kept = tank + " kept"
    key_2_rate = tank + " rate"

    if os.path.isfile(path):
        print("%s found\n\n"%path)
        df = pd.read_csv(file)
        if key_2_MC in df.keys() and key_2_kept in df.keys() and key_2_rate in df.keys():
            print("%s m keys found\n\n"%tank)
        else:
            print("%s m keys not found\nAdding new keys"%tank)
            df[key_2_MC] = ""
            df[key_2_kept] = ""
            df[key_2_rate] = ""
        if df.isin({'Component':[component]}).any().any():
            print("%s found\n\n"%component)
        else:
            print("%s not found\n\n"%component)
            new_component = [component]
            for i in range(len(df.keys())-1):
                new_component.append("")
            df.loc[df.index.max() + 1, :] = new_component
        index = df.index
        condition = df["Component"] == component
        component_indices = index[condition]
        component_indices_list = component_indices.tolist()
        df.at[component_indices_list[0],"%s MC"%tank] = MC
        df.at[component_indices_list[0],"%s kept"%tank] = kept
        df.at[component_indices_list[0],"%s rate"%tank] = rate
        df = df.set_index('Component')
        print(df)
        print("\n\nWriting dataframe to %s"%file)
        df.to_csv(file,sep=',',float_format='%.9e')
    else:
        print("%s not found\n\n"%path)
        df=pd.DataFrame()
        df['Component'] = ""
        df[key_2_MC] = ""
        df[key_2_kept] = ""
        df[key_2_rate] = ""
        new_component = [component]
        for i in range(len(df.keys())-1):
            new_component.append("")
        df.loc[df.index.max() + 1, :] = new_component
        index = df.index
        condition = df["Component"] == component
        component_indices = index[condition]
        component_indices_list = component_indices.tolist()
        df.at[component_indices_list[0],"%s MC"%tank] = MC
        df.at[component_indices_list[0],"%s kept"%tank] = kept
        df.at[component_indices_list[0],"%s rate"%tank] = rate
        df = df.set_index('Component')
        print(df)
        print("\n\nWriting dataframe to %s"%file)
        df.to_csv(file,sep=',',float_format='%.9e')

    return

results_csv(file,component,tank,kept,MC,rate)