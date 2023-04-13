sig_dict = {'fn':1,\
			'geo':2,\
			'heysham_2':3,\
			'hinkley_c':4,\
			'sizewell_b':5,\
			'gravelines':6,\
			'li9':7,\
			'n17':8,\
			'torness':9,\
			'world':10,\
			'energy':11,\
			'hartlepool_1':12,\
			'hartlepool_2':13,\
			'heysham_full':14}

def sig_choice(sig):
    if sig == "heysham_2":
        signal_components = ["heysham_2"]
        background_components = ["torness","world","hinkley_c","sizewell_b","gravelines","geo","li9","n17","fn"]
    elif sig == "heysham_full":
        signal_components = ["heysham_full"]
        background_components = ["torness","world","hinkley_c","sizewell_b","gravelines","geo","li9","n17","fn"]
    elif sig == "hartlepool":
        signal_components = ["hartlepool_1","hartlepool_2"]
        background_components = ["heysham_full","torness","world","hinkley_c","sizewell_b","gravelines","geo","li9","n17","fn"]
    elif sig == "hartlepool_1":
        signal_components = ["hartlepool_1"]
        background_components = ["hartlepool_2","heysham_full","torness","world","hinkley_c","sizewell_b","gravelines","geo","li9","n17","fn"]
    elif sig == "heytor":
        signal_components = ["heysham_2","torness"]
        background_components = ["world","hinkley_c","sizewell_b","gravelines","geo","li9","n17","fn"]
    elif sig == "pwr":
        signal_components = ["hinkley_c","sizewell_b","gravelines"]
        background_components = ["world","geo","li9","n17","fn"]
    elif sig == "gravhink":
        signal_components = ["hinkley_c","gravelines"]
        background_components = ["sizewell_b","world","geo","li9","n17","fn"]
    elif sig == "sizehink":
        signal_components = ["hinkley_c","sizewell_b"]
        background_components = ["gravelines","world","geo","li9","n17","fn"]
    return (signal_components,background_components)