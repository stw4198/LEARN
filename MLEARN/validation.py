import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uproot
from sklearn.metrics import classification_report, ConfusionMatrixDisplay, roc_curve, roc_auc_score, confusion_matrix, auc
import joblib
import sys
import glob
import readline
import shutil
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
from mlparser import mlparser

#load the validation data
input = sys.argv[1]
sig = sys.argv[2]
tank = int(sys.argv[3])
try:
  bg = sys.argv[4]
except:
  bg = "fn"

#load the model
try:
  clf_route = '%s_finder.sav'%bg
except:
  print("Can not find %s_finder.sav. Exiting."%bg)
  exit()
clf = joblib.load(clf_route)

signal_components,background_components = sig_choice(sig)
components = signal_components+background_components

filenames = []
for i in components:
  for file in glob.glob(input+"/"+i+"*.root"):
      filenames.append(file)
# for file in glob.glob(input+"/"+bg+"*.root"):
#     print("Found",file)
#     filenames.append(file)
#print(filenames)

# input+="/*.root"
# location = os.path.join(input)
# filenames = sorted(glob.glob(location))
# print(filenames)

d = {}
z = {}

#var = ['n100','n100_prev', 'innerPE', 'dt_prev_us', 'drPrevr', 'x', 'y', 'z', 'closestPMT', 'closestPMT_prev', 'good_pos', 'good_pos_prev', 'subid']
var = ['n100','n100_prev','n9','n9_prev','innerPE','dt_prev_us','drPrevr','x','y','z','closestPMT','good_pos','good_pos_prev','subid']

#var = ['n100','n100_prev','n9', 'n9_prev', 'innerPE', 'dt_prev_us', 'drPrevr', 'x', 'y', 'z', 'closestPMT', 'closestPMT_prev', 'good_pos', 'good_pos_prev', 'good_dir', 'good_dir_prev', 'subid']
#var = ['n9', 'n9_prev', 'innerPE', 'dt_prev_us', 'drPrevr', 'x', 'y', 'z', 'closestPMT', 'closestPMT_prev', 'good_pos', 'subid']
ext_var = ['mcid','mc_energy']#,'n100','n100_prev','good_pos', 'good_pos_prev', 'good_dir', 'good_dir_prev']
cond = '(n100>0)&(n100_prev>0)'#&(subid==0)'

for f in filenames:
  with uproot.open(f) as f1:
    data = f1["data"]
    b, a = f.rsplit('/',1)
    c, e = a.split('_classified.root',1)
    d[c] = data.arrays(var,cond,library='pd')
    z[c] = data.arrays(ext_var,cond,library='pd')

for x, y in d.items():
  print("Found:",x)
  if (x == 'fn'):
    y['label'] = 1
    y['source'] = sig_dict['fn']
  elif (x == 'geo'):
    y['label'] = 0
    y['source'] = sig_dict['geo']
  elif (x == 'heysham_2'):
    y['label'] = 0
    y['source'] = sig_dict['heysham_2']
  # elif (x == 'heysham_full'):
  #   y['label'] = 0
  #   y['source'] = sig_dict['heysham_full']
  # elif (x == 'hartlepool_1'):
  #   y['label'] = 0
  #   y['source'] = sig_dict['hartlepool_1']
  # elif (x == 'hartlepool_2'):
  #   y['label'] = 0
  #   y['source'] = sig_dict['hartlepool_2']
  elif (x == 'hinkley_c'):
    y['label'] = 0
    y['source'] = sig_dict['hinkley_c']
  elif (x == 'sizewell_b'):
    y['label'] = 0
    y['source'] = sig_dict['sizewell_b']
  elif (x == 'gravelines'):
    y['label'] = 0
    y['source'] = sig_dict['gravelines']
  elif (x == 'li9'):
    y['label'] = 0
    y['source'] = sig_dict['li9']
  elif (x == 'n17'):
    y['label'] = 0
    y['source'] = sig_dict['n17']
  elif (x == 'torness'):
    y['label'] = 0
    y['source'] = sig_dict['torness']
  elif (x == 'world'):
    y['label'] = 0
    y['source'] = sig_dict['world']
  # elif (x == 'energy'):
  #   y['label'] = 0
  #   y['source'] = sig_dict['energy']
  else:
    print("%s not used, please move file"%x)
    sys.exit()

#neutron model
X = pd.concat(d.values(), ignore_index=True, keys=d.keys())
y = X[['label']]
y = y.to_numpy()
y = y.flatten()
Z = pd.concat(z.values(), ignore_index=True, keys=z.keys())
#print(y)

source = X[['source']]
source = source.to_numpy()
source = source.flatten()
X = X.drop(['label', 'source'], axis=1)
#print('\n\nData:')
#print(X)
print('\n\nApplying the model ...')
pred = clf.predict(X.values)
prob = clf.predict_proba(X.values)
scores = clf.decision_function(X.values)
print ('\n\nClassification Report:')
print (classification_report(y, pred))
cm = confusion_matrix(y, pred, labels=clf.classes_)
print ('\n\nConfusion Matrix: ')
print(cm)
print ('\n\nCreating figures ...')
#plt.rcParams.update({'font.size': 15})
ConfusionMatrixDisplay.from_predictions(y, pred)
#plt.title("Fast Neutron finder Validation")
plt.xlabel("Predicted label",fontsize=15)
plt.ylabel("True label",fontsize=15)
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
plt.savefig('cm_fnfinder_validation.pdf')
#plt.show(block=False)
#plt.pause(3)
plt.clf()

signal=prob[:,1]
fpr, tpr, _ = roc_curve(y,signal)
auc = auc(fpr, tpr)
plt.plot(fpr, tpr, marker=',', label='Fast Neutrons (area = {:.2f})'.format(auc))
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc='best')
#plt.title('Fast Neutron finder Validation')
plt.savefig('roc_fnfinder_validation.pdf')
#plt.show(block=False)
#plt.pause(3)
plt.clf()

Tol_bright = ['#4477AA','#66CCEE','#228833',
                     '#CCBB44','#EE6677','#AA3377','#BBBBBB','k']

X.loc[:,'classifier'] = pred
X.loc[:,'label'] = y
X.loc[:,'scores'] = scores
ts = X.loc[(X.classifier==0) & (X.label==0)]
tb = X.loc[(X.classifier==1) & (X.label==1)]
fs = X.loc[(X.classifier==0) & (X.label==1)]
fb = X.loc[(X.classifier==1) & (X.label==0)]
#plot the events by their decision function score
plt.hist(ts.scores.values.flatten(), bins=50, label='True Other', alpha=1,hatch='\\/',color=Tol_bright[0])#,fill=True)
plt.hist(tb.scores.values.flatten(), bins=50, label='True Fast Neutron', alpha=1,hatch='O',color=Tol_bright[6])
plt.hist(fs.scores.values.flatten(), bins=50, label='False Other', alpha=1,hatch='|',color=Tol_bright[4])
plt.hist(fb.scores.values.flatten(), bins=50, label='False Fast Neutron', alpha=.75,color=Tol_bright[3])
plt.yscale('log')
#plt.title('Fast Neutron Finder scores, validation')
plt.legend(loc='best')
plt.xlabel('Decision score',fontsize=15)
plt.ylabel('Frequency',fontsize=15)
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
plt.savefig('df_fnfinder_validation.pdf')
#plt.show(block=False)
#plt.pause(3)
plt.clf()

X.loc[:,'source'] = source
X.loc[:,'prob_fn'] = prob[:,1]
X.loc[:,'prob_other'] = prob[:,0]
for i in ext_var:
  X.loc[:,i] = Z[i]
print ('Added ML data:\n\n',X, '\n\n')
print('Outputting ML data to classified_valdata.csv...\n\n')
X.to_csv('classified_valdata.csv')
mlparser('classified_valdata.csv','results_learn.csv',sig,tank)
