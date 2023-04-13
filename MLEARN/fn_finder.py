import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uproot
from sklearn.model_selection import train_test_split
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import classification_report, ConfusionMatrixDisplay, roc_curve, roc_auc_score, confusion_matrix, auc
import joblib
import sys
import os
import glob
import readline
import uproot
sys.path.insert(1, os.path.dirname(os.path.realpath(__file__))+"/../")
from sig_choice import sig_choice,sig_dict

input = sys.argv[1]#+"/*.root"#input('Enter the path to the training data directory.\n End the file path with "/*.root".\n')
sig = sys.argv[2]
try:
  bg = sys.argv[3]
except:
  bg = "fn"
filenames = []
signal_components, background_components = sig_choice(sig)
#print(signal_components)
for i in signal_components:
  for file in glob.glob(input+"/"+i+"*.root"):
      filenames.append(file)
for file in glob.glob(input+"/"+bg+"*.root"):
    filenames.append(file)
if len(filenames)<2:
  print("Not enough files found to train. Found: %s"%filenames)
  exit()
#sig_file = input+"/"+sig+"_classified.root"
#input+="/*.root"
#location = os.path.join(input)
#filenames = glob.glob(location)
#print(input,filenames)
d = {}

var = ['n100','n100_prev', 'innerPE', 'dt_prev_us', 'drPrevr', 'x', 'y', 'z', 'closestPMT', 'closestPMT_prev', 'good_pos', 'good_pos_prev', 'subid']
#var = ['n100','n100_prev','n9', 'n9_prev', 'innerPE', 'dt_prev_us', 'drPrevr', 'x', 'y', 'z', 'closestPMT', 'closestPMT_prev', 'good_pos', 'good_pos_prev', 'good_dir', 'good_dir_prev', 'subid']
#var = ['n9', 'n9_prev', 'innerPE', 'dt_prev_us', 'drPrevr', 'x', 'y', 'z', 'good_pos', 'good_pos_prev', 'subid']
#var = ['n9', 'n9_prev', 'innerPE', 'dt_prev_us', 'drPrevr', 'x', 'y', 'z', 'closestPMT', 'closestPMT_prev', 'good_pos', 'subid']
cond = '(n100>0)&(n100_prev>0)'#&(subid==0)'

for f in filenames:
  with uproot.open(f) as f1:
    data = f1["data"]
    b, a = f.rsplit('/',1)
    c, e = a.split('_classified.root',1)
    d[c] = data.arrays(var,cond,library='pd')
for x, y in d.items():
  print("Found:",x)
  if (x == bg):
    y['label'] = 1
  else:
    y['label'] = 0

#neutron model
X = pd.concat(d.values(), ignore_index=True, keys=d.keys())
ydf = X[['label']]
y = ydf.to_numpy()
y = y.flatten()
X = X.drop(['label'], axis=1)
train_X, test_X, train_y, test_y = train_test_split(X, y, stratify=y)
clf = AdaBoostClassifier(DecisionTreeClassifier(max_depth=2),
                                n_estimators=100, learning_rate=0.1)
print('Fitting classifier...')
clf.fit(train_X.values, train_y)
pred = clf.predict(test_X.values)
prob = clf.predict_proba(test_X.values)
scores = clf.decision_function(test_X.values)

print('Saving classifier...')
clf_file = bg+'_finder.sav'
joblib.dump(clf, clf_file)

cm= confusion_matrix(test_y, pred, labels=clf.classes_)
print('Confusion Matrix:')
print(cm)
print('Classification Report:')
print(classification_report(test_y, pred))
print('Creating Figures...')
ConfusionMatrixDisplay.from_predictions(test_y, pred)
plt.title("Fast Neutron finder, training")
plt.savefig('cm_%sfinder_training.pdf'%bg)
#plt.show(block=False)
#plt.pause(3)
plt.clf()

test_X.loc[:,'classifier'] = pred
rows = test_X.index
ydf.index = range(len(ydf))
labels = ydf.iloc[rows,:]
test_X.loc[:,'label'] = labels
test_X.loc[:,'scores'] = scores
signal=prob[:,1]
fpr, tpr, _ = roc_curve(labels,signal)
auc = auc(fpr, tpr)
plt.plot(fpr, tpr, marker=',', label='Fast Neutrons (area = {:.2f})'.format(auc))

plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc='best')
plt.title('Fast Neutron finder, training')
plt.savefig('roc_%sfinder_training.pdf'%bg)
#plt.show(block=False)
#plt.pause(3)
plt.clf()

ts = test_X.loc[(test_X.classifier==0) & (test_X.label==0)]
tb = test_X.loc[(test_X.classifier==1) & (test_X.label==1)]
fs = test_X.loc[(test_X.classifier==0) & (test_X.label==1)]
fb = test_X.loc[(test_X.classifier==1) & (test_X.label==0)]
#plot the events by their decision function score
plt.hist(ts.scores.values.flatten(), bins=50, label='True Other', alpha=.5)
plt.hist(tb.scores.values.flatten(), bins=50, label='True Fast Neutron', alpha=.5)
plt.hist(fs.scores.values.flatten(), bins=50, label='False Other', alpha=.5)
plt.hist(fb.scores.values.flatten(), bins=50, label='False Fast Neutron', alpha=.5)
plt.yscale('log')
plt.title('Fast Neutron Finder scores, training')
plt.legend(loc='best')
plt.xlabel('Decision scores')
plt.ylabel('Frequency (log)')
plt.savefig('df_%sfinder_training.pdf'%bg)
#plt.show(block=False)
#plt.pause(3)
plt.clf()
