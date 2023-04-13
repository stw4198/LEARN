import numpy as np
import os
import sys

folder = sys.argv[1]
component = sys.argv[2]

# path joining version for other paths
DIR = folder
nfiles = len([name for name in os.listdir(folder) if os.path.isfile(os.path.join(folder, name))])
print("There are %i files in %s" % (nfiles,folder))

filenames = next(os.walk(folder), (None, None, []))[2]

file_pdf = []
file_eval = []

for i in range(nfiles):
  if i < nfiles/2:
    file_pdf.append(filenames[i])
  elif i >= nfiles/2:
    file_eval.append(filenames[i])
    
print("There are %i files for PDFs and %i files for evaluating"%(len(file_pdf),len(file_eval)))

hadd_pdf = "hadd -f %s_pdf_data.root"%component
hadd_eval = "hadd -f %s_eval_data.root"%component

for i in range(len(file_pdf)):
    file_add = " " + folder + "/" + file_pdf[i]
    hadd_pdf += file_add
for i in range(len(file_eval)):
    file_add = " " + folder + "/" + file_eval[i]
    hadd_eval += file_add
    
os.system(hadd_pdf)
os.system(hadd_eval)
