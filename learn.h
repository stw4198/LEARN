#ifndef LEARN_H
#define LEARN_H

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TLeaf.h>
#include <sstream>
#include <fstream>
#include <iostream>

void pdf_gen(const char* file, const char* component, int nbins, int dTank, int rPMT, std::string x_win);
void merge_PDFs(const char* sig);
std::vector<double> likehood_classify(const char* infile, const char* component, std::string nx);
void likehood(const char* infile, const char* component, std::string x_win);
std::vector<const char*> input_select(int argc,char** argv);

#endif
