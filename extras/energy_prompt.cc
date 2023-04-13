#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TLeaf.h>
#include <iostream>
#include <TF1.h>
#include <TH2.h>

std::vector<double> energy(const char* infile, double e_lower, double e_upper, double dtw, int code_in){

  double ecut = e_upper;
  double ecut2 = e_lower;
  std::vector<double> energy_cut;

  // energy_cut.push_back(0);
  // energy_cut.push_back(0);

  gROOT->SetBatch(kTRUE);
  gROOT->ProcessLine( "gErrorIgnoreLevel = 1001;");

  TFile *f = new TFile(infile);
  if(!f->IsOpen()){
    printf("File %s does not exist.\n",infile);
    energy_cut.push_back(0);
    energy_cut.push_back(0);
    printf("Removed 0 out of 0 entries\n");
    return(energy_cut);
  }
  TTree *t_in = (TTree*)f->Get("data");

  const char* code_str = Form("code==%i",code_in);
  int nentries = t_in->GetEntries(code_str);
  int nentries_all = t_in->GetEntries();
  if (nentries==0){
      energy_cut.push_back(0);
      energy_cut.push_back(0);
      printf("Removed 0 out of %i entries\n",nentries);
      return(energy_cut);
  }
  //printf("There are %i entries\n",nentries);
  int subid,n100,mcid,mcid_2,closestPMT;
  double innerPE,mc_energy,code;

  int removed = 0;

  for(int i=1; i<nentries_all; i++){ //loops all values, not singles code value
      t_in->GetEntry(i-1);
      subid = t_in->GetLeaf("subid")->GetValue(0);
      mcid = t_in->GetLeaf("mcid")->GetValue(0);
      n100 = t_in->GetLeaf("n100")->GetValue(0);
      mc_energy = t_in->GetLeaf("mc_energy")->GetValue(0);
      innerPE = t_in->GetLeaf("innerPE")->GetValue(0);
      closestPMT = t_in->GetLeaf("closestPMT")->GetValue(0);
      code = t_in->GetLeaf("code")->GetValue(0);
      t_in->GetEntry(i);
      if(((n100>ecut or n100<ecut2) && subid==0 && code==code_in) || (closestPMT<dtw && code==code_in)){
          subid = t_in->GetLeaf("subid")->GetValue(0);
          mcid_2 = t_in->GetLeaf("mcid")->GetValue(0);
          removed++;
          if(subid==1 && mcid==mcid_2){
              removed++;
              i+=1;
          }
      }
  }
  printf("Removed %i out of %i entries\n",removed,nentries);
  energy_cut.clear();
  energy_cut.push_back(removed);
  energy_cut.push_back(nentries);
  return(energy_cut);
}

int main(int argc, char** argv){

  const char* infile = argv[1];
  double e_lower = std::stod(argv[2]);
  double e_upper = std::stod(argv[3]);
  double dtw = std::stod(argv[4]);
  int code_in = std::stoi(argv[5]);
  std::vector<double> energy_cut = energy(infile,e_lower,e_upper,dtw,code_in);

}
