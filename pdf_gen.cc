#include "learn.h"

void pdf_gen(const char* file, const char* component, int nbins, int dTank, int rPMT, std::string x_win){

  TFile *f = new TFile(file);
  if(!f->IsOpen()){
    printf("File %s does not exist.\n",file);
    return;
  }
  TTree *t = (TTree*)f->Get("data");
  TTree *run = (TTree*)f->Get("runSummary");
  
  TFile *like = new TFile(Form("%s_pdfs.root",component),"RECREATE");
   
  int nfiles = run->GetEntries();
  printf("%i files merged together\n",nfiles);
  
  int nevents = 0;
  int events;
  for(int i=0;i<nfiles;i++){
    run->GetEntry(i);
    events = run->GetLeaf("nEvents")->GetValue(0);
    nevents += events;
  }
  
  printf("There were %i MC events simulated\n",nevents);
  
  int nentries = t->GetEntries();
  printf("There are %i reconstructed entries in %s\n",nentries,file);
  
  //std::string x = "100";
  TString nx = "n" + x_win;
  TString nx_prev = nx + "_prev";

  int nkept = 0;
  
  TH1D* nX = new TH1D("nX","nX",nbins,0,1000);
  TH1D* nX_prev = new TH1D("nX_prev","nX_prev",nbins,0,1000);
  TH1D* dt_prev_us = new TH1D("dt_prev_us","dt_prev_us",nbins,0,10000);
  TH1D* drPrevr = new TH1D("drPrevr","drPrevr",nbins,0,dTank);
  TH1D* closestPMT = new TH1D("closestPMT","closestPMT",nbins,-499,rPMT);
  
  for(int i=0; i<nentries; i++){
    t->GetEntry(i);
    if (i%100000==0){
      printf("Generating PDFs: Event %d of %d\n",i,nentries);
    }
    if (t->GetLeaf(nx)->GetValue(0) > 0 and t->GetLeaf("closestPMT")->GetValue(0) > -499){
	    nkept++;
	    nX->Fill(t->GetLeaf(nx)->GetValue(0));
	    nX_prev->Fill(t->GetLeaf(nx_prev)->GetValue(0));
	    dt_prev_us->Fill(t->GetLeaf("dt_prev_us")->GetValue(0));
	    drPrevr->Fill(t->GetLeaf("drPrevr")->GetValue(0));
	    closestPMT->Fill(t->GetLeaf("closestPMT")->GetValue(0));
    }
  }
  
  printf("Kept %i events\n",nkept);
  
  nX->Scale(1/nX->GetEntries());
  nX->Write();
  nX_prev->Scale(1/nX_prev->GetEntries());
  nX_prev->Write();
  dt_prev_us->Scale(1/dt_prev_us->GetEntries());
  dt_prev_us->Write();
  drPrevr->Scale(1/drPrevr->GetEntries());
  drPrevr->Write();
  closestPMT->Scale(1/closestPMT->GetEntries());
  closestPMT->Write();
  like->Close();

}
