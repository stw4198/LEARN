#include "learn.h"

void likehood(const char* infile, const char* component, std::string x_win){

  //std::string x = "100";
  TString nx = "n" + x_win;
  TString nx_prev = nx + "_prev";

  TFile *f = new TFile(infile);
  if(!f->IsOpen()){
    printf("File %s does not exist.\n",infile);
    return;
  }
  TTree *t = (TTree*)f->Get("data");
  TTree *run = (TTree*)f->Get("runSummary");
  
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
  printf("There are %i reconstructed entries in %s\n",nentries,infile);

  TFile *out=new TFile(Form("%s_likelihoods.root",component),"RECREATE");

  int nbins = sqrt(nentries);

  TH1D* signal_like = new TH1D("signal_like","signal_like",nbins,0,0);
  TH1D* background_like = new TH1D("background_like","background_like",nbins,0,0);
  TH1D* ratio_like = new TH1D("ratio_like","ratio_like",nbins,0,0);

  //signal pdfs
  TFile* signal = new TFile("signal_pdfs.root");
  if(!signal->IsOpen()){
    printf("File signal_pdfs.root does not exist.\n");
    return;
  }
  TH1D* nX_signal = (TH1D*)signal->Get("nX");
  TH1D* nX_prev_signal = (TH1D*)signal->Get("nX_prev");
  TH1D* dt_prev_us_signal = (TH1D*)signal->Get("dt_prev_us");
  TH1D* drPrevr_signal = (TH1D*)signal->Get("drPrevr");
  TH1D* closestPMT_signal = (TH1D*)signal->Get("closestPMT");

  //background pdfs
  TFile * background = new TFile("background_pdfs.root");
  if(!background->IsOpen()){
    printf("File background_pdfs.root does not exist.\n");
    return;
  }
  TH1D* nX_background = (TH1D*)background->Get("nX");
  TH1D* nX_prev_background = (TH1D*)background->Get("nX_prev");
  TH1D* dt_prev_us_background = (TH1D*)background->Get("dt_prev_us");
  TH1D* drPrevr_background = (TH1D*)background->Get("drPrevr");
  TH1D* closestPMT_background = (TH1D*)background->Get("closestPMT");

  TFile *output = new TFile(Form("%s_like.root",component),"RECREATE");
  TTree *data = new TTree("data","low-energy detector triggered events");
  double r_like;
  double n100,n100_prev;
  double closestPMT;
  double drPrevr,dt_prev_us;
  double mc_energy;
  int subid;
  data->Branch("n100",&n100,"n100/D");
  data->Branch("n100_prev",&n100_prev,"n100_prev/D");
  data->Branch("closestPMT",&closestPMT,"closestPMT/D");
  data->Branch("drPrevr",&drPrevr,"drPrevr/D");
  data->Branch("dt_prev_us",&dt_prev_us,"dt_prev_us/D");
  data->Branch("subid",&subid,"subid/I");
  data->Branch("mc_energy",&mc_energy,"mc_energy/D");
  //TTree *data_next = new TTree("data_next","low-energy detector triggered events");
  double n100_next,n100_prev_next;
  double closestPMT_next;
  double drPrevr_next,dt_prev_us_next;
  double mc_energy_next;
  int subid_next;
  data->Branch("n100_next",&n100_next,"n100_next/D");
  data->Branch("n100_prev_next",&n100_prev_next,"n100_prev_next/D");
  data->Branch("closestPMT_next",&closestPMT_next,"closestPMT_next/D");
  data->Branch("drPrevr_next",&drPrevr_next,"drPrevr_next/D");
  data->Branch("dt_prev_us_next",&dt_prev_us_next,"dt_prev_us_next/D");
  data->Branch("subid_next",&subid_next,"subid_next/I");
  data->Branch("mc_energy_next",&mc_energy_next,"mc_energy_next/D");
  //TTree *data_prev = new TTree("data_prev","low-energy detector triggered events");
  double n100_previous,n100_prev_prev;
  double closestPMT_prev;
  double drPrevr_prev,dt_prev_us_prev;
  double mc_energy_prev;
  int subid_prev;
  data->Branch("n100_previous",&n100_previous,"n100_previous/D");
  data->Branch("n100_prev_prev",&n100_prev_prev,"n100_prev_prev/D");
  data->Branch("closestPMT_prev",&closestPMT_prev,"closestPMT_prev/D");
  data->Branch("drPrevr_prev",&drPrevr_prev,"drPrevr_prev/D");
  data->Branch("dt_prev_us_prev",&dt_prev_us_prev,"dt_prev_us_prev/D");
  data->Branch("subid_prev",&subid_prev,"subid_prev/I");
  data->Branch("mc_energy_prev",&mc_energy_prev,"mc_energy_prev/D");
  data->Branch("r_like",&r_like,"r_like/D");
  
  for(int i=0; i<nentries; i++){
    t->GetEntry(i);
    if (i%100000==0){
      printf("Creating likelihoods: Event %d of %d\n",i,nentries);
    }
    if (t->GetLeaf("n100")->GetValue(0) > 0 and t->GetLeaf("closestPMT")->GetValue(0) > -499) {
      
      double nX_sig_bin = nX_signal->GetXaxis()->FindBin(t->GetLeaf(nx)->GetValue(0));
      double nX_sig_prob = nX_signal->GetBinContent(nX_sig_bin);
      double nX_bg_bin = nX_background->GetXaxis()->FindBin(t->GetLeaf(nx)->GetValue(0));
      double nX_bg_prob = nX_background->GetBinContent(nX_bg_bin);
      
      double nX_prev_sig_bin = nX_prev_signal->GetXaxis()->FindBin(t->GetLeaf(nx_prev)->GetValue(0));
      double nX_prev_sig_prob = nX_prev_signal->GetBinContent(nX_prev_sig_bin);
      double nX_prev_bg_bin = nX_prev_background->GetXaxis()->FindBin(t->GetLeaf(nx_prev)->GetValue(0));
      double nX_prev_bg_prob = nX_prev_background->GetBinContent(nX_prev_bg_bin);
      
      double dt_prev_us_sig_bin = dt_prev_us_signal->GetXaxis()->FindBin(t->GetLeaf("dt_prev_us")->GetValue(0));
      double dt_prev_us_sig_prob = dt_prev_us_signal->GetBinContent(dt_prev_us_sig_bin);
      double dt_prev_us_bg_bin = dt_prev_us_background->GetXaxis()->FindBin(t->GetLeaf("dt_prev_us")->GetValue(0));
      double dt_prev_us_bg_prob = dt_prev_us_background->GetBinContent(dt_prev_us_bg_bin);
      
      double drPrevr_sig_bin = drPrevr_signal->GetXaxis()->FindBin(t->GetLeaf("drPrevr")->GetValue(0));
      double drPrevr_sig_prob = drPrevr_signal->GetBinContent(drPrevr_sig_bin);
      double drPrevr_bg_bin = drPrevr_background->GetXaxis()->FindBin(t->GetLeaf("drPrevr")->GetValue(0));
      double drPrevr_bg_prob = drPrevr_background->GetBinContent(drPrevr_bg_bin);
      
      double closestPMT_sig_bin = closestPMT_signal->GetXaxis()->FindBin(t->GetLeaf("closestPMT")->GetValue(0));
      double closestPMT_sig_prob = closestPMT_signal->GetBinContent(closestPMT_sig_bin);
      double closestPMT_bg_bin = closestPMT_background->GetXaxis()->FindBin(t->GetLeaf("closestPMT")->GetValue(0));
      double closestPMT_bg_prob = closestPMT_background->GetBinContent(closestPMT_bg_bin);
      
      double sig_like = log(nX_sig_prob*nX_prev_sig_prob*dt_prev_us_sig_prob*drPrevr_sig_prob*closestPMT_sig_prob);
      double bg_like = log(nX_bg_prob*nX_prev_bg_prob*dt_prev_us_bg_prob*drPrevr_bg_prob*closestPMT_bg_prob);
      {if( std::isinf(sig_like) == true){sig_like=0;}else{}}
      {if( std::isinf(bg_like) == true){bg_like=0;}else{}}
      /*double*/ r_like = sig_like-bg_like;
      signal_like->Fill(sig_like);
      background_like->Fill(bg_like);
      ratio_like->Fill(r_like);
      if(r_like>25 and r_like<45){
        printf("Event: r_like = %f, nX = %f, nX_prev = %f, dt_prev_us = %f, drPrevr = %f, closestPMT = %f, subid = %f\n",r_like,t->GetLeaf(nx)->GetValue(0),t->GetLeaf(nx_prev)->GetValue(0),t->GetLeaf("dt_prev_us")->GetValue(0),t->GetLeaf("drPrevr")->GetValue(0),t->GetLeaf("closestPMT")->GetValue(0),t->GetLeaf("subid")->GetValue(0));
        n100 = t->GetLeaf("n100")->GetValue(0);
        n100_prev = t->GetLeaf("n100_prev")->GetValue(0);
        closestPMT = t->GetLeaf("closestPMT")->GetValue(0);
        drPrevr = t->GetLeaf("drPrevr")->GetValue(0);
        dt_prev_us = t->GetLeaf("dt_prev_us")->GetValue(0);
        subid = t->GetLeaf("subid")->GetValue(0);
        mc_energy = t->GetLeaf("mc_energy")->GetValue(0);
        //data->Fill();
        t->GetEntry(i+1);
        printf("Next event: nX = %f, nX_prev = %f, dt_prev_us = %f, drPrevr = %f, closestPMT = %f, subid = %f\n",t->GetLeaf(nx)->GetValue(0),t->GetLeaf(nx_prev)->GetValue(0),t->GetLeaf("dt_prev_us")->GetValue(0),t->GetLeaf("drPrevr")->GetValue(0),t->GetLeaf("closestPMT")->GetValue(0),t->GetLeaf("subid")->GetValue(0));
        n100_next = t->GetLeaf("n100")->GetValue(0);
        n100_prev_next = t->GetLeaf("n100_prev")->GetValue(0);
        closestPMT_next = t->GetLeaf("closestPMT")->GetValue(0);
        drPrevr_next = t->GetLeaf("drPrevr")->GetValue(0);
        dt_prev_us_next = t->GetLeaf("dt_prev_us")->GetValue(0);
        subid_next = t->GetLeaf("subid")->GetValue(0);
        mc_energy_next = t->GetLeaf("mc_energy")->GetValue(0);
        //data_next->Fill();
        t->GetEntry(i-1);
        printf("Previous event: nX = %f, nX_prev = %f, dt_prev_us = %f, drPrevr = %f, closestPMT = %f, subid = %f\n\n",t->GetLeaf(nx)->GetValue(0),t->GetLeaf(nx_prev)->GetValue(0),t->GetLeaf("dt_prev_us")->GetValue(0),t->GetLeaf("drPrevr")->GetValue(0),t->GetLeaf("closestPMT")->GetValue(0),t->GetLeaf("subid")->GetValue(0));
        n100_previous = t->GetLeaf("n100")->GetValue(0);
        n100_prev_prev = t->GetLeaf("n100_prev")->GetValue(0);
        closestPMT_prev = t->GetLeaf("closestPMT")->GetValue(0);
        drPrevr_prev = t->GetLeaf("drPrevr")->GetValue(0);
        dt_prev_us_prev = t->GetLeaf("dt_prev_us")->GetValue(0);
        subid_prev = t->GetLeaf("subid")->GetValue(0);
        mc_energy_prev = t->GetLeaf("mc_energy")->GetValue(0);
        //data_prev->Fill();
        data->Fill();
      }
      // if(r_like>45){
      //   printf("Event: r_like = %f, nX = %f, nX_prev = %f, dt_prev_us = %f, drPrevr = %f, closestPMT = %f, subid = %f\n",r_like,t->GetLeaf(nx)->GetValue(0),t->GetLeaf(nx_prev)->GetValue(0),t->GetLeaf("dt_prev_us")->GetValue(0),t->GetLeaf("drPrevr")->GetValue(0),t->GetLeaf("closestPMT")->GetValue(0),t->GetLeaf("subid")->GetValue(0));
      //   n100 = t->GetLeaf("n100")->GetValue(0);
      //   n100_prev = t->GetLeaf("n100_prev")->GetValue(0);
      //   closestPMT = t->GetLeaf("closestPMT")->GetValue(0);
      //   drPrevr = t->GetLeaf("drPrevr")->GetValue(0);
      //   dt_prev_us = t->GetLeaf("dt_prev_us")->GetValue(0);
      //   subid = t->GetLeaf("subid")->GetValue(0);
      //   mc_energy = t->GetLeaf("mc_energy")->GetValue(0);
      //   //data->Fill();
      //   t->GetEntry(i+1);
      //   printf("Next event: nX = %f, nX_prev = %f, dt_prev_us = %f, drPrevr = %f, closestPMT = %f, subid = %f\n",t->GetLeaf(nx)->GetValue(0),t->GetLeaf(nx_prev)->GetValue(0),t->GetLeaf("dt_prev_us")->GetValue(0),t->GetLeaf("drPrevr")->GetValue(0),t->GetLeaf("closestPMT")->GetValue(0),t->GetLeaf("subid")->GetValue(0));
      //   n100_next = t->GetLeaf("n100")->GetValue(0);
      //   n100_prev_next = t->GetLeaf("n100_prev")->GetValue(0);
      //   closestPMT_next = t->GetLeaf("closestPMT")->GetValue(0);
      //   drPrevr_next = t->GetLeaf("drPrevr")->GetValue(0);
      //   dt_prev_us_next = t->GetLeaf("dt_prev_us")->GetValue(0);
      //   subid_next = t->GetLeaf("subid")->GetValue(0);
      //   mc_energy_next = t->GetLeaf("mc_energy")->GetValue(0);
      //   //data_next->Fill();
      //   t->GetEntry(i-1);
      //   printf("Previous event: nX = %f, nX_prev = %f, dt_prev_us = %f, drPrevr = %f, closestPMT = %f, subid = %f\n\n",t->GetLeaf(nx)->GetValue(0),t->GetLeaf(nx_prev)->GetValue(0),t->GetLeaf("dt_prev_us")->GetValue(0),t->GetLeaf("drPrevr")->GetValue(0),t->GetLeaf("closestPMT")->GetValue(0),t->GetLeaf("subid")->GetValue(0));
      //   n100_previous = t->GetLeaf("n100")->GetValue(0);
      //   n100_prev_prev = t->GetLeaf("n100_prev")->GetValue(0);
      //   closestPMT_prev = t->GetLeaf("closestPMT")->GetValue(0);
      //   drPrevr_prev = t->GetLeaf("drPrevr")->GetValue(0);
      //   dt_prev_us_prev = t->GetLeaf("dt_prev_us")->GetValue(0);
      //   subid_prev = t->GetLeaf("subid")->GetValue(0);
      //   mc_energy_prev = t->GetLeaf("mc_energy")->GetValue(0);
      //   //data_prev->Fill();
      //   data->Fill();
      // }
      // if(r_like<25){
      //   printf("Event: r_like = %f, nX = %f, nX_prev = %f, dt_prev_us = %f, drPrevr = %f, closestPMT = %f, subid = %f\n",r_like,t->GetLeaf(nx)->GetValue(0),t->GetLeaf(nx_prev)->GetValue(0),t->GetLeaf("dt_prev_us")->GetValue(0),t->GetLeaf("drPrevr")->GetValue(0),t->GetLeaf("closestPMT")->GetValue(0),t->GetLeaf("subid")->GetValue(0));
      //   n100 = t->GetLeaf("n100")->GetValue(0);
      //   n100_prev = t->GetLeaf("n100_prev")->GetValue(0);
      //   closestPMT = t->GetLeaf("closestPMT")->GetValue(0);
      //   drPrevr = t->GetLeaf("drPrevr")->GetValue(0);
      //   dt_prev_us = t->GetLeaf("dt_prev_us")->GetValue(0);
      //   subid = t->GetLeaf("subid")->GetValue(0);
      //   mc_energy = t->GetLeaf("mc_energy")->GetValue(0);
      //   //data->Fill();
      //   t->GetEntry(i+1);
      //   printf("Next event: nX = %f, nX_prev = %f, dt_prev_us = %f, drPrevr = %f, closestPMT = %f, subid = %f\n",t->GetLeaf(nx)->GetValue(0),t->GetLeaf(nx_prev)->GetValue(0),t->GetLeaf("dt_prev_us")->GetValue(0),t->GetLeaf("drPrevr")->GetValue(0),t->GetLeaf("closestPMT")->GetValue(0),t->GetLeaf("subid")->GetValue(0));
      //   n100_next = t->GetLeaf("n100")->GetValue(0);
      //   n100_prev_next = t->GetLeaf("n100_prev")->GetValue(0);
      //   closestPMT_next = t->GetLeaf("closestPMT")->GetValue(0);
      //   drPrevr_next = t->GetLeaf("drPrevr")->GetValue(0);
      //   dt_prev_us_next = t->GetLeaf("dt_prev_us")->GetValue(0);
      //   subid_next = t->GetLeaf("subid")->GetValue(0);
      //   mc_energy_next = t->GetLeaf("mc_energy")->GetValue(0);
      //   //data_next->Fill();
      //   t->GetEntry(i-1);
      //   printf("Previous event: nX = %f, nX_prev = %f, dt_prev_us = %f, drPrevr = %f, closestPMT = %f, subid = %f\n\n",t->GetLeaf(nx)->GetValue(0),t->GetLeaf(nx_prev)->GetValue(0),t->GetLeaf("dt_prev_us")->GetValue(0),t->GetLeaf("drPrevr")->GetValue(0),t->GetLeaf("closestPMT")->GetValue(0),t->GetLeaf("subid")->GetValue(0));
      //   n100_previous = t->GetLeaf("n100")->GetValue(0);
      //   n100_prev_prev = t->GetLeaf("n100_prev")->GetValue(0);
      //   closestPMT_prev = t->GetLeaf("closestPMT")->GetValue(0);
      //   drPrevr_prev = t->GetLeaf("drPrevr")->GetValue(0);
      //   dt_prev_us_prev = t->GetLeaf("dt_prev_us")->GetValue(0);
      //   subid_prev = t->GetLeaf("subid")->GetValue(0);
      //   mc_energy_prev = t->GetLeaf("mc_energy")->GetValue(0);
      //   //data_prev->Fill();
      //   data->Fill();
      // }
    }
  }
  
  // ratio_like->Scale(1/ratio_like->GetEntries());
  // signal_like->Scale(1/signal_like->GetEntries());
  // background_like->Scale(1/background_like->GetEntries());
  
  out->cd();
  signal_like->Write();
  background_like->Write();
  ratio_like->Write();
  out->Close();
  f->Close();
  output->cd();
  data->Write();
  // data_prev->Write();
  // data_next->Write();
  output->Close();

}
