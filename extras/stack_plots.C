#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <TROOT.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TColor.h>
#include <math.h>


void stack_plots(){

  //gROOT->SetBatch(kTRUE);

  TFile *world_2028 = new TFile("world_2028_classified_energy.root");
  TFile *hink = new TFile("hinkley_c_classified_energy.root");
  TFile *size = new TFile("sizewell_b_classified_energy.root");
  TFile *grav = new TFile("gravelines_classified_energy.root");
  TFile *hey = new TFile("heysham_2_classified_energy.root");
  TFile *tor = new TFile("torness_full_classified_energy.root");
  TFile *li9 = new TFile("li9_classified_energy.root");
  TFile *n17 = new TFile("n17_classified_energy.root");
  TFile *geo = new TFile("geo_classified_energy.root");
  
  TTree *t_w_2028 = (TTree*)world_2028->Get("data");
  TTree *t_hink = (TTree*)hink->Get("data");
  TTree *t_size = (TTree*)size->Get("data");
  TTree *t_grav = (TTree*)grav->Get("data");
  TTree *t_hey = (TTree*)hey->Get("data");
  TTree *t_tor = (TTree*)tor->Get("data");
  TTree *t_li9 = (TTree*)li9->Get("data");
  TTree *t_n17 = (TTree*)n17->Get("data");
  TTree *t_geo = (TTree*)geo->Get("data");

  TTree *run_w_2028 = (TTree*)world_2028->Get("runSummary");
  TTree *run_hink = (TTree*)hink->Get("runSummary");
  TTree *run_size = (TTree*)size->Get("runSummary");
  TTree *run_grav = (TTree*)grav->Get("runSummary");
  TTree *run_hey = (TTree*)hey->Get("runSummary");
  TTree *run_tor = (TTree*)tor->Get("runSummary");
  TTree *run_li9 = (TTree*)li9->Get("runSummary");
  TTree *run_n17 = (TTree*)n17->Get("runSummary");
  TTree *run_geo = (TTree*)geo->Get("runSummary");

  run_w_2028->GetEntry(0);
  int MCevents_w_2028 = 2*run_w_2028->GetLeaf("nevents")->GetValue();
  t_w_2028->GetEntry(0);
  int events_w_2028 = t_w_2028->GetEntries();

  run_hink->GetEntry(0);
  int MCevents_hink = 2*run_hink->GetLeaf("nevents")->GetValue();
  t_hink->GetEntry(0);
  int events_hink = t_hink->GetEntries();

  run_size->GetEntry(0);
  int MCevents_size = 2*run_size->GetLeaf("nevents")->GetValue();
  t_size->GetEntry(0);
  int events_size = t_size->GetEntries();

  run_grav->GetEntry(0);
  int MCevents_grav = 2*run_grav->GetLeaf("nevents")->GetValue();
  t_grav->GetEntry(0);
  int events_grav = t_grav->GetEntries();

  run_hey->GetEntry(0);
  int MCevents_hey = 2*run_hey->GetLeaf("nevents")->GetValue();
  t_hey->GetEntry(0);
  int events_hey = t_hey->GetEntries();

  run_tor->GetEntry(0);
  int MCevents_tor = 2*run_tor->GetLeaf("nevents")->GetValue();
  t_tor->GetEntry(0);
  int events_tor = t_tor->GetEntries();

  run_li9->GetEntry(0);
  int MCevents_li9 = 2*run_li9->GetLeaf("nevents")->GetValue();
  t_li9->GetEntry(0);
  int events_li9 = t_li9->GetEntries();

  run_n17->GetEntry(0);
  int MCevents_n17 = 2*run_n17->GetLeaf("nevents")->GetValue();
  t_n17->GetEntry(0);
  int events_n17 = t_n17->GetEntries();

  run_geo->GetEntry(0);
  int MCevents_geo = 2*run_geo->GetLeaf("nevents")->GetValue();
  t_geo->GetEntry(0);
  int events_geo = t_geo->GetEntries();

  t_w_2028->Draw("mc_energy>>hist_w_2028","subid==0");
  TH1 *hist_w_2028 = (TH1*)gDirectory->Get("hist_w_2028");

  t_hink->Draw("mc_energy>>hist_hink","subid==0");
  TH1 *hist_hink = (TH1*)gDirectory->Get("hist_hink");

  t_size->Draw("mc_energy>>hist_size","subid==0");
  TH1 *hist_size = (TH1*)gDirectory->Get("hist_size");

  t_grav->Draw("mc_energy>>hist_grav","subid==0");
  TH1 *hist_grav = (TH1*)gDirectory->Get("hist_grav");

  t_hey->Draw("mc_energy>>hist_hey","subid==0");
  TH1 *hist_hey = (TH1*)gDirectory->Get("hist_hey");

  t_tor->Draw("mc_energy>>hist_tor","subid==0");
  TH1 *hist_tor = (TH1*)gDirectory->Get("hist_tor");

  t_li9->Draw("mc_energy>>hist_li9","subid==0");
  TH1 *hist_li9 = (TH1*)gDirectory->Get("hist_li9");

  t_n17->Draw("mc_energy>>hist_n17","subid==0");
  TH1 *hist_n17 = (TH1*)gDirectory->Get("hist_n17");

  t_geo->Draw("mc_energy>>hist_geo","subid==0");
  TH1 *hist_geo = (TH1*)gDirectory->Get("hist_geo");

  hist_w_2028->Scale((1/hist_w_2028->GetEntries())*86400*0.00001474*events_w_2028/MCevents_w_2028);
  hist_hink->Scale((1/hist_hink->GetEntries())*86400*0.000002563*events_hink/MCevents_hink);
  hist_size->Scale((1/hist_size->GetEntries())*86400*0.00000133*events_size/MCevents_size);
  hist_grav->Scale((1/hist_grav->GetEntries())*86400*0.000002599*events_grav/MCevents_grav);
  hist_hey->Scale((1/hist_hey->GetEntries())*86400*0.000006894*events_hey/MCevents_hey);
  hist_tor->Scale((1/hist_tor->GetEntries())*86400*0.000004419*events_tor/MCevents_tor);
  hist_li9->Scale((1/hist_li9->GetEntries())*86400*0.0000325*0.069*events_li9/MCevents_li9);
  hist_n17->Scale((1/hist_n17->GetEntries())*86400*0.0000199*0.85*events_n17/MCevents_n17);
  hist_geo->Scale((1/hist_geo->GetEntries())*86400*0.000006364*events_geo/MCevents_geo);

  // hist_w_2028->Scale(MCevents_w_2028/(86400*0.00001474*events_w_2028));
  // hist_hink->Scale(MCevents_hink/(86400*0.000002563*events_hink));
  // hist_size->Scale(MCevents_size/(86400*0.00000133*events_size));
  // hist_grav->Scale(MCevents_grav/(86400*0.000002599*events_grav));
  // hist_hey->Scale(MCevents_hink/(86400*0.000006894*events_hey));
  // hist_tor->Scale(MCevents_size/(86400*0.000004419*events_tor));
  // hist_li9->Scale(MCevents_grav/(86400*0.0000325*0.069*events_li9));
  // hist_n17->Scale(MCevents_hink/(86400*0.0000199*0.85*events_n17));
  // hist_geo->Scale(MCevents_size/(86400*0.000006364*events_geo));

  // hist_w_2028->Scale(1/hist_w_2028->GetEntries());
  // hist_hink->Scale(1/hist_hink->GetEntries());
  // hist_size->Scale(1/hist_size->GetEntries());
  // hist_grav->Scale(1/hist_grav->GetEntries());
  // hist_hey->Scale(1/hist_hey->GetEntries());
  // hist_tor->Scale(1/hist_tor->GetEntries());
  // hist_li9->Scale(1/hist_li9->GetEntries());
  // hist_n17->Scale(1/hist_n17->GetEntries());
  // hist_geo->Scale(1/hist_geo->GetEntries());

  std::cout << "World: " << hist_w_2028->Integral() << std::endl;
  std::cout << "Hinkley Point C: " << hist_hink->Integral() << std::endl;
  std::cout << "Sizewell B: " << hist_size->Integral() << std::endl;
  std::cout << "Gravelines: " << hist_grav->Integral() << std::endl;
  std::cout << "Heysham 2: " << hist_hey->Integral()  << std::endl;
  std::cout << "Torness: " << hist_tor->Integral() << std::endl;
  std::cout << "Lithium 9: " << hist_li9->Integral() << std::endl;
  std::cout << "Nitrogen 17: " << hist_n17->Integral() << std::endl;
  std::cout << "Geoneutrinos: " << hist_geo->Integral() << std::endl;


  THStack * stack = new THStack("stack","Spectra");
  TLegend *legend = new TLegend(0.85,0.5,0.98,0.98);
  
 
  hist_w_2028->SetLineColor(1);
  hist_hink->SetLineColor(2);
  hist_size->SetLineColor(3);
  hist_grav->SetLineColor(4);
  hist_hey->SetLineColor(5);
  hist_tor->SetLineColor(6);
  hist_li9->SetLineColor(7);
  hist_n17->SetLineColor(8);
  hist_geo->SetLineColor(9);
  

  hist_w_2028->SetLineWidth(2);
  hist_hink->SetLineWidth(2);
  hist_size->SetLineWidth(2);
  hist_grav->SetLineWidth(2);
  hist_hey->SetLineWidth(2);
  hist_tor->SetLineWidth(2);
  hist_li9->SetLineWidth(2);
  hist_n17->SetLineWidth(2);
  hist_geo->SetLineWidth(2);
  
  legend->AddEntry(hist_w_2028,"World 2028");
  legend->AddEntry(hist_hink,"Hinkley C");
  legend->AddEntry(hist_size,"Sizewell B");
  legend->AddEntry(hist_grav,"Gravelines");
  legend->AddEntry(hist_hey,"Heysham 2");
  legend->AddEntry(hist_tor,"Torness");
  legend->AddEntry(hist_li9,"Li^{9}");
  legend->AddEntry(hist_n17,"N^{17}");
  legend->AddEntry(hist_geo,"Geo");
  
  stack->Add(hist_w_2028);
  stack->Add(hist_hink);
  stack->Add(hist_size);
  stack->Add(hist_grav);
  stack->Add(hist_hey);
  stack->Add(hist_tor);
  stack->Add(hist_li9);
  stack->Add(hist_n17);
  stack->Add(hist_geo);
  
  //stack->SetStats(0);

  TCanvas * c1 = new TCanvas("c1");
  // hist_w_2021->Draw("HIST");
  // hist_w_2028->Draw("HIST && same");
  stack->Draw("nostack hist");
  stack->GetXaxis()->SetTitle("E_{MC} [MeV]");
  //hist_w_2021->Draw("HIST && same");
  //stack->GetYaxis()->SetTitle("Normalised events");
  //stack->GetYaxis()->SetTitleOffset(1.3);
  legend->Draw();
  // c1->SaveAs("world_spectra.pdf");


}
