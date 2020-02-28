//root -l plotEfficiency.C()

#include "TLegend.h"
#include "TLatex.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TTree.h"
#include "TChain.h"
#include <vector>
#include <fstream>
#include <string>
#include "TROOT.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>


#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooArgusBG.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "TText.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"


using namespace RooFit;

float getXmax(TH1F* histo, float& YMax){

  float yVal = 0.;
  int xBin = 1;
  for(int iB=1; iB<histo->GetNbinsX(); ++iB){
    if(histo->GetBinContent(iB) > yVal){
      xBin = iB;
      yVal = histo->GetBinContent(iB);
      YMax = yVal;
      if(yVal > 0 && histo->GetBinContent(iB) < yVal) break;
    }
  }

  std::cout << histo->GetName() << " " << histo->GetBinCenter(xBin) << std::endl; 
  return histo->GetBinCenter(xBin);
}


void plotEfficiency(){
  gROOT->Macro("/afs/cern.ch/user/a/amartell/public/setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  gStyle->SetOptStat(0);
  //  gStyle->SetOptFit(1);

  //  gROOT->SetBatch(kTRUE);

  std::cout << " inizio ci sono " << std::endl; 


  int iColors[4] = {kRed+1, kBlue, kGreen+1, kBlack};
  int iStyle[4] = {20, 22, 23, 21};
  std::cout << " >>> now load files " << std::endl;

  TFile* inF = TFile::Open("../test/BPark_HLTEff_valid_noPU.root");
  //TFile* inF = TFile::Open("../test/BPark_HLTEff_valid_withPU.root");
  
  std::cout << " >>> file loaded " << std::endl;

  TH1F* h_ratio_pixTrk_trkP_pu[4];
  TH1F* h_ratio_pixTrk_trkP_eta[4];
  TH1F* h_ratio_pixTrk_trkP_phi[4];
  TH1F* h_ratio_pixTrk_trkP_pt[4];
  TH1F* h_trkP_Bmass = (TH1F*)(inF->Get("myOutputTest/h_trkP_Bmass"));
  TH1F* h_pixTrk_Bmass = (TH1F*)(inF->Get("myOutputTest/h_pixTrk_Bmass"));
  TH1F* h_ratio_pixTrk_trkP_Bmass = (TH1F*)(inF->Get("myOutputTest/h_ratio_pixTrk_trkP_Bmass"));

  for(int ij=0; ij<4; ++ij){
    h_ratio_pixTrk_trkP_pu[ij] = (TH1F*)(inF->Get(Form("myOutputTest/h_ratio_pixTrk_trkP_pu_%d", ij)));
    h_ratio_pixTrk_trkP_eta[ij] = (TH1F*)(inF->Get(Form("myOutputTest/h_ratio_pixTrk_trkP_eta_%d", ij)));
    h_ratio_pixTrk_trkP_phi[ij] = (TH1F*)(inF->Get(Form("myOutputTest/h_ratio_pixTrk_trkP_phi_%d", ij)));
    h_ratio_pixTrk_trkP_pt[ij] = (TH1F*)(inF->Get(Form("myOutputTest/h_ratio_pixTrk_trkP_pt_%d", ij)));   

    h_ratio_pixTrk_trkP_pu[ij]->SetLineColor(iColors[ij]);
    h_ratio_pixTrk_trkP_eta[ij]->SetLineColor(iColors[ij]);
    h_ratio_pixTrk_trkP_phi[ij]->SetLineColor(iColors[ij]);
    h_ratio_pixTrk_trkP_pt[ij]->SetLineColor(iColors[ij]);

    h_ratio_pixTrk_trkP_pu[ij]->SetLineWidth(2);
    h_ratio_pixTrk_trkP_eta[ij]->SetLineWidth(2);
    h_ratio_pixTrk_trkP_phi[ij]->SetLineWidth(2);
    h_ratio_pixTrk_trkP_pt[ij]->SetLineWidth(2);

    h_ratio_pixTrk_trkP_pu[ij]->SetMarkerColor(iColors[ij]);
    h_ratio_pixTrk_trkP_eta[ij]->SetMarkerColor(iColors[ij]);
    h_ratio_pixTrk_trkP_phi[ij]->SetMarkerColor(iColors[ij]);
    h_ratio_pixTrk_trkP_pt[ij]->SetMarkerColor(iColors[ij]);

    h_ratio_pixTrk_trkP_pu[ij]->SetMarkerStyle(iStyle[ij]);
    h_ratio_pixTrk_trkP_eta[ij]->SetMarkerStyle(iStyle[ij]);
    h_ratio_pixTrk_trkP_phi[ij]->SetMarkerStyle(iStyle[ij]);
    h_ratio_pixTrk_trkP_pt[ij]->SetMarkerStyle(iStyle[ij]);
  }
  h_ratio_pixTrk_trkP_Bmass->SetLineWidth(2);
  h_trkP_Bmass->SetLineWidth(2);
  h_pixTrk_Bmass->SetLineWidth(2);

  h_ratio_pixTrk_trkP_Bmass->SetLineColor(kBlack);
  h_trkP_Bmass->SetLineColor(kBlack);
  h_pixTrk_Bmass->SetLineColor(kBlack);

  h_ratio_pixTrk_trkP_Bmass->SetMarkerStyle(iStyle[3]);
  h_trkP_Bmass->SetMarkerStyle(iStyle[3]);
  h_pixTrk_Bmass->SetMarkerStyle(iStyle[3]);

  h_ratio_pixTrk_trkP_Bmass->SetMarkerColor(kBlack);
  h_trkP_Bmass->SetMarkerColor(kBlack);
  h_pixTrk_Bmass->SetMarkerColor(kBlack);

  std::vector<std::string> partName;
  partName.push_back("ele1");
  partName.push_back("ele2");
  partName.push_back("kaon");
  partName.push_back("Bcand");

  TLegend *legTGM = new TLegend(0.65,0.80,0.85,0.95,NULL,"brNDC");
  legTGM->SetFillStyle(0);
  legTGM->SetTextFont(42);
  legTGM->SetTextSize(0.04);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  legTGM->SetNColumns(2);
  for(int iF=0; iF<partName.size(); ++iF){
    legTGM->AddEntry(h_ratio_pixTrk_trkP_pt[iF], partName.at(iF).c_str(), "lp");
  }



  TCanvas* ch_PU = new TCanvas();
  ch_PU->cd();
  h_ratio_pixTrk_trkP_pu[0]->GetXaxis()->SetTitle("nPU");
  h_ratio_pixTrk_trkP_pu[0]->GetYaxis()->SetTitle("eff (pixTracks / gen)");
  h_ratio_pixTrk_trkP_pu[0]->GetYaxis()->SetRangeUser(0., 1.2);
  h_ratio_pixTrk_trkP_pu[0]->GetXaxis()->SetRangeUser(50., 80.);
  h_ratio_pixTrk_trkP_pu[0]->Draw("pe");
  for(int iF=1; iF<partName.size(); ++iF){
    h_ratio_pixTrk_trkP_pu[iF]->Draw("pe, same");
  }
  legTGM->Draw("same");
  ch_PU->Print("plots/h_ratio_pixTrk_trkP_PU.png", "png");
  ch_PU->Print("plots/h_ratio_pixTrk_trkP_PU.pdf", "pdf");
  ch_PU->Print("plots/h_ratio_pixTrk_trkP_PU.root", "root");

  //
  TCanvas* ch_eta = new TCanvas();
  ch_eta->cd();
  h_ratio_pixTrk_trkP_eta[0]->GetXaxis()->SetTitle("gen #eta");
  h_ratio_pixTrk_trkP_eta[0]->GetYaxis()->SetTitle("eff (pixTracks / gen)");
  h_ratio_pixTrk_trkP_eta[0]->GetYaxis()->SetRangeUser(0., 1.2);
  h_ratio_pixTrk_trkP_eta[0]->Draw("pe");
  for(int iF=1; iF<partName.size(); ++iF){
    h_ratio_pixTrk_trkP_eta[iF]->Draw("pe, same");
  }
  legTGM->Draw("same");
  ch_eta->Print("plots/h_ratio_pixTrk_trkP_eta.png", "png");
  ch_eta->Print("plots/h_ratio_pixTrk_trkP_eta.pdf", "pdf");
  ch_eta->Print("plots/h_ratio_pixTrk_trkP_eta.root", "root");

  //
  TCanvas* ch_phi = new TCanvas();
  ch_phi->cd();
  h_ratio_pixTrk_trkP_phi[0]->GetXaxis()->SetTitle("gen #phi");
  h_ratio_pixTrk_trkP_phi[0]->GetYaxis()->SetTitle("eff (pixTracks / gen)");
  h_ratio_pixTrk_trkP_phi[0]->GetYaxis()->SetRangeUser(0., 1.2);
  h_ratio_pixTrk_trkP_phi[0]->Draw("pe");
  for(int iF=1; iF<partName.size(); ++iF){
    h_ratio_pixTrk_trkP_phi[iF]->Draw("pe, same");
  }
  legTGM->Draw("same");
  ch_phi->Print("plots/h_ratio_pixTrk_trkP_phi.png", "png");
  ch_phi->Print("plots/h_ratio_pixTrk_trkP_phi.pdf", "pdf");
  ch_phi->Print("plots/h_ratio_pixTrk_trkP_phi.root", "root");

  //
  TCanvas* ch_pt = new TCanvas();
  ch_pt->cd();
  h_ratio_pixTrk_trkP_pt[0]->GetXaxis()->SetTitle("gen pt");
  h_ratio_pixTrk_trkP_pt[0]->GetYaxis()->SetTitle("eff (pixTracks / gen)");
  h_ratio_pixTrk_trkP_pt[0]->GetYaxis()->SetRangeUser(0., 1.2);
  h_ratio_pixTrk_trkP_pt[0]->GetXaxis()->SetRangeUser(0., 10);
  h_ratio_pixTrk_trkP_pt[0]->Draw("pe");
  for(int iF=1; iF<partName.size(); ++iF){
    h_ratio_pixTrk_trkP_pt[iF]->Draw("pe, same");
  }
  legTGM->Draw("same");
  ch_pt->Print("plots/h_ratio_pixTrk_trkP_pt.png", "png");
  ch_pt->Print("plots/h_ratio_pixTrk_trkP_pt.pdf", "pdf");
  ch_pt->Print("plots/h_ratio_pixTrk_trkP_pt.root", "root");


  TCanvas* ch_Bmass_eff = new TCanvas();
  ch_Bmass_eff->cd();
  h_ratio_pixTrk_trkP_Bmass->GetXaxis()->SetTitle("gen Bmass");
  h_ratio_pixTrk_trkP_Bmass->GetYaxis()->SetTitle("eff (pixTracks / gen)");
  h_ratio_pixTrk_trkP_Bmass->GetYaxis()->SetRangeUser(0., 1.2);
  h_ratio_pixTrk_trkP_Bmass->Draw("pe");
  ch_Bmass_eff->Print("plots/h_ratio_pixTrk_trkP_Bmass_eff.png", "png");
  ch_Bmass_eff->Print("plots/h_ratio_pixTrk_trkP_Bmass_eff.pdf", "pdf");
  ch_Bmass_eff->Print("plots/h_ratio_pixTrk_trkP_Bmass_eff.root", "root");


  TCanvas* ch_Bmass_gen = new TCanvas();
  ch_Bmass_gen->cd();
  h_trkP_Bmass->GetXaxis()->SetTitle("gen Bmass");
  h_trkP_Bmass->GetYaxis()->SetTitle("n.Events");
  h_trkP_Bmass->Draw("h");
  ch_Bmass_gen->Print("plots/h_ratio_pixTrk_trkP_Bmass_gen.png", "png");
  ch_Bmass_gen->Print("plots/h_ratio_pixTrk_trkP_Bmass_gen.pdf", "pdf");
  ch_Bmass_gen->Print("plots/h_ratio_pixTrk_trkP_Bmass_gen.root", "root");


  TCanvas* ch_Bmass_pixTracks = new TCanvas();
  ch_Bmass_pixTracks->cd();
  h_pixTrk_Bmass->GetXaxis()->SetTitle("pixTracks Bmass");
  h_pixTrk_Bmass->GetYaxis()->SetTitle("n.Events");
  h_pixTrk_Bmass->Draw("h");
  ch_Bmass_pixTracks->Print("plots/h_ratio_pixTrk_trkP_Bmass_pixTracks.png", "png");
  ch_Bmass_pixTracks->Print("plots/h_ratio_pixTrk_trkP_Bmass_pixTracks.pdf", "pdf");
  ch_Bmass_pixTracks->Print("plots/h_ratio_pixTrk_trkP_Bmass_pixTracks.root", "root");



  return;
}
