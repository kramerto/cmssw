// C/C++ standard libraries
#include <iostream> // for input/output prints
#include <fstream>  // for reading the ASCII files

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TTreeReader.h"

TFile *file;

TTree *tree;

// dE/dx_vs_p
TH2D *dEdx_vs_p;
TH2D *dEdx_vs_p_211;
TH2D *dEdx_vs_p_11;
TH2D *dEdx_vs_p_13;
TH2D *dEdx_vs_p_321;
TH2D *dEdx_vs_p_2212;
TCanvas *dEdx_vs_pCanvas;
TLegend *dEdx_vs_pLegend;

// lostHits
TH1D *lostHits;
TH1D *lostHits_dEdx0;
THStack *lostHitsStack = new THStack("LostHitsStack", "");
TCanvas *lostHitsCanvas;
TLegend *lostHitsLegend;

// chi2OverNdof
TH1D *chi2OverNdof;
TH1D *chi2OverNdof_unmatched;
THStack *chi2OverNdofStack = new THStack("chi2OverNdofStack", "");
TCanvas *chi2OverNdofCanvas;
TLegend *chi2OverNdofLegend;

// p_vs_N
TH2D *p_vs_N;
TGraph *p_vs_NTGraph;
TCanvas *p_vs_NCanvas;

// dEdx_vs_N
TH2D *dEdx_vs_N;
TCanvas *dEdx_vs_NCanvas;

// particleID
TH1D *particleID;
TCanvas *particleIDCanvas;

// deltaR
TH1D *deltaR;
TCanvas *deltaRCanvas;

// dZ
TH1D *dZ;
TH1D *dZ_unmatched;
TCanvas *dZCanvas;
THStack *dZStack = new THStack("dZStack", "");

void histogram_macro() {
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);

  // open root file
  file = new TFile("histodemo.root", "READ");

  // get the trees
  tree = (TTree *)file->Get("demo/Tree");

  // dEdx_vs_p

  // define the histogram
  dEdx_vs_p = new TH2D("dEdx_vs_p", "dEdx_vs_p", 1000, 0, 5, 1000, 0, 10);
  dEdx_vs_p->GetXaxis()->SetTitle("p in GeV/c");
  dEdx_vs_p->GetYaxis()->SetTitle("dE/dx Estimator (dedxHarmonic2) in MeV/cm");
  dEdx_vs_p->GetYaxis()->SetTitleOffset(0.6);
  dEdx_vs_p->GetXaxis()->SetTitleOffset(1.3);

  dEdx_vs_p_11 =
      new TH2D("dEdx_vs_p_11", "dEdx_vs_p_11", 1000, 0, 5, 1000, 0, 10);

  dEdx_vs_p_13 =
      new TH2D("dEdx_vs_p_13", "dEdx_vs_p_13", 1000, 0, 5, 1000, 0, 10);
  dEdx_vs_p_2212 =
      new TH2D("dEdx_vs_p_2212", "dEdx_vs_p_2212", 1000, 0, 5, 1000, 0, 10);

  dEdx_vs_p_321 =
      new TH2D("dEdx_vs_p_321", "dEdx_vs_p_321", 1000, 0, 5, 1000, 0, 10);

  dEdx_vs_p_211 =
      new TH2D("dEdx_vs_p_211", "dEdx_vs_p_211", 1000, 0, 5, 1000, 0, 10);

  // define the canvas
  dEdx_vs_pCanvas =
      new TCanvas("dEdx_vs_pCanvas", "dEdx_vs_pCanvas", 1020, 520);
  dEdx_vs_pCanvas->SetFillColor(kWhite);
  dEdx_vs_pCanvas->SetBorderMode(0);
  dEdx_vs_pCanvas->SetLeftMargin(0.14);
  dEdx_vs_pCanvas->SetRightMargin(0.16);
  dEdx_vs_pCanvas->SetBottomMargin(0.14);

  //  fill the histogram
  TTreeReader reader("demo/Tree", file);
  using doubleVector = std::vector<double>;
  TTreeReaderValue<doubleVector> dEdxRV(reader, "dEdx");
  TTreeReaderValue<doubleVector> pRV(reader, "p");
  TTreeReaderValue<std::vector<int>> pIDRV(reader, "particleID");
  while (reader.Next()) {
    std::vector<double> dEdxVector = *dEdxRV;
    std::vector<double> pVector = *pRV;
    std::vector<int> pIDVector = *pIDRV;
    for (unsigned int i = 0; i < dEdxVector.size(); i++) {
      if (std::abs(pIDVector[i]) == 211) {
        dEdx_vs_p_211->SetMarkerColor(kYellow + 1);
        dEdx_vs_p_211->SetFillColor(kYellow + 1);
        dEdx_vs_p_211->SetMarkerStyle(8);
        dEdx_vs_p_211->SetMarkerSize(0.4);
        dEdx_vs_p_211->Fill(pVector[i], dEdxVector[i]);
      } else if (std::abs(pIDVector[i]) == 321) {
        dEdx_vs_p_321->SetMarkerColor(kGreen + 2);
        dEdx_vs_p_321->SetFillColor(kGreen + 2);
        dEdx_vs_p_321->SetMarkerStyle(8);
        dEdx_vs_p_321->SetMarkerSize(0.4);
        dEdx_vs_p_321->Fill(pVector[i], dEdxVector[i]);
      } else if (std::abs(pIDVector[i]) == 11) {
        dEdx_vs_p_11->SetMarkerColor(kRed);
        dEdx_vs_p_11->SetFillColor(kRed);
        dEdx_vs_p_11->SetMarkerStyle(8);
        dEdx_vs_p_11->SetMarkerSize(0.4);
        dEdx_vs_p_11->Fill(pVector[i], dEdxVector[i]);
      } else if (std::abs(pIDVector[i]) == 2212) {
        dEdx_vs_p_2212->SetMarkerColor(kBlue);
        dEdx_vs_p_2212->SetFillColor(kBlue);
        dEdx_vs_p_2212->SetMarkerStyle(8);
        dEdx_vs_p_2212->SetMarkerSize(0.4);
        dEdx_vs_p_2212->Fill(pVector[i], dEdxVector[i]);
      } else if (std::abs(pIDVector[i]) == 13) {
        dEdx_vs_p_13->SetMarkerColor(kCyan);
        dEdx_vs_p_13->SetFillColor(kCyan);
        dEdx_vs_p_13->SetMarkerStyle(8);
        dEdx_vs_p_13->SetMarkerSize(0.4);
        dEdx_vs_p_13->Fill(pVector[i], dEdxVector[i]);
      } else {
        dEdx_vs_p->SetMarkerColor(kBlack);
        dEdx_vs_p->SetFillColor(kBlack);
        dEdx_vs_p->SetMarkerStyle(8);
        dEdx_vs_p->SetMarkerSize(0.4);
        dEdx_vs_p->Fill(pVector[i], dEdxVector[i]);
      }
    }
  }

  // draw the histogram
  dEdx_vs_p->Draw();
  dEdx_vs_p_211->Draw("same");
  dEdx_vs_p_321->Draw("same");
  dEdx_vs_p_2212->Draw("same");
  dEdx_vs_p_11->Draw("same");
  dEdx_vs_p_13->Draw("same");

  // legend
  dEdx_vs_pLegend = new TLegend(0.65, 0.6, 0.84, 0.9);
  dEdx_vs_pLegend->SetTextSize(0.05);
  dEdx_vs_pLegend->AddEntry(dEdx_vs_p, "unmatched", "f");
  dEdx_vs_pLegend->AddEntry(dEdx_vs_p_211, "#pi^{+/-}", "f");
  dEdx_vs_pLegend->AddEntry(dEdx_vs_p_321, "K^{+/-}", "f");
  dEdx_vs_pLegend->AddEntry(dEdx_vs_p_2212, "p^{+/-}", "f");
  dEdx_vs_pLegend->AddEntry(dEdx_vs_p_11, "e^{+/-}", "f");
  dEdx_vs_pLegend->AddEntry(dEdx_vs_p_13, "#mu^{+/-}", "f");
  dEdx_vs_pLegend->Draw();

  // save the histogram
  dEdx_vs_pCanvas->Print("dEdx_vs_p.eps");

  //  // lostHits
  //
  //  // define the histogram
  //  lostHits = new TH1D("Lost_hits", "Lost_hits", 4, 0, 4);
  //  lostHits->GetYaxis()->SetTitleOffset(0.6);
  //  lostHits->GetXaxis()->SetTitleOffset(1.3);
  //  lostHits->SetLineColor(kBlack);
  //
  //  lostHits_dEdx0 = new TH1D("Lost_hits_dEdx0", "Lost_hits_dEdx0", 4, 0, 4);
  //  lostHits_dEdx0->GetYaxis()->SetTitleOffset(0.6);
  //  lostHits_dEdx0->GetXaxis()->SetTitleOffset(1.3);
  //  lostHits_dEdx0->SetFillColor(kRed);
  //  lostHits_dEdx0->SetLineColor(kBlack);
  //
  //  // define the canvas
  //  lostHitsCanvas = new TCanvas("lostHitsCanvas", "lostHitsCanvas", 1020,
  //  520);
  //  lostHitsCanvas->SetFillColor(0);
  //  lostHitsCanvas->SetBorderMode(0);
  //  lostHitsCanvas->SetLeftMargin(0.14);
  //  lostHitsCanvas->SetRightMargin(0.16);
  //  lostHitsCanvas->SetBottomMargin(0.14);
  //
  //  // draw the histogram
  //  tree->Draw("nLostHits>>Lost_hits", "dEdx>0");
  //  tree->Draw("nLostHits>>Lost_hits_dEdx0", "dEdx==0");
  //  lostHitsStack->Add(lostHits);
  //  lostHitsStack->Add(lostHits_dEdx0);
  //  lostHitsStack->Draw();
  //  lostHitsStack->GetXaxis()->SetTitle("Number of invalid hits per track");
  //  lostHitsStack->GetYaxis()->SetTitle("Multiplicity");
  //  lostHitsStack->GetYaxis()->SetTitleOffset(1.3);
  //  lostHitsStack->GetXaxis()->SetTitleOffset(1.3);
  //
  //  // legend
  //  lostHitsLegend = new TLegend(0.84, 0.7, 0.5, 0.9);
  //  lostHitsLegend->SetTextSize(0.05);
  //  lostHitsLegend->AddEntry(lostHits, "Lost hits for dE/dx > 0", "f");
  //  lostHitsLegend->AddEntry(lostHits_dEdx0, "Lost hits for dE/dx = 0", "f");
  //  lostHitsLegend->Draw();
  //
  //  // save the histogram
  //  lostHitsCanvas->Print("lost_hits.eps");

  // chi2OverNdof

  // define the histogram
  chi2OverNdof = new TH1D("chi2OverNdof", "chi2OverNdof", 100, 0, 4);
  chi2OverNdof->GetYaxis()->SetTitleOffset(0.6);
  chi2OverNdof->GetXaxis()->SetTitleOffset(1.3);
  chi2OverNdof->SetLineColor(kBlack);

  chi2OverNdof_unmatched =
      new TH1D("chi2OverNdof_unmatched", "chi2OverNdof_unmatched", 100, 0, 4);
  chi2OverNdof_unmatched->GetYaxis()->SetTitleOffset(0.6);
  chi2OverNdof_unmatched->GetXaxis()->SetTitleOffset(1.3);
  chi2OverNdof_unmatched->SetLineColor(kRed);

  // define the canvas
  chi2OverNdofCanvas =
      new TCanvas("chi2OverNdofCanvas", "chi2OverNdofCanvas", 1020, 520);
  chi2OverNdofCanvas->SetFillColor(0);
  chi2OverNdofCanvas->SetBorderMode(0);
  chi2OverNdofCanvas->SetLeftMargin(0.14);
  chi2OverNdofCanvas->SetRightMargin(0.16);
  chi2OverNdofCanvas->SetBottomMargin(0.14);

  // draw the histogram
  tree->Draw("chi2/ndof>>chi2OverNdof",
             "(std::abs(particleID)==11)||(std::abs(particleID)==211)||("
             "std::abs(particleID)==321)||(std::abs(particleID)==2212)||(std::"
             "abs(particleID)==13)");
  tree->Draw("chi2/ndof>>chi2OverNdof_unmatched",
             "!((std::abs(particleID)==11)||(std::abs(particleID)==211)||("
             "std::abs(particleID)==321)||(std::abs(particleID)==2212)||(std::"
             "abs(particleID)==13))");
  chi2OverNdofStack->Add(chi2OverNdof);
  chi2OverNdofStack->Add(chi2OverNdof_unmatched);
  chi2OverNdofStack->Draw("nostack");
  chi2OverNdofStack->GetXaxis()->SetTitle("#chi^{2}/n_{dof}");
  chi2OverNdofStack->GetYaxis()->SetTitle("Number of tracks");
  chi2OverNdofStack->GetYaxis()->SetTitleOffset(1.3);
  chi2OverNdofStack->GetXaxis()->SetTitleOffset(1.3);

  // legend
  chi2OverNdofLegend = new TLegend(0.4, 0.7, 0.84, 0.9);
  chi2OverNdofLegend->SetTextSize(0.05);
  chi2OverNdofLegend->AddEntry(chi2OverNdof,
                               "#chi^{2}/n_{dof} for matched particles", "l");
  chi2OverNdofLegend->AddEntry(chi2OverNdof_unmatched,
                               "#chi^{2}/n_{dof} for unmatched particles", "l");
  chi2OverNdofLegend->Draw();

  // save the histogram
  chi2OverNdofCanvas->Print("chi2OverNdof.eps");

  // p_vs_N

  // define the histogram
  p_vs_N = new TH2D("p_vs_N", "p_vs_N", 40, 0, 40, 100, 0, 5);
  p_vs_N->GetYaxis()->SetTitle("p in GeV/c");
  p_vs_N->GetXaxis()->SetTitle("Number of Hits per Track");
  p_vs_N->GetYaxis()->SetTitleOffset(0.6);
  p_vs_N->GetXaxis()->SetTitleOffset(1.3);

  //  // define TGraph
  //  p_vs_NTGraph = new TGraph();
  //  p_vs_NTGraph->SetName("p_vs_NTGraph");

  // define the canvas
  p_vs_NCanvas = new TCanvas("p_vs_NCanvas", "p_vs_NCanvas", 1020, 520);
  p_vs_NCanvas->SetFillColor(0);
  p_vs_NCanvas->SetBorderMode(0);
  p_vs_NCanvas->SetLeftMargin(0.14);
  p_vs_NCanvas->SetRightMargin(0.16);
  p_vs_NCanvas->SetBottomMargin(0.14);

  // draw the histogram
  tree->Draw("p:nFoundHits>>p_vs_N");
  p_vs_N->Draw("colz");

  //  // Fill TGraph
  //
  //  TTreeReader reader("demo/Tree", file);
  //
  //  using doubleVector = std::vector<double>;
  //  TTreeReaderValue<doubleVector> foundHitsRV(reader, "nFoundHits");
  //  TTreeReaderValue<doubleVector> pRV(reader, "p");
  //  while (reader.Next()) {
  //    std::vector<double> foundHitsVector = *foundHitsRV;
  //    std::vector<double> pVector = *pRV;
  //    for (unsigned int i = 0; i < foundHitsVector.size(); i++) {
  //      p_vs_NTGraph->SetPoint(i, foundHitsVector[i], pVector[i]);
  //    }
  //  }
  //
  //  p_vs_NTGraph->Draw("A*");

  // save the histogram
  p_vs_NCanvas->Print("p_vs_N.eps");

  // dEdx_vs_N

  // define the histogram
  dEdx_vs_N = new TH2D("dEdx_vs_N", "dEdx_vs_N", 40, 0, 40, 100, 0, 10);
  dEdx_vs_N->GetYaxis()->SetTitle("dE/dx in MeV/cm");
  dEdx_vs_N->GetXaxis()->SetTitle("Number of Hits per Track");
  dEdx_vs_N->GetYaxis()->SetTitleOffset(0.6);
  dEdx_vs_N->GetXaxis()->SetTitleOffset(1.3);

  // define the canvas
  dEdx_vs_NCanvas =
      new TCanvas("dEdx_vs_NCanvas", "dEdx_vs_NCanvas", 1020, 520);
  dEdx_vs_NCanvas->SetFillColor(0);
  dEdx_vs_NCanvas->SetBorderMode(0);
  dEdx_vs_NCanvas->SetLeftMargin(0.14);
  dEdx_vs_NCanvas->SetRightMargin(0.16);
  dEdx_vs_NCanvas->SetBottomMargin(0.14);

  // draw the histogram
  tree->Draw("dEdx:nFoundHits>>dEdx_vs_N");
  dEdx_vs_N->Draw("colz");

  // save the histogram
  dEdx_vs_NCanvas->Print("dEdx_vs_N.eps");

  // particleID

  // define the histogram
  particleID = new TH1D("particleID", "particleID", 10000, -5000, 5000);
  particleID->GetYaxis()->SetTitleOffset(0.6);
  particleID->GetXaxis()->SetTitleOffset(1.3);

  // define the canvas
  particleIDCanvas =
      new TCanvas("particleIDCanvas", "particleIDCanvas", 1020, 520);
  particleIDCanvas->SetFillColor(0);
  particleIDCanvas->SetBorderMode(0);
  particleIDCanvas->SetLeftMargin(0.14);
  particleIDCanvas->SetRightMargin(0.16);
  particleIDCanvas->SetBottomMargin(0.14);

  // draw the histogram
  tree->Draw("particleID>>particleID");

  // save the histogram
  particleIDCanvas->Print("particleID.eps");

  // deltaR

  // define the histogram
  deltaR = new TH1D("deltaR", "deltaR", 100, 0, 0.03);
  deltaR->GetYaxis()->SetTitleOffset(0.6);
  deltaR->GetXaxis()->SetTitleOffset(1.3);

  // define the canvas
  deltaRCanvas = new TCanvas("deltaRCanvas", "deltaRCanvas", 1020, 520);
  deltaRCanvas->SetFillColor(0);
  deltaRCanvas->SetBorderMode(0);
  deltaRCanvas->SetLeftMargin(0.14);
  deltaRCanvas->SetRightMargin(0.16);
  deltaRCanvas->SetBottomMargin(0.14);

  // draw the histogram
  tree->Draw("deltaR>>deltaR");

  // save the histogram
  deltaRCanvas->Print("deltaR.eps");

  // dZ

  // define the histogram
  dZ = new TH1D("dZ", "dZ", 100, 0, 0.2);
  dZ->GetYaxis()->SetTitleOffset(0.6);
  dZ->GetXaxis()->SetTitleOffset(1.3);
  dZ->SetLineColor(kBlack);

  dZ_unmatched = new TH1D("dZ_unmatched", "dZ_unmatched", 100, 0, 0.2);
  dZ_unmatched->GetYaxis()->SetTitleOffset(0.6);
  dZ_unmatched->GetXaxis()->SetTitleOffset(1.3);
  dZ_unmatched->SetLineColor(kRed);

  // define the canvas
  dZCanvas = new TCanvas("dZCanvas", "dZCanvas", 1020, 520);
  dZCanvas->SetFillColor(0);
  dZCanvas->SetBorderMode(0);
  dZCanvas->SetLeftMargin(0.14);
  dZCanvas->SetRightMargin(0.16);
  dZCanvas->SetBottomMargin(0.14);

  // draw the histogram
  tree->Draw("dZ>>dZ",
             "(std::abs(particleID)==11)||(std::abs(particleID)==211)||("
             "std::abs(particleID)==321)||(std::abs(particleID)==2212)||(std::"
             "abs(particleID)==13)");
  tree->Draw("dZ>>dZ_unmatched",
             "!((std::abs(particleID)==11)||(std::abs(particleID)==211)||("
             "std::abs(particleID)==321)||(std::abs(particleID)==2212)||(std::"
             "abs(particleID)==13))");
  dZStack->Add(dZ);
  dZStack->Add(dZ_unmatched);
  dZStack->Draw("nostack");
  dZStack->GetXaxis()->SetTitle("dZ");
  dZStack->GetYaxis()->SetTitle("Number of tracks");
  dZStack->GetYaxis()->SetTitleOffset(1.3);
  dZStack->GetXaxis()->SetTitleOffset(1.3);

  // save the histogram
  dZCanvas->Print("dZ.eps");
}
