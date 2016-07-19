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
#include "TF1.h"
#include "TPaveStats.h"
#include "TProfile.h"

using doubleVector = std::vector<double>;

double fitFunction(double *var, double *par);
double landauVavilovBichsel(double *var, double *par);
double betheBloch(double *var, double *par);
void dEdx_vs_p(bool save, int pBins, double pMin, double pMax, int dEdxbins,
               double dEdxMin, double dEdxMax, std::string estimator);
void protonGraph(bool save, double pMin, double pMax, double dEdxMin,
                 double dEdxMax, std::string estimator);
void kaonGraph(bool save, double pMin, double pMax, double dEdxMin,
               double dEdxMax, std::string estimator);
void pionGraph(bool save, double pMin, double pMax, double dEdxMin,
               double dEdxMax, std::string estimator);
void bSMGraph(bool save, double pMin, double pMax, double dEdxMin,
              double dEdxMax, std::string estimator);
void massHistogram(bool save, int bins, double minMass, double maxMass,
                   std::string estimator, double dEdxCut);
void pixelVsStripsVsEstimators(bool save, double pTCut);
void dEdx_vs_eta(bool save, int etaBins, double etaMin, double etaMax,
                 int dEdxBins, double dEdxMin, double dEdxMax,
                 std::string estimator, double pTCut);
void dEdx_vs_phi(bool save, int phiBins, double phiMin, double phiMax,
                 int dEdxBins, double dEdxMin, double dEdxMax,
                 std::string estimator);
void dEdx_vs_runID(bool save, int runBins, double firstRun, double lastRun,
                   int dEdxBins, double dEdxMin, double dEdxMax,
                   std::string estimator);
void dEdxForOneROC(bool save, unsigned int Id, int bins, double dEdxMin,
                   double dEdxMax);
void dEdxForOneAPV(bool save, unsigned int Id, int bins, double dEdxMin,
                   double dEdxMax);

TF1 *protonFunction = new TF1("protonFunction", fitFunction, 0.5, 5000, 3);
TF1 *kaonFunction = new TF1("kaonFunction", fitFunction, 0.2, 5000, 3);
TF1 *pionFunction = new TF1("pionFunction", fitFunction, 0.01, 5000, 3);
TF1 *HSCPFunction = new TF1("HSCPFunction", fitFunction, 0, 5000, 3);
TF1 *HSCPBethe = new TF1("HSCPBethe", betheBloch, 0, 5000, 3);
TF1 *HSCPLandau = new TF1("HSCPLandau", landauVavilovBichsel, 0, 5000, 3);
TF1 *LVBProton = new TF1("LVBProton", landauVavilovBichsel, 0, 5000, 3);
TF1 *LVBKaon = new TF1("LVBKaon", landauVavilovBichsel, 0, 5000, 3);
TF1 *LVBPion = new TF1("LVBPion", landauVavilovBichsel, 0, 5000, 3);
TF1 *BetheBlochProton = new TF1("BetheBlochProton", betheBloch, 0, 5000, 3);

// open root file
std::string path_("/nfs/dust/cms/user/kramerto/");
std::string name_("Run2016B_SingleMuon_HSCPCalibration");
std::string nTuple_("_NTuple");
std::string root_(".root");
TFile *file_ = new TFile((path_ + name_ + nTuple_ + root_).c_str(), "READ");

// get the trees
TTree *tree_ = (TTree *)file_->Get("demo/Tree");

void histogram_macro() {

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);

  dEdx_vs_p(false, 5000, 0, 250, 1000, 0, 10, "recoHarmonic2");
  //  dEdx_vs_p(false, 5000, 15, 250, 1000, 0, 10, "myHarmonic2Strip");
  //  protonGraph(false, 0, 5, 0, 15, "recoHarmonic2"); // fit
  //  kaonGraph(false, 0, 5, 0, 15, "recoHarmonic2");
  //  pionGraph(false, 0, 5, 0, 15, "recoHarmonic2");
  //  bSMGraph(false, 0, 5, 0, 15, "recoHarmonic2");
  //    massHistogram(false, 100, 0, 3, "recoHarmonic2", 4.5);
  //  pixelVsStripsVsEstimators(false, 15);
  //  dEdx_vs_eta(false, 50, -2.3, 2.3, 1000, 0, 10, "recoHarmonic2", 15);
  //  dEdx_vs_eta(false, 115, -2.3, 2.3, 1000, 0, 10, "myHarmonic2", 15);
  //  dEdx_vs_eta(false, 50, -2.3, 2.3, 1000, 0, 10, "myHarmonic2Pixel", 15);
  //  dEdx_vs_eta(false, 50, -2.3, 2.3, 1000, 0, 10, "myHarmonic2Strip", 15);
  //  dEdx_vs_phi(false, 50, 0, 3.1416, 1000, 0, 10, "recoHarmonic2");
  //  dEdx_vs_phi(false, 100, 0, 3.1416, 1000, 0, 10, "myHarmonic2");
  //  dEdx_vs_phi(false, 50, 0, 3.1416, 1000, 0, 10, "myHarmonic2Pixel");
  //  dEdx_vs_phi(false, 50, 0, 3.1416, 1000, 0, 10, "myHarmonic2Strip");
  //  dEdx_vs_runID(false, 2210, 272760, 274970, 1000, 0, 10, "recoHarmonic2");
  //  dEdx_vs_runID(false, 2210, 272760, 274970, 1000, 0, 10, "myHarmonic2");
  //  dEdx_vs_runID(false, 2210, 272760, 274970, 1000, 0, 10,
  //  "myHarmonic2Pixel");
  //  dEdx_vs_runID(false, 2210, 272760, 274970, 1000, 0, 10,
  //  "myHarmonic2Strip");
  //  dEdxForOneROC(false, 2417004741, 50, 0, 6);
  //  dEdxForOneAPV(false, 3762756739, 50, 0, 6);
}

double fitFunction(double *var, double *par) {
  return (par[0] + (pow(par[2], 2) * par[1]) / pow(var[0], 2));
}

double landauVavilovBichsel(double *var, double *par) {
  double pi = 4 * atan(1);
  double m_e = 0.000511;  // GeV
  double I = 0.000000173; // GeV
  double Z = 14.0;
  double A = 28.085;
  double K = 0.307075;
  double ro_Si = 2.3296;
  double E = std::sqrt(pow(var[0], 2) + (pow(par[0], 2)));
  double beta = (var[0]) / (E);
  double gamma = 1 / std::sqrt(1 - pow(beta, 2));
  double W_max = 2 * m_e * pow(beta * gamma, 2) /
                 (1 + 2 * gamma * m_e / par[0] + pow(m_e / par[0], 2));
  double plasma_energy_Si = 0.000000031052829448;
  double delta = 2 * log(plasma_energy_Si / I) + 2 * log(beta * gamma) + 1;
  double xi = (K / 2) * (Z / A) * (1 * ro_Si) / pow(beta, 2);

  return (par[1] *
          (xi * (log(2 * m_e * pow(beta * gamma, 2) / I) + log(0.032 * xi / I) +
                 0.2 - pow(beta, 2) - par[2] * log(beta * gamma))));
}

double betheBloch(double *var, double *par) {
  double pi = 4 * atan(1);
  double m_e = 0.000511;  // GeV
  double I = 0.000000173; // GeV
  double Z = 14.0;
  double A = 28.085;
  double K = 0.307075;
  double ro_Si = 2.3296;
  double E = std::sqrt(pow(var[0], 2) + (pow(par[0], 2)));
  double beta = (var[0]) / (E);
  double gamma = 1 / std::sqrt(1 - pow(beta, 2));
  double W_max = 2 * m_e * pow(beta * gamma, 2) /
                 (1 + 2 * gamma * m_e / par[0] + pow(m_e / par[0], 2));
  double plasma_energy_Si = 0.000000031052829448;
  double delta = 2 * log(plasma_energy_Si / I) + 2 * log(beta * gamma) + 1;
  double xi = (K / 2) * (Z / A) * (1 * ro_Si) / pow(beta, 2);

  return (Z / A * K / pow(beta, 2) * ro_Si * par[1] *
          ((log(2 * m_e * pow(beta * gamma, 2) * W_max / (pow(I, 2))) / 2) -
           pow(beta, 2) - par[2] * log(beta * gamma)));
}

void dEdx_vs_p(bool save, int pBins, double pMin, double pMax, int dEdxbins,
               double dEdxMin, double dEdxMax, std::string estimator) {
  // Define histograms
  TH2D *dEdx_vs_p_unmatched =
      new TH2D(("dEdx_vs_p_unmatched_" + estimator).c_str(), "dE/dx vs p",
               pBins, pMin, pMax, dEdxbins, dEdxMin, dEdxMax);
  dEdx_vs_p_unmatched->GetXaxis()->SetTitle("p in GeV/c");
  dEdx_vs_p_unmatched->GetYaxis()->SetTitle(
      ("dE/dx (" + estimator + ") in MeV/cm ").c_str());
  dEdx_vs_p_unmatched->GetYaxis()->SetTitleOffset(0.6);
  dEdx_vs_p_unmatched->GetXaxis()->SetTitleOffset(1.3);
  TH2D *dEdx_vs_p_211 =
      (TH2D *)dEdx_vs_p_unmatched->Clone(("dEdx_vs_p_211" + estimator).c_str());
  TH2D *dEdx_vs_p_11 =
      (TH2D *)dEdx_vs_p_unmatched->Clone(("dEdx_vs_p_11" + estimator).c_str());
  TH2D *dEdx_vs_p_13 =
      (TH2D *)dEdx_vs_p_unmatched->Clone(("dEdx_vs_p_13" + estimator).c_str());
  TH2D *dEdx_vs_p_321 =
      (TH2D *)dEdx_vs_p_unmatched->Clone(("dEdx_vs_p_321" + estimator).c_str());
  TH2D *dEdx_vs_p_2212 = (TH2D *)dEdx_vs_p_unmatched->Clone(
      ("dEdx_vs_p_2212" + estimator).c_str());
  TH2D *dEdx_vs_p_BSM =
      (TH2D *)dEdx_vs_p_unmatched->Clone(("dEdx_vs_p_BSM" + estimator).c_str());
  TH2D *dEdx_vs_p_others =
      (TH2D *)dEdx_vs_p_unmatched->Clone("dEdx_vs_p_others");

  // Define canvas
  TCanvas *dEdx_vs_pCanvas = new TCanvas(
      ("dEdx_vs_pCanvas_" + estimator).c_str(), "dEdx_vs_pCanvas", 1020, 520);
  dEdx_vs_pCanvas->SetFillColor(kWhite);
  dEdx_vs_pCanvas->SetBorderMode(0);
  dEdx_vs_pCanvas->SetLeftMargin(0.14);
  dEdx_vs_pCanvas->SetRightMargin(0.16);
  dEdx_vs_pCanvas->SetBottomMargin(0.14);

  // Define Legend
  TLegend *dEdx_vs_pLegend = new TLegend(0.65, 0.4, 0.84, 0.9);
  dEdx_vs_pLegend->SetTextSize(0.05);

  // Fill histograms
  TTreeReader reader("demo/Tree", file_);
  TTreeReaderValue<doubleVector> dEdxRV(reader, estimator.c_str());
  TTreeReaderValue<doubleVector> pRV(reader, "p");
  TTreeReaderValue<std::vector<int>> pIDRV(reader, "matchedID");
  TTreeReaderValue<doubleVector> pTRV(reader, "pT");

  while (reader.Next()) {
    std::vector<double> dEdxVector = *dEdxRV;
    std::vector<double> pVector = *pRV;
    std::vector<int> pIDVector = *pIDRV;
    std::vector<double> pTVector = *pTRV;
    for (unsigned int i = 0; i < dEdxVector.size(); i++) {
      //      if (std::abs(pIDVector[i]) == 211) {
      //        dEdx_vs_p_211->SetMarkerColor(kYellow + 1);
      //        dEdx_vs_p_211->SetFillColor(kYellow + 1);
      //        dEdx_vs_p_211->SetMarkerStyle(8);
      //        dEdx_vs_p_211->SetMarkerSize(0.4);
      //        dEdx_vs_p_211->Fill(pVector[i], dEdxVector[i]);
      //      } else if (std::abs(pIDVector[i]) == 321) {
      //        dEdx_vs_p_321->SetMarkerColor(kGreen + 2);
      //        dEdx_vs_p_321->SetFillColor(kGreen + 2);
      //        dEdx_vs_p_321->SetMarkerStyle(8);
      //        dEdx_vs_p_321->SetMarkerSize(0.4);
      //        dEdx_vs_p_321->Fill(pVector[i], dEdxVector[i]);
      //      } else if (std::abs(pIDVector[i]) == 11) {
      //        dEdx_vs_p_11->SetMarkerColor(kRed);
      //        dEdx_vs_p_11->SetFillColor(kRed);
      //        dEdx_vs_p_11->SetMarkerStyle(8);
      //        dEdx_vs_p_11->SetMarkerSize(0.4);
      //        dEdx_vs_p_11->Fill(pVector[i], dEdxVector[i]);
      //      } else if (std::abs(pIDVector[i]) == 2212) {
      //        dEdx_vs_p_2212->SetMarkerColor(kBlue);
      //        dEdx_vs_p_2212->SetFillColor(kBlue);
      //        dEdx_vs_p_2212->SetMarkerStyle(8);
      //        dEdx_vs_p_2212->SetMarkerSize(0.4);
      //        dEdx_vs_p_2212->Fill(pVector[i], dEdxVector[i]);
      //      } else if (std::abs(pIDVector[i]) == 13) {
      //        dEdx_vs_p_13->SetMarkerColor(kCyan);
      //        dEdx_vs_p_13->SetFillColor(kCyan);
      //        dEdx_vs_p_13->SetMarkerStyle(8);
      //        dEdx_vs_p_13->SetMarkerSize(0.4);
      //        dEdx_vs_p_13->Fill(pVector[i], dEdxVector[i]);
      //      } else if (std::abs(pIDVector[i]) > 1000000) {
      //        dEdx_vs_p_BSM->SetMarkerColor(kMagenta);
      //        dEdx_vs_p_BSM->SetFillColor(kMagenta);
      //        dEdx_vs_p_BSM->SetMarkerStyle(8);
      //        dEdx_vs_p_BSM->SetMarkerSize(0.4);
      //        dEdx_vs_p_BSM->Fill(pVector[i], dEdxVector[i]);
      //      } else if (pIDVector[i] == 0) {
      //        dEdx_vs_p_unmatched->SetMarkerColor(kGray);
      //        dEdx_vs_p_unmatched->SetFillColor(kGray);
      //        dEdx_vs_p_unmatched->SetMarkerStyle(8);
      //        dEdx_vs_p_unmatched->SetMarkerSize(0.4);
      dEdx_vs_p_unmatched->Fill(pVector[i], dEdxVector[i]);
      //      } else {
      //        dEdx_vs_p_others->SetMarkerColor(kBlack);
      //        dEdx_vs_p_others->SetFillColor(kBlack);
      //        dEdx_vs_p_others->SetMarkerStyle(8);
      //        dEdx_vs_p_others->SetMarkerSize(0.4);
      //        dEdx_vs_p_others->Fill(pVector[i], dEdxVector[i]);
      //      }
    }
  }

  // Draw histogram
  dEdx_vs_p_unmatched->Draw("colz");
  //  dEdx_vs_p_others->Draw("same");
  //  dEdx_vs_p_211->Draw("same");
  //  dEdx_vs_p_321->Draw("same");
  //  dEdx_vs_p_2212->Draw("same");
  //  dEdx_vs_p_11->Draw("same");
  //  dEdx_vs_p_13->Draw("same");
  //  dEdx_vs_p_BSM->Draw("same");

  // Fill and draw legend
  //dEdx_vs_pLegend->AddEntry(dEdx_vs_p_unmatched, "unmatched", "f");
  //dEdx_vs_pLegend->AddEntry(dEdx_vs_p_211, "#pi^{+/-}", "f");
  //dEdx_vs_pLegend->AddEntry(dEdx_vs_p_321, "K^{+/-}", "f");
  //dEdx_vs_pLegend->AddEntry(dEdx_vs_p_2212, "p^{+/-}", "f");
  //dEdx_vs_pLegend->AddEntry(dEdx_vs_p_11, "e^{+/-}", "f");
  //dEdx_vs_pLegend->AddEntry(dEdx_vs_p_13, "#mu^{+/-}", "f");
  //dEdx_vs_pLegend->AddEntry(dEdx_vs_p_BSM, "BSM particles", "f");
  //dEdx_vs_pLegend->AddEntry(dEdx_vs_p_others, "others", "f");
  //dEdx_vs_pLegend->AddEntry(HSCPFunction, "protonFit", "l");
  //dEdx_vs_pLegend->AddEntry(HSCPBethe, "Bethe Bloch", "l");
  //dEdx_vs_pLegend->AddEntry(HSCPLandau, "Landau/V./B.", "l");
  //dEdx_vs_pLegend->Draw();

  // Save
  if (save)
    dEdx_vs_pCanvas->Print((name_ + estimator + "_dEdx_vs_p.eps").c_str());
}

void protonGraph(bool save, double pMin, double pMax, double dEdxMin,
                 double dEdxMax, std::string estimator) {
  TGraph *dEdx_vs_p2212Graph;
  TCanvas *dEdx_vs_p2212GraphCanvas;
  std::vector<double> p2212Fit;
  std::vector<double> dEdx2212Fit;
  std::vector<double> p2212;
  std::vector<double> dEdx2212;

  // define PROTON TGRAPH(fit) Canvas
  dEdx_vs_p2212GraphCanvas = new TCanvas("dEdx_vs_p2212GraphCanvas",
                                         "dEdx_vs_p2212GraphCanvas", 1020, 520);
  dEdx_vs_p2212GraphCanvas->SetFillColor(kWhite);
  dEdx_vs_p2212GraphCanvas->SetBorderMode(0);
  dEdx_vs_p2212GraphCanvas->SetLeftMargin(0.14);
  dEdx_vs_p2212GraphCanvas->SetRightMargin(0.16);
  dEdx_vs_p2212GraphCanvas->SetBottomMargin(0.14);

  // Fill proton TGraph
  TTreeReader reader("demo/Tree", file_);
  TTreeReaderValue<doubleVector> dEdxRV(reader, estimator.c_str());
  TTreeReaderValue<doubleVector> pRV(reader, "p");
  TTreeReaderValue<std::vector<int>> pIDRV(reader, "matchedID");

  while (reader.Next()) {
    std::vector<double> dEdxVector = *dEdxRV;
    std::vector<double> pVector = *pRV;
    std::vector<int> pIDVector = *pIDRV;
    for (unsigned int i = 0; i < dEdxVector.size(); i++) {
      if (std::abs(pIDVector[i]) == 2212) {
        dEdx2212.push_back(dEdxVector[i]);
        p2212.push_back(pVector[i]);
        if (((dEdxVector[i] > 7) && (dEdxVector[i] < 11)) ||
            ((pVector[i] > 2) && (pVector[i] < 4))) {
          dEdx2212Fit.push_back(dEdxVector[i]);
          p2212Fit.push_back(pVector[i]);
        }
      }
    }
  }

  // define proton TGraph
  dEdx_vs_p2212Graph =
      new TGraph(p2212Fit.size(), &(p2212Fit[0]), &(dEdx2212Fit[0]));

  dEdx_vs_p2212Graph->Draw("");
  protonFunction->SetParNames("C", "K", "M");
  protonFunction->FixParameter(2, 0.9382720813);

  // fit
  dEdx_vs_p2212Graph->Fit("protonFunction", "MR");

  // draw proton TGraph (and fits)
  dEdx_vs_p2212Graph = new TGraph(p2212.size(), &(p2212[0]), &(dEdx2212[0]));
  dEdx_vs_p2212Graph->SetNameTitle("dEdx_vs_p2212Graph",
                                   "dE/dx vs p (protons)");
  dEdx_vs_p2212Graph->GetXaxis()->SetRangeUser(pMin, pMax);
  dEdx_vs_p2212Graph->GetYaxis()->SetRangeUser(dEdxMin, dEdxMax);
  dEdx_vs_p2212Graph->GetXaxis()->SetTitle("p in GeV/c");
  dEdx_vs_p2212Graph->GetYaxis()->SetTitle(
      ("dE/dx Estimator (" + estimator + ") in MeV/cm ").c_str());
  dEdx_vs_p2212Graph->Draw("AP");
  protonFunction->Draw("same");
  dEdx_vs_p2212GraphCanvas->Update();
  // Save
  if (save)
    dEdx_vs_p2212GraphCanvas->Print(
        (name_ + estimator + "_dEdx_vs_p_protonTGraph.eps").c_str());
}
void kaonGraph(bool save, double pMin, double pMax, double dEdxMin,
               double dEdxMax, std::string estimator) {
  TGraph *dEdx_vs_p321Graph;
  TCanvas *dEdx_vs_p321GraphCanvas;
  std::vector<double> p321;
  std::vector<double> dEdx321;
  std::vector<double> p321Fit;
  std::vector<double> dEdx321Fit;

  // define KAON TGRAPH(fit) Canvas
  dEdx_vs_p321GraphCanvas = new TCanvas("dEdx_vs_p321GraphCanvas",
                                        "dEdx_vs_p321GraphCanvas", 1020, 520);
  dEdx_vs_p321GraphCanvas->SetFillColor(kWhite);
  dEdx_vs_p321GraphCanvas->SetBorderMode(0);
  dEdx_vs_p321GraphCanvas->SetLeftMargin(0.14);
  dEdx_vs_p321GraphCanvas->SetRightMargin(0.16);
  dEdx_vs_p321GraphCanvas->SetBottomMargin(0.14);

  // Fill kaon TGraph
  TTreeReader reader("demo/Tree", file_);
  TTreeReaderValue<doubleVector> dEdxRV(reader, estimator.c_str());
  TTreeReaderValue<doubleVector> pRV(reader, "p");
  TTreeReaderValue<std::vector<int>> pIDRV(reader, "matchedID");

  while (reader.Next()) {
    std::vector<double> dEdxVector = *dEdxRV;
    std::vector<double> pVector = *pRV;
    std::vector<int> pIDVector = *pIDRV;
    for (unsigned int i = 0; i < dEdxVector.size(); i++) {
      if (std::abs(pIDVector[i]) == 321) {
        dEdx321.push_back(dEdxVector[i]);
        p321.push_back(pVector[i]);
        //        if ((dEdxVector[i] > 4) || (pVector[i] > 0.8)) {
        //          dEdx321Fit.push_back(dEdxVector[i]);
        //          p321Fit.push_back(pVector[i]);
        //      }
      }
    }
  }

  // define kaon TGraph
  dEdx_vs_p321Graph =
      new TGraph(p321Fit.size(), &(p321Fit[0]), &(dEdx321Fit[0]));

  dEdx_vs_p321Graph->Draw("");
  kaonFunction->SetParameter(0, protonFunction->GetParameter(0));
  kaonFunction->SetParameter(1, protonFunction->GetParameter(1));
  kaonFunction->FixParameter(2, 0.493667);
  kaonFunction->Draw("same");

  dEdx_vs_p321GraphCanvas->Update();

  // fit
  //  dEdx_vs_p321Graph->Fit("kaonFunction", "MR");

  // draw proton TGraph (and fits)
  dEdx_vs_p321Graph = new TGraph(p321.size(), &(p321[0]), &(dEdx321[0]));
  dEdx_vs_p321Graph->SetNameTitle("dEdx_vs_p321Graph", "dE/dx vs p (kaons)");
  dEdx_vs_p321Graph->GetXaxis()->SetRangeUser(pMin, pMax);
  dEdx_vs_p321Graph->GetYaxis()->SetRangeUser(dEdxMin, dEdxMax);
  dEdx_vs_p321Graph->GetXaxis()->SetTitle("p in GeV/c");
  dEdx_vs_p321Graph->GetYaxis()->SetTitle(
      ("dE/dx Estimator (" + estimator + ") in MeV/cm ").c_str());
  dEdx_vs_p321Graph->Draw("AP");
  kaonFunction->Draw("same");
  dEdx_vs_p321GraphCanvas->Update();
  // Save
  if (save)
    dEdx_vs_p321GraphCanvas->Print(
        (name_ + estimator + "_dEdx_vs_p_kaonTGraph.eps").c_str());
}

void pionGraph(bool save, double pMin, double pMax, double dEdxMin,
               double dEdxMax, std::string estimator) {
  TGraph *dEdx_vs_p211Graph;
  TCanvas *dEdx_vs_p211GraphCanvas;
  std::vector<double> p211;
  std::vector<double> dEdx211;

  // define PION TGRAPH Canvas
  dEdx_vs_p211GraphCanvas = new TCanvas("dEdx_vs_p211GraphCanvas",
                                        "dEdx_vs_p211GraphCanvas", 1020, 520);
  dEdx_vs_p211GraphCanvas->SetFillColor(kWhite);
  dEdx_vs_p211GraphCanvas->SetBorderMode(0);
  dEdx_vs_p211GraphCanvas->SetLeftMargin(0.14);
  dEdx_vs_p211GraphCanvas->SetRightMargin(0.16);
  dEdx_vs_p211GraphCanvas->SetBottomMargin(0.14);

  // fill pion TGraph
  TTreeReader reader("demo/Tree", file_);
  TTreeReaderValue<doubleVector> dEdxRV(reader, estimator.c_str());
  TTreeReaderValue<doubleVector> pRV(reader, "p");
  TTreeReaderValue<std::vector<int>> pIDRV(reader, "matchedID");
  while (reader.Next()) {
    std::vector<double> dEdxVector = *dEdxRV;
    std::vector<double> pVector = *pRV;
    std::vector<int> pIDVector = *pIDRV;
    for (unsigned int i = 0; i < dEdxVector.size(); i++) {
      if (std::abs(pIDVector[i]) == 211) {
        dEdx211.push_back(dEdxVector[i]);
        p211.push_back(pVector[i]);
      }
    }
  }

  // define pion TGraph
  dEdx_vs_p211Graph = new TGraph(p211.size(), &(p211[0]), &(dEdx211[0]));
  dEdx_vs_p211Graph->SetNameTitle("dEdx_vs_p321Graph", "dE/dx vs p (pions)");
  dEdx_vs_p211Graph->GetXaxis()->SetRangeUser(pMin, pMax);
  dEdx_vs_p211Graph->GetYaxis()->SetRangeUser(dEdxMin, dEdxMax);
  dEdx_vs_p211Graph->GetXaxis()->SetTitle("p in GeV/c");
  dEdx_vs_p211Graph->GetYaxis()->SetTitle(
      ("dE/dx Estimator (" + estimator + ") in MeV/cm ").c_str());
  dEdx_vs_p211Graph->Draw("AP");
  pionFunction->SetParameter(0, protonFunction->GetParameter(0));
  pionFunction->SetParameter(1, protonFunction->GetParameter(1));
  pionFunction->FixParameter(2, 0.13957018);
  pionFunction->Draw("same");
  //  LVBPion->SetParNames("PionMass", "ScaleFactor",
  //  "densityCorrectionFactor");
  //  LVBPion->SetParameter(0, 1);
  //  LVBPion->FixParameter(1, LVBProton->GetParameter(1));
  //  LVBPion->FixParameter(2, LVBProton->GetParameter(2));
  //  LVBPion->SetNpx(1000);
  //  dEdx_vs_p211Graph->Fit("LVBPion", "MR");
  //  LVBPion->Draw("same");

  dEdx_vs_p211GraphCanvas->Update();
  if (save)
    dEdx_vs_p211GraphCanvas->Print(
        (name_ + estimator + "_dEdx_vs_p_pionTGraph.eps").c_str());
}

void bSMGraph(bool save, double pMin, double pMax, double dEdxMin,
              double dEdxMax, std::string estimator) {
  TGraph *dEdx_vs_pBSMGraph;
  TCanvas *dEdx_vs_pBSMGraphCanvas;
  std::vector<double> pBSM;
  std::vector<double> dEdxBSM;

  // define PION TGRAPH Canvas
  dEdx_vs_pBSMGraphCanvas = new TCanvas("dEdx_vs_pBSMGraphCanvas",
                                        "dEdx_vs_pBSMGraphCanvas", 1020, 520);
  dEdx_vs_pBSMGraphCanvas->SetFillColor(kWhite);
  dEdx_vs_pBSMGraphCanvas->SetBorderMode(0);
  dEdx_vs_pBSMGraphCanvas->SetLeftMargin(0.14);
  dEdx_vs_pBSMGraphCanvas->SetRightMargin(0.16);
  dEdx_vs_pBSMGraphCanvas->SetBottomMargin(0.14);

  // fill pion TGraph
  TTreeReader reader("demo/Tree", file_);
  TTreeReaderValue<doubleVector> dEdxRV(reader, estimator.c_str());
  TTreeReaderValue<doubleVector> pRV(reader, "p");
  TTreeReaderValue<std::vector<int>> pIDRV(reader, "matchedID");
  while (reader.Next()) {
    std::vector<double> dEdxVector = *dEdxRV;
    std::vector<double> pVector = *pRV;
    std::vector<int> pIDVector = *pIDRV;
    for (unsigned int i = 0; i < dEdxVector.size(); i++) {
      if (std::abs(pIDVector[i]) > 1000000) {
        dEdxBSM.push_back(dEdxVector[i]);
        pBSM.push_back(pVector[i]);
      }
    }
  }

  // define pion TGraph
  dEdx_vs_pBSMGraph = new TGraph(pBSM.size(), &(pBSM[0]), &(dEdxBSM[0]));
  dEdx_vs_pBSMGraph->SetNameTitle("dEdx_vs_p321Graph", "dE/dx vs p (pions)");
  dEdx_vs_pBSMGraph->GetXaxis()->SetRangeUser(pMin, pMax);
  dEdx_vs_pBSMGraph->GetYaxis()->SetRangeUser(dEdxMin, dEdxMax);
  dEdx_vs_pBSMGraph->GetXaxis()->SetTitle("p in GeV/c");
  dEdx_vs_pBSMGraph->GetYaxis()->SetTitle(
      ("dE/dx Estimator (" + estimator + ") in MeV/cm ").c_str());
  dEdx_vs_pBSMGraph->Draw("AP");
  pionFunction->SetParameter(0, protonFunction->GetParameter(0));
  pionFunction->SetParameter(1, protonFunction->GetParameter(1));
  pionFunction->FixParameter(2, 0.13957018);
  pionFunction->Draw("same");

  dEdx_vs_pBSMGraphCanvas->Update();
  if (save)
    dEdx_vs_pBSMGraphCanvas->Print(
        (name_ + estimator + "_dEdx_vs_p_BSMTGraph.eps").c_str());
}

void massHistogram(bool save, int bins, double minMass, double maxMass,
                   std::string estimator, double dEdxCut) {
  TH1D *massHisto;
  TH1D *massHisto2212;
  TH1D *massHisto321;
  TH1D *massHistoBSM;
  TCanvas *massHistoCanvas;
  TLegend *massHistoLegend;

  std::vector<double> p2212;
  std::vector<double> dEdx2212;
  std::vector<double> p321;
  std::vector<double> dEdx321;
  std::vector<double> pAllParticles;
  std::vector<double> dEdxAllParticles;
  std::vector<double> pBSM;
  std::vector<double> dEdxBSM;

  TTreeReader reader("demo/Tree", file_);
  TTreeReaderValue<doubleVector> dEdxRV(reader, estimator.c_str());
  TTreeReaderValue<doubleVector> pRV(reader, "p");
  TTreeReaderValue<std::vector<int>> pIDRV(reader, "matchedID");
  TTreeReaderValue<doubleVector> pTRV(reader, "pT");

  while (reader.Next()) {
    std::vector<double> dEdxVector = *dEdxRV;
    std::vector<double> pVector = *pRV;
    std::vector<int> pIDVector = *pIDRV;
    for (unsigned int i = 0; i < dEdxVector.size(); i++) {
      if (std::abs(pIDVector[i]) == 321) {
        dEdx321.push_back(dEdxVector[i]);
        p321.push_back(pVector[i]);
        dEdxAllParticles.push_back(dEdxVector[i]);
        pAllParticles.push_back(pVector[i]);
      } else if (std::abs(pIDVector[i]) == 2212) {
        dEdx2212.push_back(dEdxVector[i]);
        p2212.push_back(pVector[i]);
        dEdxAllParticles.push_back(dEdxVector[i]);
        pAllParticles.push_back(pVector[i]);
      } else if (std::abs(pIDVector[i]) > 1000000) {
        dEdxBSM.push_back(dEdxVector[i]);
        pBSM.push_back(pVector[i]);
        dEdxAllParticles.push_back(dEdxVector[i]);
        pAllParticles.push_back(pVector[i]);
      } else {
        dEdxAllParticles.push_back(dEdxVector[i]);
        pAllParticles.push_back(pVector[i]);
      }
    }
  }

  std::ostringstream strs;
  strs << dEdxCut;
  std::string dEdxCutString = strs.str();
  std::string title = "Reconstructed mass (dE/dx > ";
  title.append(dEdxCutString);
  title.append(" MeV/cm)");
  massHisto = new TH1D("massHisto", title.c_str(), bins, minMass, maxMass);
  massHisto->SetLineColor(kBlack);
  massHisto2212 =
      new TH1D("massHisto2212", title.c_str(), bins, minMass, maxMass);
  massHisto2212->SetLineColor(kBlue);
  massHisto321 =
      new TH1D("massHisto321", title.c_str(), bins, minMass, maxMass);
  massHisto321->SetLineColor(kGreen + 2);
  massHistoBSM =
      new TH1D("massHistoBSM", title.c_str(), bins, minMass, maxMass);
  massHistoBSM->SetLineColor(kMagenta);

  for (unsigned int i = 0; i < dEdxAllParticles.size(); i++) {
    if (dEdxAllParticles[i] > dEdxCut) {
      massHisto->Fill(
          pAllParticles[i] *
          std::sqrt((dEdxAllParticles[i] - protonFunction->GetParameter(0)) /
                    protonFunction->GetParameter(1)));
    }
  }

  for (unsigned int i = 0; i < dEdx2212.size(); i++) {
    if (dEdx2212[i] > dEdxCut) {
      massHisto2212->Fill(
          p2212[i] * std::sqrt((dEdx2212[i] - protonFunction->GetParameter(0)) /
                               protonFunction->GetParameter(1)));
    }
  }
  for (unsigned int i = 0; i < dEdx321.size(); i++) {
    if (dEdx321[i] > dEdxCut) {
      massHisto321->Fill(
          p321[i] * std::sqrt((dEdx321[i] - protonFunction->GetParameter(0)) /
                              protonFunction->GetParameter(1)));
    }
  }

  for (unsigned int i = 0; i < dEdxBSM.size(); i++) {
    if (dEdxBSM[i] > dEdxCut) {
      massHistoBSM->Fill(
          pBSM[i] * std::sqrt((dEdxBSM[i] - protonFunction->GetParameter(0)) /
                              protonFunction->GetParameter(1)));
    }
  }

  massHistoCanvas =
      new TCanvas("massHistoCanvas", "massHistoCanvas", 1020, 520);
  massHistoCanvas->SetFillColor(kWhite);
  massHistoCanvas->SetBorderMode(0);
  massHistoCanvas->SetLeftMargin(0.14);
  massHistoCanvas->SetRightMargin(0.16);
  massHistoCanvas->SetBottomMargin(0.14);

  massHisto->GetXaxis()->SetTitle("m in GeV/c^{2}");
  massHisto->GetYaxis()->SetTitle("N_{tracks}/bin");

  massHisto->Draw();
  massHisto2212->Draw("same");
  massHisto321->Draw("same");
  massHistoBSM->Draw("same");

  // legend
  massHistoLegend = new TLegend(0.64, 0.45, 0.84, 0.75);
  massHistoLegend->SetTextSize(0.05);
  massHistoLegend->AddEntry(massHisto321, "K^{+/-}", "l");
  massHistoLegend->AddEntry(massHisto2212, "p^{+/-}", "l");
  massHistoLegend->AddEntry(massHistoBSM, "BSM particles", "l");
  massHistoLegend->Draw();

  if (save)
    massHistoCanvas->Print((name_ + estimator + "_massHistogram.eps").c_str());
}
void pixelVsStripsVsEstimators(bool save, double pTCut) {
  std::ostringstream strs;
  strs << pTCut;
  std::string pTCutString = strs.str();
  std::string title = "dE/dx (Pixel vs Strip vs Estimators), RECO, p_{T} > " +
                      pTCutString + " GeV/c";
  TH1D *dEdxPixel = new TH1D("dEdxPixel", title.c_str(), 70, 0, 7);
  dEdxPixel->GetXaxis()->SetTitle("dE/dx (Harmonic-2) in MeV/cm");
  dEdxPixel->GetYaxis()->SetTitle("N_{tracks}/bin");
  TH1D *dEdxStrip = (TH1D *)dEdxPixel->Clone("dEdxStrip");
  TH1D *dEdxMyHarmonic2 = (TH1D *)dEdxPixel->Clone("dEdxMyHarmonic2");
  TH1D *dEdxRecoHarmonic2 = (TH1D *)dEdxPixel->Clone("dEdxRecoHarmonic2");

  TTreeReader reader("demo/Tree", file_);
  TTreeReaderValue<doubleVector> dEdxRV(reader, "recoHarmonic2");
  TTreeReaderValue<doubleVector> myHarmonic2RV(reader, "myHarmonic2");
  TTreeReaderValue<doubleVector> nFoundHitsRV(reader, "nFoundHits");
  TTreeReaderValue<std::vector<int>> pIDRV(reader, "matchedID");
  TTreeReaderValue<doubleVector> pTRV(reader, "pT");
  TTreeReaderValue<doubleVector> myHarmonic2PixelRV(reader, "myHarmonic2Pixel");
  TTreeReaderValue<doubleVector> myHarmonic2StripRV(reader, "myHarmonic2Strip");

  while (reader.Next()) {
    std::vector<double> dEdxVector = *dEdxRV;
    std::vector<double> nFoundHitsVector = *nFoundHitsRV;
    std::vector<int> pIDVector = *pIDRV;
    std::vector<double> pTVector = *pTRV;
    std::vector<double> myHarmonic2Vector = *myHarmonic2RV;
    std::vector<double> myHarmonic2PixelVector = *myHarmonic2PixelRV;
    std::vector<double> myHarmonic2StripVector = *myHarmonic2StripRV;
    for (unsigned int i = 0; i < dEdxVector.size(); i++) {
      if (pTVector[i] > pTCut && nFoundHitsVector[i] > 9) {
        if (myHarmonic2PixelVector[i] > 0)
          dEdxPixel->Fill(myHarmonic2PixelVector[i]);
        if (myHarmonic2StripVector[i] > 0)
          dEdxStrip->Fill(myHarmonic2StripVector[i]);
        if (myHarmonic2Vector[i] > 0)
          dEdxMyHarmonic2->Fill(myHarmonic2Vector[i]);
        if (dEdxVector[i] > 0)
          dEdxRecoHarmonic2->Fill(dEdxVector[i]);
      }
    }
  }

  TCanvas *dEdxPixelVsStripCanvas = new TCanvas(
      "dEdxPixelVsStripCanvas", "dEdxPixelVsStripCanvas", 1020, 520);
  dEdxPixelVsStripCanvas->SetFillColor(0);
  dEdxPixelVsStripCanvas->SetBorderMode(0);
  dEdxPixelVsStripCanvas->SetLeftMargin(0.14);
  dEdxPixelVsStripCanvas->SetRightMargin(0.16);
  dEdxPixelVsStripCanvas->SetBottomMargin(0.14);

  dEdxMyHarmonic2->SetLineColor(kAzure + 8);
  dEdxMyHarmonic2->Draw("same");
  dEdxRecoHarmonic2->SetLineColor(kBlack);
  dEdxRecoHarmonic2->Draw("same");
  dEdxStrip->SetLineColor(kRed);
  dEdxStrip->Draw("same");
  dEdxPixel->SetLineColor(kGreen + 2);
  dEdxPixel->Draw("same");

  // legend
  TLegend *dEdxPixelVsStripLegend = new TLegend(0.64, 0.45, 0.84, 0.75);
  dEdxPixelVsStripLegend->SetTextSize(0.05);
  dEdxPixelVsStripLegend->AddEntry(dEdxRecoHarmonic2, "Reco", "l");
  dEdxPixelVsStripLegend->AddEntry(dEdxMyHarmonic2, "Pixel+Strip", "l");
  dEdxPixelVsStripLegend->AddEntry(dEdxStrip, "Strip", "l");
  dEdxPixelVsStripLegend->AddEntry(dEdxPixel, "Pixel", "l");
  dEdxPixelVsStripLegend->Draw();

  if (save)
    dEdxPixelVsStripCanvas->Print(
        (name_ + "_stripsVsPixelVsEstimators.eps").c_str());
}
void dEdx_vs_eta(bool save, int etaBins, double etaMin, double etaMax,
                 int dEdxBins, double dEdxMin, double dEdxMax,
                 std::string estimator, double pTCut) {
  std::ostringstream strs;
  strs << pTCut;
  std::string pTCutString = strs.str();
  std::string title =
      "dE/dx(" + estimator + ") vs eta, p_{T} > " + pTCutString + " GeV/c";
  TH2D *dEdx_vs_eta_Histo =
      new TH2D(("dEdx_vs_eta_Histo" + estimator).c_str(), title.c_str(),
               etaBins, etaMin, etaMax, dEdxBins, dEdxMin, dEdxMax);
  dEdx_vs_eta_Histo->GetXaxis()->SetTitle("#eta_{track}");
  dEdx_vs_eta_Histo->GetYaxis()->SetTitle(
      ("dE/dx (" + estimator + ") in MeV/cm").c_str());
  TCanvas *dEdx_vs_eta_Canvas =
      new TCanvas(("dEdx_vs_eta_Canvas" + estimator).c_str(),
                  "dEdx_vs_eta_Canvas", 1020, 520);
  dEdx_vs_eta_Canvas->SetFillColor(0);
  dEdx_vs_eta_Canvas->SetBorderMode(0);
  dEdx_vs_eta_Canvas->SetLeftMargin(0.14);
  dEdx_vs_eta_Canvas->SetRightMargin(0.16);
  dEdx_vs_eta_Canvas->SetBottomMargin(0.14);
  tree_->Draw((estimator + ":eta>>dEdx_vs_eta_Histo" + estimator).c_str(),
              ("pT>" + pTCutString + "&&nFoundHits>9").c_str());
  //  tree_->Draw((estimator + ":eta>>dEdx_vs_eta_Histo" + estimator).c_str(),
  //              "p>30&&p<40");
  dEdx_vs_eta_Histo->Draw();
  TCanvas *dEdx_vs_eta_Profile_Canvas =
      new TCanvas(("dEdx_vs_eta_Profile_Canvas" + estimator).c_str(),
                  "dEdx_vs_eta_Profile_Canvas", 1020, 520);
  dEdx_vs_eta_Profile_Canvas->SetFillColor(0);
  dEdx_vs_eta_Profile_Canvas->SetBorderMode(0);
  dEdx_vs_eta_Profile_Canvas->SetLeftMargin(0.14);
  dEdx_vs_eta_Profile_Canvas->SetRightMargin(0.16);
  dEdx_vs_eta_Profile_Canvas->SetBottomMargin(0.14);
  TProfile *dEdx_vs_eta_Profile = dEdx_vs_eta_Histo->ProfileX(
      ("dEdx_vs_eta_Profile" + estimator).c_str(), 0, dEdxBins);
  dEdx_vs_eta_Profile->Draw();
  if (save) {
    dEdx_vs_eta_Canvas->Print((name_ + estimator + "_dEdx_vs_eta.eps").c_str());
    dEdx_vs_eta_Profile_Canvas->Print(
        (name_ + estimator + "_dEdx_vs_eta_Profile.eps").c_str());
  }
}
void dEdx_vs_phi(bool save, int phiBins, double phiMin, double phiMax,
                 int dEdxBins, double dEdxMin, double dEdxMax,
                 std::string estimator) {
  TH2D *dEdx_vs_phi_Histo =
      new TH2D(("dEdx_vs_phi_Histo" + estimator).c_str(),
               ("dE/dx(" + estimator + ") vs phi").c_str(), phiBins, phiMin,
               phiMax, dEdxBins, dEdxMin, dEdxMax);
  dEdx_vs_phi_Histo->GetXaxis()->SetTitle("#varphi_{track}");
  dEdx_vs_phi_Histo->GetYaxis()->SetTitle(
      ("dE/dx (" + estimator + ") in MeV/cm").c_str());
  TCanvas *dEdx_vs_phi_Canvas =
      new TCanvas(("dEdx_vs_phi_Canvas" + estimator).c_str(),
                  "dEdx_vs_phi_Canvas", 1020, 520);
  dEdx_vs_phi_Canvas->SetFillColor(0);
  dEdx_vs_phi_Canvas->SetBorderMode(0);
  dEdx_vs_phi_Canvas->SetLeftMargin(0.14);
  dEdx_vs_phi_Canvas->SetRightMargin(0.16);
  dEdx_vs_phi_Canvas->SetBottomMargin(0.14);
  tree_->Draw((estimator + ":phi>>dEdx_vs_phi_Histo" + estimator).c_str(),
              "pT>15&&nFoundHits>9");
  dEdx_vs_phi_Histo->Draw();
  TCanvas *dEdx_vs_phi_Profile_Canvas =
      new TCanvas(("dEdx_vs_phi_Profile_Canvas" + estimator).c_str(),
                  "dEdx_vs_phi_Profile_Canvas", 1020, 520);
  dEdx_vs_phi_Profile_Canvas->SetFillColor(0);
  dEdx_vs_phi_Profile_Canvas->SetBorderMode(0);
  dEdx_vs_phi_Profile_Canvas->SetLeftMargin(0.14);
  dEdx_vs_phi_Profile_Canvas->SetRightMargin(0.16);
  dEdx_vs_phi_Profile_Canvas->SetBottomMargin(0.14);
  TProfile *dEdx_vs_phi_Profile = dEdx_vs_phi_Histo->ProfileX(
      ("dEdx_vs_phi_Profile" + estimator).c_str(), 0, dEdxBins);
  dEdx_vs_phi_Profile->Draw();
  if (save) {
    dEdx_vs_phi_Canvas->Print((name_ + estimator + "_dEdx_vs_phi.eps").c_str());
    dEdx_vs_phi_Profile_Canvas->Print(
        (name_ + estimator + "_dEdx_vs_phi_Profile.eps").c_str());
  }
}
void dEdx_vs_runID(bool save, int runBins, double firstRun, double lastRun,
                   int dEdxBins, double dEdxMin, double dEdxMax,
                   std::string estimator) {
  TH2D *dEdx_vs_runID_Histo =
      new TH2D(("dEdx_vs_runID_Histo" + estimator).c_str(),
               ("dE/dx(" + estimator + ") vs Run").c_str(), runBins, firstRun,
               lastRun, dEdxBins, dEdxMin, dEdxMax);
  dEdx_vs_runID_Histo->GetXaxis()->SetTitle("Run");
  dEdx_vs_runID_Histo->GetYaxis()->SetTitle(
      ("dE/dx (" + estimator + ") in MeV/cm").c_str());
  TCanvas *dEdx_vs_runID_Canvas =
      new TCanvas(("dEdx_vs_runID_Canvas" + estimator).c_str(),
                  "dEdx_vs_runID_Canvas", 1020, 520);
  dEdx_vs_runID_Canvas->SetFillColor(0);
  dEdx_vs_runID_Canvas->SetBorderMode(0);
  dEdx_vs_runID_Canvas->SetLeftMargin(0.14);
  dEdx_vs_runID_Canvas->SetRightMargin(0.16);
  dEdx_vs_runID_Canvas->SetBottomMargin(0.14);
  tree_->Draw((estimator + ":runID>>dEdx_vs_runID_Histo" + estimator).c_str(),
              "pT>15");
  dEdx_vs_runID_Histo->Draw();
  TCanvas *dEdx_vs_runID_Profile_Canvas =
      new TCanvas(("dEdx_vs_runID_Profile_Canvas" + estimator).c_str(),
                  "dEdx_vs_runID_Profile_Canvas", 1020, 520);
  dEdx_vs_runID_Profile_Canvas->SetFillColor(0);
  dEdx_vs_runID_Profile_Canvas->SetBorderMode(0);
  dEdx_vs_runID_Profile_Canvas->SetLeftMargin(0.14);
  dEdx_vs_runID_Profile_Canvas->SetRightMargin(0.16);
  dEdx_vs_runID_Profile_Canvas->SetBottomMargin(0.14);
  TProfile *dEdx_vs_runID_Profile = dEdx_vs_runID_Histo->ProfileX(
      ("dEdx_vs_runID_Profile" + estimator).c_str(), 0, dEdxBins);
  dEdx_vs_runID_Profile->Draw();
  if (save) {
    dEdx_vs_runID_Canvas->Print(
        (name_ + estimator + "_dEdx_vs_runID.eps").c_str());
    dEdx_vs_runID_Profile_Canvas->Print(
        (name_ + estimator + "_dEdx_vs_runID_Profile.eps").c_str());
  }
}

void dEdxForOneROC(bool save, unsigned int Id, int bins, double dEdxMin,
                   double dEdxMax) {
  TH1D *dEdxForOneROC =
      new TH1D("dEdxForOneROC", "dE/dx for one ROC", bins, dEdxMin, dEdxMax);

  TLine *ROCHarmonic2Line;

  TTreeReader reader("demo/Tree", file_);
  TTreeReaderValue<doubleVector> dEdxPerPixelRV(reader, "dEdxPerPixel");
  TTreeReaderValue<std::vector<unsigned int>> ROCRV(reader, "ROCId");
  TTreeReaderValue<std::vector<unsigned int>> pixelRawIdRV(reader,
                                                           "pixelRawId");

  double ROCMean = 0;
  double ROCHarmonic2 = 0;
  unsigned int ROCCounter = 0;

  while (reader.Next()) {
    std::vector<double> dEdxPerPixelVector = *dEdxPerPixelRV;
    std::vector<unsigned int> ROCVector = *ROCRV;
    std::vector<unsigned int> pixelRawIdVector = *pixelRawIdRV;
    for (unsigned int i = 0; i < dEdxPerPixelVector.size(); i++) {
      if ((pixelRawIdVector[i] << 3 | ROCVector[i]) == Id) { // 2417004741
        dEdxForOneROC->Fill(dEdxPerPixelVector[i]);
        ROCCounter++;
        ROCMean += dEdxPerPixelVector[i];
        ROCHarmonic2 += pow(dEdxPerPixelVector[i], -2.0);
      }
    }
  }
  std::cout << "ROCCounter: " << ROCCounter
            << ", ROCMean: " << ROCMean / ROCCounter
            << ", ROCHarmonic2: " << pow(ROCHarmonic2 / ROCCounter, -0.5)
            << std::endl;
  ROCHarmonic2Line = new TLine(pow(ROCHarmonic2 / ROCCounter, -0.5), 0,
                               pow(ROCHarmonic2 / ROCCounter, -0.5),
                               dEdxForOneROC->GetMaximum());
  ROCHarmonic2Line->SetLineColor(kBlue);
  TCanvas *dEdxForOneROC_Canvas =
      new TCanvas("dEdxForOneROC_Canvas", "dEdxForOneROC_Canvas", 1020, 520);
  dEdxForOneROC_Canvas->SetFillColor(kWhite);
  dEdxForOneROC_Canvas->SetBorderMode(0);
  dEdxForOneROC_Canvas->SetLeftMargin(0.14);
  dEdxForOneROC_Canvas->SetRightMargin(0.16);
  dEdxForOneROC_Canvas->SetBottomMargin(0.14);
  dEdxForOneROC->Draw();
  ROCHarmonic2Line->Draw("same");
  if (save)
    dEdxForOneROC_Canvas->Print((name_ + "dEdxForOneROC.eps").c_str());
}
void dEdxForOneAPV(bool save, unsigned int Id, int bins, double dEdxMin,
                   double dEdxMax) {
  TH1D *dEdxForOneAPV =
      new TH1D("dEdxForOneAPV", "dE/dx for one APV", bins, dEdxMin, dEdxMax);

  TLine *APVHarmonic2Line;

  TTreeReader reader("demo/Tree", file_);
  TTreeReaderValue<doubleVector> dEdxPerStripRV(reader, "dEdxPerStrip");
  TTreeReaderValue<std::vector<unsigned int>> APVRV(reader, "APVId");
  TTreeReaderValue<std::vector<unsigned int>> stripRawIdRV(reader,
                                                           "stripRawId");

  double APVMean = 0;
  double APVHarmonic2 = 0;
  unsigned int APVCounter = 0;

  while (reader.Next()) {
    std::vector<double> dEdxPerStripVector = *dEdxPerStripRV;
    std::vector<unsigned int> APVVector = *APVRV;
    std::vector<unsigned int> stripRawIdVector = *stripRawIdRV;
    for (unsigned int i = 0; i < dEdxPerStripVector.size(); i++) {
      if ((stripRawIdVector[i] << 3 | APVVector[i]) == Id) { // 3762756739
        dEdxForOneAPV->Fill(dEdxPerStripVector[i]);
        APVCounter++;
        APVMean += dEdxPerStripVector[i];
        APVHarmonic2 += pow(dEdxPerStripVector[i], -2.0);
      }
    }
  }
  std::cout << "APVCounter: " << APVCounter
            << ", APVMean: " << APVMean / APVCounter
            << ", APVHarmonic2: " << pow(APVHarmonic2 / APVCounter, -0.5)
            << std::endl;
  APVHarmonic2Line = new TLine(pow(APVHarmonic2 / APVCounter, -0.5), 0,
                               pow(APVHarmonic2 / APVCounter, -0.5),
                               dEdxForOneAPV->GetMaximum());
  APVHarmonic2Line->SetLineColor(kBlue);
  TCanvas *dEdxForOneAPV_Canvas =
      new TCanvas("dEdxForOneAPV_Canvas", "dEdxForOneAPV_Canvas", 1020, 520);
  dEdxForOneAPV_Canvas->SetFillColor(kWhite);
  dEdxForOneAPV_Canvas->SetBorderMode(0);
  dEdxForOneAPV_Canvas->SetLeftMargin(0.14);
  dEdxForOneAPV_Canvas->SetRightMargin(0.16);
  dEdxForOneAPV_Canvas->SetBottomMargin(0.14);
  dEdxForOneAPV->Draw();
  APVHarmonic2Line->Draw("same");
  if (save)
    dEdxForOneAPV_Canvas->Print((name_ + "dEdxForOneAPV.eps").c_str());
}
