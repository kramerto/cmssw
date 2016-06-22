// C/C++ standard libraries
#include <iostream> // for input/output prints
#include <fstream>  // for reading the ASCII files
#include "TH1D.h"
#include "TGraph.h"

using doubleVector = std::vector<double>;

void gains_macro() {

  TFile *file_ =
      new TFile("Demo/DemoAnalyzer/data/Data13TeVGains_v2.root", "READ");

  std::vector<std::string> run = {
      //      "Gains_252116_to_254227", "Gains_254227_to_254437",
      //      "Gains_254437_to_254532", "Gains_254532_to_254790",
      //      "Gains_254790_to_254879", "Gains_254879_to_255031",
      //      "Gains_255031_to_256630", "Gains_256630_to_256734",
      //      "Gains_256734_to_256941", "Gains_256941_to_256957",
      //      "Gains_256957_to_257490", "Gains_257490_to_257531",
      //      "Gains_257531_to_257682", "Gains_257682_to_257823",
      //      "Gains_257823_to_258129", "Gains_258129_to_258174",
      //      "Gains_258174_to_258287", "Gains_258287_to_258702",
      //      "Gains_258702_to_258713", "Gains_258713_to_258750",
      //      "Gains_258750_to_259237", "Gains_259237_to_259352",
      //      "Gains_259352_to_259399", "Gains_259399_to_259626",
      //      "Gains_259626_to_259686", "Gains_259686_to_259809",
      //      "Gains_259809_to_259861", "Gains_259861_to_260373",
      //      "Gains_260069_to_260427", "Gains_260373_to_260069",
      //      "Gains_260427_to_260533", "Gains_260533_to_260577",
      //      "Gains_260577_to_260627",
      "Gains_272760_to_274093", "Gains_274094_to_274197",
      "Gains_274198_to_274249", "Gains_274250_to_274386",
      "Gains_274387_to_274420", "Gains_274421_to_274953",
      "Gains_274954_to_274970"};

  std::vector<double> GainFirst;
  std::vector<double> GainLast;
  std::vector<double> GainInbetween;
  std::vector<double> RunId;
  std::vector<double> DetAPVId;
  std::vector<double> Index;
  std::vector<double> GainsVector;
  std::vector<double> RunIdsVector;
  std::vector<unsigned int> DetAPVIdsVector;
  std::vector<unsigned int> IndexVector;

  for (unsigned int i = 0; i < run.size(); i++) {
    TTreeReader reader(run[i].c_str(), file_);
    TTreeReaderValue<unsigned char> APVRV(reader, "APVId");
    TTreeReaderValue<double> GainRV(reader, "Gain");
    TTreeReaderValue<unsigned int> DetIdRV(reader, "DetId");
    TTreeReaderValue<unsigned int> IndexRV(reader, "Index");

    while (reader.Next()) {
      double GainDouble = *GainRV;
      unsigned char APVChar = *APVRV;
      unsigned int DetIdInt = *DetIdRV;
      unsigned int IndexInt = *IndexRV;
      GainsVector.push_back(GainDouble);
      RunIdsVector.push_back(atoi(run[i].substr(6, 12).c_str()));
      DetAPVIdsVector.push_back(DetIdInt << 3 | APVChar);
      IndexVector.push_back(IndexInt);
    }
    GainFirst.push_back(GainsVector[0 + (i * 88624)]);
    GainLast.push_back(GainsVector[88623 + (i * 88624)]);
    GainInbetween.push_back(GainsVector[30000 + (i * 88624)]);
    RunId.push_back(RunIdsVector[0 + (i * 88624)]);
    DetAPVId.push_back(DetAPVIdsVector[0 + (i * 88624)]);
    Index.push_back(IndexVector[0 + (i * 88624)]);
  }

  TCanvas *Gain_vs_Run_Canvas =
      new TCanvas("Gain_vs_Run_Canvas", "Gain_vs_Run_Canvas", 1020, 520);
  Gain_vs_Run_Canvas->SetFillColor(kWhite);
  Gain_vs_Run_Canvas->SetBorderMode(0);
  Gain_vs_Run_Canvas->SetLeftMargin(0.14);
  Gain_vs_Run_Canvas->SetRightMargin(0.16);
  Gain_vs_Run_Canvas->SetBottomMargin(0.14);

  TGraph *Gain_vs_Run_First =
      new TGraph(GainFirst.size(), &(RunId[0]), &(GainFirst[0]));

  Gain_vs_Run_First->SetMarkerStyle(8);

  TGraph *Gain_vs_Run_Inbetween =
      new TGraph(GainInbetween.size(), &(RunId[0]), &(GainInbetween[0]));

  Gain_vs_Run_Inbetween->SetMarkerStyle(8);
  Gain_vs_Run_Inbetween->SetLineColor(kBlue);
  Gain_vs_Run_Inbetween->SetMarkerColor(kBlue);

  TGraph *Gain_vs_Run_Last =
      new TGraph(GainLast.size(), &(RunId[0]), &(GainLast[0]));
  Gain_vs_Run_Last->SetTitle("Gains for a single APV");
  Gain_vs_Run_Last->GetXaxis()->SetTitle("RunId");
  Gain_vs_Run_Last->GetYaxis()->SetTitle("Gains");
  Gain_vs_Run_Last->SetMarkerStyle(8);
  Gain_vs_Run_Last->SetLineColor(kRed);
  Gain_vs_Run_Last->SetMarkerColor(kRed);
  Gain_vs_Run_Last->GetYaxis()->SetRangeUser(0.7, 1.2);
  Gain_vs_Run_Last->Draw("APL");
  Gain_vs_Run_First->Draw("PL");
  Gain_vs_Run_Inbetween->Draw("PL");

  TCanvas *Delta_Gains_Canvas =
      new TCanvas("Delta_Gains_Canvas", "Delta_Gains_Canvas", 1020, 520);
  Gain_vs_Run_Canvas->SetFillColor(kWhite);
  Gain_vs_Run_Canvas->SetBorderMode(0);
  Gain_vs_Run_Canvas->SetLeftMargin(0.14);
  Gain_vs_Run_Canvas->SetRightMargin(0.16);
  Gain_vs_Run_Canvas->SetBottomMargin(0.14);

  TH1D *Delta_Gains[6];
  for (unsigned int j = 0; j < 6; j++) {
    Delta_Gains[j] = new TH1D(("Delta_Gains" + std::to_string(j)).c_str(),
                              "#DeltaGains for all APVs", 100, -0.2, 0.2);
    Delta_Gains[j]->SetLineColor(j + 1);
    for (unsigned int i = 0; i < /*GainsVector.size() / run.size()*/ 72784;
         i++) {
      Delta_Gains[j]->Fill(GainsVector[i] - GainsVector[i + (j + 1) * 88624]);
    }
    Delta_Gains[j]->Draw("same");
  }
}
