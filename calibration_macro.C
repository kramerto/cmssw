// C/C++ standard libraries
#include <iostream> // for input/output prints
#include <fstream>  // for reading the ASCII files
#include <unordered_map>

using doubleVector = std::vector<double>;

void calibration_macro() {
  //  TFile *Data13TeVGainsV2_ =
  //      new TFile("Demo/DemoAnalyzer/data/Data13TeVGains_v2.root", "READ");
  //
  //  std::string run_ = "Gains_274250_to_274386";
  //
  //  TTree *calibrationTree_ = (TTree *)Data13TeVGainsV2_->Get(run_.c_str());
  //
  //  std::unordered_map<unsigned int, unsigned int> *IdToIndex =
  //      new std::unordered_map<unsigned int, unsigned int>();
  //  std::map<unsigned int, unsigned int> *IndexToId =
  //      new std::map<unsigned int, unsigned int>();
  //
  //  unsigned int tree_DetId;
  //  calibrationTree_->SetBranchAddress("DetId", &tree_DetId);
  //  unsigned char tree_APVId;
  //  calibrationTree_->SetBranchAddress("APVId", &tree_APVId);
  //  unsigned int tree_Index;
  //  calibrationTree_->SetBranchAddress("Index", &tree_Index);
  //  for (unsigned int ientry = 0; ientry < calibrationTree_->GetEntries();
  //       ientry++) {
  //    calibrationTree_->GetEntry(ientry);
  //    (*IdToIndex)[tree_DetId << 3 | tree_APVId] = tree_Index;
  //    (*IndexToId)[tree_Index] = tree_DetId << 3 | tree_APVId;
  //    if (tree_DetId == 436246001)
  //      std::cout << "DetID 436246001 gefunden, APV: " << (int)tree_APVId
  //                << std::endl;
  //  }
  std::cout << "Opening File..." << std::endl;
  TFile *Run2016BSingleMuon_ =
      new TFile("/nfs/dust/cms/user/kramerto/"
                "Run2016B_SingleMuon_noCalibration_NTuple.root",
                "READ");
  std::cout << "Done" << std::endl;
  //  TTree *dataTree_ = (TTree *)Run2016BSingleMuon_->Get("demo/Tree");
  //  double dEdxPerChipSum[IdToIndex->size()];
  //  double dEdxPerChipHarmonic2[IdToIndex->size()];
  //  unsigned int counter[IdToIndex->size()];
  //  unsigned int rawId[IdToIndex->size()];
  //  unsigned int chipId[IdToIndex->size()];
  //  for (unsigned int s = 0; s < IdToIndex->size(); s++) {
  //    dEdxPerChipSum[s] = 0;
  //    dEdxPerChipHarmonic2[s] = 0;
  //    counter[s] = 0;
  //    rawId[s] = 0;
  //    chipId[s] = 0;
  //  }

  std::unordered_map<unsigned int, double> *dEdxPerChipSum =
      new std::unordered_map<unsigned int, double>();
  dEdxPerChipSum->clear();
  std::unordered_map<unsigned int, double> *dEdxPerChipHarmonic2 =
      new std::unordered_map<unsigned int, double>();
  dEdxPerChipHarmonic2->clear();
  std::unordered_map<unsigned int, unsigned int> *counter =
      new std::unordered_map<unsigned int, unsigned int>();
  counter->clear();
  std::unordered_map<unsigned int, unsigned int> *rawId =
      new std::unordered_map<unsigned int, unsigned int>();
  rawId->clear();
  std::unordered_map<unsigned int, unsigned int> *chipId =
      new std::unordered_map<unsigned int, unsigned int>();
  chipId->clear();
  std::cout << "Initialized and cleared maps" << std::endl;
  TTreeReader reader("demo/Tree", Run2016BSingleMuon_);
  TTreeReaderValue<doubleVector> dEdxPerStripRV(reader, "dEdxPerStrip");
  TTreeReaderValue<doubleVector> dEdxPerPixelRV(reader, "dEdxPerPixel");
  TTreeReaderValue<std::vector<unsigned int>> APVRV(reader, "APVId");
  TTreeReaderValue<std::vector<unsigned int>> ROCRV(reader, "ROCId");
  TTreeReaderValue<std::vector<unsigned int>> stripRawIdRV(reader,
                                                           "stripRawId");
  TTreeReaderValue<std::vector<unsigned int>> pixelRawIdRV(reader,
                                                           "pixelRawId");
  while (reader.Next()) {
    std::vector<double> dEdxPerStripVector = *dEdxPerStripRV;
    std::vector<double> dEdxPerPixelVector = *dEdxPerPixelRV;
    std::vector<unsigned int> APVVector = *APVRV;
    std::vector<unsigned int> ROCVector = *ROCRV;
    std::vector<unsigned int> stripRawIdVector = *stripRawIdRV;
    std::vector<unsigned int> pixelRawIdVector = *pixelRawIdRV;
    //    std::cout << "reader.Next()" << std::endl;
    for (unsigned int i = 0; i < dEdxPerStripVector.size(); i++) {
      //      std::cout << "Strip Loop index: " << i << std::endl;
      if (APVVector[i] != -1) {
        //        if ((*counter)[stripRawIdVector[i] << 3 | APVVector[i]] < 1) {
        //          (*counter)[stripRawIdVector[i] << 3 | APVVector[i]] = 0;
        //          (*dEdxPerChipSum)[stripRawIdVector[i] << 3 | APVVector[i]] =
        //          0;
        //        }
        (*dEdxPerChipSum)[stripRawIdVector[i] << 3 | APVVector[i]] +=
            dEdxPerStripVector[i];
        (*dEdxPerChipHarmonic2)[stripRawIdVector[i] << 3 | APVVector[i]] +=
            pow(dEdxPerStripVector[i], -2.0);
        (*counter)[stripRawIdVector[i] << 3 | APVVector[i]]++;
        (*rawId)[stripRawIdVector[i] << 3 | APVVector[i]] = stripRawIdVector[i];
        (*chipId)[stripRawIdVector[i] << 3 | APVVector[i]] = APVVector[i];
      }
    }
    for (unsigned int k = 0; k < dEdxPerPixelVector.size(); k++) {
      //      std::cout << "Pixel Loop index: " << k << std::endl;
      if (ROCVector[k] != -1) {
        //        if ((*counter)[pixelRawIdVector[k] << 3 | ROCVector[k]] < 1) {
        //          (*counter)[pixelRawIdVector[k] << 3 | ROCVector[k]] = 0;
        //          (*dEdxPerChipSum)[pixelRawIdVector[k] << 3 | ROCVector[k]] =
        //          0;
        //        }
        (*dEdxPerChipSum)[pixelRawIdVector[k] << 3 | ROCVector[k]] +=
            dEdxPerStripVector[k];
        (*dEdxPerChipHarmonic2)[pixelRawIdVector[k] << 3 | ROCVector[k]] +=
            pow(dEdxPerPixelVector[k], -2.0);
        (*counter)[pixelRawIdVector[k] << 3 | ROCVector[k]]++;
        (*rawId)[pixelRawIdVector[k] << 3 | ROCVector[k]] = pixelRawIdVector[k];
        (*chipId)[pixelRawIdVector[k] << 3 | ROCVector[k]] = ROCVector[k];
      }
    }
  }
  std::unordered_map<unsigned int, double> *calibration =
      new std::unordered_map<unsigned int, double>();
  for (auto &kv : (*dEdxPerChipHarmonic2)) {
    std::cout << "Filling calibration map" << std::endl;
    (*calibration)[kv.first] =
        pow(dEdxPerChipHarmonic2->at(kv.first) / counter->at(kv.first), -0.5) /
        3.0;
  }
  std::ofstream file("Demo/DemoAnalyzer/data/MyCalibration.txt");
  for (auto &kv : (*calibration)) {
    std::cout << "Saving calibration map to file" << std::endl;
    file << (std::to_string(kv.first) + " " + std::to_string(kv.second) + "\n")
                .c_str();
  }
  file.close();
}
