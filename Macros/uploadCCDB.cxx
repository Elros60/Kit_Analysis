// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file dqFlowAccWeights.C
/// \brief A simple macro to read produced accptance weights for Q-vector and
/// submit them to CCDB
///
/// \author Chi ZHANG, CEA-Saclay, chi.zhang@cern.ch

#include <cmath>
#include <gsl/span>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "CCDB/BasicCCDBManager.h"
#include "Framework/Logger.h"
#include <TFile.h>
#include <TProfile.h>
#include <TSystem.h>

using namespace o2;
using namespace std;

vector<string> tokenize(string input_string);

void uploadCCDB(std::string FileName,
                std::string RecoPass = "DQ_LHC23_PbPb_pass4_muon",
                std::string EventSel = "eventStandardSel8NoTFBorder") {
  const std::string &ccdbHost = "http://alice-ccdb.cern.ch";
  const std::string &objectPathSP =
      "Users/c/chizh/FlowResolution/ScalarProduct";
  const std::string &objectPathEP = "Users/c/chizh/FlowResolution/EventPlane";

  // Load weights
  TFile *file;
  try {
    file = TFile::Open(FileName.c_str());
  } catch (std::exception const &e) {
    LOG(fatal) << "Cannot open input file!";
  }

  // Send to CCDB
  if (!ccdbHost.empty()) {
    o2::ccdb::CcdbApi api;
    api.init(ccdbHost.c_str());
    TList *slist_sp = (TList *)file->Get("RunByRun_HistogramFromProfile_SP");
    TList *slist_ep = (TList *)file->Get("RunByRun_HistogramFromProfile_EP");

    LOGP(info, "Storing alignment object on {}/{}", ccdbHost, objectPathSP);
    for (auto spIt = slist_sp->First(); spIt <= slist_sp->Last(); ++spIt) {
      TH1D *hist = (TH1D *)spIt;
      hist->Print();
      vector<string> info = tokenize(hist->GetName());
      auto soreor =
          o2::ccdb::BasicCCDBManager::getRunDuration(api, stoi(info[0]));
      int64_t runStart = soreor.first;
      int64_t runEnd = soreor.second;
      map<string, string> metadata; // can be empty
      metadata["runNumber"] = info[0];
      metadata["recoPass"] = RecoPass.c_str();
      metadata["eventSelection"] = EventSel.c_str();
      LOG(info) << "Uploading SP resolution for run: " << info[0]
                << " - Timestamp from: " << runStart << " to: " << runEnd;
      try {
        api.storeAsTFileAny(hist, objectPathSP, metadata, runStart, runEnd);
      } catch (std::exception const &e) {
        LOG(fatal) << "Failed at CCDB submission!";
      }
      delete hist;
    }

    LOGP(info, "Storing alignment object on {}/{}", ccdbHost, objectPathEP);
    for (auto epIt = slist_ep->First(); epIt <= slist_ep->Last(); ++epIt) {
      TH1D *hist = (TH1D *)epIt;
      hist->Print();
      vector<string> info = tokenize(hist->GetName());
      auto soreor =
          o2::ccdb::BasicCCDBManager::getRunDuration(api, stoi(info[0]));
      // stoi(info[0])
      int64_t runStart = soreor.first;
      int64_t runEnd = soreor.second;
      map<string, string> metadata; // can be empty
      metadata["runNumber"] = info[0];
      metadata["recoPass"] = RecoPass.c_str();
      metadata["eventSelection"] = EventSel.c_str();
      LOG(info) << "Uploading EP resolution for run: " << info[0]
                << " - Timestamp from: " << runStart << " to: " << runEnd;
      try {
        api.storeAsTFileAny(hist, objectPathEP, metadata, runStart, runEnd);
      } catch (std::exception const &e) {
        LOG(fatal) << "Failed at CCDB submission!";
      }
      delete hist;
    }
  }
  file->Close();
}

vector<string> tokenize(string input_string) {
  vector<string> output_string;
  string temp;
  istringstream ss(input_string);
  while (getline(ss, temp, '_')) {
    output_string.push_back(temp);
  }
  return output_string;
}