#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include "Utilities.C"

using namespace std;
using namespace atlashi;

double DeltaPhi (double phi1, double phi2, const bool sign=0) {
  phi1 = InTwoPi(phi1);
  phi2 = InTwoPi(phi2);
  double dphi = abs(phi1 - phi2);
  while (dphi > TMath::Pi()) dphi = abs (dphi - 2*TMath::Pi());

  if (sign && InTwoPi (phi2 + dphi) == phi1)
     dphi *= -1;

  return dphi;
}

void minBiasHister(){
	gStyle->SetOptStat(0);
	TFile* f = new TFile("../pythiadata/minbias.root", "READ");
	TTree* t = (TTree*) f->Get("tree");
	TFile *thisFile = new TFile("hists.root","UPDATE");

	vector<float>* part_pt = nullptr, *part_eta = nullptr, *part_phi = nullptr;

	t->SetBranchAddress ("part_pt",   &part_pt);
	t->SetBranchAddress ("part_eta",  &part_eta);
	t->SetBranchAddress ("part_phi",  &part_phi);

	TH1F* production = new TH1F("minbias","",7,logspace(2,65,7));
	production->Sumw2();

	//make npart
	for (int iEvt = 0; iEvt < t->GetEntries(); iEvt++) {
		t->GetEntry (iEvt);
		for (unsigned i=0; i < part_pt->size(); i++) {
			production->Fill(part_pt->at(i));
		}
	}

	production->Scale(1./t->GetEntries(),"width");
	production->Scale(1./TMath::Pi());

	thisFile->Write();
	thisFile->Close();
}
