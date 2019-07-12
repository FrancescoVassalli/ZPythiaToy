#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include "Utilities.C"

using namespace std;
using namespace atlashi;

namespace {
	double InTwoPi (double phi) {
		while (phi < 0 || 2*TMath::Pi() <= phi) {
			if (phi < 0) phi += 2*TMath::Pi();
			else phi -= 2*TMath::Pi();
		}
		return phi;
	}

	double DeltaPhi (double phi1, double phi2, const bool sign=0) {
		phi1 = InTwoPi(phi1);
		phi2 = InTwoPi(phi2);
		double dphi = abs(phi1 - phi2);
		while (dphi > TMath::Pi()) dphi = abs (dphi - 2*TMath::Pi());

		if (sign && InTwoPi (phi2 + dphi) == phi1)
			dphi *= -1;

		return dphi;
	}
	unsigned plotcount=0;
}

void singleHister(){
	gStyle->SetOptStat(0);
	TFile *thisFile = new TFile("../plots/zplots.root","RECREATE");

	string name = "../pythiadata/childTaggedHisterffo1.root";
	TChain * t=new TChain("tree");
	t->Add("../pythiadata/childTaggedHisterfff1.root");
	t->Add("../pythiadata/childTaggedHisterfff2.root");

	int code;
	int id1;
	int id2;
	float x1pdf;
	float x2pdf;
	float Q;
	bool isValence1;
	bool isValence2;

	int z_n, jet_n, part_n;
	vector<float>* z_pt = nullptr, *z_eta = nullptr, *z_phi = nullptr, *z_m = nullptr;
	vector<float>* part_pt = nullptr, *part_eta = nullptr, *part_phi = nullptr;
	std::vector<bool> *part_child=nullptr;

	/*t->SetBranchAddress ("code",  &code);
	t->SetBranchAddress ("id1",   &id1);
	t->SetBranchAddress ("id2",   &id2);
	t->SetBranchAddress ("x1pdf", &x1pdf);
	t->SetBranchAddress ("x2pdf", &x2pdf);
	t->SetBranchAddress ("Q",     &Q);
*/
	t->SetBranchAddress ("z_n",   &z_n);
	t->SetBranchAddress ("z_pt",  &z_pt);
	t->SetBranchAddress ("z_eta", &z_eta);
	t->SetBranchAddress ("z_phi", &z_phi);
	t->SetBranchAddress ("z_m",   &z_m);

	t->SetBranchAddress ("part_n",    &part_n);
	t->SetBranchAddress ("part_pt",   &part_pt);
	t->SetBranchAddress ("part_eta",  &part_eta);
	t->SetBranchAddress ("part_phi",  &part_phi);
	t->SetBranchAddress ("part_child", &part_child);

	TH1F* pt = new TH1F("pt","",7,logspace(2,65,7));
	TH1F* dphi = new TH1F("dphi","",10,0,TMath::Pi());
	pt->Sumw2();
	dphi->Sumw2();
	unsigned nchild=0;
	for (int iEvt = 0; iEvt < t->GetEntries(); iEvt++)
	{
		t->GetEntry (iEvt);
		for (int i = 0; i < part_child->size(); ++i)
		{
			pt->Fill(part_pt->at(i));
			dphi->Fill(DeltaPhi(part_phi->at(i),z_phi->at(0)));
		}
	}
	TCanvas* tc = new TCanvas();
	tc->SetLogy();
	tc->SetLogx();
	pt->Scale(1./t->GetEntries(),"width");
	pt->Draw("e1");

	TCanvas* tc2 = new TCanvas();
	tc2->SetLogy();
	//tc2->SetLogx();
	dphi->Scale(1./t->GetEntries(),"width");
	dphi->Draw("e1");

	thisFile->Write();
	//thisFile->Close();
}
