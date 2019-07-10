#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

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


void childTaggedHister(){
	gStyle->SetOptStat(0);
	TChain* t = new TChain("tree");
	t->Add("ChildTagged1.root");
	t->Add("ChildTagged2.root");
	//t->Add("testout.root");
	TFile *thisFile = new TFile("zplots.root","RECREATE");

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

	t->SetBranchAddress ("code",  &code);
	t->SetBranchAddress ("id1",   &id1);
	t->SetBranchAddress ("id2",   &id2);
	t->SetBranchAddress ("x1pdf", &x1pdf);
	t->SetBranchAddress ("x2pdf", &x2pdf);
	t->SetBranchAddress ("Q",     &Q);

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

	std::vector<TH1F*> ntrack_plots;
	ntrack_plots.push_back(new TH1F("central","",15,0,45));
	ntrack_plots.push_back(new TH1F("medium","",15,0,45));
	ntrack_plots.push_back(new TH1F("outer","",15,0,45));

	std::vector<TH1F*> dphi_plots;
	dphi_plots.push_back(new TH1F("children","",15,0,TMath::Pi()));
	dphi_plots.push_back(new TH1F("strangers","",15,0,TMath::Pi()));

	unsigned totalZ=0;

	for (int iEvt = 0; iEvt < t->GetEntries(); iEvt++) {
		t->GetEntry (iEvt);
		if(z_n!=1)continue;
		totalZ+=z_n;
		for (unsigned i=0; i < part_pt->size(); i++) {
			if(DeltaPhi(part_phi->at(i),z_phi->at(0)) <TMath::Pi()/2&&DeltaPhi(part_phi->at(i),z_phi->at(0))>0){
				ntrack_plots[0]->Fill(part_pt->at(i));
			}
			else if(DeltaPhi(part_phi->at(i),z_phi->at(0))>3*TMath::Pi()/4&&DeltaPhi(part_phi->at(i),z_phi->at(0))<15*TMath::Pi()/16){
				ntrack_plots[1]->Fill(part_pt->at(i));
			}
			else if(DeltaPhi(part_phi->at(i),z_phi->at(0))<TMath::Pi()&&DeltaPhi(part_phi->at(i),z_phi->at(0))>15*TMath::Pi()/16){
				ntrack_plots[2]->Fill(part_pt->at(i));
			}
			if (part_child->at(i))
			{
				dphi_plots[0]->Fill(DeltaPhi(part_phi->at(i),z_phi->at(0)));
			}
			else{
				dphi_plots[1]->Fill(DeltaPhi(part_phi->at(i),z_phi->at(0)));
			}
		}
	}
	TCanvas* tc = new TCanvas();
	tc->SetLogy();
	tc->SetLogx();
	unsigned count=0;
	short colors[3]={kBlack,kRed,kBlue};
	TLegend* tl = new TLegend(.2,.1,.4,.4);
	for (std::vector<TH1F*>::iterator i = ntrack_plots.begin(); i != ntrack_plots.end(); ++i)
	{
		(*i)->Scale(1./totalZ);
		(*i)->GetYaxis()->SetRangeUser(10e-7,10e1);
		(*i)->SetLineColor(colors[count]);
		if (count++==0)(*i)->Draw();
		else (*i)->Draw("same");
		tl->AddEntry((*i),(*i)->GetName(),"l");
	}
	tl->Draw();

	TCanvas* tc2 = new TCanvas();
	tc2->SetLogy();
	count=0;
	TLegend* tl2 = new TLegend(.2,.1,.4,.4);
	for (std::vector<TH1F*>::iterator i = dphi_plots.begin(); i != dphi_plots.end(); ++i)
	{
		(*i)->Scale(1./totalZ);
		(*i)->SetLineColor(colors[count]);
		if (count++==0)(*i)->Draw();
		else (*i)->Draw("same");
		tl2->AddEntry((*i),(*i)->GetName(),"l");
	}
	tl2->Draw();

	TCanvas* tc3 = new TCanvas();
	TH1F *diff = (TH1F*) dphi_plots[0]->Clone("dphi_diff");
	diff->Add(dphi_plots[1],-1);
	diff->GetYaxis()->SetRangeUser(-1,1);
	diff->Draw();

	thisFile->Write();
	//thisFile->Close();
}
