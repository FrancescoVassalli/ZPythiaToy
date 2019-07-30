#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "Utilities.h"
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;
using namespace atlashi;

void ZTrackHister(){
	gStyle->SetOptStat(0);
	TChain* t= new TChain("tree");
	t->Add("../pythiadata/jet_mpi_inclusive1.root");
	t->Add("../pythiadata/jet_mpi_inclusive2.root");
	TFile *thisFile = new TFile("hists.root","UPDATE");

	int code;
	int id1;
	int id2;
	float x1pdf;
	float x2pdf;
	float Q;
	bool isValence1;
	bool isValence2;

	int z_n, part_n, lead_n, sub_n;
	vector<float> *z_pt=nullptr, *z_eta=nullptr, *z_phi=nullptr,* z_m=nullptr;
	vector<float> *l_pt=nullptr, *l_eta=nullptr, *l_phi=nullptr, *l_m=nullptr;
	vector<float> *part_pt=nullptr, *part_eta=nullptr, *part_phi=nullptr;
	vector<float> *jet_r04_pt=nullptr, *jet_r04_eta=nullptr, *jet_r04_phi=nullptr;
	std::vector<bool> *part_child=nullptr;
	float lead_pt, lead_eta, lead_phi;
	float sub_pt, sub_eta, sub_phi;

	t->SetBranchAddress ("z_n",   &z_n);
	t->SetBranchAddress ("z_pt",  &z_pt);
	t->SetBranchAddress ("z_eta", &z_eta);
	t->SetBranchAddress ("z_phi", &z_phi);
	t->SetBranchAddress ("z_m",   &z_m);

	t->SetBranchAddress ("part_n",     &part_n);
	t->SetBranchAddress ("part_pt",    &part_pt);
	t->SetBranchAddress ("part_eta",   &part_eta);
	t->SetBranchAddress ("part_phi",   &part_phi);
	t->SetBranchAddress ("part_child", &part_child);

	t->SetBranchAddress ("lead_n",   &lead_n);
	t->SetBranchAddress ("lead_pt",  &lead_pt);
	t->SetBranchAddress ("lead_eta", &lead_eta);
	t->SetBranchAddress ("lead_phi", &lead_phi);

	t->SetBranchAddress ("sublead_n",   &sub_n);
	t->SetBranchAddress ("sublead_pt",  &sub_pt);
	t->SetBranchAddress ("sublead_eta", &sub_eta);
	t->SetBranchAddress ("sublead_phi", &sub_phi);

	t->SetBranchAddress ("jet_r04_pt",  &jet_r04_pt);
	t->SetBranchAddress ("jet_r04_eta", &jet_r04_eta);
	t->SetBranchAddress ("jet_r04_phi", &jet_r04_phi);

	std::vector<TH1F*> dphi_wide_plots;
	dphi_wide_plots.push_back(new TH1F("15-25 p_{T}^{Z}_2-3.3 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_wide_plots.push_back(new TH1F("15-25 p_{T}^{Z}_3.3-5.4 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_wide_plots.push_back(new TH1F("15-25 p_{T}^{Z}_5.4-8.9 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_wide_plots.push_back(new TH1F("15-25 p_{T}^{Z}_8.9-14.6 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_wide_plots.push_back(new TH1F("15-25 p_{T}^{Z}_14.6-24.0 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_wide_plots.push_back(new TH1F("15-25 p_{T}^{Z}_24.0-39.5 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_wide_plots.push_back(new TH1F("15-25 p_{T}^{Z}_39.5-65 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_wide_plots.push_back(new TH1F("25+ p_{T}^{Z}_2-3.3 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_wide_plots.push_back(new TH1F("25+ p_{T}^{Z}_3.3-5.4 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_wide_plots.push_back(new TH1F("25+ p_{T}^{Z}_5.4-8.9 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_wide_plots.push_back(new TH1F("25+ p_{T}^{Z}_8.9-14.6 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_wide_plots.push_back(new TH1F("25+ p_{T}^{Z}_14.6-24.0 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_wide_plots.push_back(new TH1F("25+ p_{T}^{Z}_24.0-39.5 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_wide_plots.push_back(new TH1F("25+ p_{T}^{Z}_39.5-65 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));

	for (std::vector<TH1F*>::iterator i = dphi_wide_plots.begin(); i != dphi_wide_plots.end(); ++i)
	{
		(*i)->Sumw2();
	}

	unsigned wideCounter[2]={0, 0};
	//make npart
	for (int iEvt = 0; iEvt < t->GetEntries(); iEvt++) {
		t->GetEntry (iEvt);
		if(z_n!=1||z_pt->at(0)<15)continue;
		if (z_pt->at(0)<25)
		{
			wideCounter[0]=wideCounter[0]+1;
			for (int i = 0; i < part_phi->size(); ++i)
			{
				if (part_pt->at(i)<2) continue;
				float dphi=InTwoPi(DeltaPhi(part_phi->at(i),z_phi->at(0),true));
				if (dphi>3.*TMath::PiOver2()) dphi-=2*TMath::Pi();
				if(part_pt->at(i)<3.3)
				{
					dphi_wide_plots[0]->Fill(dphi);
				}
				else if(part_pt->at(i)<5.4){
					dphi_wide_plots[1]->Fill(dphi);
				}
				else if(part_pt->at(i)<8.9){
					dphi_wide_plots[2]->Fill(dphi);
				}
				else if(part_pt->at(i)<14.6){
					dphi_wide_plots[3]->Fill(dphi);
				}
				else if(part_pt->at(i)<24){
					dphi_wide_plots[4]->Fill(dphi);
				}
				else if(part_pt->at(i)<39.5){
					dphi_wide_plots[5]->Fill(dphi);
				}
				else if(part_pt->at(i)<65){
					dphi_wide_plots[6]->Fill(dphi);
				}
			}
		}//Z pT <25
		else{
			wideCounter[1]=wideCounter[1]+1;
			for (int i = 0; i < part_phi->size(); ++i)
			{
				if (part_pt->at(i)<2) continue;
				float dphi=InTwoPi(DeltaPhi(part_phi->at(i),z_phi->at(0),true));
				if (dphi>3.*TMath::PiOver2()) dphi-=2*TMath::Pi();
				if(part_pt->at(i)<3.3)
				{
					dphi_wide_plots[7]->Fill(dphi);
				}
				else if(part_pt->at(i)<5.4){
					dphi_wide_plots[8]->Fill(dphi);
				}
				else if(part_pt->at(i)<8.9){
					dphi_wide_plots[9]->Fill(dphi);
				}
				else if(part_pt->at(i)<14.6){
					dphi_wide_plots[10]->Fill(dphi);
				}
				else if(part_pt->at(i)<24){
					dphi_wide_plots[11]->Fill(dphi);
				}
				else if(part_pt->at(i)<39.5){
					dphi_wide_plots[12]->Fill(dphi);
				}
				else if(part_pt->at(i)<65){
					dphi_wide_plots[13]->Fill(dphi);
				}
			}
		}
	}
	//plot npart 
	double bins[7]={1.3,2.1,8.9-5.4,14.6-8.9,24-14.6,39.5-24,65-39.5};
	for (int i = 0; i < 7; ++i)
	{
		dphi_wide_plots[i]->Scale(1./wideCounter[0],"width");
		dphi_wide_plots[i]->Scale(bins[i]);
		dphi_wide_plots[i+7]->Scale(1./wideCounter[1],"width");
		dphi_wide_plots[i+7]->Scale(bins[i]);
	}
	thisFile->Write();
	//thisFile->Close();
}
