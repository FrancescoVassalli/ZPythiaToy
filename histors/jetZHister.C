#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "Utilities.h"
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;
using namespace atlashi;

void jetZHister(){
	gStyle->SetOptStat(0);
	TChain* t= new TChain("tree");
	t->Add("../pythiadata/jet_mpi_inclusive1.root");
	t->Add("../pythiadata/jet_mpi_inclusive2.root");
	t->Add("../pythiadata/jet_mpi_inclusive3.root");
	t->Add("../pythiadata/jet_mpi_inclusive4.root");
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

	const float myBins[11]={2,2.58,3.28,4.21,5.4,7.3,8.89,11.5,14.62,24.04,39.53};

	std::vector<TH1F*> ntrack_plots;
	ntrack_plots.push_back(new TH1F("p_{Z}^{T}>25_leading_#Delta#phi#in[0,#frac{#pi}{2}]","",10,myBins));
	ntrack_plots.push_back(new TH1F("p_{Z}^{T}>25_leading_#Delta#phi#in[#frac{3#pi}{4},#frac{15#pi}{16}]","",10,myBins));
	ntrack_plots.push_back(new TH1F("p_{Z}^{T}>25_leading_#Delta#phi#in[#frac{15#pi}{16},#pi]","",10,myBins));
	ntrack_plots.push_back(new TH1F("p_{Z}^{T}>25_sub_#Delta#phi#in[0,#frac{#pi}{2}]","",10,myBins));
	ntrack_plots.push_back(new TH1F("p_{Z}^{T}>25_sub_#Delta#phi#in[#frac{3#pi}{4},#frac{15#pi}{16}]","",10,myBins));
	ntrack_plots.push_back(new TH1F("p_{Z}^{T}>25_sub_#Delta#phi#in[#frac{15#pi}{16},#pi]","",10,myBins));
	ntrack_plots.push_back(new TH1F("p_{Z}^{T}>25_other_#Delta#phi#in[0,#frac{#pi}{2}]","",10,myBins));
	ntrack_plots.push_back(new TH1F("p_{Z}^{T}>25_other_#Delta#phi#in[#frac{3#pi}{4},#frac{15#pi}{16}]","",10,myBins));
	ntrack_plots.push_back(new TH1F("p_{Z}^{T}>25_other_#Delta#phi#in[#frac{15#pi}{16},#pi]","",10,myBins));

	std::vector<TH1F*> ntrack_lowpt_plots;
	ntrack_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_leading_#Delta#phi#in[0,#frac{#pi}{2}]","",10,myBins));
	ntrack_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_leading_#Delta#phi#in[#frac{3#pi}{4},#frac{15#pi}{16}]","",10,myBins));
	ntrack_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_leading_#Delta#phi#in[#frac{15#pi}{16},#pi]","",10,myBins));
	ntrack_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_sub_#Delta#phi#in[0,#frac{#pi}{2}]","",10,myBins));
	ntrack_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_sub_#Delta#phi#in[#frac{3#pi}{4},#frac{15#pi}{16}]","",10,myBins));
	ntrack_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_sub_#Delta#phi#in[#frac{15#pi}{16},#pi]","",10,myBins));
	ntrack_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_other_#Delta#phi#in[0,#frac{#pi}{2}]","",10,myBins));
	ntrack_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_other_#Delta#phi#in[#frac{3#pi}{4},#frac{15#pi}{16}]","",10,myBins));
	ntrack_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_other_#Delta#phi#in[#frac{15#pi}{16},#pi]","",10,myBins));

	std::vector<TH1F*> dphi_plots;
	dphi_plots.push_back(new TH1F("10-20","",100,0,2*TMath::Pi()));
	dphi_plots.push_back(new TH1F("20-40","",100,0,2*TMath::Pi()));
	dphi_plots.push_back(new TH1F("40-80","",100,0,2*TMath::Pi()));
	dphi_plots.push_back(new TH1F("80+","",100,0,  2*TMath::Pi()));

	std::vector<TH1F*> dphi_jet_plots;
	dphi_jet_plots.push_back(new TH1F("10-20_jets","",20,0,2*TMath::Pi()));
	dphi_jet_plots.push_back(new TH1F("20-40_jets","",20,0,2*TMath::Pi()));
	dphi_jet_plots.push_back(new TH1F("40-80_jets","",20,0,2*TMath::Pi()));
	dphi_jet_plots.push_back(new TH1F("80+_jets","",20,0,  2*TMath::Pi()));

	std::vector<TH1F*> dphi_lead_plots;
	dphi_lead_plots.push_back(new TH1F("p_{Z}^{T}>25_leading_2-3.3 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_plots.push_back(new TH1F("p_{Z}^{T}>25_leading_3.3-5.4 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_plots.push_back(new TH1F("p_{Z}^{T}>25_leading_5.4-8.9 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_plots.push_back(new TH1F("p_{Z}^{T}>25_leading_8.9-14.6 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_plots.push_back(new TH1F("p_{Z}^{T}>25_leading_14.6-24.0 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_plots.push_back(new TH1F("p_{Z}^{T}>25_leading_24.0-39.5 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_plots.push_back(new TH1F("p_{Z}^{T}>25_leading_39.5-65 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_plots.push_back(new TH1F("p_{Z}^{T}>25_non-lead_2-3.3 p_{T}","",20,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_plots.push_back(new TH1F("p_{Z}^{T}>25_non-lead_3.3-5.4 p_{T}","",20,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_plots.push_back(new TH1F("p_{Z}^{T}>25_non-lead_5.4-8.9 p_{T}","",20,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_plots.push_back(new TH1F("p_{Z}^{T}>25_non-lead_8.9-14.6 p_{T}","",20,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_plots.push_back(new TH1F("p_{Z}^{T}>25_non-lead_14.6-24.0 p_{T}","",20,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_plots.push_back(new TH1F("p_{Z}^{T}>25_non-lead_24.0-39.5 p_{T}","",20,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_plots.push_back(new TH1F("p_{Z}^{T}>25_non-lead_39.5-65 p_{T}","",20,-1*TMath::PiOver2(),3./2*TMath::Pi()));

	std::vector<TH1F*> dphi_lead_lowpt_plots;
	dphi_lead_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_leading_2-3.3 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_leading_3.3-5.4 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_leading_5.4-8.9 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_leading_8.9-14.6 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_leading_14.6-24.0 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_leading_24.0-39.5 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_leading_39.5-65 p_{T}","",40,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_non-lead_2-3.3 p_{T}","",20,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_non-lead_3.3-5.4 p_{T}","",20,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_non-lead_5.4-8.9 p_{T}","",20,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_non-lead_8.9-14.6 p_{T}","",20,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_non-lead_14.6-24.0 p_{T}","",20,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_non-lead_24.0-39.5 p_{T}","",20,-1*TMath::PiOver2(),3./2*TMath::Pi()));
	dphi_lead_lowpt_plots.push_back(new TH1F("15<p_{Z}^{T}<25_non-lead_39.5-65 p_{T}","",20,-1*TMath::PiOver2(),3./2*TMath::Pi()));

	//initialize the errors
	for (unsigned i =0; i<dphi_lead_plots.size(); i++)
	{
		if (i<9)
		{
			ntrack_plots[i]->Sumw2();
			ntrack_lowpt_plots[i]->Sumw2();
			if (i<4)
			{
				dphi_plots[i]->Sumw2();
				dphi_jet_plots[i]->Sumw2();
			}
		}
		dphi_lead_plots[i]->Sumw2();
		dphi_lead_lowpt_plots[i]->Sumw2();
	}

	unsigned zCounts[4]={0,0,0,0};
	unsigned leadZs=0;
	unsigned zCount1525=0;
	const int kMinJetPT=10;
	const int kMinTrackPT=2;
	//
	for (int iEvt = 0; iEvt < t->GetEntries(); iEvt++) {
		t->GetEntry (iEvt);
		if(z_n!=1||z_pt->at(0)<10)continue; //loop over all events with Z>10
		bool counted = false;
		if (z_pt->at(0)<20)
		{
			//fill the track histograms
			for (int i = 0; i < part_pt->size(); ++i)
			{
				if(part_pt->at(i)<=kMinTrackPT) dphi_plots[0]->Fill(InTwoPi(DeltaPhi(part_phi->at(i),z_phi->at(0),true)));
			}
			//fill the jet histograms
			for (unsigned i=0; i<jet_r04_phi->size();++i)
			{
				if(jet_r04_pt->at(i)>kMinJetPT) dphi_jet_plots[0]->Fill(InTwoPi(DeltaPhi(jet_r04_phi->at(i),z_phi->at(0),true)));
			}
			zCounts[0]=zCounts[0]+1;
		}
		else{  //ZpT>=20
			bool leadZcounted=false;

			for (unsigned i=0; i < part_pt->size(); i++) { //track loop
				if(part_pt->at(i)<=kMinTrackPT) continue;
				//Z going direction
				if (DeltaPhi(part_phi->at(i),z_phi->at(0))<TMath::PiOver2())
				{
					if (DeltaR(part_eta->at(i),lead_eta,part_phi->at(i),lead_phi)<.4)
					{
						ntrack_plots[0]->Fill(part_pt->at(i));
					}
					else if(DeltaR(part_eta->at(i),sub_eta,part_phi->at(i),sub_phi)<.4){
						ntrack_plots[3]->Fill(part_pt->at(i));
					}
					else{
						ntrack_plots[6]->Fill(part_pt->at(i));
					}
				}
				//hemisphere
				else if(DeltaPhi(part_phi->at(i),z_phi->at(0))<15*TMath::Pi()/16&&DeltaPhi(part_phi->at(i),z_phi->at(0))>3*TMath::Pi()/4){
					if (DeltaR(part_eta->at(i),lead_eta,part_phi->at(i),lead_phi)<.4)
					{
						ntrack_plots[1]->Fill(part_pt->at(i));
					}
					else if(DeltaR(part_eta->at(i),sub_eta,part_phi->at(i),sub_phi)<.4){
						ntrack_plots[4]->Fill(part_pt->at(i));
					}
					else{
						ntrack_plots[7]->Fill(part_pt->at(i));
					}
				}
				//away side
				else if(DeltaPhi(part_phi->at(i),z_phi->at(0))<TMath::Pi()&&DeltaPhi(part_phi->at(i),z_phi->at(0))>15*TMath::Pi()/16){
					if (DeltaR(part_eta->at(i),lead_eta,part_phi->at(i),lead_phi)<.4)
					{
						ntrack_plots[2]->Fill(part_pt->at(i));
					}
					else if(DeltaR(part_eta->at(i),sub_eta,part_phi->at(i),sub_phi)<.4){
						ntrack_plots[5]->Fill(part_pt->at(i));
					}
					else{
						ntrack_plots[8]->Fill(part_pt->at(i));
					}
				}
				if (z_pt->at(0)>25)
				{
					if(!leadZcounted){
						leadZs++;
						leadZcounted=true;
					} 
					float dphi=InTwoPi(DeltaPhi(part_phi->at(i),z_phi->at(0),true));
					if (dphi>3.*TMath::PiOver2()) dphi-=2*TMath::Pi();
					//in leading jet
					if (DeltaR(part_eta->at(i),lead_eta,part_phi->at(i),lead_phi)<.4)
					{
						if(part_pt->at(i)<3.3)
						{
							dphi_lead_plots[0]->Fill(dphi);
						}
						else if(part_pt->at(i)<5.4){
							dphi_lead_plots[1]->Fill(dphi);
						}
						else if(part_pt->at(i)<8.9){
							dphi_lead_plots[2]->Fill(dphi);
						}
						else if(part_pt->at(i)<14.6){
							dphi_lead_plots[3]->Fill(dphi);
						}
						else if(part_pt->at(i)<24){
							dphi_lead_plots[4]->Fill(dphi);
						}
						else if(part_pt->at(i)<39.5){
							dphi_lead_plots[5]->Fill(dphi);
						}
						else if(part_pt->at(i)<65){
							dphi_lead_plots[6]->Fill(dphi);
						}
					}
					else{
						if(part_pt->at(i)<3.3)
						{
							dphi_lead_plots[7]->Fill(dphi);
						}
						else if(part_pt->at(i)<5.4){
							dphi_lead_plots[8]->Fill(dphi);
						}
						else if(part_pt->at(i)<8.9){
							dphi_lead_plots[9]->Fill(dphi);
						}
						else if(part_pt->at(i)<14.6){
							dphi_lead_plots[10]->Fill(dphi);
						}
						else if(part_pt->at(i)<24){
							dphi_lead_plots[11]->Fill(dphi);
						}
						else if(part_pt->at(i)<39.5){
							dphi_lead_plots[12]->Fill(dphi);
						}
						else if(part_pt->at(i)<65){
							dphi_lead_plots[13]->Fill(dphi);
						}
					}
					if (z_pt->at(0)<40)
					{
						dphi_plots[1]->Fill(InTwoPi(DeltaPhi(part_phi->at(i),z_phi->at(0),true)));
						if(!counted) {
							zCounts[1]=zCounts[1]+1;
							counted = true;
						}
					}
					else if(z_pt->at(0)<80){
						dphi_plots[2]->Fill(InTwoPi(DeltaPhi(part_phi->at(i),z_phi->at(0),true)));
						if(!counted) {
							zCounts[2]=zCounts[2]+1;
							counted=true;
						}
					}
					else{
						dphi_plots[3]->Fill(InTwoPi(DeltaPhi(part_phi->at(i),z_phi->at(0),true)));
						if(!counted) {
							zCounts[3]=zCounts[3]+1;
							counted=true;
						}
					}
				}//pT Z >25
			}//track loop
			//fill jets for ZpT>20
			if (z_pt->at(0)<40)for (unsigned i=0; i<jet_r04_phi->size();++i)
			{
				if(jet_r04_pt->at(i)>kMinJetPT) dphi_jet_plots[1]->Fill(InTwoPi(DeltaPhi(jet_r04_phi->at(i),z_phi->at(0),true)));
			}
			else if(z_pt->at(0)<80)for (unsigned i=0; i<jet_r04_phi->size();++i)
			{
				if(jet_r04_pt->at(i)>kMinJetPT) dphi_jet_plots[2]->Fill(InTwoPi(DeltaPhi(jet_r04_phi->at(i),z_phi->at(0),true)));
			}
			else for (unsigned i=0; i<jet_r04_phi->size();++i)
			{
				if(jet_r04_pt->at(i)>kMinJetPT) dphi_jet_plots[3]->Fill(InTwoPi(DeltaPhi(jet_r04_phi->at(i),z_phi->at(0),true)));
			}
		}// Z is above 20 GeV

		//15<ZpT<25
		if(z_pt->at(0)>15&&z_pt->at(0)<25) for (unsigned i=0; i < part_pt->size(); i++)
		{
			if(part_pt->at(i)<kMinTrackPT) continue;
			zCount1525++;
			float dphi=InTwoPi(DeltaPhi(part_phi->at(i),z_phi->at(0),true));
			if (dphi>3.*TMath::PiOver2()) dphi-=2*TMath::Pi();
			//in leading jet
			if (DeltaR(part_eta->at(i),lead_eta,part_phi->at(i),lead_phi)<.4)
			{
				if(part_pt->at(i)<3.3)
				{
					dphi_lead_lowpt_plots[0]->Fill(dphi);
				}
				else if(part_pt->at(i)<5.4){
					dphi_lead_lowpt_plots[1]->Fill(dphi);
				}
				else if(part_pt->at(i)<8.9){
					dphi_lead_lowpt_plots[2]->Fill(dphi);
				}
				else if(part_pt->at(i)<14.6){
					dphi_lead_lowpt_plots[3]->Fill(dphi);
				}
				else if(part_pt->at(i)<24){
					dphi_lead_lowpt_plots[4]->Fill(dphi);
				}
				else if(part_pt->at(i)<39.5){
					dphi_lead_lowpt_plots[5]->Fill(dphi);
				}
				else if(part_pt->at(i)<65){
					dphi_lead_lowpt_plots[6]->Fill(dphi);
				}
			}//leading jet
			else{//non-leading jet
				if(part_pt->at(i)<3.3)
				{
					dphi_lead_lowpt_plots[7]->Fill(dphi);
				}
				else if(part_pt->at(i)<5.4){
					dphi_lead_lowpt_plots[8]->Fill(dphi);
				}
				else if(part_pt->at(i)<8.9){
					dphi_lead_lowpt_plots[9]->Fill(dphi);
				}
				else if(part_pt->at(i)<14.6){
					dphi_lead_lowpt_plots[10]->Fill(dphi);
				}
				else if(part_pt->at(i)<24){
					dphi_lead_lowpt_plots[11]->Fill(dphi);
				}
				else if(part_pt->at(i)<39.5){
					dphi_lead_lowpt_plots[12]->Fill(dphi);
				}
				else if(part_pt->at(i)<65){
					dphi_lead_lowpt_plots[13]->Fill(dphi);
				}
			}//non-leading jet

			if (DeltaPhi(part_phi->at(i),z_phi->at(0))<TMath::PiOver2())
			{
				if (DeltaR(part_eta->at(i),lead_eta,part_phi->at(i),lead_phi)<.4)
				{
					ntrack_lowpt_plots[0]->Fill(part_pt->at(i));
				}
				else if(DeltaR(part_eta->at(i),sub_eta,part_phi->at(i),sub_phi)<.4){
					ntrack_lowpt_plots[3]->Fill(part_pt->at(i));
				}
				else{
					ntrack_lowpt_plots[6]->Fill(part_pt->at(i));
				}
			}
			//hemisphere
			else if(DeltaPhi(part_phi->at(i),z_phi->at(0))<15*TMath::Pi()/16&&DeltaPhi(part_phi->at(i),z_phi->at(0))>3*TMath::Pi()/4){
				if (DeltaR(part_eta->at(i),lead_eta,part_phi->at(i),lead_phi)<.4)
				{
					ntrack_lowpt_plots[1]->Fill(part_pt->at(i));
				}
				else if(DeltaR(part_eta->at(i),sub_eta,part_phi->at(i),sub_phi)<.4){
					ntrack_lowpt_plots[4]->Fill(part_pt->at(i));
				}
				else{
					ntrack_lowpt_plots[7]->Fill(part_pt->at(i));
				}
			}
			//away side
			else if(DeltaPhi(part_phi->at(i),z_phi->at(0))<TMath::Pi()&&DeltaPhi(part_phi->at(i),z_phi->at(0))>15*TMath::Pi()/16){
				if (DeltaR(part_eta->at(i),lead_eta,part_phi->at(i),lead_phi)<.4)
				{
					ntrack_lowpt_plots[2]->Fill(part_pt->at(i));
				}
				else if(DeltaR(part_eta->at(i),sub_eta,part_phi->at(i),sub_phi)<.4){
					ntrack_lowpt_plots[5]->Fill(part_pt->at(i));
				}
				else{
					ntrack_lowpt_plots[8]->Fill(part_pt->at(i));
				}
			}
		}//Z 15-25
	}
	//plot npart 
	double bins[3]={2/TMath::Pi(),16./(3*TMath::Pi()),16./TMath::Pi()};
	double pTbins[7]={1.3,2.1,8.9-5.4,14.6-8.9,24-14.6,39.5-24,65-39.5};

	unsigned totalZ = 0;
	for (unsigned i : zCounts){
		totalZ+=i;
	}
	for (int i = 0; i < dphi_lead_plots.size(); ++i)
	{
		if(i<ntrack_plots.size()){
			ntrack_plots[i]->Scale(1./(totalZ-zCounts[0]),"width");
			ntrack_plots[i]->Scale(bins[i%3]);
			ntrack_lowpt_plots[i]->Scale(1./zCount1525,"width");
			ntrack_lowpt_plots[i]->Scale(bins[i%3]);
		}
		if(i<dphi_plots.size()) {
			dphi_plots[i]->Scale(1./zCounts[i],"width");
			dphi_jet_plots[i]->Scale(1./zCounts[i],"width");
		}
		dphi_lead_plots[i]->Scale(1./leadZs,"width");
		dphi_lead_plots[i]->Scale(pTbins[i%7]);
		dphi_lead_lowpt_plots[i]->Scale(1./zCount1525,"width");
		dphi_lead_lowpt_plots[i]->Scale(pTbins[i%7]);
	}
	thisFile->Write();
	//thisFile->Close();
}
