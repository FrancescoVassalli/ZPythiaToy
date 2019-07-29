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
	unsigned plotcount=0;
}

void childTaggedHister(){
	gStyle->SetOptStat(0);
	TFile *thisFile = new TFile("hists.root","UPDATE");

	string name = "../pythiadata/";
	string extention = ".root";
	std::vector<string> options;
	options.push_back("mpioff_inclusive1");
	options.push_back("jet_mpi_inclusive1");
	options.push_back("jet_mpi_inclusive2");
	options.push_back("jet_mpi_inclusive3");
	options.push_back("jet_mpi_inclusive4");
	options.push_back("fff_inclusive");

	std::vector<TChain*> chains;

	for (unsigned i=0; i<options.size();++i)
	{
		chains.push_back(new TChain("tree"));
		/*string name1 = name+options[i]+"1"+extention;
		string name2 = name+options[i]+"2"+extention;
		chains[i]->Add(name1.c_str());
		chains[i]->Add(name2.c_str());*/
		string name1 = name+options[i]+extention;
		chains[i]->Add(name1.c_str());
	}
	std::vector<string>::iterator nameit=options.begin();
	for (std::vector<TChain*>::iterator chainpointer = chains.begin(); chainpointer != chains.end(); ++chainpointer)
	{
		TChain* t=*chainpointer;

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

		ntrack_plots.push_back(new TH1F((*nameit+"_#Delta#phi#in[0,#frac{#pi}{2}]").c_str(),"",7,logspace(2,65,10)));
		ntrack_plots.push_back(new TH1F((*nameit+"_#Delta#phi#in[#frac{3#pi}{4},#frac{15#pi}{16}]").c_str(),"",7,logspace(2,65,10)));
		ntrack_plots.push_back(new TH1F((*nameit+"_#Delta#phi#in[#frac{15#pi}{16},#pi]").c_str(),"",7,logspace(2,65,10)));
		for (std::vector<TH1F*>::iterator i = ntrack_plots.begin(); i != ntrack_plots.end(); ++i)
		{
			(*i)->Sumw2();
		}

		std::vector<TH1F*> ntrackChild_plots;
		ntrackChild_plots.push_back(new TH1F((*nameit+"_child_ntrack").c_str(),"",7,logspace(2,65,7)));
		ntrackChild_plots.push_back(new TH1F((*nameit+"_forgein_ntrack").c_str(),"",7,logspace(2,65,7)));
		for (std::vector<TH1F*>::iterator i = ntrackChild_plots.begin(); i != ntrackChild_plots.end(); ++i)
		{
			(*i)->Sumw2();
		}

		std::vector<TH1F*> dphi_child_plots;
		dphi_child_plots.push_back(new TH1F((*nameit+"_child_dphi").c_str(),"",10,0,TMath::Pi()));
		dphi_child_plots.push_back(new TH1F((*nameit+"_forgein_dphi").c_str(),"",10,0,TMath::Pi()));
		for (std::vector<TH1F*>::iterator i = dphi_child_plots.begin(); i != dphi_child_plots.end(); ++i)
		{
			(*i)->Sumw2();
		}

		pair<float,float> lowptRange(2,3.3);
		pair<float,float> midptRange(5.4,8.9);
		pair<float,float> highptRange(24,39.5);
		std::vector<TH1F*> dphi_pt_plots;
		dphi_pt_plots.push_back(new TH1F((*nameit+"_lowpt").c_str(),"",10,0,TMath::Pi()));
		dphi_pt_plots.push_back(new TH1F((*nameit+"_midpt").c_str(),"",10,0,TMath::Pi()));
		dphi_pt_plots.push_back(new TH1F((*nameit+"_highpt").c_str(),"",10,0,TMath::Pi()));
		for (std::vector<TH1F*>::iterator i = dphi_pt_plots.begin(); i != dphi_pt_plots.end(); ++i)
		{
			(*i)->Sumw2();
		}

		unsigned totalZ=0;
		unsigned totalChildren=0;
		for (int iEvt = 0; iEvt < t->GetEntries(); iEvt++) {
			t->GetEntry (iEvt);
			if(z_n!=0&&z_pt->at(0)<25)continue;
			totalZ+=z_n;
			for (unsigned i=0; i < part_pt->size(); i++) {
				//fill ntrack by deltaphi plots
				if(DeltaPhi(part_phi->at(i),z_phi->at(0)) <TMath::Pi()/2&&DeltaPhi(part_phi->at(i),z_phi->at(0))>0){
					ntrack_plots[0]->Fill(part_pt->at(i));
				}
				else if(DeltaPhi(part_phi->at(i),z_phi->at(0))>3*TMath::Pi()/4&&DeltaPhi(part_phi->at(i),z_phi->at(0))<15*TMath::Pi()/16){
					ntrack_plots[1]->Fill(part_pt->at(i));
				}
				else if(DeltaPhi(part_phi->at(i),z_phi->at(0))<TMath::Pi()&&DeltaPhi(part_phi->at(i),z_phi->at(0))>15*TMath::Pi()/16){
					ntrack_plots[2]->Fill(part_pt->at(i));
				}
				//fill child based plots
				if (part_child->at(i))
				{
					dphi_child_plots[0]->Fill(DeltaPhi(part_phi->at(i),z_phi->at(0)));
					ntrackChild_plots[0]->Fill(part_pt->at(i));
					totalChildren++;
				}
				else{
					ntrackChild_plots[1]->Fill(part_pt->at(i));
					dphi_child_plots[1]->Fill(DeltaPhi(part_phi->at(i),z_phi->at(0)));
				}
				//fill dphi by pt plots 
				if (part_pt->at(i)>lowptRange.first&&part_pt->at(i)<lowptRange.second)
				{
					dphi_pt_plots[0]->Fill(DeltaPhi(part_phi->at(i),z_phi->at(0)));
				}
				if (part_pt->at(i)>midptRange.first&&part_pt->at(i)<midptRange.second)
				{
					dphi_pt_plots[1]->Fill(DeltaPhi(part_phi->at(i),z_phi->at(0)));
				}
				if (part_pt->at(i)>highptRange.first&&part_pt->at(i)<highptRange.second)
				{
					dphi_pt_plots[2]->Fill(DeltaPhi(part_phi->at(i),z_phi->at(0)));
				}
			}
		}

		cout<<"children per Z="<<(double) totalChildren/totalZ<<'\n';
		//phi bins
		double bins[3]={2/TMath::Pi(),16./(3*TMath::Pi()),16./TMath::Pi()};
		unsigned count=0;
		for (std::vector<TH1F*>::iterator i = ntrack_plots.begin(); i != ntrack_plots.end(); ++i)
		{
			(*i)->Scale(bins[count]);
			(*i)->Scale(1./totalZ,"width");
			count++;
		}

		for (std::vector<TH1F*>::iterator i = ntrackChild_plots.begin(); i != ntrackChild_plots.end(); ++i)
		{
			(*i)->Scale(1./totalZ,"width");
		}

		for (std::vector<TH1F*>::iterator i = dphi_child_plots.begin(); i != dphi_child_plots.end(); ++i)
		{
			(*i)->Scale(1./totalZ,"width");
		}

		bins[0]=lowptRange.second-lowptRange.first;
		bins[1]=midptRange.second-midptRange.first;
		bins[2]=highptRange.second-highptRange.first;
		count=0;
		for (std::vector<TH1F*>::iterator i = dphi_pt_plots.begin(); i != dphi_pt_plots.end(); ++i)
		{
			(*i)->Scale(1./totalZ,"width");
			//(*i)->Scale(bins[count++]);
		}

		thisFile->Write();
		nameit++;
	}
	

	thisFile->Close();
}


/*
//plot the npart with dphi groups
TCanvas* tc = new TCanvas();
tc->SetLogy();
tc->SetLogx();
unsigned count=0;
short colors[3]={kBlack,kRed,kBlue};
double bins[3]={2/TMath::Pi(),16./(3*TMath::Pi()),16./TMath::Pi()};
TLegend* tl = new TLegend(.8,.7,.9,.95);
for (std::vector<TH1F*>::iterator i = ntrack_plots.begin(); i != ntrack_plots.end(); ++i)
{
(*i)->Scale(1./totalZ,"width");
(*i)->Scale(bins[count]);
(*i)->GetYaxis()->SetRangeUser(10e-8,10e1);
(*i)->SetLineColor(colors[count]);
if (count++==0)(*i)->Draw("e1");
else (*i)->Draw("e1 same");
tl->AddEntry((*i),(*i)->GetName(),"l");
}
tl->Draw();
cout<<"jeff integral ="<<ntrack_plots[0]->Integral("width")+ntrack_plots[1]->Integral("width")+ntrack_plots[2]->Integral("width")<<'\n';
tc->SaveAs(("../plots/"+*nameit+"_npart_dphi.pdf").c_str());

//plot the npart with child groups
count=0;
TCanvas* tcC = new TCanvas();
tcC->SetLogy();
tcC->SetLogx();
TLegend* tlC = new TLegend(.8,.7,.9,.95);
for (std::vector<TH1F*>::iterator i = ntrackChild_plots.begin(); i != ntrackChild_plots.end(); ++i)
{
(*i)->Scale(1./totalZ,"width");
(*i)->GetYaxis()->SetRangeUser(10e-7,10e1);
(*i)->SetLineColor(colors[count]);
if (count++==0)(*i)->Draw("e1");
else (*i)->Draw("e1 same");
tlC->AddEntry((*i),(*i)->GetName(),"l");
}
double error,error2;
cout<<"children integral ="<<ntrackChild_plots[0]->IntegralAndError(1,7,error2,"width")<<"+/"<<error2<<" not="<<ntrackChild_plots[1]->IntegralAndError(1,7,error,"width")<<"+/"<<error<<'\n';
tlC->Draw();
tcC->SaveAs(("../plots/"+*nameit+"_npart_child.pdf").c_str());

//plot the dphi with 
TCanvas* tc2 = new TCanvas();
tc2->SetLogy();
count=0;
TLegend* tl2 = new TLegend(.15,.8,.25,.9);
for (std::vector<TH1F*>::iterator i = dphi_child_plots.begin(); i != dphi_child_plots.end(); ++i)
{
(*i)->Scale(1./totalZ,"width");
(*i)->SetLineColor(colors[count]);
(*i)->GetYaxis()->SetRangeUser(10e-7,10e1);
if (count++==0)(*i)->Draw("e1");
else (*i)->Draw("e1 same");
tl2->AddEntry((*i),(*i)->GetName(),"l");
}
tl2->Draw();
tc2->SaveAs(("../plots/"+*nameit+"_child_dphi.pdf").c_str());

//compare dphi for initial parton children to mpi
TCanvas* tc3 = new TCanvas();
TH1F *diff = (TH1F*) dphi_child_plots[0]->Clone("dphi_diff");
diff->Add(dphi_child_plots[1],-1);
diff->GetYaxis()->SetRangeUser(-10,10);
diff->Draw();
tc3->SaveAs(("../plots/"+*nameit+"_diff_dphi.pdf").c_str());
*/
