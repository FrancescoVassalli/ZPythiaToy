#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

void myhist(){
	TCanvas *tc = new TCanvas();
	gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    TLegend *tl=new TLegend(.4,.2,.6,.4);
    const int nBins = 53;
    double bins[nBins+1]= {-2.65,-2.55,-2.45,-2.35,-2.25,-2.15,-2.05,-1.95,-1.85,-1.75,-1.65,-1.55,-1.45,-1.35,-1.25,-1.15,-1.05,-0.95,-0.85,-0.75,-0.65,-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95,2.05,2.15,2.25,2.35,2.45,2.55,2.65};
    for(int i=0; i<=nBins;i++){
    	bins[i] = bins[i]-.05; //accidentally put the bin middles in instead of the bin lows
    }

    double errors[nBins+1] = {0.098,0.087,0.088,0.093,0.093,0.096,0.096,0.097,0.097,0.097,0.097,0.097,0.096,0.096,0.095,0.095,0.095,0.094,0.094,0.093,0.093,0.093,0.089,0.089,0.09,0.09,0.09,0.089,0.09,0.091,0.091,0.092,0.092,0.092,0.092,0.093,0.094,0.093,0.094,0.094,0.094,0.094,0.096,0.097,0.097,0.097,0.097,0.098,0.098,0.095,0.094,0.098,0.099,0.172};
    double data[nBins+1]= {8.681,8.556,8.791,8.889,8.937,8.969,8.993,9.109,9.059,9.093,9.093,9.03,9.016,8.975,8.898,8.868,8.825,8.698,8.672,8.585,8.524,8.456,8.375,8.308,8.325,8.285,8.253,8.214,8.284,8.358,8.415,8.473,8.609,8.686,8.723,8.839,8.955,8.993,9.089,9.212,9.221,9.265,9.34,9.41,9.405,9.405,9.433,9.444,9.372,9.489,9.474,9.34,9.364,9.23};
    double sys[nBins+1]= {0.637,0.622,0.618,0.615,0.611,0.607,0.6,0.596,0.587,0.582,0.577,0.566,0.555,0.548,0.54,0.535,0.527,0.515,0.511,0.503,0.494,0.49,0.483,0.479,0.477,0.473,0.471,0.467,0.473,0.479,0.486,0.492,0.503,0.509,0.514,0.526,0.537,0.543,0.553,0.569,0.57,0.583,0.593,0.605,0.614,0.625,0.635,0.653,0.658,0.699,0.701,0.706,0.713,0.727};
    TH1F *h_f = new TH1F("hf","Eta distribution 5.06 TeV charged final state 60-90 \%",nBins,bins);
    TH1F *h_fe = new TH1F("hfsys","Eta distribution 5.06 TeV charged final state 60-90 \%",nBins,bins);
    h_f->SetXTitle("Eta");
    h_f->SetYTitle("norm count");
    for(int i=0; i<=nBins;i++){
    	h_f->SetBinContent(i,data[i]);
    	h_f->SetBinError(i,errors[i]);
    	h_fe->SetBinContent(i,data[i]);
    	h_fe->SetBinError(i,sys[i]);
    }
    h_f->Sumw2();
    h_fe->Sumw2();
    //h_f->Scale(1/h_f->Integral(),"width");
    //h_fe->Scale(1/h_fe->Integral(),"width");
    h_f->SetMarkerStyle(20);
    h_f->SetLineWidth(3);
    h_fe->SetMarkerStyle(20);
    h_fe->SetLineWidth(9);
    h_fe->SetLineColor(kRed);
    tl->AddEntry(h_f, "Atlas p-Pb 60-90 \%","p");
    h_fe->Draw("P same");
    h_f->Draw("Pe same");
    tc->Print("atlas6090.pdf");


}