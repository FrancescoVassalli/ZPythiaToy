
void plotJetnPart(TFile* thisFile){
	std::vector<TH1F*> ntrack_plots;
	ntrack_plots.push_back((TH1F*) thisFile->Get("near_lead"));
	ntrack_plots.push_back((TH1F*) thisFile->Get("medium_lead"));
	ntrack_plots.push_back((TH1F*) thisFile->Get("away_lead"));
	ntrack_plots.push_back((TH1F*) thisFile->Get("near_sub"));
	ntrack_plots.push_back((TH1F*) thisFile->Get("medium_sub"));
	ntrack_plots.push_back((TH1F*) thisFile->Get("away_sub"));

	unsigned typeCount=0;
	short colors[3]={kBlack,kRed,kBlue};
	short styleTypes[6]={kFullCircle,kFullTriangleUp,kFullStar,kOpenCircle,kOpenTriangleUp,kOpenStar};
	string titles[3]={"nTracks for jet in Z going #Delta#phi#in[0,#frac{#pi}{2}]",
		"nTracks for jet in hemisphere #Delta#phi#in[#frac{3#pi}{4},#frac{15#pi}{16}]",
		"nTracks for jet in away from Z #Delta#phi#in[#frac{15#pi}{16},#pi]"};
	for (int i = 0; i < 3; ++i)
	{
		TCanvas* tc = new TCanvas();
		tc->SetLogy();
		tc->SetTicky();
		TLegend* tl = new TLegend(.7,.7,.9,.9);
		ntrack_plots[i]->SetYTitle("1/N_{evt} dN_{ch}/d#Delta#phi");
		ntrack_plots[i]->SetXTitle("#it{p}_{T} [GeV]");
		ntrack_plots[i]->SetLineColor(colors[i]);
		ntrack_plots[i]->SetMarkerColor(colors[i]);
		ntrack_plots[i]->SetMarkerStyle(styleTypes[i]);
		ntrack_plots[i+3]->SetYTitle("1/N_{evt} dN_{ch}/d#Delta#phi");
		ntrack_plots[i+3]->SetXTitle("#it{p}_{T} [GeV]");
		ntrack_plots[i+3]->SetLineColor(colors[i]);
		ntrack_plots[i+3]->SetMarkerColor(colors[i]);
		ntrack_plots[i+3]->SetMarkerStyle(styleTypes[i+3]);
		ntrack_plots[i]->Draw("e1");
		ntrack_plots[i+3]->Draw("e1 same");
		tl->AddEntry(ntrack_plots[i],"leading jet","p");
		tl->AddEntry(ntrack_plots[i+3],"non-leading jet","p");
		tl->Draw();
		myText(.4,.6,colors[i],titles[i].c_str(),.04);
		string savename = "../plots/Zjet_ntrack_dphi_";
		savename+=ntrack_plots[i]->GetName();
		savename+=".pdf";
		tc->SaveAs(savename.c_str());	
	}
}

void plotJetnPart2(TFile* thisFile){
	std::vector<TH1F*> ntrack_plots;
	ntrack_plots.push_back((TH1F*) thisFile->Get("near_lead"));
	ntrack_plots.push_back((TH1F*) thisFile->Get("medium_lead"));
	ntrack_plots.push_back((TH1F*) thisFile->Get("away_lead"));
	ntrack_plots.push_back((TH1F*) thisFile->Get("near_sub"));
	ntrack_plots.push_back((TH1F*) thisFile->Get("medium_sub"));
	ntrack_plots.push_back((TH1F*) thisFile->Get("away_sub"));

	unsigned typeCount=0;
	short colors[3]={kBlack,kRed,kBlue};
	short styleTypes[6]={kFullCircle,kFullTriangleUp,kFullStar,kOpenCircle,kOpenTriangleUp,kOpenStar};
	string titles[2]={"nTracks for leading jet",
		"nTracks not in leading jet"};
	string saves[2]={"leading","nonleading"};

	for (int i = 0; i < 2; ++i)
	{
		TCanvas* tc = new TCanvas();
		tc->SetLogy();
		tc->SetTicky();
		TLegend* tl = new TLegend(.7,.7,.9,.9);
		ntrack_plots[i]->SetYTitle("1/N_{evt} dN/d#Delta#phi");
		ntrack_plots[i]->SetXTitle("#it{p}_{T} [GeV]");
		ntrack_plots[i]->SetLineColor(colors[0]);
		ntrack_plots[i]->SetMarkerColor(colors[0]);
		ntrack_plots[i]->SetMarkerStyle(styleTypes[0]);
		ntrack_plots[i+1]->SetLineColor(colors[1]);
		ntrack_plots[i+1]->SetMarkerColor(colors[1]);
		ntrack_plots[i+1]->SetMarkerStyle(styleTypes[1]);
		ntrack_plots[i+2]->SetLineColor(colors[2]);
		ntrack_plots[i+2]->SetMarkerColor(colors[2]);
		ntrack_plots[i+2]->SetMarkerStyle(styleTypes[2]);
		ntrack_plots[i]->Draw("e1");
		ntrack_plots[i+1]->Draw("e1 same");
		ntrack_plots[i+2]->Draw("e1 same");
		tl->AddEntry(ntrack_plots[i],"Z-track #Delta#phi#in[0,#frac{#pi}{2}]","p");
		tl->AddEntry(ntrack_plots[i+1],"Z-track #Delta#phi#in[#frac{3#pi}{4},#frac{15#pi}{16}]","p");
		tl->AddEntry(ntrack_plots[i+2],"Z-track #Delta#phi#in[#frac{15#pi}{16},#pi]","p");
		tl->Draw();
		myText(.4,.6,kBlack,titles[i].c_str(),.04);
		tc->SaveAs(saves[i].c_str());	
	}
}

void plotJetDPhiTracks(TFile* thisFile){
	std::vector<TH1F*> dphi_plots;
	dphi_plots.push_back((TH1F*) thisFile->Get("10-20"));
	dphi_plots.push_back((TH1F*) thisFile->Get("20-40"));
	dphi_plots.push_back((TH1F*) thisFile->Get("40-80"));
	dphi_plots.push_back((TH1F*) thisFile->Get("80+"));

	short colors[4]={kRed,kBlue,kGreen+3,kMagenta-5};
	TCanvas* tc = new TCanvas();
	tc->Draw();
	TLegend* tl = new TLegend(1.5,.7,3.5,.9);
	for (int i = 0; i < 4; ++i)
	{
		dphi_plots[i]->SetYTitle("dN / N_{Z}d#Delta#phi");
		dphi_plots[i]->SetXTitle("Z-hadron #Delta#phi");
		dphi_plots[i]->SetLineColor(colors[i]);
		dphi_plots[i]->GetYaxis()->SetRangeUser(0,14);
		if(i==0) dphi_plots[i]->Draw("");
		else dphi_plots[i]->Draw("same");
		string name  = dphi_plots[i]->GetName();
		name+="p_{T}^{Z} [GeV]";
		myText(.2,.6+.05*i,colors[i],name.c_str(),.04);
	}
	tl->Draw();
	tc->SaveAs("../plots/Zjet_ntrack_dphi.pdf");	
}

void plotJetDPhiTracksWide(TFile* thisFile,string zType){
	std::vector<TH1F*> dphi_plots;
	dphi_plots.push_back((TH1F*) thisFile->Get((zType+" #it{p}_{T}^{Z} 2-3.3 #it{p}_{T}").c_str()));
	dphi_plots.push_back((TH1F*) thisFile->Get((zType+" #it{p}_{T}^{Z} 3.3-5.4 #it{p}_{T}").c_str()));
	dphi_plots.push_back((TH1F*) thisFile->Get((zType+" #it{p}_{T}^{Z} 5.4-8.9 #it{p}_{T}").c_str()));
	dphi_plots.push_back((TH1F*) thisFile->Get((zType+" #it{p}_{T}^{Z} 8.9-14.6 #it{p}_{T}").c_str()));
	dphi_plots.push_back((TH1F*) thisFile->Get((zType+" #it{p}_{T}^{Z} 14.6-24.0 #it{p}_{T}").c_str()));
	dphi_plots.push_back((TH1F*) thisFile->Get((zType+" #it{p}_{T}^{Z} 24.0-39.5 #it{p}_{T}").c_str()));
	dphi_plots.push_back((TH1F*) thisFile->Get((zType+" #it{p}_{T}^{Z} 39.5-65 #it{p}_{T}").c_str()));
	
	short colors[7]={kBlack,kRed,kBlue,kGreen+3,kMagenta-5,kMagenta+4,kCyan-4};
	TCanvas* tc = new TCanvas();
	tc->Draw();
	TLegend* tl = new TLegend(1.5,.7,3.5,.9);
	for (int i = 0; i < dphi_plots.size(); ++i)
	{
		dphi_plots[i]->SetYTitle("1/N_{Z} dN_{ch}/d#it{p}_{T} d#Delta#phi [GeV^{-1}]");
		dphi_plots[i]->SetXTitle("Z-hadron #Delta#phi");
		dphi_plots[i]->SetLineColor(colors[i]);
		dphi_plots[i]->GetYaxis()->SetRangeUser(0,6);
		if(i==0) dphi_plots[i]->Draw("");
		else dphi_plots[i]->Draw("same");
		string label = dphi_plots[i]->GetName();
		myText(.2,.5+.05*i,colors[i],label.substr(zType.size()+10).c_str(),.04);
	}
	tl->Draw();
	string savename = "../plots/";
	savename += zType+"_trackpT_dphi.pdf";
	tc->SaveAs(savename.c_str());
}

void plottest(TFile *thisFile){
	TH1F *test = (TH1F*) thisFile->Get("non-lead 2-3.3 p_{T}");
	TCanvas* tc = new TCanvas();
	tc->Draw();
	test->Draw("e1");
}

void plotdPhiLead(TFile *thisFile,string jetType){
	std::vector<TH1F*> dphi_plots;

	dphi_plots.push_back((TH1F*) thisFile->Get((jetType + " 2-3.3 p_{T}").c_str()));
	dphi_plots.push_back((TH1F*) thisFile->Get((jetType + " 3.3-5.4 p_{T}").c_str()));
	dphi_plots.push_back((TH1F*) thisFile->Get((jetType + " 5.4-8.9 p_{T}").c_str()));
	dphi_plots.push_back((TH1F*) thisFile->Get((jetType + " 8.9-14.6 p_{T}").c_str()));
	dphi_plots.push_back((TH1F*) thisFile->Get((jetType + " 14.6-24.0 p_{T}").c_str()));
	dphi_plots.push_back((TH1F*) thisFile->Get((jetType + " 24.0-39.5 p_{T}").c_str()));
	dphi_plots.push_back((TH1F*) thisFile->Get((jetType + " 39.5-65 p_{T}").c_str()));

  std::vector <string> dphi_plot_labels;
  dphi_plot_labels.push_back ("2 < #it{p}_{T} < 3.3 GeV");
  dphi_plot_labels.push_back ("3.3 < #it{p}_{T} < 5.4 GeV");
  dphi_plot_labels.push_back ("5.4 < #it{p}_{T} < 8.9 GeV");
  dphi_plot_labels.push_back ("8.9 < #it{p}_{T} < 14.6 GeV");
  dphi_plot_labels.push_back ("14.6 < #it{p}_{T} < 24.0 GeV");
  dphi_plot_labels.push_back ("24.0 < #it{p}_{T} < 39.5 GeV");
  dphi_plot_labels.push_back ("39.5 < #it{p}_{T} < 65 GeV");
	
	//short colors[7]={kBlack,kRed,kBlue,kGreen+3,kMagenta-5,kMagenta+4,kCyan-4};
  short colors[10] = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta, kViolet-3, kCyan+1, kOrange+1, kGreen-7, kAzure+7};
	TCanvas* tc = new TCanvas();
	tc->Draw();
	//TLegend* tl = new TLegend(1.5,.7,3.5,.9);
	for (int i = 0; i < dphi_plots.size(); ++i)
	{
		dphi_plots[i]->SetYTitle("1/N_{Z} dN_{ch}/d#it{p}_{T} d#Delta#phi [GeV^{-1}]");
		dphi_plots[i]->SetXTitle("Z-hadron #Delta#phi");
		dphi_plots[i]->SetLineColor(colors[i]);
		dphi_plots[i]->SetMarkerStyle(kFullCircle);
		dphi_plots[i]->SetMarkerColor(colors[i]);
		//dphi_plots[i]->GetYaxis()->SetRangeUser(0,8);
    if (jetType == "leading")
      dphi_plots[i]->GetYaxis ()->SetRangeUser (0, 4);
    else if (jetType == "non-lead")
      dphi_plots[i]->GetYaxis ()->SetRangeUser (0, 1);
		if(i==0) dphi_plots[i]->Draw("");
		else dphi_plots[i]->Draw("same");
		myText(.2,.9-.05*i,colors[i],dphi_plot_labels[i].c_str (),.04);
	}
	//tl->Draw();
  if (jetType == "leading")
    myText (.55, 0.88, kBlack, "Tracks assoc. with lead. jet", 0.04);
  else if (jetType == "non-lead")
    myText (.55, 0.88, kBlack, "Tracks not assoc. with lead. jet", 0.04);
	string savename = "../plots/";
	savename += jetType+"_trackpT_dphi.pdf";
	tc->SaveAs(savename.c_str());
}

void plotJetDPhi(TFile* thisFile){
	std::vector<TH1F*> dphi_plots;
	dphi_plots.push_back((TH1F*) thisFile->Get("10-20 jets"));
	dphi_plots.push_back((TH1F*) thisFile->Get("20-40 jets"));
	dphi_plots.push_back((TH1F*) thisFile->Get("40-80 jets"));
	dphi_plots.push_back((TH1F*) thisFile->Get("80+ jets"));

	short colors[4]={kRed,kBlue,kGreen+3,kMagenta-5};
	TCanvas* tc = new TCanvas();
	TLegend* tl = new TLegend(1.5,.7,3.5,.9);
	for (int i = 0; i < 4; ++i)
	{
		dphi_plots[i]->SetYTitle("1/N_{Z} dN_{jet}/d#Delta#phi");
		dphi_plots[i]->SetXTitle("Z-jet #Delta#phi");
		dphi_plots[i]->GetYaxis()->SetRangeUser(0,2.5);
		dphi_plots[i]->SetLineColor(colors[i]);
		dphi_plots[i]->SetMarkerStyle(kFullCircle);
		dphi_plots[i]->SetMarkerColor(colors[i]);
		if(i==0) dphi_plots[i]->Draw("");
		else dphi_plots[i]->Draw("same");
		string name  = dphi_plots[i]->GetName();
		name+="#it{p}_{T}^{Z} [GeV]";
		tl->AddEntry(dphi_plots[i],name.c_str(),"l");
	}
	tl->Draw();
	tc->SaveAs("../plots/Zjet_njet_dphi.pdf");
}

void plotDPhi(TFile* thisFile,std::vector<string> types){
	std::vector<TH1F*> dphi_pt_plots;
	for (std::vector<string>::iterator i = types.begin(); i != types.end(); ++i)
	{
		dphi_pt_plots.push_back((TH1F*) thisFile->Get((*i+" lowpt").c_str()));
		dphi_pt_plots.push_back((TH1F*) thisFile->Get((*i+" midpt").c_str()));
		dphi_pt_plots.push_back((TH1F*) thisFile->Get((*i+" highpt").c_str()));
	}
	TCanvas* tc = new TCanvas();
	tc->SetLogy();
	tc->SetTicky();

	unsigned count=0;
	unsigned typeCount=0;
	short colors[3]={kBlack,kRed,kBlue};
	short styleTypes[4]={kOpenCircle,kFullTriangleUp,kOpenStar,kFullDiamond};
	TLegend* tl = new TLegend(.7,.1,.9,.3);
	for (std::vector<TH1F*>::iterator i = dphi_pt_plots.begin(); i != dphi_pt_plots.end(); ++i)
	{
		(*i)->GetYaxis()->SetRangeUser(10e-5,10e1);
		(*i)->SetYTitle("1/N_{evt} dN_{ch}/d#Delta#phi");
		(*i)->SetXTitle("#Delta#phi");
		(*i)->SetLineColor(colors[count%3]);
		(*i)->SetMarkerColor(colors[count%3]);
		(*i)->SetMarkerStyle(styleTypes[typeCount]);
		if (count++==0)(*i)->Draw("e1");
		else (*i)->Draw("e1 same");
		if (count%3==0) ++typeCount;
		tl->AddEntry((*i),(*i)->GetName(),"p");
	}
	tl->Draw();
	tc->SaveAs("../plots/dphi.pdf");
}

void plotdiffpT(const std::vector<string>& types, std::vector<TH1F*> pt_plots){
	//can only subtract two data sets
	if (types.size()!=2)
	{
		cout<<"num datasets!=2 difference is not defined"<<endl;
		return;
	}

	//make the new plots
	std::vector<TH1F*> subtracted_plots;
	subtracted_plots.push_back((TH1F*)pt_plots[0]->Clone());
	subtracted_plots.push_back((TH1F*)pt_plots[1]->Clone());
	subtracted_plots.push_back((TH1F*)pt_plots[2]->Clone());

	//setup the canvas
	TCanvas *tc2 = new TCanvas();
	TLegend* tl2 = new TLegend(.65,.1,.85,.25);
	short colors[3]={kBlack,kRed,kBlue};
	short styleTypes[4]={kOpenCircle,kFullTriangleUp,kOpenStar,kFullDiamond};

	unsigned count=0;
	for (std::vector<TH1F*>::iterator i = subtracted_plots.begin(); i != subtracted_plots.end(); ++i)
	{
		(*i)->Add(pt_plots[count+3],-1); //subtract the plots
		(*i)->GetYaxis()->SetRangeUser(-.4,.1); 
		(*i)->SetYTitle("dN / N mpi:off-on");
		(*i)->SetXTitle("#it{p}_{T} [GeV]");
		(*i)->SetLineColor(colors[count%3]); //color the plots
		(*i)->SetMarkerColor(colors[count%3]);
		if (count++==0)(*i)->Draw("e1");
		else (*i)->Draw("e1 same");
		string caption((*i)->GetName());
		caption=caption.substr(caption.find(" "));
		tl2->AddEntry((*i),caption.c_str(),"p");
	}
	tl2->Draw();
	tc2->SaveAs("../plots/pTdiff.pdf");

}

void plotpT(TFile* thisFile,std::vector<string> types){
	std::vector<TH1F*> pt_plots;
	for (std::vector<string>::iterator i = types.begin(); i != types.end(); ++i)
	{
		pt_plots.push_back((TH1F*) thisFile->Get((*i+" near").c_str()));
		pt_plots.push_back((TH1F*) thisFile->Get((*i+" middle").c_str()));
		pt_plots.push_back((TH1F*) thisFile->Get((*i+" away").c_str()));
	}
	TCanvas* tc = new TCanvas();
	tc->SetLogy();
	tc->SetLogx();
	tc->SetTicky();
	TLegend* tl = new TLegend(.7,.7,.9,.95);

	unsigned count=0;
	unsigned typeCount=0;
	short colors[3]={kBlack,kRed,kBlue};
	short styleTypes[4]={kOpenCircle,kFullTriangleUp,kOpenStar,kFullDiamond};
	for (std::vector<TH1F*>::iterator i = pt_plots.begin(); i != pt_plots.end(); ++i)
	{
		(*i)->GetYaxis()->SetRangeUser(10e-11,10e1);
		(*i)->SetYTitle("1/N_{evt} dN_{ch}/d#it{p}_{T} d#Delta#phi [GeV^{-1}]");
		(*i)->SetXTitle("#it{p}_{T} [GeV]");
		(*i)->SetLineColor(colors[count%3]);
		(*i)->SetMarkerColor(colors[count%3]);
		(*i)->SetMarkerStyle(styleTypes[typeCount]);
		if (count++==0)(*i)->Draw("e1");
		else (*i)->Draw("e1 same");
		if (count%3==0) ++typeCount;
		tl->AddEntry((*i),(*i)->GetName(),"p");
	}

	TH1F* minBias = (TH1F*) thisFile->Get("minbias");
	minBias->SetLineColor(kMagenta-7);
	minBias->Draw("same");
	tl->AddEntry(minBias,"minbias","l");
	tl->Draw();
	tc->SaveAs("../plots/pT.pdf");

	plotdiffpT(types,pt_plots);
}

void child_dphi(TFile* thisFile, std::vector<string> types){
	std::vector<TH1F*> plots;
	for (std::vector<string>::iterator i = types.begin(); i != types.end(); ++i)
	{
		plots.push_back((TH1F*) thisFile->Get(((*i)+" child_dphi").c_str()));
		plots.push_back((TH1F*) thisFile->Get(((*i)+" forgein_dphi").c_str()));
	}

	TCanvas* tc = new TCanvas();
	tc->SetLogy();
	tc->SetTicky();
	TLegend* tl = new TLegend(.6,.1,.9,.35);

	unsigned count=0;
	unsigned typeCount=0;
	const int kMAXCOLOR=2;
	short colors[2]={kBlack,kRed};
	short styleTypes[4]={kOpenCircle,kFullTriangleUp,kOpenStar,kFullDiamond};
	for (std::vector<TH1F*>::iterator i = plots.begin(); i != plots.end(); ++i)
	{
		(*i)->GetYaxis()->SetRangeUser(10e-5,10e1);
		(*i)->SetYTitle("1/N_{evt} dN/d#Delta#phi");
		(*i)->SetXTitle("#Delta#phi");
		(*i)->SetLineColor(colors[count%kMAXCOLOR]);
		(*i)->SetMarkerColor(colors[count%kMAXCOLOR]);
		(*i)->SetMarkerStyle(styleTypes[typeCount]);
		if (count++==0)(*i)->Draw("e1");
		else (*i)->Draw("e1 same");
		if (count%kMAXCOLOR==0) ++typeCount;
		tl->AddEntry((*i),(*i)->GetName(),"p");
	}
	tl->Draw();
	tc->SaveAs("../plots/child_dphi.pdf");

}

void child_pT(TFile* thisFile, std::vector<string> types){
	std::vector<TH1F*> plots;
	for (std::vector<string>::iterator i = types.begin(); i != types.end(); ++i)
	{
		plots.push_back((TH1F*) thisFile->Get(((*i)+" child_ntrack").c_str()));
		plots.push_back((TH1F*) thisFile->Get(((*i)+" forgein_ntrack").c_str()));
	}

	TCanvas* tc = new TCanvas();
	tc->SetLogy();
	tc->SetTicky();
	TLegend* tl = new TLegend(.6,.75,.9,.9);

	unsigned count=0;
	unsigned typeCount=0;
	const int kMAXCOLOR=2;
	short colors[2]={kBlack,kRed};
	short styleTypes[4]={kOpenCircle,kFullTriangleUp,kOpenStar,kFullDiamond};
	for (std::vector<TH1F*>::iterator i = plots.begin(); i != plots.end(); ++i)
	{
		(*i)->GetYaxis()->SetRangeUser(10e-5,10e1);
		(*i)->SetYTitle("1/N_{evt} dN/d#Delta#phi");
		(*i)->SetXTitle("#it{p}_{T} [GeV]");
		(*i)->SetLineColor(colors[count%kMAXCOLOR]);
		(*i)->SetMarkerColor(colors[count%kMAXCOLOR]);
		(*i)->SetMarkerStyle(styleTypes[typeCount]);
		if (count++==0)(*i)->Draw("e1");
		else (*i)->Draw("e1 same");
		if (count%kMAXCOLOR==0) ++typeCount;
		tl->AddEntry((*i),(*i)->GetName(),"p");
		cout<<(*i)->GetName()<<": "<<(*i)->Integral()<<'\n';
	}
	tl->Draw();
	tc->SaveAs("../plots/child_ntrack.pdf");

}

void plotter(){
	//std::vector<string> types;
	//types.push_back("inclusive_mpioff");
	//types.push_back("inclusive_mpion");
	//types.push_back("forced_mpioff");
	//types.push_back("forced_mpion");
	TFile *thisFile = new TFile("hists.root","READ");
	//plotDPhi(thisFile,types);
	//plotpT(thisFile,types);
	//child_dphi(thisFile,types);
	//child_pT(thisFile,types);
	//plotJetnPart(thisFile);
	//plotJetnPart2(thisFile);
	//plotJetDPhiTracks(thisFile);
	//plotJetDPhi(thisFile);
	//plotJetDPhiTracksWide(thisFile,"15-25");
	//plotJetDPhiTracksWide(thisFile,"25+");
	plotdPhiLead(thisFile,"leading");
	plotdPhiLead(thisFile,"non-lead");
	//plottest(thisFile);
}
