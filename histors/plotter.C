

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

	unsigned count=0;
	short colors[3]={kBlack,kRed,kBlue};
	short styleTypes[3]={kFullCircle,kFullTriangleUp,kFullStar};
	unsigned typeCount=0;
	TLegend* tl = new TLegend(.8,.1,.9,.25);
	for (std::vector<TH1F*>::iterator i = dphi_pt_plots.begin(); i != dphi_pt_plots.end(); ++i)
	{
		(*i)->GetYaxis()->SetRangeUser(10e-5,10e1);
		(*i)->SetLineColor(colors[count%3]);
		(*i)->SetMarkerColor(colors[count%3]);
		(*i)->SetMarkerStyle(styleTypes[typeCount]);
		if (count%3==0) ++typeCount;
		if (count++==0)(*i)->Draw("e1");
		else (*i)->Draw("e1 same");
		tl->AddEntry((*i),(*i)->GetName(),"p");
	}
	tl->Draw();
	tc->SaveAs("../plots/dphi.pdf");
}

/*void plotpT(TFile* thisFile){
	std::vector<TH1F*> pt_plots;
	pt_plots.push_back((TH1F*) thisFile->Get((status+" near").c_str()));
	pt_plots.push_back((TH1F*) thisFile->Get((status+" middle").c_str()));
	pt_plots.push_back((TH1F*) thisFile->Get((status+" away").c_str()));

	TCanvas* tc = new TCanvas();
	tc->SetLogy();
	tc->SetLogx();

	unsigned count=0;
	short colors[3]={kBlack,kRed,kBlue};

	TLegend* tl = new TLegend(.8,.1,.9,.25);
	for (std::vector<TH1F*>::iterator i = dphi_pt_plots.begin(); i != dphi_pt_plots.end(); ++i)
	{
		(*i)->GetYaxis()->SetRangeUser(10e-8,10e1);
		(*i)->SetLineColor(colors[count]);
		if (count++==0)(*i)->Draw("e1");
		else (*i)->Draw("e1 same");
		tl->AddEntry((*i),(*i)->GetName(),"l");
	}
	tl->Draw();
	tc->SaveAs(("../plots/"+status+"_pt.pdf").c_str());
}*/

void plotter(){
	std::vector<string> types;
	types.push_back("fff");
	types.push_back("ffo");
	TFile *thisFile = new TFile("hists.root","READ");
	plotDPhi(thisFile,types);
}
