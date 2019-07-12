

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
	short styleTypes[4]={kFullCircle,kFullTriangleUp,kFullStar,kOpenDiamond};
	unsigned typeCount=0;
	TLegend* tl = new TLegend(.8,.1,.95,.3);
	for (std::vector<TH1F*>::iterator i = dphi_pt_plots.begin(); i != dphi_pt_plots.end(); ++i)
	{
		(*i)->GetYaxis()->SetRangeUser(10e-5,10e1);
		(*i)->SetYTitle("#frac{dN}{N}");
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

	unsigned count=0;
	short colors[3]={kBlack,kRed,kBlue};
	short styleTypes[4]={kFullCircle,kFullTriangleUp,kFullStar,kOpenDiamond};
	unsigned typeCount=0;
	TLegend* tl = new TLegend(.7,.7,.9,.95);
	for (std::vector<TH1F*>::iterator i = pt_plots.begin(); i != pt_plots.end(); ++i)
	{
		(*i)->GetYaxis()->SetRangeUser(10e-7,10e1);
		(*i)->SetYTitle("#frac{dN}{N}");
		(*i)->SetXTitle("pT");
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
}

void plotter(){
	std::vector<string> types;
	types.push_back("oof");
	types.push_back("ooo");
	types.push_back("in");
	TFile *thisFile = new TFile("hists.root","READ");
	plotDPhi(thisFile,types);
	plotpT(thisFile,types);
}
