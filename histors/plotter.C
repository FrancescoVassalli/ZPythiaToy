

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
	short styleTypes[4]={kOpenCircle,kFullTriangleUp,kOpenStar,kFullDiamond};
	unsigned typeCount=0;
	TLegend* tl = new TLegend(.7,.1,.9,.3);
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

	std::vector<TH1F*> subtracted_plots;
	subtracted_plots.push_back((TH1F*)pt_plots[0]->Clone());
	subtracted_plots.push_back((TH1F*)pt_plots[1]->Clone());
	subtracted_plots.push_back((TH1F*)pt_plots[2]->Clone());


	unsigned count=0;
	short colors[3]={kBlack,kRed,kBlue};
	short styleTypes[4]={kOpenCircle,kFullTriangleUp,kOpenStar,kFullDiamond};
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

	TCanvas *tc2 = new TCanvas();
	count=0;
	TLegend* tl2 = new TLegend(.7,.7,.9,.95);
	for (std::vector<TH1F*>::iterator i = subtracted_plots.begin(); i != subtracted_plots.end(); ++i)
	{
		(*i)->Add(pt_plots[count+3],-1);
		(*i)->SetLineColor(colors[count%3]);
		if (count++==0)(*i)->Draw("e1");
		else (*i)->Draw("e1 same");
		tl2->AddEntry((*i),(*i)->GetName(),"p");
	}
	tl2->Draw();
}

void plotter(){
	std::vector<string> types;
	//types.push_back("inclusive_mpioff");
	//types.push_back("inclusive_mpion");
	types.push_back("forced_mpioff");
	types.push_back("forced_mpion");
	TFile *thisFile = new TFile("hists.root","READ");
	plotDPhi(thisFile,types);
	plotpT(thisFile,types);
}
