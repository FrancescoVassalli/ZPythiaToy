

void plotter(){
	TFile *thisFile = new TFile("hists.root","READ");


	std::vector<TH1F*> dphi_pt_plots;
	dphi_pt_plots.push_back((TH1F*) thisFile->Get("lowpt"));
	dphi_pt_plots.push_back((TH1F*) thisFile->Get("midpt"));
	dphi_pt_plots.push_back((TH1F*) thisFile->Get("highpt"));

	TCanvas* tc = new TCanvas();
	tc->SetLogy();

	unsigned count=0;
	short colors[3]={kBlack,kRed,kBlue};

	TLegend* tl = new TLegend(.8,.7,.9,.95);
	
	for (std::vector<TH1F*>::iterator i = dphi_pt_plots.begin(); i != dphi_pt_plots.end(); ++i)
	{
		(*i)->GetYaxis()->SetRangeUser(10e-8,10e1);
		(*i)->SetLineColor(colors[count]);
		if (count++==0)(*i)->Draw("e1");
		else (*i)->Draw("e1 same");
		tl->AddEntry((*i),(*i)->GetName(),"l");
	}
	tl->Draw();

}
