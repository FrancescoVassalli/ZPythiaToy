#include "TMath.h"
#include "TVector3.h"
#include "TH1F.h"


#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <sstream>

#include "Pythia8/Pythia.h"
#include "Pythia8/FJcore.h"
#include "Pythia8/Event.h"
#include "Pythia8/Basics.h"

#include <ctime>

//#include <Utilities.h>

using namespace Pythia8;
using namespace std;

Double_t E= 2.71828182845904523536;


TLorentzVector* pToTLV(Vec4 in){
	Double_t px = (double)in.px();
	Double_t py = (double)in.py();
	Double_t pz = (double)in.pz();
	Double_t e =  (double)in.e();
	TLorentzVector *out = new TLorentzVector(px,py,pz,e);
	return out;
}

int main (int argc, char *argv[]) {

	// Generator. Process selection. LHC initialization. Histogram.
	Pythia pythia;
	pythia.readString ("Beams:eCM = 5020.");
	pythia.readString("Random::setSeed = on");
	pythia.readString("Random::seed =0");
	pythia.readString("SoftQCD:nonDiffractive = on");
	pythia.readString("SoftQCD:singleDiffractive = on");
	pythia.readString("SoftQCD:doubleDiffractive = on");
	//2GeV is suggested
	ostringstream ss; ss << "PhaseSpace:pTHatMin = " << argv[1] << ".";
	pythia.readString (ss.str ().c_str ());
	pythia.init ();

	int NEVT = atoi (argv[2]);

	TFile *f = new TFile (argv[3] , "RECREATE");

	vector<float> b_part_pt, b_part_eta, b_part_phi;

	TTree *t = new TTree("tree","an isotropic pineapple tree");

	t->Branch ("part_pt", &b_part_pt);
	t->Branch ("part_eta", &b_part_eta);
	t->Branch ("part_phi", &b_part_phi);

	double pythiaTimeInSeconds=0;

	for (int iEvent = 0; iEvent < NEVT; iEvent++) {
		clock_t startPythia = clock();
		if (!pythia.next ())
			continue;
		clock_t endPythia = clock();
		pythiaTimeInSeconds+=(endPythia-startPythia) / (double) CLOCKS_PER_SEC;

		b_part_pt.clear ();
		b_part_eta.clear ();
		b_part_phi.clear ();

		for (int i = 0; i < pythia.event.size (); i++) {

			//record track info
			if (pythia.event[i].pT() >= 2 && pythia.event[i].isCharged() && pythia.event[i].isHadron()) {
				b_part_pt.push_back (pythia.event[i].pT ());
				b_part_eta.push_back (pythia.event[i].eta ());
				b_part_phi.push_back (pythia.event[i].phi ());
			}
		}
		t->Fill();
		if (NEVT>100&&iEvent % (NEVT/100) == 0)
			std::cout << iEvent / (NEVT/100) << "\% done...\r" << std::flush;
	}

	pythia.stat();

	f->Write();
	f->Close();

	cout<<"Done with pythiatime="<<pythiaTimeInSeconds<<std::endl;
	return 0;
}
