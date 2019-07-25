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

double InTwoPi (double phi) {
		while (phi < 0 || 2*TMath::Pi() <= phi) {
			if (phi < 0) phi += 2*TMath::Pi();
			else phi -= 2*TMath::Pi();
		}
		return phi;
	}

	double DeltaPhi (double phi1, double phi2, const bool sign=0) {
		phi1 = InTwoPi(phi1);
		phi2 = InTwoPi(phi2);
		double dphi = abs(phi1 - phi2);
		while (dphi > TMath::Pi()) dphi = abs (dphi - 2*TMath::Pi());

		if (sign && InTwoPi (phi2 + dphi) == phi1)
			dphi *= -1;

		return dphi;
	}

/**
 * Returns dR between two eta, phi coordinates.
 */
double DeltaR (const double eta1, const double eta2, const double phi1, const double phi2 ) {
 const double deta = eta1 - eta2;
 const double dphi = DeltaPhi (phi1, phi2, false);
 return sqrt( pow( deta, 2 ) + pow( dphi, 2 ) );
}


void tagChildren(Pythia8::Particle parent,std::set<int>* childIndexSet,Pythia* pythia){
	childIndexSet->insert(parent.daughter1());
	childIndexSet->insert(parent.daughter2());
	if (parent.daughter1()>0)
	{
		tagChildren(pythia->event[parent.daughter1()],childIndexSet,pythia);
	}
	if (parent.daughter2()>0)
	{
		tagChildren(pythia->event[parent.daughter2()],childIndexSet,pythia);
	}
}

int main (int argc, char *argv[]) {

	if (argc < 4) {
		cout<<"args error"<<endl;
	}

	// Generator. Process selection. LHC initialization. Histogram.
	Pythia pythia;

	pythia.readString ("Beams:eCM = 5020.");
	pythia.readString("Random::setSeed = on");
	pythia.readString("Random::seed =0");
	pythia.readString("PartonLevel:FSR=off");
	pythia.readString("PartonLevel:ISR=on");
	pythia.readString("PartonLevel:MPI=on");
	if (argc>=5)
	{
		if (argv[4][0]=='i')
		{
			pythia.readString ("WeakSingleBoson:ffbar2gmZ = on");       // code 221
			pythia.readString ("WeakDoubleBoson:ffbar2gmZgmZ = on");    // code 231
			pythia.readString ("WeakDoubleBoson:ffbar2ZW = on");        // code 232
			pythia.readString ("WeakBosonAndParton:ffbar2gmZgm = on");  // code 243
			pythia.readString ("WeakBosonAndParton:fgm2gmZf = on");     // code 244
		}
		if (argv[4][1]=='o')
		{
			pythia.readString("PartonLevel:MPI=on");
		}
	}

	pythia.readString("23:onMode = off");
	pythia.readString("23:onIfAny = 11 13");

	pythia.readString ("WeakZ0:gmZmode = 2"); // set to Z's

	pythia.readString ("WeakBosonAndParton:qqbar2gmZg = on");   // code 241
	pythia.readString ("WeakBosonAndParton:qg2gmZq = on");      // code 242

	ostringstream ss; ss << "PhaseSpace:pTHatMin = " << argv[1] << ".";
	pythia.readString (ss.str ().c_str ());
	pythia.readString ("PhaseSpace:mHatMin = 60");

	pythia.init ();

	SlowJet *antikT04 = new SlowJet (-1, 0.4, 1, 5, 2, 1);
	//  SlowJet *antikT10 = new SlowJet (-1, 1.0, 1, 5, 2, 1);

	int NEVT = atoi (argv[2]);

	TFile *f = new TFile (argv[3] , "RECREATE");

	int b_code;
	int b_id1;
	int b_id2;
	float b_x1pdf;
	float b_x2pdf;
	float b_Q;
	bool b_isValence1;
	bool b_isValence2;

	int b_z_n, b_part_n, b_lead_n, b_sublead_n;
	vector<float> b_z_pt, b_z_eta, b_z_phi, b_z_m;
	vector<float> b_l_pt, b_l_eta, b_l_phi, b_l_m;
	vector<float> b_part_pt, b_part_eta, b_part_phi;
	std::vector<bool> b_part_child;
	float b_lead_pt, b_lead_eta, b_lead_phi;
	float b_sublead_pt, b_sublead_eta, b_sublead_phi;
	vector<float> b_jet_r04_pt, b_jet_r04_eta, b_jet_r04_phi, b_jet_r04_e;
	//vector<float> b_jet_r10_pt, b_jet_r10_eta, b_jet_r10_phi, b_jet_r10_e;

	TTree *t = new TTree("tree","a volumptous pineapple tree");

	t->Branch ("code",  &b_code);
	t->Branch ("id1",   &b_id1);
	t->Branch ("id2",   &b_id2);
	t->Branch ("x1pdf", &b_x1pdf);
	t->Branch ("x2pdf", &b_x2pdf);
	t->Branch ("Q",     &b_Q);

	t->Branch ("z_n",   &b_z_n);
	t->Branch ("z_pt",  &b_z_pt);
	t->Branch ("z_eta", &b_z_eta);
	t->Branch ("z_phi", &b_z_phi);
	t->Branch ("z_m",   &b_z_m);

	t->Branch ("part_n", &b_part_n);
	t->Branch ("part_pt", &b_part_pt);
	t->Branch ("part_eta", &b_part_eta);
	t->Branch ("part_phi", &b_part_phi);
	t->Branch ("part_child", &b_part_child);

	t->Branch ("lead_n",   &b_lead_n);
	t->Branch ("lead_pt",  &b_lead_pt);
	t->Branch ("lead_eta", &b_lead_eta);
	t->Branch ("lead_phi", &b_lead_phi);

	t->Branch ("sublead_n",   &b_sublead_n);
	t->Branch ("sublead_pt",  &b_sublead_pt);
	t->Branch ("sublead_eta", &b_sublead_eta);
	t->Branch ("sublead_phi", &b_sublead_phi);
	/* t->Branch ("l_pt",  &b_l_pt);
		 t->Branch ("l_eta", &b_l_eta);
		 t->Branch ("l_phi", &b_l_phi);
		 t->Branch ("l_m",   &b_l_m);*/

	//t->Branch ("jet_r04_n",   &b_jet_r04_n);
	t->Branch ("jet_r04_pt",  &b_jet_r04_pt);
	t->Branch ("jet_r04_eta", &b_jet_r04_eta);
	t->Branch ("jet_r04_phi", &b_jet_r04_phi);
	t->Branch ("jet_r04_e",   &b_jet_r04_e);

		 /*t->Branch ("jet_r10_n",   &b_jet_r10_n);
		 t->Branch ("jet_r10_pt",  &b_jet_r10_pt);
		 t->Branch ("jet_r10_eta", &b_jet_r10_eta);
		 t->Branch ("jet_r10_phi", &b_jet_r10_phi);
		 t->Branch ("jet_r10_e",   &b_jet_r10_e);*/

	double pythiaTimeInSeconds=0;
	unsigned nFinalChildren=0;
	unsigned totalGluon=0;

	for (int iEvent = 0; iEvent < NEVT; iEvent++) {
		clock_t startPythia = clock();
		if (!pythia.next ())
			continue;
		clock_t endPythia = clock();
		pythiaTimeInSeconds+=(endPythia-startPythia) / (double) CLOCKS_PER_SEC;

		b_z_n = 0;
		b_z_pt.clear ();
		b_z_eta.clear ();
		b_z_phi.clear ();
		b_z_m.clear ();

		b_part_n = 0;
		b_part_pt.clear ();
		b_part_eta.clear ();
		b_part_phi.clear ();
		b_part_child.clear();

		b_l_phi.clear();
		b_l_eta.clear();
		b_l_pt.clear();
		b_l_m.clear();


		/*get the lepton info using event record 
			set<int> ZChildIndicies;
			tagChildren(pythia.event[5],&ZChildIndicies, &pythia);
			set<int>::reverse_iterator rit=ZChildIndicies.rbegin();
			b_l_phi.push_back(pythia.event[*rit].phi());
			b_l_eta.push_back(pythia.event[*rit].eta());
			b_l_pt.push_back(pythia.event[*rit].pT());
			b_l_m.push_back(pythia.event[*rit].m());
			++rit;
			b_l_phi.push_back(pythia.event[*rit].phi());
			b_l_eta.push_back(pythia.event[*rit].eta());
			b_l_pt.push_back(pythia.event[*rit].pT());
			b_l_m.push_back(pythia.event[*rit].m());*/

		set<int> partonChildIndicies;
		tagChildren(pythia.event[6],&partonChildIndicies, &pythia);
		/*for (std::set<int>::iterator i = partonChildIndicies.begin(); i != partonChildIndicies.end(); ++i)
			{
			if(pythia.event[(*i)].isFinal()&&pythia.event[(*i)].pT()>2&&pythia.event[(*i)].isHadron()&& TMath::Abs(pythia.event[(*i)].eta())<2.5) nFinalChildren++;
			}*/

		//cout<<"children of "<<pythia.event[6].id()<<"\n";
		/*for (std::set<int>::iterator it = ZChildIndicies.begin(); it != ZChildIndicies.end(); ++it)
			{
			cout<<(*it)<<'\n';
			}*/
		unsigned tempChildCount = 0;
		TLorentzVector *l1=nullptr,*l2=nullptr;
		for (int i = 0; i < pythia.event.size (); i++) {
			if (pythia.event[i].pT() >= 2 && pythia.event[i].isFinal() && TMath::Abs(pythia.event[i].eta())<2.5) {
				//reconstruct z
				if (pythia.event[i].isLepton())
				{
					if (!l1)
					{
						l1 = new TLorentzVector();
						l1->SetPtEtaPhiM(pythia.event[i].pT (), pythia.event[i].eta (), pythia.event[i].phi (), pythia.event[i].m ());
						for (int j = i+1; j <pythia.event.size () ; j++) {

							if (pythia.event[j].isFinal()&&pythia.event[i].id () == -pythia.event[j].id ()) {
								l2 = new TLorentzVector();
								l2->SetPtEtaPhiM(pythia.event[j].pT (), pythia.event[j].eta (), pythia.event[j].phi (), pythia.event[j].m ());
								if ((*l1+*l2).M () >40&& (*l1+*l2).Pt ()>10)
								{
									b_z_pt.push_back ((*l1+*l2).Pt ());
									b_z_eta.push_back ((*l1+*l2).Eta ());
									b_z_phi.push_back ((*l1+*l2).Phi ());
									b_z_m.push_back ((*l1+*l2).M ());
									b_z_n++;
									b_l_pt.push_back(l1->Pt());
									b_l_pt.push_back(l2->Pt());
									b_l_eta.push_back(l1->Eta());
									b_l_eta.push_back(l2->Eta());
									b_l_phi.push_back(l1->Phi());
									b_l_phi.push_back(l2->Phi());
								}
								delete l2;
								l2 = nullptr;
							}
						}
						delete l1;
						l1=nullptr;
					}
				}
				//record tracks
				else 
				{
					b_part_pt.push_back (pythia.event[i].pT ());
					b_part_eta.push_back (pythia.event[i].eta ());
					b_part_phi.push_back (pythia.event[i].phi ());
					b_part_child.push_back(partonChildIndicies.count(i));
					b_part_n++;
				}
				
				if (partonChildIndicies.count(i))
				{
					tempChildCount++;
				}
			}
			/*record the info for the final Z using event record
			if (pythia.event[i].pT()>=25&& abs(pythia.event[i].id ()) == 23 
					&& (abs(pythia.event[pythia.event[i].daughter1()].id())==11 
						|| abs(pythia.event[pythia.event[i].daughter1()].id())==13)) {
				if (b_z_pt.size()>=1)
				{
					cout<<i<<":\n \t"<<pythia.event[i].pT ()<<"\n \t"<<pythia.event[i].m ()<<'\n';
				}
				b_z_pt.push_back (pythia.event[i].pT ());
				b_z_eta.push_back (pythia.event[i].eta ());
				b_z_phi.push_back (pythia.event[i].phi ());
				b_z_m.push_back (pythia.event[i].m ());
				b_z_n++;
			}*/
		}
		if (b_z_n !=1) {//check there is Z in event
			iEvent--;
			continue;
		}
		nFinalChildren+=tempChildCount;
		if (pythia.event[6].id()==21)
		{
			totalGluon++;
		}

		antikT04->analyze (pythia.event);
		//antikT10->analyze (pythia.event);
		b_lead_pt = -1.;
		b_sublead_pt=-1.;
		b_jet_r04_eta.clear();
		b_jet_r04_pt.clear();
		b_jet_r04_phi.clear();
		for (int i = 0; i < antikT04->sizeJet (); i++) {
			bool matchesLepton = false;
			for (int j = 0; !matchesLepton && j < b_l_eta.size(); j++) {
				matchesLepton = DeltaR (antikT04->p (i).eta (), b_l_eta.at (j), antikT04->phi (i), b_l_phi.at (j)) < 0.2;
			}
			if (matchesLepton) continue;

			b_jet_r04_pt.push_back( antikT04->pT (i));
			b_jet_r04_eta.push_back( antikT04->p (i).eta ());
			b_jet_r04_phi.push_back( antikT04->phi (i));

			if (b_lead_pt<0)
			{
				b_lead_pt = antikT04->pT (i);
				b_lead_eta = antikT04->p (i).eta ();
				b_lead_phi = antikT04->phi (i);
				b_lead_n = antikT04->constituents(i).size();
			}
			else if(antikT04->pT (i)>b_sublead_pt){
				b_sublead_pt = antikT04->pT (i);
				b_sublead_eta = antikT04->p (i).eta ();
				b_sublead_phi = antikT04->phi (i);
				b_sublead_n = antikT04->constituents(i).size();
			}
		}
		if (b_sublead_pt>b_lead_pt)
		{
			float pTtemp, etaTemp, phiTemp;
			int nTemp;
			pTtemp= b_sublead_pt;
			etaTemp = b_sublead_eta;
			phiTemp = b_sublead_phi;
			nTemp = b_sublead_n;
			b_sublead_n = b_lead_n;
			b_sublead_pt = b_lead_pt;
			b_sublead_phi = b_lead_phi;
			b_sublead_eta = b_lead_eta;
			b_lead_eta = etaTemp;
			b_lead_phi = phiTemp;
			b_lead_pt = pTtemp;
			b_lead_n = nTemp;
		}

		b_code = pythia.info.code ();
		b_id1 = pythia.info.id1pdf ();
		b_id2 = pythia.info.id2pdf ();
		b_x1pdf = pythia.info.x1pdf ();
		b_x2pdf = pythia.info.x2pdf ();
		b_Q =  pythia.info.QFac ();

		b_isValence1 = pythia.info.isValence1 ();
		b_isValence2 = pythia.info.isValence2 ();
		t->Fill();

		/*
			for (int i = 0; i < pythia.event.size (); i++) {

			if (!pythia.event[i].isFinal()) continue; // check if in final state

			if (abs (pythia.event[i].id ()) != 11 && abs (pythia.event[i].id ()) != 13) continue; // check if electron or muon, resp.

			b_l_pt.push_back (pythia.event[i].pT ());
			b_l_eta.push_back (pythia.event[i].eta ());
			b_l_phi.push_back (pythia.event[i].phi ());
			b_l_m.push_back (pythia.event[i].m ());
			b_l_n++;
			}

			b_jet_r04_n = 0;
			b_jet_r04_pt.clear ();
			b_jet_r04_eta.clear ();
			b_jet_r04_phi.clear ();
			b_jet_r04_e.clear ();

			for (int i = 0; i < antikT04->sizeJet (); i++) {

			bool matchesLepton = false;
			for (int j = 0; !matchesLepton && j < b_l_n; j++) {
			matchesLepton = DeltaR (antikT04->p (i).eta (), b_l_eta.at (j), antikT04->phi (i), b_l_phi.at (j)) < 0.2;
			}
			if (matchesLepton) continue;

			b_jet_r04_pt.push_back (antikT04->pT (i));
			b_jet_r04_eta.push_back (antikT04->p (i).eta ());
			b_jet_r04_phi.push_back (antikT04->phi (i));
			b_jet_r04_e.push_back (antikT04->p (i).e ());
			b_jet_r04_n++;
			}

			b_jet_r10_n = 0;
			b_jet_r10_pt.clear ();
			b_jet_r10_eta.clear ();
			b_jet_r10_phi.clear ();
			b_jet_r10_e.clear ();

			for (int i = 0; i < antikT10->sizeJet (); i++) {

			bool matchesLepton = false;
			for (int j = 0; !matchesLepton && j < b_l_n; j++) {
			matchesLepton = DeltaR (antikT10->p (i).eta (), b_l_eta.at (j), antikT10->phi (i), b_l_phi.at (j)) < 0.2;
			}
			if (matchesLepton) continue;

			b_jet_r10_pt.push_back (antikT10->pT (i));
			b_jet_r10_eta.push_back (antikT10->p (i).eta ());
			b_jet_r10_phi.push_back (antikT10->phi (i));
			b_jet_r10_e.push_back (antikT10->p (i).e ());
			b_jet_r10_n++;
			}*/


		if (NEVT>100&&iEvent % (NEVT/100) == 0)
			std::cout << iEvent / (NEVT/100) << "\% done...\r" << std::flush;
	}

	pythia.stat();

	f->Write();
	f->Close();

	cout<<"Children ="<<(double)nFinalChildren<<'\n';
	cout<<"gluon fraction="<<(double)totalGluon/NEVT<<'\n';
	cout<<"Done with pythiatime="<<pythiaTimeInSeconds<<std::endl;
	return 0;
}
