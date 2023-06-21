#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;
#include "../pbmodel/F1F209Wrapper.hh"

void hms_foil_rates_carbon_inelastic(TString basename="temp",double cur=70){
   if (basename=="temp") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="worksim/"+basename+".root";
   TString outputhist;
   outputhist="worksim/hists/"+basename+"_hist.root";
 TObjArray HList(0);
     TString outputpdf;
     //   outputpdf="plots/"+basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");

  //  gSystem->Load("pbmodel/libF1F209.so");
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("h1");
Float_t         hsxfp; // position at focal plane ,+X is pointing down
 Float_t         hsyfp; // X x Y = Z so +Y pointing central ray left
 Float_t         hsxpfp; // dx/dz at focal plane
 Float_t         hsypfp; //  dy/dz at focal plane
 Float_t         hsztari; // thrown position along the beam direction
 Float_t         hsytari;  //thrown  horizontal position X x Y = Z so +Y pointing central ray left at plane perpendicular to SHMS at z=0
 Float_t         hsdeltai; // thrown  100*(p - pc)/pc with pc = central SHMS momentum 
 Float_t         hsyptari; // thrown target dy/dz horizontal slope
 Float_t         hsxptari; // thrown target dx/dz vertical slope
 Float_t         hsztar; // reconstructed position along the beam directio
 Float_t         hsytar; //reconstructed horizontal position
 Float_t         hsdelta;//reconstructed
   Float_t         hsyptar;//reconstructed
   Float_t         hsxptar;//reconstructed
   Float_t         hsxtari;// from thrown kinematics , the calculated vertical position at plane perpendicular to SHMS at z=0
   Float_t wfac,xn;
   Float_t ys,xs;
   Float_t event_type;
   Float_t sig;
   Float_t beam_e;
   Float_t p_spec_b;
   Float_t th_spec;
   // Set branch addresses.
   //   tsimc->SetBranchAddress("xsieve",&xs);
   //  tsimc->SetBranchAddress("ysieve",&ys);
   tsimc->SetBranchAddress("hsxfp",&hsxfp);
   tsimc->SetBranchAddress("hsyfp",&hsyfp);
   tsimc->SetBranchAddress("hsxpfp",&hsxpfp);
   tsimc->SetBranchAddress("hsypfp",&hsypfp);
   tsimc->SetBranchAddress("hsytari",&hsytari);
   tsimc->SetBranchAddress("ztari",&hsztari);
   tsimc->SetBranchAddress("hsdeltai",&hsdeltai);
   tsimc->SetBranchAddress("hsyptari",&hsyptari);
   tsimc->SetBranchAddress("hsxptari",&hsxptari);
   tsimc->SetBranchAddress("hsytar",&hsytar);
   tsimc->SetBranchAddress("hsdelta",&hsdelta);
   tsimc->SetBranchAddress("hsyptar",&hsyptar);
   tsimc->SetBranchAddress("hsxptar",&hsxptar);
   tsimc->SetBranchAddress("wfac",&wfac);// wfac is domega*denergy/n_thrown rad*MeV
   tsimc->SetBranchAddress("beam_e",&beam_e);
   tsimc->SetBranchAddress("p_spec",&p_spec_b);
   tsimc->SetBranchAddress("th_spec",&th_spec);

//define simulation histograms
TH1F *hytar = new TH1F("hytar", "hytar", 100, -5., 5.);
TH1F *hWw = new TH1F("hWw", "hWw", 100, 0, 5.0);
hWw->GetXaxis()->SetTitle("W [GeV]");
hWw->GetYaxis()->SetTitle("counts/sec");
TH1F *hWQ2 = new TH1F("hWQ2", "hWQ2", 100, 0, 7);
hWQ2->GetXaxis()->SetTitle("Q^2 [GeV^2]");
hWQ2->GetYaxis()->SetTitle("counts/sec");
TH1F *hWom = new TH1F("hWom", "hWom",100, 0, 5.0);
hWom->GetXaxis()->SetTitle("Ei - Ef [GeV]");
hWom->GetYaxis()->SetTitle("counts/sec");
TH1F *hth = new TH1F("hth", "hth",160, 0, 40);
hth->GetXaxis()->SetTitle("theta [deg]");
hth->GetYaxis()->SetTitle("counts/sec");
     Double_t x_r=40.;
   Double_t y_r= 40.;
 Double_t xp_r=.1;
 Double_t yp_r=.04;
    TH2F* h_xs_ys;
   h_xs_ys = new TH2F("h_xs_ys",";x_sieve;y_sieve",100,-15,15,100,-15,15);
   TH2F* h_xfp_yfp;
   h_xfp_yfp = new TH2F("h_xfp_yfp",";x_fp;y_fp",100,-x_r,x_r,100,-y_r,y_r);
   TH2F* h_xfp_ypfp;
   h_xfp_ypfp = new TH2F("h_xfp_ypfp",";x_fp;y_pfp",100,-x_r,x_r,100,-yp_r,yp_r);
   TH2F* h_xfp_xpfp;
   h_xfp_xpfp = new TH2F("h_xfp_xpfp",";x_fp;xp_fp",100,-x_r,x_r,100,-xp_r,xp_r);
   TH2F* h_yfp_ypfp;
   h_yfp_ypfp = new TH2F("h_yfp_ypfp",";y_fp;yp_fp",100,-y_r,y_r,100,-yp_r,yp_r);
   TH2F* h_yfp_xpfp;
   h_yfp_xpfp = new TH2F("h_yfp_xpfp",";y_fp;xp_fp",100,-y_r,y_r,100,-xp_r,xp_r);
   TH2F* h_xpfp_ypfp;
   h_xpfp_ypfp = new TH2F("h_xpfp_ypfp",";xp_fp;yp_fp",100,-xp_r,xp_r,100,-yp_r,yp_r);
   TH2F* h_xpfp_delta;
   h_xpfp_delta = new TH2F("h_xpfp_delta","Inelastic ;delta;xp_fp",100,-10,10,100,-0.1,+0.1);

   TH2F* h_ytar_delta;
   h_ytar_delta = new TH2F("h_ytar_delta","Inelastic ;delta;ytar",100,-10,10,100,-10.,+10.);
   TH2F* h_ytar_yptar;
   h_ytar_yptar = new TH2F("h_ytar_yptar",";yptar;ytar",100,-.05,.05,100,-10.,+10.);
   TH2F* h_ztar_yptar_all;
   h_ztar_yptar_all = new TH2F("h_ztar_yptar_all",";ztar:yptar",100,-25.,+25.,100,-.05,.05);
   TH1F* h_ytar;
   h_ytar =  new TH1F("h_ytar"," ;ytar",100,-5.,5.);


Double_t W;
Double_t W2;
Double_t Q2;
Double_t nu;
Double_t thetaDeg;





Double_t ts;
Double_t p_spec;
double Ei;//Beam energy //GeV
//Get the beam energy, spectrometer angle and momentum from the first entry of the ntuple.
tsimc->GetEntry(0);
Ei = beam_e/1000.0;
p_spec = p_spec_b/1000.0;
ts = th_spec;

Long64_t nentries = tsimc->GetEntries(); 

//double ts = 7.5; //central spectrometer angle //deg
/*
double p_spec = 3.0; //spectrometer momentum //GeV 
double Ei = 6.418; //Beam energy //GeV
double ts = 25.; //central spectrometer angle //deg
*/
double Z = 6.0;
double A = 12.0;
double Mp = 0.93825; //proton's mass //GeV
 double Mc = 12.0107*931.5/1000.;
double N_A = 6.02*1e+23; // Avogadro's number
double Q_E = 1.60*1e-19; // Conversion from charge to electron
double deg2rad = 3.14159/180.;
double cos_ts = cos(ts * deg2rad);
double sin_ts = sin(ts * deg2rad);
double car_density = 2.2; // density of carbon g/cm3
double mass_tar = 12. * 931.5; //mass of the target
//double cur; //current uA
double thick; // target thickness g/cm2
double run_time =515.;// secondes;
double lumin ;// luminosity per ub 
double Ef;
double theta; 
double sig_inelastic;
 double ex_calc;
 double sig_elas_calc;
F1F209Wrapper pF1F209;
 Float_t cfac ;
    Float_t weight;
//
  cfac=1.;
  thick= 0.044; // foil thickness g/cm2 in multifoil
  thick= 0.1749; // 0.5% single carbon
lumin= thick*cur/A*N_A/Q_E*1e-36;// lumin 1/ub for cur uA

Double_t rate;
TText *t;

	for (int i = 0; i < nentries; i++) {
    tsimc->GetEntry(i);
    
    // Define kinematics
    Ef = p_spec * (1.0 + 0.01*hsdelta); //scattered electron energy //GeV
    nu = Ei - Ef; //GeV
    W2=0;
    h_ytar->Fill(hsytari);
    h_ztar_yptar_all->Fill(hsztari,hsyptari);      
    if (nu >0) {
      theta = TMath::ACos((cos_ts - hsyptar * sin_ts) / TMath::Sqrt( 1. + hsxptar * hsxptar + hsyptar * hsyptar )); // polar 			scattering angle relative to the beam line //rad
      thetaDeg = theta / deg2rad;
      Q2 = 4.0 * Ei * Ef * (TMath::Sin(theta / 2.0) * TMath::Sin(theta / 2.0)); //GeV^2
      nu = Ei - Ef; //GeV
      W2 = -Q2 + Mp * Mp + 2.0 * Mp * nu; // GeV^2
      W=0.;
      if (W2 > 0) W = TMath::Sqrt(W2); //GeV
    }
    //call the CS function
    if (W2 > 0) {
      sig_inelastic = pF1F209.GetXS(Z, A, Ei, Ef, theta); // ub/MeV-sr
      //Fill histograms
      // wfac is domega*denergy/n_thrown rad*MeV
      // lumin 1/ub for cur uA
      weight=sig_inelastic*lumin*wfac*cfac;
      hytar->Fill(hsytar,weight);
      
      hWw->Fill(W,weight);
      hWQ2->Fill(Q2,weight);
      hWom->Fill(nu,weight);
      hth->Fill(thetaDeg,weight);
      h_xs_ys->Fill(ys,xs,weight);
      h_xfp_yfp->Fill(hsxfp,hsyfp,weight);
      h_xfp_ypfp->Fill(hsxfp,hsypfp,weight);
      h_xfp_xpfp->Fill(hsxfp,hsxpfp,weight);
      h_yfp_ypfp->Fill(hsyfp,hsypfp,weight);
      h_yfp_xpfp->Fill(hsyfp,hsxpfp,weight);
      h_xpfp_ypfp->Fill(hsxpfp,hsypfp,weight);
      h_xpfp_delta->Fill(hsdelta,hsxpfp,weight);
		}
	}
  
  cout << " theta_spec = " << ts << " p_spec = " << p_spec << endl;
  //
  rate = hWQ2->Integral();
  t = new TText(.5,.5,Form("Inelastic MC rate: %f Hz",rate));
  t->SetTextAlign(22);

	if (rate > 0) {
    TCanvas *c = new TCanvas("c", "c", 800, 1200);
    c->Divide(2,2);
    c->cd(1);
    hWw->Draw();
    c->cd(2);
    hWQ2->Draw();
    cout << "Inelastic  MC rate = " << rate << " Hz" << endl;
    c->cd(3);
    hth->Draw();
    c->cd(4);
    hWom->Draw();
    outputpdf="inelastic_carbon/"+basename+"_kin.pdf";
    c->Print(outputpdf+"(");

    c->Clear();
    t->Draw();
    c->Print(outputpdf+")");
    delete c;
    //hom_dat->Draw("same");
    //
    //
    TCanvas *cfp = new TCanvas("cfp","Focal plane ",1400,900);
    cfp->Divide(2,3);
    gStyle->SetGridStyle(1);
    cfp->cd(1);
    gPad->SetLogz();
    gPad->SetGridx();
    gPad->SetGridy();
    h_xfp_yfp->Draw("colz");
    HList.Add(h_xfp_yfp);
    title->Draw();
    cfp->cd(2);
    h_xfp_xpfp->Draw("colz");
    HList.Add(h_xfp_xpfp);
    cfp->cd(3);
    h_xpfp_ypfp->Draw("colz");
    HList.Add(h_xpfp_ypfp);
    cfp->cd(4);
    h_yfp_ypfp->Draw("colz");
    HList.Add(h_yfp_ypfp);
    cfp->cd(5);
    h_yfp_xpfp->Draw("colz");
    HList.Add(h_yfp_xpfp);
    cfp->cd(6);
    h_xfp_ypfp->Draw("colz");
    HList.Add(h_xfp_ypfp);
    outputpdf="inelastic_carbon/"+basename+".pdf";
    cfp->Print(outputpdf+"(");

    cfp->Clear();
    t->Draw();
    cfp->Print(outputpdf+")");
    delete cfp;
	}
  //
  //
  TCanvas *cytar = new TCanvas("cytar", "cytar", 800, 1200);
  cytar->Divide(1,1);
  cytar->cd(1);
  hytar->Draw();
  outputpdf="inelastic_carbon/"+basename+"_ytar.pdf";
  cytar->Print(outputpdf+"(");
  cytar->cd(2);
  cytar->Clear();
  t->Draw();
  cytar->Print(outputpdf+")");
  delete cytar;
  //
  TFile hsimc(outputhist,"recreate");
  HList.Write();
  cout << " Plotted histograms put in root file = " << outputhist << endl;
}