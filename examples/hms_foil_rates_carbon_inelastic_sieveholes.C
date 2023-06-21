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

void hms_foil_rates_carbon_inelastic_sieveholes(TString basename="temp",double cur=70){
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
   outputhist="worksim/sieveHolesHists/"+basename+"_sieve_hist.root";
   TObjArray HList(0);
   TString outputpdf1="inelastic_carbon/sieveData/"+basename+"_sieve_kin.pdf";
   TString outputpdf2="inelastic_carbon/sieveData/"+basename+"_sieve.pdf";;
   TString outputpdf3="inelastic_carbon/sieveData/"+basename+"_sieve_ytar.pdf";
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
   Float_t xsnum;
   Float_t ysnum;

   
   // Set branch addresses.
   tsimc->SetBranchAddress("xc_sieve",&xs);
   tsimc->SetBranchAddress("yc_sieve",&ys);
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
   tsimc->SetBranchAddress("xs_num",&xsnum);
   tsimc->SetBranchAddress("ys_num",&ysnum);

   //Declare maximum indexes for sieve location.
   int xmax = 9;
   int ymax = 9;
   Double_t x_r=40.;
   Double_t y_r= 40.;
   Double_t xp_r=.1;
   Double_t yp_r=.04;
   TH1F *hytar[xmax][ymax];
   TH1F *hWw[xmax][ymax];
   TH1F *hWQ2[xmax][ymax];
   TH1F* hWom[xmax][ymax];
   TH1F* hth[xmax][ymax];

   TH2F* h_xs_ys[xmax][ymax];
   TH2F* h_xfp_yfp[xmax][ymax];
   TH2F* h_xfp_ypfp[xmax][ymax];
   TH2F* h_xfp_xpfp[xmax][ymax];
   TH2F* h_yfp_ypfp[xmax][ymax];
   TH2F* h_yfp_xpfp[xmax][ymax];
   TH2F* h_xpfp_ypfp[xmax][ymax];
   TH2F* h_xpfp_delta[xmax][ymax];

   TH2F* h_ytar_delta[xmax][ymax];
   TH2F* h_ytar_yptar[xmax][ymax];
   TH2F* h_ztar_yptar_all[xmax][ymax];
   TH1F* h_ytar[xmax][ymax];

   //define simulation histograms
   for(int x = 0; x < xmax; x++){
    for(int y = 0; y < ymax; y++){
        hytar[x][y] = new TH1F(Form("hytar_x%i_y%i", x-4, y-4), Form("hytar_x%i_y%i", x-4, y-4), 100, -5., 5.);
        
        hWw[x][y]= new TH1F(Form("hWw_x%i_y%i",x-4,y-4), Form("hWw_x%d_y%d",x-4,y-4), 100, 0, 5.0);
        hWw[x][y]->GetXaxis()->SetTitle("W [GeV]");
        hWw[x][y]->GetYaxis()->SetTitle("counts/sec");
        
        hWQ2[x][y] = new TH1F(Form("hWQ2_x%i_y%i",x-4,y-4), Form("hWQ2_x%i_y%i",x-4,y-4), 100, 0, 7);
        hWQ2[x][y]->GetXaxis()->SetTitle("Q^2 [GeV^2]");
        hWQ2[x][y]->GetYaxis()->SetTitle("counts/sec");
        
        hWom[x][y] = new TH1F(Form("hWom_x%i_y%i",x-4,y-4), Form("hWom_x%i_y%i",x-4,y-4),100, 0, 5.0);
        hWom[x][y]->GetXaxis()->SetTitle("Ei - Ef [GeV]");
        hWom[x][y]->GetYaxis()->SetTitle("counts/sec");
        
        hth[x][y] = new TH1F(Form("hth_x%i_y%i",x-4,y-4), Form("hth_x%i_y%i",x-4,y-4),160, 0, 40);
        hth[x][y]->GetXaxis()->SetTitle("theta [deg]");
        hth[x][y]->GetYaxis()->SetTitle("counts/sec");

        
        h_xs_ys[x][y] = new TH2F(Form("h_xs_ys_x%i_y%i",x-4,y-4),";x_sieve;y_sieve",100,-15,15,100,-15,15);
        h_xfp_yfp[x][y] = new TH2F(Form("h_xfp_yfp_x%i_y%i",x-4,y-4),";x_fp;y_fp",100,-x_r,x_r,100,-y_r,y_r);
        h_xfp_ypfp[x][y] = new TH2F(Form("h_xfp_ypfp_x%i_y%i",x-4,y-4),";x_fp;y_pfp",100,-x_r,x_r,100,-yp_r,yp_r);
        h_xfp_xpfp[x][y] = new TH2F(Form("h_xfp_xpfp_x%i_y%i",x-4,y-4),";x_fp;xp_fp",100,-x_r,x_r,100,-xp_r,xp_r);
        h_yfp_ypfp[x][y] = new TH2F(Form("h_yfp_ypfp_x%i_y%i",x-4,y-4),";y_fp;yp_fp",100,-y_r,y_r,100,-yp_r,yp_r);
        h_yfp_xpfp[x][y] = new TH2F(Form("h_yfp_xpfp_x%i_y%i",x-4,y-4),";y_fp;xp_fp",100,-y_r,y_r,100,-xp_r,xp_r);
        h_xpfp_ypfp[x][y] = new TH2F(Form("h_xpfp_ypfp_x%i_y%i",x-4,y-4),";xp_fp;yp_fp",100,-xp_r,xp_r,100,-yp_r,yp_r);
        h_xpfp_delta[x][y] = new TH2F(Form("h_xpfp_delta_x%i_y%i",x-4,y-4),"Inelastic ;delta;xp_fp",100,-10,10,100,-0.1,+0.1);

        h_ytar_delta[x][y] = new TH2F(Form("h_ytar_delta_x%i_y%i",x-4,y-4),"Inelastic ;delta;ytar",100,-10,10,100,-10.,+10.);
        h_ytar_yptar[x][y] = new TH2F(Form("h_ytar_yptar_x%i_y%i",x-4,y-4),";yptar;ytar",100,-.05,.05,100,-10.,+10.);
        h_ztar_yptar_all[x][y] = new TH2F(Form("h_ztar_yptar_all_x%i_y%i",x-4,y-4),";ztar:yptar",100,-25.,+25.,100,-.05,.05);
        h_ytar[x][y] =  new TH1F(Form("h_ytar_x%i_y%i",x-4,y-4)," ;ytar",100,-5.,5.);
    }
   }
   
   
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
   Float_t cfac;
   Float_t weight;
   //
   cfac=1.;
   thick= 0.044; // foil thickness g/cm2 in multifoil
   thick= 0.1749; // 0.5% single carbon
   lumin= thick*cur/A*N_A/Q_E*1e-36;// lumin 1/ub for cur uA
   
   Double_t rate[xmax][ymax];
   Double_t time[xmax][ymax];
   TText *t;
   TText *t2;
   TText *t3;
   
   for (int i = 0; i < nentries; i++) {
        tsimc->GetEntry(i);
        for(int x2 = 0; x2 < xmax; x2++){
            for(int y2 = 0; y2 < ymax; y2++){
                if(xsnum == x2-4 && ysnum == y2-4){
                    // Define kinematics
                    Ef = p_spec * (1.0 + 0.01*hsdelta); //scattered electron energy //GeV
                    nu = Ei - Ef; //GeV
                    W2=0;
                    h_ytar[x2][y2]->Fill(hsytari);
                    h_ztar_yptar_all[x2][y2]->Fill(hsztari,hsyptari);      
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
                        hytar[x2][y2]->Fill(hsytar,weight);
                        
                        hWw[x2][y2]->Fill(W,weight);
                        hWQ2[x2][y2]->Fill(Q2,weight);
                        hWom[x2][y2]->Fill(nu,weight);
                        hth[x2][y2]->Fill(thetaDeg,weight);
                        h_xs_ys[x2][y2]->Fill(ys,xs,weight);
                        h_xfp_yfp[x2][y2]->Fill(hsxfp,hsyfp,weight);
                        h_xfp_ypfp[x2][y2]->Fill(hsxfp,hsypfp,weight);
                        h_xfp_xpfp[x2][y2]->Fill(hsxfp,hsxpfp,weight);
                        h_yfp_ypfp[x2][y2]->Fill(hsyfp,hsypfp,weight);
                        h_yfp_xpfp[x2][y2]->Fill(hsyfp,hsxpfp,weight);
                        h_xpfp_ypfp[x2][y2]->Fill(hsxpfp,hsypfp,weight);
                        h_xpfp_delta[x2][y2]->Fill(hsdelta,hsxpfp,weight);
                    }
                }
            }
        }
    }
    
    cout << " theta_spec = " << ts << " p_spec = " << p_spec << endl;
    //
    TCanvas *c = new TCanvas("c", "c", 800, 1200);
    TCanvas *cfp = new TCanvas("cfp","Focal plane ",1400,900);
    TCanvas *cytar = new TCanvas("cytar", "cytar", 800, 1200);
    int check = 0;

    for(int x3 = 0; x3 < xmax; x3++){
        for(int y3 = 0; y3 < ymax; y3++){
            
            rate[x3][y3] = hWQ2[x3][y3]->Integral();
            if(rate[x3][y3] != 0){
                time[x3][y3] = 200.0/rate[x3][y3];
            }
            else{
                time[x3][y3] = 0;
            }
            t = new TText(0.5,.5,Form("Inelastic MC rate: %f Hz",rate[x3][y3]));
            t->SetTextAlign(22);
            t2 = new TText(0.5, 0.6, Form("for xsnum = %i and ysnum = %i", x3-4, y3-4));
            t2->SetTextAlign(22);
            t3 = new TText(0.5,.4,Form("Time to 200 Events: %f s", time[x3][y3]));
            t3->SetTextAlign(22);

            if (rate[x3][y3] > 0) {
                
                check++;
                
                c->Clear();
                c->Divide(2,2);
                c->cd(1);
                hWw[x3][y3]->Draw();
                
                c->cd(2);
                hWQ2[x3][y3]->Draw();
                cout << "Inelastic  MC rate = " << rate[x3][y3] << " Hz for xsnum = " << x3-4 << " and ysnum = " << y3-3 << endl;
                cout << "Time to 200 events = " << time[x3][y3] << " s" << endl;
                c->cd(3);
                hth[x3][y3]->Draw();
                c->cd(4);
                hWom[x3][y3]->Draw();
                
                if(check == 1){//x3==0&&y3==0){
                    c->SaveAs(outputpdf1+"(");
                }
                else{
                    c->SaveAs(outputpdf1);
                }

                c->Clear();
                t3->Draw();
                t->Draw("SAME");
                t2->Draw("SAME");
                c->SaveAs(outputpdf1);

                //hom_dat->Draw("same");
                //
                //
                cfp->Clear();
                cfp->Divide(2,3);
                gStyle->SetGridStyle(1);
                cfp->cd(1);
                gPad->SetLogz();
                gPad->SetGridx();
                gPad->SetGridy();
                h_xfp_yfp[x3][y3]->Draw("colz");
                HList.Add(h_xfp_yfp[x3][y3]);
                title->Draw();
                cfp->cd(2);
                h_xfp_xpfp[x3][y3]->Draw("colz");
                HList.Add(h_xfp_xpfp[x3][y3]);
                cfp->cd(3);
                h_xpfp_ypfp[x3][y3]->Draw("colz");
                HList.Add(h_xpfp_ypfp[x3][y3]);
                cfp->cd(4);
                h_yfp_ypfp[x3][y3]->Draw("colz");
                HList.Add(h_yfp_ypfp[x3][y3]);
                cfp->cd(5);
                h_yfp_xpfp[x3][y3]->Draw("colz");
                HList.Add(h_yfp_xpfp[x3][y3]);
                cfp->cd(6);
                h_xfp_ypfp[x3][y3]->Draw("colz");
                HList.Add(h_xfp_ypfp[x3][y3]);
                
                
                if(check == 1){//x3==0 && y3==0){
                    cfp->SaveAs(outputpdf2+"(");
                }
                else{
                    cfp->SaveAs(outputpdf2);
                }
                
                
                cfp->Clear();
                t3->Draw();
                t->Draw("SAME");
                t2->Draw("SAME");
                cfp->SaveAs(outputpdf2);
                
            }
            //
            //
            cytar->Clear();
            cytar->Divide(1,1);
            cytar->cd(1);
            hytar[x3][y3]->Draw();

            if(x3==0 && y3==0){
                cytar->SaveAs(outputpdf3+"(");
            }
            else{
                cytar->SaveAs(outputpdf3);
            }
            cytar->cd(2);
            cytar->Clear();
            t3->Draw();
            t->Draw("SAME");
            t2->Draw("SAME");
            cytar->SaveAs(outputpdf3);
            //
        }
    }

    Double_t lowRate = 1000000;
    Double_t lowTime;
    int lowX;
    int lowY;

    for(int x5 = 0; x5 < xmax; x5++){
        for(int y5 = 0; y5 < ymax; y5++){
            if(rate[x5][y5] != 0 && rate[x5][y5] - lowRate < 0){
                lowRate = rate[x5][y5];
                lowTime = time[x5][y5];
                lowX = x5 - 4;
                lowY = y5 - 4;
            }
        }
    }

    TText* t4 = new TText(0.5,.5,Form("Lowest Inelastic MC rate: %f Hz",lowRate));
    t4->SetTextAlign(22);
    TText* t5 = new TText(0.5, 0.6, Form("Found in xsnum = %i and ysnum = %i", lowX, lowY));
    t5->SetTextAlign(22);
    TText* t6 = new TText(0.5,.4,Form("Low Time to 200 Events: %f s", lowTime));
    t6->SetTextAlign(22);

    if(check > 0){
        c->Clear();
        t6->Draw();
        t4->Draw("SAME");
        t5->Draw("SAME");
        c->SaveAs(outputpdf1+")");
        
        cfp->Clear();
        t6->Draw();
        t4->Draw("SAME");
        t5->Draw("SAME");
        cfp->SaveAs(outputpdf2+")");
    }

    cytar->Clear();
    t6->Draw();
    t4->Draw("SAME");
    t5->Draw("SAME");
    cytar->SaveAs(outputpdf3+")");


    delete c;
    delete cfp;
    delete cytar;
    TFile hsimc(outputhist,"recreate");
    HList.Write();
    cout << " Plotted histograms put in root file = " << outputhist << endl;
    
    for(int x4 = 0; x4 < xmax; x4++){
        for(int y4 = 0; y4 < ymax; y4++){
            delete hytar[x4][y4];
            delete hWw[x4][y4];
            delete hWQ2[x4][y4];
            delete hWom[x4][y4];
            delete hth[x4][y4];
            delete h_xs_ys[x4][y4];
            delete h_xfp_yfp[x4][y4];
            delete h_xfp_ypfp[x4][y4];
            delete h_xfp_xpfp[x4][y4];
            delete h_yfp_ypfp[x4][y4];
            delete h_yfp_xpfp[x4][y4];
            delete h_xpfp_ypfp[x4][y4];
            delete h_xpfp_delta[x4][y4];
            delete h_ytar_delta[x4][y4];
            delete h_ytar_yptar[x4][y4];
            delete h_ztar_yptar_all[x4][y4];
            delete h_ytar[x4][y4];
        }
    }
}