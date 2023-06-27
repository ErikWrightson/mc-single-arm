/**
 * Originally adapted from code by Mark Jones (jones) and modified to allow carbon multifoil use
 * with a sieve in place, and intaking different currents.
 * @author Erik Wrightson (wrightso)
 * @created June 2023
 * @version 06.27.2023
 */

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

/**
 * This macro runs gets the rate of events per sieve hole on a carbon optics target in the HMS Spectrometer in
 * Hall-C at JLab with the optics sieve in place. Must be linked with the object file libF1F209.so found in the pbmodel directory as it uses a
 * wrapper class to get the inelastic cross-section  of each event. This will also plot various kinematics for events
 * passing through each hole originating from each foil.
 *
 * @note There is no delta cut being made and the cuts for the sieve holes are done asssuming
 * assuming perfect hole ID using the initial target values for each event. Conversion for using with real data
 * should be done by making cuts on reconstructed variables.
 *
 * @param basename - name of the root file that is assumed to be in the worksim directory
 * @param cur - current to calculate the rates at in uA
 * @param foilSep - seperation of the non-z=0 foils in cm
 * @param threeFoil - flag for a threefoil target (true = there is a z=0 foil)
 */
void hms_foil_rates_carbon_inelastic_MSH(TString basename="temp",double cur=70, double foilSep=8, bool threeFoil=false){

    //Set the number of foils for this analysis. (Default to two.)
    int nfoil = 2;
    if(threeFoil){
        nfoil = 3;
    }

    //Prompt the user for the name of the root file if not already provided.
   if (basename=="temp") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }

   //Set global options.
   gStyle->SetPalette(1,0);
   gStyle->SetOptStat(1);
   gStyle->SetOptFit(11);
   gStyle->SetTitleOffset(1.,"Y");
   gStyle->SetTitleOffset(.7,"X");
   gStyle->SetLabelSize(0.04,"XY");
   gStyle->SetTitleSize(0.06,"XY");
   gStyle->SetPadLeftMargin(0.12);

   //Declare paths input and paths for storing output.
   TString inputroot;
   inputroot="worksim/"+basename+".root";
   TString outputhist;
   outputhist="worksim/multifoilHolesHists/"+basename+"_MSH_hist.root";
   TObjArray HList(0);
   TString outputpdf1="inelastic_carbon/multifoilHoles/"+basename+"_MSH_kin.pdf";
   TString outputpdf2="inelastic_carbon/multifoilHoles/"+basename+"_MSH.pdf";;
   TString outputpdf3="inelastic_carbon/multifoilHoles/"+basename+"_MSH_ytar.pdf";
   TString htitle=basename;
   TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
   //  gSystem->Load("pbmodel/libF1F209.so");
   
   //Open the rootfile and retrieve its tree.
   TFile *fsimc = new TFile(inputroot); 
   TTree *tsimc = (TTree*) fsimc->Get("h1");

   //Declare variables to store each entry from the hms ntuple.   
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

   //Declare all histogram arrays with enough room for each foil and all sieve holes to seperate by hole.
   TH1F *hytar[nfoil][xmax][ymax];
   TH1F *hWw[nfoil][xmax][ymax];
   TH1F *hWQ2[nfoil][xmax][ymax];
   TH1F* hWom[nfoil][xmax][ymax];
   TH1F* hth[nfoil][xmax][ymax];

   TH2F* h_xs_ys[nfoil][xmax][ymax];
   TH2F* h_xfp_yfp[nfoil][xmax][ymax];
   TH2F* h_xfp_ypfp[nfoil][xmax][ymax];
   TH2F* h_xfp_xpfp[nfoil][xmax][ymax];
   TH2F* h_yfp_ypfp[nfoil][xmax][ymax];
   TH2F* h_yfp_xpfp[nfoil][xmax][ymax];
   TH2F* h_xpfp_ypfp[nfoil][xmax][ymax];
   TH2F* h_xpfp_delta[nfoil][xmax][ymax];

   TH2F* h_ytar_delta[nfoil][xmax][ymax];
   TH2F* h_ytar_yptar[nfoil][xmax][ymax];
   TH2F* h_ztar_yptar_all[nfoil][xmax][ymax];
   TH1F* h_ytar[nfoil][xmax][ymax];

   //Instantiate histograms.
   for(int n = 0; n < nfoil; n++){
        for(int x = 0; x < xmax; x++){
            for(int y = 0; y < ymax; y++){
                hytar[n][x][y] = new TH1F(Form("hytar_f%i_x%i_y%i",n,x-4, y-4), Form("hytar_f%i_x%i_y%i",n,x-4, y-4), 100, -5., 5.);
                
                hWw[n][x][y]= new TH1F(Form("hWw_f%i_x%i_y%i",n,x-4,y-4), Form("hWw_f%i_x%d_y%d",n,x-4,y-4), 100, 0, 5.0);
                hWw[n][x][y]->GetXaxis()->SetTitle("W [GeV]");
                hWw[n][x][y]->GetYaxis()->SetTitle("counts/sec");
                
                hWQ2[n][x][y] = new TH1F(Form("hWQ2_f%i_x%i_y%i",n,x-4,y-4), Form("hWQ2_f%i_x%i_y%i",n,x-4,y-4), 100, 0, 7);
                hWQ2[n][x][y]->GetXaxis()->SetTitle("Q^2 [GeV^2]");
                hWQ2[n][x][y]->GetYaxis()->SetTitle("counts/sec");
                
                hWom[n][x][y] = new TH1F(Form("hWom_f%i_x%i_y%i",n,x-4,y-4), Form("hWom_f%i_x%i_y%i",n,x-4,y-4),100, 0, 5.0);
                hWom[n][x][y]->GetXaxis()->SetTitle("Ei - Ef [GeV]");
                hWom[n][x][y]->GetYaxis()->SetTitle("counts/sec");
                
                hth[n][x][y] = new TH1F(Form("hth_f%i_x%i_y%i",n,x-4,y-4), Form("hth_f%i_x%i_y%i",n,x-4,y-4),160, 0, 40);
                hth[n][x][y]->GetXaxis()->SetTitle("theta [deg]");
                hth[n][x][y]->GetYaxis()->SetTitle("counts/sec");
                
                h_xs_ys[n][x][y] = new TH2F(Form("h_xs_ys_f%i_x%i_y%i",n,x-4,y-4),";x_sieve;y_sieve",100,-15,15,100,-15,15);
                h_xfp_yfp[n][x][y] = new TH2F(Form("h_xfp_yfp_f%i_x%i_y%i",n,x-4,y-4),";x_fp;y_fp",100,-x_r,x_r,100,-y_r,y_r);
                h_xfp_ypfp[n][x][y] = new TH2F(Form("h_xfp_ypfp_f%i_x%i_y%i",n,x-4,y-4),";x_fp;y_pfp",100,-x_r,x_r,100,-yp_r,yp_r);
                h_xfp_xpfp[n][x][y] = new TH2F(Form("h_xfp_xpfp_f%i_x%i_y%i",n,x-4,y-4),";x_fp;xp_fp",100,-x_r,x_r,100,-xp_r,xp_r);
                h_yfp_ypfp[n][x][y] = new TH2F(Form("h_yfp_ypfp_f%i_x%i_y%i",n,x-4,y-4),";y_fp;yp_fp",100,-y_r,y_r,100,-yp_r,yp_r);
                h_yfp_xpfp[n][x][y] = new TH2F(Form("h_yfp_xpfp_f%i_x%i_y%i",n,x-4,y-4),";y_fp;xp_fp",100,-y_r,y_r,100,-xp_r,xp_r);
                h_xpfp_ypfp[n][x][y] = new TH2F(Form("h_xpfp_ypfp_f%i_x%i_y%i",n,x-4,y-4),";xp_fp;yp_fp",100,-xp_r,xp_r,100,-yp_r,yp_r);
                h_xpfp_delta[n][x][y] = new TH2F(Form("h_xpfp_delta_f%i_x%i_y%i",n,x-4,y-4),"Inelastic ;delta;xp_fp",100,-10,10,100,-0.1,+0.1);
                
                h_ytar_delta[n][x][y] = new TH2F(Form("h_ytar_delta_f%i_x%i_y%i",n,x-4,y-4),"Inelastic ;delta;ytar",100,-10,10,100,-10.,+10.);
                h_ytar_yptar[n][x][y] = new TH2F(Form("h_ytar_yptar_f%i_x%i_y%i",n,x-4,y-4),";yptar;ytar",100,-.05,.05,100,-10.,+10.);
                h_ztar_yptar_all[n][x][y] = new TH2F(Form("h_ztar_yptar_all_f%i_x%i_y%i",n,x-4,y-4),";ztar:yptar",100,-25.,+25.,100,-.05,.05);
                h_ytar[n][x][y] =  new TH1F(Form("h_ytar_f%i_x%i_y%i",n,x-4,y-4)," ;ytar",100,-5.,5.);
            }
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
   
   //Get the beam energy, spectrometer angle and momentum from the first entry of the ntuple (assumed to be the same throughout).
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
   //Included by linking with files in pbmodel folder.
   F1F209Wrapper pF1F209;
   Float_t cfac;
   Float_t weight;
   //
   cfac=1.;
   thick= 0.044; // foil thickness g/cm2 in multifoil
   //thick= 0.1749; // 0.5% single carbon
   lumin= thick*cur/A*N_A/Q_E*1e-36;// lumin 1/ub for cur uA, 1e-36 corrects for change of unit size
   
   //Store the rate for each sieve hole sourced from each foil.
   Double_t rate[nfoil][xmax][ymax];
   Double_t time[nfoil][xmax][ymax];
   //Use TText variables to print out any single data points that may be useful.
   TText *t;
   TText *t2;
   TText *t3;
   
   for (int i = 0; i < nentries; i++) {
        tsimc->GetEntry(i);
        for(int x2 = 0; x2 < xmax; x2++){
            for(int y2 = 0; y2 < ymax; y2++){
                
                //Sort for each whole that each event came through. (THIS USES THE INITIAL VARIABLES NOT THE RECONSTRUCTED)
                //Ignore holes that are not actually there (refer to sieve schematic).
                //Central hole is smaller than other but not implemented since central hole is never really a limiting factor.
                if(xsnum == x2-4 && ysnum == y2-4 && !(xsnum == -1 && ysnum == 1) && !(xsnum == 1 && ysnum == -1)){
                    
                    // Define kinematics
                    Ef = p_spec * (1.0 + 0.01*hsdelta); //scattered electron energy //GeV
                    nu = Ei - Ef; //GeV
                    W2=0; //W^2
                    //if final energy is less than initial calculcate the other kinematic values for this event.
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
                    if (W2 > 0){
                        //Get inelastic cross section for this event.
                        sig_inelastic = pF1F209.GetXS(Z, A, Ei, Ef, theta); // ub/MeV-sr
                        
                        //Weight the data by the cross-section times the luminosity and weighting factor.
                        // wfac is domega*denergy/n_thrown rad*MeV
                        // lumin 1/ub for cur uA
                        weight=sig_inelastic*lumin*wfac*cfac;
                    }
                    
                    //Check if this event originated from the -z foil.
                    //USES INITIAL VALUES. ONLY VALID ASSUMING PERFECT IDENTIFICATION OF SOURCE (or using the mc-single-arm simulator).
                    //USE RECONSTRUCTED VALUES FOR CUTTING REAL DATA.
                    if(hsztari >= -1*foilSep-(thick/2) && hsztari <= -1*foilSep+(thick/2)){
                        h_ytar[0][x2][y2]->Fill(hsytari);
                        h_ztar_yptar_all[0][x2][y2]->Fill(hsztari,hsyptari);

                        if (W2 > 0) {
                            //Fill histograms
                            hytar[0][x2][y2]->Fill(hsytar,weight);
                            
                            hWw[0][x2][y2]->Fill(W,weight);
                            hWQ2[0][x2][y2]->Fill(Q2,weight);
                            hWom[0][x2][y2]->Fill(nu,weight);
                            hth[0][x2][y2]->Fill(thetaDeg,weight);
                            h_xs_ys[0][x2][y2]->Fill(ys,xs,weight);
                            h_xfp_yfp[0][x2][y2]->Fill(hsxfp,hsyfp,weight);
                            h_xfp_ypfp[0][x2][y2]->Fill(hsxfp,hsypfp,weight);
                            h_xfp_xpfp[0][x2][y2]->Fill(hsxfp,hsxpfp,weight);
                            h_yfp_ypfp[0][x2][y2]->Fill(hsyfp,hsypfp,weight);
                            h_yfp_xpfp[0][x2][y2]->Fill(hsyfp,hsxpfp,weight);
                            h_xpfp_ypfp[0][x2][y2]->Fill(hsxpfp,hsypfp,weight);
                            h_xpfp_delta[0][x2][y2]->Fill(hsdelta,hsxpfp,weight);
                        }
                    }
                    else if(hsztari >= foilSep-(thick/2) && hsztari <= foilSep+(thick/2)){
                        //Check if this event originated from the +z foil.
                        //USES INITIAL VALUES. ONLY VALID ASSUMING PERFECT IDENTIFICATION OF SOURCE (or using the mc-single-arm simulator)
                        //USE RECONSTRUCTED VALUES FOR CUTTING REAL DATA.
                        h_ytar[nfoil-1][x2][y2]->Fill(hsytari);
                        h_ztar_yptar_all[nfoil-1][x2][y2]->Fill(hsztari,hsyptari);

                        if (W2 > 0) {                            
                            //Fill histograms
                            hytar[nfoil-1][x2][y2]->Fill(hsytar,weight);
                            
                            hWw[nfoil-1][x2][y2]->Fill(W,weight);
                            hWQ2[nfoil-1][x2][y2]->Fill(Q2,weight);
                            hWom[nfoil-1][x2][y2]->Fill(nu,weight);
                            hth[nfoil-1][x2][y2]->Fill(thetaDeg,weight);
                            h_xs_ys[nfoil-1][x2][y2]->Fill(ys,xs,weight);
                            h_xfp_yfp[nfoil-1][x2][y2]->Fill(hsxfp,hsyfp,weight);
                            h_xfp_ypfp[nfoil-1][x2][y2]->Fill(hsxfp,hsypfp,weight);
                            h_xfp_xpfp[nfoil-1][x2][y2]->Fill(hsxfp,hsxpfp,weight);
                            h_yfp_ypfp[nfoil-1][x2][y2]->Fill(hsyfp,hsypfp,weight);
                            h_yfp_xpfp[nfoil-1][x2][y2]->Fill(hsyfp,hsxpfp,weight);
                            h_xpfp_ypfp[nfoil-1][x2][y2]->Fill(hsxpfp,hsypfp,weight);
                            h_xpfp_delta[nfoil-1][x2][y2]->Fill(hsdelta,hsxpfp,weight);
                        }
                    }
                    else if(threeFoil && hsztari >= -(thick/2) && hsztari <= (thick/2)){
                        //Check if this event originated from the z=0 foil in the case of three foils.
                        //USES INITIAL VALUES. ONLY VALID ASSUMING PERFECT IDENTIFICATION OF SOURCE (or using the mc-single-arm simulator).
                        //USE RECONSTRUCTED VALUES FOR CUTTING REAL DATA.
                        
                        h_ytar[1][x2][y2]->Fill(hsytari);
                        h_ztar_yptar_all[1][x2][y2]->Fill(hsztari,hsyptari);

                        if (W2 > 0) {
                            //Fill histograms
                            hytar[1][x2][y2]->Fill(hsytar,weight);
                            
                            hWw[1][x2][y2]->Fill(W,weight);
                            hWQ2[1][x2][y2]->Fill(Q2,weight);
                            hWom[1][x2][y2]->Fill(nu,weight);
                            hth[1][x2][y2]->Fill(thetaDeg,weight);
                            h_xs_ys[1][x2][y2]->Fill(ys,xs,weight);
                            h_xfp_yfp[1][x2][y2]->Fill(hsxfp,hsyfp,weight);
                            h_xfp_ypfp[1][x2][y2]->Fill(hsxfp,hsypfp,weight);
                            h_xfp_xpfp[1][x2][y2]->Fill(hsxfp,hsxpfp,weight);
                            h_yfp_ypfp[1][x2][y2]->Fill(hsyfp,hsypfp,weight);
                            h_yfp_xpfp[1][x2][y2]->Fill(hsyfp,hsxpfp,weight);
                            h_xpfp_ypfp[1][x2][y2]->Fill(hsxpfp,hsypfp,weight);
                            h_xpfp_delta[1][x2][y2]->Fill(hsdelta,hsxpfp,weight);
                        }
                    }
                }
            }
        }
    }
    
    //Output the spectrometer angle and central momentum as a check that the data has been taken correctly.
    cout << " theta_spec = " << ts << " p_spec = " << p_spec << endl;
    //
    //Declare the canvases used to create the pdfs printed.
    TCanvas *c = new TCanvas("c", "c", 800, 1200);
    TCanvas *cfp = new TCanvas("cfp","Focal plane ",1400,900);
    TCanvas *cytar = new TCanvas("cytar", "cytar", 800, 1200);
    
    //Count how many holes have non-zero rates.
    int check = 0;

    for(int n2 = 0; n2 < nfoil; n2++){
        for(int x3 = 0; x3 < xmax; x3++){
            for(int y3 = 0; y3 < ymax; y3++){
                //Integrate over the weighted Q^2 to get the rate for each hole.
                rate[n2][x3][y3] = hWQ2[n2][x3][y3]->Integral();

                if(rate[n2][x3][y3] != 0){
                    //Get the time it takes to get to 200 events at this calculated rate.
                    time[n2][x3][y3] = 200.0/rate[n2][x3][y3];
                }
                else{
                    //If the rate is 0 avoid a divide by zero error and just 0 out the answer.
                    time[n2][x3][y3] = 0;
                }

                //Write out the rates and location for this hole.
                t = new TText(0.5,.5,Form("Inelastic MC rate: %f Hz",rate[n2][x3][y3]));
                t->SetTextAlign(22);
                t2 = new TText(0.5, 0.6, Form("for f-%i xsnum = %i and ysnum = %i",n2,x3-4, y3-4));
                t2->SetTextAlign(22);
                t3 = new TText(0.5,.4,Form("Time to 200 Events: %f s", time[n2][x3][y3]));
                t3->SetTextAlign(22);
                
                //Print the basic histograms for each hole.
                if (rate[n2][x3][y3] > 0) {
                    check++;
                    
                    c->Clear();
                    c->Divide(2,2);
                    c->cd(1);
                    hWw[n2][x3][y3]->Draw();
                    
                    c->cd(2);
                    hWQ2[n2][x3][y3]->Draw();
                    cout << "Inelastic  MC rate = " << rate[n2][x3][y3] << " Hz for foil = " << n2 << ", xsnum = " << x3-4 << ", and ysnum = " << y3-3 << endl;
                    cout << "Time to 200 events = " << time[n2][x3][y3] << " s" << endl;
                    
                    c->cd(3);
                    hth[n2][x3][y3]->Draw();
                    c->cd(4);
                    hWom[n2][x3][y3]->Draw();
                    
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
                    h_xfp_yfp[n2][x3][y3]->Draw("colz");
                    HList.Add(h_xfp_yfp[n2][x3][y3]);
                    title->Draw();
                    cfp->cd(2);
                    h_xfp_xpfp[n2][x3][y3]->Draw("colz");
                    HList.Add(h_xfp_xpfp[n2][x3][y3]);
                    cfp->cd(3);
                    h_xpfp_ypfp[n2][x3][y3]->Draw("colz");
                    HList.Add(h_xpfp_ypfp[n2][x3][y3]);
                    cfp->cd(4);
                    h_yfp_ypfp[n2][x3][y3]->Draw("colz");
                    HList.Add(h_yfp_ypfp[n2][x3][y3]);
                    cfp->cd(5);
                    h_yfp_xpfp[n2][x3][y3]->Draw("colz");
                    HList.Add(h_yfp_xpfp[n2][x3][y3]);
                    cfp->cd(6);
                    h_xfp_ypfp[n2][x3][y3]->Draw("colz");
                    HList.Add(h_xfp_ypfp[n2][x3][y3]);
                    
                    //If this is the first time through open the pdf.
                    if(check == 1){
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
                hytar[n2][x3][y3]->Draw();
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
    }
    

    Double_t lowRate = 1000000;
    Double_t lowTime;
    Double_t sumRate = 0;
    int lowX;
    int lowY;
    int lowF;

    //Declare histograms to store the rates of events sourced from each foil.
    TH2F *h_rates_foils[nfoil];
    if(threeFoil){
        h_rates_foils[0] = new TH2F("h_rates_f0","Rates from f=0;xsnum;ysnum",30,-5,5,30,-5,5);
        h_rates_foils[1] = new TH2F("h_rates_f1","Rates from f=1;xsnum;ysnum",30,-5,5,30,-5,5);
        h_rates_foils[nfoil-1] = new TH2F("h_rates_f2","Rates from f=2;xsnum;ysnum",30,-5,5,30,-5,5);
    }
    else{
        h_rates_foils[0] = new TH2F("h_rates_f0","Rates from f=0;xsnum;ysnum",30,-5,5,30,-5,5);
        h_rates_foils[nfoil-1] = new TH2F("h_rates_f1","Rates from f=1;xsnum;ysnum",30,-5,5,30,-5,5);
    }
    TH2F *h_rates_all = new TH2F("h_rates_all","Combined Rates from all foils;xsnum;ysnum",30,-5,5,30,-5,5);
    
    //Get the lowest non-zero rate for hole, its location and fill the rates histograms.
    for(int n4 = 0; n4 < nfoil; n4++){
        for(int x5 = 0; x5 < xmax; x5++){
            for(int y5 = 0; y5 < ymax; y5++){
                sumRate += rate[n4][x5][y5];
                if(rate[n4][x5][y5] != 0 && rate[n4][x5][y5] - lowRate < 0){
                    lowRate = rate[n4][x5][y5];
                    lowTime = time[n4][x5][y5];
                    lowX = x5 - 4;
                    lowY = y5 - 4;
                    lowF = n4;
                }
                h_rates_foils[n4]->Fill(x5-4,y5-4,rate[n4][x5][y5]);
            }
        }
    }
    
    int zeroCount = 0;
    int timeUnder3_0 = 0;
    int timeUnder3_1 = 0;
    int timeUnder3_2 = 0;
    int comboTimeUnder3 = 0;

    //Get the number of holes with a rate of 0 from all foils, the number of holes that would reach 200 counts within 3 hours,
    //and the combined rates for each hole.
    for(int x6 = 0; x6 < xmax; x6++){
        for (int y6 = 0; y6 < ymax; y6++){
            if(threeFoil){
                if(rate[0][x6][y6] == 0 && rate[1][x6][y6] == 0 && rate[nfoil-1][x6][y6]){
                    zeroCount++;
                }
                if(time[0][x6][y6] != 0 && time[0][x6][y6] < 10800) timeUnder3_0++;
                if(time[1][x6][y6] != 0 && time[1][x6][y6] < 10800) timeUnder3_1++;
                if(time[nfoil-1][x6][y6] != 0 &&time[nfoil-1][x6][y6] < 10800) timeUnder3_2++;
                if(200.0/(rate[0][x6][y6]+rate[1][x6][y6]+rate[nfoil-1][x6][y6]) < 10800) comboTimeUnder3++;
                h_rates_all->Fill(x6-4,y6-4,rate[0][x6][y6]+rate[1][x6][y6]+rate[nfoil-1][x6][y6]);
            }
            else{
                if(rate[0][x6][y6] == 0 && rate[nfoil-1][x6][y6]){
                    zeroCount++;
                }
                if(time[0][x6][y6] != 0 && time[0][x6][y6] < 10800) timeUnder3_0++;
                if(time[nfoil-1][x6][y6] != 0 && time[nfoil-1][x6][y6] < 10800) timeUnder3_1++;
                if(200.0/(rate[0][x6][y6]+rate[nfoil-1][x6][y6]) < 10800) comboTimeUnder3++;
                h_rates_all->Fill(x6-4,y6-4,rate[0][x6][y6]+rate[nfoil-1][x6][y6]);

            }
        }
    }

    gStyle->SetPaintTextFormat("4.3f");
    //Print the rates corresponding to each hole on each of the pdf files. Do this for each foil and the combined rates.
    if (check > 0){
        
        c->Clear();
        h_rates_foils[0]->Draw("COLZ");
        h_rates_foils[0]->Draw("SAMETEXT");
        c->SaveAs(outputpdf1);
        if(threeFoil){
            c->Clear();
            h_rates_foils[1]->Draw("COLZ");
            h_rates_foils[1]->Draw("SAMETEXT");
            c->SaveAs(outputpdf1);
        }
        c->Clear();
        h_rates_foils[nfoil-1]->Draw("COLZ");
        h_rates_foils[nfoil-1]->Draw("SAMETEXT");
        c->SaveAs(outputpdf1);

        c->Clear();
        h_rates_all->Draw("COLZ");
        h_rates_all->Draw("SAMETEXT");
        c->SaveAs(outputpdf1);

        cfp->Clear();
        h_rates_foils[0]->Draw("COLZ");
        h_rates_foils[0]->Draw("SAMETEXT");
        cfp->SaveAs(outputpdf2);
        if(threeFoil){
            cfp->Clear();
            h_rates_foils[1]->Draw("COLZ");
            h_rates_foils[1]->Draw("SAMETEXT");
            cfp->SaveAs(outputpdf2);
        }
        cfp->Clear();
        h_rates_foils[nfoil-1]->Draw("COLZ");
        h_rates_foils[nfoil-1]->Draw("SAMETEXT");
        cfp->SaveAs(outputpdf2);

        cfp->Clear();
        h_rates_all->Draw("COLZ");
        h_rates_all->Draw("SAMETEXT");
        cfp->SaveAs(outputpdf2);
    }

    cytar->Clear();
    h_rates_foils[0]->Draw("COLZ");
    h_rates_foils[0]->Draw("SAMETEXT");
    cytar->SaveAs(outputpdf3);
    if(threeFoil){
        cytar->Clear();
        h_rates_foils[1]->Draw("COLZ");
        h_rates_foils[1]->Draw("SAMETEXT");
        cytar->SaveAs(outputpdf3);
    }
    cytar->Clear();
    h_rates_foils[nfoil-1]->Draw("COLZ");
    h_rates_foils[nfoil-1]->Draw("SAMETEXT");
    cytar->SaveAs(outputpdf3);

    cytar->Clear();
    h_rates_all->Draw("COLZ");
    h_rates_all->Draw("SAMETEXT");  
    HList.Add(h_rates_all);
    cytar->SaveAs(outputpdf3); 

    int hrs = (int)(lowTime/3600);
    int min = (int)((lowTime-hrs*3600)/60);
    Double_t sec = lowTime-hrs*3600-min*60;

    TText* t4 = new TText(0.5,.5,Form("Lowest Inelastic MC rate: %f Hz",lowRate));
    t4->SetTextAlign(22);
    t4->SetY(2);

    TText* t5 = new TText(0.5, 0.6, Form("Found in f-%i xsnum = %i and ysnum = %i", lowF, lowX, lowY));
    t5->SetTextAlign(22);
    t5->SetY(3);

    TText* t6 = new TText(0.5,.4,Form("Low Time to 200 Events: %f s", lowTime));
    t6->SetTextAlign(22);
    t6->SetY(1);

    TText* t7 = new TText(0.5,0.2,Form("Summed rate: %f Hz", sumRate));
    t7->SetTextAlign(22);
    t7->SetY(-2);

    TText* t8 = new TText(0.1,0.3,Form("Low Time to 200 Events: %i hr %i min %2.2fs",hrs, min, sec ));
    t8->SetTextAlign(22);
    t8->SetY(0.5);
    
    TText* t9 = new TText(0.5,0.1,Form("Number of holes with 0 counts: %i holes", zeroCount));
    t9->SetTextAlign(22);
    t9->SetY(-2.5);

    TText* t10 = new TText(0.5,0.9,Form("Number of holes with 200 in 3hrs f=%i: %i holes", 0, timeUnder3_0));
    t10->SetTextAlign(22);
    t10->SetY(5.5);

    TText* t11 = new TText(0.5,0.8,Form("Number of holes with 200 in 3hrs f=%i: %i holes", 1, timeUnder3_1));
    t11->SetTextAlign(22);
    t11->SetY(4.5);

    TText* t12 = new TText(0.5,0.8,Form("Number of holes with 200 in 3hrs f=%i: %i holes", 2, timeUnder3_2));
    t12->SetTextAlign(22);
    t12->SetY(3.5);
    
    TText* t13 = new TText(0.5,0.8,Form("Holes with 200 in 3hrs combined rate: %i holes", comboTimeUnder3));
    t13->SetTextAlign(22);
    t13->SetY(-3.5);

    //Print any useful data out to the final page of each pdf
    if(check > 0){
        c->Clear();
        t6->Draw();
        t4->Draw("SAME");
        t5->Draw("SAME");
        t8->Draw("SAME");
        t7->Draw("SAME");
        t9->Draw("SAME");
        t10->Draw("SAME");
        t11->Draw("SAME");
        if(threeFoil) t12->Draw("SAME");
        t13->Draw("SAME");
        c->SaveAs(outputpdf1+")");
        
        cfp->Clear();
        t6->Draw();
        t4->Draw("SAME");
        t5->Draw("SAME");
        t8->Draw("SAME");
        t7->Draw("SAME");
        t9->Draw("SAME");
        t10->Draw("SAME");
        t11->Draw("SAME");
        if(threeFoil) t12->Draw("SAME");
        t13->Draw("SAME");
        cfp->SaveAs(outputpdf2+")");
    }

    cytar->Clear();
    t6->Draw();
    t4->Draw("SAME");
    t5->Draw("SAME");
    t8->Draw("SAME");
    t7->Draw("SAME");
    t9->Draw("SAME");
    t10->Draw("SAME");
    t11->Draw("SAME");
    if(threeFoil) t12->Draw("SAME");
    t13->Draw("SAME");
    cytar->SaveAs(outputpdf3+")");


    delete c;
    delete cfp;
    delete cytar;

    //Save all the filled histograms.
    TFile hsimc(outputhist,"recreate");
    HList.Write();
    cout << " Plotted histograms put in root file = " << outputhist << endl;
    
    //Delete pointers.
    delete h_rates_all;
    delete h_rates_foils[0];
    if(threeFoil) delete h_rates_foils[1];
    delete h_rates_foils[nfoil-1];

    for(int n3 = 0; n3 < nfoil; n3++){
        for(int x4 = 0; x4 < xmax; x4++){
            for(int y4 = 0; y4 < ymax; y4++){
                delete hytar[n3][x4][y4];
                delete hWw[n3][x4][y4];
                delete hWQ2[n3][x4][y4];
                delete hWom[n3][x4][y4];
                delete hth[n3][x4][y4];
                delete h_xs_ys[n3][x4][y4];
                delete h_xfp_yfp[n3][x4][y4];
                delete h_xfp_ypfp[n3][x4][y4];
                delete h_xfp_xpfp[n3][x4][y4];
                delete h_yfp_ypfp[n3][x4][y4];
                delete h_yfp_xpfp[n3][x4][y4];
                delete h_xpfp_ypfp[n3][x4][y4];
                delete h_xpfp_delta[n3][x4][y4];
                delete h_ytar_delta[n3][x4][y4];
                delete h_ytar_yptar[n3][x4][y4];
                delete h_ztar_yptar_all[n3][x4][y4];
                delete h_ytar[n3][x4][y4];
            }
        }
    }
}