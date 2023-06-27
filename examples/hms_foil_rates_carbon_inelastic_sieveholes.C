/**
 * Originally adapted from code by Mark Jones (jones) and modified to allow carbon
 * single foil use with a sieve in place, and intaking variable currents.
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
 * This macro runs gets the rate of events per sieve hole from scattering on a carbon optics target in the HMS Spectrometer
 * in Hall-C at JLab with the optics sieve in place. Must be linked with the object file libF1F209.so found in the pbmodel directory as it uses a
 * wrapper class to get the inelastic cross-section  of each event. This will also plot various kinematics for events
 * passing through each sieve hole. 
 *
 * @note There is no delta cut being made and the cuts for the sieve holes are done asssuming assuming perfect hole
 * ID using the initial target values for each event. Conversion for use with real data should be done by making cuts
 * on reconstructed variables.
 *
 * @param basename - name of the root file that is assumed to be in the worksim directory
 * @param cur - current to calculate the rates at in uA
 */
void hms_foil_rates_carbon_inelastic_sieveholes(TString basename="temp",double cur=70){
    
    //Get the name of a root file located in the worksim directory if one was not provided.
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

    //Declare paths to input and output files.
    TString inputroot;
    inputroot="worksim/"+basename+".root";
    TString outputhist;
    outputhist="worksim/sieveHolesHists/"+basename+"_sieve_hist.root";
    TObjArray HList(0);
    TString outputpdf1="inelastic_carbon/sieveData/"+basename+"_sieve_kin.pdf";
    TString outputpdf2="inelastic_carbon/sieveData/"+basename+"_sieve.pdf";;
    TString outputpdf3="inelastic_carbon/sieveData/"+basename+"_sieve_ytar.pdf";
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

    //Declare histograms for each sieve hole.
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
    
    //Instantiate histograms
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
    double thick; // target thickness g/cm2
    double run_time =515.;// secondes;
    double lumin ;// luminosity per ub 
    double Ef;
    double theta; 
    double sig_inelastic;
    double ex_calc;
    double sig_elas_calc;

    //Included by linking with object file in pbmodel folder.
    F1F209Wrapper pF1F209;
    Float_t cfac;
    Float_t weight;
    //
    cfac=1.;
    //thick= 0.044; // foil thickness g/cm2 in multifoil
    thick= 0.1749; // 0.5% single carbon
    lumin= thick*cur/A*N_A/Q_E*1e-36;// lumin 1/ub for cur uA
    
    //Store the rate for each sieve hole.
    Double_t rate[xmax][ymax];
    Double_t time[xmax][ymax];
    TText *t;
    TText *t2;
    TText *t3;
    
    for (int i = 0; i < nentries; i++) {
        tsimc->GetEntry(i);
        for(int x2 = 0; x2 < xmax; x2++){
            for(int y2 = 0; y2 < ymax; y2++){
                
                //Cut on the sieve num and store this entry in the corresponding matrix location.
                //Ignore holes that are not actually there (refer to sieve schematic).
                //Central hole is smaller than other but not implemented since central hole is never really a limiting factor.
                //USE RECONSTRUCTED SIEVE VARIABLES FOR REAL DATA ANALYSIS
                if(xsnum == x2-4 && ysnum == y2-4 && !(xsnum == -1 && ysnum == 1) && !(xsnum == 1 && ysnum == -1)){
                    
                    // Define kinematics
                    Ef = p_spec * (1.0 + 0.01*hsdelta); //scattered electron energy //GeV
                    nu = Ei - Ef; //GeV
                    W2=0;
                    h_ytar[x2][y2]->Fill(hsytari);
                    h_ztar_yptar_all[x2][y2]->Fill(hsztari,hsyptari);  
                    
                    //If final energy is less than initial calculcate the other kinematic values for this event.    
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
                        
                        //Get inelastic cross section for this event.
                        sig_inelastic = pF1F209.GetXS(Z, A, Ei, Ef, theta); // ub/MeV-sr
                        
                        
                        // wfac is domega*denergy/n_thrown rad*MeV
                        // lumin 1/ub for cur uA
                        //Weight the data by the cross-section times the luminosity and weighting factor.
                        weight=sig_inelastic*lumin*wfac*cfac;

                        //Fill histograms
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
    
    //Output the spectrometer angle and central momentum as a check that the data has been taken correctly.
    cout << " theta_spec = " << ts << " p_spec = " << p_spec << endl;
    //
    //Declare the canvases used to create the pdfs printed.
    TCanvas *c = new TCanvas("c", "c", 800, 1200);
    TCanvas *cfp = new TCanvas("cfp","Focal plane ",1400,900);
    TCanvas *cytar = new TCanvas("cytar", "cytar", 800, 1200);

    //Check how many holes have non-zero rates.
    int check = 0;

    for(int x3 = 0; x3 < xmax; x3++){
        for(int y3 = 0; y3 < ymax; y3++){

            //Integrate over the weighted Q^2 to get the rate for each hole.
            rate[x3][y3] = hWQ2[x3][y3]->Integral();

            //If this hole has a non-zero rate then get the amount of time needed to get 200 events
            //through this hole. 200 was the rate determined to be the minimum useful amount for an
            //optics calibration. Else just set the time to 0 to avoid a div by 0 error.
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

            //If the rate is higher than 0 plot the different kinematics for each hole.
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
    Double_t sumRate = 0;
    int lowX;
    int lowY;
    int zeroCount = 0;
    int timeUnder3 = 0;

    //Declare a histogram to hold the rates for each hole.
    TH2F *h_rates = new TH2F("h_rates",";xsnum;ysnum",30,-5,5,30,-5,5);

    //Get the summed rate from all the holes, the rate from the hole with the lowest rate,
    //its location, and its time to 200 counts. Fill the rates histogram.
    for(int x5 = 0; x5 < xmax; x5++){
        for(int y5 = 0; y5 < ymax; y5++){
            sumRate += rate[x5][y5];
            if(rate[x5][y5] != 0 && rate[x5][y5] - lowRate < 0){
                lowRate = rate[x5][y5];
                lowTime = time[x5][y5];
                lowX = x5 - 4;
                lowY = y5 - 4;
            }
            if(rate[x5][y5] == 0){
                zeroCount++;
            }
            if(time[x5][y5] != 0 &&time[x5][y5] < 10800){
                timeUnder3++;
            }
            h_rates->Fill(x5-4,y5-4,rate[x5][y5]);
        }
    }

    //Convert the lowTime in seconds to hrs:min:sec
    int hrs = (int)(lowTime/3600);
    int min = (int)((lowTime-hrs*3600)/60);
    Double_t sec = lowTime-hrs*3600-min*60;
    
    //Print the rates histograms.
    gStyle->SetPaintTextFormat("4.3f");
    if (check > 0){
        c->Clear();
        h_rates->Draw("COLZ");
        h_rates->Draw("SAMETEXT");
        c->SaveAs(outputpdf1);

        cfp->Clear();
        h_rates->Draw("COLZ");
        h_rates->Draw("SAMETEXT");
        cfp->SaveAs(outputpdf2);
    }
    cytar->Clear();
    h_rates->Draw("COLZ");
    h_rates->Draw("SAMETEXT");  
    HList.Add(h_rates);
    cytar->SaveAs(outputpdf3); 

    TText* t4 = new TText(0.5,0.5,Form("Lowest Inelastic MC rate: %f Hz",lowRate));
    t4->SetTextAlign(22);
    t4->SetY(2);

    TText* t5 = new TText(0.5, 0.6, Form("Found in xsnum = %i and ysnum = %i", lowX, lowY));
    t5->SetTextAlign(22);
    t5->SetY(3);

    TText* t6 = new TText(0.5,.4,Form("Low Time to 200 Events: %f s", lowTime));
    t6->SetTextAlign(22);
    t6->SetY(1);

    TText* t7 = new TText(0.5,0.2,Form("Summed rate: %f Hz", sumRate));
    t7->SetTextAlign(22);
    t7->SetY(-2);

    TText* t8 = new TText(0.1,0.3,Form("Low Time to 200 Events: %i hr %i min %2.2f s",hrs, min, sec ));
    t8->SetTextAlign(22);
    t8->SetY(0.5);

    TText* t9 = new TText(0.5,0.1,Form("Number of holes with 0 counts: %i holes", zeroCount));
    t9->SetTextAlign(22);
    t9->SetY(-2.5);

    TText* t10 = new TText(0.1,0.7,Form("Number of holes with 200 in 3hrs: %i holes", timeUnder3));
    t10->SetTextAlign(22);
    t10->SetY(4);

    //Print some useful values (i.e. the lowest rate, the number of 0-count holes, etc.) to the final page of the pdfs.
    if(check > 0){
        c->Clear();
        t6->Draw();
        t4->Draw("SAME");
        t5->Draw("SAME");
        t8->Draw("SAME");
        t7->Draw("SAME");
        t9->Draw("SAME");
        t10->Draw("SAME");
        c->SaveAs(outputpdf1+")");

        cfp->Clear();
        t6->Draw();
        t4->Draw("SAME");
        t5->Draw("SAME");
        t8->Draw("SAME");
        t7->Draw("SAME");
        t9->SetY(0.1);
        t9->Draw("SAME");
        t10->Draw("SAME");
        cfp->SaveAs(outputpdf2+")");
    } 

    cytar->Clear();
    t6->Draw();
    t4->Draw("SAME");
    t5->Draw("SAME");
    t8->Draw("SAME");
    t7->Draw("SAME");
    t9->SetY(0.1);
    t9->Draw("SAME");
    t10->Draw("SAME");
    cytar->SaveAs(outputpdf3+")");


    delete c;
    delete cfp;
    delete cytar;

    //Store all histograms to a root file.
    TFile hsimc(outputhist,"recreate");
    HList.Write();
    cout << " Plotted histograms put in root file = " << outputhist << endl;
    
    //Delete pointers.
    delete h_rates;
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