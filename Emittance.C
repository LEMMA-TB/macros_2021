
// usage:
// tested with ROOT/6.12.06
// root [0] .L Emittance.C++
// root [1] Emittance()

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <iomanip>

#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TAxis.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TRandom3.h"
#include "TPaveText.h"

using namespace std;

// function that returns emittance
Double_t getemittance(vector<Double_t> xv, vector<Double_t> xpv){

  UInt_t centering=1; // default: 1, only other reasonable value: 0

  UInt_t n_points = xv.size();

  Double_t emittance=0.;

  Double_t x=0.;
  Double_t x2=0.;
  Double_t xp=0.;
  Double_t xp2=0.;
  Double_t xxp=0.;

  for(UInt_t i=0;i<xv.size();i++){

    x+=(1000000*xv[i]);                // [nm]
    x2+=(1000000*xv[i]*1000000*xv[i]); // [nm x nm]
    xp+=(xpv[i]);                      // [rad]
    xp2+=(xpv[i]*xpv[i]);              // [rad x rad]
    xxp+=(1000000*xv[i]*xpv[i]);       // [nm x rad]

  }

  x*=1./float(n_points);
  x2*=1./float(n_points);
  xp*=1./float(n_points);
  xp2*=1./float(n_points);
  xxp*=1./float(n_points);

  if( centering==0 ){
    emittance=TMath::Sqrt(x2*xp2-xxp*xxp); // [nm x rad]
  }else{
    emittance=((x2-x*x)*(xp2-xp*xp)-(xxp-x*xp)*(xxp-x*xp));
    emittance=TMath::Sqrt(emittance);  
  }

  return emittance;

}

//Extrapolate track function: track points are refitted
Double_t extrapolate_track_x(TString inputFileName, Double_t z0, Double_t x_pos_mum[12], Double_t x_pos_mum_err[12], Double_t z_x_pos_mum[12], Double_t& x_ext, Double_t& x_ext_err, Double_t& dx_on_dz_ext, Double_t& dx_on_dz_ext_err){

  //Magnetic field (box)
  Double_t zM=0.,B=0.;
  if ( inputFileName.Contains("aug18") ){
    zM=17193-846;
    B=1.7476;
  }
  if( inputFileName.Contains("sep18") || inputFileName.Contains("Next") ){
    zM=10.*(1573.79+0.5*(1773.79-1573.79)-82.78);
    B=2.01;
    if( inputFileName.Contains("Next") ) B=2.0; // rounded value implemented in Geant4
  }
  Double_t z1=zM-1000.;
  Double_t z2=zM+1000.;
  if( B==0. ) cout << "B undefined in extrapolate_track_x !" << endl;

  Double_t chi2 = 99999;

  //3 cases: z > z_magnet_out=z2, z < z_magnet_entry=z1, z1<z<z2 

  if (z0>z2) {

    Int_t npoints=0;

    for (Int_t k=0; k<12; k++) {
 
      if (z_x_pos_mum[k]>z2) npoints++;

    }

    if (npoints < 2) return chi2;
   
    Double_t z_x_pos_mum_err[12];
    for (Int_t k=0; k<12; k++) {z_x_pos_mum_err[k]=0;}

    //Change TGraph with TGraphErrors to use errors
    //TGraphErrors* graph = new TGraphErrors(npoints, z_x_pos_mum, x_pos_mum, z_x_pos_mum_err, x_pos_mum_err);

    TGraph* graph = new TGraph(npoints, z_x_pos_mum, x_pos_mum);

    //Linear fit using points after magnet

    TF1* line = new TF1("line", "[0]+[1]*(x-[2])");

    Double_t theta_min = 0.03;
    Double_t theta_max = 0.07;

    line->FixParameter(2,z0);


    graph->Fit(line,"Q");


    chi2 =  line->GetChisquare()/line->GetNDF();

    x_ext = line->GetParameter(0);
    x_ext_err = line->GetParError(0);

    dx_on_dz_ext = line->GetParameter(1);
    dx_on_dz_ext_err = line->GetParError(1);

    delete graph;
    delete line;

  } else if (z0 < z1)
    {

      Int_t npoints=0;

      for (Int_t k=0; k<12; k++) {
 
	if ((z_x_pos_mum[k]<-0.001 || z_x_pos_mum[k]>0.001)  && z_x_pos_mum[k] < 22700) npoints++;

      }    

      if (npoints < 2) return chi2;

      Double_t z_x_pos_mum_err[12];
      for (Int_t k=0; k<12; k++) {z_x_pos_mum_err[k]=0;}

      //Change TGraph with TGraphErrors to use errors
      //TGraphErrors* graph = new TGraphErrors(npoints, z_x_pos_mum, x_pos_mum, z_x_pos_mum_err, x_pos_mum_err);

      TGraph* graph = new TGraph(npoints, z_x_pos_mum, x_pos_mum);

      //Fit to the complete track: line + parabula + line, 3 free parameters

      TF1* trajectory = new TF1("trajectory", "(x<[3])*([0]+[1]*(x-[3]))+(x>[4])*([0]-[2]*([4]-[3])*([4]-[3])+([1]+2*[2]*([4]-[3]))*(x-[3]))+(x>[3])*(x<[4])*([0]+[1]*(x-[3])+[2]*(x-[3])*(x-[3]))");

      trajectory->FixParameter(3,z1);
      trajectory->FixParameter(4,z2);
  
      graph->Fit(trajectory,"Q");

      Double_t R = 1./(2.*trajectory->GetParameter(2));
      Double_t p = -B/(1e+9/TMath::C()) * R;
  
      chi2 = trajectory->GetChisquare()/trajectory->GetNDF();

      x_ext = trajectory->Eval(z0);
      x_ext_err = trajectory->GetParError(0);

      dx_on_dz_ext = trajectory->GetParameter(1);
      dx_on_dz_ext_err = trajectory->GetParError(1);


      delete graph;
      delete trajectory;

    } else 
    {
      cout<< "Extrapolation not defined" << endl;
    }

  return chi2;

}

void Emittance(TString inputFileName="reco.root",Double_t z_ref=3811.5){

  TH1::SetDefaultSumw2(true);

  cout << "processing: " << inputFileName << endl;
  
  Double_t chi2Si5MuM;
  Double_t x_pos_mum[12];
  Double_t x_pos_mum_err[12];
  Double_t z_x_pos_mum[12];
  Double_t x_pos_DT_mum[8];
  Double_t z_pos_DT_mum[8];
  Double_t p_mum;
  Double_t p_mup;
  Double_t chi2Si5MuP;
  Double_t x_pos_mup[12];
  Double_t x_pos_mup_err[12];
  Double_t z_x_pos_mup[12];
  Double_t x_pos_DT_mup[8];
  Double_t z_pos_DT_mup[8];
  Double_t vtx_x;
  Double_t vtx_z;
  Double_t vtx_chi2;
  Int_t    subdet[100];
  Int_t    itrack[100];
  Double_t xh[100];
  Double_t yh[100];
  Double_t zh[100];
  Int_t    nhits;
  Double_t Calo_EnDep[25];
  Int_t    event_type;
  Double_t gen_vtx_mum[7];  // used for MC only
  Double_t gen_vtx_mup[7];  // used for MC only
  Double_t gen_ref_mum[6];  // used for MC only
  Double_t gen_ref_mup[6];  // used for MC only
  Double_t gen_pos_mum[12]; // used for MC only
  Double_t gen_pos_mup[12]; // used for MC only

  TFile* inputFile = new TFile(inputFileName);
  TTree* inputTree = (TTree*)inputFile->Get("lemma");

  inputTree->SetBranchAddress("chi2Si5MuM",	&chi2Si5MuM);
  inputTree->SetBranchAddress("x_pos_mum",      &x_pos_mum[0]); 
  inputTree->SetBranchAddress("x_pos_mum_err",  &x_pos_mum_err[0]);
  inputTree->SetBranchAddress("z_x_pos_mum",    &z_x_pos_mum[0]);
  inputTree->SetBranchAddress("x_pos_DT_mum",   &x_pos_DT_mum[0]);
  inputTree->SetBranchAddress("z_pos_DT_mum",   &z_pos_DT_mum[0]);
  inputTree->SetBranchAddress("p_mum",          &p_mum);	     
  inputTree->SetBranchAddress("p_mup",          &p_mup);	     
  inputTree->SetBranchAddress("chi2Si5MuP",     &chi2Si5MuP);
  inputTree->SetBranchAddress("x_pos_mup",      &x_pos_mup[0]); 
  inputTree->SetBranchAddress("x_pos_mup_err",  &x_pos_mup_err[0]);
  inputTree->SetBranchAddress("z_x_pos_mup",    &z_x_pos_mup[0]);
  inputTree->SetBranchAddress("x_pos_DT_mup",   &x_pos_DT_mup[0]);
  inputTree->SetBranchAddress("z_pos_DT_mup",   &z_pos_DT_mup[0]);
  inputTree->SetBranchAddress("vtx_x",          &vtx_x);
  inputTree->SetBranchAddress("vtx_z",          &vtx_z);
  inputTree->SetBranchAddress("vtx_chi2",       &vtx_chi2);
  inputTree->SetBranchAddress("subdet",         &subdet[0]);   
  inputTree->SetBranchAddress("itrack",         &itrack[0]);   
  inputTree->SetBranchAddress("xh",             &xh[0]);	     
  inputTree->SetBranchAddress("yh",             &yh[0]);	     
  inputTree->SetBranchAddress("zh",             &zh[0]);	     
  inputTree->SetBranchAddress("nhits",          &nhits);	     
  inputTree->SetBranchAddress("Calo_EnDep",     &Calo_EnDep[0]);
  inputTree->SetBranchAddress("event_type",     &event_type);   
  inputTree->SetBranchAddress("gen_vtx_mum", &gen_vtx_mum[0]);
  inputTree->SetBranchAddress("gen_vtx_mup", &gen_vtx_mup[0]);
  inputTree->SetBranchAddress("gen_ref_mum", &gen_ref_mum[0]);
  inputTree->SetBranchAddress("gen_ref_mup", &gen_ref_mup[0]); 
  inputTree->SetBranchAddress("gen_pos_mum", &gen_pos_mum[0]);
  inputTree->SetBranchAddress("gen_pos_mup", &gen_pos_mup[0]);

  TH1F* phTest = new TH1F("hTest","extrapolated mu+ position",100,-30.,30.);
  TH1F* phDeltaTest = new TH1F("hDeltaTest"," ",100,-1.,1.);

  Double_t h_min_x_emitt=-0.3,h_max_x_emitt=0.3,h_min_xprime_emitt=-0.002,h_max_xprime_emitt=0.002;
  // gen level
  // TH2D emittance plots
  TH2D* hist2D_emittance_x_mup_MC = new TH2D("hist2D_emittance_x_mup_MC","hist2D_emittance_x_mup_MC",100, h_min_x_emitt, h_max_x_emitt, 100, h_min_xprime_emitt, h_max_xprime_emitt);
  TH2D* hist2D_emittance_x_mum_MC = new TH2D("hist2D_emittance_x_mum_MC","hist2D_emittance_x_mum_MC",100, h_min_x_emitt, h_max_x_emitt, 100, h_min_xprime_emitt, h_max_xprime_emitt);
  // TH1D emittance plots
  TH1D* hist1D_emittance_x_mup_MC       = new TH1D("hist1D_emittance_x_mup_MC",      "hist1D_emittance_x_mup_MC",      100, h_min_x_emitt,      h_max_x_emitt);
  TH1D* hist1D_emittance_x_prime_mup_MC = new TH1D("hist1D_emittance_x_prime_mup_MC","hist1D_emittance_x_prime_mup_MC",100, h_min_xprime_emitt, h_max_xprime_emitt);
  TH1D* hist1D_emittance_x_mum_MC       = new TH1D("hist1D_emittance_x_mum_MC",      "hist1D_emittance_x_mum_MC",      100, h_min_x_emitt,      h_max_x_emitt);
  TH1D* hist1D_emittance_x_prime_mum_MC = new TH1D("hist1D_emittance_x_prime_mum_MC","hist1D_emittance_x_prime_mum_MC",100, h_min_xprime_emitt, h_max_xprime_emitt);
  // rec level
  // TH2D emittance plots
  TH2D* hist2D_emittance_x_mup_recMC = new TH2D("hist2D_emittance_x_mup_recMC","hist2D_emittance_x_mup_recMC",100, h_min_x_emitt, h_max_x_emitt, 100, h_min_xprime_emitt, h_max_xprime_emitt);
  TH2D* hist2D_emittance_x_mum_recMC = new TH2D("hist2D_emittance_x_mum_recMC","hist2D_emittance_x_mum_recMC",100, h_min_x_emitt, h_max_x_emitt, 100, h_min_xprime_emitt, h_max_xprime_emitt);
  // TH1D emittance plots
  TH1D* hist1D_emittance_x_mup_recMC       = new TH1D("hist1D_emittance_x_mup_recMC",      "hist1D_emittance_x_mup_recMC",      100, h_min_x_emitt,      h_max_x_emitt);
  TH1D* hist1D_emittance_x_prime_mup_recMC = new TH1D("hist1D_emittance_x_prime_mup_recMC","hist1D_emittance_x_prime_mup_recMC",100, h_min_xprime_emitt, h_max_xprime_emitt);
  TH1D* hist1D_emittance_x_mum_recMC       = new TH1D("hist1D_emittance_x_mum_recMC",      "hist1D_emittance_x_mum_recMC",      100, h_min_x_emitt,      h_max_x_emitt);
  TH1D* hist1D_emittance_x_prime_mum_recMC = new TH1D("hist1D_emittance_x_prime_mum_recMC","hist1D_emittance_x_prime_mum_recMC",100, h_min_xprime_emitt, h_max_xprime_emitt);

  vector<double> vec_emittance_x_mup_MC;
  vector<double> vec_emittance_xprime_mup_MC;
  vector<double> vec_emittance_x_mum_MC;
  vector<double> vec_emittance_xprime_mum_MC;
  //
  vector<double> vec_emittance_x_mup_recMC;
  vector<double> vec_emittance_xprime_mup_recMC;
  vector<double> vec_emittance_x_mum_recMC;
  vector<double> vec_emittance_xprime_mum_recMC;

  // --- loop over tree entries 
  Long64_t entries = inputTree->GetEntries();
  for(Long64_t z=0; z<entries; z++){

    inputTree->GetEntry(z);

    // --- condition for candidate events
    if( p_mup > 0. && p_mum > 0. ) {

      if( chi2Si5MuP<0. || chi2Si5MuP>100. ) continue;
      if( chi2Si5MuM<0. || chi2Si5MuM>100. ) continue;
      if( vtx_chi2>100. ) continue;

      Double_t x_atZref_eplus=0.;
      Double_t x_prime_atZref_eplus=0.;
      // --- incoming e+
      // px of e+ = Cx mu- * En of mu- + Cx mu+ * En of mu+ = px of mu- + px of mu+
      Double_t px_eplus = gen_vtx_mum[3]*gen_vtx_mum[6] + gen_vtx_mup[3]*gen_vtx_mup[6];
      Double_t py_eplus = gen_vtx_mum[4]*gen_vtx_mum[6] + gen_vtx_mup[4]*gen_vtx_mup[6];
      Double_t pz_eplus = gen_vtx_mum[5]*gen_vtx_mum[6] + gen_vtx_mup[5]*gen_vtx_mup[6];
      // x  of e+ (extrapolation on reference plane)  = x_vtx - (z_vtx - z_ref)*(px_e+/pz_e+)
      // N.B. gen_vtx_mup[0] = gen_vtx_mum[0] and gen_vtx_mup[2] = gen_vtx_mum[2] 
      x_atZref_eplus = gen_vtx_mup[0] - (gen_vtx_mup[2] - z_ref)*(px_eplus/pz_eplus);
      // x' of e+ (extrapolation on reference plane)  = px / p  (the direction remain the same as at vtx)
      x_prime_atZref_eplus = px_eplus / sqrt(px_eplus*px_eplus + py_eplus*py_eplus + pz_eplus*pz_eplus);

      Double_t x_atZref_mup=gen_ref_mup[0];
      Double_t x_prime_atZref_mup=(gen_ref_mup[3]/gen_ref_mup[5]);

      Double_t x_emittance_mup = x_atZref_mup - x_atZref_eplus;
      Double_t x_prime_emittance_mup = x_prime_atZref_mup - x_prime_atZref_eplus;
      //
      hist2D_emittance_x_mup_MC->Fill(x_emittance_mup, x_prime_emittance_mup);
      hist1D_emittance_x_mup_MC->Fill(x_emittance_mup);
      hist1D_emittance_x_prime_mup_MC->Fill(x_prime_emittance_mup);
      //
      vec_emittance_x_mup_MC     .push_back(x_emittance_mup);
      vec_emittance_xprime_mup_MC.push_back(x_prime_emittance_mup);

      // --- use extrapolate_track_x
      // z_ref=vtx_z; // consistency test for phDeltaTest
      Double_t chi2_mup=9999;
      Double_t x_ext_mup=-9999,x_ext_err_mup=-9999;
      Double_t dx_on_dz_ext_mup=-9999,dx_on_dz_ext_err_mup=-9999;
      if( false ){
        cout << "vtx: " << vtx_x << " " << vtx_z << endl;
        cout << "hits:" << endl;
        for(Int_t i=0;i<12;i++){
            cout << i << " " << x_pos_mup[i] << " " << x_pos_mup_err[i] << " " << z_x_pos_mup[i] << endl;
	}
      }
      chi2_mup = extrapolate_track_x("Next",z_ref, x_pos_mup, x_pos_mup_err, z_x_pos_mup, 
		  		     x_ext_mup, x_ext_err_mup, dx_on_dz_ext_mup, dx_on_dz_ext_err_mup);
      phTest->Fill(x_ext_mup);
      phDeltaTest->Fill(x_ext_mup-vtx_x);

      x_atZref_mup=x_ext_mup;
      x_prime_atZref_mup=dx_on_dz_ext_mup;

      x_emittance_mup = x_atZref_mup - x_atZref_eplus;
      x_prime_emittance_mup = x_prime_atZref_mup - x_prime_atZref_eplus;
      //
      hist2D_emittance_x_mup_recMC->Fill(x_emittance_mup, x_prime_emittance_mup);
      hist1D_emittance_x_mup_recMC->Fill(x_emittance_mup);
      hist1D_emittance_x_prime_mup_recMC->Fill(x_prime_emittance_mup);
      //
      vec_emittance_x_mup_recMC     .push_back(x_emittance_mup);
      vec_emittance_xprime_mup_recMC.push_back(x_prime_emittance_mup);


    } // end if (p_mup > 0. && p_mum > 0.)

  } // end loop over tree entries 

  Double_t emittance_mup_MC = getemittance(vec_emittance_x_mup_MC, vec_emittance_xprime_mup_MC);
  cout << "emittance_mup_MC = " << emittance_mup_MC << endl;
  Double_t emittance_mup_recMC = getemittance(vec_emittance_x_mup_recMC, vec_emittance_xprime_mup_recMC);
  cout << "emittance_mup_recMC = " << emittance_mup_recMC << endl;

  // --- plot histos

  gStyle->SetOptStat(1);

  TCanvas *c_1DEmit_mup_MC = new TCanvas("c_1DEmit_mup_MC"); c_1DEmit_mup_MC->Divide(1,2);
  c_1DEmit_mup_MC->cd(1); hist1D_emittance_x_mup_MC->Draw("box");
  c_1DEmit_mup_MC->cd(2); hist1D_emittance_x_prime_mup_MC->Draw("box");
  TCanvas *c_2DEmit_mup_MC = new TCanvas("c_2DEmit_mup_MC"); c_2DEmit_mup_MC->cd();
  hist2D_emittance_x_mup_MC->Draw("box");

  TCanvas *c_1DEmit_mup_recMC = new TCanvas("c_1DEmit_mup_recMC"); c_1DEmit_mup_recMC->Divide(1,2);
  c_1DEmit_mup_recMC->cd(1); hist1D_emittance_x_mup_recMC->Draw("box");
  c_1DEmit_mup_recMC->cd(2); hist1D_emittance_x_prime_mup_recMC->Draw("box");
  TCanvas *c_2DEmit_mup_recMC = new TCanvas("c_2DEmit_mup_recMC"); c_2DEmit_mup_recMC->cd();
  hist2D_emittance_x_mup_recMC->Draw("box");

  TCanvas *c_Test = new TCanvas("c_Test"); c_Test->Divide(1,2);
  c_Test->cd(1); phTest->Draw();
  c_Test->cd(2); phDeltaTest->Draw();

  cout<< " Plots done ! " << endl; 

}
