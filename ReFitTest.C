
// usage:
// tested with ROOT/6.12.06
// root [0] .L ReFitTest.C++
// root [1] ReFitTest()

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

void ReFitTest(TString inputFileName="reco.root"){

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
  Double_t gen_pos_mum[12]; // used for MC only
  Double_t gen_pos_mup[12]; // used for MC only
  Double_t gen_vtx_mum[7];  // used for MC only
  Double_t gen_vtx_mup[7];  // used for MC only

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
  inputTree->SetBranchAddress("gen_pos_mum", &gen_pos_mum[0]);
  inputTree->SetBranchAddress("gen_pos_mup", &gen_pos_mup[0]); 
  inputTree->SetBranchAddress("gen_vtx_mum", &gen_vtx_mum[0]);
  inputTree->SetBranchAddress("gen_vtx_mup", &gen_vtx_mup[0]); 

  TH1F* phTest = new TH1F("hTest","extrapolated mu+ position",100,-30.,30.);
  TH1F* phDeltaTest = new TH1F("hDeltaTest"," ",100,-1.,1.);

  // --- loop over tree entries 
  Long64_t entries = inputTree->GetEntries();
  for(Long64_t z=0; z<entries; z++){

    inputTree->GetEntry(z);

    // --- condition for candidate events
    if( p_mup > 0. && p_mum > 0. ) {
    if( chi2Si5MuP>500. ) continue;
    if( chi2Si5MuM>500. ) continue;
    if( vtx_chi2==9999. ) continue;

      // --- use extrapolate_track_x
      // Double_t z_ref=10.*(441.63+16.3+6.-82.78); // [mm] target end position
      Double_t z_ref=vtx_z;
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

    } // end if (p_mup > 0. && p_mum > 0.)

  } // end loop over tree entries 

  // --- plot histos

  gStyle->SetOptStat(1);

  TCanvas *c_Test = new TCanvas("c_Test"); c_Test->Divide(1,2);
  c_Test->cd(1); phTest->Draw();
  c_Test->cd(2); phDeltaTest->Draw();

  cout<< " Plots done ! " << endl; 

}
