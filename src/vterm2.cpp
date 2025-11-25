#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

struct Params {
  double g;
  double m;
  double air_k;
};

//ODEs

double f_ri(double x, const vector<double> &y, void *params=0){ 
  (void) x;
  return y[1];
}

double f_vi(double x, const vector<double> &y, void *params=0){ 
  (void) x;
  Params *p = (Params*)params;
  return -p->air_k * sqrt(y[1]*y[1] + y[3]*y[3]) * y[1] / p->m;
}

double f_rj(double x, const vector<double> &y, void *params=0){  
  (void) x;
  return y[3];
}

double f_vj(double x, const vector<double> &y, void *params=0){  
  (void) x;
  Params *p = (Params*)params;
  return -p->air_k * sqrt(y[1]*y[1] + y[3]*y[3]) * y[3] / p->m - p->g;
}


int main(int argc, char **argv){

  Params pars;
  pars.g = 9.81;
  pars.m = 1.0;
  pars.air_k = 0.1;

  void *p_par = (void*) &pars;

  TApplication theApp("App", &argc, argv);

  UInt_t dh = gClient->GetDisplayHeight()/2;
  UInt_t dw = 1.1*dh;

  vector<pfunc_t> v_fun(4);
  v_fun[0]=f_ri;
  v_fun[1]=f_vi;
  v_fun[2]=f_rj;
  v_fun[3]=f_vj;


  TGraph *gVterm = new TGraph();
  gVterm->SetName("terminal_velocity_vs_mass");
  gVterm->SetTitle("Terminal Velocity vs Mass; Mass [kg]; Terminal Velocity |v| [m/s]");

  int point = 0;

  cout << "\nMass [kg]    V_terminal [m/s]\n";
  cout << "-----------------------------\n";

  for (double mass = 1.0; mass <= 10.0; mass += 1.0) {

    pars.m = mass;

    vector<double> y(4);

    // Start high and let it fall straight down
    y[0] = 0.0;    // x position irrelevant
    y[1] = 0.0;    // vx = 0
    y[2] = 100.0;  // start high
    y[3] = 0.0;    // vy = 0

    double t = 0.0;
    double tmax = 50.0;
    int nsteps = 2000;
    double h = (tmax - t) / nsteps;

    double prev_vy = 0;
    double vterm = 0;

    // Manual stepping to detect terminal velocity
    for (int i = 0; i < nsteps; i++) {

      y = RK4StepN(v_fun, y, t, h, p_par);
      t += h;

      double vy = y[3];

      // If velocity stops changing significantly → terminal velocity
      if (fabs(vy - prev_vy) < 1e-4 && i > 100) {
        vterm = fabs(vy);
        break;
      }

      prev_vy = vy;
    }

    cout << mass << "          " << vterm << endl;

    gVterm->SetPoint(point++, mass, vterm);
  }

  TCanvas *c = new TCanvas("c", "Terminal Velocity vs Mass", dw, dh);

  gVterm->SetLineWidth(2);
  gVterm->SetMarkerStyle(20);
  gVterm->Draw("ALP");

  //so yet again, saved to a ROOT file in the spirit of RKnDemo.cpp but I needed another script to take a look
  //look at test2.cpp to check out this plot
  TFile *tf = new TFile("vterm2.root", "RECREATE");
  gVterm->Write();
  tf->Close();

  theApp.SetIdleTimer(15, ".q");
  theApp.Run();

  return 0;
}
