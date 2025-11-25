#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <iostream>

using namespace std;

void test2() {

  TFile *f = TFile::Open("~/compphys/ode3/outputs/vterm2.root");

  TGraph *g = (TGraph*) f->Get("terminal_velocity_vs_mass");
  //check to make sure it isn't empty 
  cout << "Graph loaded. Number of points = " << g->GetN() << endl;

  TCanvas *c = new TCanvas("c", "Terminal Velocity vs Mass", 900, 600);

  g->SetLineWidth(2);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.0);
  g->Draw("ALP");

  //Analytic curve vt = sqrt(mg / k) 
  double gval = 9.81;
  double k = 0.1;

  TF1 *analytic = new TF1("analytic", Form("sqrt(x*%f/%f)", gval, k), 1, 10);

  analytic->SetLineColor(kRed);
  analytic->SetLineStyle(2);
  analytic->SetLineWidth(2);
  analytic->Draw("same");

  TLegend *leg = new TLegend(0.15, 0.7, 0.45, 0.88);
  leg->AddEntry(g, "RK4 Simulation", "lp");
  leg->AddEntry(analytic, "Analytic: v_{t} = #sqrt{mg/k}", "l");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  c->SaveAs("terminal_velocity_vs_mass.png");
}
