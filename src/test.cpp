#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>

using namespace std;

void test() {

  TFile *f = TFile::Open("~/compphys/ode3/outputs/vterm.root");

  const int Ngraphs = 4;
  const char* graphNames[Ngraphs] = {
    "energy_50steps",
    "energy_100steps",
    "energy_200steps",
    "energy_400steps"
  };

  TGraph* graphs[Ngraphs];

  for (int i = 0; i < Ngraphs; i++) {
    graphs[i] = (TGraph*) f->Get(graphNames[i]);

    //once again, I had some issues with JSROOT so checking the number of points is a good idea
    cout << "Loaded " << graphNames[i] 
         << " with " << graphs[i]->GetN() 
         << " points" << endl;
  }

  TCanvas *c = new TCanvas("c", "Energy vs Time - Step Size Comparison", 1000, 700);
  int colors[4] = {2, 4, 6, 8};  // red, blue, magenta, green


  graphs[0]->SetLineColor(colors[0]);
  graphs[0]->SetLineWidth(2);
  graphs[0]->SetTitle("Energy vs Time; t [s]; Energy [J]");
  graphs[0]->Draw("AL");

  for (int i = 1; i < Ngraphs; i++) {
    graphs[i]->SetLineColor(colors[i]);
    graphs[i]->SetLineWidth(2);
    graphs[i]->Draw("L SAME");
  }

  TLegend *leg = new TLegend(0.15, 0.7, 0.4, 0.88);
  leg->SetBorderSize(0);

  leg->AddEntry(graphs[0], "50 steps",  "l");
  leg->AddEntry(graphs[1], "100 steps", "l");
  leg->AddEntry(graphs[2], "200 steps", "l");
  leg->AddEntry(graphs[3], "400 steps", "l");

  leg->Draw();

  c->SaveAs("energy_step_comparison.png");
}
