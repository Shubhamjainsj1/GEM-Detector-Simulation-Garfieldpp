#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <fstream>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/TrackHeed.hh"

using namespace Garfield;

int main(int argc, char* argv[]) {
  TApplication app("app", &argc, argv);

  // ------------------ GAS ------------------
  MediumMagboltz gas("ar", 80., "co2", 20.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.Initialise(true);
  gas.EnablePenningTransfer(0.51, 0., "ar");
  gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");

  // ------------------ FIELD ------------------
  ComponentAnsys123 fm;
  fm.Initialise("ELIST.lis", "NLIST.lis", "MPLIST.lis", "PRNSOL.lis", "mm");
  fm.EnableMirrorPeriodicityX();
  fm.EnableMirrorPeriodicityY();
  fm.SetGas(&gas);

  constexpr double pitch = 0.014;

  // ------------------ VISUALIZATION ------------------
  ViewField fieldView(&fm);
  ViewFEMesh meshView(&fm);

  TCanvas* cf = new TCanvas("cf", "", 600, 600);
  fieldView.SetPlane(0, -1, 0, 0, 0, 0);
  fieldView.SetArea(-0.5 * pitch, -0.02, 0.5 * pitch, 0.02);
  fieldView.SetCanvas(cf);
  fieldView.PlotContour();

  meshView.SetArea(-0.5 * pitch, -0.02, 0.5 * pitch, 0.02);
  meshView.SetCanvas(cf);
  meshView.SetPlane(0, -1, 0, 0, 0, 0);
  meshView.SetFillMesh(true);
  meshView.Plot(true);

  // ------------------ SENSOR ------------------
  Sensor sensor(&fm);
  sensor.SetArea(-5 * pitch, -5 * pitch, -0.02,
                  5 * pitch,  5 * pitch,  0.03);

  // ------------------ TRACK ------------------
  TrackHeed track;
  track.SetParticle("mu-");
  track.SetMomentum(3.e9);
  track.SetSensor(&sensor);
  track.Initialise(&gas);

  // ------------------ AVALANCHE ------------------
  AvalancheMicroscopic aval(&sensor);

  // ------------------ DRIFT VIEW ------------------
  ViewDrift driftView;

  // ------------------ HISTOGRAM ------------------
  TH1F* hGain = new TH1F("gain", "Gain Distribution", 100, 0, 500);

  // ------------------ DATA ------------------
  std::vector<double> gainPerEvent;
  std::ofstream fout("gain.txt");

  // ------------------ SIMULATION ------------------
  int nEvents = 50;
  double totalGain = 0;
  double totalGain2 = 0;

  for (int i = 0; i < nEvents; ++i) {

    // Plot only first event (performance fix)
    if (i == 0) {
      aval.EnablePlotting(&driftView, 5);
    }

    // Randomize track position (optional but better)
    double x0 = (rand() / (double)RAND_MAX - 0.5) * pitch;
    double y0 = (rand() / (double)RAND_MAX - 0.5) * pitch;

    track.NewTrack(x0, y0, 0.02, 0., 0., -1., 0.);

    double xe, ye, ze, te;
    int ne, ni, np;
    double ec, extra;

    double eventGain = 0;

    while (track.GetCluster(xe, ye, ze, te, ne, ni, np, ec, extra)) {

      for (int j = 0; j < ne; ++j) {

        aval.AvalancheElectron(xe, ye, ze, te, 0., 0., 0.);

        // ✅ CORRECT WORKING GAIN
        int nEndpoints = aval.GetNumberOfElectronEndpoints();
        eventGain += nEndpoints;
      }
    }

    // Store + print
    gainPerEvent.push_back(eventGain);
    fout << i << " " << eventGain << std::endl;

    std::cout << "Event " << i << " Gain = " << eventGain << std::endl;

    if (eventGain == 0) {
      std::cout << "Zero gain event detected!" << std::endl;
    }

    hGain->Fill(eventGain);

    totalGain += eventGain;
    totalGain2 += eventGain * eventGain;
  }

  fout.close();

  // ------------------ RESULTS ------------------
  double avgGain = totalGain / nEvents;
  double variance = (totalGain2 / nEvents) - (avgGain * avgGain);
  double error = sqrt(variance);

  std::cout << "===================================" << std::endl;
  std::cout << "Average Gain = " << avgGain << std::endl;
  std::cout << "Error (RMS) = " << error << std::endl;
  std::cout << "===================================" << std::endl;

  // ------------------ HISTOGRAM ------------------
  TCanvas* cGain = new TCanvas("cGain", "Gain Histogram", 600, 600);
  hGain->GetXaxis()->SetTitle("Gain");
  hGain->GetYaxis()->SetTitle("Counts");
  hGain->Draw();

  // ------------------ DRIFT VISUALIZATION ------------------
  TCanvas* cd = new TCanvas();
  meshView.SetCanvas(cd);
  meshView.SetViewDrift(&driftView);
  meshView.Plot(true, false);

  app.Run();
}