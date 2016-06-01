// Copyright [2016] Mikhail Mikhasenko

#ifndef __MSTRUCTURE_HOLDER_H__
#define __MSTRUCTURE_HOLDER_H__

#include <deflib.h>
#include <pair>
#include <iostream>
#include <vector>

#include "TGraphErrors.h"
#include "TH1D.h"

class MStructureHolder {
 public:
  MStructureHolder();

 public:
  void AddWave(const TH1D& intensity, cd (*amp)(int, double, const double*), double Npar,
               const std::vector<std::pair<double, double> > &phasespace);
  void AddInterference(const TH1D& interference, int, int, double (*interf)(cd, cd));

  void SetWaveModel(int i, cd (*amp)(int, double, const double*), double Npar);

  double GetWaveLowRange(int i) {return _waverange[i].first; }
  double GetWaveUpRange (int i) {return _waverange[i].second;}
  double GetInterfLowRange(int i) {return _interfrange[i].first; }
  double GetInterfUpRange (int i) {return _interfrange[i].second;}
  void SetWaveRange(int i, double mleft, double mright);
  void SetInterfRange(int i, double mleft, double mright);

  double GetChi2() const;
  double GetWaveChi2(int w, double mleft = 0., double mright = 0.) const;
  double GetInterferenceChi2(int l, double mleft = 0, double mright = -1.) const;

  void JustToPlot(int w, bool status = true);
  void JustToPlot(int i, int w, bool status = true);

  TGraphErrors *GetWavePlot(int w, int flag, double mleft = 0, double mright = 0, TGraphErrors *gr = 0) const;
  TGraphErrors *GetInterfPlot(int l, double mleft = 0, double mright = 0, TGraphErrors *gr = 0) const;
  TGraphErrors *GetInterfPlotWithError(int l, double mleft, double mright,
                                       std::vector<double*> pars, void (*setPar)(const double *) );
  TGraphErrors *GetPolarPlot(int w, int flag, double mleft, double mright, TGraphErrors *gr = 0, int Npoints = 1000);

  void reGenerateData();

 public:
  std::vector<std::vector<double> > pars;

 private:
  bool can_calculate_together;
  std::vector<bool> _justToplot_wave;
  std::vector<bool> _justToplot_interf;

 private:
  // intensity
  std::vector<TH1D*> _intensity;
  std::vector<TH1D*> _ointensity;
  std::vector<cd (*)(int, double s, const double* pars)> _amplitude;
  std::vector<std::pair<double, double> > _waverange;
  // interference
  std::vector<std::pair<TH1D*, std::pair<int, int> > > _interference;
  std::vector<TH1D*> _ointerference;
  std::vector<double (*)(cd, cd)> _calcinterf;
  std::vector<std::pair<double, double> > _interfrange;
  // phase space
  std::vector<std::vector<std::pair<double, double> > *> _phasespace;

 public:
  static double getphsp(double s, const std::vector<std::pair<double, double> > &phsp);

  void Print();

};

#endif
