#include "MStructureHolder.h"
#include "TRandom.h"
#include "TString.h"

MStructureHolder::MStructureHolder() :
  can_calculate_together(true) {
}

void MStructureHolder::AddWave(TH1D& intensity, cd (*amp)(int, double, const double*), double Npar,
			       std::vector<std::pair<double,double> > &phasespace) {

  _intensity.push_back(&intensity);
  _ointensity.push_back(0);
  _amplitude.push_back(amp);
  _phasespace.push_back(&phasespace);
  _justToplot_wave.push_back(false);
  std::vector<double>_wpar(Npar); pars.push_back(_wpar);  
  //if(N1!=N2||N2!=N3) std::cerr<<"Warning<MStr.Hold.::AddWave> : '!='"<<srd::endl;
  
  const int Nbins = intensity.GetXaxis()->GetNbins();
  _waverange.push_back(std::make_pair<double,double>(intensity.GetXaxis()->GetBinLowEdge(1),
						    intensity.GetXaxis()->GetBinUpEdge (Nbins)
						    ));
  //check if it has the same amount of bins and the ranges are the same
  if(
     Nbins != _intensity[0]->GetXaxis()->GetNbins() ||
     intensity.GetXaxis()->GetBinLowEdge(1) != _intensity[0]->GetXaxis()->GetBinLowEdge(1) || 
     intensity.GetXaxis()->GetBinUpEdge(Nbins) != _intensity[0]->GetXaxis()->GetBinUpEdge(Nbins) 
     ) {
    std::cout << "can_calculate_together " << intensity.GetName() << std::endl;
    can_calculate_together=false;
  }

}

void MStructureHolder::SetWaveModel(int w, cd (*amp)(int, double, const double*), double Npar) {
  if(w<0||w>=pars.size()) {std::cerr<<"Error<SetWaveModel> : w range"<<std::endl; return;}
  _amplitude[w] = amp;
  std::vector<double>_wpar(Npar); pars[w] = _wpar;
}


void MStructureHolder::AddInterference(TH1D& interference, int i, int j,
				       double (*calcinterf)(cd,cd)) {
  _interference.push_back(std::make_pair<TH1D*, std::pair<int,int> >
			  (&interference, std::make_pair<int,int>(int(i),int(j))));
  _ointerference.push_back(0);
  _calcinterf.push_back(calcinterf);
  _justToplot_interf.push_back(false);
  //if(N1!=N2) std::cerr<<"Warning<MStr.Hold.::AddInterf.> : '!='"<<srd::endl;

  const int Nbins = interference.GetXaxis()->GetNbins();
  _interfrange.push_back(std::make_pair<double,double>(interference.GetXaxis()->GetBinLowEdge(1),
						       interference.GetXaxis()->GetBinUpEdge (Nbins)
						       ));
  //check if it has the same amount of bins and the ranges are the same
  if(
     Nbins != _interference[0].first->GetXaxis()->GetNbins() ||
     interference.GetXaxis()->GetBinLowEdge(   1) != _interference[0].first->GetXaxis()->GetBinLowEdge(1) || 
     interference.GetXaxis()->GetBinUpEdge(Nbins) != _interference[0].first->GetXaxis()->GetBinUpEdge(Nbins)
     ) {
    std::cout << "can_calculate_together " << interference.GetName() << std::endl;
    can_calculate_together=false;
  }
}


double MStructureHolder::GetWaveChi2(int w, double mleft, double mright) const {
  double _mleft  = (mleft==mright) ? _waverange[w].first  : mleft;
  double _mright = (mleft==mright) ? _waverange[w].second : mright;

  if(w<0||w>=_intensity.size()) {std::cerr<<"Error<GetWaveChi2> : w range"<<std::endl; return 0;}
  TH1D &h = *_intensity[w];
  const int Nbins = h.GetXaxis()->GetNbins();
  double chi2 = 0;
  for(int i=1;i<=Nbins;i++) {
    double M = h.GetBinCenter(i);
    if(M<_mleft||M>_mright) continue;
    double s = M*M;
    cd amp = _amplitude[w](0,M*M,pars[w].data());
    double ph = getphsp(s,*_phasespace[w]);
    double intens = norm(amp)*ph;
    double err = h.GetBinError(i);
//    std::cout << "M = " << M << ", A = " << amp 
//	      << ", C = " << std::setprecision(10) << h.GetBinContent(i) 
//	      << ", E = " << std::setprecision(10) << h.GetBinError(i)
//	      << ", H = " << std::setprecision(10) << ph
//	      << ", P = " << std::setprecision(10) << pow(intens - h.GetBinContent(i),2)/(err*err) << std::endl;
    if(err==0) {std::cerr<<"Error<GetWaveChi2> : err == 0 "<<std::endl; err = 1;}
    chi2 += pow(intens - h.GetBinContent(i),2)/(err*err);
  }
  return chi2;
}

double MStructureHolder::GetInterferenceChi2(int l, double mleft, double mright) const {
  double _mleft  = (mleft>=mright) ? _interfrange[l].first  : mleft;
  double _mright = (mleft>=mright) ? _interfrange[l].second : mright;

  int w1 = _interference[l].second.first;
  int w2 = _interference[l].second.second;
  TH1D &h = *_interference[l].first;
  const int Nbins = h.GetXaxis()->GetNbins();
  double chi2 = 0;
  for(int i=1;i<=Nbins;i++) {
    double M = h.GetBinCenter(i);
    if(M<_mleft||M>_mright) continue;
    double s = M*M;
    cd amp1 = _amplitude[w1](0,M*M,pars[w1].data());
    cd amp2 = _amplitude[w2](0,M*M,pars[w2].data());
    double interf = _calcinterf[l](amp1,amp2);
    double err = h.GetBinError(i);
    if(err==0) {std::cerr<<"Error<GetInterferenceChi2> : err == 0 "<<std::endl; err = 1;}
    double dphi = interf - h.GetBinContent(i);
    if(dphi>M_PI) dphi-=2*M_PI; if(dphi<-M_PI) dphi+=2*M_PI;
    double value = pow(dphi,2)/(err*err);
    chi2 += value;
  }
  return chi2;  
}

double MStructureHolder::getphsp(double s, std::vector<std::pair<double,double> > &table) {
  const int N = table.size();
  const double lft = table[0].first;
  const double rht = table[N-1].first;
  const double Mstep = table[1].first - lft;
  const int Nsteps = (sqrt(s) - lft)/Mstep;
  if(Nsteps<0 || Nsteps>=N-1) {std::cerr<<"Error!! in getvalue! sqrt(s) = "<<sqrt(s)<<std::endl; return 0;}
  const double value = table[Nsteps].second + 
    ( table[Nsteps+1].second - table[Nsteps].second ) / 
    ( table[Nsteps+1].first  - table[Nsteps].first  ) * (sqrt(s) - table[Nsteps].first); 
  return value;
}

void MStructureHolder::SetWaveRange(int w, double mleft, double mright) {
  if(w<0||w>=_intensity.size()) {std::cerr<<"Error<SetWaveRange> : w range"<<std::endl; return;}
  _waverange[w].first = mleft;
  _waverange[w].second = mright;
}
void MStructureHolder::SetInterfRange(int w, double mleft, double mright) {
  if(w<0||w>=_interference.size()) {std::cerr<<"Error<SetWaveRange> : w range"<<std::endl; return;}
  _interfrange[w].first = mleft;
  _interfrange[w].second = mright;  
}

double MStructureHolder::GetChi2() const {
  if(!can_calculate_together) {
    std::cerr << "------------ can_not_calculate_together!----------------" << std::endl;
    return 0;
  } else {
    double chi2 = 0;
    TH1D &h = *_intensity[0];
    const int Nbins = h.GetXaxis()->GetNbins();
    const int NInterf = _interference.size();
    const int Namps = _amplitude.size();
    for(int i=1;i<=Nbins;i++) {
      double M = h.GetBinCenter(i);
      double s = M*M;
      //first calculate waves
      cd amp[Namps]; 
      for(int w=0;w<Namps;w++) {
	// check if it is possible to avoid rest of calculations
	bool not_needed_for_interf = true;
	for(int wi=0;wi<NInterf;wi++) 
	  if(!_justToplot_interf[wi] && 
	     (w==_interference[wi].second.first || w==_interference[wi].second.second) &&
	     (M > _interfrange[wi].first || M<_interfrange[wi].second )) { 
	    not_needed_for_interf = false; break;}
	if(not_needed_for_interf && (_justToplot_wave[w] || M < _waverange[w].first || M > _waverange[w].second)) continue;
	// calculate amplitude
	amp[w] = _amplitude[w](0,s,pars[w].data());
	// check if it is possible to avoid rest of calculations
	if(_justToplot_wave[w]) continue;	
	if(M < _waverange[w].first || M > _waverange[w].second) continue;
	// calculate rest
	double intens  = norm(amp[w])*getphsp(s,*_phasespace[w]);
	double content = _intensity[w]->GetBinContent(i); 
	double err = _intensity[w]->GetBinError(i);
	if(err==0) {std::cerr<<"Error<GetChi2> : err == 0 "<<std::endl; err = 1;}
	//std::cout<<"intens: "<<M << ", " << intens << ", " << content << ", err = "<<err<<std::endl;
	chi2 += pow(intens - content,2)/(err*err);
      }
      for(int w=0;w<NInterf;w++) {
	//std::cout << chi2 << std::endl;
	if(_justToplot_interf[w]) {/*std::cout<<"skipping interf "<<w<<std::endl;*/ continue;}
	if(M < _interfrange[w].first || M > _interfrange[w].second) continue;
	const int w1 = _interference[w].second.first;
	const int w2 = _interference[w].second.second;
	cd amp1 = amp[ w1 ]; //if(amp1==SPECIAL_NUMBER) {std::cout << "Error: chech our code!" << std::endl; amp1=_amplitude[ w1 ]( 0,s,pars[ w1 ].data());}
	cd amp2 = amp[ w2 ]; //if(amp2==SPECIAL_NUMBER) {std::cout << "Error: chech our code!" << std::endl; amp2=_amplitude[ w2 ]( 0,s,pars[ w2 ].data());}
	double interf = _calcinterf[w](amp1,amp2);
	double content = _interference[w].first->GetBinContent(i);
	double err = _interference[w].first->GetBinError(i);
	//std::cout<<"interf: "<<M << ", " << interf << ", " << amp1 << ", " << amp2 << content <<std::endl;
	if(err==0) {std::cerr<<"Error<GetChi2> : err == 0 "<<std::endl; err = 1;}
	double diff = interf - content;
	if(diff>M_PI) diff-=2*M_PI; if(diff<-M_PI) diff+=2*M_PI;
	chi2 += diff*diff/(err*err);
      }
    }
    //std::cout << chi2 << std::endl;
    return chi2;
  }
}

void MStructureHolder::JustToPlot(int w, bool status) {
  if(w<0||w>=_intensity    .size()) {std::cerr<<"Error<JustToPlot> : w range"<<std::endl; return;}
  _justToplot_wave[w] = status;
}
void MStructureHolder::JustToPlot(int i, int w, bool status) {
  if(w<0||w>=_interference.size()) {std::cerr<<"Error<JustToPlot> : w range"<<std::endl; return;}
  _justToplot_interf[w] = status;
}

void MStructureHolder::reGenerateData() {
  for(int w=0;w<_intensity.size();w++) {
    if(_ointensity[w]==0) _ointensity[w] = (TH1D*)_intensity[w]->Clone(TString::Format("orig_%s",_intensity[w]->GetName())); 
    TH1D &h  = * _intensity[w];
    TH1D &h0 = *_ointensity[w];
    const int Nbins = h.GetXaxis()->GetNbins();
    for(int j=1;j<=Nbins;j++) {
      double content = h0.GetBinContent(j);
      double error = h0.GetBinError(j);
      if(error != 0.0) h.SetBinContent(j, gRandom->Gaus(content,error));
    }
  }
  std::cout << "Done with intensities" << std::endl;
  for(int w=0;w<_interference.size();w++) {
    if(_ointerference[w]==0) _ointerference[w] = (TH1D*)_interference[w].first->Clone(TString::Format("orig_%s",_interference[w].first->GetName())); 
    TH1D &h  = * _interference[w].first;
    TH1D &h0 = *_ointerference[w];
    const int Nbins = h.GetXaxis()->GetNbins();
    for(int j=1;j<=Nbins;j++) {
      double content = h0.GetBinContent(j);
      double error = h0.GetBinError(j);
      if(error != 0.0) h.SetBinContent(j, gRandom->Gaus(content,error));
    }
  }
  std::cout << "Done with interderences" << std::endl;
}


void MStructureHolder::Print() {
  std::cout << "Intensities: < ";
  for(int w=0;w<_intensity.size();w++) 
    std::cout << ((_justToplot_wave[w]) ? "(" : "") 
	      << _intensity[w]->GetName()
	      << ((_justToplot_wave[w]) ? ")" : "") << " ";
  std::cout << ">" << std::endl;
  std::cout << "Phases: < ";
  for(int w=0;w<_interference.size();w++) 
    std::cout << ((_justToplot_interf[w]) ? "(" : "") 
	      << _interference[w].first->GetName()
	      << ((_justToplot_interf[w]) ? ")" : "") << " ";
  std::cout << ">" << std::endl;
}

/*...........................................................................................*/
/*...........................................................................................*/
/*.....................................Plots In Style........................................*/
/*...........................................................................................*/
/*...........................................................................................*/

TGraphErrors *MStructureHolder::GetWavePlot(int w, int flag, double mleft, double mright, 
					    TGraphErrors *gr)  const {
 
  double _mleft  = (mleft>=mright) ? _waverange[w].first  : mleft;
  double _mright = (mleft>=mright) ? _waverange[w].second : mright;

  TGraphErrors *_gr = (gr!=0) ? gr : new TGraphErrors(1000);
  const int Npoints = _gr->GetN();

  for(int i=0;i<Npoints;i++) {
    double M = _mleft + i*(_mright-_mleft)/(Npoints-1);
    cd amp = _amplitude[w](flag,M*M,pars[w].data());
    double intens = norm(amp)*getphsp(M*M,*_phasespace[w]);
    _gr->GetX()[i] = M;
    _gr->GetY()[i] = intens;
  }
  return _gr;
}

TGraphErrors *MStructureHolder::GetPolarPlot(int w, int flag, double mleft, double mright, 
					     TGraphErrors *gr, int Npoints) {
  double _mleft  = (mleft>=mright) ? _waverange[w].first  : mleft;
  double _mright = (mleft>=mright) ? _waverange[w].second : mright;

  TGraphErrors *_gr = (gr!=0) ? gr : new TGraphErrors(Npoints);

  for(int i=0;i<Npoints;i++) {
    double M = _mleft + i*(_mright-_mleft)/(Npoints-1);
    cd amp = _amplitude[w](flag,M*M,pars[w].data());
    _gr->GetX()[i] = real(amp);
    _gr->GetY()[i] = imag(amp);
  }
  return _gr;
}

TGraphErrors *MStructureHolder::GetInterfPlot(int l, double mleft, double mright, 
					      TGraphErrors *gr) const {
 
  double _mleft  = (mleft==mright) ? _interfrange[l].first  : mleft;
  double _mright = (mleft==mright) ? _interfrange[l].second : mright;

  TGraphErrors *_gr = (gr!=0) ? gr : new TGraphErrors(1000);
  const int Npoints = _gr->GetN();

  int w1 = _interference[l].second.first;
  int w2 = _interference[l].second.second;
  for(int i=0;i<Npoints;i++) {
    double M = _mleft + i*(_mright-_mleft)/(Npoints-1);
    cd amp1 = _amplitude[w1](0,M*M,pars[w1].data());
    cd amp2 = _amplitude[w2](0,M*M,pars[w2].data());
    double interf = _calcinterf[l](amp1,amp2);
    _gr->GetX()[i] = M;
    _gr->GetY()[i] = interf;
  }
  return _gr;
}

TGraphErrors *MStructureHolder::GetInterfPlotWithError(int l, double mleft, double mright,
						       std::vector<double*> inpars, void (*setPar)(const double *) 
						       ) {

  double _mleft  = (mleft==mright) ? _interfrange[l].first  : mleft;
  double _mright = (mleft==mright) ? _interfrange[l].second : mright;
  
  TGraphErrors *_gr = new TGraphErrors(100);
  const int Npoints = _gr->GetN();

  int w1 = _interference[l].second.first;
  int w2 = _interference[l].second.second;
  for(int i=0;i<Npoints;i++) {
    double M = _mleft + i*(_mright-_mleft)/(Npoints-1);
    std::vector<double> value(inpars.size());
    for(int j=0;j<inpars.size();j++) {
      setPar(inpars[j]);
      cd amp1 = _amplitude[w1](0,M*M,pars[w1].data());
      cd amp2 = _amplitude[w2](0,M*M,pars[w2].data());
      double interf = _calcinterf[l](amp1,amp2);
      value[j] = interf;
    }
    double mean = 0; for(int j=0;j<inpars.size();j++) mean += value[j]/inpars.size(); 
    double sigma2 = 0; for(int j=0;j<inpars.size();j++) sigma2 += pow(value[j]-mean,2)/inpars.size(); 
    std::cout << mean << " +- "<< sqrt(sigma2) << std::endl;
    _gr->GetX()[i] = M;
    _gr->GetY()[i] = mean;
    _gr->GetEY()[i] = sqrt(sigma2);
  }
  return _gr;
}
