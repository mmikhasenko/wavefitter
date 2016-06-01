

cd Afull(cd s, vector<cd> apar, vector<cd> npar, 
	 double s1, double s2) {

  //Constract D
  cd dsum =0;
  if(s1!=gs1) { gs1=s1; build_integrals_table(*garr, s1); }
  for(int i=0;i<npar.size();i++) dsum += npar[i]*vtable(real(s),(*garr)[i]);
  cd Denomin = 1. - s / (2*M_PI) * dsum;

  //Constract N, Alpha and A.
  cd Numir = expand(s,s1,npar);
  cd Aprod = expand(s,s2,apar);
  cd Afull = Aprod*Numir/Denomin;

  return Afull;
}

void build_integrals_table(vector<vector<pair<double,cd> > > &data, double s1) {
  data.clear();
  for(int i=0;i<NPAR_FIT;i++) {
    vector<pair<double,cd> > vj;
    double step = 1.0*(VFUNC_RB - VFUNC_LB)/(VFUNC_NPOINT-1);
    for(int j=0;j<VFUNC_NPOINT;j++) {
      double s = VFUNC_LB + step*j;
      cd value = vf(s,i,s1,pow(RHO_MASS-PI_MASS,2),pow(RHO_MASS+PI_MASS,2));
      vj.push_back(make_pair<double,cd>(double(s),cd(value)));
    }
    data.push_back(vj);
  }
}


//////////////////////////////////////////////////////////////////////////////                  
///  //   /  ///      //      ///     ///      ///      ///  ///////     /////
///  //      /////  ////  ///////  //////  //  ///  //  ///  ///////   ///////
///  //  /   /////  ////    /////  /  ///    /////      ///  /////////   /////
///  //  //  /////  ////      ///     ///  /   ///  //  ///      ///     /////
//////////////////////////////////////////////////////////////////////////////

cd vint_sub(cd sp, cd s, double k, double s0, double sl, double sr) {
  cd rho_sp2 = (sp-sl)*(sp-sr)/(sp*sp);
  cd rho_s2  = (s -sl)*(s -sr)/(s*s);
  return (pow(w(sp,s0),k)*sqrt(rho_sp2)-pow(w(s,s0),k)*sqrt(rho_s2))/(sp*(sp-s));
}

double re_vint_sub(double *x, double *par) {
  cd spr(1/x[0],0.), sr(par[0],par[1]);
  return real(vint_sub(spr,sr,par[2],par[3],par[4],par[5]))/(x[0]*x[0]);
}
double im_vint_sub(double *x, double *par) {
  cd spr(1/x[0],0.), sr(par[0],par[1]);
  return imag(vint_sub(spr,sr,par[2],par[3],par[4],par[5]))/(x[0]*x[0]);
}

cd integral(cd s, double k, double s0, double sl, double sr) {

  TF1 fre("fre", re_vint_sub, 0, 1./sr,6);
  TF1 fim("fim", im_vint_sub, 0, 1./sr,6);
  double pars[] = {real(s),imag(s),k,s0,sl,sr};
  fre.SetParameters(pars);
  fim.SetParameters(pars);

  ROOT::Math::WrappedTF1 wre(fre);  ROOT::Math::GaussIntegrator igre;
  ROOT::Math::WrappedTF1 wim(fim);  ROOT::Math::GaussIntegrator igim;
  
  igre.SetFunction(wre); igre.SetRelTolerance(0.01);
  igim.SetFunction(wim); igim.SetRelTolerance(0.01);

  cd res(igre.Integral(0, 1./sr),igim.Integral(0, 1./sr));
  return res;
}

cd vf(cd s, double k, double s0, double sl, double sr) {
  //cout << "integral" << endl;
  cd rho_s2  = (s -sl)*(s -sr)/(s*s);
  cd mult = pow(w(s,s0),k)*sqrt(rho_s2)/s;
  cd unit(0,1);
  cd alog = 1.-(s+unit*EPSILON)/sr;
  cd first_int = integral(s,k,s0,sl,sr);
  return first_int-mult*log(alog);
}
