//Misha Mikhasenko
//14.07.2015

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

#include "MIsobar.h"
#include "./deflib.h"
#define SCALEX_FOR_CUT 2.
#define SCALEY_FOR_CUT 10.0

MIsobar::MIsobar(double Mi, double G0i,
                 double m1i, double m2i, double m3i,
                 int Li, double Ri):
  M(Mi), G0(G0i), m1(m1i), m2(m2i), m3(m3i), L(Li), R(Ri) {;}

double MIsobar::rho(double s, double s12, double m3_2) const {
  return (sqrt(s) > sqrt(s12)+sqrt(m3_2)) ? 1./(8*M_PI)*sqrt( LAMBDA(s, s12, m3_2)/(s*s) ) : 0;
}
cd     MIsobar::rho(cd     s, cd     s12, double m3_2) const {
  return 1./(8*M_PI)*sqrt( LAMBDA(s, s12, m3_2)/(s*s) );
}
cd     MIsobar::rhoPi(cd     s, cd     s12, double m3_2) const {
  return 1./(8*M_PI)*sqrtPi( LAMBDA(s, s12, m3_2)/(s*s) );
}

double MIsobar::U(double s) const {
  double G = G0 / rho(M*M, m1*m1, m2*m2)  *  rho(s, m1*m1, m2*m2);
  if (L > 0) {
    const double p_2 = LAMBDA(s, m1*m1, m2*m2)/(4.*s);
    const double p0_2 = LAMBDA(M*M, m1*m1, m2*m2)/(4.*M*M);
    G *=  pow(p_2/p0_2 * (1.+R*R*p0_2)/(1.+R*R*p_2), L);
  }
  return 2.*M*G/(pow(M*M-s, 2)+pow(M*G, 2));
}
cd     MIsobar::U(cd     s) const {
  cd G = G0/rho(M*M, m1*m1, m2*m2)  *  rhoPi(s, m1*m1, m2*m2);
  if (L > 0) {
    const cd p_2 = LAMBDA(s, m1*m1, m2*m2)/(4.*s);
    const double p0_2 = LAMBDA(M*M, m1*m1, m2*m2)/(4.*M*M);
    G *=  pow(p_2/p0_2 * (1.+R*R*p0_2)/(1.+R*R*p_2), L);
  }
  return -2.0*imag( 1./(M*M-s+cd(0, M)*G) );
  // return 2.*imag( 1./(M*M-s-cd(0, M)*G) );
}

//////////////////////////////////////////////////////////////////////////////
///  //   /  ///      //      ///     ///      ///      ///  ///////     /////
///  //      /////  ////  ///////  //////  //  ///  //  ///  ///////   ///////
///  //  /   /////  ////    /////  /  ///    /////      ///  /////////   /////
///  //  //  /////  ////      ///     ///  /   ///  //  ///      ///     /////
//////////////////////////////////////////////////////////////////////////////

double MIsobar::dph(double *x, double *par) const  {
  double s = par[0];
  double th1 = pow(m1+m2, 2);
  double th2 = pow(sqrt(s)-m3, 2);
  double s12 = th1+x[0]*(th2-th1);
  double Ul = U(s12);
  double rhol = rho(s, s12, m3*m3);
  return (th2-th1)/(2.*M_PI)*Ul*rhol;
}

// double MIsobar::rdph(double *x, double *par) const {
//  cd s(par[0],par[1]);
//  cd th1 = pow(m1+m2,2);
//  cd th2 = pow(sqrt(s)-m3,2);
//  cd s12 = th1 + x[0]*(th2-th1);
//  cd Ul = U(s12);
//  cd rhol = rho(s,s12,m3*m3);
//  return real( (th2-th1) / (2.*M_PI) * Ul * rhol );
//}
//
//double MIsobar::idph(double *x, double *par) const {
//  cd s(par[0],par[1]);
//  cd th1 = pow(m1+m2,2);
//  cd th2 = pow(sqrt(s)-m3,2);
//  cd s12 = th1 + x[0]*(th2-th1);
//  cd Ul = U(s12);
//  cd rhol = rho(s,s12,m3*m3);
//  //std::cout << "s12 = "<<s12<<", U = "<<Ul<<", rhol = "<<rhol << std::endl;
//  return imag( (th2-th1) / (2.*M_PI) * Ul * rhol );
//}

double MIsobar::rho3(double s) const {

  TF1 fdph("fdrho",this,&MIsobar::dph, 0,1 ,1,"MIsobar","dph");
  fdph.SetParameter(0,s);
  
  ROOT::Math::WrappedTF1 wph(fdph);  ROOT::Math::GaussIntegrator igph;
  igph.SetFunction(wph); igph.SetRelTolerance(1e-8);

  return igph.Integral(0,1);
}

cd MIsobar::rho3(cd s) const {

  //std::cout << "Limits " << real(pow(sqrtPi(s)-m3,2)) << " +i" << imag(pow(sqrtPi(s)-m3,2)) << std::endl;
  double par[] = {real(s), imag(s)};
  TF1 frdph("fdrho_re",this,&MIsobar::rdph, 0,1, 2,"MIsobar","rdph"); frdph.SetParameters(par);
  TF1 fidph("fdrho_im",this,&MIsobar::idph, 0,1, 2,"MIsobar","idph"); fidph.SetParameters(par);
  
  ROOT::Math::WrappedTF1 wph_re(frdph);  ROOT::Math::GaussIntegrator igph_re; igph_re.SetFunction(wph_re); igph_re.SetRelTolerance(1e-8);
  ROOT::Math::WrappedTF1 wph_im(fidph);  ROOT::Math::GaussIntegrator igph_im; igph_im.SetFunction(wph_im); igph_im.SetRelTolerance(1e-8);

  double real_part = igph_re.Integral(0,1);
  double imag_part = igph_im.Integral(0,1);

  cd value(real_part,imag_part);

  return value;

}

/* Triangle way to integrate */
double MIsobar::rdph(double *x, double *par) const {
  cd s(par[0],par[1]);
  cd th1 = pow(m1+m2,2);
  cd th2 = pow(sqrt(s)-m3,2);
  cd thM(pow(m1+m2,2) + (real(th2)-pow(m1+m2,2))/SCALEX_FOR_CUT,
	 imag(th2)/SCALEY_FOR_CUT);//(th2+th1)/2.0;

  if(x[0]<=1./2.) {
    cd s12 = th1 + 2*x[0]*(thM-th1);
    cd Ul = U(s12);
    cd rhol = rho(s,s12,m3*m3);
    return real( 2.*(thM-th1) / (2.*M_PI) * Ul * rhol );
  } else if(x[0]<=1.) { 
    cd s12 = thM + 2*(x[0]-1./2.)*(th2-thM);
    cd Ul = U(s12);
    cd rhol = rho(s,s12,m3*m3);
    return real( 2.*(th2-thM) / (2.*M_PI) * Ul * rhol );    
  } else {
    std::cerr<<"Error: inconsistency in MIsobar::rdph"<<std::endl;
    return 0.;
  }  
}

double MIsobar::idph(double *x, double *par) const {
  cd s(par[0],par[1]);
  cd th1 = pow(m1+m2,2);
  cd th2 = pow(sqrt(s)-m3,2);
  cd thM(pow(m1+m2,2) + (real(th2)-pow(m1+m2,2))/SCALEX_FOR_CUT,
	 imag(th2)/SCALEY_FOR_CUT);//(th2+th1)/2.0;

  if(x[0]<=1./2.) {
    cd s12 = th1 + 2*x[0]*(thM-th1);
    cd Ul = U(s12);
    cd rhol = rho(s,s12,m3*m3);
    return imag( 2.*(thM-th1) / (2.*M_PI) * Ul * rhol );
  } else if(x[0]<=1.0) { 
    cd s12 = thM + 2*(x[0]-1./2.)*(th2-thM);
    cd Ul = U(s12);
    cd rhol = rho(s,s12,m3*m3);
    return imag( 2.*(th2-thM) / (2.*M_PI) * Ul * rhol );    
  } else {
    std::cerr<<"Error: inconsistency in MIsobar::rdph"<<std::endl;
    return 0.;
  }
}
