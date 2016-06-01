//Misha Mikhasenko
//14.07.2015

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

#include "MIsobar.h"
#include "deflib.h"


MIsobar::MIsobar(double Mi, double G0i,
		   double m1i, double m2i, double m3i):
M(Mi), G0(G0i), m1(m1i), m2(m2i), m3(m3i) {;}
  
double MIsobar::rho(double s, double s12, double m3_2) const {
  return 1./(8*M_PI)*sqrt( LAMBDA(s,s12,m3_2)/(s*s) );
}
cd     MIsobar::rho(cd     s, cd     s12, double m3_2) const {
  return 1./(8*M_PI)*sqrtPi( LAMBDA(s,s12,m3_2)/(s*s) );
}

double MIsobar::U(double s) const {
  return 2.*M*G0/(pow(M*M-s,2)+pow(M*G0,2));
}
cd     MIsobar::U(cd     s) const {
  cd unit(0,1);
  return 2.*imag( 1./(M*M-s-unit*M*G0) );
}

//////////////////////////////////////////////////////////////////////////////                  
///  //   /  ///      //      ///     ///      ///      ///  ///////     /////
///  //      /////  ////  ///////  //////  //  ///  //  ///  ///////   ///////
///  //  /   /////  ////    /////  /  ///    /////      ///  /////////   /////
///  //  //  /////  ////      ///     ///  /   ///  //  ///      ///     /////
//////////////////////////////////////////////////////////////////////////////

double MIsobar::dph(double *x, double *par) const  {
  double s = par[0];
  double s12 = x[0];
  double Ul = U(s12);
  double rhol = rho(s,s12,m3*m3);
  return 1./(2.*M_PI)*Ul*rhol;
}

double MIsobar::rdph1(double *x, double *par) const {
  cd s(par[0],par[1]);
  cd s12(pow(m1+m2,2),x[0]);
  cd Ul = U(s12);
  cd rhol = rho(s,s12,m3*m3);
  return real( 1./(2.*M_PI)*Ul*rhol );
}
double MIsobar::idph1(double *x, double *par) const {
  cd s(par[0],par[1]);
  cd s12(pow(m1+m2,2),x[0]);
  cd Ul = U(s12);
  cd rhol = rho(s,s12,m3*m3);
  return imag( 1./(2.*M_PI)*Ul*rhol );
}
double MIsobar::rdph2(double *x, double *par) const {
  cd s(par[0],par[1]);
  cd s12( x[0], imag(pow(sqrtPi(s)-m3,2)) ); //I am not sure
  cd Ul = U(s12);
  cd rhol = rho(s,s12,m3*m3);
  return real( 1./(2.*M_PI)*Ul*rhol );
}
double MIsobar::idph2(double *x, double *par) const {
  cd s(par[0],par[1]);
  cd s12( x[0], imag(pow(sqrtPi(s)-m3,2)) ); //I am not sure
  cd Ul = U(s12);
  cd rhol = rho(s,s12,m3*m3);
  return imag( 1./(2.*M_PI)*Ul*rhol );
}

double MIsobar::rho3(double s) const {

  TF1 fdph("fdrho",this,&MIsobar::dph, pow(m1+m2,2), pow(sqrt(s)-m3,2),1,"MIsobar","dph");
  fdph.SetParameter(0,s);
  
  ROOT::Math::WrappedTF1 wph(fdph);  ROOT::Math::GaussIntegrator igph;
  igph.SetFunction(wph); igph.SetRelTolerance(0.001);

  return igph.Integral(pow(m1+m2,2),pow(sqrt(s)-m3,2));
}

cd MIsobar::rho3(cd s) const {

  //std::cout << "Limits " << real(pow(sqrtPi(s)-m3,2)) << " +i" << imag(pow(sqrtPi(s)-m3,2)) << std::endl;
  double par[] = {real(s), imag(s)};
  TF1 frdph1("fdrho1_re",this,&MIsobar::rdph1,            0, imag(pow(sqrtPi(s)-m3,2)),2,"MIsobar","rdph1"); frdph1.SetParameters(par);
  TF1 fidph1("fdrho1_im",this,&MIsobar::idph1,            0, imag(pow(sqrtPi(s)-m3,2)),2,"MIsobar","idph1"); fidph1.SetParameters(par);
  TF1 frdph2("fdrho2_re",this,&MIsobar::rdph2, pow(m1+m2,2), real(pow(sqrtPi(s)-m3,2)),2,"MIsobar","rdph2"); frdph2.SetParameters(par);
  TF1 fidph2("fdrho2_im",this,&MIsobar::idph2, pow(m1+m2,2), real(pow(sqrtPi(s)-m3,2)),2,"MIsobar","idph2"); fidph2.SetParameters(par);
  
  ROOT::Math::WrappedTF1 wph1_re(frdph1);  ROOT::Math::GaussIntegrator igph1_re; igph1_re.SetFunction(wph1_re); igph1_re.SetRelTolerance(0.001);
  ROOT::Math::WrappedTF1 wph1_im(fidph1);  ROOT::Math::GaussIntegrator igph1_im; igph1_im.SetFunction(wph1_im); igph1_im.SetRelTolerance(0.001);
  ROOT::Math::WrappedTF1 wph2_re(frdph2);  ROOT::Math::GaussIntegrator igph2_re; igph2_re.SetFunction(wph2_re); igph2_re.SetRelTolerance(0.001);
  ROOT::Math::WrappedTF1 wph2_im(fidph2);  ROOT::Math::GaussIntegrator igph2_im; igph2_im.SetFunction(wph2_im); igph2_im.SetRelTolerance(0.001);

  double r1 = igph1_re.Integral(0,imag(pow(sqrtPi(s)-m3,2)));
  double r2 = igph2_re.Integral(pow(m1+m2,2), real(pow(sqrtPi(s)-m3,2)));
  double i1 = igph1_im.Integral(0,imag(pow(sqrtPi(s)-m3,2)));
  double i2 = igph2_im.Integral(pow(m1+m2,2), real(pow(sqrtPi(s)-m3,2)));

  cd value(r1+r2,i1+i2);

  return value;

}
