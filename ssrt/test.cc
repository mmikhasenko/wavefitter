#include <iostream>
#include <complex>
#include <fstream>
#include <initializer_list>
#include <MIsobar.h>
#include <NoverD.h>
#include <constants.h>
#include <deflib.h>

#include <TROOT.h>
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>


using namespace std; 
typedef std::complex<double> cd;

TMultiGraph *mg(TGraph &g1,TGraph &g2);
TMultiGraph *mgm(std::initializer_list<TGraph*> a_args);
void Print(int Nbx, int Nby,cd *fs_2b, TH2D &t2d,ofstream &f);

cd rho(cd s);
cd rho3(cd s);
int main() {
  
  vector<double> vtest(5),vtest2(2);
  double dtest[7] = {1.1,2.3,3.3,4.4,5.6,8.6,9.7};
  vtest  = vector<double>(dtest,dtest+5);
  vtest2 = vector<double>(dtest+5,dtest+7);
  cout << vtest[4] << endl;
  cout << vtest2[1] << endl;
  return 0;

  gROOT->ProcessLine(".x ~/Documents/root-scripts/cv_n.C");

  MIsobar RhoPi(RHO_MASS,RHO_WIDTH,PI_MASS,PI_MASS,PI_MASS);
  const int N=100;
  double m_l=3*PI_MASS, m_r = 2;
  double step = (m_r-m_l)/(N-1);
  double x[N],y[N],yc[N];
  for(int i=0;i<N;i++) {
    x[i] = m_l+i*step;
    y [i] = real(RhoPi.rho3(x[i]*x[i]));
    yc[i] = (x[i]>RHO_MASS+PI_MASS) ? real(rho(x[i]*x[i])) : 0;
  }
  
  TGraph gr (N,x,y);  gr .SetLineWidth(2); gr.SetTitle("Phase space;#sqrt(s)");    
  TGraph grc(N,x,yc); grc.SetLineWidth(2); grc.SetLineColor(kRed); 
  TCanvas c1("c1","can",0,0,500,500);
  gr.Draw("alp"); grc.Draw("same");
  c1.SaveAs("/tmp/rho.png");
  
  MIsobar aRhoPi(RHO_MASS,RHO_WIDTH,PI_MASS,PI_MASS,PI_MASS,1,5);
  cout << aRhoPi.U(cd(0.5,0.0153)) << "  " << aRhoPi.U(cd(0.5,-0.0153)) << endl;

  //NoverD
  NoverD nd1(3,-16.,rho,pow(RHO_MASS+PI_MASS,2),1000,300);
  NoverD nd2(3,-16.,rho3,pow(3*PI_MASS,2)      ,1000,300);
//  nd1.npar[0]=8; nd1.npar[1]=200.; nd1.npar[2]=1000;
//  nd2.npar[0]=8; nd2.npar[1]=200.; nd2.npar[2]=1000;
  nd1.npar[0]=8*M_PI*647.04123; nd1.npar[1]=8*M_PI*2109.0199; nd1.npar[2]=8*M_PI*1719.2750;
  nd2.npar[0]=8*M_PI*647.04123; nd2.npar[1]=8*M_PI*2109.0199; nd2.npar[2]=8*M_PI*1719.2750;

  cout << "*************************************************************************" << endl;
  cout << "************************Calculation of real axis*************************" << endl;
  const int Nm=1000;
  double m1[Nm], m2[Nm];
  double reA1[Nm], reA2[Nm], imA1[Nm], imA2[Nm];
  double aA1[Nm], aA2[Nm], phA1[Nm], phA2[Nm];
  for(int i=0;i<Nm;i++) {
    //
    m1[i] = RHO_MASS+PI_MASS+(2.0-RHO_MASS-PI_MASS)/(Nm-1)*i;
    //Two body phase space;
    cd v1 = nd1.A(m1[i]*m1[i]);
    reA1[i] = real(v1); imA1[i] = imag(v1);
    aA1[i] = norm(v1); phA1[i] = atan2(imag(v1),real(v1));
    //Three body phase space
    m2[i] = 3*PI_MASS+(2.0-3*PI_MASS)/(Nm-1)*i;
    cd v2 = nd2.A(m2[i]*m2[i]);
    reA2[i] = real(v2); imA2[i] = imag(v2);
    aA2[i] = norm(v2); phA2[i] = atan2(imag(v2),real(v2));
  }
  TGraph g1re(Nm,m1,reA1); g1re.SetLineWidth(1.3); g1re.SetTitle("Real part of T;m(GeV)"); g1re.SetLineColor(kBlue);  g1re.SetLineStyle(2);
  TGraph g1im(Nm,m1,imA1); g1im.SetLineWidth(1.3); g1im.SetTitle("Imag part of T;m(GeV)"); g1im.SetLineColor(kGreen); g1im.SetLineStyle(2);                       
  TGraph g2re(Nm,m2,reA2); g2re.SetLineWidth(1.3); g2re.SetTitle("Real part of T;m(GeV)"); g2re.SetLineColor(kBlue);  g2re.SetLineStyle(1);
  TGraph g2im(Nm,m2,imA2); g2im.SetLineWidth(1.3); g2im.SetTitle("Imag part of T;m(GeV)"); g2im.SetLineColor(kGreen); g2im.SetLineStyle(1);
  //for(int i=0;i<Nm;i++) std::cout << reA2[i] << endl;
  TGraph g1a (Nm,m1, aA1); g1a .SetLineWidth(1.3); g1a .SetTitle("Intensity;m(GeV)");                                 g1a .SetLineStyle(2);
  TGraph g1ph(Nm,m1,phA1); g1ph.SetLineWidth(1.3); g1ph.SetTitle("Phase;m(GeV)");          g1ph.SetLineColor(kRed);   g1ph.SetLineStyle(2);                  
  TGraph g2a (Nm,m2, aA2); g2a .SetLineWidth(1.3); g2a .SetTitle("Intensity;m(GeV)");                                 g2a .SetLineStyle(1);
  TGraph g2ph(Nm,m2,phA2); g2ph.SetLineWidth(1.3); g2ph.SetTitle("Phase;m(GeV)");          g2ph.SetLineColor(kRed);   g2ph.SetLineStyle(1);

  TCanvas c2("c2","can2",0,0,1500,500);
  c2.Divide(3,1);
  c2.cd(1); mgm({&g1re,&g2re,&g1im,&g2im})->Draw("al");
  c2.cd(2); mgm({&g1a ,&g2a })->Draw("al"); 
  c2.cd(3); mgm({&g1ph,&g2ph})->Draw("al");
  c2.SaveAs("/tmp/func.png");
  c2.SaveAs("/tmp/func.pdf");

  //cout << "rho(1.5,-0.4) = " <<rho3(cd(1.5,-0.4)) << endl;

  cout << "*************************************************************************" << endl;
  cout << "************************Calculation of sheets****************************" << endl;
  const int Nbx=100, Nby=100;
  const double lx=0.2, rx=2.2, ly=-1.0, ry=1.0;
  TH2D t2d("t2d","t2d",Nbx,lx,rx,Nby,ly,ry); //t2d.Print("range");
  cd fs_2b[Nbx][Nby],ss_2b[Nbx][Nby],fs_q2b[Nbx][Nby],ss_q2b[Nbx][Nby];
  cd rho_fs[Nbx][Nby], rho_ss[Nbx][Nby];
  for(int i=1;i<=Nbx;i++) 
    for(int j=1;j<=Nby;j++) {
      cd s(t2d.GetXaxis()->GetBinCenter(i),t2d.GetYaxis()->GetBinCenter(j));
      rho_fs[i-1][j-1]=rho3(s);
      //if(i!=15||(j!=5&&j!=16)) 
      //continue;
      cout << "s = " <<s << ":" << endl;
      cd DI_2b   = nd1.DI (s);           cout << "------------->  DI_2b = "   <<   DI_2b << endl;  fs_2b[i-1][j-1]=  DI_2b;
      cd DII_2b  = DI_2b  - nd1.Disc(s); cout << "------------->  DII_2b = "  <<  DII_2b << endl;  ss_2b[i-1][j-1]= DII_2b;
      cd DI_q2b  = nd2.DI (s);           cout << "------------->  DI_q2b = "  <<  DI_q2b << endl; fs_q2b[i-1][j-1]= DI_q2b;
      cd DII_q2b = DI_q2b - nd2.Disc(s); cout << "------------->  DII_q2b = " << DII_q2b << endl; ss_q2b[i-1][j-1]=DII_q2b;      
    }
  
  ofstream tfout;
  tfout.open ("/tmp/2.5.out");
  tfout << "-----------> FIRST SHEET OF 2B" << endl;
  Print(Nbx,Nby,&fs_2b[0][0],t2d,tfout);
  tfout << "-----------> SECOND SHEET OF 2B" << endl;
  Print(Nbx,Nby,&ss_2b[0][0],t2d,tfout);
  tfout << "-----------> FIRST SHEET OF Q2B" << endl;
  Print(Nbx,Nby,&fs_q2b[0][0],t2d,tfout);
  tfout << "-----------> SECOND SHEET OF Q2B" << endl;
  Print(Nbx,Nby,&ss_q2b[0][0],t2d,tfout);
  tfout << "-----------> SECOND SHEET PHASE SPACE" << endl;
  Print(Nbx,Nby,&rho_fs[0][0],t2d,tfout);
  tfout.close();
  
  TH2D splane_fs_2b ( "splane_fs_2b","S-plane. First sheet 2-body"       ,Nbx,lx,rx,Nby,ly,ry);
  TH2D splane_ss_2b ( "splane_ss_2b","S-plane. Second sheet 2-body"      ,Nbx,lx,rx,Nby,ly,ry);
  TH2D splane_fs_q2b("splane_fs_q2b","S-plane. First sheet quasi-2-body" ,Nbx,lx,rx,Nby,ly,ry);
  TH2D splane_ss_q2b("splane_ss_q2b","S-plane. Second sheet quasi-2-body",Nbx,lx,rx,Nby,ly,ry);
  TH2D splane_rho_q2b("splane_rho_q2b","S-plane. Just quasi-2-body ph.sp.",Nbx,lx,rx,Nby,ly,ry);
  for(int i=1;i<=Nbx;i++)
    for(int j=1;j<=Nby;j++) {      
      int bin = t2d.GetBin(i,j);
      splane_rho_q2b.SetBinContent(bin, abs( rho_fs[i-1][j-1] ));    
      splane_fs_2b  .SetBinContent(bin, abs(  fs_2b[i-1][j-1] ));
      splane_ss_2b  .SetBinContent(bin, abs(  ss_2b[i-1][j-1] ));
      splane_fs_q2b .SetBinContent(bin, abs( fs_q2b[i-1][j-1] ));
      splane_ss_q2b .SetBinContent(bin, abs( ss_q2b[i-1][j-1] ));    
    }

  TCanvas c3("c3","can3",0,0,1000,1000);
  c3.Divide(2,2);

  c3.cd(1); splane_fs_2b .SetStats(kFALSE); splane_fs_2b .Draw("colz");
  c3.cd(2); splane_ss_2b .SetStats(kFALSE); splane_ss_2b .Draw("colz");
  c3.cd(3); splane_fs_q2b.SetStats(kFALSE); splane_fs_q2b.Draw("colz");
  c3.cd(4); splane_ss_q2b.SetStats(kFALSE); splane_ss_q2b.Draw("colz");
  c3.SaveAs("/tmp/structure.png");
  c3.SaveAs("/tmp/structure.pdf");

  TCanvas c4("c4","can4",0,0,700,700);
  splane_rho_q2b.GetZaxis()->SetRangeUser(0,0.03);
  splane_rho_q2b.SetStats(kFALSE); splane_rho_q2b .Draw("colz");
  c4.SaveAs("/tmp/rho.pdf");

  TFile fout("/tmp/sheets.root","recreate");

  splane_fs_2b .Write();
  splane_ss_2b .Write();
  splane_fs_q2b.Write();
  splane_ss_q2b.Write();
  c3.Write();
  splane_rho_q2b.Write();
  c4.Write();
  fout.Close();


  return 1;
}

cd rho(cd s) {
  if(real(s)>1000 && fabs(imag(s))<1e-5) return 1./(8*M_PI)*(1.0-(pow(RHO_MASS,2)+pow(PI_MASS,2))/s);
  cd phspc2 = (s-pow(RHO_MASS+PI_MASS,2))*(s-pow(RHO_MASS-PI_MASS,2))/(s*s);
  return 1./(8*M_PI)*sqrtPi(phspc2);
}

cd rho3(cd s) {
  if(real(s)>1000 && fabs(imag(s))<1e-5) return 1./(8*M_PI)*(1.0-(pow(RHO_MASS,2)+pow(PI_MASS,2))/s);
  MIsobar RhoPi(RHO_MASS,RHO_WIDTH,PI_MASS,PI_MASS,PI_MASS,5,5.0);
  cd value = RhoPi.rho3(s);
  //cout << "s = " << s << ", -> " << value << endl;
  return value;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

TMultiGraph *mg(TGraph &g1,TGraph &g2) {
  TMultiGraph *_mg = new TMultiGraph();
  _mg->Add(&g1);
  _mg->Add(&g2);
  return _mg;
}

TMultiGraph *mgm(std::initializer_list<TGraph*> a_args) {
  TMultiGraph *_mg = new TMultiGraph();
  for (auto i: a_args) _mg->Add(i);
  return _mg;
}

void Print(int Nbx, int Nby,cd *arr, TH2D &t2d, ofstream &f) { 
  for(int i=1;i<=Nbx;i++)
    for(int j=1;j<=Nby;j++) {
      int bin = t2d.GetBin(i,j);
      cd s(t2d.GetXaxis()->GetBinCenter(i),t2d.GetYaxis()->GetBinCenter(j));
      f << real(s) << " "
	<< imag(s) << " "
	<< real(arr[(i-1)*Nby+(j-1)]) << " "
	<< imag(arr[(i-1)*Nby+(j-1)]) << endl;
    }
}
 
