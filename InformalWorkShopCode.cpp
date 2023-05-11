/*
 *  InformalWorkShopCode.cpp
 *  
 *
 *  Created by obinna umeh on 28/02/2016.
 *  Copyright 2016 University of Cape Town. All rights reserved.
 *
 */

#include "InformalWorkShopCode.h"


#include "MyHIBias.h"

#include <boost/math/special_functions/legendre.hpp>


#include <boost/bind.hpp>
#include "MyBessel.cpp"

using boost::bind;
using boost::cref;


const double QMIN = 1e-4;
const double QMAX = 1e4;
const double XMAX = 0.999999;

const  double zmin = 1e-6;


double fnlini = 0.0;
double zini1 = 0.0;
double zini2 = 0.0;

LinearPower::LinearPower()
{
    //SetPowerSpectrum();
}

LinearPower::~LinearPower()
{
    //SetPowerSpectrum();
}

double T_b(double z)
{
	HIBias HIB;
	return HIB.Tb_Sp(z);
}

double Alpha( double k,double z)
{
	
	Cosmology C;
	double g0_overg_infty = 3.0/4.0;
	double c_light =  C.H_0*3000.0/C.h;
    // double c_light =  C.H_0*3000.0;
    const  double c_light_km_s =  299792;
	double numerator = 2.0*k*k*c_light*c_light*C.D_growth_Tabulated(z)*C.Tk_tab(k);
	double denominator = 3.0*C.Omega_m*pow(C.H_0,2) ;
	
	double tmp= g0_overg_infty*numerator/denominator;
	
	//cout <<"\t"<<"Alpha="<<tmp<<endl;
	return tmp;
}


double b_one(double fnl, double k, double z)
{
    HIBias HIB;
    Cosmology C;
    double delta_c = 1.686;
    double HIb10 = HIB.HI_b_ST_10(z);
	double HIb01 = HIB.HI_b_ST_01(fnl,z); 
    double DeltaB = 2.0* fnl *delta_c*(HIb10-1.0)/Alpha(k, z);
	// double tmp  = HIb10+ (HIb01/alpha_z(k, z));
    double tmp = HIb10 + DeltaB;
	// cout<<"\tAlpha="<<Alpha(k, z)<<"\t" <<"\t b_1=" <<"\t"<< tmp <<"\tT(k)="<<C.T_EH(k)<<endl;
    //return HIb10+ fnl *delta_c*(HIb10-1.0);
	return  tmp;
}
double Omz_in(double z)
{
    Cosmology C;
    double tmpz = C.Omega_m*pow((1.0+z),3)/(pow(C.E(z),2));
    return tmpz;
}
double CH2(double z)
{
    Cosmology C;
    //double c_light = C.H_0*3000.0;
     const  double c_light_km_s =  299792.0;
    // double c_light =  C.H_0*3000.0;
    return C.H(z)/(c_light_km_s*(1.0+z));
}
double theta(double z)
{
	double clightms = 299792458.0;
	double nu_21 = 1420.4*pow(10,6);//Hz
	double D_dish  = 15.0;//in meters
	double nu = nu_21/(1.0+z);
	double tmp  =1.22* clightms/(D_dish*nu)/M_PI;///Divided by Pi for dimensionless. 
	return tmp;
}

double WindowFun2(double z,double R, double k)
{
	double x =k*R/(1.0+z); 
	
	//double x = k*chi(z)*theta(z)/(sqrt(8.0*log(2.0)));
	
	
	if(x<1.0e-4)
	{
		return 1.0;
	}
	//note that W(Rk) is sensitive to this upper limit.
	else if(x>5000.)
	{
		return 0.0;
	}
	else
	{
		//return pow(3/(x*x*x),2)*pow(sin(x)-x*cos(x),2);
		double tmp =3.0*(sin(x)/pow(x,3)-cos(x)/pow(x,2));
		return tmp; 
	}
}



double f22(double z)
{
	Cosmology C;
	double tmpz = C.Omega_m*pow((1.0+z),3)/(pow(C.E(z),2));
	return pow(tmpz,4.0/7.0);
}

double LinearPower::fg_out(double z)
{
	return f22(z);
}

double Hpr_in(double z)
{
	//Cosmology C;
	double delta = 1.0e-4;
	double x1 = z -delta;
	double x2 = z + delta ;
	double y1 = CH2(x1);
	double y2 = CH2(x2);
	double deri = (y2-y1)/(x2-x1);
	double tmp = -(1.0 + z)*CH2(z)*deri;
	return tmp;
         //return -pow(CH(z)/2,2);///using CH' = -CH^2
}



double chi(double z)
{
     Cosmology C;
    double zmin1 = 1e-6;
     //double tmpchi = 0;
     double tmpchi = C.ComovingDistance(zmin1,z);
     return tmpchi;
}
double dchi(double z)
{
    Constants K;
    Cosmology C;
    double c_light =  C.H_0*3000.0;
    double coef =  ((K.c/1.0e3)/C.H_0);
    double tmp = coef/(C.E(z));
    return tmp;
	
	
	/*
	 double delta = 1.0e-5;///Need to be small for reasonal accuracy
	 double x1 = z -delta;
	 double x2 = z + delta ;
	 double y1 = chi(x1);
	 double y2 =chi(x2);
	 return (y2-y1)/(x2-x1);*/
}

double Delta_ISW_intrgnd(double z)
{
	Cosmology C;
	double Dz =   C.D_growth_Tabulated(z);
	double tmp1  = dchi(z)* pow(CH2(z),3)*Omz_in(z)*(f22(z)-1.0)*Dz;
	//cout << "\t"<<"ISW_integrand"<<"\t"<<tmp1<<endl;
	return tmp1;
}


double Delta_ISW(double z)
{ 
	
	Cosmology C;
	double Dz =    C.D_growth_Tabulated(z);
	const double EPSREL = 1e-10;
	const double EPSAB = 1e-10;
	double tmp2 = Integrate(bind(Delta_ISW_intrgnd,_1),zmin,z,EPSREL, EPSAB)/(Dz*pow(CH2(z),2));
	//cout << "\t"<<"ISW_int"<<"\t"<<tmp2<<endl;
	return tmp2;
}



double Sntgrd(double z, double k)
{
    Cosmology C;
    HIBias HIB;
    double Tb = T_b(z); 
    double HIb20 = HIB.HI_b_ST_20(z);
    double Dz = C.D_growth_Tabulated(z);
   //WindowFun(double theta,double z, double k)
    double coef = pow( Tb*HIb20,2)*pow(Dz,4)*k*k*C.P_k_EH_Tabulated(k)*C.P_k_EH_Tabulated(k)/(4.0*pow(M_PI,2));
    return coef;

}


double Sigmavintgrd(double z, double k)
{
    Cosmology C;
    double Dz = C.D_growth_Tabulated(z);
   //WindowFun(double theta,double z, double k)
    double coef = pow(Dz,2)*C.P_k_EH_Tabulated(k)/(6*pow(M_PI,2));

     return coef;

}

double Sigmaintgrd(double z, double k)
{
    Cosmology C;
    double Dz = C.D_growth_Tabulated(z);
   //WindowFun(double theta,double z, double k)
    double coef = pow(Dz,2)*k*k*C.P_k_EH_Tabulated(k)/(2*pow(M_PI,2));

     return coef;

}


double Mysigma1(double z)
{
      Cosmology C;
     //double R = 0.0;
      const double EPSREL = 1e-2;
      const double EPSAB = 1e-2;
     double upp = 2.0/(2.0 +C.n_s);
      double knl = 0.2*pow((1.0+z),upp);///in Mpc/h
       double QMIN = 1e-4;
       double QMAX = knl;
       // double QMAX = 1e4;

      double tmp1 = Integrate(bind(Sigmaintgrd,z,_1),QMIN,QMAX,EPSREL,EPSAB);
  // cout<<"\t i"<<"\t"<<"z"<<"\t"<<z <<"\t"<<"Sig_tab"<<"\t"<<"\t"<<tmp1 <<endl;
       return tmp1;
}
double MysigmaR(double z, double R)
{
      Cosmology C;
     //double R = 0.0;
      const double EPSREL = 1e-2;
      const double EPSAB = 1e-2;
     double upp = 2.0/(2.0 +C.n_s);
      double knl = 0.2*pow((1.0+z),upp);///in Mpc/h
       double QMIN = 1e-4;
       double QMAX = R;
       // double QMAX = 1e4;

      double tmp1 = Integrate(bind(Sigmaintgrd,z,_1),QMIN,QMAX,EPSREL,EPSAB);
   //cout<<"\t i"<<"\t"<<"z"<<"\t"<<z <<"\t"<<"Sig_tab"<<"\t"<<"\t"<<tmp1 <<endl;
       return tmp1;
}
double Mysigmav(double z)
{
      Cosmology C;
     //double R = 0.0;
      const double EPSREL = 1e-2;
      const double EPSAB = 1e-2;
     double upp = 2.0/(2.0 +C.n_s);
      double knl = 0.2*pow((1.0+z),upp);///in Mpc/h
       double QMIN = 1e-4;
       double QMAX = knl;

      double tmp1 = Integrate(bind(Sigmavintgrd,z,_1),QMIN,QMAX,EPSREL,EPSAB);
   //cout<<"\t i"<<"\t"<<"z"<<"\t"<<z <<"\t"<<"SigV_tab"<<"\t"<<"\t"<<tmp1 <<endl;
       return tmp1;
}


double NonlinearShotnoise(double z)
{
      const double EPSREL = 1e-4;
      const double EPSAB = 1e-4;
       double QMIN = 1e-4;
       double QMAX = 10e4;
      double tmp2 = Integrate(bind(Sntgrd,z,_1),QMIN,QMAX,EPSREL,EPSAB);
       cout<<"\t i"<<"\t"<<"z"<<"\t"<<z <<"\t"<<"ShotNoise_tab"<<"\t"<<"\t"<<tmp2 <<endl;
       return tmp2;
}

double Mysigmaderivative(double z)
{
        double delta = 1.0e-4;
	double x1 = z -delta;
	double x2 = z + delta ;
	double y1 = Mysigma1(x1);
	double y2 = Mysigma1(x2);
	double deri = (y2-y1)/(x2-x1);
	double tmp = deri;
	return tmp;
}

double LinearPower::HproverH(double z)
{
	//Cosmology C;
	double delta = 1.0e-4;
	double x1 = z -delta;
	double x2 = z + delta ;
	double y1 = CH2(x1);
	double y2 = CH2(x2);
	double deri = (y2-y1)/(x2-x1);
	double tmp = -(1.0 + z)*deri/CH2(z);
	return tmp;
         //return -pow(CH(z)/2,2);///using CH' = -CH^2
}



int k_bins2 =100;
std::vector <double> z_Dz_tab2(k_bins2);
std::vector <double> Lens_tab(k_bins2);
std::vector <double> f_tab(k_bins2);
std::vector <double> ISW_tab(k_bins2);
std::vector <double> OM_tab(k_bins2);
std::vector <double> jz_tab(k_bins2);
std::vector <double> Tb_tab2(k_bins2);

std::vector <double> k_arrayin(k_bins2);


std::vector <double> Hpr_tab(k_bins2);
 std::vector <double> Sig_tab1(k_bins2);
 std::vector <double> Sigv_tab1(k_bins2);
 std::vector <double> Sigderi_tab1(k_bins2);

 std::vector <double> Sn_tab(k_bins2);



void LinearPower::InitialzieSpeed4(void)
{
	   HIBias HIB;
	//const int k_bins2 = 100;
	double  EPSREL = 1e-8;
	double z_Dz_min2;
	double z_Dz_max2;
	z_Dz_min2 = log(1e-3);
	z_Dz_max2 = log(10.0);
	const double QMIN = 1e-4;
	const double QMAX = 1e5;
	
	
	
	
#pragma omp parallel private(thread_number)
	{
#pragma omp for schedule(static) nowait
		for(int i = 0;i< k_bins2 ;i++)
		{
			z_Dz_tab2[i] = ((z_Dz_max2-z_Dz_min2)*((double) i)/((double)  k_bins2-1) + z_Dz_min2 );
			//k_arrayin[i] = exp( ((log(QMAX)-log(QMIN))*((double) i)/((double) k_bins2-1) + log(QMIN) ));
			Hpr_tab[i] = Hpr_in(exp(z_Dz_tab2[i]));
			f_tab[i]  = f22(exp(z_Dz_tab2[i]));
			//ISW_tab[i]  =   Delta_ISW(exp(z_Dz_tab2[i]));
                        Sig_tab1[i] =  Mysigma1(exp(z_Dz_tab2[i]));
                        //Sigv_tab1[i] =  Mysigmav(exp(z_Dz_tab2[i]));
                        Tb_tab2[i] = HIB.T_bST(exp(z_Dz_tab2[i]));
			OM_tab[i]  =   Omz_in(exp(z_Dz_tab2[i]));//Mysigma2(int x,double fnl,double R,double z)
                         Sigderi_tab1[i]= Mysigmaderivative(exp(z_Dz_tab2[i])); //NonlinearShotnoise(double z)

			Sn_tab[i] = NonlinearShotnoise(exp(z_Dz_tab2[i]));
		}
		
	}
}


double f2(double z)
{
	Spline<double, double> CubicSpline_f(z_Dz_tab2,f_tab);
	return CubicSpline_f.interpolate(log(z));
   	
}

double Omz(double z)
{
	Spline<double, double> CubicSpline_OM(z_Dz_tab2,OM_tab);
	return CubicSpline_OM.interpolate(log(z));
   	
}

double Hpr(double z)
{
	Spline<double, double> CubicSpline_Hpr(z_Dz_tab2,Hpr_tab);
	return CubicSpline_Hpr.interpolate(log(z));
   	
}

double Tb_spline(double z)
{
       Spline<double, double> CubicSpline_Tb(z_Dz_tab2,Tb_tab2);
       double tmp =  CubicSpline_Tb.interpolate(log(z));
 // cout <<"t_b="<<tmp<<endl;
     return tmp*1e3;
   	
}

double Sigma_out(double z)
{
	Spline<double, double> CubicSpline_sigma(z_Dz_tab2,Sig_tab1);
	return CubicSpline_sigma.interpolate(log(z));
   	
}
double Sigmav_out(double z)
{
	//Spline<double, double> CubicSpline_sigmav(z_Dz_tab2,Sigv_tab1);
	return 0.0;//CubicSpline_sigmav.interpolate(log(z));
   	
}
double Sigmaderi_out(double z)
{
	Spline<double, double> CubicSpline_sigma(z_Dz_tab2,Sigderi_tab1);
	return CubicSpline_sigma.interpolate(log(z));
   	
}

double ShotNoise(double z)
{
	Spline<double, double> CubicSpline_Sn(z_Dz_tab2,Sn_tab);

	return CubicSpline_Sn.interpolate(log(z));
   	
}


double LinearPower::ShotNoiseout(double z)
{
 return ShotNoise(z);
}

double LinearPower::sigmasqout(double z)
{
 return NonlinearShotnoise(z);//Sigma_out(z);
}

double Delta_GRISW_intrgnd(double z)
{
	Cosmology C;
	double Dz =   C.D_growth_Tabulated(z);
	double tmp1  =  dchi(z)*pow(CH2(z),2)*Omz_in(z)*Dz;
	return tmp1;
}



double Delta_GRISW(double z)
{ 
	
	Cosmology C;
	double Dz =    C.D_growth_Tabulated(z);
	const double EPSREL = 1e-8;
	const double EPSAB = 1e-8;
	double tmp2 = Integrate(bind(Delta_GRISW_intrgnd,_1),zmin,z,EPSREL, EPSAB);
	double tmp = tmp2;
	return tmp;
}
double bevo(double z)
{
	HIBias HIB;
	double delta_c = 1.686;
	//double tmp = delta_c*(HIB.HI_b_ST_10(z)-1.0);
        double DeltaT = HIB.HI_b_ST_20(z)*Sigmav_out(z);
        double tmp1 = (1.0 +  DeltaT);
        double tmp = HIB.be_Sp(z) - CH2(z)*(1.0+z)*Sigmaderi_out(z)/tmp1;
	return tmp;
}


double CurlyB(double Q, double z)
{
	double tmp = -f2(z)*(2.0-bevo(z)/CH2(z) + Hpr_in(z)/pow(CH2(z),2));
	return tmp;
}


double CurlyBmatter(double Q, double z)
{
	double tmp = -f2(z)*(2.0 + Hpr_in(z)/pow(CH2(z),2));
	return tmp;
}

double CurlyA(double Q, double z)
{
	double first = f2(z)*(3.0-bevo(z)/CH2(z) - 3.0*Omz_in(z)/2.0) ;
	double second  = -(3.0*Omz_in(z)/2.0);/// -3.0*ISW(z));
	double mult_term =(2.0-bevo(z)/CH2(z) + Hpr_in(z)/pow(CH2(z),2));
	double tmp1 = (first + second*mult_term);
	return tmp1;
}
double kernelZ1(double fnl, double k, double z , double mu)
{
	double Q = 1.0;
	double tmp  = b_one(fnl, k, z)+ mu*mu*f2(z)+ CurlyA(Q,z)*pow(CH2(z)/k,2);
	return tmp;
	
}



double EvoCurlyB(double be, double z)
{
	double tmp = -f2(z)*(2.0-be/CH2(z)+ Hpr_in(z)/pow(CH2(z),2));
	return tmp;
}



double EvokernelZ1(double fnl,double k, double z , double mu)//EvokernelZ1
{ 
       HIBias HIB;
	 
	double HIb10 = HIB.HI_b_ST_10(z); //b_one(fnl, k, z)
      
	double tmp  = b_one(fnl, k, z)  + mu*mu*f2(z)+ CurlyA(1.0,z)*pow(CH2(z)/k,2);
	return tmp;
	
}


double Omegampr(double z)
{
       double delta = 1.0e-4;
	double x1 = z -delta;
	double x2 = z + delta ;
	double y1 = Omz(x1);
	double y2 = Omz(x2);
	double deri = (y2-y1)/(x2-x1);
	double tmp = -(1.0 + z)*CH2(z)*deri;
	return tmp;
}
double gamma4(double z)
{

   double Q = 1.0;
      double first = 2.0*f2(z)-1.0-bevo(z)/CH2(z) + 4.0*Q  + 2.0*Hpr(z)/pow(CH2(z),2);
	double second  = 3.0*Omegampr(z)/(2.0*CH2(z));
	double mult_term = (3.0*Omz(z)/2.0)*first;
	double tmp1 = (mult_term + second);
	return tmp1;
}
double gamma5(double z)
{
double tmp = - (bevo(z)/CH2(z)-3)*f2(z);
return tmp;
}

double CurlyAmatter(double Q, double z)
{
       
	double first = f2(z)*(3.0 - 3.0*Omz_in(z)/2.0) ;
	double second  = -(3.0*Omz_in(z)/2.0);/// -3.0*ISW(z));
	double mult_term =(2.0 + Hpr_in(z)/pow(CH2(z),2));
	double tmp1 = (first+ second*mult_term);
	return tmp1;
}

double EvokernelZ1matter(double fnl,double k, double z , double mu)
{ 
       HIBias HIB;
	 
	//double HIb10 = HIB.HI_b_ST_10(z); 
      
	double tmp  = 1.0 + mu*mu*f2(z)+ CurlyAmatter(1.0,z)*pow(CH2(z)/k,2);
	return tmp;
	
}



double LinearPower::LinearPower_HIHI(double fnl, double k, double z1,double mu)
{
	Cosmology C;
	double Tb = T_b(z1);  
	//cout <<"T_b="<<"\t"<<Tb<<endl;
     double tmp  = pow(Tb*C.D_growth_Tabulated(z1)*(b_one(fnl, k, z1)+ mu*mu*f2(z1)),2)*C.P_k_EH_Tabulated(k);
	//double tmp = pow(Tb*(b_one(fnl, k, z1) + mu*mu*f2(z1)),2)*C.P_EH(k,z1);
	return  tmp;
}

double LinearPower::Adhocbias(double fnl, double k, double z1,double mu)
{
  return b_one(fnl, k, z1);
}


double LinearPower::LinearPower_HIHIRSD(double fnl, double k, double z1,double mu)
{
	Cosmology C;
	HIBias HIB;
	double Tb = T_b(z1); 
	//double HIb10 = HIB.HI_b_ST_10(z1); 
        double Dz = C.D_growth_Tabulated(z1);
	double pkm =  pow(Dz,2)*C.P_k_EH_Tabulated(k);
	// double fnlin  = fnl -5.0/3.0; 
    // double tmp = pow(Tb*(HIb10+mu*mu* f2(z1)),2)*C.P_EH(k,z1);EvokernelZ1(double fnl,double k, double z , double mu)
       double tmp =  pow(EvokernelZ1(fnl, k, z1 , mu),2)+pow(mu*CurlyB(1.0, z1)*(CH2(z1)/k),2);
	return  pow(Tb,2)*tmp*pkm;
}


double LinearPower::LinearPower_HIHIRSDmatter(double fnl, double k, double z1,double mu)
{
	Cosmology C;
	HIBias HIB;
	double Tb = T_b(z1); 
	double Dz = C.D_growth_Tabulated(z1);
      double power_multi= pow(Dz,2)*C.P_k_EH_Tabulated(k);
      double tmp =  pow(EvokernelZ1matter(fnl, k, z1 , mu),2)+pow(mu*CurlyBmatter(1.0, z1)*(CH2(z1)/k),2);
	return  Tb*Tb*tmp*power_multi;
}


double ThirdOrderP(double fnl, double k, double z1,double mu)
{
	Cosmology C;
	HIBias HIB;
	double Tb = T_b(z1); 
	double HIb10 = HIB.HI_b_ST_10(z1); 
	// double fnlin  = fnl -5.0/3.0; 
    // double tmp = pow(Tb*(HIb10+mu*mu* f2(z1)),2)*C.P_EH(k,z1);
    double tmp =  kernelZ1(fnl, k, z1 , mu);
	return  tmp*C.P_EH(k,z1);
}


double LinearPower::CompareKaizer(double fnl, double k, double z1,double mu)
{
	//AllPower_HIHIRSD(int x1, int x2,double fnl, double k, double z,double mu)
	//double tmp = (LinearPower_HIHIRSD(fnl, k, z1,1.0)-LinearPower_HIHI(fnl, k, z1))/LinearPower_HIHI(fnl, k, z1);
	double tmp = (LinearPower_HIHIRSD(fnl, k, z1,mu)-LinearPower_HIHI(fnl, k, z1,mu))/LinearPower_HIHI(fnl, k, z1,mu);
	return tmp;
}


double LinearPower::LinearPower_HIHI_intgrand(int l, double fnl, double k, double z,double mu)
{
	
	return (2*l+1)*LinearPower_HIHIRSD(fnl,k, z,mu)*boost::math::legendre_p(l, 0, mu)/2.0;
	
}

double LinearPower::LinearPower_HIHI_intgrand2(int l, double fnl, double k, double z,double mu)
{
	
	return (2*l+1)*LinearPower_HIHI(fnl,  k, z,mu)*boost::math::legendre_p(l, 0, mu)/2.0;
	
}
double LinearPower::LinearPower_HIHI_intgrandmatter(int l, double fnl, double k, double z,double mu)
{
	
	//return (2*l+1)*LinearPower_HIHIRSDmatter(fnl,k, z,mu)*boost::math::legendre_p(l, 0, mu)/2.0; 
     return LinearPower_HIHIRSDmatter(fnl,k, z,mu)/2.0;
	
}




double LinearPower::LinearHIHI_Average( int l, double fnl, double k, double z)
{
	const double EPSREL = 1e-5;
	const double EPSABS = 1e-5;
	const  double mu_min = -1.0;
	const double mu_max = 1.0;
	double tmp = Integrate(bind(&LinearPower::LinearPower_HIHI_intgrand2,this,l,fnl,k,z,_1),mu_min,mu_max,EPSREL,EPSABS);
	return tmp;
}
double LinearPower::LinearHIHI_AverageRSD( int l, double fnl, double k, double z)
{
	const double EPSREL = 1e-7;
	const double EPSABS = 1e-8;
	const  double mu_min = -1.0;
	const double mu_max = 1.0;
	double tmp = Integrate(bind(&LinearPower::LinearPower_HIHI_intgrand,this,l,fnl,k,z,_1),mu_min,mu_max,EPSREL,EPSABS);
	return tmp;
}

double LinearPower::LinearHIHI_AverageRSDmatter(int l, double fnl, double k, double z)
{
	const double EPSREL = 1e-7;
	const double EPSABS = 1e-8;
	const  double mu_min = -1.0;
	const double mu_max = 1.0;
	double tmp = Integrate(bind(&LinearPower::LinearPower_HIHI_intgrandmatter,this,l,fnl,k,z,_1),mu_min,mu_max,EPSREL,EPSABS);
	return tmp;
}

double LinearPower::k_H(double z)
{
	Cosmology C;
	double c_light = C.H_0*3000.0;//in  km s^-1
	const  double c_light_km_s =  299792.0;
	double k =  (1.0/c_light_km_s )*C.H(z)/(1.0+z);
	//double k = C.H(z)/(1.0+z);
	return k;
}

double LinearPower::LinearPower_mm(double k, double z)
{   
	Cosmology C;
	double Dz = C.D_growth_Tabulated(z);
	return pow(Dz,2)*C.P_k_EH_Tabulated(k);
}






double WindowFun(double z, double k)
{
Cosmology C;
double upp = 2.0/(2.0 +C.n_s);
//double R = 10.0*pow((1.0+z),-upp);
double R = 5.0*pow((1.0+z),-upp);///Smith et all 2002
///double upp = 2.0/(2.0 +C.n_s);
 //     double knl = 0.2*pow((1.0+z),upp);///in h/Mpc
double x = k*R;///(1.0+z);
if(x<1.0e-5)
{
		return 1.0;
}

	//note that W(Rk) is sensitive to this upper limit.
else if(x>50000.)
{
		return 0.0;
}
else 
{
	//double top_hat = 3.0*(sin(x)/pow(x,3)-cos(x)/pow(x,2));// pow(3/(x*x*x),2)*pow(sin(x)-x*cos(x),2);
   double Guassian = exp(-x*x/2.0);
//cout <<"\t"<<"top_hat"<<"\t"<< top_hat<<endl;	
        return  Guassian;//top_hat;//
}

}





double IntegrandP13(double z, double k, double r)
{
       Cosmology C;
       double Dz = C.D_growth_Tabulated(z);
       double powermm = pow(Dz,2)*C.P_k_EH_Tabulated(k*r);
       double powermmy = pow(Dz,2)*C.P_k_EH_Tabulated(k);
       double numerator = 3.0*pow((pow(r,2)-1.0),3)/pow(r,2);
       double tmp1 = 12/pow2(r) - 158 + 100*pow2(r) - 42*pow4(r) + 3/pow3(r) * pow3(r*r - 1) * (7*r*r + 2) * log((1+r)/fabs(1-r));
      /// double tmp1 = 12.0/(r*r) - 158.0 + 100.0*r*r - 42.0*pow(r,4) +  numerator*(7.0*r*r + 2.0)*log(abs((1.0 + r)/(1.0 - r)));
       double prefactor = pow(k,3)*powermmy/(252.0*pow(2.0*M_PI,2));
       double tmp2 = powermm*tmp1*prefactor;
       return tmp2;
}

double IntegrandP22(double z, double k, double r, double mu)
{
      Cosmology C;
        double y = sqrt(r*r-2.0*mu *r + 1.0);
       double Dz = C.D_growth_Tabulated(z);
       double powermmr = pow(Dz,2)*C.P_k_EH_Tabulated(k*r);
      double powermmy = pow(Dz,2)*C.P_k_EH_Tabulated(k*y);
      double tmp1 = pow((3.0* r + 7.0 * mu - 10.0 * r*mu*mu),2)/pow((1.0+ r*r - 2.0* r*mu),2);
    double prefactor = pow(k,3)/(98.0*pow(2.0*M_PI,2));
     double  tmp2 = powermmr*powermmy*prefactor;
    return tmp2;

}

double F2P(double k, double k1, double mu1)
{


   double r = k1/k;
     double y = sqrt(r*r-2.0*mu1 *r + 1.0);
    // double k2 = k*y;
    // double k1mod = k*r;
   double first_term = 5.0/7.0;
   double second_term = (r/y + y/r)*(r*mu1-r*r)/(2.0*r*y);
   double third_term = (2.0/7.0)*pow((r*mu1-r*r)/(r*y),2);
   double tmp =  first_term+second_term+third_term;
   return 2.0*tmp;
}

double PowerSpectrumIntegrand(double z, double k,double k1,double mu)
    {
     Cosmology C;
      double r = k1/k;
     double y = sqrt(r*r-2.0*mu *r + 1.0);
     double Dz = C.D_growth_Tabulated(z);
     double power_multi= pow(Dz,4)*C.P_k_EH_Tabulated(k*r)*C.P_k_EH_Tabulated(k*y);
     double coeff =  k1*k1/(4.0*pow(M_PI,2));
     double tmp = coeff*F2P(k,k*r,mu)*F2P(k,k*r,mu)*power_multi;
    /// cout<<"\t Integrand="<<tmp1<<endl;
    return tmp;
}



double  LinearPower::PowerSpectrum22(double k,double z)
{

     Cosmology C;
      const double EPSREL = 1e-5;
     const double EPSABS = 1e-5;
     // const double QMIN = 1e-4;
      //const double QMAX = 1e4;
     // const double XMAX = 0.98;
     double a[2] = {QMIN, -1};
     double b[2] = {QMAX, XMAX};
     double tmp= Integrate<2>(bind(PowerSpectrumIntegrand,z,k,_1,_2),a, b, EPSREL,EPSABS);
    // double tmp = Integrate<2>(bind(IntegrandP22,z,k,_1,_2),a, b,EPSREL,EPSABS);
 //cout<<"Integrating="<<tmp<<endl;
    return tmp;
    
}


static double f22_dd(double k, double q, double x)
 {
  Cosmology C;//C.P_k_EH_Tabulated(k1)

    double r = q/k;
   double d = 1 + r*r - 2*r*x;
    if(d < 1e-5)
        return 0;
    else
        return C.P_k_EH_Tabulated(q)*C.P_k_EH_Tabulated(k*sqrt(d))*pow2(3*r + 7*x - 10*r*x*x)/pow2(d);
}

/* P_{\delta\delta}^{(22)} */
double LinearPower::PowerSpectrum222(double k, double z) 
 {
      Cosmology C;
  const double EPSREL = 1e-5;
     const double EPSABS = 1e-5;
     double Dz = C.D_growth_Tabulated(z);
    double a[2] = {QMIN, -1 };
    double b[2] = {QMAX, XMAX };
    return pow4(Dz)*k*k/(98*4*M_PI*M_PI)*Integrate<2>(bind(f22_dd, k, _1, _2), a, b, EPSREL, EPSABS*C.P_k_EH_Tabulated(k));
}















double  LinearPower::PowerSpectrum13(double k,double z)
{

     Cosmology C;
      const double EPSREL = 1e-5;
     const double EPSABS = 1e-5;
     // const double QMIN = 1e-5;
     // const double QMAX = 1e-1;
     // const double XMAX = 0.98;
   //  double a[2] = {QMIN,-1.0};
   // double  b[2] = {QMAX,0.99};
    
     double tmp = Integrate(bind(IntegrandP13,z,k,_1),QMIN,QMAX,EPSREL,C.P_k_EH_Tabulated(k)*EPSABS);
// cout<<"Integrating="<<tmp<<endl;
    return tmp;
    
}






///////////////////////Copter tools////////////////////////////////////////////



static double f13_dd(double k, double q) 
{
     Cosmology C;   
    double r = q/k;
    double s;
   
      
      
    if(r < 1e-2)
        s = -168 + (928./5.)*pow2(r) - (4512./35.)*pow4(r) + (416./21.)*pow6(r);
    else if(fabs(r-1) < 1e-10)
        s = -88 + 8*(r-1);
    else if(r > 100)
        s = -488./5. + (96./5.)/pow2(r) - (160./21.)/pow4(r) - (1376./1155.)/pow6(r);
    else
        s = 12/pow2(r) - 158 + 100*pow2(r) - 42*pow4(r) + 3/pow3(r) * pow3(r*r - 1) * (7*r*r + 2) * log((1+r)/fabs(1-r));

    return C.P_k_EH_Tabulated(q) * s;
}

/* P_{\delta\delta}^{(13)} */
double LinearPower::P13_dd(double k,double z) 
const {
///ThirdOrderP(double fnl, double k, double z1,double mu)
   // Cosmology C;
   //  double Dz = C.D_growth_Tabulated(z);
    //const double EPSREL = 1e-5;
    // const double EPSABS = 1e-5;
   // return -pow4(Dz)*k*k/(252*4*M_PI*M_PI) * C.P_k_EH_Tabulated(k) * Integrate<ExpSub>(bind(f13_dd, k, _1), QMIN, QMAX,EPSREL, 1e-4*C.P_k_EH_Tabulated(k));

 //return pow4(Dz)*k*k/(252*4*M_PI*M_PI) *  C.P_k_EH_Tabulated(k) * Integrate<ExpSub>(bind(f13_dd, k, _1), QMIN, QMAX,EPSREL, C.P_k_EH_Tabulated(k)*EPSABS);
return P13_ddint(k,z);
}
double P13_ddint(double k,double z) 
 {
///ThirdOrderP(double fnl, double k, double z1,double mu)
    Cosmology C;
     double Dz = C.D_growth_Tabulated(z);
    const double EPSREL = 1e-5;
     const double EPSABS = 1e-5;
   // return -pow4(Dz)*k*k/(252*4*M_PI*M_PI) * C.P_k_EH_Tabulated(k) * Integrate<ExpSub>(bind(f13_dd, k, _1), QMIN, QMAX,EPSREL, 1e-4*C.P_k_EH_Tabulated(k));

 double tmp = pow4(Dz)*k*k/(252*4*M_PI*M_PI)*C.P_k_EH_Tabulated(k)*Integrate(bind(f13_dd, k, _1), QMIN, QMAX,EPSREL, C.P_k_EH_Tabulated(k)*EPSABS);

 return  tmp;
}






//sigma_R(double R)


double SecondOrderKernelP(int x,double z, double k,double k1,double mu)
{     HIBias HIB;
     double HIb10 = HIB.HI_b_ST_10(z); 
     double HIb20 = HIB.HI_b_ST_20(z);
double r = k1/k;
     double y = sqrt(r*r-2.0*mu *r + 1.0);
if(x == 1)
{
return HIb10*F2P(k,k*r,mu);//F2P(k,k*r,mu);F2P(double k, double k1, double mu1)
}
else if(x ==2)
{
return HIb20;
}
else
{
double tmp  = HIb10*F2P(k,k*r,mu) + HIb20;
return tmp;
}
}


double S_kernelP(double z, double k, double k1, double mu)
{
    double r = k1/k;
     double y = sqrt(r*r-2.0*mu *r + 1.0);
    Cosmology C;
  // double beta = -3.0*Omz_in(z);//*pow(CH2(z),2);
   double tidal = (r*mu-r*r)*(r*mu-r*r)/((r*y)*(r*y)) - 1.0/3.0;
    // double tidal  = ( pow(k1*(r*mu - r*r),2) - 5.0 *pow( k1*k*y,2)/9.0)/(Alpha(k*y, z)*Alpha(k1, z));
       return tidal;

}


 double TidalBias(double z)
 {  
   HIBias HIB;
    double HIb10 = HIB.HI_b_ST_10(z);
    double tmp = -2.0*(HIb10-1.0)/7.0;
    return tmp;
 }

double SecondOrderKernelPtidal(int x,double z, double k,double k1,double mu)
{     HIBias HIB;
     double HIb10 = HIB.HI_b_ST_10(z); 
     double HIb20 = HIB.HI_b_ST_20(z);
double r = k1/k;
     double y = sqrt(r*r-2.0*mu *r + 1.0);
if(x == 1)
{
return HIb10*F2P(k,k*r,mu);//F2P(k,k*r,mu);
}
else if(x ==2)
{
return HIb20;
}
else if(x ==3)
{
return TidalBias(z)*S_kernelP(z, k, k1, mu);
}
else
{
double tmp  = HIb10*F2P(k,k*r,mu) + HIb20;// + TidalBias(z)*S_kernelP(z, k, k1, mu);
return tmp;
}
}


double PowerSpectrumIntegrandHItidal(int x1,int x2, double z, double k,double k1,double mu)
    {
     Cosmology C;
      double r = k1/k;
     double y = sqrt(r*r-2.0*mu *r + 1.0);
     double Dz = C.D_growth_Tabulated(z);
     double power_multi= pow(Dz,4)*C.P_k_EH_Tabulated(k*r)*C.P_k_EH_Tabulated(k*y);
     double coeff =  k1*k1/(4.0*pow(M_PI,2));
     double tmp = coeff*SecondOrderKernelPtidal(x1,z, k,k*r,mu)*SecondOrderKernelPtidal(x2,z, k,k*r,mu)*power_multi;
    /// cout<<"\t Integrand="<<tmp1<<endl;
    return tmp;
}



double PowerSpectrumIntegrandHI(int x1,int x2, double z, double k,double k1,double mu)
    {
     Cosmology C;
      double r = k1/k;
     double y = sqrt(r*r-2.0*mu *r + 1.0);
     double Dz = C.D_growth_Tabulated(z);
     double power_multi= pow(Dz,4)*C.P_k_EH_Tabulated(k*r)*C.P_k_EH_Tabulated(k*y);
     double coeff =  k1*k1/(4.0*pow(M_PI,2));
     double tmp = coeff*SecondOrderKernelP(x1,z, k,k*r,mu)*SecondOrderKernelP(x2,z, k,k*r,mu)*power_multi;
    /// cout<<"\t Integrand="<<tmp1<<endl;
    return tmp;
}


double  LinearPower::Non_LinearPower_HIHItidal(int x1,int x2, double k,double z1)
{

      Cosmology C;
      const double EPSREL = 1e-5;
     const double EPSABS = 1e-5;
    double a[2] = {QMIN, -1 };
    double b[2] = {QMAX, XMAX };
  
     double tmp = Integrate<2>(bind(PowerSpectrumIntegrandHItidal,x1,x2,z1,k,_1,_2),a, b,EPSREL,EPSABS);
    return tmp;
  
}


 double  LinearPower::Non_LinearPower_HIHI(int x1,int x2, double k,double z1)
{

      Cosmology C;
      const double EPSREL = 1e-5;
     const double EPSABS = 1e-5;
    double a[2] = {QMIN, -1 };
    double b[2] = {QMAX, XMAX };
  
     double tmp = Integrate<2>(bind(PowerSpectrumIntegrandHI,x1,x2,z1,k,_1,_2),a, b,EPSREL,EPSABS);
    return tmp;
  
}

double SecondOrderCrosstermintegrang(int x,double z1,double z2, double k,double k1,double mu)
{
       Cosmology C;
      double r = k1/k;
     double y = sqrt(r*r-2.0*mu *r + 1.0);
     double Dz1 = C.D_growth_Tabulated(z1);
     double Dz2 = C.D_growth_Tabulated(z2);
     double power_multi= pow(Dz1*Dz2,2)*C.P_k_EH_Tabulated(k*r)*C.P_k_EH_Tabulated(k*y);
     double coeff =  k1*k1/(4.0*pow(M_PI,2));

if(x == 1)
{
double tmp1 = (SecondOrderKernelP(10,z1, k, k1,mu)*f2(z2)*GRSD_2(k, k1, mu) + SecondOrderKernelP(10,z2, k, k1,mu)*f2(z1)*GRSD_2( k, k1, mu));
return coeff*power_multi*tmp1;///*WindowFun(z1, k);//;
}
else if(x ==2)
{
double tmp2 = f2(z1)*GRSD_2(k, k1, mu)*f2(z2)*GRSD_2(k,k1, mu);
return coeff*power_multi*tmp2;///*WindowFun(z1, k);//;
}
else
{
cout<<"Bad choice"<<endl;
return 0.0;
}
}


 double  Non_LinearPower_CrossTerms(int x1, double k,double z1,double z2)
{

      Cosmology C;
      const double EPSREL = 1e-5;
     const double EPSABS = 1e-5;
    double a[2] = {QMIN, -1 };
    double b[2] = {QMAX, XMAX };
  
     double tmp = Integrate<2>(bind(SecondOrderCrosstermintegrang,x1,z1,z2,k,_1,_2),a, b,EPSREL,EPSABS);
    return tmp;
  
}

double LinearPower::HIPowerSpectrum(int x, double k,double z) 
{

     Cosmology C;
      HIBias HIB;
     // double upp = -2.0/(2.0 +C.n_s);
     // double R = 5.0*pow((1.0+z),upp);;///in Mpc/h
     double HIb10 = HIB.HI_b_ST_10(z); 
     double HIb20 = HIB.HI_b_ST_20(z);
     double HIb30 = HIB.HI_b_ST_30(z);
     double Dz = C.D_growth_Tabulated(z);
     double pkm  = pow(Dz,2)*C.P_k_EH_Tabulated(k);

     double p13 =  HIb10*HIb10*P13_dd(k,z)  + HIb10*pkm*(HIb30 + 68.0*HIb20/21.0)*Sigma_out(z);


if(x == 1)
{
return pow(T_b(z),2)*pow(HIb10*Dz,2)*C.P_k_EH_Tabulated(k);
}
else if(x ==2) 	
{
///return pow(T_b(z),2)*(Non_LinearPower_HIHI(10,10, k,z)-Non_LinearPower_HIHI(2,2, 10e-3, z));//-ShotNoise(z);//*WindowFun(z, k);
return pow(T_b(z),2)*(Non_LinearPower_HIHI(10,10, k,z));//-ShotNoise(z);//*WindowFun(z, k);
}
else if(x ==3)
{
return pow(T_b(z),2)*p13;//*WindowFun(z, k);
}
else
{
double tmp  = pow(HIb10*Dz,2)*C.P_k_EH_Tabulated(k) + (Non_LinearPower_HIHI(10,10,  k,z)/2.0 + p13);//*WindowFun(z, k);
//return pow(T_b(z),2)*(tmp-Non_LinearPower_HIHI(2,2, 10e-3, z)/2.0);// -ShotNoise(z);
return pow(T_b(z),2)*tmp;
}

}

double LinearPower::HIPowerSpectrumR(int x, double k,double z,double R)
{
	
     Cosmology C;	
      HIBias HIB;
     // double upp = -2.0/(2.0 +C.n_s);
     // double R = 5.0*pow((1.0+z),upp);;///in Mpc/h
     double HIb10 = HIB.HI_b_ST_10(z); 
     double HIb20 = HIB.HI_b_ST_20(z);
     double HIb30 = HIB.HI_b_ST_30(z);
     double Dz = C.D_growth_Tabulated(z);
     double pkm  = pow(Dz,2)*C.P_k_EH_Tabulated(k);

     double p13 =  HIb10*HIb10*P13_dd(k,z)  + HIb10*pkm*(HIb30 + 68.0*HIb20/21.0)*MysigmaR(z,R);


if(x == 1)
{
return pow(T_b(z),2)*pow(HIb10*Dz,2)*C.P_k_EH_Tabulated(k);
}
else if(x ==2)
{
return pow(T_b(z),2)*Non_LinearPower_HIHI(10,10, k,z);//*WindowFun(z, k);
}
else if(x ==3)
{
return pow(T_b(z),2)*p13;//*WindowFun(z, k);
}
else
{
double tmp  = pow(HIb10*Dz,2)*C.P_k_EH_Tabulated(k) + (Non_LinearPower_HIHI(10,10,  k,z)/2.0 + p13);//*WindowFun(z, k);
return pow(T_b(z),2)*tmp;
}

}


double LinearPower::HIPowerSpectrum_NG(int x, double fnl, double k,double z) 
{

     Cosmology C;
      HIBias HIB;
      double upp = -2.0/(2.0 +C.n_s);
      double R = 5.0*pow((1.0+z),upp);;///in Mpc/h
     double HIb10 = HIB.HI_b_ST_10(z); 
     double HIb20 = HIB.HI_b_ST_20(z);
     double HIb30 = HIB.HI_b_ST_30(z);
     double Dz = C.D_growth_Tabulated(z);
     double pkm  = pow(Dz,2)*C.P_k_EH_Tabulated(k);

     double p13 =  b_one(fnl, k,z)*(HIb10*P13_dd(k,z)  + pkm*(HIb30 + 68.0*HIb20/21.0)*Sigma_out(z));


if(x == 1)
{
return pow(b_one(fnl, k,z)*Dz,2)*C.P_k_EH_Tabulated(k);
}
else if(x ==2)
{
return Non_LinearPower_HIHItidal(10,10, k,z)-Non_LinearPower_HIHI(2,2, 10e-3, z);//*WindowFun(z, k);
}
else if(x ==3)
{
return p13;//*WindowFun(z, k);
}
else
{
double tmp  = pow(b_one(fnl, k,z)*Dz,2)*C.P_k_EH_Tabulated(k) + (Non_LinearPower_HIHItidal(10,10,  k,z)/2.0 + p13);//*WindowFun(z, k);
return tmp-Non_LinearPower_HIHI(2,2, 10e-3, z)/2.0;
}

}


double LinearPower::MatterPowerSpectrum(int x, double k,double z) 
{
if(x == 1)
{
return LinearPower_mm(k,z);
}
else if(x == 2)
{
return PowerSpectrum22(k,z);//*WindowFun(z, k);
}
else if(x == 3)
{
return P13_dd(k,z);//*WindowFun(z, k);
}
else
{
double tmp  = LinearPower_mm(k,z) + (PowerSpectrum22(k,z)/2.0 + P13_dd(k,z));//*WindowFun(z, k);
return tmp;
}
   


}



double LinearPower::DeltaPowerSpectrum(int x, double k,double z)
{
double Pm = LinearPower_mm(k,z) + PowerSpectrum22(k,z) + P13_dd(k,z);
if(x == 1)
{
return (HIPowerSpectrum(10, k, z) - HIPowerSpectrum(1, k, z))/HIPowerSpectrum(1, k, z) ;
}
else if(x ==2)
{
return (MatterPowerSpectrum(10, k, z) - MatterPowerSpectrum(1, k, z))/MatterPowerSpectrum(1, k, z) ;
}
else 
return 0.0;
}


double LinearPower::beffrealspace(double fnl,double k,double z)
{
 double tmp = sqrt(HIPowerSpectrum_NG(10,fnl,  k, z) /MatterPowerSpectrum(10, k,z) );
return tmp;
}




/////////////////////////////////////////////////Redshift space ////////////////////////////////////////////////////

/*
 double FME_2( double k1, double k2, double muk)
{    double r = q/k;
 // double muk = mu*mu1+ sqrt(1.0-pow(mu1,2))*sqrt(1.0-pow(mu,2));
    double y = sqrt(r*r-2.0*muk *r + 1.0);
    double first_term = 5./7.0;
    double second_term  =  ((r*muk-r*r)/(r*y))*(r/y +y/r)/2.0; 
    double third_term =  2.0*pow((r*muk-r*r)/(r*y),2)/7.0;
    double tmp  = first_term+ second_term+third_term;
    return 2.0*tmp;
}
*/


 double FRSD_2( double k, double q, double muk)
{    double r = q/k;
 // double muk = mu*mu1+ sqrt(1.0-pow(mu1,2))*sqrt(1.0-pow(mu,2));
    double y = sqrt(r*r-2.0*muk *r + 1.0);
    double first_term = 5./7.0;
    double second_term  =  ((r*muk-r*r)/(r*y))*(r/y +y/r)/2.0; 
    double third_term =  2.0*pow((r*muk-r*r)/(r*y),2)/7.0;
    double tmp  = first_term+ second_term+third_term;
    return 2.0*tmp;
}


double GRSD_2( double k, double q, double muk)
{
     double r = q/k;
//double muk = mu*mu1+ sqrt(1.0-pow(mu1,2))*sqrt(1.0-pow(mu,2));
    double y = sqrt(r*r-2.0*muk *r + 1.0);
    double first_term = 3./7.0;
    double second_term  =  ((r*muk-r*r)/(r*y))*(r/y +y/r)/2.0; 
    double third_term =  4.0*pow((r*muk-r*r)/(r*y),2)/7.0;
    double tmp  = first_term+ second_term+third_term;
    return 2.0*tmp;
}
double KNRSD(double fnl,double z,double k, double k1, double mu, double muk)
{

    HIBias HIB;
    double HIb10 = HIB.HI_b_ST_10(z);
      double r = k1/k;
     double y = sqrt(r*r-2.0*muk*r + 1.0);
    double  mu1 = muk *mu + sqrt((1.0-muk*muk)*(1.0-mu*mu));
     double mu2 = (mu- r*mu1)/(y);

    double biasr =  b_one(fnl, k*r, z);
  double biasy =  b_one(fnl, k*y, z);//


    double tmp1 = biasy*mu1*mu1 + biasr*mu2*mu2 + mu1*mu2*(biasr*(r/y) + biasy*(y/r));
    double tmp2 = 2.0*mu1*mu1*mu2*mu2  + mu1*mu2*((mu1*mu1*r/y) + (mu2*mu2*y/r));

   double tmp = f2(z)*tmp1 + pow(f2(z),2)*tmp2;

  return tmp;

  // double tmptest = 2.0*mu1*mu1*mu2*mu2*pow(f2(z),2);

   // return tmptest;
}


double KNRSDmatter(double fnl,double z,double k, double k1, double mu, double muk)
{

    
      double r = k1/k;
     double y = sqrt(r*r-2.0*muk *r + 1.0);
    double  mu1 = muk *mu + sqrt((-1.0+muk*muk)*(-1.0+mu*mu));
     double mu2 = (mu- r*mu1)/(y);

  //  double biasr =  b_one(fnl, k*r, z);
  //double biasy =  b_one(fnl, k*y, z);//


    double tmp1 = mu1*mu1 + mu2*mu2 + mu1*mu2*((r/y) + (y/r));
    double tmp2 = 2.0*mu1*mu1*mu2*mu2  + mu1*mu2*((mu1*mu1*r/y) + (mu2*mu2*y/r));

   double tmp = f2(z)*tmp1 + pow(f2(z),2)*tmp2;

  return tmp;

  // double tmptest = 2.0*mu1*mu1*mu2*mu2*pow(f2(z),2);

   //s return tmptest;
}


double fnlterm(double fnl,double z,double k1, double k2,double k3)
{
  
     return 2.0*fnl*Alpha(k3,z)/(Alpha(k1,z)*Alpha(k2,z));//AlphaAlpha2(z,k1,k2);   
}


static double f22_dt(double k, double q, double x)
 {

     Cosmology C;
    double r = q/k;
    double d = 1 + r*r - 2*r*x;
    if(d < 1e-5)
        return 0;
    else
        return C.P_k_EH_Tabulated(q) * C.P_k_EH_Tabulated(k*sqrt(d)) * (3*r + 7*x - 10*r*x*x)*(7*x - r - 6*r*x*x) / pow2(d);
}

/* P_{\delta\theta}^{(22)} */
double P22_dt(double k, double z) 
 {
    Cosmology C;  
   const double EPSREL = 1e-5;
     const double EPSABS = 1e-5;
    double Dz = C.D_growth_Tabulated(z);
    double a[2] = { QMIN, -1 };
    double b[2] = { QMAX, XMAX };
    return pow(Dz,4)*k*k/(98*4*M_PI*M_PI) * Integrate<2>(bind(f22_dt, k, _1, _2), a, b, EPSREL, 1e-4*C.P_k_EH_Tabulated(k));
}


static double f22_tt(double k,double q,double x)
 {
    Cosmology C;
    double r = q/k;
    double d = 1 + r*r - 2*r*x;
    if(d < 1e-5)
        return 0;
    else
        return C.P_k_EH_Tabulated(q) * C.P_k_EH_Tabulated(k*sqrt(d)) * pow2(7*x - r - 6*r*x*x) / pow2(d);
}

/* P_{\theta\theta}^{(22)} */
double P22_tt(double k, double z)  
{
    Cosmology C;  
   const double EPSREL = 1e-5;
     const double EPSABS = 1e-5;
     double Dz = C.D_growth_Tabulated(z);
    double a[2] = { QMIN, -1 };
    double b[2] = { QMAX, XMAX };
    return pow(Dz,4)*k*k/(98*4*M_PI*M_PI) * Integrate<2>(bind(f22_tt, k, _1, _2), a, b, EPSREL, 1e-4*C.P_k_EH_Tabulated(k));
}

double LinearPower::P22_ttout(double k, double z)  
{
return P22_tt(k,z);
}


double  f13_tt(double k, double q) 
{
    Cosmology C;
    double r = q/k;
    double s;
    if(r < 1e-2)
        s = -56 - (32./5.)*pow2(r) - (96./7.)*pow4(r) + (352./105.)*pow6(r);
    else if(fabs(r-1) < 1e-10)
        s = -72 - 40*(r-1);
    else if(r > 100)
        s = -504./5. + (1248./35.)/pow2(r) - (608./105.)/pow4(r) - (160./231.)/pow6(r);
    else
        s = 12/pow2(r) - 82 + 4*pow2(r) - 6*pow4(r) + 3/pow3(r) * pow3(r*r - 1) * (r*r + 2) * log((1+r)/fabs(1-r));

    return C.P_k_EH_Tabulated(q) * s;
}

/* P_{\theta\theta}^{(13)} */
    double P13_tt(double  k,double z)
 {
   Cosmology C;  
   const double EPSREL = 1e-5;
     const double EPSABS = 1e-5;
     double Dz = C.D_growth_Tabulated(z);
    return pow4(Dz)*pow2(k)/(84*4*M_PI*M_PI) * C.P_k_EH_Tabulated(k) * Integrate(bind(f13_tt, k, _1), QMIN, QMAX, EPSREL, 1e-4*C.P_k_EH_Tabulated(k));
}

double f13_dt(double k, double q) 
{
   Cosmology C;  
   double r = q/k;
    double s;
    if(r < 1e-2)
        s = -168 + (416./5.)*pow2(r) - (2976./35.)*pow4(r) + (224./15.)*pow6(r);
    else if(fabs(r-1) < 1e-10)
        s = -152 - 56*(r-1);
    else if(r > 100)
        s = -200 + (2208./35.)/pow2(r) - (1312./105.)/pow4(r) - (1888./1155.)/pow6(r);
    else
        s = 24/pow2(r) - 202 + 56*pow2(r) - 30*pow4(r) + 3/pow3(r) * pow3(r*r - 1) * (5*r*r + 4) * log((1+r)/fabs(1-r));

    return C.P_k_EH_Tabulated(q) * s;
}

/* P_{\delta\theta}^{(13)} */
double P13_dt(double k, double z) 
 {
   Cosmology C;  
   const double EPSREL = 1e-5;
     const double EPSABS = 1e-5;
     double Dz = C.D_growth_Tabulated(z);
    return pow4(Dz)*k*k/(252*4*M_PI*M_PI) * C.P_k_EH_Tabulated(k) * Integrate<ExpSub>(bind(f13_dt, k, _1), QMIN, QMAX, EPSREL, 1e-4*C.P_k_EH_Tabulated(k));
}



/* P_{ab}(k) */
double LinearPower::AllPk(int x, double  k, double z) 
 {

        Cosmology C;
        double Dz = C.D_growth_Tabulated(z);
      double power_multi= Dz*Dz*C.P_k_EH_Tabulated(k);
       double Dz4 =  pow(Dz,4);
    if(x ==1) 
{
      
            return power_multi + P13_dd(k,z) + PowerSpectrum22(k,z)/2.0;
}
else if(x ==2)
{
      
            return power_multi + P13_dt(k,z) + P22_dt(k,z);
}
else if(x == 3)
{
            return power_multi + P13_tt(k,z) + P22_tt(k,z);
}
else 
  {    
            cout<<"\t"<<"invalid choice , choose x < 4: 1= matter-matter, 2= matter-velocity, 3=velocity-velocity\n"<<endl;
            return 0;
    }
}

double LinearPower::ScoccimarroPk(double  k, double z, double mu) 
{
       HIBias HIB;
     double HIb10 = HIB.HI_b_ST_10(z); 
    double tmp =  pow(HIb10,2)*AllPk(1,  k,z)  + 2.0* f2(z)*mu*mu* HIb10*AllPk(2,  k,z) +  pow(mu*mu*f2(z),2)*AllPk(3,  k,z);
   return pow(T_b(z),2)*tmp;

}

double LinearPower::ScoccimarroPkAverage(double  k, double z) 
{
       HIBias HIB;
     double HIb10 = HIB.HI_b_ST_10(z); 
    double tmp =  pow(HIb10,2)*AllPk(1,  k,z)  +  f2(z)* HIb10*AllPk(2,  k,z)/3.0 +  pow(f2(z),2)*AllPk(3,k,z);
   return pow(T_b(z),2)*tmp;

}


	
double RenormSecondOrderKernelPRSD(int x,double fnl, double z, double k,double k1,double mu,double muk)
{    
       HIBias HIB;
     double HIb10 = HIB.HI_b_ST_10(z); 
     double HIb20 = HIB.HI_b_ST_20(z); 
double  r = k1/k;
     double y = sqrt(r*r-2.0*muk *r + 1.0);
     double  mu1 = muk *mu + sqrt((1.0-muk*muk)*(1.0-mu*mu));
     double mu2 = (k*mu- k1*mu1)/(k*y);
   double bias2 = HIb20;//b_two_integrand(fnl, z ,k1, k*y, mu);
if(x == 1)
{
return HIb10*F2P(k,k*r,muk);//FRSD_2(k, k1,muk) ;//GRF2P(fnl,z, k, k1, mu);
}
else if (x == 2)
{
return GRSD_2(k, k1, muk)*mu*mu*f2(z);
}
else if (x == 3)
{
return HIb10*fnlterm(fnl,z,k1, k*y,k);//HIb10*GRF2P_mat(fnl, z, k,  k1, muk);//
}
else if(x == 4)
{
return HIb20;//b_two_integrand(fnl, z ,k1, k*y);
}
else if(x == 5)
{
return KNRSD(fnl, z, k,  k1, mu, muk); 
}
else
{
double tmp  = HIb10*F2P(k,k*r,muk)  + HIb20 + pow(mu,2)*f2(z)*GRSD_2(k, k1, muk) + KNRSD(fnl, z, k,  k1,  mu, muk);//+ HIb10*fnlterm(fnl,z,k1, k*y,k)
return tmp;
}
}

double RenormPowerSpectrumIntegrandRSD(int x1,int x2,double fnl, double z, double k,double mu,double k1,double muk)
    {
     Cosmology C;
      double  r =  k1/k;
      double y = sqrt(r*r-2.0*muk *r + 1.0);
      double Dz = C.D_growth_Tabulated(z);
      double power_multi= pow(Dz,4)*C.P_k_EH_Tabulated(k*r)*C.P_k_EH_Tabulated(k*y);
      double coeff =  k1*k1/(4.0*pow(M_PI,2));
      double tmp2 =  RenormSecondOrderKernelPRSD(x1, fnl, z, k,k*r,mu, muk)*RenormSecondOrderKernelPRSD(x2, fnl, z, k,k*r,mu, muk); 
     // double tmp = coeff*SecondOrderKernelP(x1,z, k,k*r,muk)*SecondOrderKernelP(x2,z, k,k*r,muk)*power_multi;
      double tmp3 = coeff*tmp2*power_multi;
      //cout<<"\t"<<"Integrating="<<"\t"<<tmp3<<endl;
       return tmp3;
}

double FoG(double k, double z, double mu)
{
double tmp = pow(f2(z),2)*k*k*pow(mu,2)*Sigmav_out(z)/2.0;
return exp(-tmp);

}

 double RenormNon_LinearPower_HIHIRSD22(int x1,int x2,double fnl, double k,double  z,double mu)
{
      Cosmology C;
      const double EPSREL = 1e-3;
     const double EPSABS = 1e-3;
      //const double QMIN = 1e-4;
    //  double upp = 2.0/(2.0 +C.n_s);
     // double knl = 0.2*pow((1.0+z),upp);///in h/Mpc//1.0/Sigmav_out(z);//
      // double QMAX = knl;
      // const double QMAX = 1e4;
     // const double XMAX = 0.98;
     double a[2] = {QMIN,-1.0};
    double  b[2] = {QMAX, XMAX};
    //double fnlin  = fnl -5.0/3.0; 
// const double kin = 71.94;

//RenormPowerSpectrumIntegrandRSD(int x1,int x2,double fnl, double z, double k,double mu,double r,double muk)
 double tmpfinal = Integrate<2>(bind(RenormPowerSpectrumIntegrandRSD,x1,x2,fnl, z,k,mu,_1,_2),a, b,EPSREL,EPSABS);//LOS_tab(z);/
 ///cout<<"\t"<<"Integrating="<<"\t"<<tmpfinal <<endl;
   return pow(T_b(z),2)*tmpfinal;
  // return tmpfinal;
    
}

double LinearPower::RenormNon_LinearPower_HIHIRSD22out(int x1,int x2,double fnl, double k,double  z,double mu)
{
return RenormNon_LinearPower_HIHIRSD22(x1,x2,fnl,k,z, mu)/2.0;
}




//////////////////////////////////////////////////////////////////////matter/////////////////////////////////


double RenormSecondOrderKernelPRSDmatter(int x,double fnl, double z, double k,double k1,double mu,double muk)
{    
     
double  r = k1/k;
     double y = sqrt(r*r-2.0*muk *r + 1.0);
     double  mu1 = muk *mu + sqrt((1.0-muk*muk)*(1.0-mu*mu));
     double mu2 = (k*mu- k1*mu1)/(k*y);
   double bias2 = 0;//b_two_integrand(fnl, z ,k1, k*y, mu);
if(x == 1)
{
return FRSD_2(k, k1,muk) ;//GRF2P(fnl,z, k, k1, mu);
}
else if (x == 2)
{
return mu*mu*f2(z)*GRSD_2(k, k1, muk);
}
else if (x == 3)
{
return fnlterm(fnl,z,k1, k*y,k);//HIb10*GRF2P_mat(fnl, z, k,  k1, muk);//
}
else if(x == 4)
{
return 0.0;//b_two_integrand(fnl, z ,k1, k*y);
}
else if(x == 5)
{
return KNRSDmatter(fnl, z, k,  k1, mu, muk); 
}
else
{
double tmp  = FRSD_2(k, k1,muk)+pow(mu,2)*f2(z)*GRSD_2(k, k1, muk)  + KNRSDmatter(fnl, z, k,  k1,  mu, muk);//+ HIb10*fnlterm(fnl,z,k1, k*y,k)
return tmp;
}
}


double RenormPowerSpectrumIntegrandRSDmatter(int x1,int x2,double fnl, double z, double k,double mu,double k1,double muk)
    {
     Cosmology C;
      double  r =  k1/k;
      double y = sqrt(r*r-2.0*muk *r + 1.0);
      double Dz = C.D_growth_Tabulated(z);
      double power_multi= Dz*Dz*Dz*Dz*C.P_k_EH_Tabulated(k*r)*C.P_k_EH_Tabulated(k*y);
      double coeff =  k1*k1/(4.0*pow(M_PI,2));
      double tmp2 =  RenormSecondOrderKernelPRSDmatter(x1, fnl, z,  k,k1,mu, muk)*RenormSecondOrderKernelPRSDmatter(x2, fnl, z,  k,k1,mu, muk); 
      double tmp3 = coeff*tmp2*power_multi;
      //cout<<"\t"<<"Integrating="<<"\t"<<tmp3<<endl;
       return tmp3;
}



 double RenormNon_LinearPower_HIHIRSD22matter(int x1,int x2,double fnl, double k,double  z,double mu)
{
      Cosmology C;
      const double EPSREL = 1e-2;
     const double EPSABS = 1e-2;
      const double QMIN = 1e-4;
      double upp = 2.0/(2.0 +C.n_s);
      double knl = 0.2*pow((1.0+z),upp);///in h/Mpc//1.0/Sigmav_out(z);//
      // double QMAX = knl;
       double QMAX = 1e4;
      const double XMAX = 0.98;
     double a[2] = {QMIN,-1.0};
    double  b[2] = {QMAX, XMAX};
    //double fnlin  = fnl -5.0/3.0; 
// const double kin = 71.94;


 double tmpfinal = Integrate<2>(bind(RenormPowerSpectrumIntegrandRSDmatter,x1,x2,fnl, z,k,mu,_1,_2),a, b,EPSREL,EPSABS);//LOS_tab(z);/
 ///cout<<"\t"<<"Integrating="<<"\t"<<tmpfinal <<endl;
   return pow(T_b(z),2)*tmpfinal;
    
}















double ThirdorderRSDintegrdmatter(double z, double k,double mu, double k1)
{
    // HIBias HIB;
     Cosmology C;
      double r = k1/k;
    ///  double y = sqrt(r*r-2.0*muk *r + 1.0);
      double HIb10 = 1.0;//HIB.HI_b_ST_10(z); 
      double HIb20 = 0.0;//HIB.HI_b_ST_20(z); 
     // double  mu1 = muk *mu + sqrt((1.0-muk*muk)*(1.0-mu*mu));
      double Dz = C.D_growth_Tabulated(z);
      double power_multi= Dz*Dz*C.P_k_EH_Tabulated(k*r);
      double coeff =  k1*k1/(4.0*pow(M_PI,2));
/*
      double firstterm = -12.0*pow(f2(z),3)*pow4(mu1)*pow2(mu)-6.0*HIb10*pow2(mu1)*pow2(f2(z))*(pow2(mu1) + pow2(mu)) + 6.0* pow2(mu1)*pow2(mu)*f2(z)*(2.0*f2(z) - 1.0);

      double secondterm = HIb10*f2(z)*F2Simp(k, k1, muk)*((pow2(mu)/r) - mu1*mu + 2.0*pow2(mu1)+pow2(mu)) + HIb20*f2(z)*((pow2(mu)/r) - mu1*mu + 6.0*pow2(mu1)+3.0*pow2(mu));

// double secondterm = HIb10*f2(z)*FF2(k1, -k, muk)*((pow2(mu)/r) - mu1*mu - 2.0*pow2(mu1)-pow2(mu));// + HIb20*f2(z)*((pow2(mu)/r) - mu1*mu - 6.0*pow2(mu1)-3.0*pow2(mu));

     double thirdterm = (HIb10 + pow2(mu)*f2(z))*pow2(f2(z))*(2.0*pow2(mu1)*(3.0*pow2(mu1)-pow2(mu))-pow2(mu1)*pow2(mu));
 
    double fourtterm1 = f2(z)*(HIb10 + pow2(mu)*f2(z))*(pow2(mu) - mu1*mu*r) ;

    double fourtterm2 = -3.0*HIb10*f2(z)*(2.0*mu1*mu*r - pow2(mu) - pow2(mu1) *pow2(r) ) ;

    double fourtterm31 = -mu1*pow3(mu) + 4.0*pow2(mu)*pow2(mu1) - 2.0* (mu1*mu + pow2(mu))*pow2(mu);

    double fourtterm32 =  pow4(mu)/r + r*pow2(mu1)*(8.0*mu1*mu + pow2(mu))+ 2.0*r*pow2(mu)*(2.0*mu1*mu + pow2(mu1));

    double fourtterm33 = -pow2(r)*(2.0*pow2(mu1)*pow2(mu) + pow2(mu1)*(mu1*mu + 4.0* pow2(mu1)) ) ;

    double fourtterm3 =  pow2(f2(z))*(fourtterm31 + fourtterm32 + fourtterm33);

    double fourtterm = G2Simp(k, k1, muk)*(fourtterm1 + fourtterm2 + fourtterm3)/pow2(y);

      double tmp = firstterm + secondterm  + thirdterm + fourtterm;

      return power_multi *coeff* tmp/3.0;
*/

double Xt1 = 9.0* pow(r,4)- 24.0*pow(r,2) + 19.0;
double Xt2 = -9.0* pow(r,7) + 33.0 * pow(r,5)  + 33.0* pow(r,3) - 9.0*r;
double Xt3 = -27.0*pow(r,6)+ 63.0*pow(r,4)  - 109.0* pow(r,2) + 9.0;
double logterm =log((1.0+r)/fabs(1.0-r));//log((1+r)/fabs(1-r)

double B1101 = 1.0/2.0;
double B1110 = (-2.0*(Xt1) + 9.0*(pow((r*r -1.0),3)/r)*logterm )/84.0;
double B1210 = -1.0/3.0;
double B1200 = -(2.0*(Xt2) + 9.0*pow((r*r -1.0),4) *logterm)/(336.0*pow(r,3));
double B2220 = (2.0*r*(Xt3) + 9.0*(3.0*r*r + 1.0)*pow((r*r -1.0),3) *logterm)/(336.0*pow(r,3));
double B2300 = -1.0/3.0;

double tmp1 = mu*mu*f2(z)*(HIb10*B1110+HIb20*B1110)+ pow(mu,4)* f2(z)*f2(z)*(HIb10*HIb10*B2220 + f2(z)*B2300) + pow(mu*f2(z),2)*(HIb10*B1210 + B1200);

 double tmp  =  power_multi*coeff*tmp1;
  return tmp;

}

 double CurlyItermmatter(double k,double  z, double mu)
{
    Cosmology C;  
    const double EPSREL = 1e-2;
     const double EPSABS = 1e-2;
      const double QMIN = 1e-4;
     double upp = 2.0/(2.0 +C.n_s);
      double knl = 0.2*pow((1.0+z),upp);///in h/Mpc//1.0/Sigmav_out(z);//
       double QMAX = knl;
      const double XMAX = 0.98;
    // double a[2] = {QMIN,-1.0};
    //double  b[2] = {QMAX, XMAX};
 double tmpfinal = Integrate(bind(ThirdorderRSDintegrdmatter,z,k,mu,_1),QMIN,QMAX,EPSREL,EPSABS);
//cout<<"\t"<<"Curly term="<<"\t"<<tmpfinal <<endl;
   return tmpfinal;
    
}

double RenormNon_LinearPower_HIHI13matter(double fnl, double k,double  z,double mu)
{    
      Cosmology C;
      //HIBias HIB;
     double HIb10 = 1.0;//HIB.HI_b_ST_10(z); 
     double HIb20 = 0.0;//HIB.HI_b_ST_20(z); 
     double HIb30 = 0.0;//HIB.HI_b_ST_30(z);
    // double upp = -2.0/(2.0 +C.n_s);
    // double R = 5.0*pow((1.0+z),upp);///in Mpc/h
     double Dz = C.D_growth_Tabulated(z);
     double pkm  = pow(Dz,2)*C.P_k_EH_Tabulated(k);
     double firstterm =  EvokernelZ1matter(fnl, k, z,mu)*pkm*(HIb30 + 68.0*HIb20/21.0)*Sigma_out(z);//pow(Dz,2)**C.sigma_R(R)/2.0;//;
     double secondterm =  EvokernelZ1matter(fnl, k, z,mu)*(HIb10*P13_ddint(k,z)  + pow2(mu)*f2(z)*P13_tt(k,z));//P13_dt(k,z)P13_tt(
     double thirdterm = EvokernelZ1matter(fnl, k, z,mu)*pkm*CurlyItermmatter(k, z,mu);
     double tmp = firstterm + secondterm + thirdterm;
  return pow(T_b(z),2)*tmp;
}



double LinearPower::TPowerSpectrummatter(int x,double fnl, double k,double z,double mu) 
{
//WindowFun(z, k)*
if(x == 1)
{
return LinearPower_HIHIRSDmatter(fnl,k,z,mu);//*FoG(k,z,mu);
}
else if(x ==2)
{
return RenormNon_LinearPower_HIHIRSD22matter(10,10,fnl, k,z,mu);//*FoG(k,z,mu);//WindowFun(z, k)*
}
else if(x ==3)
{
return RenormNon_LinearPower_HIHI13matter(fnl, k, z,mu);///*FoG(k,z,mu);//WindowFun(z, k)*
}
else
{
double tmp  = LinearPower_HIHIRSDmatter(fnl,k,z,mu) + (RenormNon_LinearPower_HIHIRSD22matter(10,10,fnl, k,z,mu)/2.0 + RenormNon_LinearPower_HIHI13matter(fnl, k, z,mu));
return tmp;//*FoG(k,z,mu);
}

}




double  LinearPower::TPowerSpectrum_Averagematter(int x, double fnl, double k, double z)
{
 Cosmology C;
 const double EPSREL = 1e-3;
  const double EPSABS = 1e-3;
const double mu_min = -0.998;
const double mu_max = 0.998;
double tmp = Integrate(bind(&LinearPower::TPowerSpectrummatter,this,x,fnl,k,z,_1),mu_min,mu_max,EPSREL,EPSABS)/2.0;
return tmp;
}



double LinearPower::beffHImatter(int x, double fnl, double k,double z)
{
if(x == 1)//LinearHIHI_AverageRSDmatter(int l, double fnl, double k, double z)
   ///LinearHIHI_AverageRSD( int l, double fnl, double k, double z)
{
double tmp = sqrt(LinearHIHI_AverageRSD(0,fnl, k, z)/LinearHIHI_AverageRSDmatter(0,fnl, k, z));
return tmp;
}
else
{
 double tmp = sqrt(TPowerSpectrum_Average(10,fnl, k, z)/TPowerSpectrum_Averagematter(10,fnl, k, z));
return tmp;
}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double k1dotk2(double k1, double k2, double mu)
{
    return k1*k2*mu;
}

double FF2(double k1, double k2, double mu12)
{
   double first_term = 5.0/7.0;
   double second_term = (k1/k2 + k2/k1)*k1dotk2(k1, k2, mu12)/(2.0*k1*k2);
   double third_term = (2.0/7.0)*pow((k1dotk2(k1, k2,mu12)/(k1*k2)),2);
   double tmp =  first_term+second_term+third_term;
   return 2.0*tmp;
}
double GG2(double k1, double k2, double mu12)
{
   double first_term = 3.0/7.0;
   double second_term = (k1/k2 + k2/k1)*k1dotk2(k1, k2, mu12)/(2.0*k1*k2);
   double third_term = (4.0/7.0)*pow((k1dotk2(k1, k2,mu12)/(k1*k2)),2);
   double tmp =  first_term+second_term+third_term;
   return 2.0*tmp;
}


double F2Simp(double k, double k1, double mu1)
{
   double r = k1/k;
   double first_term = 5.0/7.0;
   double second_term = -(r + 1.0/r)*(mu1)/(2.0);
   double third_term = (2.0/7.0)*pow(mu1,2);
   double tmp =  first_term+second_term+third_term;
   return 2.0*tmp;
}

double G2Simp(double k, double k1, double mu1)
{
     double r = k1/k;
    double first_term = 3.0/7.0;
    double second_term  =  -(r + 1.0/r)*(mu1)/(2.0); 
    double third_term =  4.0*pow(mu1,2)/7.0;
    double tmp  = first_term+ second_term+third_term;
    return 2.0*tmp;
}

double ThirdorderRSDintegrd2(double z, double k,double mu, double k1, double muk)
{
     HIBias HIB;
     Cosmology C;
      double r = k1/k;
      double  mu1 = mu1 *muk + sqrt((1.0-mu*mu)*(1.0-muk*muk));
      double y = sqrt(r*r-2.0*muk *r + 1.0);
      double HIb10 = HIB.HI_b_ST_10(z); 
      double HIb20 = HIB.HI_b_ST_20(z); 
      
      double Dz = C.D_growth_Tabulated(z);
      double power_multi= Dz*Dz*C.P_k_EH_Tabulated(k*r);
      double coeff =  k1*k1/(4.0*pow(M_PI,2));
      double firstterm =  mu*mu*f2(z)*HIb20*k*(1.0/k+2.0/k1);
      double secondterm = 2.0*mu*mu*f2(z)*k*(HIb10*F2Simp(k, k1, muk) + f2(z)*G2Simp(k, k1, muk)*pow2(mu1*k1-mu*k)/pow2(k*y))/k1;
     double thirdterm = -2.0*mu*f2(z)*k*((HIb10 + pow(mu,2)*f2(z)) + 2.0*(HIb10 + pow2(mu1)*f2(z)))*G2Simp(k, k1, muk)*(mu1*k1-mu*k)/pow2(k*y);
    double fourtterm = pow(mu*f2(z)*k,2)*(HIb10 + pow(mu,2)*f2(z))*pow(mu1,2)*(1.0/2.0 + k1/k)/pow2(k1);
      double tmp = firstterm + secondterm  + thirdterm - fourtterm;
      return power_multi *coeff* tmp;



/*
double Xt1 = 9.0* pow(r,4)- 24.0*pow(r,2) + 19.0;
double Xt2 = -9.0* pow(r,7) + 33.0 * pow(r,5)  + 33.0* pow(r,3) - 9.0*r;
double Xt3 = -27.0*pow(r,6)+ 63.0*pow(r,4)  - 109.0* pow(r,2) + 9.0;
double logterm =log((1.0+r)/fabs(1.0-r));//log((1+r)/fabs(1-r)

double B1101 = 1.0/6.0;
double B1110 = (-2.0*(Xt1) + 9.0*(pow((r*r -1.0),3)/r)*logterm )/84.0;
double B1210 = -1.0/3.0;
double B1200 = -(2.0*(Xt2) + 9.0*pow((r*r -1.0),4) *logterm)/(336.0*pow(r,3));
double B2220 = (2.0*r*(Xt3) + 9.0*(3.0*r*r + 1.0)*pow((r*r -1.0),3) *logterm)/(336.0*pow(r,3));
double B2300 = -1.0/3.0;

double tmp1 = mu*mu*f2(z)*(HIb10*B1110+HIb20*B1110)+ pow(mu,4)* f2(z)*f2(z)*(HIb10*HIb10*B2220 + f2(z)*B2300) + pow(mu*f2(z),2)*(HIb10*B1210 + B1200);

 double tmp  =  power_multi*coeff*tmp1;
  return tmp;
*/

}


 double CurlyIterm2(double k,double  z, double mu)
{
    Cosmology C;  
    const double EPSREL = 1e-2;
     const double EPSABS = 1e-2;
      const double QMIN = 1e-4;
     double upp = 2.0/(2.0 +C.n_s);
      double knl = 0.2*pow((1.0+z),upp);///in h/Mpc//1.0/Sigmav_out(z);//
       double QMAX = knl;
       //const double QMAX = 1e4;
      const double XMAX = 0.998;
     double a[2] = {QMIN,-1.0};
    double  b[2] = {QMAX, XMAX};
 double tmpfinal = Integrate<2>(bind(ThirdorderRSDintegrd2,z,k,mu,_1,_2),a, b,EPSREL,EPSABS);
///cout<<"\t"<<"Curly term="<<"\t"<<tmpfinal <<endl;
   return tmpfinal;
    
}

double ThirdorderRSDintegrd1(double z, double k,double mu, double k1)
{
     HIBias HIB;
     Cosmology C;
      double r = k1/k;
    ///  double y = sqrt(r*r-2.0*muk *r + 1.0);
      double HIb10 = HIB.HI_b_ST_10(z); 
      double HIb20 = HIB.HI_b_ST_20(z); 
     // double  mu1 = muk *mu + sqrt((1.0-muk*muk)*(1.0-mu*mu));
      double Dz = C.D_growth_Tabulated(z);
      double power_multi= Dz*Dz*C.P_k_EH_Tabulated(k*r);
      double coeff =  k1*k1/(4.0*pow(M_PI,2));


   
double Xt1 = 9.0* pow(r,4)- 24.0*pow(r,2) + 19.0;
double Xt2 = -9.0* pow(r,7) + 33.0 * pow(r,5)  + 33.0* pow(r,3) - 9.0*r;
double Xt3 = -27.0*pow(r,6)+ 63.0*pow(r,4)  - 109.0* pow(r,2) + 9.0;
double logterm =log((1.0+r)/fabs(1.0-r));//log((1+r)/fabs(1-r)

double B1101 = 1.0/2.0;
double B1110 = (-2.0*(Xt1) + 9.0*(pow((r*r -1.0),3)/r)*logterm )/84.0;
double B1210 = -1.0/3.0;
double B1200 = -(2.0*(Xt2) + 9.0*pow((r*r -1.0),4) *logterm)/(336.0*pow(r,3));
double B2220 = (2.0*r*(Xt3) + 9.0*(3.0*r*r + 1.0)*pow((r*r -1.0),3) *logterm)/(336.0*pow(r,3));
double B2300 = -1.0/3.0;

double tmp1 = mu*mu*f2(z)*(HIb10*B1110+HIb20*B1110)+ pow(mu,4)* f2(z)*f2(z)*(HIb10*HIb10*B2220 + f2(z)*B2300) + pow(mu*f2(z),2)*(HIb10*B1210 + B1200);

 double tmp  =  power_multi*coeff*tmp1;
  return tmp;

}

 double CurlyIterm1(double k,double  z, double mu)
{
    Cosmology C;  
    const double EPSREL = 1e-2;
     const double EPSABS = 1e-2;
      const double QMIN = 1e-4;
     double upp = 2.0/(2.0 + C.n_s);
      double knl = 0.2*pow((1.0+z),upp);///in h/Mpc//1.0/Sigmav_out(z);//
       double QMAX = knl;
      // const double QMAX = 1e4;
      const double XMAX = 0.98;
   //  double a[2] = {QMIN,-1.0};
   // double  b[2] = {QMAX, XMAX};
 double tmpfinal = Integrate(bind(ThirdorderRSDintegrd1,z,k,mu,_1),QMIN, QMAX,EPSREL,EPSABS*C.P_k_EH_Tabulated(k));
///cout<<"\t"<<"Curly term="<<"\t"<<tmpfinal <<endl;
   return tmpfinal;
    
}

double LinearPower::ExportCurlyIterm1( double k,double z,double mu) 
{
return CurlyIterm1(k,z,mu);
}

double LinearPower::ExportCurlyIterm2( double k,double z,double mu) 
{
return CurlyIterm2(k,z,mu);
}

double RenormNon_LinearPower_HIHI13(double fnl, double k,double  z,double mu)
{    
      Cosmology C;
      HIBias HIB;
     double HIb10 = HIB.HI_b_ST_10(z); 
     double HIb20 = HIB.HI_b_ST_20(z); 
     double HIb30 = HIB.HI_b_ST_30(z);
    // double upp = -2.0/(2.0 +C.n_s);
    // double R = 5.0*pow((1.0+z),upp);///in Mpc/h double p13 =  HIb10*HIb10*P13_dd(k,z)  + HIb10*pkm*(HIb30 + 68.0*HIb20/21.0)*Sigma_out(z);
     double Dz = C.D_growth_Tabulated(z);
     double pkm  = pow(Dz,2)*C.P_k_EH_Tabulated(k);
     double firstterm =  pkm*(HIb30 + 68.0*HIb20/21.0)*Sigma_out(z);//pow(Dz,2)**C.sigma_R(R)/2.0;//;
     double secondterm =  (HIb10*P13_ddint(k,z)  + pow2(mu)*f2(z)*P13_tt(k,z));//P13_dt(k,z)P13_tt(
     double thirdterm = pkm*CurlyIterm1(k, z,mu);//EvokernelZ1(fnl, k, z,mu)//+mu*CurlyB(1.0, z)*(CH2(z)/k)
     double tmp = (EvokernelZ1(fnl, k, z,mu))*(firstterm + secondterm + thirdterm);
  return pow(T_b(z),2)*tmp;
}


double LinearPower::TPowerSpectrum(int x,double fnl, double k,double z,double mu) 
{
//WindowFun(z, k)*
if(x == 1)
{
return LinearPower_HIHIRSD(fnl,k,z,mu);//*FoG(k,z,mu);
}
else if(x ==2)
{
return RenormNon_LinearPower_HIHIRSD22(10,10,fnl, k,z,mu)-RenormNon_LinearPower_HIHIRSD22(4,4,fnl, 1e-3,z,mu);//*FoG(k,z,mu);//WindowFun(z, k)*
}
else if(x ==3)
{
return RenormNon_LinearPower_HIHI13(fnl, k, z,mu);///*FoG(k,z,mu);//WindowFun(z, k)*
}
else
{
double tmp  = LinearPower_HIHIRSD(fnl,k,z,mu) + ((RenormNon_LinearPower_HIHIRSD22(10,10,fnl, k,z,mu)-RenormNon_LinearPower_HIHIRSD22(4,4,fnl, 1e-3,z,mu))/2.0 + RenormNon_LinearPower_HIHI13(fnl, k, z,mu));
return tmp;//*FoG(k,z,mu);
}

}




double  LinearPower::TPowerSpectrum_Average(int x, double fnl, double k, double z)
{
 Cosmology C;
 const double EPSREL = 1e-2;
  const double EPSABS = 1e-2;
const double mu_min = -0.998;
const double mu_max = 0.998;
double tmp = Integrate(bind(&LinearPower::TPowerSpectrum,this,x,fnl,k,z,_1),mu_min,mu_max,EPSREL,EPSABS)/2.0;
return tmp;
}

double LinearPower::beffredshift(double fnl, double k,double z)
{
 double tmp = sqrt(TPowerSpectrum_Average(10,fnl, k, z)/(pow(T_b(z),2)*MatterPowerSpectrum(10, k,z)));
return tmp;
}

double LinearPower::beffredshift2(double fnl, double k,double z,double mu)
{
 double tmp = sqrt(TPowerSpectrum(10,fnl, k, z,mu)/(pow(T_b(z),2)*MatterPowerSpectrum(10, k,z)));
return tmp;
}


double LinearPower::beffrealspace2(int x, double z)
{
    HIBias HIB;
     double HIb10 = HIB.HI_b_ST_10(z); 
     double HIb20 = HIB.HI_b_ST_20(z); 
     double HIb30 = HIB.HI_b_ST_30(z);
if(x == 1)
{
return HIb10;
}
else if(x ==2)
{
return HIb20;
}
else if(x ==3)
{
return HIb30;
}
else
{
     double tmp =  HIb10  + (HIb30 + 68.0*HIb20/21.0)*Sigma_out(z)/2.0;
     return tmp;
}
}



/*
double LinearPower::k_H(double z)
{
Cosmology C;
double c_light = C.H_0*3000.0/C.h;//in  km s^-1
const  double c_light_km_s =  299792.0;
double k =  (C.H_0/c_light_km_s )*C.H(z)/(C.H_0*(1.0+z));
//double k = C.H(z)/(1.0+z);
return k/C.h;
}
*/
double LinearPower::k_BAO(double z)
{
///in h/Mpc
return (1.0 +z)/105.0;
}

///////////////////////////////////////////Angular Power ///////////////////////////////////////




int k_bins4 =1000;
//std::vector <double> z_Dz_tab3(k_bins4);
 std::vector <double> Pz1_tab(k_bins4);
 std::vector <double> Pz2_tab(k_bins4);
 std::vector <double> Pz3_tab(k_bins4);
 std::vector <double> Pz4_tab(k_bins4);
 std::vector <double> Pz5_tab(k_bins4);
 std::vector <double> Pz6_tab(k_bins4);


std::vector <double> PHI1_tab(k_bins4);
 std::vector <double> PHI2_tab(k_bins4);
 std::vector <double> PHI3_tab(k_bins4);
 std::vector <double> PHI4_tab(k_bins4);

 std::vector <double> PHIG1_tab(k_bins4);
 std::vector <double> PG2_tab(k_bins4);

std::vector <double> theta_array(k_bins4);

 std::vector <double> k_array(k_bins4);

//double zini1;
////
void LinearPower::InitialzieSpeedM(double z)
{

const  double kmin2 = 0.0001;
 const  double kmax2 = 1e2;

///fnlini = fnl;
zini1 =z;
//zini2 = z2;

#pragma omp parallel private(thread_number)
{
 #pragma omp for schedule(static) nowait
for(int i = 0;i< k_bins4;i++)
{
k_array[i] = exp( ((log(kmax2)-log(kmin2))*((double) i)/((double) k_bins4-1) + log(kmin2) ));


Pz1_tab[i] = MatterPowerSpectrum(1, k_array[i],z);
Pz2_tab[i] = MatterPowerSpectrum(2, k_array[i],z);//*WindowFun(z, k_array[i]);
Pz3_tab[i] = MatterPowerSpectrum(3, k_array[i],z);//*WindowFun(z, k_array[i]);
Pz4_tab[i] = MatterPowerSpectrum(10,k_array[i],z);//*WindowFun(z, k_array[i]);

Pz5_tab[i] = P13_tt(k_array[i],z);//*WindowFun(z, k_array[i]);

/*
PHI1_tab[i] = HIPowerSpectrum(1, k_array[i],z);
PHI2_tab[i] = HIPowerSpectrum(2, k_array[i],z);
PHI3_tab[i] = HIPowerSpectrum(3, k_array[i],z);
PHI4_tab[i] = HIPowerSpectrum(10,k_array[i],z);
*/

cout<<"Matter"<<" "<<i<<"\t"<<k_array[i]<<"\t"<<Pz1_tab[i]<<"\t"<<Pz2_tab[i]<<"\t"<<Pz3_tab[i]<<"\t"<<Pz4_tab[i]<<"\t"<<Pz5_tab[i]<<endl;
}

}
}



void LinearPower::InitialzieSpeedHI(double z)
{

const  double kmin2 = 0.0001;
 const  double kmax2 = 1e4;

///fnlini = fnl;
zini1 =z;
//zini2 = z2;

#pragma omp parallel private(thread_number)
{
 #pragma omp for schedule(static) nowait
for(int i = 0;i< k_bins4;i++)
{
k_array[i] = exp( ((log(kmax2)-log(kmin2))*((double) i)/((double) k_bins4-1) + log(kmin2) ));

/*
Pz1_tab[i] = MatterPowerSpectrum(1, k_array[i],z);
Pz2_tab[i] = MatterPowerSpectrum(2, k_array[i],z);
Pz3_tab[i] = MatterPowerSpectrum(3, k_array[i],z);
Pz4_tab[i] = MatterPowerSpectrum(10,k_array[i],z);
Pz6_tab[i] = P13_tt(k_array[i],z);
*/


PHI1_tab[i] = HIPowerSpectrum(1, k_array[i],z);
PHI2_tab[i] = HIPowerSpectrum(2, k_array[i],z);//*WindowFun(z, k_array[i]);
PHI3_tab[i] = HIPowerSpectrum(3, k_array[i],z);//*WindowFun(z, k_array[i]);
PHI4_tab[i] = HIPowerSpectrum(10,k_array[i],z);//*WindowFun(z, k_array[i]);

PHIG1_tab[i] = Non_LinearPower_CrossTerms(1,k_array[i],z,z);//*WindowFun(z, k_array[i]);
PG2_tab[i] = Non_LinearPower_CrossTerms(2,k_array[i],z,z);//*WindowFun(z, k_array[i]);

cout<<"Tracer"<<" "<<i<<"\t"<<k_array[i]<<"\t"<<PHI1_tab[i]<<"\t"<<PHI2_tab[i]<<"\t"<<PHI3_tab[i]<<"\t"<<PHI4_tab[i]<<"\t"<<PHIG1_tab[i]<<"\t"<<PG2_tab[i]<<endl;
}

}
}


double MyPk(int x,double k)
{

//((age >= 18) && (age <= 35)) 
if(x == 1)
{

      Spline<double, double> CubicSpline_pk1(k_array,Pz1_tab);
      return CubicSpline_pk1.interpolate(k);    
   	
}
else if(x == 2)
{

      Spline<double, double> CubicSpline_pk2(k_array,Pz2_tab);
      return CubicSpline_pk2.interpolate(k);    
   	
}
else if(x == 3)
{

      Spline<double, double> CubicSpline_pk3(k_array,Pz3_tab);
      return CubicSpline_pk3.interpolate(k);    
   	
}
else if(x == 4)
{

      Spline<double, double> CubicSpline_pkt13(k_array,Pz5_tab);
      return CubicSpline_pkt13.interpolate(k);    
   	
}
else 
{

      Spline<double, double> CubicSpline_pkall(k_array,Pz4_tab);
      return CubicSpline_pkall.interpolate(k);    
   	
}

}

double LinearPower::MyPk_tab(int x,double k)
{
return MyPk(x, k);
}






double F_All_intrgndMatter(int x,int l,double z, double k)
{
     Cosmology C;
   double x1 =  k*chi(z);
   double x2 =  k*chi(z);
     //double power_mm = MyPk(x, k);
    double norm = l*(l+1.0)/(2.0*M_PI);
     double test2 = jl(l,x1)*jl(l,x2);
     //;//WindowFun(z, k)
if(x == 1)
{
 double tmp2 = k*k*MyPk(1, k)*test2;
return  norm*tmp2;  
}
else if(x ==2)
{
 double tmp2 = k*k*MyPk(2, k)*test2*WindowFun(z, k);
return  norm*tmp2;  
}
else if(x == 3)
{
 double tmp2 = k*k*MyPk(3, k)*test2*WindowFun(z, k);
return  norm*tmp2;  
}
else
{
    double tmp2 = MyPk(1, k) +  (MyPk(2, k)/2.0 +  MyPk(3, k))*WindowFun(z, k);
return  norm*k*k*tmp2*test2;  
}
}



double LinearPower::AngularPowerMatter(int x, int l, double z)
{
  
  const double EPSREL = 1e-6;
  const double EPSABS = 1e-6;

     
     
const  double  ln_k_min1 = 1e-4;// log((double) DLNKMIN);  
const  double  ln_k_max2 = 90;//log((double) DLNKMAX);

// const double  ln_k_min1 =  log((double) DLNKMIN);  
 //const double    ln_k_max2 = log((double) DLNKMAX);
 
 double tmp = 2.0*Integrate(bind(F_All_intrgndMatter,x,l,z,_1),ln_k_min1,ln_k_max2,EPSREL,EPSABS)/M_PI;
 return tmp;

}



////////////////////////////////////////////////////////////////////////Scale Dependence//////////////////////////////////////////


double F_All_intrgndMatterScale(int x,int l,double R,double z, double k)
{
     Cosmology C;
   double x1 =  k*chi(z);
   double x2 =  k*chi(z);
     //double power_mm = MyPk(x, k);
    double norm = l*(l+1.0)/(2.0*M_PI);
     double test2 = jl(l,x1)*jl(l,x2);
     double paramx = k*R;
     double Guassian = exp(-paramx*paramx/2.0);
if(x == 1)
{
 double tmp2 = k*k*MyPk(1, k)*test2;
return  norm*tmp2;  
}
else if(x ==2)
{
 double tmp2 = k*k*MyPk(2, k)*test2*Guassian;
return  norm*tmp2;  
}
else if(x == 3)
{
 double tmp2 = k*k*MyPk(3, k)*test2*Guassian;
return  norm*tmp2;  
}
else if(x == 4)
{
 double tmp2 = k*k*(MyPk(2, k)/2.0 + MyPk(3, k))*test2*Guassian;
return  norm*tmp2;  
}
else
{
    double tmp2 = MyPk(1, k) +  (MyPk(2, k)/2.0 +  MyPk(3, k))*Guassian;
return  norm*k*k*tmp2*test2;  
}
}



double LinearPower::AngularPowerMatterScale(int x, int l,double R, double z)
{
  
  const double EPSREL = 1e-8;
  const double EPSABS = 1e-8;
 const double EPSRELA = 1e-10;
  const double EPSABSA = 1e-10;
     
     
const  double  ln_k_min1 = 1e-4;// log((double) DLNKMIN);  
const  double  ln_k_max2 = 98;//log((double) DLNKMAX);

// const double  ln_k_min1 =  log((double) DLNKMIN);  
 //const double    ln_k_max2 = log((double) DLNKMAX);
 
if(x == 4)
{
double tmp = 2.0*Integrate(bind(F_All_intrgndMatterScale,x,l,R,z,_1),ln_k_min1,ln_k_max2,EPSRELA,EPSABSA)/M_PI;
 return tmp;

}
else
{
 double tmp = 2.0*Integrate(bind(F_All_intrgndMatterScale,x,l,R,z,_1),ln_k_min1,ln_k_max2,EPSREL,EPSABS)/M_PI;
 return tmp;
}

}




double MyHIPk(int x,double k)
{

//((age >= 18) && (age <= 35)) 
if(x == 1)
{

      Spline<double, double> CubicSpline_HIpk1(k_array,PHI1_tab);
      return CubicSpline_HIpk1.interpolate(k);    
   	
}
else if(x == 2)
{

      Spline<double, double> CubicSpline_HIpk2(k_array,PHI2_tab);
      return CubicSpline_HIpk2.interpolate(k);    
   	
}
else if(x == 3)
{

      Spline<double, double> CubicSpline_HIpk3(k_array,PHI3_tab);
      return CubicSpline_HIpk3.interpolate(k);    
   	
}


else if(x == 4)
{

      Spline<double, double> CubicSpline_Cross(k_array,PHIG1_tab);
      return CubicSpline_Cross.interpolate(k);    
   	
}
else if(x == 5)
{

      Spline<double, double> CubicSpline_G2(k_array,PG2_tab);
      return CubicSpline_G2.interpolate(k);    
   	
}


else 
{

      Spline<double, double> CubicSpline_HIpkall(k_array,PHI4_tab);
      return CubicSpline_HIpkall.interpolate(k);    
   	
}

}








double F_All_intrgndHI(int x,int l,double z, double k)
{
     Cosmology C;
   double x1 =  k*chi(z);
   double x2 =  k*chi(z);
     //double power_mm = MyHIPk(x, k);
    double norm = l*(l+1.0)/(2.0*M_PI);
     double test2 = jl(l,x1)*jl(l,x2);
    // double tmp2 =  k*k*MyHIPk(x, k)*test2;//

if(x == 1)
{
 double tmp2 = k*k*MyHIPk(1, k)*test2;
return  norm*tmp2;  
}
else if(x ==2)
{
 double tmp2 = k*k*MyHIPk(2, k)*test2;//*WindowFun(z, k);
return  norm*tmp2;  
}
else if(x == 3)
{
 double tmp2 = k*k*MyHIPk(3, k)*test2;//*WindowFun(z, k);
return  norm*tmp2;  
}
else if(x == 4)
{
 double tmp2 = k*k*(MyHIPk(2, k)/2.0 + MyHIPk(3, k))*test2;//*WindowFun(z, k);
return  norm*tmp2;  
}
else
{
    double tmp2 = MyHIPk(1, k) +  (MyHIPk(2, k)/2.0 +  MyHIPk(3, k));//*WindowFun(z, k);
return  norm*k*k*tmp2*test2;  
}
 
    
}


double LinearPower::AngularPowerHI(int x, int l, double z)
{
  
  const double EPSREL = 1e-6;
  const double EPSABS = 1e-6;

     
     
const  double  ln_k_min1 = 1e-4;// log((double) DLNKMIN);  
const  double  ln_k_max2 = 90;//log((double) DLNKMAX);

// const double  ln_k_min1 =  log((double) DLNKMIN);  
 //const double    ln_k_max2 = log((double) DLNKMAX);
 

double tmp = 2.0*Integrate(bind(F_All_intrgndHI,x,l,z,_1),ln_k_min1,ln_k_max2,EPSREL,EPSABS)/M_PI;
 return tmp;



}








double F_All_intrgndHIScale(int x,int l,double R,double z, double k)
{
     Cosmology C;
   double x1 =  k*chi(z);
   double x2 =  k*chi(z);
     //double power_mm = MyHIPk(x, k);
    double norm = l*(l+1.0)/(2.0*M_PI);
     double test2 = jl(l,x1)*jl(l,x2);
    // double tmp2 =  k*k*MyHIPk(x, k)*test2;//
double paramx = k*R;
double Guassian = exp(-paramx*paramx/2.0);

if(x == 1)
{
 double tmp2 = k*k*MyHIPk(1, k)*test2;
return  norm*tmp2;  
}
else if(x ==2)
{
 double tmp2 = k*k*MyHIPk(2, k)*test2*Guassian;
return  norm*tmp2;  
}
else if(x == 3)
{
 double tmp2 = k*k*MyHIPk(3, k)*test2*Guassian;
return  norm*tmp2;  
}
else if(x == 4)
{
 double tmp2 = k*k*(MyHIPk(2, k)/2.0 + MyHIPk(3, k))*test2*Guassian;
return  norm*tmp2;  
}
else
{
    double tmp2 = MyHIPk(1, k) +  (MyHIPk(2, k)/2.0 +  MyHIPk(3, k))*Guassian;
return  norm*k*k*tmp2*test2;  
}
 
    
}


double LinearPower::AngularPowerHIScale(int x, int l, double R, double z)
{
  
  const double EPSREL = 1e-6;
  const double EPSABS = 1e-6;

      const double EPSRELA = 1e-8;
  const double EPSABSA = 1e-8;
     
const  double  ln_k_min1 = 1e-4;// log((double) DLNKMIN);  
const  double  ln_k_max2 = 98;//log((double) DLNKMAX);

// const double  ln_k_min1 =  log((double) DLNKMIN);  
 //const double    ln_k_max2 = log((double) DLNKMAX);
 

if(x == 4)
{
 double tmp = 2.0*Integrate(bind(F_All_intrgndHIScale,x,l,R,z,_1),ln_k_min1,ln_k_max2,EPSRELA,EPSABSA)/M_PI;
 return tmp;
}
else
{
 double tmp = 2.0*Integrate(bind(F_All_intrgndHIScale,x,l,R,z,_1),ln_k_min1,ln_k_max2,EPSREL,EPSABS)/M_PI;
 return tmp;
}
}








/////////////////////////////////////////////////////////Relativistic angular power spectrum ////////////////



/*
double bevo(double z)
{
HIBias HIB;
double delta_c = 1.686;
//double tmp = (1.0+z)*CH2(z)*delta_c*(HIB.HI_b_ST_10(z)-1.0);
//return tmp;
return HIB.be_Sp(z);
}
*/

double SelectionWindow(double zstar, double z)
{

double zbin = 0.03;
double chibinstar  =  chi(zstar + zbin/2.0 ) - chi(zstar -zbin/2.0);

double norm = sqrt(2.0*M_PI*pow(chibinstar/2.0,2));

double tmp = exp(-pow((chi(z) - chi(zstar)),2)/(2.0*pow(chibinstar/2.0,2)))/norm;

return tmp;
}

double Doppler(int l,double k, double z)
{
    double tmp1 = -f2(z)*(2.0- bevo(z)/CH2(z)-Hpr(z)/CH2(z));

    double tmp =  tmp1*jlpr(l, k*chi(z))*CH2(z)/k;
  //cout<<"\t"<<"l"<<"\t"<<l<<"\t"<<"k"<<"\t"<<k<<"\t"<< "Doppler"<<"\t"<<tmp<<endl;
   return tmp;
}


double Delta_bias(double z)
{
      HIBias HIB;
      Cosmology C;
      double HIb10 = HIB.HI_b_ST_10(z);
      return  HIb10;
}

double Delta_kazier(int l, double k, double z)
{
     return -jlpr2(l, k*chi(z)) *f2(z);
}


double Deltafnl(double fnl, double k, double z)
{
        HIBias HIB;
        Cosmology C;
        double delta_c = 1.686;
         double HIb10 = HIB.HI_b_ST_10(z);
	 double HIb01 = HIB.HI_b_ST_01(fnl,z); 
         
	double g0_overg_infty = 3.0/4.0;
	double c_light =  C.H_0*3000.0/C.h;

        double numerator = 2.0*c_light*c_light*C.D_growth_Tabulated(z)*C.Tk_tab(k) * pow(CH2(z),2);
	double denominator = 3.0*C.Omega_m*pow(C.H_0,2) ;
        double alph = g0_overg_infty*numerator/denominator;
        
         double DeltaB = 2.0*fnl*delta_c*(HIb10-1.0)/alph;
        double tmp =  DeltaB;
	
	return  tmp;
}


double Delta_ISW_intrgnd(int l, double k,double z)
{
   Cosmology C;
   double Dz =   C.D_growth_Tabulated(z);
   double sph_bessel2 = jl(l,k*chi(z));
   double tmp1  =  pow(CH2(z),2)*Omz(z)*(f2(z)-1.0)*Dz *sph_bessel2;
   return tmp1;
}

double Delta_ISW(int l,double k,double z)
{ 

  Cosmology C;
  double zmin = 1e-3;
  double Dz =    C.D_growth_Tabulated(z);
  const double EPSREL = 1e-10;
   const double EPSAB = 1e-10;
  double sph_bessel2 = jl(l,k*chi(z));
  double divisor = Dz*pow(CH2(z),2);
  double Coef = -3.0*(2.0 - bevo(z)/CH2(z)+Hpr(z)/pow(CH2(z),2));

  double tmp2 = Integrate(bind(Delta_ISW_intrgnd,l,k,_1),zmin,z,EPSREL, EPSAB)/divisor;

  double tmp = Coef*tmp2;
  return tmp;
}

double Delta_SW(double z)
{
  double Coef = (3.0 - bevo(z)/CH2(z) + Hpr(z)/pow(CH2(z),2));



  double tmp = -3.0* Omz(z)*Coef/2.0;
  return tmp;
}

double Delta_TD(double z)
{
   double first = f2(z)*((-3.0 + bevo(z)/CH2(z)) - 3.0*Omz(z)/2.0) + 3.0*Omz(z)/2.0;
   return first;
}

/*
double RenormCurlyA(double fnl,double R,double k, double z)
{
//Renormbevo(double koff, double z)
 
 double second  = (3.0*Omz(z)/2.0 + 3.0*ISW(z));
  double mult_term =(2.0-Renormbevo_in(fnl,R, z)/CH2(z)-(1.0+z)*Hpr_in(z)/CH2(z));
 double tmp1 = -(first+ second*mult_term);
return tmp1;
}
*/

double Delta_together(int l, double fnl, double k, double z)
{
double sph_bessel2 = jl(l,k*chi(z));
double tmp = Delta_bias(z) +  (Deltafnl(fnl, k, z)   +  Delta_SW(z)  + Delta_TD(z))*pow(CH2(z)/k,2);// ) 
return tmp*sph_bessel2;
}

double Delta_projection(int l, double fnl, double k, double z)
{
double sph_bessel2 = jl(l,k*chi(z));
double tmp = Delta_ISW(l, k,z)*pow(CH2(z)/k,2)+  (Deltafnl(fnl, k, z)   +  Delta_SW(z)  + Delta_TD(z))*pow(CH2(z)/k,2);// ) 
return tmp*sph_bessel2;
}

double Delta_R(int x,int l, double fnl, double k, double z)
{
double sph_bessel2 = jl(l,k*chi(z));

 if(x == 1)
{
return Delta_bias(z)*sph_bessel2;
}
else if(x ==2)
{
return Delta_kazier(l, k, z);
}
else if(x == 3)
{
return Doppler(l, k, z);
}
else if(x ==4)
{
return Deltafnl(fnl, k,z)*sph_bessel2*pow(CH2(z)/k,2);
}
else if(x == 5)
{
return Delta_ISW(l,k,z)*pow(CH2(z)/k,2);
}
else if(x == 6)
{
return Delta_SW(z)*sph_bessel2*pow(CH2(z)/k,2);
}
else if(x== 7)
{
return Delta_TD(z)*sph_bessel2*pow(CH2(z)/k,2);
}
else if(x == 8)
{
return Delta_projection(l,fnl,k,z);
}
else if(x == 9)
{
return Delta_kazier(l, k, z) + Delta_bias(z)*sph_bessel2;
}
else 
{
 return Delta_kazier(l, k, z) + Doppler(l,k,z)   +  Delta_together(l,fnl,k,z);

}
}


double Transferfunintgrand(int x, int l, double fnl,double k,double zstar,double z)
{

     Cosmology C;
    
      double Dz1 = C.D_growth_Tabulated(z);
      double Tb1 = Tb_spline(z);
    
      double tmp  =  Dz1*Tb1*Delta_R(x,l,fnl, k,z);//*SelectionWindow(zstar,z)*dchi(z);
   //  cout<<"\t"<<"k"<<"\t"<<k<<"\t"<< "TF"<<"\t"<<tmp<<endl;
      return tmp;


}
double Transferfuncfirstorder(int x, int l, double fnl,double k,double z)
{
    const double EPSREL = 1e-5;
    const double EPSABS = 1e-5;
    const double zmin2 = 1e-3;
     const double zmax2 = 9.0;
     //double tmp1 = Integrate(bind(transferfunintgrand,x,l,fnl,k,z,_1),zmin2,zmax2,EPSREL,EPSABS);
      double tmp1  = Transferfunintgrand(x, l, fnl, k, z, z);
    return tmp1;

}
double F_All_intrgndTb_int(int x,int l, double fnl,double zstar1, double zstar2, double k)
{
    
     ///Cosmology C;
  
    // double power_mm = C.P_k_EH_Tabulated(k);
    // double norm = l*(l+1.0)/(2.0*M_PI);
    // double tmppower =  k*k*C.P_k_EH_Tabulated(k)*norm;//


   double tmp1 = Transferfuncfirstorder(x, l, fnl, k,zstar1);
   double tmp2 = Transferfuncfirstorder(x, l, fnl, k, zstar2);
   double tmp = tmp1*tmp2;//*tmppower;

///cout<<"\t"<<"k"<<"\t"<<k<<"\t"<< "integrand"<<"\t"<<tmp<<endl;
    return tmp;
}




const int myz_tab = 101;
std::vector<std::vector<double> > tftabA1(myz_tab,std::vector<double>(myz_tab));
 std::vector<std::vector<double> > tftabA2(myz_tab,std::vector<double>(myz_tab));
std::vector<std::vector<double> > tftabA3(myz_tab,std::vector<double>(myz_tab));
 std::vector<std::vector<double> > tftabA4(myz_tab,std::vector<double>(myz_tab));


std::vector<std::vector<double> > tftabA5(myz_tab,std::vector<double>(myz_tab));
 std::vector<std::vector<double> > tftabA6(myz_tab,std::vector<double>(myz_tab));
std::vector<std::vector<double> > tftabA7(myz_tab,std::vector<double>(myz_tab));
 std::vector<std::vector<double> > tftabA8(myz_tab,std::vector<double>(myz_tab));


 std::vector<std::vector<double> > tftabA9(myz_tab,std::vector<double>(myz_tab));
 std::vector<std::vector<double> > tftabA10(myz_tab,std::vector<double>(myz_tab));
 std::vector<std::vector<double> > tftabA11(myz_tab,std::vector<double>(myz_tab));



 std::vector<std::vector<double> > Nontftab1(myz_tab,std::vector<double>(myz_tab));
 std::vector<std::vector<double> > Nontftab2(myz_tab,std::vector<double>(myz_tab));
 std::vector<std::vector<double> > Nontftab3(myz_tab,std::vector<double>(myz_tab));

std::vector<std::vector<double> > WL_tab(myz_tab,std::vector<double>(myz_tab));

 std::vector <double> k_array2(myz_tab);
 std::vector <double> l_array(myz_tab);
 std::vector <double> l_arrayint(myz_tab);
  int lmax = 200;

void LinearPower::InitialzieTF_tabA(double fnl, double z1, double z2)
{

const  double kmin2 = 1e-4;
 const  double kmax2 = 1e2;
 //int lmax = 120;

  int  lmin = lB[1];
 //const  int lmax = lB[40];

fnlini = fnl;
zini1 = z1;
zini2 = z2;

//#pragma omp parallel private(thread_number)
{
 //#pragma omp for schedule(static) nowait
#pragma omp parallel for schedule(dynamic,1) collapse(2)
for(int i = 0;i< myz_tab;i++)
{
   for(int j = 0; j < myz_tab; j++)
    {
//l_array[i] =  lmin + i*(lmax-lmin)/myz_tab;//exp( ((log(lmax)-log(lmin))*((int) j)/((int) myz_tab-1) + log(lmin) ));//
l_arrayint[i] = lmin + i*(lmax-lmin)/myz_tab;//exp( ((log(lmax)-log(lmin))*((double) j)/((double) myz_tab-1) + log(lmin) ));//
k_array2[j] = exp( ((log(kmax2)-log(kmin2))*((double) j)/((double) myz_tab-1) + log(kmin2) ));//lB[i];
tftabA1[i][j] =  F_All_intrgndTb_int(1,l_arrayint[i], fnl, zini1,zini2, k_array2[j]);
tftabA2[i][j] =  F_All_intrgndTb_int(2,l_arrayint[i], fnl, zini1,zini2, k_array2[j]);
tftabA3[i][j] =  F_All_intrgndTb_int(3,l_arrayint[i], fnl, zini1,zini2, k_array2[j]);

tftabA4[i][j] =  F_All_intrgndTb_int(4,l_arrayint[i], fnl, zini1,zini2, k_array2[j]);
tftabA5[i][j] =  F_All_intrgndTb_int(5,l_arrayint[i], fnl, zini1,zini2, k_array2[j]);
tftabA6[i][j] =  F_All_intrgndTb_int(6,l_arrayint[i], fnl, zini1,zini2, k_array2[j]);

tftabA7[i][j] =  F_All_intrgndTb_int(7,l_arrayint[i], fnl, zini1,zini2, k_array2[j]);
tftabA8[i][j] =  F_All_intrgndTb_int(8,l_arrayint[i], fnl, zini1,zini2, k_array2[j]);
tftabA9[i][j] =  F_All_intrgndTb_int(9,l_arrayint[i], fnl, zini1,zini2, k_array2[j]);
tftabA10[i][j] =  F_All_intrgndTb_int(10,l_arrayint[i], fnl, zini1,zini2, k_array2[j]);





cout<<l_arrayint[i]<<"\t"<<k_array2[j]<<"\t"<<tftabA1[i][j]<<"\t"<<tftabA2[i][j]<<"\t"<<tftabA3[i][j]<<"\t"<<tftabA4[i][j]<<"\t"<<tftabA5[i][j]<<"\t"<<tftabA6[i][j]<<"\t"<<tftabA7[i][j]<<"\t"<<tftabA8[i][j]<<"\t"<<tftabA9[i][j]<<"\t"<<tftabA10[i][j]<<endl;
}
}


}
}





double MyFirstOrderPkintgrand(int x,double l, double k)
{

//((age >= 18) && (age <= 35)) 
if(x == 1)
{

     // double Mybicubicinterpolate(char label , std::vector <double> X,std::vector <double> Y, std::vector<std::vector<double> > F,double x, double y)
return Mybicubicinterpolate(vec[101] , l_arrayint,k_array2,tftabA1,l,k);    
   	
}
else if(x == 2)
{

     return Mybicubicinterpolate(vec[102] , l_arrayint,k_array2,tftabA2,l,k);   
   	
}
else if(x == 3)
{

     return Mybicubicinterpolate(vec[103] , l_arrayint,k_array2,tftabA3,l,k);   
   	
}
else if(x == 4)
{

     return Mybicubicinterpolate(vec[104] , l_arrayint,k_array2,tftabA4,l,k);   
   	
}
else if(x == 5)
{

     return Mybicubicinterpolate(vec[105] , l_arrayint,k_array2,tftabA5,l,k);   
   	
}
else if(x == 6)
{

     return Mybicubicinterpolate(vec[106] , l_arrayint,k_array2,tftabA6,l,k);   
   	
}
else if(x == 7)
{

     return Mybicubicinterpolate(vec[107] , l_arrayint,k_array2,tftabA7,l,k);   
   	
}
else if(x == 8)
{

     return Mybicubicinterpolate(vec[106] , l_arrayint,k_array2,tftabA8,l,k);   
   	
}
else if(x == 9)
{

     return Mybicubicinterpolate(vec[107] , l_arrayint,k_array2,tftabA9,l,k);   
   	
}
else 
{

     return Mybicubicinterpolate(vec[108] , l_arrayint,k_array2,tftabA10,l,k);   
   	
}

}




double F_All_intrgndTb(int x,int l, double fnl,double zstar1, double zstar2, double k)
{
    
     Cosmology C;
  
    // double power_mm = C.P_k_EH_Tabulated(k);
     double norm = l*(l+1.0)/(2.0*M_PI);
     double tmppower =  k*k*C.P_k_EH_Tabulated(k)*norm;//

   double tmpz1 = Transferfuncfirstorder(x, l, fnl, k,zstar1);
   double tmpz2 = Transferfuncfirstorder(x, l, fnl, k, zstar2);
   double tmp1 = tmpz1*tmpz2;//*tmppower;
   double tmp = tmp1*tmppower;

///cout<<"\t"<<"k"<<"\t"<<k<<"\t"<< "integrand"<<"\t"<<tmp<<endl;
    return tmp;
}


double F_All_intrgndTb2(int x,int l, double k)
{
    
     Cosmology C;
  
    // double power_mm = C.P_k_EH_Tabulated(k);
     double norm = l*(l+1.0)/(2.0*M_PI);
     double tmppower =  k*k*C.P_k_EH_Tabulated(k)*norm;//
   double tmp1 = MyFirstOrderPkintgrand(x,l, k);
   double tmp = tmp1*tmppower;

///cout<<"\t"<<"k"<<"\t"<<k<<"\t"<< "integrand"<<"\t"<<tmp<<endl;
    return tmp;
}

double LinearPower::AngularPowerTb(int x, int l,double fnl, double z1, double z2)
{
  
  const double EPSREL = 1e-7;
  const double EPSABS = 1e-7;

  const double EPSREL2 = 1e-10;
  const double EPSABS2 = 1e-10;

  const double EPSREL3 = 1e-20;
  const double EPSABS3 = 1e-20;

     
     
const  double  ln_k_min1 = 1e-4;// log((double) DLNKMIN);  
const  double  ln_k_max2 = 70;//log((double) DLNKMAX);

// const double  ln_k_min1 =  log((double) DLNKMIN);  
 //const double    ln_k_max2 = log((double) DLNKMAX);//F_All_intrgndTb(int x,int l,double zstar1, double zstar2, double k)
if(x == 1)
{
 
 double tmp = 2.0*Integrate(bind(F_All_intrgndTb,x,l,fnl,z1,z2,_1),ln_k_min1,ln_k_max2,EPSREL,EPSABS)/M_PI;
 return tmp;
}
if(x == 2)
{
 
 double tmp = 2.0*Integrate(bind(F_All_intrgndTb,x,l,fnl,z1,z2,_1),ln_k_min1,ln_k_max2,EPSREL,EPSABS)/M_PI;
 return tmp;
}


else if(x == 3)
{
double tmp = 2.0*Integrate(bind(F_All_intrgndTb,x,l,fnl,z1,z2,_1),ln_k_min1,ln_k_max2,EPSREL3,EPSABS3)/M_PI;
 return tmp;
}

else if(x == 4)
{
double tmp = 2.0*Integrate(bind(F_All_intrgndTb,x,l,fnl,z1,z2,_1),ln_k_min1,ln_k_max2,EPSREL3,EPSABS3)/M_PI;
 return tmp;
}

else if(x == 8)
{
double tmp = 2.0*Integrate(bind(F_All_intrgndTb,x,l,fnl,z1,z2,_1),ln_k_min1,ln_k_max2,EPSREL2,EPSABS2)/M_PI;
 return tmp;
}
else if(x == 9)
{
double tmp = 2.0*Integrate(bind(F_All_intrgndTb,x,l,fnl,z1,z2,_1),ln_k_min1,ln_k_max2,EPSREL,EPSABS)/M_PI;
 return tmp;
}

else if(x == 10)
{

double tmp = 2.0*Integrate(bind(F_All_intrgndTb,x,l,fnl,z1, z2,_1),ln_k_min1,ln_k_max2,EPSREL,EPSABS)/M_PI;
 return tmp;
}

else
{

double tmp = 2.0*Integrate(bind(F_All_intrgndTb,x,l,fnl,z1,z2,_1),ln_k_min1,ln_k_max2,EPSREL3,EPSABS3)/M_PI;
 return tmp;
}


}

//


double LinearPower::AngularPowerTb_tab(int x, int l,double fnl, double z1, double z2)
{
  
  const double EPSREL = 1e-17;
  const double EPSABS = 1e-17;

     
     
const  double  ln_k_min1 = 1e-4;// log((double) DLNKMIN);  
const  double  ln_k_max2 = 50;//log((double) DLNKMAX);

// const double  ln_k_min1 =  log((double) DLNKMIN);  MyFirstOrderPkintgrand(int x,double l, double k)
 //const double    ln_k_max2 = log((double) DLNKMAX);//F_All_intrgndTb(int x,int l,double zstar1, double zstar2, double k)
 
 double tmp = 2.0*Integrate(bind(F_All_intrgndTb2,x,l,_1),ln_k_min1,ln_k_max2,EPSREL,EPSABS)/M_PI;
 return tmp;

}



////////////////////////////////////Nonlinear part///////////////////////////////

double SecondOrderharmonicintgrand(int l, double zstar, double k,double k1,double mu,double z)
{
      Cosmology C;
    
      double Dz1 = C.D_growth_Tabulated(z);
      double Tb1 = Tb_spline(z);
      double sph_bessel2 = jl(l,k*chi(z));


      double tmp  =  pow(Dz1,2)*Tb1*(SecondOrderKernelP(5,z, k,k1,mu)*sph_bessel2 -f2(z)*GRSD_2(k, k1, mu)*jlpr2(l,k*chi(z)));//*SelectionWindow(zstar,z)*dchi(z);
   //  cout<<"\t"<<"k"<<"\t"<<k<<"\t"<< "TF"<<"\t"<<tmp<<endl;
      return tmp; 

}
double  SecondOrderharmonic( int l, double z, double k,double k1,double mu)
{
   const double EPSREL = 1e-5;
    const double EPSABS = 1e-5;
    const double zmin2 = 1e-3;
     const double zmax2 = 9.0;
     
  //double tmp1 = Integrate(bind(SecondOrderharmonicintgrand,x1,l,z,k,k1,mu,_1),zmin2,zmax2,EPSREL,EPSABS);
   double tmp1 = SecondOrderharmonicintgrand(l,z,k,k1,mu,z);
   return tmp1;
   
}

double HarmonicPowerSpectrumIntegrandHI(int l, double z1, double z2, double k,double r,double mu)
    {
     Cosmology C;
      //double r = k1/k;
     double y = sqrt(r*r-2.0*mu *r + 1.0);
    // double Dz1 = C.D_growth_Tabulated(z1);
     //double Dz2 = C.D_growth_Tabulated(z2);
     double power_multi= C.P_k_EH_Tabulated(k*r)*C.P_k_EH_Tabulated(k*y);
     double coeff =  k*r*k*r*k/(2.0*pow(M_PI,2));
     double tmp = coeff*SecondOrderharmonic(l, z1, k,k*r,mu)*SecondOrderharmonic(l,z2,k,k*r,mu)*power_multi/2.0;
    /// cout<<"\t Integrand="<<tmp1<<endl;
    return tmp;
}




 double  NonLinearPowerHIHI22(int l, double k,double z1, double z2)
{

      Cosmology C;
    //  double sph_bessel2 = jl(l,k*chi((z1+z2)/2.0));
  
      const double EPSREL = 1e-5;
     const double EPSABS = 1e-5;
    double a[2] = {QMIN, -1 };
    double b[2] = {QMAX, XMAX };
  
    // double tmp = Integrate<2>(bind(HarmonicPowerSpectrumIntegrandHI,l,z1,z2,k,_1,_2),a, b,EPSREL,EPSABS*C.P_k_EH_Tabulated(k));
   // return tmp;


  /// Cosmology C;
    
      //double Dz1 = C.D_growth_Tabulated(z1);
     //double Dz2 = C.D_growth_Tabulated(z2);
    
     double Tb1 = Tb_spline(z1);
       double Tb2 = Tb_spline(z2);
      double sph_bessel1 = jl(l,k*chi(z1));
      double sph_bessel2 = jl(l,k*chi(z2));//pow(Dz1*Dz2,2)*
    
      double tmp  =  MyHIPk(2, k)*sph_bessel1*sph_bessel2 - MyHIPk(4, k)*(sph_bessel2*jlpr2(l,k*chi(z1)) + sph_bessel1*jlpr2(l,k*chi(z2)))/2.0 + MyHIPk(5, k)*jlpr2(l,k*chi(z1))*jlpr2(l,k*chi(z2));//*SelectionWindow(zstar,z)*dchi(z);
      return Tb2*Tb1*tmp;
  
}


double  Thirdorderloopintegrand( int l, double k,double zstar, double z)
{

      Cosmology C;
      HIBias HIB;
      //double upp = -2.0/(2.0 +C.n_s);
      //double R = 5.0*pow((1.0+z),upp);///in Mpc/h
     double Dz = C.D_growth_Tabulated(z);
     double Tb1 = Tb_spline(z);
     double HIb10 = HIB.HI_b_ST_10(z); 
     double HIb20 = HIB.HI_b_ST_20(z);
     double HIb30 = HIB.HI_b_ST_30(z);
     double sph_bessel2 = jl(l,k*chi(z));
     double pkm  = pow(Dz,2)*C.P_k_EH_Tabulated(k);
//
     double p13 =  (HIb10*MyPk(3, k)  + pkm*(HIb30 + 68.0*HIb20/21.0)*Sigma_out(z))*sph_bessel2 -f2(z)*MyPk(4, k)*jlpr2(l,k*chi(z));
     double tmp =   Tb1*p13/Dz;//*SelectionWindow(zstar,z)*dchi(z);
     return tmp;
}

double  Thirdorderloop(int l, double k, double z)
{
   const double EPSREL = 1e-5;
    const double EPSABS = 1e-5;
    const double zmin2 = 1e-3;
     const double zmax2 = 9.0;
     
  //double tmp1 = Integrate(bind(Thirdorderloopintegrand,l,k,z,_1),zmin2,zmax2,EPSREL,EPSABS);
   double tmp1 = Thirdorderloopintegrand(l, k, z,z);
   return tmp1;
   
}

double HarmonicHIPowerSpectrumloop(int x,int l, double fnl, double z1, double z2, double k) 
{

     Cosmology C;
      HIBias HIB;
      
 double p22 =  NonLinearPowerHIHI22(l, k,z1,z2);//Transferfuncfirstorder(int x, int l, double fnl,double k,double z)
double p13 = (Transferfuncfirstorder(9,l,fnl,k,z1)*Thirdorderloop(l,k,z2) + Transferfuncfirstorder(9,l,fnl,k,z2)*Thirdorderloop(l,k,  z1))/2.0;
   ///double z =  (z1 +z2)/2.0;

if(x == 1)
{
return p22*WindowFun(z1, k);
}
else if(x ==2)
{
return p13*WindowFun(z1, k);
}
else
{
     double tmp  = (p22/2.0 + p13)*WindowFun(z1, k);
     return tmp;
}

}


double F_Nonlinear_intrgndTb(int x,int l, double fnl,double z1, double z2, double k)
{    
     double norm = k*k*l*(l+1.0)/(2.0*M_PI);
     double tmp =  HarmonicHIPowerSpectrumloop(x,l, fnl,z1, z2, k) ;
///cout<<"\t"<<"k"<<"\t"<<k<<"\t"<< "integrand"<<"\t"<<tmp<<endl;
    return tmp*norm;
}


 std::vector <double> k_arrayNL(myz_tab);
 //std::vector <double> l_array(myz_tab);
 std::vector <double> l_arrayNL(myz_tab);

void LinearPower::InitialzieTF_tabA2(double fnl, double z1, double z2)
{

const  double kmin2 = 1e-4;
 const  double kmax2 = 1e2;


  int  lmin = lB[1];
 //const  int lmax = lB[40];

fnlini = fnl;
zini1 = z1;
zini2 = z2;

//#pragma omp parallel private(thread_number)
{
 //#pragma omp for schedule(static) nowait
#pragma omp parallel for schedule(dynamic,1) collapse(2)
for(int i = 0;i< myz_tab;i++)
{
   for(int j = 0; j < myz_tab; j++)
    {
//l_array[i] =  lmin + i*(lmax-lmin)/myz_tab;//exp( ((log(lmax)-log(lmin))*((int) j)/((int) myz_tab-1) + log(lmin) ));//
l_arrayNL[i] = lmin + i*(lmax-lmin)/myz_tab;//exp( ((log(lmax)-log(lmin))*((double) j)/((double) myz_tab-1) + log(lmin) ));//
k_arrayNL[j] = exp( ((log(kmax2)-log(kmin2))*((double) j)/((double) myz_tab-1) + log(kmin2) ));//lB[i];

Nontftab1[i][j] =   F_Nonlinear_intrgndTb(1,l_arrayNL[i],fnl,zini1,zini2, k_arrayNL[j]);
Nontftab2[i][j] =   F_Nonlinear_intrgndTb(2,l_arrayNL[i],fnl,zini1,zini2, k_arrayNL[j]);
Nontftab3[i][j] =   F_Nonlinear_intrgndTb(3,l_arrayNL[i],fnl,zini1,zini2, k_arrayNL[j]);
//double F_Nonlinear_intrgndTb(int x,int l, double fnl,double zstar1, double zstar2, double k)

cout<<l_arrayNL[i]<<"\t"<<k_arrayNL[j]<<"\t"<<"\t"<<Nontftab1[i][j]<<"\t"<<Nontftab2[i][j]<<"\t"<<Nontftab3[i][j]<<endl;
}
}


}
}


double MyNonlinearOrderPkintgrand(int x,double l, double k)
{
if(x == 1)
{
return Mybicubicinterpolate("A9", l_arrayint,k_array2,Nontftab1,l,k);
}
else if(x == 2)
{
return Mybicubicinterpolate("A10", l_arrayint,k_array2,Nontftab2,l,k);
}
else
{
return Mybicubicinterpolate("A11", l_arrayint,k_array2,Nontftab3,l,k);
}
}

double F_Nonlinear_intrgndTb_in(int x,int l, double k)
{    
     ///double norm = k*k*l*(l+1.0)/(2.0*M_PI);
     double tmp =  MyNonlinearOrderPkintgrand(x,l, k) ;
///cout<<"\t"<<"k"<<"\t"<<k<<"\t"<< "integrand"<<"\t"<<tmp<<endl;
    return tmp;///*norm;
}
double LinearPower::AngularNonlinearPowerTb_tab(int x, int l,double fnl, double z1, double z2)
{
  
  const double EPSREL = 1e-12;
  const double EPSABS = 1e-12;

     
     
const  double  ln_k_min1 = 1e-4;// log((double) DLNKMIN);  
const  double  ln_k_max2 = 50;//log((double) DLNKMAX);

// const double  ln_k_min1 =  log((double) DLNKMIN);  MyFirstOrderPkintgrand(int x,double l, double k)
 //const double    ln_k_max2 = log((double) DLNKMAX);//F_All_intrgndTb(int x,int l,double zstar1, double zstar2, double k)
 
 double tmp = 2.0*Integrate(bind(F_Nonlinear_intrgndTb_in,x,l,_1),ln_k_min1,ln_k_max2,EPSREL,EPSABS)/M_PI;
 return tmp;

}


double LinearPower::AngularNonlinearPowerTb(int x, int l,double fnl, double z1, double z2)
{
  
  const double EPSREL = 1e-8;
  const double EPSABS = 1e-8;

     
     
const  double  ln_k_min1 = 1e-4;// log((double) DLNKMIN);  
const  double  ln_k_max2 = 50;//log((double) DLNKMAX);

// const double  ln_k_min1 =  log((double) DLNKMIN);  MyFirstOrderPkintgrand(int x,double l, double k)
 //const double    ln_k_max2 = log((double) DLNKMAX);//F_All_intrgndTb(int x,int l,double zstar1, double zstar2, double k)
 //F_Nonlinear_intrgndTb(int x,int l, double fnl,double z1, double z2, double k)
 double tmp = 2.0*Integrate(bind(F_Nonlinear_intrgndTb,x,l,fnl,z1,z2,_1),ln_k_min1,ln_k_max2,EPSREL,EPSABS)/M_PI;
 return tmp;

}




