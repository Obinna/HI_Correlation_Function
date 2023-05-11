/*
 *  RealSpacePower.h
 *
 *
 *  Created by Obinna Umeh on 11/10/2020.
 *  Copyright 2020 University of Portsmouth. All rights reserved.
 *
 */

#include "RealSpacePower.h"

#include "MyHIBias.h"

#include <boost/math/special_functions/legendre.hpp>


#include <boost/bind.hpp>
//#include "MyBessel.cpp"

using boost::bind;
using boost::cref;


const double QMIN = 1e-4;
const double QMAX = 1e4;
const double XMAX = 0.999999;

const  double zmin = 1e-6;


double fnlini = 0.0;
double zini1 = 0.0;
double zini2 = 0.0;

RealSpacePower::RealSpacePower()
{
    //SetPowerSpectrum();
}

RealSpacePower::~RealSpacePower()
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
    double x =k*R;///(1.0+z);
	
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




double f22(double z)
{
	Cosmology C;
	double tmpz = C.Omega_m*pow((1.0+z),3)/(pow(C.E(z),2));
	return pow(tmpz,4.0/7.0);
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




double Sntgrd(double z, double k)
{
    Cosmology C;
    HIBias HIB;
    //double Tb = T_b(z);
    double HIb20 = HIB.HI_b_ST_20(z);
    double Dz = C.D_growth_Tabulated(z);
   //WindowFun(double theta,double z, double k)
    double coef = pow( HIb20,2)*pow(Dz,4)*pow(k*C.P_k_EH_Tabulated(k),2)/(4.0*pow(M_PI,2));
    return coef;

}
double SntgrdR(double z,double R, double k)
{
    Cosmology C;
    HIBias HIB;
    //double Tb = T_b(z);
    double HIb20 = HIB.HI_b_ST_20(z);
    double Dz = C.D_growth_Tabulated(z);
   //WindowFun(double theta,double z, double k)
    double coef = pow(WindowFun2(z,R, k)* HIb20,2)*pow(Dz,4)*pow(k*C.P_k_EH_Tabulated(k),2)/(4.0*pow(M_PI,2));
    return coef;

}

double RealSpacePower::NonlinearShotnoise(double z)
{
      const double EPSREL = 1e-10;
      const double EPSAB = 1e-10;
       double QMIN = 1e-4;
       double QMAX = 10e4;
      double tmp2 = Integrate(bind(Sntgrd,z,_1),QMIN,QMAX,EPSREL,EPSAB);
       //cout<<"\t i"<<"\t"<<"z"<<"\t"<<z <<"\t"<<"ShotNoise_tab"<<"\t"<<"\t"<<tmp2 <<endl;
       return tmp2;
}

double RealSpacePower::NonlinearShotnoiseR(double z, double R)
{
      const double EPSREL = 1e-20;
      const double EPSAB = 1e-20;
       double QMIN = 1e-4;
       double QMAX = 10e4;
      double tmp2 = Integrate(bind(SntgrdR,z,R,_1),QMIN,QMAX,EPSREL,EPSAB);
       //cout<<"\t i"<<"\t"<<"z"<<"\t"<<z <<"\t"<<"ShotNoise_tab"<<"\t"<<"\t"<<tmp2 <<endl;
       return tmp2;
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

/*
double RealSpacePower::HproverH(double z)
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
*/


int k_bins2 =50;
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



void RealSpacePower::InitialzieSpeed4(void)
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
            //Sigderi_tab1[i]= Mysigmaderivative(exp(z_Dz_tab2[i])); //NonlinearShotnoise(double z)

			//Sn_tab[i] = NonlinearShotnoise(exp(z_Dz_tab2[i]));
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
/*
double Sigmav_out(double z)
{
	//Spline<double, double> CubicSpline_sigmav(z_Dz_tab2,Sigv_tab1);
	return 0.0;//CubicSpline_sigmav.interpolate(log(z));
   	
}

double Sigmaderi_out(double z)
{
	Spline<double, double> CubicSpline_sigma(z_Dz_tab2,Sigderi_tab1);
	return CubicSpline_sigma.interpolate(log(z));
   	
}*/

/*
double ShotNoisePK(double z)
{
	Spline<double, double> CubicSpline_Sn(z_Dz_tab2,Sn_tab);

	return CubicSpline_Sn.interpolate(log(z));
   	
}
*/
/*
double RealSpacePower::ShotNoiseout(double z)
{
 return ShotNoisePK(z);
}

double RealSpacePower::sigmasqout(double z)
{
 return NonlinearShotnoise(z);//Sigma_out(z);
}

*/






double RealSpacePower::k_H(double z)
{
	Cosmology C;
	double c_light = C.H_0*3000.0;//in  km s^-1
	const  double c_light_km_s =  299792.0;
	double k =  (1.0/c_light_km_s )*C.H(z)/(1.0+z);
	//double k = C.H(z)/(1.0+z);
	return k;
}








/////////////////////////////////////////Matter power psectrum///////////////////////////////////////////////



double RealSpacePower::MatterPowerSpectrum(int x, double k,double z)
{
   if(x == 1)
{
    return LinearPower_mm(k,z);
}
    else if(x == 2)
{
   return PowerSpectrum222(k,z);
}
else if(x == 3)
{
     return P13_ddint(k,z);
}
else
{
    double tmp  = LinearPower_mm(k,z) + (PowerSpectrum222(k,z) + P13_ddint(k,z));
return tmp;
}

}


double RealSpacePower::MatterPowerSpectrumR(int x, double k,double z, double R)
{
   if(x == 1)
{
    return LinearPower_mm(k,z)*pow(WindowFun2(z,R, k),2);
}
    else if(x == 2)
{
   return PowerSpectrum222(k,z)*pow(WindowFun2(z, R, k),2);
}
else if(x == 3)
{
     return P13_ddint(k,z)*pow(WindowFun2(z, R, k),2);
}
else
{
    double tmp  = LinearPower_mm(k,z) + (PowerSpectrum222(k,z) + P13_ddint(k,z));
return tmp*pow(WindowFun2(z,R, k),2);
}

}



double LinearPower_mm( double k, double z1)
{
    Cosmology C;
    
     double tmp  = pow(C.D_growth_Tabulated(z1),2)*C.P_k_EH_Tabulated(k);
    return  tmp;
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
double PowerSpectrum222(double k, double z)
 {
      Cosmology C;
  const double EPSREL = 1e-5;
     const double EPSABS = 1e-5;
     double Dz = C.D_growth_Tabulated(z);
    double a[2] = {QMIN, -1 };
    double b[2] = {QMAX, XMAX };
    return pow4(Dz)*k*k/(98*4*M_PI*M_PI)*Integrate<2>(bind(f22_dd, k, _1, _2), a, b, EPSREL, EPSABS*C.P_k_EH_Tabulated(k));
}



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













//////HI power spectrum///////////////////////////////////////////////////////////////




double RealSpacePower::HIPowerSpectrum(int x, double k,double z)
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

     double p13 =  HIb10*HIb10*P13_ddint(k,z)  + HIb10*pkm*(HIb30 + 68.0*HIb20/21.0)*Sigma_out(z);


     if(x == 1)
       {
         return  pow(HIb10*Dz,2)*C.P_k_EH_Tabulated(k);
       }
else if(x ==2)
{
///return pow(T_b(z),2)*(Non_RealSpacePower_HIHI(10,10, k,z)-Non_RealSpacePower_HIHI(2,2, 10e-3, z));//-ShotNoise(z);//*WindowFun(z, k);
return Non_RealSpacePower_HIHI(k,z);//-ShotNoise(z);//*WindowFun(z, k);
}
else if(x ==3)
{
return p13;//*WindowFun(z, k);
}
else
{
double tmp  = pow(HIb10*Dz,2)*C.P_k_EH_Tabulated(k) + (Non_RealSpacePower_HIHI(k,z)/2.0 + p13);
    
return tmp;
}

}




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



double PowerSpectrumIntegrandHI( double z, double k,double k1,double mu)
    {
     Cosmology C;
      double r = k1/k;
     double y = sqrt(r*r-2.0*mu *r + 1.0);
     double Dz = C.D_growth_Tabulated(z);
     double power_multi= pow(Dz,4)*C.P_k_EH_Tabulated(k*r)*C.P_k_EH_Tabulated(k*y);
     double coeff =  k1*k1/(4.0*pow(M_PI,2));
     double tmp = coeff*SecondOrderKernelP(10,z, k,k*r,mu)*SecondOrderKernelP(10,z, k,k*r,mu)*power_multi;
    /// cout<<"\t Integrand="<<tmp1<<endl;
    return tmp;
}




 double  Non_RealSpacePower_HIHI( double k,double z1)
{
    Cosmology C;
    const double EPSREL = 1e-5;
    const double EPSABS = 1e-5;
    double a[2] = {QMIN, -1 };
    double b[2] = {QMAX, XMAX };
    double tmp = Integrate<2>(bind(PowerSpectrumIntegrandHI,z1,k,_1,_2),a, b,EPSREL,EPSABS);
    return tmp;
  
}




/////////////////////////////////////////////Windowned HI power spectrum ///////////////////////////////////////////
double RealSpacePower::HIPowerSpectrumR(int x, double k,double z,double R)
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

     double p13 =  HIb10*HIb10*pow(WindowFun2(z,R, k),2)*P13_ddint(k,z)
    + WindowFun2(z,R, k)*HIb10*pkm*(WindowFun2(z,R, k)*HIb30 + 68.0*HIb20/21.0)*MysigmaR(z,R);
//pow(WindowFun2(z,R, k),2)

if(x == 1)
{
return pow(HIb10*Dz*WindowFun2(z,R, k),2)*C.P_k_EH_Tabulated(k);
}
else if(x ==2)
{
return Non_RealSpacePower_HIHIWindow(k,z,R);//*WindowFun(z, k);
}
else if(x ==3)
{
return p13;//*WindowFun(z, k);
}
else
{
double tmp  = pow(HIb10*Dz*WindowFun2(z,R, k),2)*C.P_k_EH_Tabulated(k) + (Non_RealSpacePower_HIHIWindow(k,z,R)/2.0 + p13);
return tmp;
}

}



double SecondOrderKernelPWindow(int x,double R, double z, double k,double k1,double mu)
{     HIBias HIB;
     double HIb10 = HIB.HI_b_ST_10(z);
     double HIb20 = HIB.HI_b_ST_20(z);
     double r = k1/k;
     double y = sqrt(r*r-2.0*mu *r + 1.0);//pow(WindowFun2(z,R, k),2)
if(x == 1)
{
return HIb10*F2P(k,k*r,mu)*WindowFun2(z,R, k);//F2P(k,k*r,mu);F2P(double k, double k1, double mu1)
}
else if(x ==2)
{
return HIb20*WindowFun2(z,R, k*r)*WindowFun2(z,R, k*y);
}
else
{
double tmp  = HIb10*WindowFun2(z,R, k)*F2P(k,k*r,mu) + HIb20*WindowFun2(z,R, k*r)*WindowFun2(z,R, k*y);
return tmp;
}
}

double PowerSpectrumIntegrandHIWindow(double R, double z, double k,double k1,double mu)
    {
     Cosmology C;
      double r = k1/k;
     double y = sqrt(r*r-2.0*mu *r + 1.0);
     double Dz = C.D_growth_Tabulated(z);
     double power_multi= pow(Dz,4)*C.P_k_EH_Tabulated(k*r)*C.P_k_EH_Tabulated(k*y);
     double coeff =  k1*k1/(4.0*pow(M_PI,2));
     double tmp = coeff*SecondOrderKernelPWindow(10,R,z, k,k*r,mu)*SecondOrderKernelPWindow(10,R,z, k,k*r,mu)*power_multi;
    /// cout<<"\t Integrand="<<tmp1<<endl;
    return tmp;
}




 double  Non_RealSpacePower_HIHIWindow(double k,double z1, double R)
{
    Cosmology C;
    const double EPSREL = 1e-5;
    const double EPSABS = 1e-5;
    double a[2] = {QMIN, -1 };
    double b[2] = {QMAX, XMAX };
    double tmp = Integrate<2>(bind(PowerSpectrumIntegrandHIWindow,R,z1,k,_1,_2),a, b,EPSREL,EPSABS);
    return tmp;
  
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



double  PowerSpectrum22(double k,double z)
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









double PowerSpectrum13(double k,double z)
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



/*
double RealSpacePower::DeltaPowerSpectrum(int x, double k,double z)
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


double RealSpacePower::beffrealspace(double fnl,double k,double z)
{
 double tmp = sqrt(HIPowerSpectrum_NG(10,fnl,  k, z) /MatterPowerSpectrum(10, k,z) );
return tmp;
}

*/



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

