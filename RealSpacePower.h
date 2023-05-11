/*
 *  RealSpacePower.h
 *  
 *
 *  Created by Obinna Umeh on 11/10/2020.
 *  Copyright 2020 University of Portsmouth. All rights reserved.
 *
 */

//#ifndef  _RealSpacePower
//#define  _RealSpacePower



class RealSpacePower
	{
    public:
		
		
		
        double HIPowerSpectrum(int x, double k,double z);
        double HIPowerSpectrumR(int x, double k,double z, double R);
        
        double MatterPowerSpectrum(int x, double k,double z) ;
        double MatterPowerSpectrumR(int x, double k,double z, double R) ;
		
	

       // void InitialzieSpeedM(double z);
       // void InitialzieSpeedHI(double z);
       // double sigmasqout(double z);
        double NonlinearShotnoise(double z);
        double NonlinearShotnoiseR(double z, double R);

       // double ShotNoiseout(double z);

        double k_H(double z);
        double k_BAO(double z);
        double DeltaPowerSpectrum(int x, double k,double z);
        
        double beffrealspace(double fnl,double k,double z);
  
		
		void InitialzieSpeed4(void);
        RealSpacePower();// { cout << "\n\t Creating Linear power nspectrum";}
		// void SetPowerSpectrum(void);
		~RealSpacePower();// { cout << "\n\t Destroying  Linear power nspectrum";}
		
		
		
		
	};


double F2P(double k, double k1, double mu1);
double LinearPower_mm( double k, double z1);
//static double f22_dd(double k, double q, double x);
double PowerSpectrum222(double k, double z);
double P13_ddint(double k,double z) ;

double SecondOrderKernelP(int x,double z, double k,double k1,double mu);
double PowerSpectrumIntegrandHI(double z, double k,double k1,double mu);
double  Non_RealSpacePower_HIHI(double k,double z1);


double SecondOrderKernelPWindow(int x,double R, double z, double k,double k1,double mu);
double PowerSpectrumIntegrandHIWindow(int x1,double R, double z, double k,double k1,double mu);
double  Non_RealSpacePower_HIHIWindow( double k,double z1, double R);

double IntegrandP22(double z, double k, double r, double mu);
double IntegrandP13(double z, double k, double r);

double  NonLinearPowerHIHI22( int l, double k,double z1, double z2);
double PowerSpectrumIntegrandHI(int x1,int x2, double z, double k,double k1,double mu);
double  Non_RealSpacePower_HIHI(int x1,int x2, double k,double z1);


double  PowerSpectrum13(double k,double z);
double  PowerSpectrum22(double k,double z);

static inline double pow2(double x) { return x*x; }
static inline double pow3(double x) { return x*x*x; }
static inline double pow4(double x) { return pow2(pow2(x)); }
static inline double pow5(double x) { return x*pow4(x); }
static inline double pow6(double x) { return pow3(pow2(x)); }
static inline double pow7(double x) { return x*pow6(x); }
static inline double pow8(double x) { return pow2(pow4(x)); }

double F2Simp(double k, double k1, double mu1);
double G2Simp(double k, double k1, double mu1);
double FF2(double k1, double k2, double mu12);
double GG2(double k1, double k2, double mu12);
















	
