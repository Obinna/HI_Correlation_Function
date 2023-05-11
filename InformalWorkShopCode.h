/*
 *  InformalWorkShopCode.h
 *  
 *
 *  Created by obinna umeh on 28/02/2016.
 *  Copyright 2016 University of the Western Cape. All rights reserved.
 *
 */

//#ifndef  _InformalWorkShopCode
//#define  _InformalWorkShopCode



class LinearPower
	{
    public:
		
		
		double LinearPower_HIHI(double fnl, double k, double z1,double mu);
		
		double CompareKaizer(double fnl, double k, double z1,double mu);
		
		double LinearPower_mm(double k, double z);
               double  PowerSpectrum13(double k,double z);
               double  PowerSpectrum22(double k,double z);
              double PowerSpectrum222(double k,double z);
		 double P13_dd( double k,double z) const;
		double DM(double z);
		double u_perp(double z, double k, double mu);
		double u_para(double z, double k, double mu);
		double nuofz(double z);
		
		double Power_Noise(int choose_exp,double S_area,double tobs,double N_dish, double D_dish,double Dnu, double nu,double w, double u);
		double PowerSpectrum_tot(int choose_exp,double fnl,double S_area,double tobs,double N_dish,double D_dish,double Dnu, double nu,double u_para, double u_perp,double mu);
		double Cov_PowerSpectrum(int choose_exp,double fnl,double S_area,double tobs,double N_dish,double D_dish, double nu,double u_para1, double u_para2,double u_perp1,double u_perp2, double mu);
		
		double LinearPower_HIHIRSD(double fnl, double k, double z1,double mu);
		double  LinearHIHI_Average( int l, double fnl, double k, double z);
		double LinearHIHI_AverageRSD( int l, double fnl, double k, double z);
		double DeltaHIHI_Average(int x1, int x2, int l, double fnl, double k, double z);
		///double k_H(double z);
		
		
		double LinearPower_HIHI_intgrand(int l, double fnl, double k, double z,double mu);
		double LinearPower_HIHI_intgrand2(int l, double fnl, double k, double z,double mu);
		
		 double HIPowerSpectrum(int x, double k,double z);
		 double HIPowerSpectrumR(int x, double k,double z,double R);
                double HIPowerSpectrum_NG(int x, double fnl, double k,double z) ;
                  double Non_LinearPower_HIHI(int x1,int x2, double k,double z1);
                 double DeltaPowerSpectrum(int x, double k,double z);
             double Non_LinearPower_HIHItidal(int x1,int x2, double k,double z1);

                 double MatterPowerSpectrum(int x, double k,double z) ;
double RenormNon_LinearPower_HIHIRSD22out(int x1,int x2,double fnl, double k,double  z,double mu);
double P22_ttout(double k, double z) ;

 double Adhocbias(double fnl, double k, double z1,double mu);
double fg_out(double z);

               void InitialzieSpeedM(double z);
             void InitialzieSpeedHI(double z);
           double sigmasqout(double z);

 double ShotNoiseout(double z);

 double k_H(double z);
     double k_BAO(double z);

    double beffredshift(double fnl, double k,double z);
   double beffrealspace(double fnl,double k,double z);
double ScoccimarroPk(double  k, double z, double mu);
double ScoccimarroPkAverage(double  k, double z);
double AllPk(int x, double  k, double z);
double ExportCurlyIterm1(double k,double z,double mu);
double ExportCurlyIterm2(double k,double z,double mu);
double beffredshift2(double fnl, double k,double z,double mu);
double beffrealspace2(int x, double z);
    double TPowerSpectrum(int x,double fnl, double k,double z,double mu);
     double TPowerSpectrum_Average(int x, double fnl, double k, double z);




//////////////////////////////Matter power spectrum/////////////////////////////////////////


double LinearPower_HIHIRSDmatter(double fnl, double k, double z1,double mu);
double LinearPower_HIHI_intgrandmatter(int l, double fnl, double k, double z,double mu);
double LinearHIHI_AverageRSDmatter(int l, double fnl, double k, double z);
double beffHImatter(int x, double fnl, double k,double z);
double TPowerSpectrummatter(int x,double fnl, double k,double z,double mu) ;
double TPowerSpectrum_Averagematter(int x, double fnl, double k, double z);

  







////////////////////////////////////////Angular power spectrum/////////////////////////////////

               double AngularPowerMatter(int x, int l, double z);
               double AngularPowerHI(int x, int l, double z);
              double AngularPowerTb(int x, int l,double fnl, double z1, double z2);
              double AngularPowerTb_tab(int x, int l,double fnl, double z1, double z2);
             void InitialzieTF_tabA(double fnl, double z1, double z2);
             void InitialzieTF_tabA2(double fnl, double z1, double z2);
               double MyPk_tab(int x,double k);
                double AngularNonlinearPowerTb_tab(int x, int l,double fnl, double z1, double z2);
		double AngularNonlinearPowerTb(int x, int l,double fnl, double z1, double z2);
		double HproverH(double z);
		double AngularPowerHIScale(int x, int l, double R, double z);
              double AngularPowerMatterScale(int x, int l,double R, double z);
		
		void InitialzieSpeed4(void);
		LinearPower();// { cout << "\n\t Creating Linear power nspectrum";}
		// void SetPowerSpectrum(void);
		~LinearPower();// { cout << "\n\t Destroying  Linear power nspectrum";}
		
		
		
		
	};





double IntegrandP22(double z, double k, double r, double mu);
double IntegrandP13(double z, double k, double r);
double F_All_intrgndTb(int x,int l,double fnl,double zstar1, double zstar2, double k);
double HarmonicHIPowerSpectrumloop(int x,int l, double fnl, double z1, double z2, double k) ;
double  Thirdorderloop(int l, double k, double z);
 double  NonLinearPowerHIHI22( int l, double k,double z1, double z2);
double F_Nonlinear_intrgndTb(int x,int l, double fnl,double zstar1, double zstar2, double k);
double F_Nonlinear_intrgndTb2(int x,int l, double fnl,double z1, double z2, double k);
double P13_ddint(double k,double z) ;
double F_All_intrgndTb_int(int x,int l, double fnl,double zstar1, double zstar2, double k);

double GRSD_2( double k, double q, double muk);

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
















	
