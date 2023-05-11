/*! \file cosmology.c
 *  \brief Definitions for the Cosmology class that contains
 *         useful cosmological calculations.
 */
#include<iostream>
#include<stdio.h>
#include<math.h>
//#include<cmath>
//#include <vector>
#include"MyCosmology.h"
#include"Copter_Jeremy/rk_int.h"
#include"Copter_Jeremy/routines.h"

#include <boost/bind.hpp>
#include <boost/math/special_functions/bessel.hpp>
using boost::cref;
using boost::bind;
#include "Copter_Jeremy/Quadrature.h"
#include "Copter_Jeremy/Quadrature.cpp"



#include "Copter_Jeremy/MySpline.cpp"


using namespace std;
#define DLNKMIN 1.0e-6 //real space top hat min k table value
#define DLNKMAX 1.0e5 //real space top hat max k table value


int thread_number =8;
 int     n_ln_z =150;
 int    n_W_x =150;
 int    n_ln_k =150;

	    //tabulated log redshift
std::vector <double> z_Dz_tab(n_ln_k);
std::vector <double> Dz_tab(n_ln_k);  //tabulated growth function
std::vector <double> x_W_x(n_ln_k);
std::vector <double> W2_x(n_ln_k); //tabulated window function
std::vector <double> k_EH(n_ln_k);
std::vector <double> P_k_EH(n_ln_k);
std::vector <double> R_EH(n_ln_k);
std::vector <double> M_EH(n_ln_k);
std::vector <double> Delta2_k_EH(n_ln_k);
std::vector <double> sigma_EH_tab(n_ln_k);
std::vector <double> Tk_EH_tab(n_ln_k);

std::vector <double> PowM_EH(n_ln_k);



//int    n_W_x;
//int    n_ln_k;
/*
vector <double> x_W_x(n_ln_k); //tabulated window function
vector <double> W2_x(n_ln_k);
vector <double> k_EH(n_ln_k);
vector <double> M_EH(n_ln_k);
vector <double> R_EH(n_ln_k);
vector <double> P_k_EH(n_ln_k);
vector <double> Delta2_k_EH(n_ln_k);
vector <double> sigma_EH_tab(n_ln_k);
*/



/*! \fn Cosmology(void)
 *  \brief Constructor for the Cosmology class. */
Cosmology::Cosmology(void)
{
    InitializeCosmology();
}

void Cosmology::SetUpTools(void)
{
 InitializeSigma_EH();
    InitializeDz();
  TabulateSigma_EH();
}

/*! \fn InitializeCosmology(void)
 *  \brief Initializer for the Cosmology class. */
void Cosmology::InitializeCosmology(void)
{
    // Hubble constant at z=0 in h km/s/Mpc
    H_0 = 100.0;

    // Inverse of Hubble constant in h^-1 Gyr
    H_inv = 1 / (H_0 * 1.0e3 / mpc * year_in_sec * 1.0e9);

    //critical density in h^2 Msun / Mpc^3
    //Dodelson eqn 1.3
    rho_crit = 3.0*H_0*H_0/(8*pi*G_cosmo*1.0e-3); 

    //critical density in h^2 g / cm^3
    rho_c    = rho_crit * msun_cgs / pow(mpc_cgs,3);


    //number of samples in tabulated power spectra
   // n_ln_k   = 1000;

    //largest scale for tabulated power spectrum 
    ln_k_min = log((double) DLNKMIN); 

    //smallest scale for tabulated power spectrum 
    ln_k_max = log((double) DLNKMAX);

    //set sigma tabulation flags and normalization
    //initialize_sigma_flag      = 1;
   // initialize_sigma_EH_flag   = 1;
    //tabulated_sigma_EH_flag    = 1;
    sigma_EH_renormalization   = 1.0;

    //set the empty lookup table pointers to null
    //SetSigmaMemoryToNull();


    //Pick a Cosmology
    //SetWMAP9Cosmology();
    //SetPrelimCosmology();
    //SetZentner2007Cosmology();
    //SetMillenniumCosmology();
    SetPLANCKCosmology();

    //Initialize sigma and DZ
  //InitializeSigma_EH();
   // InitializeDz();
   

}

/*! \fn Cosmology(void)
 *  \brief Destructor for the Cosmology class. */
Cosmology::~Cosmology(void)
{

    //Nothing to be done!

 /*   //if the growth factor lookup table has
    //been initialize, free its memory
    if(initialize_Dz_flag==1)
    {
        free(Dz_tab);
	free(z_Dz_tab);
	gsl_spline_free(spline_D_z);
	gsl_interp_accel_free(acc_D_z);
    }

    //if the EH98 sigma lookup table has been initialized,
    //free its memory
    if(initialize_sigma_EH_flag==1)
    {
	gsl_spline_free(spline_sigma_EH);
	gsl_interp_accel_free(acc_sigma_EH);
	gsl_spline_free(spline_sigma_M_EH);
	gsl_interp_accel_free(acc_sigma_M_EH);
    }

*/
    //if the sigma lookup table has been initialized,
    //free its memory
/*
    if(initialize_sigma_flag==1)
    {
	FreeSigmaMemory();
    }*/

}


///Cosmology member functions below here

/*! \fn double ABTonFlux(double AB)
*  \brief AB Magnitude to Flux, where F_nu is in ergs/s/cm^2/Hz*/
double Cosmology::ABToFlux(double AB)
{
	//to F_nu
	return pow(10.0, -0.4*(AB+48.6));
}

/*! \fn double ABTonJy(double AB)
*  \brief AB Magnitude to Nanojansky, where nJy is in nJy */
double Cosmology::ABTonJy(double AB)
{
	//to nJy
	return pow(10.0, -0.4*(AB-31.4));
}


/*! \fn double ComovingDistance(double z_min, double z_max)
 *  \brief Comoving distance in h^-1 Mpc */
double Cosmology::ComovingDistance(double z_min, double z_max)
{
	//comoving distance
	//Dodelson equation 2.89
	//in h^-1 Mpc
	//Checked against Cosmology Calculator

	double fp[4];
	double dz_int;
	double eps = 1.0e-6;

	if(z_min==z_max)
		return 0;

	dz_int = 0.01*(z_max-z_min);

	fp[0] = Omega_m;
	fp[1] = Omega_l;
	fp[2] = Omega_k;
	fp[3] = Omega_r;

	//in h^-1 Mpc
	return ((c/1.0e3)/H_0)*integrate(ComovingDistanceIntegrand,fp,4,z_min,z_max,dz_int,eps);
}

/*! \fn double ComovingVolume(double z_min, double z_max)
 *  \brief Comoving volume in h^-3 Mpc^3 */
double Cosmology::ComovingVolume(double z_min, double z_max)
{
        //Checked against Cosmology Calculator
    double r31 = pow(ComovingDistance(0,z_min),3.);
    double r32 = pow(ComovingDistance(0,z_max),3.);
    return 4.0*M_PI/3.0 * (r32 - r31);
}


/*! \fn double D_growth(double z)
 *  \brief Linear growth factor as a function of redshift*/
double Cosmology::D_growth(double z)
{
	//linear growth function

	return D_growth_integrated(z)/D_growth_0;
}

/*! \fn double D_growth(double z)
 *  \brief Unnormalized linear growth factor as a function of redshift*/
double Cosmology::D_growth_integrated(double z)
{
	//Growth factor

	double fp[4];
	double dz_int;
	double eps = 1.0e-9;
	double z_max = log(1.0e20);

	if(z>=exp(z_max))
		return 0;

	if(z<1.0e-9)
		z = 1.0e-9;

	dz_int = 0.001*(z_max-log(z));

	fp[0] = Omega_m;
	fp[1] = Omega_l;
	fp[2] = Omega_k;
	fp[3] = Omega_r;

	return h*H(z)*pow(H_inv/h,3)*integrate(D_growth_integrand,fp,2,log(z),z_max,dz_int,eps);
}

/*! \fn double DistanceModulus(double z)
 *  \brief Distance modulus to redshift z in delta mag for a flat F_nu source*/
double Cosmology::DistanceModulus(double z)
{
    double dl = 1.0e6*LuminosityDistance(z)/h;    //luminosity distance in pc
    double K_correction = -2.5*log10(1+z);        //K correction, see Hogg 1999, eq 27

    return 5.0*(log10(dl)-1.0) + K_correction;
}

/*! \fn double E(double z)
 *  \brief Dimensionless Hubble parameter at redshift z */
double Cosmology::E(double z)
{
	//returns unit-free
	//Hubble parameter 
	//at redshift z 

	return sqrt( Omega_m*pow(1.+z,3) + Omega_k*pow(1+z,2) + Omega_l + Omega_r*pow(1+z,4));
}

/*! \fn double FluxToAB(double F_nu)
 *  \brief Flux to AB Magnitude, where F_nu is in erg/s/cm^2/Hz */
double Cosmology::FluxToAB(double F_nu)
{
	//returns AB mags
	return -2.5*log10(F_nu) - 48.60;
}


/*! \fn double H(double z)
 *  \brief Hubble parameter at redshift z */
double Cosmology::H(double z)
{
	//returns 
	//Hubble parameter 
	//at redshift z 
	//in h km/s / Mpc
	return H_0*E(z);
}

/*! \fn double rho_m(double z)
 *  \brief physical matter density 
 * h^2 Msun/Mpc^3
 * at redshift z */
double Cosmology::rho_m(double z)
{
	//Matter density 
	//Dodelson equation 2.72
	//z is the redshift
	//returns physical h^2 Msun/Mpc^3
	return Omega_m*rho_crit*pow((1.0+z),3);
}

/*! \fn void InitializeDz(void)
 *  \brief Initialize linear growth function lookup tables */
void   Cosmology::InitializeDz(void)
{
	double z_Dz_min;
	double z_Dz_max;

	//initialize_Dz_flag = 1;
	//n_ln_z    = 1000;
	
	//z_Dz_tab = calloc_double_array(n_ln_z);
	//Dz_tab   = calloc_double_array(n_ln_z);

	//tabulate linear growth function
	z_Dz_min = log(1.0e-8);
	z_Dz_max = log(1000.0);
//#pragma omp parallel num_threads(8)
       {
    //   #pragma omp parallel for
	for(int i=0;i<n_ln_z;i++)
	{
		z_Dz_tab[i] = (  (z_Dz_max-z_Dz_min)*((double) i)/((double) n_ln_z-1) + z_Dz_min );
		Dz_tab[i]   = D_growth(exp(z_Dz_tab[i]));

         //cout<<"i="<<i<<"\t"<<"z="<<exp(z_Dz_tab[i]) <<"\t"<<"D(z)="<<Dz_tab[i] <<endl;
	}
      }
	


	//print the linear growth function
	//PrintfDz();
}

/*! \fn void InitializeSigma(void)
*  \brief Initialize the sigma lookup table memory */
void   Cosmology::InitializeSigma(void)
{
	//remember that we've initialized sigma
	initialize_sigma_flag = 1;


	//allocate memory
	//k_EH = calloc_double_array(n_ln_k);
	//R_EH = calloc_double_array(n_ln_k);
	//M_EH = calloc_double_array(n_ln_k);
	//P_k_EH		= calloc_double_array(n_ln_k);
	//Delta2_k_EH	= calloc_double_array(n_ln_k);
	//sigma_EH_tab	= calloc_double_array(n_ln_k);


	//initialize top hat window function

	n_W_x	    = n_ln_k;
	//W2_x	    = calloc_double_array(n_W_x);
	//x_W_x	    = calloc_double_array(n_W_x);
	x_W_x_min   = log(1.0e-8);
	x_W_x_max   = log(12.0);
//#pragma omp parallel num_threads(8)
       {
   //    #pragma omp parallel for

	for(int i=0;i<n_W_x;i++)
	{
		x_W_x[i] = exp(  (x_W_x_max-x_W_x_min)*((double) i)/((double) (n_W_x-1)) + x_W_x_min );
		W2_x[i]  = W2_RTH(x_W_x[i]);
	}
	//acc_W_x    = gsl_interp_accel_alloc();
	//spline_W_x = gsl_spline_alloc(gsl_interp_cspline,n_ln_k);
	//gsl_spline_init(spline_W_x,x_W_x,W2_x,n_W_x);	
}
}

/*! \fn void InitializeSigma_EH(void)
 *  \brief Initialize the sigma for EH98 power spectrum */
void   Cosmology::InitializeSigma_EH(void)
{
	//printf("Initializing sigma (%d)...\n",initialize_sigma_flag);
	//fflush(stdout);
	/*if(initialize_sigma_flag!=1)
		InitializeSigma();

	initialize_sigma_EH_flag = 1;*/
   //cout << "\t Now entering  the loop for tabulation of power spectrum \n \t\t in order to speed-up computation of \n \t dark matter variance in readiness for Halo bias"<<endl;
	//tabulate the power spectrum
//#pragma omp parallel num_threads(8)
  {
    //   #pragma omp parallel for
	for(int i=0;i<n_ln_k;i++)
	{
		//h Mpc^-1
              k_EH[i]= ( (ln_k_max-ln_k_min)*((double) i)/((double) n_ln_k-1) + ln_k_min );
		
		//h^-1 Mpc
		R_EH[i]       = 1./exp(k_EH[i]);


		// h^-3 Mpc^3
		P_k_EH[i]     = P_EH(exp(k_EH[i]),0);
                Tk_EH_tab[i]  = T_EH(exp(k_EH[i]));
		//unit free
		Delta2_k_EH[i] = pow(exp(k_EH[i]),3)*P_k_EH[i]/(2*pi*pi); 

//printf("%d %e %e %e %e\n",i,exp(k_EH[i]),R_EH[i],P_k_EH[i],Delta2_k_EH[i]);
	}

	}


	//renormalize to sigma 8
	sigma_EH_renormalization = sigma_8/sigma_R_EH_Tabulated(8.0);
//cout <<"Sigma_normalization" << sigma_EH_renormalization<<endl;
	//renormalize the power spectrum
//#pragma omp parallel num_threads(8)

       {
 //     #pragma omp parallel for
	for(int i=0;i<n_ln_k;i++)
	{
		Delta2_k_EH[i] *= sigma_EH_renormalization*sigma_EH_renormalization;
		P_k_EH[i]      *= sigma_EH_renormalization*sigma_EH_renormalization;
 ///printf("%d %e %e\n",i,P_k_EH[i],Delta2_k_EH[i]);
	}
}
}	

	//tabulate sigma
	//TabulateSigma_EH();

/*! \fn double LuminosityDistance(double z)
 *  \brief Luminosity distance to redshift z in h^-1 Mpc */
double Cosmology::LuminosityDistance(double z)
{
	//physical h^-1 Mpc
	//Dodelson eq 2.50
	//Checked against Cosmology Calculator
	return ComovingDistance(0,z)*(1+z);
}
/*! \fn double AngularDistance(double z)
 *  \brief Angular distance to redshift z in h^-1 Mpc */
double Cosmology::AngularDistance(double z)
{
	//physical h^-1 Mpc
	//Dodelson eq 2.50
	//Checked against Cosmology Calculator
	return ComovingDistance(0,z)/(1+z);
}

/*! \fn double nJyToAB(double nJy)
 *  \brief Nanojanskys to AB Magnitude, where nJy is in nJy */
double Cosmology::nJyToAB(double nJy)
{
	//returns AB mags
	return 31.4-2.5*log10(nJy);
}

/*! \fn	void PrintfDz(void)
*  \brief Print the linear growth function */
void Cosmology::PrintfDz(void)
{
    double log10_z_min =  -6.0;
    double log10_z_max =  3.0;
    int nz = 1000;
    double z;

    printf("%d\n",nz);
//#pragma omp parallel num_threads(8)
       {
    //   #pragma omp parallel for
    for(int i=0;i<nz;i++)
    {
	z = pow(10, (log10_z_max - log10_z_min)*((double) i)/((double) (nz-1)) + log10_z_min);
	//printf("%e\t%e\n",z,D_growth(z));
    }
}
}

void Cosmology::SetPrelimCosmology(void)
{
    double hh = 0.7;
	sprintf(cosmology_model_name,"Prelim\n");
	//SetCosmology(0.70,0.0,1.,0.04,0,4.15e-5/(0.705*0.705),0.9646,0.810,2.725);
	//SetCosmology(hh,0.0,1.,0.04,0,4.15e-5/(hh*hh),0.9646,0.810,2.725);
	SetCosmology(hh,0.0,1.,0.04,0,0,0.9646,0.810,2.725);
}


/*! \fn void SetWMAP9Cosmology(void)
 *  \brief Set cosmological parameters to the WMAP9 best fit. */
void Cosmology::SetWMAP9Cosmology(void)
{
	//WMAP9 cosmology with h=0.705, Omega_L = 0.728, Omega_m = 0.272,
	//Omega_b = 0.0448468
	//see Hinshaw et al. 2012
	sprintf(cosmology_model_name,"WMAP9\n");
	SetCosmology(0.705,0.728,0.272,0.0448468,0,4.15e-5/(0.705*0.705),0.9646,0.810,2.725);
}

/*! \fn void SetZentner2007Cosmology(void)
 *  \brief Set cosmological parameters to the Zentner 2007 review. */
void Cosmology::SetZentner2007Cosmology(void)
{
	//WMAP9 cosmology with h=0.7, Omega_L = 0.7, Omega_m = 0.3,
	//Omega_b*h^2 = 0.022
	//sigma_8 = 0.93
	sprintf(cosmology_model_name,"Zentner 2007\n");
	SetCosmology(0.7,0.7,0.3,0.022/pow(0.7,2),0,4.15e-5/(0.7*0.7),1.0,0.93,2.725);
}

/*! \fn void SetMillenniumCosmology(void)
 *  \brief Set cosmological parameters to the Millennium Simulation. */
void Cosmology::SetMillenniumCosmology(void)
{
	//h=0.73, Omega_L = 0.75, Omega_m = 0.25,
	//Omega_b = 0.045
	//n_s = 1
	//sigma_8 = 0.9
	sprintf(cosmology_model_name,"Millennium Simulation\n");
	SetCosmology(0.73,0.75,0.25,0.045,0,4.15e-5/(0.7*0.7),1.0,0.9,2.725);
}


void Cosmology::SetPLANCKCosmology(void)
{
	//PLANCK cosmology with h=0.678, Omega_L = 1-0.308, Omega_m = 0.308,
	//Omega_b = 0.0448468,Omega_k_in = 0,Omega_r = 4.15e-5/(0.701*0.701),n_s = 0.9608,sigma_8 = 0.826,T_cmb = 2.728,Omega_HI = 4e-4, b10 = 1.39, b20 = -0.34
	//see Ade et al. 2013
	sprintf(cosmology_model_name,"PLANCK\n");
	SetCosmology(0.678,0.692,0.308,0.0448468,0,4.15e-5/(0.705*0.705),0.9608,0.826,2.725,4e-4,1.39,-0.34, pow(10,9), pow(10,14));
}

/*
void Cosmology::SetPLANCKCosmology(void)
{
	//PLANCK cosmology with h=0.673, Omega_L = 0.728, Omega_m = 0.272,
	//Omega_b = 0.0448468,Omega_k_in = 0,Omega_r = 4.15e-5/(0.701*0.701),n_s = 0.9608,sigma_8 = 0.826,T_cmb = 2.728,Omega_HI = 4e-4, b10 = 1.39, b20 = -0.34
	//see Ade et al. 2013
	sprintf(cosmology_model_name,"PLANCK\n");
	SetCosmology(0.673,0.0,0.95,0.0448468,0,4.15e-5/(0.705*0.705),0.9608,0.826,2.725,4e-4,1.39,-0.34, pow(10,9), pow(10,14));
}
*/


/*! \fn void SetCosmology(double h_in, double Omega_l_in, double Omega_m_in, double Omega_b_in, double Omega_k_in, double Omega_r_in, double n_s_in, double sigma_8_in, double T_cmb_in)
 *  \brief Set cosmological parameters */
void Cosmology::SetCosmology(double h_in, double Omega_l_in, double Omega_m_in, double Omega_b_in, double Omega_k_in, double Omega_r_in, double n_s_in, double sigma_8_in, double T_cmb_in, double Omega_HI_in,double b10_in,double b20_in, double Mmin_in, double Mmax_in)
{
	//default Hubble parameter
	h        = h_in;

	//default Omegas
	Omega_l = Omega_l_in;
	Omega_m = Omega_m_in;
	Omega_b = Omega_b_in;
	Omega_c = Omega_m - Omega_b;
	Omega_k = Omega_k_in;
	Omega_r = Omega_r_in;
        Omega_HI  = Omega_HI_in;

	//spectral slope
	n_s       = n_s_in; 
       // bias parameters
        b10  = b10_in;
        b20  = b20_in;

    ///Mass limit
     Mmin =  Mmin_in;
     Mmax =  Mmax_in;


	//power spectrum normalization
	sigma_8   = sigma_8_in;

	//normalize linear growth factor
	D_growth_0   = 1;
	D_growth_0 = D_growth(0);


	//Transfer function

	//CMB Temperature (default from fixen et al. 1996)
	T_cmb     = T_cmb_in;
	Theta_cmb = T_cmb/2.7;

	//epoch of matter-radiation equality
	z_eq        = z_eq_MatterToRadiation();

	//particle horizon at matter-radiation equality
	k_eq        = k_eq_MatterToRadiation();

	//Epoch when baryons are released from Compton drag
	z_d         = z_d_ComptonDrag();

	//Sound horizon at drag epoch
	s_drag	    = s_DragEpoch();

	//Fitting parameters for EH transfer function
	k_silk	    = k_SilkDamping();

	//EH98 transfer function normalization
	
	//cdm transfer function
	alpha_c_EH   = alpha_c_EH_compute();
	beta_c_EH    = beta_c_EH_compute();

	//baryon transfer function
	alpha_b_EH   = alpha_b_EH_compute();
	beta_b_EH    = beta_b_EH_compute();
	beta_node_EH = beta_node_EH_compute();

	//power spectrum normalization
	delta_H = delta_H_BW();

	//initialize power spectrum
	//InitializeSigma_EH();
       //TabulateSigma_EH();
}

/*! \fn void ShowCosmology(void)
 *  \brief Print cosmology info. */
void Cosmology::ShowCosmology(void)
{
	printf("*************************************\n");
	printf("Model    = %s\n",cosmology_model_name);
	printf("h        = % 5.4e \n",h);
	printf("H_0      = % 5.4e h     km/s/Mpc\n",H_0);
	printf("H_inv    = % 5.4e h^-1  Gyr\n",H_inv);
	printf("rho_crit = % 5.4e h^2   Msun/Mpc^3\n",rho_crit);
	printf("rho_c    = % 5.4e h^2   g/cm^3\n",rho_c);
	printf("\n");
	printf("Omegas  :\n");
	printf("Omega_b = %e\n",Omega_b);
	printf("Omega_m = %e\n",Omega_m);
	printf("Omega_l = %e\n",Omega_l);
	printf("Omega_k = %e\n",Omega_k);
	printf("Omega_r = %e\n",Omega_r);
        printf("\tOmega_HI = %e\n",Omega_HI);
        printf("\t\nBias Parameters :\n");
        printf("\tb10 = %e\n",b10);
        printf("\tb20 = %e\n",b20);
	printf("\nPower Spectrum Info :\n");
	printf("n_s       = %e\n",n_s);
	printf("sigma_8   = %e\n",sigma_8);
	printf("Theta_cmb = %e\n",Theta_cmb);
	printf("z_eq      = %e\n",z_eq);
	printf("k_eq      = %e\n",k_eq);
	printf("z_d       = %e\n",z_d);
	printf("s_d       = %e\n",s_drag);
	printf("\n");
	printf("*************************************\n\n");
}


 double Cosmology::Omega_m_z(double z)
{
     double tmpz = Omega_m*pow((1+z),3)/(pow((Omega_m*pow((1+z),3)+ (1-Omega_m)),1/2));
     return tmpz;
}

double Cosmology::delta_c_z(double z)
{
   double tmp = (3*pow((12 *M_PI),2/3)/20)*(1+0.0123*log10(Omega_m_z(z)));
    //double tmp = 1.686*D_growth_integrated(0.000001)/D_growth_integrated(z);
    return 1.686;
}

/*
double Cosmology::T_b(double z)
{
         //in micro K  I removed 0.001 at the denominator
       return 566.0*h* Omega_HI*pow((1.0+z),2)/(0.003*E(z)*1000);
}
*/

/*! \fn void TabulateSigma_EH(void)
 *  \brief Tabulate sigma for the EH98 power spectrum */
void Cosmology::TabulateSigma_EH(void)
{
	//tabulate M vs. R and sigma vs. R
//#pragma omp parallel private(thread_number)
{
 //#pragma omp for schedule(static) nowait
	for(int i=0;i<n_ln_k;i++)
	{
		//h^-1 Msun in reverse

		//Mass and R_EH are consistent
		//sigma is not
		M_EH[i]   = log10(4.0*pi*rho_m(0)*pow(R_EH[n_ln_k-1-i],3)/3.0);
        PowM_EH[i]  = 4.0*pi*rho_m(0)*pow(R_EH[n_ln_k-1-i],3)/3.0;
		sigma_EH_tab[i] = sigma_R_EH_Tabulated(R_EH[n_ln_k-1-i]);
//cout<<"\t Mass="<<"\t"<< PowM_EH[i]<<"\tlog10 of M"<<"\t"<<M_EH[i]<<"\t R"<<"\t"<< R_EH[i]<<endl;
//printf("%d %e %e %e\n",i,M_EH[i] ,R_EH[i],sigma_EH_tab[i]);
	}
}
	
}
double Cosmology::UniversalAge(double a)
{
	//universal age
	//in h^-1 Gyr

	double fp[4];
	double da_int;
	double eps = 1.0e-9;
	double a_min = 1.0e-14;

	if(a<=a_min)
		return 0;

	da_int = 0.01*fabs(log(a)-log(a_min));

	fp[0] = Omega_m;
	fp[1] = Omega_l;
	fp[2] = Omega_k;
	fp[3] = Omega_r;

	//in h^-1 Gyr
	return H_inv*integrate(UniversalAgeIntegrand,fp,4,log(a),log(a_min),da_int,eps);
}
double UniversalAgeIntegrand(double a, void *params)
{
	double *fp = (double *) params;

	double Omega_m = fp[0];
	double Omega_l = fp[1];
	double Omega_k = fp[2];
	double Omega_r = fp[3];

	a = exp(a);
	return -a/(a*sqrt( Omega_m*pow(a,-3) + Omega_k*pow(a,-2) + Omega_l + Omega_r*pow(a,-4)));
}

double Cosmology::ConformalTime(double z)
{
	//conformal time 
	//Dodelson equation 2.90
	//in h^-1 Gyr

	double fp[4];
	double da_int;
	double eps = 1.0e-9;
	double a_min = 1.0e-11;
	double a = 1./(1+z);
	da_int = 0.01*(log(a)-log(a_min));
	//da_int = 0.01*(log(a) - log(a_min));
	//dt_int = 0.01*(log(t) - log(t_min));

	fp[0] = Omega_m;
	fp[1] = Omega_l;
	fp[2] = Omega_k;
	fp[3] = Omega_r;

	//in h^-1 Gyr
	return H_inv*integrate(ConformalTimeIntegrand,fp,4,log(a_min),log(a),da_int,eps);
}
double ConformalTimeIntegrand(double lna, void *params)
{
	double *fp = (double *) params;
	double a = exp(lna);
	double Omega_m = fp[0];
	double Omega_l = fp[1];
	double Omega_k = fp[2];
	double Omega_r = fp[3];
	return 1./(a*sqrt( Omega_m*pow(a,-3) + Omega_k*pow(a,-2) + Omega_l + Omega_r*pow(a,-4)));
}

//Power spectrum routines here

/*! \fn double MatterToRadiationRatio(void)
 *  \brief proportion of matter to radiation
 *  EH98 section 2 */
double Cosmology::MatterToRadiationRatio(void)
{
	//Section 2 of Eisenstein & Hu 1998
	return Omega_m*h*h*pow(Theta_cmb,-4);
}
/*! \fn double BaryonToRadiationRatio(void)
 *  \brief proportion of baryon to radiation
 *  EH98 section 2 */
double Cosmology::BaryonToRadiationRatio(void)
{
	//Section 2 of Eisenstein & Hu 1998
	return Omega_b*h*h*pow(Theta_cmb,-4);
}


/*! \fn double z_eq_MatterToRadiation(void)
 *  \brief Redshift of matter-radiation equality 
 *  EH98 eqn 2 */
double Cosmology::z_eq_MatterToRadiation(void)
{
	//redshift of matter-radiation equality
	//Equation 2 of Eisenstein and Hu 1998
	return 2.50e4 * MatterToRadiationRatio();
}
/*! \fn double k_eq_MatterToRadiation(void)
 *  \brief Particle horizon at matter-radiation equality */
double Cosmology::k_eq_MatterToRadiation(void)
{
	//scale of particle horizon at the
	//redshift of matter-radiation equality
	//in Mpc^-1

	return sqrt(2*Omega_m*(h*H_0/(c/1.0e3))*(h*H_0/(c/1.0e3))*z_eq_MatterToRadiation());

	//approximate
	//return 7.46e-2 * Omega_m*h*h/(Theta_cmb*Theta_cmb);
}
/*! \fn double z_d_ComptonDrag(void)
 *  \brief Epoch when baryons are released from Compton drag */
double Cosmology::z_d_ComptonDrag(void)
{
	//epoch when baryons are released
	//from Compton Drag
	//eqn 4 of Eisenstein and Hu 1998
	//see also Hu and Sugiyama 1996
	double ohh = Omega_m*h*h;
	double b1 = 0.313*pow(ohh, -0.419)*(1 + 0.607*pow(ohh,0.674));
	double b2 = 0.238*pow(ohh,  0.223);
	double A;
	double B;

	A = 1291. *pow(ohh,0.251)/(1 + 0.659*pow(ohh,0.828));
	B = (1+b1*pow(Omega_b*h*h,b2));
	return A*B;
}
/*! \fn double R_BaryonToPhotonMomentumDensity(double z)
 *  \brief Ratio of Baryon to Photon Momentum Density */
double Cosmology::R_BaryonToPhotonMomentumDensity(double z)
{
	//eqn 5 of EH98
	return 31.5*BaryonToRadiationRatio()*(1.0e3/z);
}

/*! \fn double s_DragEpoch(void)
 *  \brief Sound horizon at the drag epoch */
double Cosmology::s_DragEpoch(void)
{
	//sound horizon at the drag epoch
	//in Mpc
	//eq 6 of eh98
	double R_eq = R_BaryonToPhotonMomentumDensity(z_eq);
	double R_d  = R_BaryonToPhotonMomentumDensity(z_d);
	return (2./(3.*k_eq))*sqrt(6./R_eq)*log( (sqrt(1.+R_d)+sqrt(R_d+R_eq)) / (1.+sqrt(R_eq)) );
}

/*! \fn double delta_H_BW(void)
 *  \brief Set the power spectrum normalization */
double Cosmology::delta_H_BW(void)
{
	//power spectrum normalization

	double ntilde = n_s-1;
	double dH;
	if( Omega_l ==0)
	{
		//for open from Bunn & White 1997, eqn A2 of Eisentstein & Hu 1998
		dH = 1.95e-5 * pow(Omega_m,-0.35 -0.19*log(Omega_m)-0.17*ntilde)*exp(-ntilde-0.14*ntilde*ntilde);
	}else{
		//for flat from Bunn & White 1997, eqn A2 of Eisentstein & Hu 1998
		dH = 1.94e-5 * pow(Omega_m,-0.785-0.05*log(Omega_m))*exp(-0.95*ntilde-0.169*ntilde*ntilde);
	}
	return dH;
}

/*! \fn double k_SilkDamping(void)
 *  \brief Scale where damping due to diffusion
 *	    of photons past baryons becomes important */
double Cosmology::k_SilkDamping(void)
{
	//k_silk (silk 1968)
	//damping due to diffusion
	//of photons past baryons
	//eqn 7 of EH98
	//in Mpc^-1
	double ks = 1.6*pow(Omega_b*h*h,0.52)*pow(Omega_m*h*h,0.73)*(1 + pow(10.4*Omega_m*h*h,-0.95)); 
	return ks;
}

/*! \fn double alpha_c_EH_compute(void)
 *  \brief Fitting parameter for Eisenstein & Hu
 *  transfer function */
double Cosmology::alpha_c_EH_compute(void)
{
	//eqn 11 of EH98
	double fb  = Omega_b/Omega_m;
	double ohh = Omega_m*h*h;
	double a1  = pow(46.9*ohh,0.670)*(1 + pow(32.1*ohh, -0.532));
	double a2  = pow(12.0*ohh,0.424)*(1 + pow(45.0*ohh, -0.582));
	return pow(a1,-fb)*pow(a2,-pow(fb,3));
}
/*! \fn double beta_c_EH_compute(void)
 *  \brief Fitting parameter for Eisenstein & Hu
 *  transfer function */
double Cosmology::beta_c_EH_compute(void)
{
	//eqn 12 of EH98
	double b1 = 0.944/(1+pow(458*Omega_m*h*h,-0.708));
	double b2 = pow(0.395*Omega_m*h*h,-0.0266);
	return 1./(1+ b1*( pow(Omega_c/Omega_m, b2) -1 ));
}


/*! \fn double alpha_b_EH_compute(void)
 *  \brief Baryon transfer function amplitude at ks>>1 */
double Cosmology::alpha_b_EH_compute(void)
{
	//eqn 14 of EH98
	double R_d  = R_BaryonToPhotonMomentumDensity(z_d);

	return 2.07*k_eq*s_drag*pow(1+R_d,-0.75)*alpha_b_EH_G( (1+z_eq)/(1+z_d) );
}
/*! \fn double alpha_b_EH_G(double y)
 *  \brief Baryon transfer function functionality at ks>>1 */
double Cosmology::alpha_b_EH_G(double y)
{
	//eqn 15 of EH98
	return y*(-6*sqrt(1+y) + (2 + 3*y)*log( ( sqrt(1+y) +1 ) / ( sqrt(1+y) - 1) ) );
}
/*! \fn double beta_b_EH_compute(void)
 *  \brief Fractional scale where velocity contributions to
 *  the baryon transfer function become important */
double Cosmology::beta_b_EH_compute(void)
{
	//eqn 24 of EH98
	//checked against wayne
	return 0.5 + (Omega_b/Omega_m) + (3 - 2*(Omega_b/Omega_m))*sqrt( pow(17.2*Omega_m*h*h,2) + 1);
}
/*! \fn double s_tilde_EH(double k)
 *  \brief Phenomenological shift in baryon acoustic scale */
double Cosmology::s_tilde_EH(double k)
{
	//eqn 22 of eh98
	return s_drag/pow( 1 + pow(beta_node_EH/(k*s_drag),3) , 0.3333);
}
/*! \fn double beta_node_EH_compute(void)
 *  \brief Shift in acoustic nodes to higher k*/
double Cosmology::beta_node_EH_compute(void)
{
	//eqn 23 of EH98
	//checked against wayne's code
	return 8.41*pow(Omega_m*h*h, 0.435);
}

/*! \fn double T_EH(double k)
 *  \brief EH98 Transfer function */
double Cosmology::T_EH(double k)
{
	//EH98 transfer function fitting, eq 16
	//k is in h Mpc^-1
	//note conversion to Mpc^-1

	return (Omega_b/Omega_m)*T_b_EH(k*h) + (Omega_c/Omega_m)*T_c_EH(k*h);

}

/*! \fn double T_0_tilde_EH(double k, double alpha, double beta)
 *  \brief Base line EH transfer function */
double Cosmology::T_0_tilde_EH(double k, double alpha, double beta)
{
	//eqn 19 of eh 98
	//note that k is in Mpc

	//eqn 9 of eh 98
	double q = k/(13.41*k_eq);
	
	//eqn 20 of eh 98
	double C = (14.2/alpha) + (386./(1+69.9*pow(q,1.08)));

	//eqn 19 of eh98
	double num = log(e + 1.8*beta*q);

	return num / ( num + C*q*q);
}

/*! \fn double T_b_EH(double k)
 *  \brief Baryonic contribution of the EH transfer function */
double Cosmology::T_b_EH(double k)
{
	//k is in Mpc
	//eqn 21 of EH98
	double stilde = s_tilde_EH(k);
	double j_0    = sin(k*stilde)/(k*stilde);
	double xx     = k*s_drag;
	return j_0*( T_0_tilde_EH(k,1,1)/(1+pow(xx/5.2,2))  + alpha_b_EH*exp(-pow(k/k_silk,1.4))/(1+pow(beta_b_EH/xx,3)));
}
/*! \fn double T_b_EH(double k)
 *  \brief CDM contribution of the EH transfer function */
double Cosmology::T_c_EH(double k)
{
	//k is in Mpc
	double f = 1./( 1. + pow(k*s_drag/5.4, 4));

	return f*T_0_tilde_EH(k,1.,beta_c_EH) + (1-f)*T_0_tilde_EH(k,alpha_c_EH,beta_c_EH);
}

/*! \fn double P_EH(double k, double z)
 *  \brief EH98 cdm+baryons power spectrum */
double Cosmology::P_EH(double k, double z)
{
	//EH98 cdm+baryons power spectrum
	//dodelson eqn 7.9 + 7.71
	//k is in h Mpc^-1
	double A = 2.*pi*pi*delta_H*delta_H*pow(k,n_s)/pow(H_0/(c/1.0e3),n_s+3);
	double Tk;
	double Dz;

	if(k==0)	
		return 0;
	Tk = T_EH(k);
	Dz = D_growth(z);
	return A*pow(Tk*Dz,2)*pow(sigma_EH_renormalization,2);
}

/*! \fn double Delta2_EH(double k, double z)
 *  \brief Unit-free CDM+baryon power spectum */
double Cosmology::Delta2_EH(double k, double z)
{
	//EH98 unit-free cdm+baryon power spectrum
	//dodelson eqn 7.9 + 7.71 * k^3 / 2pi
	//k is in h Mpc^-1
	return pow(k,3)*P_EH(k,z)/(2*pi*pi);
}

/*! \fn double sigma(double M)
 *  \brief RMS density fluctuations on a smoothing
 *  mass scale M */
double Cosmology::sigma(double M)
{
	//RMS density flucuation of 
	//smoothing mass scale M
	return sigma_EH(M);
}

/*! \fn double sigma_EH(double M)
 *  \brief RMS density fluctuations on a smoothing
 *  mass scale M for the EH98 power spectrum*/
double Cosmology::sigma_EH(double M)
{
	//RMS density flucuation of 
	//smoothing mass scale M
	//if(initialize_sigma_EH_flag!=1)
	//	InitializeSigma_EH();
	return sigma_M_EH_Tabulated(M);
}

/*! \fn void SetSigmaMemoryToNull(void)
*  \brief Set pointers associated with sigma
*  lookup tables to NULL. */
void Cosmology::SetSigmaMemoryToNull(void)
{
	/*
      Delta2_k_EH    = NULL;	
	P_k_EH         = NULL;	
	k_EH           = NULL;
	R_EH           = NULL;
	M_EH           = NULL;
	sigma_EH_tab   = NULL;
*/
}

/*! \fn void FreeSigmaMemory(void)
 *  \brief Free memory associated with sigma
 *  lookup tables. */
void Cosmology::FreeSigmaMemory(void)
{
	//free lookup table memory
	/*
       delete[] k_EH;
	delete[] R_EH;
	delete[] M_EH;
	delete[] P_k_EH;
	delete[] sigma_EH_tab;
	delete[] Delta2_k_EH;
*/

	//set pointers to null
	SetSigmaMemoryToNull();
}

/*! \fn void Print_sigma_M_EH(void)
 *  \brief Print sigma as a function of M */
void Cosmology::Print_sigma_M_EH(void)
{
	/*if(initialize_sigma_EH_flag!=1)
		InitializeSigma_EH();
	printf("%d\n",n_ln_k);
	for(int i=n_ln_k-1;i>=0;i--)*/
	{
		//printf("%5d  %e  %e  %e\n",i,sigma_R_EH_Tabulated(R_EH[i]),pow(10,M_EH[n_ln_k-1-i]),R_EH[i]);

	}
}
/*! \fn void Print_Delta2_k_EH(void)
 *  \brief  Print Delta2 as a function of k */
void   Cosmology::Print_Delta2_k_EH(void)
{
	double k, Delta;
	/*if(initialize_sigma_EH_flag!=1)
		InitializeSigma_EH();
	printf("%d\n",n_ln_k);*/
	for(int i=0;i<n_ln_k;i++)
	{
		//k = exp(k_EH[i]);
		//Delta = Delta2_k_EH_Tabulated(k);
		k = exp(k_EH[i]);
		Delta = Delta2_k_EH[i];
		//printf("%e %e\n",k,Delta);
	}
}


/////////////////////////////////////////
// Member public functions for halo mass functions
/////////////////////////////////////////
/*! \fn double nu(double M, double z)
*  \brief Peak rarity*/
/*
double Cosmology::nu(double M, double z)
{
  //return peak rarity
  return nu_Mz(M,&z);
}*/

/*! \fn double dndm_sheth_tormen(double M, double z)
*  \brief Press-Schechter Halo Mass Function at redshift z */
/*
double Cosmology::dndm_press_schechter(double M, double z)
{
  //mass is in h^-1 Msun
  double fp[2];
 
  //store redshift
  fp[0] = z; 

  //store matter density in h^2 Msun Mpc^-3
  fp[1] = rho_m(z);

  //return dndM for Press-Schechter
  return dndm_PS(M,fp);
}*/


/*! \fn double dndm_sheth_tormen(double M, double z)
*  \brief Sheth-Tormen Halo Mass Function at redshift z */
/*
double Cosmology::dndm_sheth_tormen(double M, double z)
{
  //mass is in h^-1 Msun
  double fp[2];
 
  //store redshift
  fp[0] = z; 

  //store matter density in h^2 Msun Mpc^-3
  fp[1] = rho_m(z);

  //return dndM for Sheth-Tormen
  return dndm_ST(M,fp);
}*/

/////////////////////////////////////////
// Member public functions for halo bias functions
/////////////////////////////////////////
/*! \fn double b_MW(double M, double z)
*  \brief Mo & White 1996 Bias Function */
/*
double Cosmology::b_mo_white(double M, double z)
{
  //mass is in h^-1 Msun
  double fp[1];
 
  //store redshift
  fp[0] = z; 

  //return bias from Mo & White 1996
  return b_MW(M,fp);
}*/

/*! \fn double b_tinker(double M, double z)
*  \brief Tinker et al. 2010 Bias Function */
/*
double Cosmology::b_tinker(double M, double z)
{
  //mass is in h^-1 Msun
  double fp[1];
 
  //store redshift
  fp[0] = z; 

  //return bias from Tinker et al. 2010
  return b_T(M,fp);
}*/



/////////////////////////////////////////
// Non member public functions below here
/////////////////////////////////////////

/*! \fn double ComovingDistanceIntegrand(double z, void *params)
 *  \brief Comoving distance integrand.*/
double ComovingDistanceIntegrand(double z, void *params)
{
    double *fp = (double *) params;

    double Omega_m = fp[0];
    double Omega_l = fp[1];
    double Omega_k = fp[2];
    double Omega_r = fp[3];

    //Eqn 15 of Hogg 1999
    // dz / E(z)

    return 1./sqrt( Omega_m*pow(1.+z,3) + Omega_k*pow(1+z,2) + Omega_l + Omega_r*pow(1+z,4));
}



/*! \fn double D_growth_integrand(double lnz, void *params)
 *  \brief Integrand for the linear growth factor */
double D_growth_integrand(double lnz, void *params)
{
	//integrand of D(z)
	double *fp = (double *) params;
	double z = exp(lnz);
	double Omega_m = fp[0];
	double Omega_l = fp[1];
	double Omega_k = fp[2];
	double Omega_r = fp[3];
	double H = sqrt(Omega_m*pow(1+z,3) + Omega_k*pow(1+z,2) + Omega_l + Omega_r*pow(1+z,4));

	return z*(1+z)/pow(H,3);
}


/*! \fn double D_growth_Tabulated(double z)
 *  \brief Tabulated linear growth factor */
double Cosmology::D_growth_Tabulated(double z)
{
	double Dz;

	/*if(Dz_tab==NULL)
	{
		printf("Must allocate Dz first before using D_growth_Tabulated!\n");
		exit(-1);
	}*/
	if(z<1.0e-9)
		return 1.0;
	if(log(z)<z_Dz_tab[0])
		return Dz_tab[0];
	if(log(z)>z_Dz_tab[n_ln_z-1])
		return Dz_tab[n_ln_z-1];


	//Dz = gsl_spline_eval(spline_D_z,log(z),acc_D_z);
       Spline<double, double> CubicSpline_D_z(z_Dz_tab, Dz_tab);
       Dz =  CubicSpline_D_z.interpolate(log(z));

	return Dz;
}



/*! \fn double sigma2_R_EH_integrand(double lnk, void *params)
 *  \brief k-space integrand for RMS density fluctuations on a size scale R */
double sigma2_R_EH_integrand(double lnk, void *params)
{
	double *fp = (double *) params;
	double k = exp(lnk);
	double R = fp[0];
	double x = R*k;

	//return W2_RTH_Tabulated(x)*Delta2_k_EH_Tabulated(k);

#ifdef KTOPHAT
	return W_KTH(k,R)*W_KTH(k,R)*Delta2_k_EH_Tabulated(k);
#else
	return W2_RTH(x)*Delta2_k_EH_Tabulated(k);
#endif

}

/*! \fn double Delta2_k_EH_Tabulated(double k)
 *  \brief Dimensionless smoothed power spectrum for EH98 */
double Delta2_k_EH_Tabulated(double k)
{
	double Dk;

	/*if(Delta2_k_EH==NULL)
	{
		printf("Must allocate memory first before using D_k_EH_Tabulated!\n");
		exit(-1);
	}*/
	if(log(k)<k_EH[0])
		return Delta2_k_EH[0];
	if(log(k)>k_EH[n_ln_k-1])
		return Delta2_k_EH[n_ln_k-1];

         Spline<double, double> CubicSpline_Delta2_k_EH(k_EH, Delta2_k_EH);
         Dk = CubicSpline_Delta2_k_EH.interpolate(log(k));
	//Dk = gsl_spline_eval(spline_sigma_EH,log(k),acc_sigma_EH);

	return Dk;
}

/*! \fn double P_k_EH_Tabulated(double k)
 *  \brief P(k) tabulated for EH08 power spectrum */
double Cosmology::P_k_EH_Tabulated(double k)
{
	double Pk;

	/*if(P_k_EH==NULL)
	{
		printf("Must allocate memory first before using P_k_EH_Tabulated!\n");
		exit(-1);
	}*/
	if(log(k)<k_EH[0])
		return P_k_EH[0];
	if(log(k)>k_EH[n_ln_k-1])
		return P_k_EH[n_ln_k-1];

        Spline<double, double> CubicSpline_P_EH(k_EH, P_k_EH);
        Pk = CubicSpline_P_EH.interpolate(log(k));
	//Pk = gsl_spline_eval(spline_P_EH,log(k),acc_P_EH);

	return Pk;
}


double Cosmology::Tk_tab(double k)
{
return Tk_tab_in(k);
}

double Tk_tab_in(double k)
{
       double Tk;
         

    if(log(k)<k_EH[0])
		return Tk_EH_tab[0];
	if(log(k)>k_EH[n_ln_k-1])
		return Tk_EH_tab[n_ln_k-1];	

          Spline<double, double> CubicSpline_T_EH(k_EH, Tk_EH_tab);
       Tk = CubicSpline_T_EH.interpolate(log(k));
      
	return Tk;
         
}


/*! \fn double sigma_k_EH_Tabulated(double k) 
 *  \brief sigma(k) tabulated for EH08 power spectrum */
double sigma_k_EH_Tabulated(double k)
{
	//k is in h Mpc^-1
	return sigma_R_EH_Tabulated(1./k);
}
//double Renormsigma_EK(double k)
/*! \fn double W2_RTH(double x)
 *  \brief Square of the Fourier transform of the top hat window function */
double W2_RTH(double x)
{
	//fourier transform of the top hat
	//window function
	if(x<1.0e-4)
		return 1.0;

	//note that W(Rk) is sensitive to this upper limit.
	if(x>5000.)
		return 0.0;
	return pow(3/(x*x*x),2)*pow(sin(x)-x*cos(x),2);
}

/*! \fn	double sigma_M_EH_Tabulated(double M)
 *  \brief sigma as a function of M for EH98 power spectrum */
//Cosmology::
double sigma_M_EH_Tabulated(double M)
{
	//M is in h^-1 Msun

	double sigma;

	/*if(sigma_EH_tab==NULL)
	{
		printf("Must allocate memory first before using sigma_M_EH_Tabulated!\n");
		exit(-1);
	}*/
	if(M<pow(10,M_EH[0]))
		return sigma_EH_tab[0];
	if(M>pow(10,M_EH[n_ln_k-1]))
		return sigma_EH_tab[n_ln_k-1];

        Spline<double, double> CubicSpline_sigma_M_EH(M_EH, sigma_EH_tab);
        sigma = CubicSpline_sigma_M_EH.interpolate(log10(M));
	//sigma = gsl_spline_eval(spline_sigma_M_EH,log10(M),acc_sigma_M_EH);


	return sigma;
	
}


/*  mass scale M for the EH98 power spectrum*/
double Cosmology::sigma_R(double R)
{
	//RMS density flucuation of 
	//smoothing mass scale M
	//if(initialize_sigma_EH_flag!=1)
	//	InitializeSigma_EH();
	return sigma_R_EH_Tabulated(R);
}




/*! \fn double sigma_R_EH_Tabulated(double R) 
 *  \brief RMS density fluctuations on a size scale R */
double sigma_R_EH_Tabulated(double R)
{
	//R is in h^-1 Mpc
       const double EPSREL = 1e-2;
	double sigma2;
	double fp[1];

	double ln_k_min = log((double) DLNKMIN); 
	double ln_k_max = log((double) DLNKMAX);
	double dinit;
	double eps = 1.0e-7;

	dinit = 0.01*(ln_k_max-ln_k_min);

	fp[0] = R;

	

	sigma2 = integrate(sigma2_R_EH_integrand,fp,1,ln_k_min,ln_k_max,dinit,eps);

 //sigma2 = Integrate(bind(sigma2_R_EH_integrand,R,_1),ln_k_min,ln_k_max, EPSREL);
//cout <<"Sigma Integration ="<<sigma2<<endl;
	return sqrt(sigma2);
	
}

///////////////////////////////////////
//peak rarity below here
///////////////////////////////////////

/*! \fn double nu_Mz(double M, void *params)
 *  \brief Peak rarity nu = delta_c / [D(z)*sigma(M)] */
/*
double nu_Mz(double M, void *params)
{
  double *fp     = (double *) params;
  double z       = fp[0];
  double delta_c = 1.686;
  double Dz      = D_growth_Tabulated(z);
  double sig     = sigma_M_EH_Tabulated(M);

  //printf("M %e z %e Dz %e sig %e\n",M,z,Dz,sig);

  //nu = delta_c / sigma(M,z)
  return delta_c / (Dz*sig);
}*/

/*! \fn double dnudM(double M, void *params)
 *  \brief Derivative of nu(M) wrt M */
/*
double dnudM(double M, void *params)
{
  //nu = delta_c/[D(z) * sigma(M)]
  //dnudM = -delta_c/[D(z) * sigma(M)^2] x dsigmadM
  //dnudM = -[nu/sigma(M)] x dsigmadM
  return -nu_Mz(M,params)*dsigmadM(M)/sigma_M_EH_Tabulated(M);
}*/

/*! \fn double dsigmadM(double M)
 *  \brief derivative of sigma(M) wrt M */
/*
double dsigmadM(double M)
{
	//dumb derivative
	double sigma_A = sigma_M_EH_Tabulated(1.01*M);
	double sigma_B = sigma_M_EH_Tabulated(0.99*M);
	return (sigma_A-sigma_B)/(0.02*M);
}*/

///////////////////////////////////////
// first crossing distributions below here
///////////////////////////////////////

/*! \fn double f_nu_PS(double nu, void *params)
 *  \brief First-Crossing Distribution for Press-Schechter */
/*
double f_nu_PS(double nu, void *params)
{
  //equation 9 of Mo & White 2002
  return sqrt(2/M_PI)*exp(-0.5*nu*nu);
}*/


/*! \fn double f_nu_ST(double nu, void *params)
 *  \brief First-Crossing Distribution for Sheth-Tormen*/
/*
double f_nu_ST(double nu, void *params)
{
  //equation 14 of Mo & White 2002
  double a   = 0.707;
  double nup = sqrt(a)*nu;
  double A   = 0.322;
  double q   = 0.3;
  return A*(1 + pow(nup,-2*q))*sqrt(2/M_PI)*exp(-0.5*nup*nup);
}*/

///////////////////////////////////////
// Halo mass functions below here
///////////////////////////////////////

/*! \fn double dndm_PS(double M, void *params)
 *  \brief Halo Mass Function for Press-Schecher*/
/*
double dndm_PS(double M, void *params)
{
  //M is in h^-1 Msun
  double *fp    = (double *) params;
  double z      = fp[0]; //redshift
  double rho_m  = fp[1]; //mean background density in h^2 Msun / Mpc^3
  double nu     = nu_Mz(M,params);
  double dnu_dM = dnudM(M,params);

  //number of halos per h^3 Mpc^-3 per dM
  return (rho_m/M)*dnu_dM*f_nu_PS(nu,params);
}*/

/*! \fn double dndm_ST(double M, void *params)
 *  \brief Halo Mass Function for Sheth-Tormen*/
/*
double dndm_ST(double M, void *params)
{
  //M is in h^-1 Msun
  double *fp    = (double *) params;
  double z      = fp[0]; //redshift
  double rho_m  = fp[1]; //mean background density in h^2 Msun / Mpc^3
  double nu     = nu_Mz(M,params);
  double dnu_dM = dnudM(M,params);

  //printf("M %e z %e rho_m %e nu %e dnu_dM %e\n",M,z,rho_m,nu,dnu_dM);

  //number of halos per h^3 Mpc^-3 per dM
  return (rho_m/M)*dnu_dM*f_nu_ST(nu,params);
}*/

/*! \fn double b_MW(double M, void *params)
 *  \brief Halo Bias Function for Press-Schecher*/
/*
double b_MW(double M, void *params)
{
  //M is in h^-1 Msun
  double *fp    = (double *) params;
  double nu     = nu_Mz(M,params);
  double delta_c = 1.686;

  //halo bias
  return 1.0 + (nu*nu - 1.0)/delta_c;
}*/

/*! \fn double b_T(double M, void *params)
 *  \brief Halo Bias Function for Tinker mass function*/
/*
double b_T(double M, void *params)
{
  //M is in h^-1 Msun
  double *fp    = (double *) params;
  double z      = fp[0];
  double nu     = nu_Mz(M,params);
  double delta_c = 1.686;
  double alpha_0 = 0.368;
  double beta_0  = 0.589;
  double gamma_0 = 0.864;
  double phi_0   = -0.729;
  double eta_0   = -0.243;
  double beta = beta_0*pow(1+z,0.2);
  double phi  = phi_0*pow(1+z,-0.08);
  double eta  = eta_0*pow(1+z,0.27);
  double gamma = gamma_0*pow(1+z,-0.01);

  //halo bias
  return 1.0 + ( gamma*nu*nu - (1+2*eta) + 2*phi/(1+pow(beta*nu,2*phi)))/delta_c;
}
*/

 

/*
double Cosmology::M_EH_sigma_Tabulated(double sig)
{
	//M is in h^-1 Msun

	double Mass;
        Spline<double, double> CubicSpline_Sigma_EH(sigma_EH_tab,PowM_EH);
        Mass = CubicSpline_Sigma_EH.interpolate(sig);
	return Mass;	
}*/











