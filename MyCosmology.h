/*! \file cosmology.h
 *  \brief Declaration for the Cosmology class that contains
 *         useful cosmological calculations.
 */
#ifndef  _COSMOLOGY
#define  _COSMOLOGY
//#include <vector>
#include<gsl/gsl_spline.h>
#include"Copter_Jeremy/constants.h"

/*! \class Cosmology
 *  \brief Class containing useful cosmology calculations. */
class Cosmology: public Constants
{
    public:
        //name of the cosmological model
        char cosmology_model_name[200];

        //Hubble constant at z=0 in h km/s/Mpc
        double H_0;

        //Hubble parameter h
        double h;

        //inv Hubble constant at z=0 in h^-1 Gyr
        double H_inv;

        //critical density in h^2 Msun / Mpc^3
        double rho_crit;

        //critical density in h^2 g / cm^3
        double rho_c;

        //Omegas

        //Omega_matter at z=0
        double Omega_m;
 
        //Omega_CDM at z=0
        double Omega_c;

        //Omega_baryon at z=0
        double Omega_b;

        //Omega_dark_energy at z=0
        double Omega_l;

        //Omega_curvature at z=0
        double Omega_k;

        //Omega_radiation at z=0
        double Omega_r;


        // Omega HI
        double Omega_HI;
	//linear growth factor at z=0
	double D_growth_0;

	//baryon fraction
	double f_baryon;

	//spectral index
	double n_s;

	//rms fluctuations on 8 h^-1 Mpc scales
	double sigma_8;

	//CMB temperature
	double T_cmb;
	double Theta_cmb;

	//linear spherical collapse overdensity
	double delta_c;

	//number of samples in tabulated power spectra
	//int n_ln_k;

	//largest scale for tabulated power spectrum 
	double ln_k_min;

	//smallest scale for tabulated power spectrum 
	double ln_k_max;

	//smallest R*k in W(R,k) table
	double x_W_x_min;

	//largest R*k in W(R,k) table
	double x_W_x_max;
        // Bias world
        // linear bias
        double b10;
        // Non-linear bias
        double b20;
    /// Lower Limit of the halo mass to host HI
    double Mmin;
  /// upper Limit of the halo mass to host HI
   double Mmax;
 
       


        //constructor
        Cosmology(void);

	/*! \fn void InitializeCosmology(void)
	*  \brief Initialize the cosmology class */
	void InitializeCosmology(void);

	/*! \fn void InitializeDz(void)
	*  \brief Initialize linear growth function lookup tables */
	 void InitializeDz(void);
          //void InitializeSigma_EH(void);
          //void TabulateSigma_EH(void);
        //destructor
        ~Cosmology(void);
        
        //print cosmology
        void ShowCosmology(void);



        //initialize cosmology

	//void SetCosmology(double h_in = 0.701, double Omega_l_in = 0.721, double Omega_m_in = (1-0.721), double Omega_b_in = 0.0462, double Omega_k_in = 0, double Omega_r_in = 4.15e-5/(0.701*0.701),double n_s_in = 0.960, double sigma_8_in = 0.817, double T_cmb_in = 2.728);

void SetCosmology(double h_in = 0.678, double Omega_l_in = 1-0.308, double Omega_m_in = 0.308, double Omega_b_in = 0.0462, double Omega_k_in = 0, double Omega_r_in = 4.15e-5/(0.701*0.701),double n_s_in = 0.960, double sigma_8_in = 0.817, double T_cmb_in = 2.728, double Omega_HI_in = 4e-4,double b10_in = 1.39,double b20_in = -0.34, double Mmin = pow(10,9), double Mmax = pow(10,14));		

	// Omega_L = 0.728, Omega_m=0.272, Omega_b=0.0448, h=0.705, sigma_8 = 0.810, n=0.9646
        void SetWMAP9Cosmology(void);
        void SetPrelimCosmology(void);
        void SetZentner2007Cosmology(void);
        void SetMillenniumCosmology(void);
       void SetPLANCKCosmology(void);


	//Cosmology member functions below here


	/*! \fn double ABTonJy(double AB)
 	*  \brief AB Magnitude to Nanojansky, where nJy is in nJy */
	double ABTonJy(double AB);

	/*! \fn double ABToFlux(double AB)
 	*  \brief AB Magnitude to Flux, where F_nu is in ergs/s/cm^2/Hz*/
	double ABToFlux(double AB);

	//linear growth function
	double D_growth(double z);
	double D_growth_integrated(double z);

	//distance modulus
	//to redshift z in delta mags for a flat F_nu source
	double DistanceModulus(double z);
    
	//unit-free Hubble parameter H(z)/H_0 at
	//redshift z
	double E(double z);

	/*! \fn double rho_m(double z)
	*  \brief physical matter density 
	* h^2 Msun/Mpc^3
	* at redshift z */
	double rho_m(double z);

	/*! \fn double FluxToAB(double F_nu)
 	*  \brief Flux to AB Magnitude, where F_nu is in erg/s/cm^2/Hz */
	double FluxToAB(double F_nu);

        //Comoving Distance in h^-1 Mpc
        //between redshift z_min and z_max
        double ComovingDistance(double z_min, double z_max);

        //Comoving Volume in h^-3 Mpc^3
        //between redshift z_min and z_max
        double ComovingVolume(double z_min, double z_max);

	//Hubble parameter at
	//redshift z
	double H(double z);

        //Luminosity Distance
        //to redshift in h^-1 Mpc
	double LuminosityDistance(double z);

        //Angular Distance
        //to redshift in h^-1 Mpc
	double AngularDistance(double z);



 

     ///For collapse model
        double Omega_m_z(double z);
        double delta_c_z(double z);

     ///HI Brightness temperature
        //double  T_b(double z);

	//initialization flags below here

	// initialize linear growth lookup table
	int initialize_Dz_flag;

	/*! \fn double nJyToAB(double nJy)
 	*  \brief Nanojanskys to AB Magnitude, where nJy is in nJy */
	double nJyToAB(double nJy);

	//printing routines below here
	/*! \fn	void PrintfDz(void)
	*  \brief Print the linear growth function */
	void PrintfDz(void);




	//Routines for calculating the power spectrum

	/*! \fn double delta_H_BW(void);
	*  \brief power spectrum normalization
	*  from Bunn & White 1997 */
	double delta_H_BW(void);
	double delta_H;

	//matter to photon ratio
	double MatterToRadiationRatio(void);

	//baryon to photon ratio
	double BaryonToRadiationRatio(void);

	//Redshift of matter-radiation equality
	double z_eq;
	double z_eq_MatterToRadiation(void);

	//Scale of particle horizon at matter-radiation equality
	double k_eq;
	double k_eq_MatterToRadiation(void);

	//Epoch when baryons are released from Compton drag
	double z_d;
	double z_d_ComptonDrag(void);

	//Ratio of Baryon to Photon Momentum Density
	double R_BaryonToPhotonMomentumDensity(double z);

	//Sound horizon at the drag epoch in Mpc
	double s_drag;
	double s_DragEpoch(void);

	//Scale of Silk damping
	double k_silk;
	double k_SilkDamping(void);

	//transfer function fitting parameters
	//cdm transfer function
	double alpha_c_EH;
	double alpha_c_EH_compute(void);
	double beta_c_EH;
	double beta_c_EH_compute(void);

	//baryon transfer function
	double alpha_b_EH;
	double alpha_b_EH_compute(void);
	double alpha_b_EH_G(double y);

	double beta_b_EH;
	double beta_b_EH_compute(void);
	double beta_node_EH;
	double beta_node_EH_compute(void);
	double s_tilde_EH(double k);

	double T_EH(double k);
	double T_0_tilde_EH(double k, double alpha, double beta);
	double T_c_EH(double k);
	double T_b_EH(double k);

	/*! \fn double P_EH(double k, double z)
	*  \brief EH98 cdm+baryons power spectrum */
	double P_EH(double k, double z);

	/*! \fn double Delta2_EH(double k, double z)
	*  \brief Unit-free CDM+baryon power spectum */
	double Delta2_EH(double k, double z);

	/*! \fn double sigma_EH(double M)
	*  \brief RMS density fluctuations on a smoothing
	*  mass scale M for the EH98 power spectrum*/
	double sigma_EH(double M);

	/*! \fn double sigma(double M)
	*  \brief RMS density fluctuations on a smoothing
	*  mass scale M */
	double sigma(double M);

double sigma_R(double R);
     

	//Universal Age and Conformal Time
	double UniversalAge(double a);
	double ConformalTime(double z);

	//flags and normalization for sigma
	int initialize_sigma_flag;
	int initialize_sigma_EH_flag;
	int tabulated_sigma_EH_flag;
	double sigma_EH_renormalization;

	/*! \fn void InitializeSigma(void)
	*  \brief Initialize the sigma lookup table memory */
	void InitializeSigma(void);

	/*! \fn void SetSigmaMemoryToNull(void)
	*  \brief Set pointers associated with sigma
	*  lookup tables to NULL. */
	void SetSigmaMemoryToNull(void);

	/*! \fn void FreeSigmaMemory(void)
	*  \brief Free memory associated with sigma
	*  lookup tables. */
	void FreeSigmaMemory(void);

	/*! \fn void InitializeSigma_EH(void)
	*  \brief Initialize the sigma for EH98 power spectrum */
	void InitializeSigma_EH(void);

	/*! \fn void TabulateSigma_EH(void)
	*  \brief Tabulate sigma for the EH98 power spectrum */
	void TabulateSigma_EH(void);

	/*! \fn void Print_sigma_M_EH(void)
	*  \brief Print sigma as a function of M */
	void Print_sigma_M_EH(void);

	/*! \fn void Print_Delta2_k_EH(void)
	*  \brief  Print Delta2 as a function of k */
	void Print_Delta2_k_EH(void);

	/*! \fn double dndm_press_schechter(double M, double z)
	 *  \brief Press-Schechter Halo Mass Function at redshift z */
	//double dndm_press_schechter(double M, double z);

	/*! \fn double dndm_sheth_tormen(double M, double z)
	 *  \brief Sheth-Tormen Halo Mass Function at redshift z */
	//double dndm_sheth_tormen(double M, double z);

	/*! \fn double b_MW(double M, double z)
	 *  \brief Mo & White 1996 Bias Function */
	//double b_mo_white(double M, double z);

	/*! \fn double b_tinker(double M, double z)
	 *  \brief Tinker et al. 2010 Bias Function */
	//double b_tinker(double M, double z);

	/*! \fn double nu(double M, double z)
	 *  \brief Peak rarity*/
	//double nu(double M, double z);

      

double M_EH_sigma_Tabulated(double sig);


 ///Tabs

double Tk_tab(double k);



/*! \fn double sigma_k_EH_Tabulated(double k) 
 *  \brief sigma(k) tabulated for EH08 power spectrum */
double sigma_k_EH_Tabulated(double k);





/*! \fn double P_k_EH_Tabulated(double k)
 *  \brief P(k) tabulated for EH08 power spectrum */
double P_k_EH_Tabulated(double k);
/*! \fn double D_growth_Tabulated(double z)
 *  \brief Tabulated linear growth factor */
double D_growth_Tabulated(double z);

void SetUpTools(void);

};

/*! \fn	double sigma_M_EH_Tabulated(double M)
 *  \brief sigma as a function of M for EH98 power spectrum */
double sigma_M_EH_Tabulated(double M);
/*! \fn double Delta2_k_EH_Tabulated(double k)
 *  \brief Dimensionless smoothed power spectrum for EH98 */
double Delta2_k_EH_Tabulated(double k);
/*! \fn double sigma2_R_EH_integrand(double lnk, void *params)
 *  \brief k-space integrand for RMS density fluctuations on a size scale R */
double sigma2_R_EH_integrand(double lnk, void *params);
/*! \fn double ComovingDistanceIntegrand(double ln_z, void *params)
 *  \brief Integrand for the comoving distance calculation */
double ComovingDistanceIntegrand(double ln_z, void *params);

double ConformalTimeIntegrand(double lnz, void *params);


/*! \fn double D_growth_integrand(double lnz, void *params)
 *  \brief Integrand for the linear growth factor */
double D_growth_integrand(double lnz, void *params);



double UniversalAgeIntegrand(double a, void *params);





//Top Hat Window
double W2_RTH(double x);
double W_RTH(double x);


double Tk_tab_in(double k);
/*! \fn double sigma_R_EH_Tabulated(double R) 
 *  \brief RMS density fluctuations on a size scale R */
double sigma_R_EH_Tabulated(double R);




#endif //_COSMOLOGY
