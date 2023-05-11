
#ifndef _ToBIAS
#define _ToBIAS
//#include<gsl/gsl_spline.h>
//#include"constants.h"


class ToBias
{
    public:

    //   Bias(void);

	/*! \fn double dndm_press_schechter(double M, double z)
	 *  \brief Press-Schechter Halo Mass Function at redshift z */
	double n_m_press_schechter(double M, double z);

	/*! \fn double dndm_sheth_tormen(double M, double z)
	 *  \brief Sheth-Tormen Halo Mass Function at redshift z */
	double n_m_sheth_tormen(double M, double z);

         /*! \fn double dndm_LV(double M, double z)
	 *  \brief LV et al non-Gaussian Halo Mass Function at redshift z */
	double n_m_LV_et_al(double fnl,double M, double z);


////Mass derivatives of mass functions



    /*! \fn double dndm_press_schechter(double M, double z)
	 *  \brief Press-Schechter Halo Mass Function at redshift z */
	double dndm_press_schechter(double M, double z);

	/*! \fn double dndm_sheth_tormen(double M, double z)
	 *  \brief Sheth-Tormen Halo Mass Function at redshift z */
	double dndm_sheth_tormen(double M, double z);

         /*! \fn double dndm_LV(double M, double z)
	 *  \brief LV et al non-Gaussian Halo Mass Function at redshift z */

       /*! \fn double dndm_LV(double M, double z)
	 *  \brief LV et al non-Gaussian Halo Mass Function at redshift z */
	double dndm_LV(double fnl,double M, double z);
        
//////Press-Schechter Biass function 

	/*! \fn double b_PS(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
	double E_b_PS_10(double M, double z);
     
 
       /*! \fn double b_PS(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
	double E_b_PS_01(double fnl, double M, double z);

  


       /*! \fn double b_PS(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
	double E_b_PS_02(double fnl,double M, double z);
     
 
       /*! \fn double b_PS(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
	double E_b_PS_20( double M, double z);

         
        /*! \fn double b_PS(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
	double E_b_PS_11(double fnl,double M, double z);



/// Sheth and Torman Eulerican bias function

	/*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
	double E_b_ST_10(double M, double z);


       /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Bias Function */
	double E_b_ST_01(double fnl,double M, double z);
   
 

        /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
	double E_b_ST_20(double M, double z);


       /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Bias Function */
	double E_b_ST_02(double fnl,double M, double z);

       
         /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Bias Function */
	double E_b_ST_11(double fnl,double M, double z);
	

double E_b_ST_30(double M, double z);


///Lagrangian frame///
double b30_STL(double M, double z);
double b20_STL(double M, double z);
double b10_STL(double M, double z);



   ///// LV et al Bias function 

     /// Sheth and Torman  Halo Eulerican bias function

	/*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
	double E_b_LV_10(double fnl,double M, double z);


       /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Bias Function */
	double E_b_LV_01(double fnl,double M, double z);
   
 

        /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
	double E_b_LV_20(double fnl,double M, double z);


       /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Bias Function */
	double E_b_LV_02(double fnl,double M, double z);

       
         /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Bias Function */
	double E_b_LV_11(double fnl,double M, double z);

        /*! \fn double nu(double M, double z)
	 *  \brief Peak rarity*/
	double nu(double M, double z);

       /*! \fn double b_eff_mo_white(double z)*/
      // double b_eff_PS(double z);

         /*! \fn double b_eff_tinker(double z)*/
      // double b_eff_ST(double z);
      double Mofnu(double nu, double z);
      void Initialize_nu_tab( void);

 //double  OmegaHIST(double z);
//double OmegaHILV(double fnl, double z);
// Skewness paramter
      double C_3(double fnl, double M, double z);
 
   // double E_b_PS_20(double M, double z);
~ToBias();

 // private:


        

};

  ///Bias
 double b_PS_20(double M, double z);
 //Mass from nu

double Mofnu_int(double nu, double z);
/*! Effective bias intgrand for PS*/
//
   double domin_b_eff_PS_intgrnd(double z, double M);
 double Num_b_eff_PS_intgrnd(double z, double M);

/*! Effective bias intgrand for ST*/
double domin_b_eff_ST_intgrnd(double z, double M);
double Num_b_eff_ST_intgrnd(double z, double M);

double dlnnudlnM(double M,double z);
double b_eff_ST_result(double z);

//linear spherical collapse overdensity
	double b_LV_10(double fnl,double M, double z);


       /*! \fn double dsigmadM(double M)
 *  \brief derivative of sigma(M) wrt M */
   double dsigmadM(double M);

    /*Non-Gaussian skwness parameter*/
     double  sigma_S_3(double fnl, double M,double z);

    /*derivative Non-Gaussian skwness parameter*/
    double  dsigma_S_3dnu(double fnl, double M,double z);

    //Non-gaussian extension of PS by LV
    double Rnu(double fnl, double M, double z);
  
   //Derivative of Non-gaussian extension of PS
    double dRnudnu(double fnl,double M, double z);
    //Second Derivative of Non-gaussian extension of PS
    double ddRnuddnu_R(double fnl, double M, double z);

/*! \fn double nu_Mz(double M, void *params)
 *  \brief Peak rarity nu = delta_c / [D(z)*sigma(M)] */
     double nu_Mz(double M, double z);

/*! \fn double dnudM(double M, void *params)
 *  \brief Derivative of nu(M) wrt M */
    double dnudM(double M, double z);


  //ddnuddlnM
    double ddnuddM(double M, double z);

/*! \fn double dsigmadM(double M)
 *  \brief derivative of sigma(M) wrt M */
    double dsigmadM(double M);

/*! \fn double f_nu_PS(double nu, void *params)
 *  \brief First-Crossing Distribution for Press-Schechter */
    double nu_f_nu_PS(double nu);

/*! \fn double dndm_PS(double M, void *params)
 *  \brief Number of Collapsed objects for Press-Schecher*/
    double n_m_PS(double M,double z);

/*! \fn double f_nu_ST(double nu, void *params)
 *  \brief First-Crossing Distribution for Sheth-Tormen*/
    double nu_f_nu_ST(double nu);

/*! \fn double dndm_ST(double M, void *params)
 *  \brief Number of Collapsed objects for Sheth-Tormen*/
double n_m_ST(double M, double z);

  double n_m_ST_int(double M, double z);

/*! \fn double dndm_ST(double M, void *params)
 *  \brief Number of Collapsed objects for LV et al*/
   //double n_m_LV(double fnl, double M, double z);


/*! \fn double b_MW(double M, void *params)
 *  \brief Halo Bias Function for Press-Schecher*/
double b_PS(double M, double z);


  /*! \fn double b_ST(double M, void *params)
 *  \brief Halo Bias Function for ST function*/
double b_ST_10(double M,double z);

/*! \fn double b_T(double M, void *params)
 *  \brief Halo Bias Function for ST function*/
double b_ST_20(double M, double z);

 double nu_Mz_z(double M, double z);
    
	/*! \fn double dndm_press_schechter(double M, double z)
	 *  \brief Press-Schechter Halo Mass Function at redshift z */
	double dndm_press_schechter_in(double M, double z);

	/*! \fn double dndm_sheth_tormen(double M, double z)
	 *  \brief Sheth-Tormen Halo Mass Function at redshift z */
	double dndm_sheth_tormen_in(double M, double z);

         /*! \fn double dndm_LV(double M, double z)
	 *  \brief LV et al non-Gaussian Halo Mass Function at redshift z */
	double dndm_LV_in(double fnl,double M, double z);


double n_m_press_schechter_int(double M, double z);

double n_m_ST_int(double M, double z);
double n_m_LV_int(double fnl, double M, double z);


double C_3int(double fnl, double nu,double z);


//Top Hat Window
//double W2_RTH(double x);
//double W_RTH(double x);




#endif //ToBIAS
