#ifndef _HIBIAS
#define _HIBIAS


class HIBias
{
    public:
   //// calculates Omega_HI given a Halo mass function


         double T_bST(double z);
        
         double Tb_Sp(double z);
         double be_Sp(double z);
         double OmegaHIST( double z);
         double ShotNoise( double z);


     //* HI mass in terms of Halo mass//
        double HImass( double z, double M);

       //////Press-Schechter Biass function

	/*! \fn double b_PS(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
	    double HI_b_PS_10( double z);
     
 
       /*! \fn double b_PS(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
        double HI_b_PS_01(double fnl, double z);

  


       /*! \fn double b_PS(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
	    double HI_b_PS_02(double fnl, double z);
     
 
       /*! \fn double b_PS(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
        double HI_b_PS_20(  double z);

         
        /*! \fn double b_PS(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
	    double HI_b_PS_11(double fnl, double z);



/// Sheth and Torman Eulerican bias function

	/*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
        double HI_b_ST_10(double z);


       /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Bias Function */
	    double HI_b_ST_01(double fnl, double z);
   
 

        /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
	   double HI_b_ST_20( double z);


       /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Bias Function */
	   double HI_b_ST_02(double fnl, double z);

       
         /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Bias Function */
	   double HI_b_ST_11(double fnl, double z);
	
      double HI_b_ST_30(double z);



   ///// LV et al Bias function 

     /// Sheth and Torman  Halo Eulerican bias function

	/*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
	   double HI_b_LV_10(double fnl, double z);


       /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Bias Function */
	   double HI_b_LV_01(double fnl, double z);
   
 

        /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Halo Bias Function */
	   double HI_b_LV_20(double fnl,double z);


       /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Bias Function */
	   double HI_b_LV_02(double fnl, double z);

       
         /*! \fn double b_ST(double M, double z)
	 *  \brief Eulerian Bias Function */
	    double HI_b_LV_11(double fnl, double z);
    
    void InitializeAverageloop(void);

       ~HIBias();




};

double CH_interim(double z);
double bevo_int(double z);
double OmegaHIST_int(double z);
    double HI_b_ST_20_in(double z);
   double HI_b_ST_10_in(double z);
      void InitializeST_10_Average(double z);
    double HI_b_ST_10_int(double z);
     void InitializeST_20_Average(double z);
    double HI_b_ST_20_int(double z);
  double HImass_int( double z, double M);

double HI_b_ST_30_int(double z);
     double ShotNoise_intgrnd(double z, double M);
       double Num_b_ST_intgrnd_10(double z, double M);


       double Num_b_ST_intgrnd_20(double z, double M);
double Num_b_ST_intgrnd_30(double z, double M);

      double Deno_b_ST_intgrnd_10(double z, double M);

      double Deno_b_ST_intgrnd_20(double z, double M);

double  T_bST_int(double z);
double ShotNoise_intgrnd(double z, double M);
double n_HILV( double z);
double n_HILVpr(double z);
double b_eHIST(double z);

#endif //HIBIAS
