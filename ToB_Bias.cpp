#include"ToB_Bias.h"

#include"MyCosmology.h"

using namespace std;


/*
ToBias::ToBias(void)
{
    //SetPowerSpectrum();
}
*/
ToBias::~ToBias(void)
{
    //SetPowerSpectrum();
}

const double Mmin = 5;
 const double Mmax = 16;
const double Delta_m = 3;
  const double M_i = 14;
   double Av_M_min = 12;
   double Av_M_max = 15;

  

//const double Mmin = 10e12;
// const double Mmax = 10e16;


//const double M_i = 10e14;
 //  double M_min = M_i- 10e5/2;
 //  double M_max = M_i+ 10e5/2;


/////////////////////////////////////////
// Member public functions for halo mass functions
/////////////////////////////////////////
/*! \fn double nu(double M, double z)
*  \brief Peak rarity*/
double ToBias::nu(double M, double z)
{
   Cosmology C;
  //Tabulate T;
  //return peak rarity
//cout <<"Sigma(M)="<<T.sigma_EH(5)<<endl;
  return nu_Mz(M, z);
  // return C.D_growth_Tabulated(z);
}


double nu_Mz_z(double M, double z)
{
  //Tabulate T;
  //return peak rarity
//cout <<"Sigma(M)="<<T.sigma_EH(5)<<endl;
  return nu_Mz(M,z);
   //return M*z;
}

int nn_Mass =500;
 std::vector <double> MMass_tab(nn_Mass);
 std::vector <double> z_tab(nn_Mass);
std::vector <double> nu_tab(nn_Mass);

/*
void ToBias::Initialize_nu_tab( void)
{ 
   double zmin1 = 0.000001;
   for(int i =0; i < nn_Mass;i++)
   {
   MMass_tab[i] =pow(10, Mmin + i*(Mmax-Mmin)/nn_Mass);
  // MMass_tab[i] =Mmin + i*(Mmax-Mmin)/nn_Mass;
   nu_tab[i] = nu_Mz_z(MMass_tab[i],zmin1);
   }
}*/

/*
double ToBias::Mofnu(double nu, double z)
{
        return Mofnu_int(nu, z);
}


double Mofnu_int(double nu, double z)
{

        Cosmology C;
	double M_nu_z;        
        double Dz      = C.D_growth_Tabulated(z);
        double sig = C.delta_c_z(z)/(Dz*sqrt(nu));
        M_nu_z =  C.M_EH_sigma_Tabulated( sig);

	return M_nu_z;
}
*/

/*
double nu_func(int n, int m,int o)
{
	return nu_Mz_z(n,m,o);
}*/
///////////////////////////////////////
//peak rarity below here
///////////////////////////////////////

/*! \fn double nu_Mz(double M, void *params)
  //equation 3.1 (Toby)
 *  \brief Peak rarity nu = delta_c / [D(z)*sigma(M)] */


double nu_Mz(double M, double z)
{

   //Tabulate T;
   Cosmology C;
  //double *fp     = (double *) params;
  //double z       = fp[0];
  double delta_c = C.delta_c_z(z);// 1.686;
  double Dz      = C.D_growth_Tabulated(z);
  double sig     = C.sigma_EH(M);

  //printf("M %e z %e Dz %e sig %e\n",M,z,Dz,sig);

  //nu = pow(delta_c / sigma(M,z),2)
  double tmp =  delta_c/(Dz*sig);
 //cout<<"nu="<<tmp*tmp<<endl;
  return tmp*tmp;
  //return sig;
}

double diff1(double(*f)(double) , double x0)
{
 double delta = 1.0e-5;///Need to be small for reasonal accuracy
double x1 = x0 -delta;
double x2 = x0 + delta ;
double y1 = f(x1);
double y2 =f(x2);
return (y2-y1)/(x2-x1);
}
/*

double dnudM(double M, void *params)
{
  //nu = delta_c/[D(z) * sigma(M)]
  //dnudM = -delta_c/[D(z) * sigma(M)^2] x dsigmadM
  //dnudM = -[nu/sigma(M)] x dsigmadM
  return -2.0*nu_Mz(M,params)*dsigmadM(M)/sigma_M_EH_Tabulated(M);
}

*/
/*! \fn double dsigmadM(double M)
 *  \brief derivative of sigma(M) wrt M */
/*
double dsigmadM(double M)
{
	//dumb derivative
      Cosmology C;
	double sigma_A = C.sigma_EH(1.01*M);
	double sigma_B = C.sigma_EH(0.99*M);
	return (sigma_A-sigma_B)/(0.02*M);
}

double ddsigmadMsq(double M)
{
       double sigma_A =dsigmadM(1.01*M);
	double sigma_B = dsigmadM(0.99*M);
	return (sigma_A-sigma_B)/(0.02*M);
}
*/
//////////////////Press-Schechter Halo Mass Function tools.//////////////////////////////////////

/*! \fn double dndm_sheth_tormen(double M, double z)	
*  \brief Press-Schechter number of collapsed objects at redshift z */
double ToBias::n_m_press_schechter(double M, double z)
{
  //return n(m) for Press-Schechter
  return n_m_press_schechter_int(M,z);
}

double n_m_press_schechter_int(double M, double z)
{

   // Cosmology C;
  //mass is in h^-1 Msun
  //double fp[2];
  
  //store redshift
  //fp[0] = z; 

  //store matter density in h^2 Msun Mpc^-3
 // fp[1] = C.rho_m(z);

  //return dndM for Press-Schechter
  return n_m_PS(M,z);
}

///////////////////////////////////////
// Halo mass functions below here
///////////////////////////////////////

/*! \fn double dndm_PS(double M, void *params)
 *  \brief Halo Mass Function for Press-Schecher*/
double n_m_PS(double M, double z)
{
  //M is in h^-1 Msun
  Cosmology C;
  //double *fp    = (double *) params;
  //double z      = fp[0]; //redshift
  //double rho_m  = fp[1]; //mean background density in h^2 Msun / Mpc^3
  double nu     = nu_Mz(M,z);
  double M_X = pow(10,M);
 // double dlnnu_dlnM = M*dnudM(M,params)/nu;

  //number of halos per h^3 Mpc^-3 per dM
  return (C.rho_m(z)/(M_X*M_X))*dlnnudlnM(M,z)*nu_f_nu_PS(nu);
}


/*! \fn double f_nu_PS(double nu, void *params)
 *  \brief First-Crossing Distribution for Press-Schechter */
double nu_f_nu_PS(double nu)
{  
  
  //equation 3.3 of Toby
  // double nu     = nu_Mz(M,params);
  return (sqrt(nu/(2.0*M_PI))*exp(-0.5*nu));
}



//////////////////Sheth-Tormen mass function  tools    /////////     tools/////////////////////////////////////////////////////////
//////////////////////////////////////////////////////

/*! \fn double dndm_sheth_tormen(double M, double z)
*  \brief Sheth-Tormen number of collapsed objects at redshift z */
double ToBias::n_m_sheth_tormen(double M, double z)
{
   
  //return n(M) for Sheth-Tormen
  return n_m_ST_int(M,z);
}



///for internal use only

double n_m_ST_int(double M, double z)
{
    //Cosmology C;
  //mass is in h^-1 Msun
  //double fp[2];
 
  //store redshift
  //fp[0] = z; 

  //store matter density in h^2 Msun Mpc^-3
  //fp[1] = C.rho_m(z);

  //return dndM for Sheth-Tormen
  return n_m_ST(M,z);
}

/*! \fn double dndm_ST(double M, void *params)
 *  \brief Halo Mass Function for Sheth-Tormen*/
double n_m_ST(double M, double z)
{
  //M is in h^-1 Msun
  Cosmology C;
  //double *fp    = (double *) params;
 // double z      = fp[0]; //redshift
 // double rho_m  = fp[1]; //mean background density in h^2 Msun / Mpc^3
  double nu     = nu_Mz(M,z);
  //double dlnnu_dlnM = M*dnudM(M,params)/nu;

  //printf("M %e z %e rho_m %e nu %e dnu_dM %e\n",M,z,rho_m,nu,dnu_dM);

  //number of halos per h^3 Mpc^-3 per dM
 double tmp =  (C.rho_m(z)/(M*M))*nu_f_nu_ST(nu)*dlnnudlnM(M,z);//;
//cout <<"n_mST="<<C.rho_m(z)<<endl;
   return tmp;
}

/*! \fn double f_nu_ST(double nu, void *params)
 *  \brief First-Crossing Distribution for Sheth-Tormen*/
double nu_f_nu_ST(double nu)
{
  //equation 3.4 of Toby
  double p  = 0.3;
  double q   = 0.707;
  double A   = 0.322; 
  double nuq = q*nu;
  double tmp = A*(1.0 + pow(nuq,-p))*sqrt(nuq/(2.0*M_PI))*exp(-0.5*nuq);
// cout <<"nufST="<<tmp<<endl;
  return tmp;
}

 //dlnnudlnM
double dlnnudlnM(double M,double z)
{
    //Tabulate T;
    // nu();
    double delta = 1.0e-4;
    double M1 = M - delta ;
    double M2 = M + delta ;
    double coef = M/nu_Mz_z(M,z);
    double lnM1 = log(M) - delta ;
    double lnM2 = log(M) + delta ;
    double nu1 = log(nu_Mz_z(exp(lnM1),z));
    double nu2 = log(nu_Mz_z(exp(lnM2),z));
    double deri =  (nu2-nu1)/(lnM2-lnM1);
    double tmp  = deri;
  //cout <<"derivative="<<"\t"<<tmp<<endl;
   return tmp;
  //nu = delta_c/[D(z) * sigma(M)]
  //dnudM = -delta_c/[D(z) * sigma(M)^2] x dsigmadM
  //dnudM = -[nu/sigma(M)] x dsigmadM
  //return -nu_Mz(M,params)*dsigmadM(M)/T.sigma_M_EH_Tabulated(M);
}
double dnudM(double M, double z)
{
    //Tabulate T;
    // nu();
    double delta = 1e-5;
    double M1 = log(M) - delta ;
    double M2 = log(M) + delta ;
    double nu1 = nu_Mz(exp(M1),z);
    double nu2 = nu_Mz(exp(M2),z);
     double tmp =  (nu2-nu1)/(exp(M2)-exp(M1));
     return tmp;
  //nu = delta_c/[D(z) * sigma(M)]
  //dnudM = -delta_c/[D(z) * sigma(M)^2] x dsigmadM
  //dnudM = -[nu/sigma(M)] x dsigmadM
  //return -nu_Mz(M,params)*dsigmadM(M)/T.sigma_M_EH_Tabulated(M);
}


double dnudM_z(double M,double z)
{
 Cosmology C;
  //mass is in h^-1 Msun
  //double fp[2];
 
  //store redshift
  //fp[0] = z; 

  //store matter density in h^2 Msun Mpc^-3
  //fp[1] = C.rho_m(z);
 double tmp = dnudM( M, z);
//cout <<"Diff="<<"\t"<<tmp<<endl;
return tmp;
}



///////////// LV et al Mass function tools////////////////////////////////////////




//////////////////Sheth-Tormen mass function  tools    /////////     tools/////////////////////////////////////////////////////////
//////////////////////////////////////////////////////

/*! \fn double dndm_sheth_tormen(double M, double z)
*  \brief LV number of collapsed objects at redshift z */
double ToBias::n_m_LV_et_al(double fnl,double M, double z)
{
   
  //return n(M) for LV et al
     return n_m_LV_int(fnl,M,z);
}

/*! \fn double dndm_ST(double M, void *params)
 *  \brief Number of Collapsed objects for LV et al*/

double n_m_LV_int(double fnl, double M, double z)
{
   // double *fp    = (double *) params;
    double n_m_ST_in  = n_m_ST_int(M, z);
    double R_nu  = Rnu(fnl,M,z);
    return n_m_ST_in*R_nu;
}


//Non-gaussian extension of PS
double Rnu(double fnl,double nu, double z)
{
    double sigma_S_3_in = sigma_S_3(fnl,nu,z);
    double dsigma_S_3dnu_in = dsigma_S_3dnu(fnl,nu,z);
    double R_nu =  1.0 + sigma_S_3_in*(pow(nu,3/2)-3*pow(nu,1/2))/(6.0)-dsigma_S_3dnu_in*(pow(nu,3/2)-pow(nu,1/2))/3;
//cout << "R=" <<R_nu<<endl;
     return R_nu;
}

double ToBias::C_3(double fnl, double nu,double z)
{
 
    return C_3int( fnl, nu, z);

}

double C_3int(double fnl, double nu,double z)
{
 //Cosmology C;
  //double REH =  pow((3.0*M/(4.0*pi*rho_m(0))),1.0/3.0);
  //double sigmaR = sigma_R_EH_Tabulated(REH);
    Cosmology C;
    double delta_c = C.delta_c_z(z);
    double Dz    = C.D_growth_Tabulated(z);
    double sig = delta_c/( Dz*sqrt(nu));
   double S3 = fnl*(3.20e-4)*pow(sig,3.141)/sig;//pow(sig,3);
//cout <<"S3=" <<S3<<endl;
   return  S3;

//return fnl*(6.6e-4)*(1.0-log(0.016*M))*D_growth_Tabulated(z)*sigmaR;//sigma_M_EH_Tabulated(M);

}

/*Non-Gaussian skwness parameter*/
     double  sigma_S_3(double fnl,double nu,double z)
 {
    
    return  C_3int(fnl,nu,z);
}



    /*derivative Non-Gaussian skwness parameter*/
    double  dsigma_S_3dnu(double fnl,double nu,double z)
{
    double delta = 1.0e-2;
    double M1 = nu - delta ;
    double M2 = nu + delta ;
    double C_3_1 =sigma_S_3(fnl,M1,z);
    double C_3_2 = sigma_S_3(fnl,M2,z);
    double tmp = (C_3_2-C_3_1)/(M2-M1);
//cout<<"diff="<<tmp<<endl;
   return  tmp;
}




     double ToBias::dndm_press_schechter(double M, double z)
{
   double delta = 1.0e-5;
    double M1 = M - delta ;
    double M2 =M + delta ;
    //double lnM1 = log(M) - delta ;
   // double lnM2 =log(M) + delta ;
     //double A = log(n_m_press_schechter_int(M1,z));
    //double B = log(n_m_press_schechter_int(M2,z));
    //double tmp = n_m_press_schechter_int(M,z)*(B-A)/((lnM2-lnM1)*M);
      double A = n_m_press_schechter_int(M1,z);
    double B = n_m_press_schechter_int(M2,z);
    double tmp = (B-A)/(M2-M1);
   return tmp;
}

	
	double ToBias::dndm_sheth_tormen(double M, double z)
{
   double delta = 1.0e-5;
    double M1 = M - delta ;
    double M2 = M + delta ;
    //double lnM1 = log(M) - delta ;
   // double lnM2 = log(M) + delta ;
  // double C_factor = log(10.0)*M;
   //double A =log(n_m_ST_int(M1,z));
    //double B =log(n_m_ST_int(M2,z));
   // double tmp = n_m_ST_int(M,z)*(B-A)/(M*(lnM2-lnM1));
   double A = n_m_ST_int(M1,z);
    double B = n_m_ST_int(M2,z);
    double tmp = (B-A)/(M2-M1);
   return tmp;
}

	
	double dndm_sheth_tormen_in(double M, double z)
{
   double delta = 1.0e-5;
    double M1 = M - delta ;
    double M2 = M + delta ;
    //double lnM1 = log(M) - delta ;
   // double lnM2 = log(M) + delta ;
  // double C_factor = log(10.0)*M;
   //double A =log(n_m_ST_int(M1,z));
    //double B =log(n_m_ST_int(M2,z));
   // double tmp = n_m_ST_int(M,z)*(B-A)/(M*(lnM2-lnM1));
   double A = n_m_ST_int(M1,z);
    double B = n_m_ST_int(M2,z);
    double tmp = (B-A)/(M2-M1);
   return tmp;
}
         
 double ToBias::dndm_LV(double fnl,double M,double z)
{
  double delta = 1.0e-5;
    double M1 = M - delta ;
    double M2 = M + delta ;
   //  double lnM1 = log(M) - delta ;
   // double lnM2 = log(M) + delta ;
     //double A =log(n_m_LV_int(fnl,M1,z));
    //double B = log(n_m_LV_int(fnl,M2,z));
    //double tmp =n_m_LV_int(fnl,M,z) *(B-A)/(M*(lnM2-lnM1));
   double A =n_m_LV_int(fnl,M1,z);
    double B = n_m_LV_int(fnl,M2,z);
    double tmp =(B-A)/(M2-M1);
   return tmp;
}


////Return Lagrangian bias


/*! \fn double b_MW(double M, void *params)
 *  \brief Halo Bias Function for Press-Schecher*/
double b_PS(double M, double z)
{
  //M is in h^-1 Msun

  double nu     = nu_Mz_z(M,z);
  //double delta_c = 1.686;
   Cosmology C;
  double delta_c = C.delta_c_z(z);// 1.686;

  //halo bias
  return  (nu - 1.0)/delta_c;
}
double ToBias::E_b_PS_10(double M, double z)
{
  return 1+b_PS(M,z);
}

double b_PS_20(double M, double z)
{
    Cosmology C;
  double delta_c = C.delta_c_z(z);// 1.686;
   double nu     = nu_Mz_z(M,z);
   //double delta_c = 1.686;
   double first_term = (nu*nu-2*nu-1.0)/(delta_c*delta_c);
   double second_term = -(nu-1)/(delta_c*delta_c);
   return first_term+second_term;
}



double ToBias::E_b_PS_20(double M, double z)
{
return 8*b_PS(M, z)/21 + b_PS_20(M, z);
}



////Shet Torman
double b_ST_10(double M,double  z)
{
 
   Cosmology C;
  double delta_c = C.delta_c_z(z);// 1.686;
  double nu     = nu_Mz_z(M,z);
  // double delta_c = 1.686;
  double p  = 0.3;
  double q   = 0.707;
  double A   = 0.322;
  double nuq = q*nu;
   double dndnu =-(nuq-1)/(2*nu) - p/(nu*(1+pow(nuq,p)));
   // double first_term = -(nuq-1)/delta_c;
   // double second_term = -(2.0*p)/(delta_c*(1+pow(nuq,p)));

   //Lagrangian Halo bias

  // return first_term+second_term;
  return -(2*nu/delta_c)* dndnu;
}


double b_ST_20(double M,double  z)
{
  //M is in h^-1 Msun
 // double *fp    = (double *) params;
 // double z      = fp[0];
    Cosmology C;
  double delta_c = C.delta_c_z(z);// 1.686;
   double nu     = nu_Mz_z(M,z);
  // double delta_c = 1.686;
   double p  = 0.3;
   double q   = 0.707;
   double A   = 0.322; 
   double nuq = q*nu;
   double ddnddnv =(pow(p,2)+ nuq*p)/(nu*nu*(1.0+pow(nuq,q)))
                  +(pow(nuq,2)-2*nuq-1)/(4*nu*nu);
   double dndnu =-(nuq-1)/(2*nu) - p/(nu*(1+pow(nuq,p)));

   double tmp = (4*nu*nu*(ddnddnv))/(delta_c*delta_c)
                + (2*nu*(dndnu))/(delta_c*delta_c);
  return tmp;
}


////Public functions for the Eulerian Halo bias in Shet-Torman model

double ToBias::E_b_ST_10(double M, double z)
{
double tmp = 1+b_ST_10(M, z);
//cout<<"b_LV10="<<"\t"<<tmp<<endl;
return tmp;
}

double ToBias::E_b_ST_20(double M, double z)
{
return 8.0*b_ST_10(M, z)/21.0 + b_ST_20(M, z);
}

double ToBias::E_b_ST_01(double fnl,double M, double z)
{
   Cosmology C;
  double delta_c = C.delta_c_z(z);// 1.686;
 //double delta_c = 1.686;
return 2* fnl*delta_c*b_ST_10(M, z);
}

double ToBias::E_b_ST_02(double fnl,double M, double z)
{
   Cosmology C;
  double delta_c = C.delta_c_z(z);// 1.686;
//double delta_c = 1.686;
// Here we set gnl to zero
return 4*pow(fnl,2)*delta_c*(delta_c*b_ST_20(M, z)-2*b_ST_10(M, z));
}

double ToBias::E_b_ST_11(double fnl,double M, double z)
{
  Cosmology C;
  double delta_c = C.delta_c_z(z);// 1.686;
//double delta_c = 1.686;
return 2* fnl* (delta_c* b_ST_20(M, z)-b_ST_10(M, z))+2* fnl*delta_c*b_ST_10(M, z);
}
///LV et al











/*
double b_LV_10(double fnl, double M,double  z)
{
 
  double nu     = nu_Mz_z(M,z);
   double delta_c = 1.686;
  double p  = 0.3;
  double q   = 0.707;
  double A   = 0.322;
  double nuq = q*nu;

   //Eulerian Halo bias
  return   (-(2*nu)/delta_c)*(-(nuq-1)/2*nu - p/(nu*(1+pow(nuq,p))+dRnudnu_R(fnl,M,z);
}
*/

//Derivative of Non-gaussian extension of PS
double dRnudnu_R(double fnl, double nu, double z)
{
   return  dRnudnu(fnl, nu, z)/Rnu(fnl,nu,z);
}
//Derivative of Non-gaussian extension of PS
double dRnudnu(double fnl, double nu, double z)
{
    double delta = 1.0e-5;
    double nu1 = nu - delta ;
    double nu2 = nu + delta ;
    double Rnu1 =Rnu(fnl,nu1,z);
    double Rnu2 = Rnu(fnl,nu2,z);
    double tmp = (Rnu2-Rnu1)/((nu2-nu1));
   return  tmp;
}






double b_LV_10(double fnl, double M,double  z)
{
     Cosmology C;
   double delta_c = C.delta_c_z(z);// 1.686;

    double p  = 0.3;
    double q   = 0.707;
    double A   = 0.322; 
   // double delta_c = 1.686;
    double nu     = nu_Mz_z(M,z);
    double nuq = q*nu;
   // double first_term = (nuq-1)/delta_c;
    //double second_term = (2*p)/(delta_c*(1+pow(nuq,p)));
    //double third_term  = -(2*nu/delta_c)*dRnudnu_R(fnl,nu,z);
    //double tmp  = first_term+second_term +third_term;
  double dndnu =-(nuq-1)/(2*nu) - p/(nu*(1+pow(nuq,p)))+dRnudnu_R(fnl,nu,z);
  return -(2*nu/delta_c)* dndnu;
    //return  tmp;
}




double b_LV_20(double fnl, double M,double  z)
{
    Cosmology C;
   double delta_c = C.delta_c_z(z);// 1.686;

    double p  = 0.3;
    double q   = 0.707;
    double A   = 0.322; 
    //double delta_c = 1.686;
    double nu     = nu_Mz_z(M,z);
    double nuq = q*nu;
   // double first_term = (4*(p*p+nuq*p)-(nuq-1)*(1+pow(nuq,p))-2*p)/(delta_c*delta_c*(1+pow(nuq,p)));
    //double second_term = (pow(nuq,2)-2*nuq-1)/(delta_c*delta_c);
    //double Third_term  = (8*nu/(delta_c*delta_c))*((nuq-1)/2+ p/(1+pow(nuq,p)))*dRnudnu_R(fnl,nu,z);
    //double Fourth_term  = ddRnuddnu_R( fnl,nu,z)+(2*nu/(delta_c*delta_c))*dRnudnu(fnl,nu, z);
   //double tmp= first_term +second_term-Third_term+Fourth_term ;

     double ddnddnv_ST =(pow(p,2)+ nuq*p)/(nu*nu*(1.0+pow(nuq,q)))
                  +(pow(nuq,2)-2*nuq-1)/(4*nu*nu);
   double dndnu_ST =-(nuq-1)/(2*nu) - p/(nu*(1+pow(nuq,p)));
    
   double ddnddnv = ddnddnv_ST+ 2*dndnu_ST*dRnudnu_R(fnl,nu,z)+ddRnuddnu_R( fnl,nu,z);
   double dndnu  =  dndnu_ST + dRnudnu_R(fnl,nu,z);

   double tmp = (4*nu*nu*(ddnddnv))/(delta_c*delta_c)
                + (2*nu*(dndnu))/(delta_c*delta_c);
 //double tmp = (4*nu*nu*(ddnddnv_ST))/(delta_c*delta_c)
      //          + (2*nu*(dndnu_ST))/(delta_c*delta_c);
  return tmp;
   
}




double ddnuddM_z(double M, double z)
{
    double delta = 1.0e-5;
    double M1 = M - delta ;
    double M2 = M + delta ;
    double dnudM1 = dnudM_z(M1,z);
    double dnudM2 = dnudM_z(M2,z);
     return (dnudM2-dnudM1)/(M2-M1);
}


//Second Derivative of Non-gaussian extension of PS
double ddRnuddnu_R(double fnl,double nu, double z)
{
    double delta = 1.0e-5;
    double M1 = nu - delta ;
    double M2 = nu + delta ;
    double ddRnu1 =dRnudnu(fnl,M1,z);
    double ddRnu2 = dRnudnu(fnl,M2,z);
    double dRnu1 =Rnu(fnl,M1,z);
    double dRnu2 = Rnu(fnl,M2,z);
    //double dnuDM = dnudM_z(M, z);
   //double tmp1 = (dRnu2-dRnu1)/(M2-M1);
    double tmp2 = (ddRnu2-ddRnu1)/(M2-M1);
   // double tmp3 =  tmp2*pow(dnuDM,-2)+tmp1*pow(ddnuddM_z(M,z),-1);
    return tmp2;///Rnu(fnl,M1,z);
}




double ToBias::E_b_LV_10(double fnl,double M, double z)
{
return 1.0+b_LV_10(fnl,M, z);
//return 1+b_ST_10(M, z);
}

double ToBias::E_b_LV_20(double fnl,double M, double z)
{
return 8.0*b_LV_10(fnl,M, z)/21.0 + b_LV_20(fnl,M, z);
}

double ToBias::E_b_LV_01(double fnl,double M, double z)
{
   Cosmology C;
  double delta_c = C.delta_c_z(z);// 1.686;
  // double Mass = C.M_EH_sigma_Tabulated(delta_c/pow(nu,0.5));
    double nu  =   nu_Mz_z(M,z);
   
  
  double tmp = 2* fnl*delta_c*b_LV_10(fnl,M, z) + (2.0*fnl*C_3(fnl,nu,z)*(pow(nu,3.0/2.0)-3.0*pow(nu,0.5)))/(3.0*Rnu(fnl,nu,z));
return tmp;
}

double ToBias::E_b_LV_02(double fnl,double M, double z)
{
    Cosmology C;
  double delta_c = C.delta_c_z(z);// 1.686;
  // double delta_c = 1.686;
// Here we set gnl to zero
return 4*pow(fnl,2)*delta_c*(delta_c*b_LV_20(fnl,M, z)-2*b_LV_10(fnl,M, z));
}

double ToBias::E_b_LV_11(double fnl,double M, double z)
{
    Cosmology C;
  double delta_c = C.delta_c_z(z);// 1.686;
  //double delta_c = 1.686;
return 2* fnl* (delta_c* b_LV_20(fnl,M, z)-b_LV_10(fnl,M, z))+2* fnl*delta_c*b_LV_10(fnl,M, z);
}



///Computation of Omega_HI




/////////////////////////////////////////////////////Third order bias ////////////////////////////////////////////




double ToBias::b10_STL(double M, double z)
{
return  b_ST_10(M,z);

}
double ToBias::b20_STL(double M, double z)
{
return  b_ST_20(M,z);
}


double ToBias::b30_STL(double M, double z)
{
 Cosmology C;
    double delta_c = C.delta_c_z(z);// 1.686;
    double p  = 0.3;
    double q   = 0.707;
    double A   = 0.322; 
    //double delta_c = 1.686;
    double nu     = nu_Mz_z(M,z);
    double nuq = q*nu;
    double dndnu =-(nuq-1.0)/(2.0*nu) - p/(nu*(1.0+pow(nuq,p)));
   // double ddnddv =(pow(p,2)+ nuq*p)/(nu*nu*(1.0+pow(nuq,q)))+(pow(nuq,2)-2.0*nuq-1.0)/(4.0*nu*nu);
   // double dddndddv = -(8.0*pow(p,3) + 12.0* pow(p,2)*(1.0+ nuq) + p*(6.0*nuq*nuq - 2.0))/(8.0*pow(nu,3)*(1.0+pow(nuq,p)) )+ (3.0 + 3.0*nuq + 3.0*nuq*nuq - pow(nuq,3))/(8.0*pow(nu,3));

double ddnddv = (4*pow(p,2) + 4*nu*p*q + (-1 - 2*nu*q + pow(nu,2)*pow(q,2))*
      (1 + pow(nu*q,p)))/(4.*pow(nu,2)*(1 + pow(nu*q,p)));

 double dddndddv = -(8*pow(p,3) + 12*pow(p,2)*(1 + nu*q) + 
      p*(-2 + 6*pow(nu,2)*pow(q,2)) + 
      (-3 - 3*nu*q - 3*pow(nu,2)*pow(q,2) + pow(nu,3)*pow(q,3))*
       (1 + pow(nu*q,p)))/(8.*pow(nu,3)*(1 + pow(nu*q,p)));

double tmp1 = pow(nu,3)*dddndddv + 3.0* pow(nu,2)*ddnddv + nu*dndnu;
double tmp   = -8.0* tmp1/pow(delta_c,3);
return tmp;


}


double ToBias::E_b_ST_30(double M, double z)
{
const double a1 =1.0;
const double a2 = -17.0/21.0;
const double a3 = 341.0/567.0;

double tmp1 = 6.0*(a2+ a3)*b_ST_10( M, z) + 3.0*(pow(a1,2)+ 2.0*a1*a2)*b_ST_20(M,z) + pow(a1,3)*b30_STL(M, z);

return tmp1;
}
















