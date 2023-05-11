




#include"MyCosmology.h"
#include"ToB_Bias.h"
#include"MyHIBias.h"

using namespace std;

/*
HIBias::HIBias()
{
    //SetPowerSpectrum();
}*/
HIBias::~HIBias()
 {
    //SetPowerSpectrum();
 }


const double Delta_HI =pow(10, 12);
const double HIM_i = pow(10,14);
double HIM_min = pow(10,8);
double HIM_max = pow(10,13);
int n_Mass = 100;
const double zminH = 1e-3;
const double zmaxH = 6;
const double EPSREL = 1e-2;
//******************************************************
//Some constant parmeters
//********************************************************
// The mass here are Halo mass, it has been changed to improve accuracy 
  
   



///Speed up HI_b_ST_10_int
    int z_bins =30;
    std::vector <double> z_B_ST_10(z_bins);
     std::vector <double> HI_b_ST_10_tab(z_bins);
    std::vector <double> z_B_ST_20(z_bins);
     std::vector <double> HI_b_ST_20_tab(z_bins);
    std::vector <double> Tb_tab(z_bins);
   std::vector <double> HI_b_ST_30_tab(z_bins);
  std::vector <double> be_tab(z_bins);
   std::vector <double> SHotNoise_tab(z_bins);


///********************************************************************************
//Engine room for the ST model
//******************************************************************************************
     std::vector <double> Mass_tab(n_Mass);
      //std::vector <double> Mass_tab2(n_Mass);
     std::vector <double> Numerator_ST_10(n_Mass);
     std::vector <double> Denominator_ST_10(n_Mass);

   std::vector <double> Numerator_ST_20(n_Mass);
     std::vector <double> Denominator_ST_20(n_Mass);


//********************************************************************
//Initailize
///************************************************

/*
void HIBias::InitializeAverageloop(void)
{
//InitializeST_10_Average();
//InitializeAverageloop();
}*/


//Following Bagla et al 2010
double HIBias::HImass(double z, double M)
{
return HImass_int(z, M);
}

double HImass_int( double z, double M)
{
  ////Define the fitting paramters

       double A;double d;
       double ctwo; 
      double cone = pow(10,11);
       double b = 2.65;
       double XHI_gal = 0.15;

       Cosmology C;
      double OmegaM = C.Omega_m;
      double Omegab = C.Omega_b;
      double tmp1 = pow((1.0+M/cone),b);
      double tmp2 = pow((1.0+M/ctwo),d);

    double another_norm = 100;

    return 100.0*pow(M,0.6);//Mass_HI;


}

double Hmass_min(double z)
{
    Cosmology C;
//double tmp = pow(10,9)/pow((1.0+z),1.5);
    double tmp =  1e8;
    return tmp;
//return C.Mmin;
}

double Hmass_max(double z)
{
    Cosmology C;
//double tmp = 8.0*pow(10,14)/(27.0*pow((1.0+z),1.5));
//double tmp = 8.0*C.Mmax/(27.0*pow((1.0+z),1.5));
     double tmp  = 1e14;
     return tmp;
//return C.Mmax;
}
double CH_interim(double z)
{
    Cosmology C;
  // double c_light =  C.H_0*3000.0/C.h;
     const double c_light =  C.H_0*3000.0;
    const  double c_light_km_s =  299792.0;
    return C.H(z)/(c_light_km_s*(1.0+z));
}






double CurlyHoverH0(double z)
{
Cosmology C;
double a = 1.0/(1.0 + z);
double tmp = a*C.E(z);
return tmp;
}
double  HIBias::T_bST(double z)
{
      Cosmology C;
         //in mill K I removed 1 at the denominator//In K now
     // double tmp =  566.0*C.h*OmegaHIST_int(z)*pow((1.0+z),2)/(0.003* CurlyHoverH0(z)*pow(10,6));
      //return tmp;
    return T_bST_int(z);
}
double  T_bST_int(double z)
{
      Cosmology C;
         //in mill K I removed 1 at the denominator//In K now
      double tmp =  566.0*C.h*OmegaHIST_int(z)*pow((1.0+z),2)/(0.003* CurlyHoverH0(z)*pow(10,6));
      return tmp;
}




  double HIBias::OmegaHIST( double z)
{
return OmegaHIST_int(z);
}



double OmegaHIST_int(double z)
{
 Cosmology C;
  //cout<<"calling ST mass function"<<"\t"<<Deno_b_ST_intgrnd_10(z, 10e10)<<"\t"<<"rho_crit="<<C.rho_crit<<endl;
    const double EPSREL = 1e-2;
      const double EPSABS = 1e-2;
   
    double rhoHI = Integrate(bind(Deno_b_ST_intgrnd_10,z, _1),Hmass_min(z),Hmass_max(z),EPSREL,EPSABS);
    double tmp = rhoHI/(pow((1.0+z),3)*C.rho_crit);
      return tmp;
}



////Public functions for the Eulerian HI bias in Shet-Torman model


double HIBias::HI_b_ST_10( double z)
{
return HI_b_ST_10_int(z);
}

double HIBias::HI_b_ST_20( double z)
{
return HI_b_ST_20_int(z);
}

double HIBias::HI_b_ST_01(double fnl, double z)
{
   Cosmology C;
  double delta_c = C.delta_c_z(z);// 1.686;
  return 2.0* fnl*delta_c*(HI_b_ST_10_int(z)-1.0);
}

double HIBias::HI_b_ST_02(double fnl, double z)
{

   Cosmology C;
  double delta_c = C.delta_c_z(z);// 1.686;
//double delta_c = 1.686;
// Here we set gnl to zero
return 4.0*pow(fnl,2)*delta_c*(delta_c*HI_b_ST_20_int( z)-2.0*HI_b_ST_10_int(z));
}

double HIBias::HI_b_ST_11(double fnl, double z)
{
   Cosmology C;
  double delta_c = C.delta_c_z(z);// 1.686;
//double delta_c = 1.686;
return 2.0* fnl* (delta_c* HI_b_ST_20_int( z)-HI_b_ST_10_int( z))+2.0* fnl*delta_c*HI_b_ST_10_int(z);
}

double HIBias::HI_b_ST_30(double z)
{
return HI_b_ST_30_int(z);
}



double NoiseNumerator(double z, double M)
   {

     //mass is in h^-1 Msun
      ToBias TB;
     //double dnudm =   M*TB.dndm_sheth_tormen(M, z);
      double dnudm =   pow(HImass_int(z, M),2)*TB.n_m_sheth_tormen(M, z);
  //cout <<"dndm="<<"\t"<<TB.n_m_sheth_tormen(M, z) <<endl;
      return dnudm;
}

/*
double ShotNoiseHIST_int(double z)
{
 Cosmology C;
  //cout<<"calling ST mass function"<<"\t"<<Deno_b_ST_intgrnd_10(z, 10e10)<<"\t"<<"rho_crit="<<C.rho_crit<<endl;
    const double EPSREL = 1e-2;
      const double EPSABS = 1e-2;
   
    double rhoHI = Integrate(bind(Deno_b_ST_intgrnd_10,z, _1),Hmass_min(z),Hmass_max(z),EPSREL,EPSABS);
    double numerator = Integrate(bind(NoiseNumerator,z, _1),Hmass_min(z),Hmass_max(z),EPSREL,EPSABS);
    double tmp = pow(T_bST_int(z)/rhoHI,2)*numerator;
      return tmp;
}
*/
/*
double HIBias::ShotNoisePk(double z)
{
return ShotNoise_int(z);
}
*/
///**********************************************************************************************************
// Public function for HI Bias in LV model
//***********************************************



 
 



   double HI_b_ST_10_in(double z)
{   
    const double EPSREL = 1e-2; 
   const double EPSABS = 1e-2; 
  

     double tmp1 = Integrate(bind(Num_b_ST_intgrnd_10,z,_1),Hmass_min(z),Hmass_max(z),EPSREL,EPSABS);
     double tmp2 = Integrate(bind(Deno_b_ST_intgrnd_10,z,_1),Hmass_min(z),Hmass_max(z),EPSREL,EPSABS);
    
     double tmp = tmp1/tmp2;
    return tmp;
}

double HI_b_ST_20_in(double z)
{ 

     const double EPSREL = 1e-2;
    const double EPSABS = 1e-2;
     double tmp1 = Integrate(bind(Num_b_ST_intgrnd_20,z,_1),Hmass_min(z),Hmass_max(z),EPSREL,EPSABS);
     double tmp2 = Integrate(bind(Deno_b_ST_intgrnd_20,z,_1), Hmass_min(z),Hmass_max(z),EPSREL,EPSABS);

     double tmp = tmp1/tmp2;
    return tmp;
}

double HI_b_ST_30_in(double z)
{ 

     const double EPSREL = 1e-2;
    const double EPSABS = 1e-2;
     double tmp1 = Integrate(bind(Num_b_ST_intgrnd_30,z,_1),Hmass_min(z),Hmass_max(z),EPSREL,EPSABS);
     double tmp2 = Integrate(bind(Deno_b_ST_intgrnd_20,z,_1), Hmass_min(z),Hmass_max(z),EPSREL,EPSABS);

     double tmp = tmp1/tmp2;
    return tmp;
}


double ShotNoise_in(double z)
{
    const double EPSREL = 1e-2;
   const double EPSABS = 1e-2;
    double tmp1 = Integrate(bind(ShotNoise_intgrnd,z,_1),Hmass_min(z),Hmass_max(z),EPSREL,EPSABS);
    double tmp2 = Integrate(bind(Deno_b_ST_intgrnd_20,z,_1), Hmass_min(z),Hmass_max(z),EPSREL,EPSABS);

    double tmp = tmp1/(pow(tmp2,2));
   return tmp;
}



double b_eHIST(double z)
{
double tmp =  CH_interim(z)*(3.0 - (1.0+z)*n_HILVpr( z)/n_HILV(z));
return tmp;
//return n_HILVpr(fnl, z);
}



double n_HILVpr(double z)
{

 double delta = 1e-2;
    double z1 = z - delta ;
    double z2 = z + delta ;
    double nu1 = n_HILV(z1);
    double nu2 = n_HILV( z2);
     double tmp =  (nu2-nu1)/(z2-z1);
     return tmp;
}
double n_HILV( double z)
{
   Constants K;
   Cosmology C;
   double c_light =  C.H_0*3000.0;
   //const  double c_light_km_s =  299792.0/C.h;
   double nu_21 = 1420.0*1e6;// Hz
   double E_21 = K.h_planck*nu_21; //J
   double CurlyA = (3.0 *pow(K.h_planck*K.c,3)* K.A_10)/(32.0*K.pi*pow(E_21,2));
   double tmp  = pow((1.0 + z),3)* CH_interim(z)*T_bST_int(z)/CurlyA;
   return tmp;
}

void HIBias::InitializeAverageloop(void)
{
//cout<<"\t \t Pre-Computing HI bias\n\t .........////////....../////......./////\n \t\t Please wait,\n \t there are averaging integrals invloved....."<<endl;
//InitializeST_10_Average();
    int localthread_number =4;
 #pragma omp parallel private(localthread_number)
{
 #pragma omp for schedule(static) nowait
 for(int i=0;i<z_bins;i++)
{
//z_B_ST_10[i] = zminH + i*(zmaxH-zminH)/z_bins;
   z_B_ST_10[i] = zminH + i*(zmaxH-zminH)/z_bins;
//z_B_ST_10[i] = exp( (  (log(zmaxH)-log(zminH))*((double) i)/((double) z_bins-1) + log(zminH) ));
    HI_b_ST_10_tab[i] = HI_b_ST_10_in(z_B_ST_10[i]);
    HI_b_ST_20_tab[i]  = HI_b_ST_20_in(z_B_ST_10[i]);
    HI_b_ST_30_tab[i]  = HI_b_ST_30_in(z_B_ST_10[i]);
    SHotNoise_tab[i] = ShotNoise_in(z_B_ST_10[i]);
    Tb_tab[i] = T_bST(z_B_ST_10[i]);
    be_tab[i] =b_eHIST(z_B_ST_10[i]);

cout<<"i="<<i<<"\t z="<<z_B_ST_10[i]<<"\t b_10 ="<<HI_b_ST_10_tab[i]<<"\tb_20 ="<<HI_b_ST_20_tab[i]<<"\tb_30 ="<<HI_b_ST_30_tab[i]<<"\tTb="<<Tb_tab[i]<<"\t b_e="<<be_tab[i]<<"\t shot-noise="<<SHotNoise_tab[i]<<"\t be/H="<<be_tab[i]/CH_interim(z_B_ST_10[i])<<endl;

}
}
}

   double HI_b_ST_10_int(double z)
{    
    Cosmology C;
   Spline<double, double> CubicSpline_ST_10(z_B_ST_10,HI_b_ST_10_tab);
     double tmp1 = 0.0;
   tmp1 = CubicSpline_ST_10.interpolate(z);
   
   // tmp1 = C.b10;
//cout<<"\t b10 in use is from MyCosmology.cpp and it is equal to "<<C.b10<<endl;
   return tmp1;
 // return 0.85*sqrt(1+z);
}

   double HI_b_ST_20_int(double z)
{    
    Cosmology C;
  
    Spline<double, double> CubicSpline_ST_20(z_B_ST_10,HI_b_ST_20_tab);
      double tmp2 =0.0;
     tmp2 = CubicSpline_ST_20.interpolate(z);
    // tmp2 = C.b20;
//cout<<"\t b20 in use is from MyCosmology.cpp and it is equal to "<<tmp2<<endl;
    return tmp2;
   //return -0.29*sqrt(1+z);
}
    
   double HI_b_ST_30_int(double z)
{    
    Cosmology C;
  
    Spline<double, double> CubicSpline_ST_30(z_B_ST_10,HI_b_ST_30_tab);
      double tmp2 =0.0;
     tmp2 = CubicSpline_ST_30.interpolate(z);
    // tmp2 = C.b20;
//cout<<"\t b20 in use is from MyCosmology.cpp and it is equal to "<<tmp2<<endl;
    return tmp2;
    // return 0.85*sqrt(1+z) - 0.3*(1+z); 
}



double HIBias::Tb_Sp(double z)
{
       Spline<double, double> CubicSpline_Tb(z_B_ST_10,Tb_tab);
       double tmp =  CubicSpline_Tb.interpolate(z);
 // cout <<"t_b="<<tmp<<endl;
     return tmp;
   	
}
double HIBias::be_Sp(double z)
{
       Spline<double, double> CubicSpline_be(z_B_ST_10,be_tab);
       double tmp =  CubicSpline_be.interpolate(z);
 // cout <<"t_b="<<tmp<<endl;
     return tmp;
   	
}
                       
double HIBias::ShotNoise(double z)
{
          Spline<double, double> CubicSpline_SN(z_B_ST_10,SHotNoise_tab);
        double tmp =  CubicSpline_SN.interpolate(z);
                        // cout <<"t_b="<<tmp<<endl;
        return tmp;
                              
}





///***************************************Integrand tools***************************************************************







///ST integrand


 double Num_b_ST_intgrnd_10(double z, double M)
   {

     //mass is in h^-1 Msun
     ToBias TB;

     //double dnudm =   TB. E_b_ST_10(M,z)*M*TB.dndm_sheth_tormen(M, z);
     double dnudm =   TB.E_b_ST_10(M,z)*HImass_int(z, M)*TB.n_m_sheth_tormen(M, z);
  
      return dnudm;
}

   double Num_b_ST_intgrnd_20(double z, double M)
   {
      //mass is in h^-1 Msun
      ToBias TB;
     //double dnudm = TB. E_b_ST_20(M,z)*M*TB.dndm_sheth_tormen(M, z);
      double dnudm = TB.E_b_ST_20(M,z)*HImass_int(z, M)*TB.n_m_sheth_tormen(M, z);
    
 return dnudm;
}

  double Num_b_ST_intgrnd_30(double z, double M)
   {
      //mass is in h^-1 Msun
      ToBias TB;
     //double dnudm = TB. E_b_ST_20(M,z)*M*TB.dndm_sheth_tormen(M, z);
      double dnudm = TB.E_b_ST_30(M,z)*HImass_int(z, M)*TB.n_m_sheth_tormen(M, z);
     return dnudm;
}



double Deno_b_ST_intgrnd_10(double z, double M)
   {

     //mass is in h^-1 Msun
      ToBias TB;
     //double dnudm =   M*TB.dndm_sheth_tormen(M, z);
      double dnudm =   HImass_int(z, M)*TB.n_m_sheth_tormen(M, z);
  //cout <<"dndm="<<"\t"<<TB.n_m_sheth_tormen(M, z) <<endl;
      return dnudm;
}

   double Deno_b_ST_intgrnd_20(double z, double M)
   {
      //mass is in h^-1 Msun  
      ToBias TB;
     //double dnudm = M*TB.dndm_sheth_tormen(M, z);
     double dnudm = HImass_int(z, M)*TB.n_m_sheth_tormen(M, z);
   // cout <<"\tHImass="<<"\t"<<HImass_int(z, M)<<endl;
 return dnudm;
}

double ShotNoise_intgrnd(double z, double M)
{
   //mass is in h^-1 Msun
   ToBias TB;
  //double dnudm = M*TB.dndm_sheth_tormen(M, z);
  double dnudm = HImass_int(z, M)*HImass_int(z, M)*TB.n_m_sheth_tormen(M, z);
 
return dnudm;
}



