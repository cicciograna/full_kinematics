#ifndef FULL_KINEMATICS_HH
#define FULL_KINEMATICS_HH

# include "TMath.h"
# include "TMatrixD.h"
# include "TF1.h"
# include "TRandom3.h"
# include "TRandom.h"
# include <fstream>
# include <iostream>
# include <iomanip>

# include <gmp.h>     // hi-precision library
# include <mpfr.h>
# include "mpreal.h"
// # include <mpf2mpfr.h>

using std::cout;
using std::endl;
using std::ofstream;
using std::fixed;
using std::setprecision;
using mpfr::mpreal;

struct mpreal_matrix{
  mpreal r0c0, r0c1, r0c2;
  mpreal r1c0, r1c1, r1c2;
  mpreal r2c0, r2c1, r2c2;
};

struct mpreal_4vector{
  mpreal x, y, z, e;
};

// Forward declaration
double GenerateNuPz(double neutrino_mass, double xenon_mass, double cesium_mass);
mpreal GenerateNuPzGMP(mpreal mnu, mpreal mxe, mpreal mcs);
double GenerateEnergyFromPM(double particle_pz, double particle_mass);
mpreal GenerateEnergyFromPMGMP(mpreal pz, mpreal mass);

mpreal RandomTheta();
mpreal RandomValidXeStarTheta();
mpreal RandomPhi();
mpreal RandomValidXrayTheta();
mpreal RandomValidXrayPhi();
mpreal RandomThetaWithErr(double error);
mpreal RandomPhiWithErr(double error);
mpreal RandomValidAugerTheta();
mpreal RandomThetaScatter();
mpreal RandomThetaScatter(mpreal velocity);

// ***** RANDOM GENERATOR *****
TRandom3 *r = new TRandom3(0);

// ***** NUMBER OF ITERATIONS *****
// const double iterations = 1e2;

// ***** FIXED PARAMETER *****

double c_mmusd = 299792.458;    // mm/us
mpreal c_mmus  = 299792.458;    // mm/us

mpreal mcs, mnu;

mpreal mxe_CF;
mpreal mxe_star;
mpreal mxe_defc;
mpreal mxe_phot;
mpreal mxe_corr;


// **************************
// ***** PARTICLE CLASS *****
// **************************
class FK_Particle{

public:
  // CONSTRUCTORS
  FK_Particle(){ };
  FK_Particle(double px, double py, double pz, double energy){
    p4[0] = px;
    p4[1] = py;
    p4[2] = pz;
    p4[3] = energy;
    SetMomentumGMP();
    SetMassGMP();
  };

  // "SET" METHODS
  void SetEnergy(double en)     {       // Sets FK_Particle energy
    p4[3] = en;
    gen = en;
  }

  // Methods which set the x,y,z; they update the magnitude of the momentum and the direction theta/phi
  void SetPx    (double px)      {p4[0] = px; gpx = px; UpdateThetaPhi();}
  void SetPy    (double py)      {p4[1] = py; gpy = py; UpdateThetaPhi();}
  void SetPz    (double pz)      {p4[2] = pz; gpz = pz; UpdateThetaPhi();}


  // Sets the vectorial FK_Particle of a 4-vector - used only in "matmult" function
  void SetP     (double qdm[3]) {
    for(int i=0;i<3;i++){p4[i] = qdm[i];}
    gpx = qdm[0];
    gpy = qdm[1];
    gpz = qdm[2];
    UpdateThetaPhi();}

  // Sets the vectorial FK_Particle of a 4-momentum one by one; updates the magnitude of the momentum and the direction theta/phi
  void SetP    (double px, double py, double pz){
    p4[0] = px; gpx = px;
    p4[1] = py; gpy = py;
    p4[2] = pz; gpz = pz;
    UpdateThetaPhi();
  }

  // Sets the 4 components of a 4-momentum one by one; updates the magnitude of the momentum and the direction theta/phi
  void SetP4    (double px, double py, double pz, double en){
    p4[0] = px;
    p4[1] = py;
    p4[2] = pz;
    p4[3] = en;
    gpx = px;
    gpy = py;
    gpz = pz;
    gen = en;
    SetMomentumGMP();
    SetMassGMP();
    UpdateThetaPhi();
  }

  // Sets the 4 components of a 4-momentum from a 4-vector; updates the magnitude of the momentum and the direction theta/phi
  void SetP4    (double qdm4[4]){
    for(int i=0;i<4;i++){
      p4[i] = qdm4[i];
    }
    gpx = qdm4[0];
    gpy = qdm4[1];
    gpz = qdm4[2];
    gen = qdm4[3];
    SetMomentumGMP();
    SetMassGMP();
    UpdateThetaPhi();
  }

  void SetP4GMP  (mpreal gmp_lpx, mpreal gmp_lpy, mpreal gmp_lpz, mpreal gmp_lenergy){
    gpx = gmp_lpx    ;
    gpy = gmp_lpy    ;
    gpz = gmp_lpz    ;
    gen = gmp_lenergy;
    gpt = sqrt(gpx*gpx + gpy*gpy + gpz*gpz);
    gmass = sqrt(gen*gen - gpt*gpt);
//     SetMomentumGMP();
//     SetMassGMP();
    p4[0] = gmp_lpx.toDouble();
    p4[1] = gmp_lpy.toDouble();
    p4[2] = gmp_lpz.toDouble();
    p4[3] = gmp_lenergy.toDouble();
    UpdateThetaPhi();
  }

  void SetTheta (double t)      {theta = t;}    // Sets theta
  void SetPhi   (double p)      {phi = p;}      // Sets phi

  void SetMomentumGMP(){ gpt = sqrt(gpx*gpx + gpy*gpy + gpz*gpz); }


  void SetMassGMP()      { gmass = sqrt(gen*gen - gpt*gpt); }

  // "GET" METHODS
  double GetPComp(int i)    {return p4[i];} // Returns the i-th component of a 4-momentum
  mpreal GetPxGMP()         {return gpx;}
  mpreal GetPyGMP()         {return gpy;}
  mpreal GetPzGMP()         {return gpz;}
  mpreal GetEnergyGMP()     {return gen;}
  double GetPxDouble()      {return p4[0];}
  double GetPyDouble()      {return p4[1];}
  double GetPzDouble()      {return p4[2];}
  double GetParrayDouble()  {return *p4;}
  double GetEnergyDouble()  {return p4[3];}
  double GetThetaDouble()   {return theta;}
  double GetPhiDouble()     {return phi;}
  double GetThetaGMP()      {return theta;}
  double GetPhiGMP()        {return phi;}
  mpreal GetMomentumGMP()   {return gpt;}
  double GetMomentumDouble(){return TMath::Sqrt( p4[0]*p4[0] + p4[1]*p4[1] + p4[2]*p4[2] );}
  mpreal GetMassGMP()       {return sqrt(gen*gen - gpt*gpt);}
  double GetMassDouble()    {return p4[3] * TMath::Sqrt(1 - GetMomentumDouble()/p4[3] * GetMomentumDouble()/p4[3]);}
  double GetVtotDouble()    {double vtot = (GetMomentumDouble() * c_mmusd)/ TMath::Sqrt( TMath::Power(GetMomentumDouble(),2) + TMath::Power(GetMassDouble(),2) );
                             if(GetMassDouble() == 0. || TMath::IsNaN(vtot) ){return c_mmusd;}
                             else return vtot;
                            }

  mpreal GetVtotGMP()       {mpreal vtot = (gpt*c_mmus)/sqrt(gpt*gpt + gmass*gmass);
                             if(gmass == 0. || mpfr::isnan(vtot) ){return c_mmus;}
                             else return vtot;
                            }


  mpreal GetVxGMP()         {return GetVtotGMP() * sin( gtheta ) * cos( gphi );}
  mpreal GetVyGMP()         {return GetVtotGMP() * sin( gtheta ) * sin( gphi );}
  mpreal GetVzGMP()         {return GetVtotGMP() * cos( gtheta );}

  double GetVxDouble()      {return GetVtotDouble() * sin( theta ) * cos( phi );}
  double GetVyDouble()      {return GetVtotDouble() * sin( theta ) * sin( phi );}
  double GetVzDouble()      {return GetVtotDouble() * cos( theta );}

  mpreal GetKEGMP()         {return (gen - gmass);}
  double GetKEDouble()      {mpreal ene = (gen - gmass); return ene.toDouble();}

  // Prints all the info for a FK_Particle
  void   PrintFourVector()    {
    std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(12);
    cout << "px = " << gpx << " keV/c\t; py = " << gpy << " keV/c\t; pz = " << gpz << " keV/c\t; energy = " << gen << " keV\t; ptot = " << gpt << " keV/c\t; mass = " << gmass << " keV/c²\t" << endl;
  }

  // Prints all the info for a FK_Particle
  void   PrintExpandedInfo()    {
    PrintFourVector();
    std::cout << "vx = " << GetVxGMP() << " mm/us\t; vy = " << GetVyGMP() << " mm/us\t; vz = " << GetVzGMP() << " mm/us\t; kinetic energy = " << GetKEGMP() << " keV\t; vtot = " << GetVtotGMP() << " mm/us\t" << std::endl;
  }

private:
  double p4[4] = {0.,0.,0.,0.}; // 4-momentum (px, py, pz, E)
  double theta; // theta angle
  double phi;   // phi angle
  double tof;   // time of flight
  mpreal gtheta; // theta angle
  mpreal gphi;   // phi angle
  mpreal gtof;   // time of flight
  mpreal gpx;
  mpreal gpy;
  mpreal gpz;
  mpreal gpt;
  mpreal gen;
  mpreal gke;
  mpreal gvx;
  mpreal gvy;
  mpreal gvz;
  mpreal gvt;
  mpreal gmass;

  // Theta-Phi updater - invoked at each change of momentum
  void UpdateThetaPhi(){
    gtheta = acos( gpz / gpt );
    gphi = atan2(gpy,gpx);
    gphi < 0 ? gphi += 2*TMath::Pi() : gphi;
    theta = gtheta.toDouble();
    phi = gphi.toDouble();
  }
};

// nu_pz = √[ mCs⁴ + (mnu² - mXe²)² - 2·mCs² · (mnu² + mXe²) ] / 2·mCs
mpreal GenerateNuPzGMP(mpreal mnu, mpreal mxe, mpreal mcs){ return ( sqrt( pow(mcs,4) + pow(mnu*mnu - mxe*mxe,2) - 2*mcs*mcs * (mnu*mnu + mxe*mxe) ) / (2*mcs) ); }

// en = √(pz² + mass²)
mpreal GenerateEnergyFromPMGMP(mpreal pz, mpreal mass){ return sqrt(pz*pz + mass*mass); }

FK_Particle *GenerateAugerGMP(mpreal auger_energy, mpreal auger_one_p){
  mpreal auger_one_theta = RandomTheta();
  mpreal auger_one_phi   = RandomPhi();
  FK_Particle *Auger_one = new FK_Particle();
  mpreal px, py, pz, en;
  px = auger_one_p * sin( auger_one_theta ) * cos( auger_one_phi ) ;
  py = auger_one_p * sin( auger_one_theta ) * sin( auger_one_phi ) ;
  pz = auger_one_p * cos( auger_one_theta )                               ;
  en = auger_energy;

  Auger_one->SetP4GMP( px, py, pz, en );

  return Auger_one;
}

FK_Particle *GenerateValidAugerGMP(mpreal auger_energy, mpreal auger_one_p){
  mpreal auger_one_theta = RandomValidAugerTheta();
  mpreal auger_one_phi   = RandomPhi();
  FK_Particle *Auger_one = new FK_Particle();
  mpreal px, py, pz, en;
  px = auger_one_p * sin( auger_one_theta ) * cos( auger_one_phi ) ;
  py = auger_one_p * sin( auger_one_theta ) * sin( auger_one_phi ) ;
  pz = auger_one_p * cos( auger_one_theta )                               ;
  en = auger_energy;

  Auger_one->SetP4GMP( px, py, pz, en );

  return Auger_one;
}

FK_Particle *GenerateXrayGMP(mpreal xray_energy, mpreal xray_p){
  mpreal xray_theta = RandomTheta();
  mpreal xray_phi   = RandomPhi();
  FK_Particle  *Xray = new FK_Particle();
  mpreal px, py, pz, en;
  px = xray_p * sin( xray_theta ) * cos( xray_phi ) ;
  py = xray_p * sin( xray_theta ) * sin( xray_phi ) ;
  pz = xray_p * cos( xray_theta )                          ;
  en = xray_energy;

  Xray->SetP4GMP( px, py, pz, en );

  return Xray;
}

FK_Particle *GenerateValidXrayGMP(mpreal xray_energy, mpreal xray_p){
  mpreal xray_theta = RandomValidXrayTheta();
  mpreal xray_phi   = RandomValidXrayPhi();
  FK_Particle  *Xray = new FK_Particle();
  mpreal px, py, pz, en;
  px = xray_p * sin( xray_theta ) * cos( xray_phi ) ;
  py = xray_p * sin( xray_theta ) * sin( xray_phi ) ;
  pz = xray_p * cos( xray_theta )                          ;
  en = xray_energy;

  Xray->SetP4GMP( px, py, pz, en );

  return Xray;
}


// ***** OTHER FUNCTIONS *****

mpreal RandomTheta()                      {return acos(2 * r->Uniform() - 1);                                   }
mpreal RandomValidXeStarTheta()           {return acos( sqrt(2)/2 + (1 - sqrt(2)/2) * r->Uniform() );           }
mpreal RandomPhi()                        {return (2 * TMath::Pi() * r->Uniform());                             }
mpreal RandomThetaWithErr(double error)   {return (acos(2 * r->Uniform() - 1) + r->Uniform(-1,1)*error);        }
mpreal RandomPhiWithErr(double error)     {return ((2 * TMath::Pi() * r->Uniform()) + r->Uniform(-1,1)*error);  }

mpreal RandomValidXrayTheta()             {
                                            bool valid = false;
                                            mpreal theTheta;
                                            while(valid == false){
                                              theTheta = acos(2 * r->Uniform() - 1);
//                                               if( (theTheta >= 0.92 && theTheta <= 1.51) ||
//                                                   (theTheta >= 1.63 && theTheta <= 2.22)){valid = true;}
                                              if( (theTheta >= 1.26 && theTheta <= 1.51) ||
                                                  (theTheta >= 1.63 && theTheta <= 1.88)){valid = true;}
                                            }
                                            return theTheta;
                                          }
mpreal RandomValidXrayPhi()               {
                                            bool valid = false;
                                            mpreal thePhi;
                                            while(valid == false){
                                              thePhi = 2 * TMath::Pi() * r->Uniform();
                                              if( (thePhi >= 1.22 && thePhi <= 1.93) ||
                                                  (thePhi >= 4.36 && thePhi <= 5.07)){valid = true;}
                                            }
                                            return thePhi;
                                          }

mpreal RandomValidAugerTheta()            {return acos(r->Uniform());                                           }

mpreal RandomThetaScatter(mpreal velocity)  {
  // TF1 for atom-ion scattering
  gRandom->SetSeed(0);
  mpreal mass = 2.1737066e-25;

  TF1 *ai_scatter = new TF1( "ai_scatter","( ( TMath::Sqrt( (3*TMath::Pi()*[0]) / (64.*[1]) ) * TMath::Power(1/(x+0.003),5./2.) + TMath::Sqrt( (1*TMath::Pi()*[0]) / (768.*[1]) ) * TMath::Power(1/(x+0.003),1./2.) + TMath::Sqrt( (49*TMath::Pi()*[0]) / (2764800.*[1]) ) * TMath::Power(1/(x+0.003),-3./2.) ) * TMath::Sin(x+0.003) ) + 1/(4. * [2] * TMath::Sqrt(8*[0]/[3])) ", 0, TMath::Pi() );
//  TF1 *ai_scatter = new TF1( "ai_scatter","TMath::Sqrt(3*TMath::Pi())/8. * TMath::Sqrt([0]/([1]*x)) * 1/(x)", 0, TMath::Pi() );
  mpreal C4 = 1.37055E-56;
  mpreal Erel = 0.5 * mass/2. * velocity*velocity;
//  ai_scatter->SetParameters(C4.toDouble() , Erel.toDouble());
  ai_scatter->SetParameters(C4.toDouble() , Erel.toDouble(), velocity.toDouble(), mass.toDouble()/2.);
  ai_scatter->SetNpx(100000);
  return ai_scatter->GetRandom();
}

mpreal *crossprod(mpreal *array1, mpreal *array2){
  mpreal *thearray = new mpreal[3];

  // Cross product
  thearray[0] = array1[1]*array2[2] - array1[2]*array2[1];
  thearray[1] = array1[2]*array2[0] - array1[0]*array2[2];
  thearray[2] = array1[0]*array2[1] - array1[1]*array2[0];

  mpreal sumsq = sqrt(thearray[0]*thearray[0] + thearray[1]*thearray[1] + thearray[2]*thearray[2]);

  for(int i=0;i<3;i++){thearray[i] /= sumsq;} // Normalization of the adaus

  return thearray;

}

void PrintGMPMatrix(mpreal_matrix *thematrix){
  cout << thematrix->r0c0 << " " << thematrix->r0c1 << " " << thematrix->r0c2 << endl;
  cout << thematrix->r1c0 << " " << thematrix->r1c1 << " " << thematrix->r1c2 << endl;
  cout << thematrix->r2c0 << " " << thematrix->r2c1 << " " << thematrix->r2c2 << endl;
}

mpreal_matrix TransposeGMPMatrix(mpreal_matrix *in){
  mpreal_matrix out;
  out.r0c0 = in->r0c0; out.r0c1 = in->r1c0; out.r0c2 = in->r2c0;
  out.r1c0 = in->r0c1; out.r1c1 = in->r1c1; out.r1c2 = in->r2c1;
  out.r2c0 = in->r0c2; out.r2c1 = in->r1c2; out.r2c2 = in->r2c2;
  return out;
}

mpreal_matrix InvertGMPMatrix(mpreal_matrix *in){
  mpreal_matrix out;

  // Determinant of the input matrix
  mpreal det =  (in->r0c0 * in->r1c1 * in->r2c2) +
                (in->r0c2 * in->r1c0 * in->r2c1) +
                (in->r0c1 * in->r1c2 * in->r2c0) -
                (in->r0c2 * in->r1c1 * in->r2c0) -
                (in->r0c1 * in->r1c0 * in->r2c2) -
                (in->r0c0 * in->r1c2 * in->r2c1);

  // Elements of the output matrix calculated as det(minor) / det
  out.r0c0 = ( (in->r1c1 * in->r2c2) - (in->r1c2 * in->r2c1) ) / det;
  out.r0c1 = ( (in->r1c0 * in->r2c2) - (in->r1c2 * in->r2c0) ) / det;
  out.r0c2 = ( (in->r1c0 * in->r2c1) - (in->r1c1 * in->r2c0) ) / det;

  out.r1c0 = ( (in->r0c1 * in->r2c2) - (in->r0c2 * in->r2c1) ) / det;
  out.r1c1 = ( (in->r0c0 * in->r2c2) - (in->r0c2 * in->r2c0) ) / det;
  out.r1c2 = ( (in->r0c0 * in->r2c1) - (in->r0c1 * in->r2c0) ) / det;

  out.r2c0 = ( (in->r0c1 * in->r1c2) - (in->r0c2 * in->r1c1) ) / det;
  out.r2c1 = ( (in->r0c0 * in->r1c2) - (in->r0c2 * in->r1c0) ) / det;
  out.r2c2 = ( (in->r0c0 * in->r1c1) - (in->r0c1 * in->r1c0) ) / det;

  return out;
}

mpreal_4vector GMP_matmult(FK_Particle &FK_Particle, mpreal_matrix *mat){
  mpreal_4vector out;
  mpreal gmp_lpx, gmp_lpy, gmp_lpz, gmp_len;

  gmp_lpx = FK_Particle.GetPxGMP(); gmp_lpy = FK_Particle.GetPyGMP(); gmp_lpz = FK_Particle.GetPzGMP(); gmp_len = FK_Particle.GetEnergyGMP();

  // Multiplication between Matrix' Rows and FK_Particle P4 spatial components
  out.x = (mat->r0c0*gmp_lpx) + (mat->r0c1*gmp_lpy) + (mat->r0c2*gmp_lpz);
  out.y = (mat->r1c0*gmp_lpx) + (mat->r1c1*gmp_lpy) + (mat->r1c2*gmp_lpz);
  out.z = (mat->r2c0*gmp_lpx) + (mat->r2c1*gmp_lpy) + (mat->r2c2*gmp_lpz);
  out.e = gmp_len;

  return out;
}

void GMP_rotate4vector(mpreal_matrix *rotmatrix, FK_Particle &theoriginal, FK_Particle &therotated){
  mpreal_4vector matresult = GMP_matmult(theoriginal, rotmatrix);
  therotated.SetP4GMP(matresult.x, matresult.y, matresult.z, matresult.e);
}


#endif
