///////////////////////////////////////////////////////////////////////////////////
// full_kinematics_mpreal.cc (v5.0) - Author: Francesco Granato
// 2017-10-06
//
// 5.0  - Changed single calls to SIMION, decay by decay, to two calls, for 4 and 5 particles events
// 4.1  - Added WRITE_MOMENTA_TO_FILE to export momenta
// 4.0  - Converted every quantity to mpreal hi-precision variables
// 3.0  - Introduced particle vectors
// 2.0  - Added GenerateAuger() and GenerateXray() functions
// 1.0  - Initial draft
///////////////////////////////////////////////////////////////////////////////////

// includes for kinematics
# include "full_kinematics_mpreal.hh"
# include "TF1.h"
# include "TMath.h"
# include <iostream>
# include <fstream>
# include <stdio.h>

// includes to run external program (SIMION)
# include <sstream>
# include <cstdlib>
# include <cstring>

// includes for ROOT reconstruction
# include "TROOT.h"
# include "TSystem.h"
# include "TSystemDirectory.h"
# include <vector>
# include "TString.h"
# include "TFile.h"
# include "TTree.h"
# include "TH1D.h"
# include "TF1.h"

# include "no-single-double_scatter.hh" // defines if scattering is used
# include "hires_lores.hh"  // defines if we're using HiRes-LoAcc or LoRes-HiAcc tune

// -----------------------------------------------

// # define RANDOMIZE_CS_POSITION

// # define WRITE_MOMENTA_TO_FILE
# define WRITE_FLY2
# define SIMION
// # define HADD_FILES
# define TFULL_WRITER


// PRINT DEFINITIONS
// # define PRINT_CS
// # define PRINT_XE_STAR_ORIG
// # define PRINT_NU_ORIG
// # define PRINT_XE_STAR_ROTATED
// # define PRINT_XE_STAR_ROTATED_SCATTERED
// # define PRINT_SCATTER_STUFF
// # define PRINT_NU_ROTATED
// # define PRINT_AUGER_ONE
// # define PRINT_X_RAY
// # define PRINT_XE_PLUS
// # define PRINT_RECONSTRUCTED_XE_PLUS_MASS
// # define PRINT_RECONSTRUCTED_NU_MASS
// # define PRINT_CS_RANDOM_XYZ
// # define PRINT_NEW_AXIS
// # define PRINT_ROT_MATRIX
// # define PRINT_TRANSVERSE_ROT_MATRIX
// # define PRINT_INVERTED_ROT_MATRIX

using namespace std;

void full_kinematics(double startentry, double stopentry){
  gRandom->SetSeed(0);
  gROOT->ProcessLine("#include <vector>"); // To include <vector> in ROOT
  mpfr_set_default_prec(256);

  gSystem->cd("/home/tug26830/full_kinematics/");

  # ifdef WRITE_FLY2
  TString fly2_dir    = TString::Format("fly2_%06.f-%06.f",         startentry, stopentry);
  gROOT->ProcessLine( TString::Format(".! mkdir -p /home/tug26830/work/fk_output/%s", fly2_dir.Data())    );
  ofstream *out_fly2 = new ofstream(); out_fly2->open(TString::Format("/home/tug26830/work/fk_output/%s/out.fly2", fly2_dir.Data())); *out_fly2 << fixed << setprecision(20);
  //     ofstream *out_fly2 = new ofstream(); out_fly2->open(TString::Format("fly2/out_%s_%d.fly2",e_type.Data(),i+1)); *out_fly2 << fixed << setprecision(20);

  int nparticles;
  ofstream *out_nparts = new ofstream(); out_nparts->open(TString::Format("/home/tug26830/work/fk_output/%s/out_nparts.txt", fly2_dir.Data()));
  *out_fly2 << "particles {\n";
  # endif

  int nscat = 0;
  # ifdef SINGLE_ION_SCATTER
  nscat = 1;
  # endif
  # ifdef DOUBLE_ION_SCATTER
  nscat = 2;
  # endif

  # ifdef SIMION
  TString input_dir   = TString::Format("input_txt_%06.f-%06.f",    startentry, stopentry);
  TString output_dir  = TString::Format("output_root_%06.f-%06.f",  startentry, stopentry);
  TString momenta_dir = TString::Format("momenta_%06.f-%06.f",      startentry, stopentry);
  gROOT->ProcessLine( TString::Format(".! mkdir -p /home/tug26830/work/fk_output/%s", input_dir.Data())   );
  gROOT->ProcessLine( TString::Format(".! mkdir -p /home/tug26830/work/fk_output/%s", output_dir.Data())  );

  TSystemDirectory *dir = new TSystemDirectory( output_dir.Data(), output_dir.Data() ); int nfiles;

  # ifdef WRITE_MOMENTA_TO_FILE
  gROOT->ProcessLine( TString::Format(".! mkdir -p /home/tug26830/work/fk_output/%s", momenta_dir.Data()) );
  # endif
  # endif

  # ifdef WRITE_MOMENTA_TO_FILE
  gROOT->ProcessLine( TString::Format(".! rm /home/tug26830/work/fk_output/%s/*", momenta_dir.Data()) ); // */
  # endif

  // *** FK_Particle parameters ***
  // Everything is in keV
  mpreal CF                 = 931494.095;                               // AMU -> keV

  // Electrons
  mpreal auger_mass         =    510.998910;                            // electron mass
  mpreal O1_energy          =      0.0233;                              // Xenon O1 binding energy - source
  mpreal O23_energy         =      0.0128;                              // Xenon O2,3 binding energy (average between the two)

  // Neutrino
  mpreal nu_mass            =      0.;                                  // neutrino mass

  // X-Ray
  mpreal xray_en            =     34.415;                               // X-ray energy (KM2 transition)
  mpreal xray_p             =      xray_en;                             // X-ray energy (KM2 transition)
  mpreal xray_kin           =      xray_en;                             // X-ray energy (KM2 transition)
  mpreal xray_mass          =      0.;                                  // X-ray mass

  // Cesium
  mpreal cs_mass_defect     = -88058.9;                                 // Cs mass defect
  mpreal cs_mass            =    131. * CF + cs_mass_defect;            // Cs mass (ground state)
  mpreal cs_qv              =    354.752;                               // Cs -> Xe + nu Q-value
  mpreal cs_mass_err        =      5.343e-6 * CF;                       // error on Cs mass
  mpreal cs_mass_def_err    =      4.977;                               // error on Cs mass defect
  mpreal cs_init_xyz[3]     = {   0.,
                                  0.,
                               1090.};                                  // initial Cs position (mm - SIMION coords)
  mpreal cs_source_size     =     1.;                                   // source radius (mm)

  // Xenon
  mpreal xe_mass_defect     = -88413.6;                                 // Xe mass defect
  mpreal xe_photoabs_en     =     34.565;                               // K edge for Xe
  mpreal xe_mass            =    131. * CF + xe_mass_defect;            // Xe mass (ground state)
  mpreal xe_mass_corr       =      0.0038939;                           // correction to Xe mass after K-capture - energy of Cs P-shell (taken from http://www.chembio.uoguelph.ca/educmat/atomdata/bindener/grp1num.htm )
  mpreal xe_ioniz_en        =      0.01213;                             // ionization energy
  mpreal xe_mass_err        =      0.236e-6 * CF;                       // error on Xe mass
  mpreal xe_mass_def_err    =      0.;                                  // error on Xe mass defect

  // Version modified by the factor "expected nparts / observed nparts"
  // TF1 *ai_scatter = new TF1( "ai_scatter","( ( ( TMath::Sqrt( (3*TMath::Pi()*[0]) / (64.*[1]) ) * TMath::Power(1/(x+0.003),5./2.) + TMath::Sqrt( (1*TMath::Pi()*[0]) / (768.*[1]) ) * TMath::Power(1/(x+0.003),1./2.) + TMath::Sqrt( (49*TMath::Pi()*[0]) / (2764800.*[1]) ) * TMath::Power(1/(x+0.003),-3./2.) ) * TMath::Sin(x+0.003) ) + 1/(4. * [2]) * TMath::Sqrt(8*[0]/[3]) )*0.00407936", 0, TMath::Pi() );

  // Unmodified version
  //   TF1 *ai_scatter = new TF1( "ai_scatter","( ( TMath::Sqrt( (3*TMath::Pi()*[0]) / (64.*[1]) ) * TMath::Power(1/(x+0.003),5./2.) + TMath::Sqrt( (1*TMath::Pi()*[0]) / (768.*[1]) ) * TMath::Power(1/(x+0.003),1./2.) + TMath::Sqrt( (49*TMath::Pi()*[0]) / (2764800.*[1]) ) * TMath::Power(1/(x+0.003),-3./2.) ) * TMath::Sin(x+0.003) ) + 1/(4. * [2]) * TMath::Sqrt(8*[0]/[3])", 0, TMath::Pi() );

  // Analytical version
  // TF1 *ai_scatter = new TF1( "ai_scatter","TMath::Sqrt(3*TMath::Pi())/8. * TMath::Sqrt([0]/([1]*(x+0.003))) * 1/((x+0.003))", 0, TMath::Pi() );

  // Analytical version - renormalized
  // TF1 *ai_scatter = new TF1( "ai_scatter","(TMath::Sqrt(3*TMath::Pi())/8. * TMath::Sqrt([0]/([1]*(x+0.003))) * 1/((x+0.003)))*0.00359243", 0, TMath::Pi() );

  // Analytical version - piecewise function from quantum and classical cross section + LANGEVIN
  TF1 *ai_quant = new TF1( "ai_quant","(TMath::Sqrt(3*TMath::Pi())/8. * TMath::Sqrt([0]/[1]) * 1/(TMath::Power(0.003,2)*TMath::Sin(0.003))) * TMath::Sqrt(x)*TMath::Sin(x)", 0, TMath::Pi() );
  TF1 *ai_class = new TF1( "ai_class","(TMath::Sqrt(3*TMath::Pi())/8. * TMath::Sqrt([0]/[1]) * 1/(TMath::Power(x,1.5)))", 0, TMath::Pi() );
  TF1 *ai_lange = new TF1( "ai_lange","0.25*TMath::Sqrt(4*[0]/[1])*TMath::Sin(x)", 0, TMath::Pi() );
  TF1 *ai_scatter = new TF1( "ai_scatter","(x<0.003)*(ai_quant + ai_lange) + (x>=0.003)*(ai_class + ai_lange)", 0, TMath::Pi() );

  // very low cut
//  TF1 *ai_scatter = new TF1( "ai_scatter","(x<0.00001)*(ai_quant + ai_lange) + (x>=0.00001)*(ai_class + ai_lange)", 0, TMath::Pi() );

  // Only classical + Langevin, to generate events only over a certain minimum
//  double minsmith = 0.003; // == 1 (keV/c2)2
//  double minsmith = 0.0008; // == 0.07 (keV/c2)2
//  TF1 *ai_scatter = new TF1( "ai_scatter", "(ai_class + ai_lange)", minsmith, TMath::Pi() );

//   ai_scatter->SetNpx(3500);
  ai_scatter->SetNpx(10000);
  mpreal C4 = 6.8508667E-57;
  mpreal mass = 2.1737066e-25;

  // Statement for setting the precision of output numbers
  std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(20);

  mpreal mxe_star = 131. * CF + xe_mass_defect + xe_photoabs_en - xe_mass_corr;

  // VARIABLES DEFINITION

  double process_probability;
  int process_selector = 0;
  int nauger    = 0;
  double temp_auger_energy;

  vector<mpreal> *Augers_momen_GMP = new vector<mpreal>;
  vector<mpreal> *Augers_binds_GMP = new vector<mpreal>;

  TString e_type;
  TString e_name;

  stringstream stream;

  // The Cs is standing still in the origin of the system; the origin is randomly taken
  // Cs
  FK_Particle *Cs = new FK_Particle();
  Cs->SetP4GMP(0.,0.,0.,cs_mass);

  # ifdef PRINT_CS
  cout << "Cs" << endl;
  Cs->PrintExpandedInfo();
  cout << endl;
  # endif

  // The geometry is closed, so the neutrino will always have a definite energy and momentum
  mpreal nupz       = GenerateNuPzGMP(nu_mass, mxe_star, cs_mass);    // pz of Xe* and nu
  mpreal xe_orig_en = GenerateEnergyFromPMGMP(nupz, mxe_star);        // resulting energy for Xe*
  mpreal nu_orig_en = GenerateEnergyFromPMGMP(nupz, nu_mass);         // resulting energy for nu

  // The initial decay of Cs will leave the Xe* and the neutrino
  FK_Particle *Xe_star_orig = new FK_Particle();
  FK_Particle *Nu_orig      = new FK_Particle();

  Xe_star_orig->SetP4GMP(0, 0,  nupz, xe_orig_en );
  Nu_orig     ->SetP4GMP(0, 0, -nupz, nu_orig_en );

  # ifdef PRINT_XE_STAR_ORIG
  cout << "Xe*_orig" << endl;
  Xe_star_orig->PrintExpandedInfo();
  cout << endl;
  # endif

  # ifdef PRINT_NU_ORIG
  cout << "Nu_orig" << endl;
  Nu_orig->PrintExpandedInfo();
  cout << endl;
  # endif

  // This ends the first, "closed" part of the reaction

  // Second part, Auger and x-ray generation

  int niter = (int)stopentry - (int)startentry;

  int step_delta, step;
  // Loop on events in the chain
  if(niter > 10.){
    step_delta = -2; // smaller values (even < 0) = faster counter; larger values = slower counter
    step = TMath::Power(10,step_delta + TMath::Floor(TMath::Log10(niter)));    // Event counter
  }
  for(int i=startentry; i<stopentry; ++i){
//     if(niter > 10.){ if((i-(int)startentry)%step == 0){cout << (i-(int)startentry) << "/" << niter << endl;} }


    // Random process selection
    process_probability = r->Uniform();
    if(process_probability >= 0.00 && process_probability < 0.03) { process_selector = 0; }   // 100 eV Auger               - prob. 0.03
    if(process_probability >= 0.03 && process_probability < 0.09) { process_selector = 1; }   // 111 eV Auger               - prob. 0.06
    if(process_probability >= 0.09 && process_probability < 0.25) { process_selector = 2; }   // 122 eV Auger               - prob. 0.16
    if(process_probability >= 0.25 && process_probability < 0.34) { process_selector = 3; }   //  54 eV Auger + 34 eV Auger - prob. 0.09
    if(process_probability >= 0.34 && process_probability < 0.52) { process_selector = 4; }   //  54 eV Auger + 45 eV Auger - prob. 0.18
    if(process_probability >= 0.52 && process_probability < 0.74) { process_selector = 5; }   //  65 eV Auger + 45 eV Auger - prob. 0.22
    if(process_probability >= 0.74 && process_probability < 0.88) { process_selector = 6; }   //  65 eV Auger + 34 eV Auger - prob. 0.14
    if(process_probability >= 0.88 && process_probability < 1.00) { process_selector = 7; }   //  65 eV Auger + 23 eV Auger - prob. 0.12

    // MANUAL OVERRIDE FOR THE PROCESS SELECTION
//     process_selector = 0;

    # ifdef WRITE_FLY2
    // This defines how many particles in the fly2 belong to the decay
    if(process_selector <= 2){nparticles = 4+nscat;}
    else{nparticles = 5+nscat;}

    *out_nparts << nparticles << "\n";
    # endif

    if      (process_selector == 0){    // 100 eV Auger               - prob. 0.03
      e_type = "sa100";
      e_name = "100 eV Auger";
      nauger    = 1;

      temp_auger_energy = 100.;

      Augers_momen_GMP->push_back(sqrt(temp_auger_energy*temp_auger_energy + 2*auger_mass*temp_auger_energy*1e3)*1e-3);
      Augers_binds_GMP->push_back(O23_energy);
    }
    else if (process_selector == 1){    // 111 eV Auger               - prob. 0.06
      e_type = "sa111";
      e_name = "111 eV Auger";
      nauger    = 1;

      temp_auger_energy = 111.;

      Augers_momen_GMP->push_back(sqrt(temp_auger_energy*temp_auger_energy + 2*auger_mass*temp_auger_energy*1e3)*1e-3);
      Augers_binds_GMP->push_back(O23_energy);
    }
    else if (process_selector == 2){    // 122 eV Auger               - prob. 0.16
      e_type = "sa122";
      e_name = "122 eV Auger";
      nauger    = 1;

      temp_auger_energy = 122.;

      Augers_momen_GMP->push_back(sqrt(temp_auger_energy*temp_auger_energy + 2*auger_mass*temp_auger_energy*1e3)*1e-3);
      Augers_binds_GMP->push_back(O23_energy);
    }
    else if (process_selector == 3){    //  54 eV Auger + 34 eV Auger - prob. 0.09
      e_type = "da54_34";
      e_name = "54 eV + 34 eV Augers";
      nauger    = 2;

      temp_auger_energy = 54.;

      Augers_momen_GMP->push_back(sqrt(temp_auger_energy*temp_auger_energy + 2*auger_mass*temp_auger_energy*1e3)*1e-3);
      Augers_binds_GMP->push_back(O23_energy);


      temp_auger_energy = 34.;

      Augers_momen_GMP->push_back(sqrt(temp_auger_energy*temp_auger_energy + 2*auger_mass*temp_auger_energy*1e3)*1e-3);
      Augers_binds_GMP->push_back(O23_energy);
    }
    else if (process_selector == 4){    //  54 eV Auger + 45 eV Auger - prob. 0.18
      e_type = "da54_45";
      e_name = "54 eV + 45 eV Augers";
      nauger    = 2;

      temp_auger_energy = 54.;

      Augers_momen_GMP->push_back(sqrt(temp_auger_energy*temp_auger_energy + 2*auger_mass*temp_auger_energy*1e3)*1e-3);
      Augers_binds_GMP->push_back(O23_energy);


      temp_auger_energy = 45.;

      Augers_momen_GMP->push_back(sqrt(temp_auger_energy*temp_auger_energy + 2*auger_mass*temp_auger_energy*1e3)*1e-3);
      Augers_binds_GMP->push_back(O23_energy);
    }
    else if (process_selector == 5){    //  65 eV Auger + 45 eV Auger - prob. 0.22
      e_type = "da65_45";
      e_name = "65 eV + 45 eV Augers";
      nauger    = 2;

      temp_auger_energy = 65.;

      Augers_momen_GMP->push_back(sqrt(temp_auger_energy*temp_auger_energy + 2*auger_mass*temp_auger_energy*1e3)*1e-3);
      Augers_binds_GMP->push_back(O23_energy);


      temp_auger_energy = 45.;

      Augers_momen_GMP->push_back(sqrt(temp_auger_energy*temp_auger_energy + 2*auger_mass*temp_auger_energy*1e3)*1e-3);
      Augers_binds_GMP->push_back(O23_energy);
    }
    else if (process_selector == 6){    //  65 eV Auger + 34 eV Auger - prob. 0.14
      e_type = "da65_34";
      e_name = "65 eV + 34 eV Augers";
      nauger    = 2;

      temp_auger_energy = 65.;

      Augers_momen_GMP->push_back(sqrt(temp_auger_energy*temp_auger_energy + 2*auger_mass*temp_auger_energy*1e3)*1e-3);
      Augers_binds_GMP->push_back(O23_energy);


      temp_auger_energy = 34.;

      Augers_momen_GMP->push_back(sqrt(temp_auger_energy*temp_auger_energy + 2*auger_mass*temp_auger_energy*1e3)*1e-3);
      Augers_binds_GMP->push_back(O23_energy);
    }
    else if (process_selector == 7){    //  65 eV Auger + 23 eV Auger - prob. 0.12
      e_type = "da65_23";
      e_name = "65 eV + 23 eV Augers";
      nauger    = 2;

      temp_auger_energy = 65.;

      Augers_momen_GMP->push_back(sqrt(temp_auger_energy*temp_auger_energy + 2*auger_mass*temp_auger_energy*1e3)*1e-3);
      Augers_binds_GMP->push_back(O23_energy);


      temp_auger_energy = 23.;

      Augers_momen_GMP->push_back(sqrt(temp_auger_energy*temp_auger_energy + 2*auger_mass*temp_auger_energy*1e3)*1e-3);
      Augers_binds_GMP->push_back(O23_energy);
    }

    // The Cs initial position is randomized
    mpreal cs_xyz[3];
    for(int i=0;i<3;++i){cs_xyz[i] = cs_init_xyz[i];}

    # ifdef RANDOMIZE_CS_POSITION
    // Gaussian distribution of the three components of the initial position
    TF1 *xgaus = new TF1("xgaus","gaus",cs_init_xyz[0].toDouble()-10, cs_init_xyz[0].toDouble()+10); xgaus->SetNpx(1000);
    TF1 *ygaus = new TF1("ygaus","gaus",cs_init_xyz[1].toDouble()-10, cs_init_xyz[1].toDouble()+10); ygaus->SetNpx(1000);
    TF1 *zgaus = new TF1("zgaus","gaus",cs_init_xyz[2].toDouble()-10, cs_init_xyz[2].toDouble()+10); zgaus->SetNpx(1000);
    double sigma  = 0.682;  // sigma that for a 1e9 atoms gives a peak density of 2e8/mm³
    xgaus->SetParameters(1, cs_init_xyz[0].toDouble(), sigma);
    ygaus->SetParameters(1, cs_init_xyz[1].toDouble(), sigma);
    zgaus->SetParameters(1, cs_init_xyz[2].toDouble(), sigma);
    cs_xyz[0] = xgaus->GetRandom();
    cs_xyz[1] = ygaus->GetRandom();
    cs_xyz[2] = zgaus->GetRandom();

   // old uniform spherical distribution
//    mpreal cs_theta               = RandomTheta();
//    mpreal cs_phi                 = RandomPhi();
//    mpreal cs_radius              = cs_source_size * TMath::Power(r->Uniform(), 1./3.);
//    cs_xyz[0] += cs_radius * sin( cs_theta ) * cos( cs_phi );
//    cs_xyz[1] += cs_radius * sin( cs_theta ) * sin( cs_phi );
//    cs_xyz[2] += cs_radius * cos( cs_theta );
    # endif

    # ifdef PRINT_CS_RANDOM_XYZ
    cout << "Cs random position:" << endl;
    cout << cs_xyz[0] << " " << cs_xyz[1] << " " << cs_xyz[2] << endl;
    cout << endl;
    # endif

    // Now one can generate a random system of axes to apply a random rotation to the vectors
    // One generates two random versors, than calculates the cross product between them;
    // then calculates the cross product between the first one and the result

    // *** Random coordinate system generation ***

    mpreal A1_theta = RandomTheta();
    mpreal   A1_phi = RandomPhi();
    mpreal  D_theta = RandomTheta(); // Dummy vector for Axis2 generation
    mpreal    D_phi = RandomPhi();   // via cross product with Axis1

    mpreal A1[3] = { sin(A1_theta) * cos(A1_phi), sin(A1_theta) * sin(A1_phi), cos(A1_theta) };
    mpreal  D[3] = { sin(D_theta)  * cos(D_phi) , sin(D_theta)  * sin(D_phi) , cos(D_theta)  };

    mpreal *A2 = crossprod(A1,D );      // Axis2
    mpreal *A3 = crossprod(A1,A2);      // Axis3

    # ifdef PRINT_NEW_AXIS
    cout << "New rotated coordinate system" << endl;
    cout << " Axis 1 - i = " << A1[0] << " ; j = " << A1[1] << " ; k = " << A1[2] << endl;
    cout << " Axis 2 - i = " << A2[0] << " ; j = " << A2[1] << " ; k = " << A2[2] << endl;
    cout << " Axis 3 - i = " << A3[0] << " ; j = " << A3[1] << " ; k = " << A3[2] << endl;

    cout << endl << "Orthogonality test" << endl;
    cout << " A1·A2 = " << A1[0]*A2[0] + A1[1]*A2[1] + A1[2]*A2[2] << endl;
    cout << " A1·A3 = " << A1[0]*A3[0] + A1[1]*A3[1] + A1[2]*A3[2] << endl;
    cout << " A2·A3 = " << A2[0]*A3[0] + A2[1]*A3[1] + A2[2]*A3[2] << endl;
    cout << endl;
    # endif

    mpreal_matrix GMP_trans;

    GMP_trans.r0c0 = A1[0]; GMP_trans.r0c1 = A1[1]; GMP_trans.r0c2 = A1[2];
    GMP_trans.r1c0 = A2[0]; GMP_trans.r1c1 = A2[1]; GMP_trans.r1c2 = A2[2];
    GMP_trans.r2c0 = A3[0]; GMP_trans.r2c1 = A3[1]; GMP_trans.r2c2 = A3[2];

    # ifdef PRINT_ROT_MATRIX
    cout << "Rotation matrix from ex, ey, ez to rotated system" << endl;
    PrintGMPMatrix(&GMP_trans);
    cout << endl;
    # endif

    mpreal_matrix GMP_transT = TransposeGMPMatrix(&GMP_trans);

    # ifdef PRINT_TRANSVERSE_ROT_MATRIX
    cout << "Transposed rotation matrix from ex, ey, ez to rotated system" << endl;
    PrintGMPMatrix(&GMP_transT);
    cout << endl;
    # endif

    GMP_transT = InvertGMPMatrix(&GMP_transT);

    # ifdef PRINT_INVERTED_ROT_MATRIX
    cout << "Inverted rotation matrix from ex, ey, ez to rotated system" << endl;
    PrintGMPMatrix(&GMP_transT);
    cout << endl;
    # endif

    FK_Particle *Xe_star = new FK_Particle();
    FK_Particle *Nu      = new FK_Particle();

    GMP_rotate4vector(&GMP_transT, *Xe_star_orig, *Xe_star);
    GMP_rotate4vector(&GMP_transT, *Nu_orig,      *Nu);

    # ifdef PRINT_XE_STAR_ROTATED
    cout << "Xe*" << endl;
    Xe_star->PrintExpandedInfo();
    cout << endl;
    # endif

    # ifdef PRINT_NU_ROTATED
    cout << "Nu" << endl;
    Nu->PrintExpandedInfo();
    cout << endl;
    # endif

    // An Xray and an Auger are always generated
    // The Auger has random generation angles
    vector<FK_Particle*> *Augers = new vector<FK_Particle*>;
    for(int i=0;i<nauger;++i){
      FK_Particle *Auger_one = GenerateAugerGMP( sqrt(Augers_momen_GMP->at(i)*Augers_momen_GMP->at(i) + auger_mass*auger_mass), Augers_momen_GMP->at(i));
      Augers->push_back(Auger_one);
    }

    # ifdef PRINT_AUGER_ONE
    cout << "Augers" << endl;
    for(int i=0;i<nauger;++i){
      Augers->at(i)->PrintExpandedInfo();
      cout << endl;
    }
    # endif

    // The X-ray has random generation angles
    FK_Particle *Xray = GenerateValidXrayGMP( xray_en, xray_p );
    # ifdef PRINT_X_RAY
    cout << "Xray" << endl;
    Xray->PrintExpandedInfo();
    cout << endl;
    # endif



    // The Xe* goes into Xe+, and its components are calculated from the ones of the Xe*, the Auger and the X-ray (conservation of 4-momentum)
    FK_Particle *Xe_plus = new FK_Particle();
    mpreal total_Augers_px     = 0;
    mpreal total_Augers_py     = 0;
    mpreal total_Augers_pz     = 0;
    mpreal total_Augers_energy = 0;

    for(int i=0;i<nauger;++i){
      total_Augers_px     += Augers->at(i)->GetPxGMP();
      total_Augers_py     += Augers->at(i)->GetPyGMP();
      total_Augers_pz     += Augers->at(i)->GetPzGMP();
      total_Augers_energy += Augers->at(i)->GetEnergyGMP();
    }

    Xe_plus->SetP4GMP( Xe_star->GetPxGMP()     - total_Augers_px     - Xray->GetPxGMP(),
                       Xe_star->GetPyGMP()     - total_Augers_py     - Xray->GetPyGMP(),
                       Xe_star->GetPzGMP()     - total_Augers_pz     - Xray->GetPzGMP(),
                       Xe_star->GetEnergyGMP() - total_Augers_energy - Xray->GetEnergyGMP() );

    # ifdef PRINT_XE_PLUS
    cout << "Xe+" << endl;
    Xe_plus->PrintExpandedInfo();
    cout << endl;
    # endif

    # ifdef PRINT_RECONSTRUCTED_XE_PLUS_MASS
    if(nauger == 1){
      xe_mass_plus = xe_mass - auger_mass + Augers_bindings->at(0);
      cout << "Calculated  mass of the Xe+ = " << zmass << endl;
      cout << "Theoretical mass of the Xe+ = " << xe_mass_plus << endl;
      cout << "Delta = " << zmass - xe_mass_plus << endl;
    }
    if(nauger == 2){
      xe_mass_plusplus = xe_mass - 2*auger_mass + Augers_bindings->at(0) + Augers_bindings->at(1);
      cout << "Calculated  mass of the Xe++ = " << zmass << endl;
      cout << "Theoretical mass of the Xe++ = " << xe_mass_plusplus << endl;
      cout << "Delta = " << zmass - xe_mass_plusplus << endl;
    }
    cout << endl;
    # endif

    // mν² = En² - Px² - Py² - Pz²
    mpreal GMP_En_sq = pow(cs_mass - Xe_plus->GetEnergyGMP() - total_Augers_energy - Xray->GetEnergyGMP(), 2);
    mpreal GMP_Px_sq = pow(          Xe_plus->GetPxGMP()     + total_Augers_px     + Xray->GetPxGMP(),     2);
    mpreal GMP_Py_sq = pow(          Xe_plus->GetPyGMP()     + total_Augers_py     + Xray->GetPyGMP(),     2);
    mpreal GMP_Pz_sq = pow(          Xe_plus->GetPzGMP()     + total_Augers_pz     + Xray->GetPzGMP(),     2);

    mpreal GMP_Rec_mnu2 = GMP_En_sq - GMP_Px_sq - GMP_Py_sq - GMP_Pz_sq;
    mpreal GMP_Rec_mnu  = sqrt(GMP_Rec_mnu2);

    # ifdef PRINT_RECONSTRUCTED_NU_MASS
    // Reconstructed neutrino mass
    cout << "Reconstructed nu mass² (keV²) = ";
    cout << GMP_Rec_mnu2 << endl;
    cout << "Reconstructed nu mass (keV) = ";
    cout << GMP_Rec_mnu << endl;
    # endif

    mpreal xe_modplab = Xe_plus->GetMomentumGMP();
    mpreal xe_plab[3] = { Xe_plus->GetPxGMP(), Xe_plus->GetPyGMP(), Xe_plus->GetPzGMP() };
    mpreal xe_energy = Xe_plus->GetEnergyGMP();

    # ifdef SINGLE_ION_SCATTER
    FK_Particle *Xe_sing = new FK_Particle();
    mpreal plab_to_pcm = cs_mass / (cs_mass + Xe_plus->GetMassGMP());
    mpreal pcm_to_plab = Xe_plus->GetMassGMP() / (cs_mass + Xe_plus->GetMassGMP());

    // Conversion from lab frame to COM frame
    mpreal xe_pcm[3]  = { Xe_plus->GetPxGMP() * plab_to_pcm, Xe_plus->GetPyGMP() * plab_to_pcm, Xe_plus->GetPzGMP() * plab_to_pcm };
    mpreal xe_modpcm  = sqrt( xe_pcm[0]*xe_pcm[0] + xe_pcm[1]*xe_pcm[1] + xe_pcm[2]*xe_pcm[2] ) ;

    // NOTE: the ai_scatter gives us theta angles, and we want that these angles
    // determine the (very small) aperture of a scattering cone in the direction of the
    // nonscattered Xenon.
    // To this end, we have to move from the COM to a reference frame in which the
    // z-axis is collinear to the nonscattered Xenon. Doing so requires finding the
    // matrix M so that:
    // P' = M x P
    // where P=(px,py,pz) is in the COM, and P'=(0,0,|P|) is in the new collinear
    // frame.
    // xhat, yhat and zhat construct this matrix. zhat is constructed so that the
    // z-axis is collinear to P.
    // Once the theta_scatter and phi_scatter angles have been generated, we have to
    // move back to the COM frame (and then to the lab frame). To do so, one should
    // use:
    // P = M^-1 x P'
    // So we should invert the matrix M. BUT! The matrix M is UNITARY, and the inverse
    // of a unitary matrix is its TRANSPOSED!
    // So when building the three components of xe_pcm_sc[3], we perform the matrix
    // product between the transposed M and P'(theta_scatter, phi_scatter)!

    // Building the M matrix rows
    mpreal zhat[3]    = { xe_pcm[0]/xe_modpcm, xe_pcm[1]/xe_modpcm, xe_pcm[2]/xe_modpcm };
    mpreal modzhat    = sqrt( zhat[0]*zhat[0] + zhat[1]*zhat[1] + zhat[2]*zhat[2] );
    mpreal xhat[3]    = { -zhat[2]/sqrt(zhat[0]*zhat[0]+zhat[2]*zhat[2]), 0., zhat[0]/sqrt(zhat[0]*zhat[0]+zhat[2]*zhat[2]) };
    mpreal *yhat      = crossprod(zhat, xhat);

    mpreal velocity = Xe_plus->GetVtotGMP()*1000.;
    mpreal Erel = 0.5 * mass/2. * velocity*velocity;
    // ai_scatter->SetParameters(C4.toDouble() , Erel.toDouble(), C4.toDouble() , Erel.toDouble());
    ai_scatter->SetParameters(C4.toDouble() , Erel.toDouble(), C4.toDouble() , Erel.toDouble(), C4.toDouble() , Erel.toDouble(), C4.toDouble() , Erel.toDouble());
    mpreal theta_scatter = ai_scatter->GetRandom();
    mpreal   phi_scatter = RandomPhi();

    // Here there's the matrix product M^T x P'. Remember that M^-1 == M^T
    mpreal xe_pcm_sc[3] = {
      xe_modpcm * ( xhat[0] * sin(theta_scatter)*cos(phi_scatter) + yhat[0] * sin(theta_scatter)*sin(phi_scatter) + zhat[0]*cos(theta_scatter) ),
      xe_modpcm * ( xhat[1] * sin(theta_scatter)*cos(phi_scatter) + yhat[1] * sin(theta_scatter)*sin(phi_scatter) + zhat[1]*cos(theta_scatter) ),
      xe_modpcm * ( xhat[2] * sin(theta_scatter)*cos(phi_scatter) + yhat[2] * sin(theta_scatter)*sin(phi_scatter) + zhat[2]*cos(theta_scatter) )
    };

    mpreal xe_modpcm_sc = sqrt( xe_pcm_sc[0]*xe_pcm_sc[0] + xe_pcm_sc[1]*xe_pcm_sc[1] +xe_pcm_sc[2]*xe_pcm_sc[2] );

    // Conversion from COM frame to lab frame
    mpreal xe_plab_sc[3] = {
      xe_pcm_sc[0] + pcm_to_plab*xe_plab[0],
      xe_pcm_sc[1] + pcm_to_plab*xe_plab[1],
      xe_pcm_sc[2] + pcm_to_plab*xe_plab[2]
    };
    mpreal xe_modplab_sc = sqrt( xe_plab_sc[0]*xe_plab_sc[0] + xe_plab_sc[1]*xe_plab_sc[1] +xe_plab_sc[2]*xe_plab_sc[2] );

    Xe_sing->SetP4GMP( xe_plab_sc[0], xe_plab_sc[1], xe_plab_sc[2], sqrt(xe_modplab_sc*xe_modplab_sc + Xe_plus->GetMassGMP()*Xe_plus->GetMassGMP()) );

    # ifdef PRINT_SCATTER_STUFF
    cout << endl;
    cout << "xe_plab" << endl;
    cout << "px = " << Xe_plus->GetPxGMP()  << " keV/c\t; py = " << Xe_plus->GetPyGMP()  << " keV/c\t; pz = " << Xe_plus->GetPzGMP()  << " keV/c\t; ptot = " << Xe_plus->GetMomentumGMP() << endl;
    cout << "theta_lab = " << acos(Xe_plus->GetPzGMP()/Xe_plus->GetMomentumGMP()) << " ; phi_lab = " << atan2(Xe_plus->GetPyGMP(),Xe_plus->GetPxGMP()) << endl;
    cout << endl;

    cout << "xe_pcm" << endl;
    cout << "px = " << xe_pcm[0]  << " keV/c\t; py = " << xe_pcm[1]  << " keV/c\t; pz = " << xe_pcm[2]  << " keV/c\t; ptot = " << xe_modpcm << endl;
    cout << "theta_cm = " << acos(xe_pcm[2]/xe_modpcm) << " ; phi_cm = " << atan2(xe_pcm[1],xe_pcm[0]) << endl;
    cout << endl;

//     cout << "zhat" << endl;
//     cout << "px = " << zhat[0]  << " keV/c\t; py = " << zhat[1]  << " keV/c\t; pz = " << zhat[2]  << " keV/c\t; ptot = " << modzhat << endl;
//     cout << endl;
//
//     cout << "xhat" << endl;
//     cout << "px = " << xhat[0]  << " keV/c\t; py = " << xhat[1]  << " keV/c\t; pz = " << xhat[2]  << " keV/c\t" << endl;
//     cout << endl;
//
//     cout << "yhat" << endl;
//     cout << "px = " << yhat[0]  << " keV/c\t; py = " << yhat[1]  << " keV/c\t; pz = " << yhat[2]  << " keV/c\t" << endl;
//     cout << endl;
//
//     cout << "Orthogonality check:" << endl;
//     cout << "zhat·xhat: zhat[0]*xhat[0] + zhat[1]*xhat[1] + zhat[2]*xhat[2] = \t" << zhat[0]*xhat[0] + zhat[1]*xhat[1] + zhat[2]*xhat[2] << endl;
//     cout << "zhat·yhat: zhat[0]*yhat[0] + zhat[1]*yhat[1] + zhat[2]*yhat[2] = \t" << zhat[0]*yhat[0] + zhat[1]*yhat[1] + zhat[2]*yhat[2] << endl;
//     cout << "xhat·yhat: xhat[0]*yhat[0] + xhat[1]*yhat[1] + xhat[2]*yhat[2] = \t" << xhat[0]*yhat[0] + xhat[1]*yhat[1] + xhat[2]*yhat[2] << endl;
//     cout << endl;

    cout << "theta_sc = " << theta_scatter << " ; phi_sc = " << phi_scatter << endl << endl;

    cout << "xe_pcm_sc" << endl;
    cout << "px = " << xe_pcm_sc[0]  << " keV/c\t; py = " << xe_pcm_sc[1]  << " keV/c\t; pz = " << xe_pcm_sc[2]  << " keV/c\t; ptot = " << xe_modpcm_sc << endl;
    cout << "theta_cm_sc = " << acos(xe_pcm_sc[2]/xe_modpcm_sc) << " ; phi_cm = " << atan2(xe_pcm_sc[1],xe_pcm_sc[0]) << endl;
    cout << endl;

    cout << "Xe* scattered" << endl;
    cout << "px = " << xe_plab_sc[0]  << " keV/c\t; py = " << xe_plab_sc[1]  << " keV/c\t; pz = " << xe_plab_sc[2]  << " keV/c\t; ptot = " << xe_modplab_sc << endl;
    cout << "theta_lab_sc = " << acos(xe_plab_sc[2]/xe_modplab_sc) << " ; phi_lab = " << atan2(xe_plab_sc[1],xe_plab_sc[0]) << endl;
    cout << endl;

    cout << "Deltaptot = " << xe_modplab - xe_modplab_sc << endl << endl;

    cout << "**********************************************************************************" << endl << endl;
    # endif
    # endif

    # ifdef DOUBLE_ION_SCATTER
    FK_Particle *Xe_dobl = new FK_Particle();
    plab_to_pcm = cs_mass / (cs_mass + Xe_sing->GetMassGMP());
    pcm_to_plab = Xe_sing->GetMassGMP() / (cs_mass + Xe_sing->GetMassGMP());

    xe_pcm[0]  = Xe_sing->GetPxGMP() * plab_to_pcm;
    xe_pcm[1]  = Xe_sing->GetPyGMP() * plab_to_pcm;
    xe_pcm[2]  = Xe_sing->GetPzGMP() * plab_to_pcm;
    xe_modpcm  = sqrt( xe_pcm[0]*xe_pcm[0] + xe_pcm[1]*xe_pcm[1] + xe_pcm[2]*xe_pcm[2] ) ;

    zhat[0]    = xe_pcm[0]/xe_modpcm;
    zhat[1]    = xe_pcm[1]/xe_modpcm;
    zhat[2]    = xe_pcm[2]/xe_modpcm;
    modzhat    = sqrt( zhat[0]*zhat[0] + zhat[1]*zhat[1] + zhat[2]*zhat[2] );
    xhat[0]    = -zhat[2]/modzhat;
    xhat[1]    = 0.;
    xhat[2]    = zhat[0]/modzhat;
    mpreal *dyhat      = crossprod(zhat, xhat);

    velocity = Xe_plus->GetVtotGMP()*1000.;
    Erel = 0.5 * mass/2. * velocity*velocity;
    ai_scatter->SetParameters(C4.toDouble() , Erel.toDouble(), velocity.toDouble(), mass.toDouble()/2.);
    theta_scatter = ai_scatter->GetRandom();
//     theta_scatter = RandomThetaScatter(Xe_sing->GetVtotGMP()*1000.);
    phi_scatter = RandomPhi();

    xe_pcm_sc[0] =   xe_modpcm * ( xhat[0] * sin(theta_scatter)*cos(phi_scatter) + dyhat[0] * sin(theta_scatter)*sin(phi_scatter) + zhat[0]*cos(theta_scatter) );
    xe_pcm_sc[1] =   xe_modpcm * ( xhat[1] * sin(theta_scatter)*cos(phi_scatter) + dyhat[1] * sin(theta_scatter)*sin(phi_scatter) + zhat[1]*cos(theta_scatter) );
    xe_pcm_sc[2] =   xe_modpcm * ( xhat[2] * sin(theta_scatter)*cos(phi_scatter) + dyhat[2] * sin(theta_scatter)*sin(phi_scatter) + zhat[2]*cos(theta_scatter) );

    xe_modpcm_sc = sqrt( xe_pcm_sc[0]*xe_pcm_sc[0] + xe_pcm_sc[1]*xe_pcm_sc[1] +xe_pcm_sc[2]*xe_pcm_sc[2] );

    xe_plab_sc[0] =  xe_pcm_sc[0] + pcm_to_plab*xe_plab[0];
    xe_plab_sc[1] =  xe_pcm_sc[1] + pcm_to_plab*xe_plab[1];
    xe_plab_sc[2] =  xe_pcm_sc[2] + pcm_to_plab*xe_plab[2];
    xe_modplab_sc = sqrt( xe_plab_sc[0]*xe_plab_sc[0] + xe_plab_sc[1]*xe_plab_sc[1] +xe_plab_sc[2]*xe_plab_sc[2] );

    Xe_dobl->SetP4GMP( xe_plab_sc[0], xe_plab_sc[1], xe_plab_sc[2], sqrt(xe_modplab_sc*xe_modplab_sc + Xe_sing->GetMassGMP()*Xe_sing->GetMassGMP()) );

    # ifdef PRINT_SCATTER_STUFF
    cout << endl;
    cout << "xe_plab" << endl;
    cout << "px = " << Xe_sing->GetPxGMP()  << " keV/c\t; py = " << Xe_sing->GetPyGMP()  << " keV/c\t; pz = " << Xe_sing->GetPzGMP()  << " keV/c\t; ptot = " << Xe_sing->GetMomentumGMP() << endl;
    cout << "theta_lab = " << acos()
    cout << endl;

    cout << "xe_pcm" << endl;
    cout << "px = " << xe_pcm[0]  << " keV/c\t; py = " << xe_pcm[1]  << " keV/c\t; pz = " << xe_pcm[2]  << " keV/c\t; ptot = " << xe_modpcm << endl;
    cout << endl;

    cout << "zhat" << endl;
    cout << "px = " << zhat[0]  << " keV/c\t; py = " << zhat[1]  << " keV/c\t; pz = " << zhat[2]  << " keV/c\t; ptot = " << modzhat << endl;
    cout << endl;

    cout << "xhat" << endl;
    cout << "px = " << xhat[0]  << " keV/c\t; py = " << xhat[1]  << " keV/c\t; pz = " << xhat[2]  << " keV/c\t" << endl;
    cout << endl;

    cout << "yhat" << endl;
    cout << "px = " << dyhat[0]  << " keV/c\t; py = " << dyhat[1]  << " keV/c\t; pz = " << dyhat[2]  << " keV/c\t" << endl;
    cout << endl;

    cout << "Orthogonality check:" << endl;
    cout << "zhat·xhat: zhat[0]*xhat[0] + zhat[1]*xhat[1] + zhat[2]*xhat[2] = \t" << zhat[0]*xhat[0] + zhat[1]*xhat[1] + zhat[2]*xhat[2] << endl;
    cout << "zhat·yhat: zhat[0]*yhat[0] + zhat[1]*yhat[1] + zhat[2]*yhat[2] = \t" << zhat[0]*dyhat[0] + zhat[1]*dyhat[1] + zhat[2]*dyhat[2] << endl;
    cout << "xhat·yhat: xhat[0]*yhat[0] + xhat[1]*yhat[1] + xhat[2]*yhat[2] = \t" << xhat[0]*dyhat[0] + xhat[1]*dyhat[1] + xhat[2]*dyhat[2] << endl;
    cout << endl;

    cout << "xe_pcm_sc" << endl;
    cout << "px = " << xe_pcm_sc[0]  << " keV/c\t; py = " << xe_pcm_sc[1]  << " keV/c\t; pz = " << xe_pcm_sc[2]  << " keV/c\t; ptot = " << xe_modpcm_sc << endl;
    cout << endl;

    cout << "Xe* scattered" << endl;
    cout << "px = " << xe_plab_sc[0]  << " keV/c\t; py = " << xe_plab_sc[1]  << " keV/c\t; pz = " << xe_plab_sc[2]  << " keV/c\t; ptot = " << xe_modplab_sc << endl;
    cout << endl;

    cout << "Deltaptot = " << xe_modplab - xe_modplab_sc << endl;
    # endif
    # endif

    # ifdef WRITE_MOMENTA_TO_FILE
    // Writes a txt file containing Cs initial position, and momenta and mass for every particle
    ofstream *out_momenta = new ofstream(); out_momenta->open( TString::Format("/home/tug26830/work/fk_output/momenta_%06.f-%06.f/mom_%06d.txt", startentry, stopentry, i+1) ); *out_momenta << fixed << setprecision(20);
    *out_momenta << cs_xyz [0] << "\t" << cs_xyz [1] << "\t" << cs_xyz [2] << "\n";
    *out_momenta << Xe_plus->GetPxGMP() << "\t" << Xe_plus->GetPyGMP() << "\t" << Xe_plus->GetPzGMP() << "\t" << Xe_plus->GetMomentumGMP() << "\t" << Xe_plus->GetMassGMP() << "\t" << Xe_plus->GetEnergyGMP() << "\n";
    # ifdef SINGLE_ION_SCATTER
    *out_momenta << Xe_sing->GetPxGMP() << "\t" << Xe_sing->GetPyGMP() << "\t" << Xe_sing->GetPzGMP() << "\t" << Xe_sing->GetMomentumGMP() << "\t" << Xe_sing->GetMassGMP() << "\t" << Xe_sing->GetEnergyGMP() << "\n";
    # endif
    # ifdef DOUBLE_ION_SCATTER
    *out_momenta << Xe_dobl->GetPxGMP() << "\t" << Xe_dobl->GetPyGMP() << "\t" << Xe_dobl->GetPzGMP() << "\t" << Xe_dobl->GetMomentumGMP() << "\t" << Xe_dobl->GetMassGMP() << "\t" << Xe_dobl->GetEnergyGMP() << "\n";
    # endif
    *out_momenta << Xray->GetPxGMP() << "\t" << Xray->GetPyGMP() << "\t" << Xray->GetPzGMP() << "\t" << Xray->GetMomentumGMP() << "\t" << Xray->GetMassGMP() << "\t" << Xray->GetEnergyGMP() << "\n";
    *out_momenta << Nu->GetPxGMP() << "\t" << Nu->GetPyGMP() << "\t" << Nu->GetPzGMP() << "\t" << Nu->GetMomentumGMP() << "\t" << Nu->GetMassGMP() << "\t" << Nu->GetEnergyGMP() << "\n";
    for(int i=0;i<nauger;++i){
    *out_momenta << Augers->at(i)->GetPxGMP() << "\t" << Augers->at(i)->GetPyGMP() << "\t" << Augers->at(i)->GetPzGMP() << "\t" << Augers->at(i)->GetMomentumGMP() << "\t" << Augers->at(i)->GetMassGMP() << "\t" << Augers->at(i)->GetEnergyGMP() << "\n";
    }
    out_momenta->close();
    # endif

    # ifdef WRITE_FLY2
    *out_fly2 << "  standard_beam { ";
    *out_fly2 << "    mass = " << Xe_plus->GetMassGMP()/CF << ", ";
    if(nauger == 1){
    *out_fly2 << "    charge = 1, ";
    }
    if(nauger == 2){
    *out_fly2 << "    charge = 2, ";
    }
    *out_fly2 << "    x = " << cs_xyz[0] << " , y = " << cs_xyz[1] << " , z = " << cs_xyz[2] << ", ";
    *out_fly2 << "    direction = vector(" << Xe_plus->GetPxGMP()/Xe_plus->GetMomentumGMP() << " , " << Xe_plus->GetPyGMP()/Xe_plus->GetMomentumGMP() << " , " << Xe_plus->GetPzGMP()/Xe_plus->GetMomentumGMP() << "), ";
    *out_fly2 << "    momentum = " << Xe_plus->GetMomentumGMP()*1e-3 << ", ";
    *out_fly2 << "    color = 0, ";
    *out_fly2 << "  },\n";
    # ifdef SINGLE_ION_SCATTER
    *out_fly2 << "  standard_beam { ";
    *out_fly2 << "    mass = " << Xe_sing->GetMassGMP()/CF << ", ";
    if(nauger == 1){
    *out_fly2 << "    charge = 1, ";
    }
    if(nauger == 2){
    *out_fly2 << "    charge = 2, ";
    }
    *out_fly2 << "    x = " << cs_xyz[0] << " , y = " << cs_xyz[1] << " , z = " << cs_xyz[2] << ", ";
    *out_fly2 << "    direction = vector(" << Xe_sing->GetPxGMP()/Xe_sing->GetMomentumGMP() << " , " << Xe_sing->GetPyGMP()/Xe_sing->GetMomentumGMP() << " , " << Xe_sing->GetPzGMP()/Xe_sing->GetMomentumGMP() << "), ";
    *out_fly2 << "    momentum = " << Xe_sing->GetMomentumGMP()*1e-3 << ", ";
    *out_fly2 << "    color = 8, ";
    *out_fly2 << "  },\n";
    # endif
    # ifdef DOUBLE_ION_SCATTER
    *out_fly2 << "  standard_beam { ";
    *out_fly2 << "    mass = " << Xe_dobl->GetMassGMP()/CF << ", ";
    if(nauger == 1){
    *out_fly2 << "    charge = 1, ";
    }
    if(nauger == 2){
    *out_fly2 << "    charge = 2, ";
    }
    *out_fly2 << "    x = " << cs_xyz[0] << " , y = " << cs_xyz[1] << " , z = " << cs_xyz[2] << ", ";
    *out_fly2 << "    direction = vector(" << Xe_dobl->GetPxGMP()/Xe_dobl->GetMomentumGMP() << " , " << Xe_dobl->GetPyGMP()/Xe_dobl->GetMomentumGMP() << " , " << Xe_dobl->GetPzGMP()/Xe_dobl->GetMomentumGMP() << "), ";
    *out_fly2 << "    momentum = " << Xe_dobl->GetMomentumGMP()*1e-3 << ", ";
    *out_fly2 << "    color = 0, ";
    *out_fly2 << "  },\n";
    # endif
    *out_fly2 << "  standard_beam { ";
    *out_fly2 << "    mass = 1E-10, "; // manually changed to a very low number so that SIMION flies it
    *out_fly2 << "    charge = 0, ";
    *out_fly2 << "    x = " << cs_xyz[0] << " , y = " << cs_xyz[1] << " , z = " << cs_xyz[2] << ", ";
    *out_fly2 << "    direction = vector(" << Xray->GetPxGMP() / Xray->GetMomentumGMP() << " , " << Xray->GetPyGMP() / Xray->GetMomentumGMP() << " , " << Xray->GetPzGMP() / Xray->GetMomentumGMP() << "), ";
    *out_fly2 << "    momentum = " << Xray->GetMomentumDouble()*1e-3 << ", ";
    *out_fly2 << "    color = 1, ";
    *out_fly2 << "  },\n";
    *out_fly2 << "  standard_beam { ";
    if(nu_mass > 0.){
      *out_fly2 << "    mass = " << nu_mass/CF << ", ";
    }
    else{
      *out_fly2 << "    mass = 1E-10, "; // manually changed to a very low number so that SIMION flies it
    }
    *out_fly2 << "    charge = 0, ";
    *out_fly2 << "    x = " << cs_xyz[0] << " , y = " << cs_xyz[1] << " , z = " << cs_xyz[2] << ", ";
    *out_fly2 << "    direction = vector(" << Nu->GetPxGMP() / Nu->GetMomentumGMP() << " , " << Nu->GetPyGMP() / Nu->GetMomentumGMP() << " , " << Nu->GetPzGMP() / Nu->GetMomentumGMP() << "), ";
    *out_fly2 << "    momentum = " << Nu->GetMomentumGMP()*1e-3 << ", ";
    *out_fly2 << "    color = 2, ";
    *out_fly2 << "  },\n";
    // Augers
    for(int i=0;i<nauger;++i){
      *out_fly2 << "  standard_beam { ";
      *out_fly2 << "    mass = 0.00054857990946, ";
      *out_fly2 << "    charge = -1, ";
      *out_fly2 << "    x = " << cs_xyz[0] << " , y = " << cs_xyz[1] << " , z = " << cs_xyz[2] << ", ";
      *out_fly2 << "    direction = vector(" << Augers->at(i)->GetPxGMP() / Augers->at(i)->GetMomentumGMP() << " , " << Augers->at(i)->GetPyGMP() / Augers->at(i)->GetMomentumGMP() << " , " << Augers->at(i)->GetPzGMP() / Augers->at(i)->GetMomentumGMP() << "), ";
      *out_fly2 << "    momentum = " << Augers->at(i)->GetMomentumGMP()*1e-3 << ", ";
      *out_fly2 << "  },\n";
//       if(i == 0 && nauger > 1){ *out_fly2 << ",\n"; }
//       else{*out_fly2 << "\n";}
    }
    # endif

    // This hadds files, so that no enormous number of output files is produced
    # ifdef HADD_FILES
    nfiles = dir->GetListOfFiles()->GetSize();
    if(nfiles > 100.){

      stream.str(std::string());  // used to clear the stringstream after each use
      stream << "hadd /home/tug26830/work/fk_output/" << output_dir.Data() << "/temp.tmp";
      stream << " " << "/home/tug26830/work/fk_output/" << output_dir.Data() << "/kinem_*";
      stream << " " << ">>";
      stream << " " << "/dev/null";
      system( stream.str().c_str() );

      stream.str(std::string());  // used to clear the stringstream after each use
      stream << " " << "rm";
      stream << " " << "/home/tug26830/work/fk_output/" << output_dir.Data() << "/*.root"; // */
      system( stream.str().c_str() );

      stream.str(std::string());  // used to clear the stringstream after each use
      stream << " " << "mv /home/tug26830/work/fk_output/" << output_dir.Data() << "/temp.tmp ";
      stream << " " << TString::Format("/home/tug26830/work/fk_output/%s/kinem_%06.f-%06.f_temp.root", output_dir.Data(), startentry, stopentry);
      system( stream.str().c_str() );
    }
    # endif

    Augers_momen_GMP->clear();
    Augers_binds_GMP->clear();

    # ifdef RANDOMIZE_CS_POSITION
    delete xgaus;
    delete ygaus;
    delete zgaus;
    # endif
    delete Xe_star;
    delete Nu;
    delete Augers;
    delete Xe_plus;
    # ifdef SINGLE_ION_SCATTER
    delete Xe_sing;
    # endif
    # ifdef DOUBLE_ION_SCATTER
    delete Xe_dobl;
    # endif
    # ifdef WRITE_MOMENTA_TO_FILE
    delete out_momenta;
    # endif
  }

  # ifdef WRITE_FLY2
  *out_fly2 << "} ";
  out_fly2->close();

  out_nparts->close();
  # endif

  // This invokes the SIMION routine
  # ifdef SIMION
  gROOT->ProcessLine( TString::Format(".! rm /home/tug26830/work/fk_output/%s/*", input_dir.Data()) ); // */
  stream.str(std::string());  // used to clear the stringstream after each use

//   stream << "/home/garnet/.wine/drive_c/SIMION-8.1/simion_orig";    // local computer
//   stream << "/home/tug26830/SIMION-8.1/simion";                     // compute.temple.edu
  stream << "/home/tug26830/bin/simion";                     // owlsnest.hpc.temple.edu
  stream << " " << "--nogui";
  stream << " " << "fly";
  stream << " " << "--particles=/home/tug26830/work/fk_output/" << fly2_dir.Data() << "/out.fly2";
  //     stream << " " << "--particles=./fly2/out_" << e_type.Data() << "_" << i+1 << ".fly2";
  stream << " " << "--recording-output=/home/tug26830/work/fk_output/" << input_dir.Data() << "/kinem.txt";
  stream << " " << "--retain-trajectories=0";

  # ifdef LORES_HIACC_4100V
  stream << " " << "--recording=./1090mm_Asym8G_RealLoRes_HiAccept_4100Vmm/1090mm_Asym8G_RealLoRes_HiAccept_4100Vmm.rec"; // LoRes HiAcceptance tune record file
  stream << " " << "1090mm_Asym8G_RealLoRes_HiAccept_4100Vmm.iob";                                                        // LoRes HiAcceptance tune .iob file
  # endif

  # ifdef HIRES_LOACC_1236V
  stream << " " << "--recording=./1090mm_Asym8G_RealHiRes_LoAccept_1236Vmm/1090mm_Asym8G_RealHiRes_LoAccept_1236Vmm.rec"; // HiRes LoAcceptance tune record file
  stream << " " << "1090mm_Asym8G_RealHiRes_LoAccept_1236Vmm.iob";                                                        // HiRes LoAcceptance tune .iob file
  # endif

  stream << " " << ">>";
  stream << " " << "/dev/null";
  system( stream.str().c_str() );
  # endif

  # ifdef HADD_FILES
  stream.str(std::string());  // used to clear the stringstream after each use
  stream << "hadd /home/tug26830/work/fk_output/" << output_dir.Data() << "/temp.tmp";
  stream << " " << "/home/tug26830/work/fk_output/" << output_dir.Data() << "/kinem_*";
  stream << " " << ">>";
  stream << " " << "/dev/null";
  system( stream.str().c_str() );

  stream.str(std::string());  // used to clear the stringstream after each use
  stream << " " << "rm";
  stream << " " << "/home/tug26830/work/fk_output/" << output_dir.Data() << "/*.root"; // */
  system( stream.str().c_str() );

  stream.str(std::string());  // used to clear the stringstream after each use
  stream << " " << "mv /home/tug26830/work/fk_output/" << output_dir.Data() << "/temp.tmp ";
  stream << " " << TString::Format("/home/tug26830/work/fk_output/%s/kinem_%06.f-%06.f.root", output_dir.Data(), startentry, stopentry);
  system( stream.str().c_str() );
  # endif


  // This invokes the tfull_writer routine
  # ifdef TFULL_WRITER
  stream.str(std::string());  // used to clear the stringstream after each use
  stream << "./tfull_writer_1090_nocomment_onlypasscuts_file_multiple.cc.exe";
  stream << " " << "/home/tug26830/work/fk_output/" << input_dir.Data() << "/kinem.txt";
  stream << " " << TString::Format("%06.f", startentry) << " " << TString::Format("%06.f", stopentry);
  stream << " " << TString::Format("%06d", (int)stopentry);
  system( stream.str().c_str() );
  # endif
}

#ifndef __CINT__

int main(int argc, char **argv) {
  if ( argc == 3 ) {
    full_kinematics( atoi(argv[1]), atoi(argv[2]) );
  }
  else {
    cout << "Usage:" << endl;
    cout << "./full_kinematics.cc.exe startentry stopentry [generates (startentry-stopentry) decays]" << endl;
    return 0;
  }


  return 0;
}
#endif /* __CINT __ */
