# include <vector>
# include "TROOT.h"
# include "TString.h"
# include "TFile.h"
# include "TF1.h"
# include "TTree.h"
# include "TMath.h"
# include "TGraph2D.h"
# include "TGraph.h"
# include "TChain.h"
# include "TPolyLine3D.h"
# include "TPolyMarker3D.h"
# include "TPolyMarker.h"
# include "TCanvas.h"
# include "TGLViewer.h"
# include "TAxis.h"
# include "TPad.h"
# include "TEllipse.h"
# include "TH1F.h"
# include "TRandom3.h"
# include <iostream>
# include <fstream>

// General variables

double CX_SR_LOW(3.), CX_SR_HIGH(60.);
// double CX_TOF(0.165); // 1400
// double CX_TOF(0.42); // 12 V
// double CX_TOF(0.38); // 1538 V
double CX_TOF(0.425); // 1206 V
// double CX_TOF(0.35); // 1206 V
double CX_NODE(0.01);
//   double CX_TOF(0.120); // 846 V

// double nodes[] = {0.0446439, 0.0892878, 0.133932, 0.178576, 0.22322 , 0.267864, 0.312507, 0.357151, 0.401795, 0.446439, 0.491083, 0.535727, 0.580371, 0.625015, 0.669659, 0.714303, 0.758947, 0.803591, 0.848234, 0.892878, 0.937522}; // 2000mm, E=0.12
double nodes[] = {0.0441362, 0.0882725, 0.132409, 0.176545, 0.220681, 0.264817, 0.308954, 0.35309 , 0.397226, 0.441362, 0.485499, 0.529635, 0.573771, 0.617907, 0.662044, 0.70618 , 0.750316, 0.794452, 0.838589, 0.882725, 0.926861}; // 2000mm, E=0.1206, real B
// double nodes[] = {0.0446548, 0.0893097, 0.133965, 0.178619, 0.223274, 0.267929, 0.312584, 0.357239, 0.401894, 0.446548, 0.491203, 0.535858, 0.580513, 0.625168, 0.669823, 0.714477, 0.759132, 0.803787, 0.848442, 0.893097, 0.937752}; // 2000mm, E=0.12, real B
// double nodes[] = {0.0444327, 0.0888654, 0.133298, 0.177731, 0.222163, 0.266596, 0.311029, 0.355461, 0.399894, 0.444327, 0.488759, 0.533192, 0.577625, 0.622057, 0.66649, 0.710923, 0.755356, 0.799788, 0.844221, 0.888654, 0.933086}; // 1209mm, E=0.12

// PARTICLE STRUCTURE

double CF = 931494.095; // Conversion factor amu -> keV

double c = 299792458; // m/s

double Bfield = 8.094;   // Gauss
// double Bfield = 6.0;   // Gauss
double auger_mass = 510.998910;   // keV/c²
double omega =  8.98755* 1e3 * Bfield / auger_mass;    // us^-1

double DEsplat_radius_MCP = 0.05;   // mm = 50 um
double DEtof_MCP = 0.001;   // us = 1 ns

// Origin position (mm)
double origin_x = 990.;
double origin_y = 0.;
double origin_z = 0.;

double vessel_inner_radius = 300.;

double MCP_radius = 60.;

double ion_MCP_position = -91.;

double src_inner_radius = 70.;
double src_outer_radius = src_inner_radius + 25.;

double ele_inner_radius = 171.;
double ele_outer_radius = ele_inner_radius + 75.;
double ele_MCP_position = origin_x + 2000.;


// TFile *isos = new TFile("isos_1090218113.root");
// TFile *isos = new TFile("isos_1269141.root");
TFile *isos = new TFile("isos_253941.root");

TGraph2D *iso_init_prad     = (TGraph2D*)isos->Get("iso_init_prad");
TGraph2D *iso_init_plong    = (TGraph2D*)isos->Get("iso_init_plong");
TGraph2D *iso_init_azimuth  = (TGraph2D*)isos->Get("iso_init_azimuth");

// ION SECTION

// ion
TFile *ion_isos1 = new TFile("ion_data1.root");
TGraph2D *iso_f_plong1  = (TGraph2D*)ion_isos1->Get("plon_iso_curve");

TFile *ion_isos2 = new TFile("ion_data2.root");
TGraph2D *iso_f_plong2  = (TGraph2D*)ion_isos2->Get("plon_iso_curve");

// ion trend line parameters

// ion charge = 1e
const double SLOP_TOF_TO_MCP1          = 36.9505918635/1000  ;                 // microsecond/(keV/c)
const double TOF_INTERCEPTION_TO_MCP1  = 199.8486561762  ;

const double SLOP_TD_TO_MCP1           = 502.9849488062/1000  ;                // millimeter/(keV/c)
const double TD_TO_MCP_INTERCEPTION1   = 0.1556559785  ;
const double SLOP_TD_OPS_MCP1          = 501.5807497755/1000  ;                // millimeter/(keV/c)
const double TD_OPS_MCP_INTERCEPTION1  = 0.3776078339 ;
const double TOF_CRITERION1            = TOF_INTERCEPTION_TO_MCP1 ;

// ion charge = 2e
const double SLOP_TOF_TO_MCP2          = 23.200989991623/1000  ;                 // microsecond/(keV/c)
const double TOF_INTERCEPTION_TO_MCP2  = 141.30564661252  ;

const double SLOP_TD_TO_MCP2           = 449.27100024715/1000  ;                // millimeter/(keV/c)
const double TD_TO_MCP_INTERCEPTION2   = 0.090053233216185  ;
const double SLOP_TD_OPS_MCP2          = 450.48953465832/1000  ;                // millimeter/(keV/c)
const double TD_OPS_MCP_INTERCEPTION2  = 0.15222402685268 ;
const double TOF_CRITERION2            = TOF_INTERCEPTION_TO_MCP2 ;

// TF1 *scatter = new TF1("scatter","(x<0.035)*[0]*[1] + (x>0.035)*[1]*([2]/TMath::Power(x,2.5) + [3]/TMath::Power(x,1.5) + [0])", 0, TMath::Pi());
// TF1 *scatter = new TF1("scatter","[0]*([1]/TMath::Power(x+[2],2.5) + [3]/TMath::Power(x+[2],1.5) + [4])", 0, TMath::Pi());

TRandom3 *ra = new TRandom3(0);

// Definition of the "particle" class
class particle{

public:
  double  ionn;           // ion number
  double  tof;            // tof [us]
  double  mass;           // mass of the particle [amu]
  double  charge;         // charge of the particle [e]
  double  init_x;         // initial x position [mm]
  double  init_y;         // initial y position [mm]
  double  init_z;         // initial z position [mm]
  double  init_azm;       // initial azimuth [deg]
  double  init_elev;      // initial angle with the X axis [deg] [NOTE: IT'S DIFFERENT FROM ELEVATION! GOES FROM 0° TO 90°!]
  double  init_dircos;    // initial director cosine
  double  init_vx;        // initial vx (X-component of the velocity) [mm/us]
  double  init_vy;        // initial vy (Y-component of the velocity) [mm/us]
  double  init_vz;        // initial vz (Z-component of the velocity) [mm/us]
  double  init_vtot;      // initial total velocity [mm/us]
  double  init_vrad;      // initial vrad (radial component of the velocity = √(vy²+vz²) ) [mm/us]
  double  init_vlong;     // initial vlong (longitudinal component of the velocity = vx) [mm/us]
  double  init_px;        // initial px (X-component of the momentum) [keV/c]
  double  init_py;        // initial py (Y-component of the momentum) [keV/c]
  double  init_pz;        // initial pz (Z-component of the momentum) [keV/c]
  double  init_ptot;      // initial total momentum [keV/x]
  double  init_prad;      // initial prad (radial component of the momentum = √(py²+pz²) ) [keV/c]
  double  init_plong;     // initial plong (longitudinal component of the momentum = px) [keV/c]
  double  init_ke;        // initial kinetic energy (keV)
  double  splat_x;        // splat x position [mm]
  double  splat_y;        // splat y position [mm]
  double  splat_z;        // splat z position [mm]
  double  splat_azm;      // splat azimuth [deg]
  double  splat_elev;     // splat angle with the X axis [deg] [NOTE: IT'S DIFFERENT FROM ELEVATION! GOES FROM 0° TO 90°!]
  double  splat_dircos;   // splat director cosine
  double  splat_vx;       // splat vx (X-component of the velocity) [mm/us]
  double  splat_vy;       // splat vy (Y-component of the velocity) [mm/us]
  double  splat_vz;       // splat vz (Z-component of the velocity) [mm/us]
  double  splat_vtot;     // splat total velocity [mm/us]
  double  splat_vrad;     // splat vrad (radial component of the velocity = √(vy²+vz²) ) [mm/us]
  double  splat_vlong;    // splat vlong (longitudinal component of the velocity = vx) [mm/us]
  double  splat_px;       // splat px (X-component of the momentum) [keV/c]
  double  splat_py;       // splat py (Y-component of the momentum) [keV/c]
  double  splat_pz;       // splat pz (Z-component of the momentum) [keV/c]
  double  splat_ptot;     // splat total momentum [keV/x]
  double  splat_prad;     // splat prad (radial component of the momentum = √(py²+pz²) ) [keV/c]
  double  splat_plong;    // splat plong (longitudinal component of the momentum = px) [keV/c]
  double  splat_ke;       // splat kinetic energy (keV)
  double  splat_radius;   // splat radius (= √(y²+z²) ) [mm]
  long    run_id;         // univocal identifier for the run
  bool    ontarget;       // flag to mark that the particle hit the electron MCP
  bool    onmcpplane;     // flag to mark that the particle hit the electron MCP plane

  // These vectors contain the represented quantity at each time step in the particle flight
  // Only useful for the fully recorded data
  vector<double> *vec_tof     = new vector<double>;   // tof [us]
  vector<double> *vec_x       = new vector<double>;   // x position [mm]
  vector<double> *vec_y       = new vector<double>;   // y position [mm]
  vector<double> *vec_z       = new vector<double>;   // z position [mm]
  vector<double> *vec_azm     = new vector<double>;   // azimuth [deg]
  vector<double> *vec_elev    = new vector<double>;   // elevation [deg]
  vector<double> *vec_dircos  = new vector<double>;   // director cosine
  vector<double> *vec_vx      = new vector<double>;   // vx (X-component of the velocity) [mm/us]
  vector<double> *vec_vy      = new vector<double>;   // vy (Y-component of the velocity) [mm/us]
  vector<double> *vec_vz      = new vector<double>;   // vz (Z-component of the velocity) [mm/us]
  vector<double> *vec_vtot    = new vector<double>;   // total velocity [mm/us]
  vector<double> *vec_vrad    = new vector<double>;   // vrad (radial component of the velocity = √(vy²+vz²) ) [mm/us]
  vector<double> *vec_vlong   = new vector<double>;   // vlong (longitudinal component of the velocity = vx) [mm/us]
  vector<double> *vec_px      = new vector<double>;   // px (X-component of the momentum) [keV/c]
  vector<double> *vec_py      = new vector<double>;   // py (Y-component of the momentum) [keV/c]
  vector<double> *vec_pz      = new vector<double>;   // pz (Z-component of the momentum) [keV/c]
  vector<double> *vec_ptot    = new vector<double>;   // total momentum [keV/x]
  vector<double> *vec_prad    = new vector<double>;   // prad (radial component of the momentum = √(py²+pz²) ) [keV/c]
  vector<double> *vec_plong   = new vector<double>;   // plong (longitudinal component of the momentum = px) [keV/c]
  vector<double> *vec_ke      = new vector<double>;   // kinetic energy [keV]
  vector<double> *vec_radius  = new vector<double>;   // radius [mm]]
  vector<double> *vec_volt    = new vector<double>;   // voltage [V]
  vector<double> *vec_e       = new vector<double>;   // total electric field [V/mm]
  vector<double> *vec_ex      = new vector<double>;   // ex (X-component of the electric field) [V/mm]
  vector<double> *vec_ey      = new vector<double>;   // ey (Y-component of the electric field) [V/mm]
  vector<double> *vec_ez      = new vector<double>;   // ez (Z-component of the electric field) [V/mm]
  vector<double> *vec_b       = new vector<double>;   // total magnetic field [Gauss]
  vector<double> *vec_bx      = new vector<double>;   // ex (X-component of the magnetic field) [Gauss]
  vector<double> *vec_by      = new vector<double>;   // ey (Y-component of the magnetic field) [Gauss]
  vector<double> *vec_bz      = new vector<double>;   // ez (Z-component of the magnetic field) [Gauss]



  // init_elevation() [deg] - returns the initial elevation from 0 to 180
  double init_elevation() {return TMath::ACos(init_dircos) * TMath::RadToDeg();}


  // init_azimuth() [deg] - returns the initial azimuth from 0 to 360
  double init_azimuth()       {
    double val = TMath::ATan2(init_vz,init_vy) * TMath::RadToDeg();
    if(val < 0.){val += 360.;}  // comment for -pi to +pi
    return val;
  }

  // splat_azimuth() [deg] - returns the splat azimuth from 0 to 360
  double splat_azimuth()      {
    double val = TMath::ATan2(splat_z,splat_y) * TMath::RadToDeg();
    if(val < 0.){val += 360.;}  // comment for -pi to +pi
    return val;

  }

  // Reconstructed plong [keV/c]
  double calc_plong() {
    double output;
    if(mass > 1 && charge ==  1 && tof > TOF_CRITERION1)  {output = iso_f_plong1->Interpolate(tof,splat_radius);          }
    if(mass > 1 && charge ==  1 && tof < TOF_CRITERION1)  {output = (tof - TOF_INTERCEPTION_TO_MCP1)/SLOP_TOF_TO_MCP1;    }
    if(mass > 1 && charge ==  2 && tof > TOF_CRITERION2)  {output = iso_f_plong2->Interpolate(tof,splat_radius);          }
    if(mass > 1 && charge ==  2 && tof < TOF_CRITERION2)  {output = (tof - TOF_INTERCEPTION_TO_MCP2)/SLOP_TOF_TO_MCP2;    }
    if(mass < 1 && charge == -1)                          {output = iso_init_plong->Interpolate(tof,splat_radius);        }
    return output;
  }

  double calc_prad() {
    double output;
    if(mass > 1 && charge ==  1 && tof > TOF_CRITERION1)  {output = (splat_radius - TD_OPS_MCP_INTERCEPTION1)/SLOP_TD_OPS_MCP1; }
    if(mass > 1 && charge ==  1 && tof < TOF_CRITERION1)  {output = (splat_radius - TD_TO_MCP_INTERCEPTION1)/SLOP_TD_TO_MCP1  ; }
    if(mass > 1 && charge ==  2 && tof > TOF_CRITERION2)  {output = (splat_radius - TD_OPS_MCP_INTERCEPTION2)/SLOP_TD_OPS_MCP2; }
    if(mass > 1 && charge ==  2 && tof < TOF_CRITERION2)  {output = (splat_radius - TD_TO_MCP_INTERCEPTION2)/SLOP_TD_TO_MCP2  ; }
    if(mass < 1 && charge == -1)                          {output = iso_init_prad->Interpolate(tof,splat_radius);               }
    return output;
  }

  // Reconstructed azimuth [deg]
  double calc_azimuth()   {
    double output;
    if(mass > 1){output = splat_azimuth(); }
    if(mass < 1){output = iso_init_azimuth->Interpolate((omega*tof - TMath::Floor(omega*tof/(2*TMath::Pi())) * 2*TMath::Pi())/2.,splat_azimuth()); }
    return output;
  }

  // Reconstructed ptot [keV/c]
  double calc_ptot()      {return TMath::Sqrt( TMath::Power(calc_prad(),2) + TMath::Power(calc_plong(),2) );}

  // Reconstructed energy [keV]
  double calc_energy()    {return 1e3 * (TMath::Sqrt( TMath::Power(auger_mass,2) + TMath::Power(calc_ptot(), 2) ) - auger_mass ); }

  // Reconstructed azimuth (phi) [deg] [Ullrich formula]
  double calc_phi()       {
    double val = ( splat_azimuth()*TMath::DegToRad() + (omega*tof - TMath::Floor(omega*tof/(2*TMath::Pi())) * 2*TMath::Pi())/2. ) * TMath::RadToDeg();
    if(val < 0.){val +=360.;}
    return val;
  }

  // Uncertainty on plong [keV/c]
  double delta_plong()        {
    vector<double> *pr = new vector<double>;
    for(double i=-1.0;i<1.1;i+=0.1){
      for(double j=-1.0;j<1.1;j+=0.1){
        pr->push_back( iso_init_plong->Interpolate( tof+i*DEtof_MCP,splat_radius+j*DEsplat_radius_MCP ) );
      }
    }

    return (TMath::RMS(pr->begin(),pr->end()));
  }

  // Uncertainty on prad [keV/c]
  double delta_prad()        {
    vector<double> *pr = new vector<double>;
//     double sinomegat = TMath::Abs (TMath::Sin(omega*tof/2.) );
    double sinomegat = TMath::Sin(omega*tof/2.);
    double deltat = omega/2. * TMath::Cos(omega*tof/2.);

    for(double i=-1.0;i<1.1;i+=0.1){
      for(double j=-1.0;j<1.1;j+=0.1){
        pr->push_back( iso_init_prad->Interpolate( tof+i*DEtof_MCP, splat_radius+j*DEsplat_radius_MCP ) );
      }
    }

    return (TMath::RMS(pr->begin(),pr->end()));
  }

  // Uncertainty on ptot [keV/c]
  double delta_ptot()         {return (TMath::Sqrt( TMath::Power(calc_prad()*delta_prad(),2) + TMath::Power(calc_plong()*delta_plong(),2) ) / calc_ptot());}

  // Uncertainty on phi [deg]
  double delta_phi()          {return (TMath::Sqrt( TMath::Power(DEsplat_radius_MCP/splat_radius,2) + TMath::Power(0.5*omega*DEtof_MCP,2) )) * TMath::RadToDeg(); }

  // Uncertainty on splat_azimuth [deg]
  double delta_splat_azimuth(){return (TMath::Sqrt(  TMath::Power(splat_z*DEsplat_radius_MCP,2) + TMath::Power(splat_y*DEsplat_radius_MCP,2)     ) / (splat_radius*splat_radius)) * TMath::RadToDeg(); }

  // Uncertainty on prad [keV/c]
  double delta_calc_azimuth()        {
    vector<double> *pr = new vector<double>;
    double omegat = (omega*tof - TMath::Floor(omega*tof/(2*TMath::Pi())) * 2*TMath::Pi())/2.;

    for(double i=-1.0;i<1.1;i+=0.1){
      for(double j=-1.0;j<1.1;j+=0.1){
        pr->push_back( iso_init_azimuth->Interpolate( omegat+omega*i*DEtof_MCP/2, splat_azimuth()+j*delta_splat_azimuth() ) );
      }
    }

    return (TMath::RMS(pr->begin(),pr->end()));
  }

  // Uncertainty on the reconstructed energy [keV]
  double calc_delta_energy()  {return 1e3 * ((calc_ptot() * delta_ptot()) / (TMath::Sqrt( TMath::Power(auger_mass,2) + TMath::Power(calc_ptot(), 2) ) ) ); }

  // Calculated initial value of px, py and pz, from calculated azimuth
  double calc_px()  {return  calc_plong();}
  double calc_py()  {return (calc_prad() * TMath::Cos(calc_azimuth() * TMath::DegToRad()) );}
  double calc_pz()  {return (calc_prad() * TMath::Sin(calc_azimuth() * TMath::DegToRad()) );}

  bool cut_checker(bool onmcp = true, bool sprad = true, bool timef = true, bool nodec = true)
  {
    if(charge == -1){
      if(onmcp){ if(!(onmcpplane == true)){cout << "Auger not on MCP plane" << endl; return false;} }
      if(sprad){ if(!(splat_radius >= CX_SR_LOW && splat_radius <= CX_SR_HIGH)){cout << "Auger splat_radius = " << splat_radius << endl; return false;} }
      if(timef){ if(!(tof < CX_TOF)){cout << "TOF = " << tof << " [> CX_TOF = " << CX_TOF << "]" << endl; return false;} }
      if(nodec){ for(int nnodes=0; nnodes<sizeof(nodes)/sizeof(double); ++nnodes){ if(tof >= nodes[nnodes]-CX_NODE && tof <= nodes[nnodes]+CX_NODE){cout << "TOF = " << tof << " cut by CX_NODE" << endl; return false;} } }
    }
    if(charge >= 1){
      if(onmcp){ if(!(onmcpplane == true)){cout << "Ion not on MCP plane" << endl; return false;} }
      if(sprad){ if(!(splat_radius >= CX_SR_LOW && splat_radius <= CX_SR_HIGH)){cout << "ion splat_radius = " << splat_radius << endl; return false;} }
    }
    return true;
  }

  bool silent_cut_checker(bool onmcp = true, bool sprad = true, bool timef = true, bool nodec = true)
  {
    if(charge == -1){
      if(onmcp){ if(!(onmcpplane == true)){return false;} }
      if(sprad){ if(!(splat_radius >= CX_SR_LOW && splat_radius <= CX_SR_HIGH)){return false;} }
      if(timef){ if(!(tof < CX_TOF)){return false;} }
      if(nodec){ for(int nnodes=0; nnodes<sizeof(nodes)/sizeof(double); ++nnodes){ if(tof >= nodes[nnodes]-CX_NODE && tof <= nodes[nnodes]+CX_NODE){return false;} } }
    }
    if(charge >= 1){
      if(onmcp){ if(!(onmcpplane == true)){return false;} }
      if(sprad){ if(!(splat_radius >= CX_SR_LOW && splat_radius <= CX_SR_HIGH)){return false;} }
    }
    return true;
  }


  // QUANTITIES WITH ERRORS

  double generate_tof_werr()           { return ra->Gaus(tof,           DEtof_MCP         ); }
  double generate_splat_y_werr()       { return ra->Gaus(splat_y,       DEsplat_radius_MCP); }
  double generate_splat_z_werr()       { return ra->Gaus(splat_z,       DEsplat_radius_MCP); }
  double generate_splat_radius_werr()  { return ra->Gaus(splat_radius,  DEsplat_radius_MCP); }

  double tof_werr;
  double splat_y_werr;
  double splat_z_werr;
  double splat_radius_werr;

  // splat_azimuth() [deg] - returns the splat azimuth from 0 to 360
  double splat_azimuth_werr()      {
    double val = TMath::ATan2(splat_z_werr, splat_y_werr) * TMath::RadToDeg();
    if(val < 0.){val += 360.;}  // comment for -pi to +pi
    return val;

  }

  // Reconstructed plong [keV/c]
  double calc_plong_werr()     {return iso_init_plong->Interpolate(tof_werr, splat_radius_werr); }

  // Reconstructed plong [keV/c]
  double calc_prad_werr()      {return iso_init_prad->Interpolate(tof_werr, splat_radius_werr); }

  // Reconstructed ptot [keV/c]
  double calc_ptot_werr()      {return TMath::Sqrt( TMath::Power(calc_prad_werr(),2) + TMath::Power(calc_plong_werr(),2) );}

  // Reconstructed azimuth [deg]
  double calc_azimuth_werr()   {return iso_init_azimuth->Interpolate((omega*(tof_werr) - TMath::Floor(omega*(tof_werr)/(2*TMath::Pi())) * 2*TMath::Pi())/2.,splat_azimuth_werr()); }

  // Calculated initial value of px, py and pz, from calculated azimuth
  double calc_px_werr()  {return  calc_plong_werr();}
  double calc_py_werr()  {return (calc_prad_werr() * TMath::Cos(calc_azimuth_werr() * TMath::DegToRad()) );}
  double calc_pz_werr()  {return (calc_prad_werr() * TMath::Sin(calc_azimuth_werr() * TMath::DegToRad()) );}

};

// TChain that contains all the produced .root data files
TChain *so = new TChain("tends");

// Temporary graph and polyline used in the drawing functions
TPolyLine3D *gtrack_work;
TGraph   *gfront_work;
TGraph   *gfullfront_work;

// OpenGL viewer, used for the visualization routines
TGLViewer *v;

particle *xe  = new particle;
particle *xr  = new particle;
particle *nu  = new particle;
particle *e0  = new particle;
particle *e1  = new particle;

// Tree initialization routine
void tree_initializer(){

  so->Add("output_root/*.root");
  so->SetScanField(15);

  so->SetBranchAddress("xenon"  , &xe);
  so->SetBranchAddress("xray"   , &xr);
  so->SetBranchAddress("nu"     , &nu);
  so->SetBranchAddress("auger0" , &e0);
  so->SetBranchAddress("auger1" , &e1);

//   scatter->SetParameters(TMath::Sqrt(1.37055E-56),0.2713,0.035,0.08637,0.07182);
//   scatter->SetNpx(1000);
}

// Draws a 3D world marking the position of source, the two MCPs and endplanes
void DrawTrackWorld(){
  int full_x[8] = {   0,    0,    0,    0,    2200,   2200,   2200,   2200};
  int full_y[8] = {-300, -300,  300,  300,    -300,   -300,    300,    300};
  int full_z[8] = {-300,  300, -300,  300,    -300,    300,   -300,    300};
//   if(!(gROOT->FindObject("ctrack_full"))){TCanvas *ctrack_full = new TCanvas("ctrack_full","ctrack_full",800,800);}
  TCanvas *ctrack_full = new TCanvas("ctrack_full","ctrack_full",800,800);
  TGraph2D *gtrack_full = new TGraph2D(8,full_x,full_y,full_z);
  gtrack_full->SetTitle(so->GetTitle());
  gtrack_full->GetXaxis()->SetTitle("x (mm)");
  gtrack_full->GetYaxis()->SetTitle("y (mm)");
  gtrack_full->GetZaxis()->SetTitle("z (mm)");
  gtrack_full->GetXaxis()->SetTitleOffset(1.35);
  gtrack_full->GetYaxis()->SetTitleOffset(1.35);
  gtrack_full->GetZaxis()->SetTitleOffset(1.35);gtrack_full->Draw("P");

  // TRACK CIRCLES
  int slices = 180;
  double step = 2*TMath::Pi()/slices;
  TPolyLine3D *elec_plane3d = new TPolyLine3D(slices);
  TPolyLine3D *ion_plane3d = new TPolyLine3D(slices);
  TPolyLine3D *elec_mcp3d = new TPolyLine3D(slices);
  TPolyLine3D *ion_mcp3d  = new TPolyLine3D(slices);
  for(int i=0;i<slices+1;++i){
    elec_plane3d->SetPoint(i, ele_MCP_position+10., ele_inner_radius * TMath::Cos(i * step),  ele_inner_radius * TMath::Sin(i * step));
    ion_plane3d ->SetPoint(i, ion_MCP_position-10., src_outer_radius * TMath::Cos(i * step),  src_outer_radius * TMath::Sin(i * step));
    elec_mcp3d  ->SetPoint(i, ele_MCP_position,     MCP_radius * TMath::Cos(i * step),        MCP_radius * TMath::Sin(i * step));
    ion_mcp3d   ->SetPoint(i, ion_MCP_position,     MCP_radius * TMath::Cos(i * step),        MCP_radius * TMath::Sin(i * step));
  }
  elec_plane3d->SetLineColor(kBlack);
  elec_plane3d->SetLineWidth(5);
  ion_plane3d->SetLineColor(kBlack);
  ion_plane3d->SetLineWidth(5);
  elec_mcp3d->SetLineColor(kRed);
  elec_mcp3d->SetLineWidth(2);
  ion_mcp3d->SetLineColor(kRed);
  ion_mcp3d->SetLineWidth(2);
  elec_plane3d->Draw("SAME");
  ion_plane3d->Draw("SAME");
  elec_mcp3d->Draw("SAME");
  ion_mcp3d->Draw("SAME");

  // TRACK SOURCE
  TPolyMarker3D *gtrack_source = new TPolyMarker3D();
  gtrack_source->SetPoint(0,origin_x,0,0);
  gtrack_source->SetMarkerColor(kRed);
  gtrack_source->SetMarkerStyle(20);
  gtrack_source->Draw("SAME");
}

// Draws the splat plane on the electron side, as well as source, electron MCP and electron endplane
void DrawFrontWorld(){
  TCanvas *cfront_full = new TCanvas("cfront_full","cfront_full",800,800);
  cfront_full->SetGrid(1,1);
  TPad *padfront = new TPad("padfront","",0,0,1,1);
  TH1F *frafront = padfront->DrawFrame(-vessel_inner_radius-20, -vessel_inner_radius-20, vessel_inner_radius+20, vessel_inner_radius+20,"Front view (electrons);y (mm);z (mm)");
//   frafront->SetTitle("Front view (electrons)");
//   frafront->GetXaxis()->SetTitle("y (mm)");
//   frafront->GetYaxis()->SetTitle("z (mm)");
  frafront->GetXaxis()->SetTitleOffset(1.35);
  frafront->GetYaxis()->SetTitleOffset(1.35);

  // FRONT CIRCLES
  TEllipse *mcp     = new TEllipse(0, 0, MCP_radius);
  TEllipse *rings   = new TEllipse(0, 0, ele_inner_radius + (ele_outer_radius-ele_inner_radius)/2.);
  TEllipse *vessel  = new TEllipse(0, 0, vessel_inner_radius);
  mcp   ->SetLineColor(kRed);
  mcp   ->SetLineWidth(2);
  mcp   ->SetFillColor(0);
  mcp   ->SetFillStyle(0);
  rings ->SetLineColor(kGray);
  rings ->SetLineWidth(ele_outer_radius-ele_inner_radius);
  rings ->SetFillColor(0);
  rings ->SetFillStyle(0);
  vessel->SetLineColor(kBlack);
  vessel->SetLineWidth(5);
  vessel->SetFillColor(0);
  vessel->SetFillStyle(0);

  mcp   ->Draw("SAME");
  rings ->Draw("SAME");
  vessel->Draw("SAME");

  // FRONT SOURCE
  TPolyMarker *gfront_source = new TPolyMarker();
  gfront_source->SetPoint(0,0,0);
  gfront_source->SetMarkerColor(kRed);
  gfront_source->SetMarkerStyle(20);
  gfront_source->Draw("SAME");
}

// Draws the splat plane on the electron side, as well as source, electron MCP and electron endplane
void DrawFrontZoomWorld(){
  TCanvas *cfront_zoom = new TCanvas("cfront_zoom","cfront_zoom",800,800);
  cfront_zoom->SetGrid(1,1);
  TPad *padfront = new TPad("padfront","",0,0,1,1);
  TH1F *frafront = padfront->DrawFrame(-MCP_radius-10, -MCP_radius-10, MCP_radius+10, MCP_radius+10,"MCP view (electrons);y (mm);z (mm)");
//   frafront->SetTitle("Front view (electrons)");
//   frafront->GetXaxis()->SetTitle("y (mm)");
//   frafront->GetYaxis()->SetTitle("z (mm)");
  frafront->GetXaxis()->SetTitleOffset(1.35);
  frafront->GetYaxis()->SetTitleOffset(1.35);

  // FRONT CIRCLES
  TEllipse *mcp     = new TEllipse(0, 0, MCP_radius);
  mcp   ->SetLineColor(kRed);
  mcp   ->SetLineWidth(2);
  mcp   ->SetFillColor(0);
  mcp   ->SetFillStyle(0);
  mcp   ->Draw("SAME");

  // FRONT SOURCE
  TPolyMarker *gfront_source = new TPolyMarker();
  gfront_source->SetPoint(0,0,0);
  gfront_source->SetMarkerColor(kRed);
  gfront_source->SetMarkerStyle(20);
  gfront_source->Draw("SAME");
}

void ProduceFly2(){
  cout << "particles {\n";
  cout << "  -- Xe_plus\n";
  cout << "  standard_beam {\n";
  cout << "    mass = " << xe->mass/CF << ",\n";
  if(e1->init_ke < 1){
    cout << "    charge = 1,\n";
  }
  if(e1->init_ke >= 1){
    cout << "    charge = 2,\n";
  }
  cout << "    x = " << xe->init_x << " , y = " << xe->init_y << " , z = " << xe->init_z << ",\n";
  cout << "    direction = vector(" << xe->init_px/xe->init_ptot << " , " << xe->init_py/xe->init_ptot  << " , " << xe->init_pz/xe->init_ptot  << "),\n";
  cout << "    momentum = " << xe->init_ptot*1e-3 << ",\n";
  cout << "  },\n";
  cout << "  -- Xray\n";
  cout << "  standard_beam {\n";
  cout << "    mass = 1E-10,\n"; // manually changed to a very low number so that SIMION flies it
  cout << "    charge = 0,\n";
  cout << "    x = " << xr->init_x << " , y = " << xr->init_y << " , z = " << xr->init_z << ",\n";
  cout << "    direction = vector(" << xr->init_px/xr->init_ptot << " , " << xr->init_py/xr->init_ptot  << " , " << xr->init_pz/xr->init_ptot  << "),\n";
  cout << "    momentum = " << xr->init_ptot*1e-3 << ",\n";
  cout << "    color = 1,\n";
  cout << "  },\n";
  cout << "  -- Nu\n";
  cout << "  standard_beam {\n";
  cout << "    mass = 1E-10,\n"; // manually changed to a very low number so that SIMION flies it
  cout << "    charge = 0,\n";
  cout << "    x = " << nu->init_x << " , y = " << nu->init_y << " , z = " << nu->init_z << ",\n";
  cout << "    direction = vector(" << nu->init_px/nu->init_ptot << " , " << nu->init_py/nu->init_ptot  << " , " << nu->init_pz/nu->init_ptot  << "),\n";
  cout << "    momentum = " << nu->init_ptot*1e-3 << ",\n";
  cout << "    color = 2,\n";
  cout << "  },\n";
  // Augers
  cout << "  -- Auger\n";
  cout << "  standard_beam {\n";
  cout << "    mass = 0.00054857990946,\n";
  cout << "    charge = -1,\n";
  cout << "    x = " << e0->init_x << " , y = " << e0->init_y << " , z = " << e0->init_z << ",\n";
  cout << "    direction = vector(" << e0->init_px/e0->init_ptot << " , " << e0->init_py/e0->init_ptot  << " , " << e0->init_pz/e0->init_ptot  << "),\n";
  cout << "    momentum = " << e0->init_ptot*1e-3 << ",\n";
  cout << "  }";
  if(e1->init_ke > 1){
    cout << ",\n";
    cout << "  standard_beam {\n";
    cout << "    mass = 0.00054857990946,\n";
    cout << "    charge = -1,\n";
    cout << "    x = " << e1->init_x << " , y = " << e1->init_y << " , z = " << e1->init_z << ",\n";
    cout << "    direction = vector(" << e1->init_px/e1->init_ptot << " , " << e1->init_py/e1->init_ptot  << " , " << e1->init_pz/e1->init_ptot  << "),\n";
    cout << "    momentum = " << e1->init_ptot*1e-3 << ",\n";
    cout << "  }\n";
  }
  cout << "}\n";
}

void ProduceFly2File(){
  ofstream *outfly = new ofstream(); outfly->open( TString::Format("./fly2/output.fly2") );
  *outfly << "particles {\n";
  *outfly << "  -- Xe_plus\n";
  *outfly << "  standard_beam {\n";
  *outfly << "    mass = " << xe->mass << ",\n";
  if(e1->init_ke < 1){
    *outfly << "    charge = 1,\n";
  }
  if(e1->init_ke >= 1){
    *outfly << "    charge = 2,\n";
  }
  *outfly << "    x = " << xe->init_x << " , y = " << xe->init_y << " , z = " << xe->init_z << ",\n";
  *outfly << "    direction = vector(" << xe->init_px/xe->init_ptot << " , " << xe->init_py/xe->init_ptot  << " , " << xe->init_pz/xe->init_ptot  << "),\n";
  *outfly << "    momentum = " << xe->init_ptot*1e-3 << ",\n";
  *outfly << "  },\n";
  *outfly << "  -- Xray\n";
  *outfly << "  standard_beam {\n";
  *outfly << "    mass = 1E-10,\n"; // manually changed to a very low number so that SIMION flies it
  *outfly << "    charge = 0,\n";
  *outfly << "    x = " << xr->init_x << " , y = " << xr->init_y << " , z = " << xr->init_z << ",\n";
  *outfly << "    direction = vector(" << xr->init_px/xr->init_ptot << " , " << xr->init_py/xr->init_ptot  << " , " << xr->init_pz/xr->init_ptot  << "),\n";
  *outfly << "    momentum = " << xr->init_ptot*1e-3 << ",\n";
  *outfly << "    color = 1,\n";
  *outfly << "  },\n";
  *outfly << "  -- Nu\n";
  *outfly << "  standard_beam {\n";
  *outfly << "    mass = 1E-10,\n"; // manually changed to a very low number so that SIMION flies it
  *outfly << "    charge = 0,\n";
  *outfly << "    x = " << nu->init_x << " , y = " << nu->init_y << " , z = " << nu->init_z << ",\n";
  *outfly << "    direction = vector(" << nu->init_px/nu->init_ptot << " , " << nu->init_py/nu->init_ptot  << " , " << nu->init_pz/nu->init_ptot  << "),\n";
  *outfly << "    momentum = " << nu->init_ptot*1e-3 << ",\n";
  *outfly << "    color = 2,\n";
  *outfly << "  },\n";
  *outfly << "  -- Auger\n";
  *outfly << "  standard_beam {\n";
  *outfly << "    mass = 0.00054857990946,\n";
  *outfly << "    charge = -1,\n";
  *outfly << "    x = " << e0->init_x << " , y = " << e0->init_y << " , z = " << e0->init_z << ",\n";
  *outfly << "    direction = vector(" << e0->init_px/e0->init_ptot << " , " << e0->init_py/e0->init_ptot  << " , " << e0->init_pz/e0->init_ptot  << "),\n";
  *outfly << "    momentum = " << e0->init_ptot*1e-3 << ",\n";
  *outfly << "  }";
  if(e1->init_ke > 1){
    *outfly << ",\n";
    *outfly << "  standard_beam {\n";
    *outfly << "    mass = 0.00054857990946,\n";
    *outfly << "    charge = -1,\n";
    *outfly << "    x = " << e1->init_x << " , y = " << e1->init_y << " , z = " << e1->init_z << ",\n";
    *outfly << "    direction = vector(" << e1->init_px/e1->init_ptot << " , " << e1->init_py/e1->init_ptot  << " , " << e1->init_pz/e1->init_ptot  << "),\n";
    *outfly << "    momentum = " << e1->init_ptot*1e-3 << ",\n";
    *outfly << "  }\n";
  }
  *outfly << "}\n";
  outfly->close();
}