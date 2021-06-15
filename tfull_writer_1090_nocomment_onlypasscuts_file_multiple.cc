///////////////////////////////////////////////////////////////////////////////////
// tfull_writer_2000_nocomment.c (v5.0) - Author: Francesco Granato
// 2017-09-07
//
// The program reads the SIMION output, removing the header and adding a custom
// last line for exiting the loop: its arguments determine if SO or Auger electrons
// are selected, and their energy. Two trees are created, one for simple init and
// splat variables, the other for all the time steps. The SIMION output is read
// line by line, and its values stored in the associated variables; then the
// program writes the two trees to a file, "friends" the "tfull" tree to the
// init/splat one, and closes the ROOT file.
//
// 5.0  - Added univocal integer to identify same-batch runs
// 4.0  - Removed writing of last lines (the input .txt won't be reusable)
// 3.0  - Added voltage (volt) and E field (vec_e/ex/ey/ez) to the variables
// 2.0  - Changed input from energy and particle tipe to filename
// 1.2  - Converted dircos from radians to degrees
// 1.1  - Added parts to write dummy initial line and remove dummy final line
// 1.0  - Initial draft
///////////////////////////////////////////////////////////////////////////////////

# include "TROOT.h"
# include "TSystem.h"
# include <iostream>
# include <fstream>
# include <stdlib.h>
# include <vector>
# include "TString.h"
# include "TRandom.h"
# include "TFile.h"
# include "TTree.h"
# include "TMath.h"
# include "TObject.h"
# include "/home/tug26830/full_kinematics/particle_class.hh"
# include <gmp.h>     // hi-precision library
# include <mpfr.h>
# include "mpreal.h"

# include "no-single-double_scatter.hh" // defines if scattering is used
# include "hires_lores.hh"  // defines if we're using HiRes-LoAcc or LoRes-HiAcc tune

// # define ONLY_ION_PASSING_CUTS
// # define ONLY_AUGERS_PASSING_CUTS

using namespace std;
using mpfr::mpreal;

int tfull_writer(std::string nameofthefile, int startrun, int stoprun, long run_id){

  gSystem->cd("/home/tug26830/full_kinematics/");

  // General variables
  bool xe_nos_cuts_passed = false;
  bool xe_sng_cuts_passed = false;
  bool xe_dbl_cuts_passed = false;
  bool xe_all_cuts_passed = false;
  bool e0_cuts_passed = false;
  bool e1_exists = false;
  bool e1_cuts_passed = false;

  bool e1_passer;
  bool cuts_passed;


  double CX_SR_LOW(3.), CX_SR_HIGH(60.);
  # ifdef HIRES_LOACC_1236V
  double CX_TOF(0.425); // 1236 V
  double CX_NODE(0.010); // 1236 V
  # endif

  # ifdef LORES_HIACC_4100V
  double CX_TOF(0.232); // 41 V
  double CX_NODE(0.006); // 41 V
  # endif


  double nodes[] = {0.0441362, 0.0882725, 0.132409, 0.176545, 0.220681, 0.264817, 0.308954, 0.35309 , 0.397226, 0.441362, 0.485499, 0.529635, 0.573771, 0.617907, 0.662044, 0.70618 , 0.750316, 0.794452, 0.838589, 0.882725, 0.926861}; // 2000mm, E=0.1206, real B

 // To include <vector> in ROOT
  gROOT->ProcessLine("#include <vector>");

  gSystem->Load("/home/tug26830/full_kinematics/particle_class_hh.so");

  TString filename = nameofthefile.c_str();

  // SIMION OUTPUT READER
  TString ifstream_name = TString::Format("%s", filename.Data());
//   TString ifstream_name = TString::Format("%s", infile.Data());
//   cout << "Reading file " << ifstream_name.Data() << endl;
  ifstream data_file(ifstream_name.Data());
  if(!data_file){/*cout << "File " << ifstream_name.Data() << " not found. Aborting." << endl;*/ return 0;}
//   gROOT->ProcessLine(Form(".! echo \"$(tail -n +2 %s)\" > %s", ifstream_name.Data(), ifstream_name.Data())); // removes first line of .txt input file OBSOLETE
  gROOT->ProcessLine(Form(".! echo \"1.1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\" >> %s", ifstream_name.Data()));

  filename.Remove(0, filename.Last('/')+1 );
  filename.Remove(filename.Last('.'),filename.Length() );

  // Number of particles reader
  ifstream in_nparticles(TString::Format("/home/tug26830/work/fk_output/fly2_%06d-%06d/out_nparts.txt", startrun,stoprun));
  int nparticles;

  // VARIABLES INITIALIZATION

  double origin_x = 0.;
  double origin_y = 0.;
  double origin_z = 1090.;

  double MCP_radius = 60.;

  double ion_MCP_position = 159;

  double src_inner_radius = 70.;
  double src_outer_radius = src_inner_radius + 25.;

  double xray_radius = 400.;

  double ele_inner_radius = 171.;
  double ele_outer_radius = ele_inner_radius + 75.;
  double ele_MCP_position = origin_z + 2000.;

  int particle_type = 0;
  int nscat = 0;

  particle xenon ;
  # ifdef SINGLE_ION_SCATTER
  nscat = 1;
  particle xe_sng;
  # endif
  # ifdef DOUBLE_ION_SCATTER
  nscat = 2;
  particle xe_dbl;
  # endif
  particle xray  ;
  particle nu    ;
  particle auger0;
  particle auger1;



  // tends variables
  double ionn;
  double init_x;
  double init_y;
  double init_z;
  double init_azm;
  double init_elev;
  double init_dircos;
  double init_vx;
  double init_vy;
  double init_vz;
  double init_vtot;
  double init_vrad;
  double init_vlong;
  double init_px;
  double init_py;
  double init_pz;
  double init_ptot;
  double init_prad;
  double init_plong;
  double init_ke;
  double splat_tof;
  double splat_x;
  double splat_y;
  double splat_z;
  double splat_azm;
  double splat_elev;
  double splat_dircos;
  double splat_vx;
  double splat_vy;
  double splat_vz;
  double splat_vtot;
  double splat_vrad;
  double splat_vlong;
  double splat_px;
  double splat_py;
  double splat_pz;
  double splat_ptot;
  double splat_prad;
  double splat_plong;
  double splat_ke;
  double splat_radius;
  bool   ontarget;
  bool   onmcpplane;
//   TString *comment;

  // tfull variables
  double tof;
  double mass;
  double charge;
  double x;
  double y;
  double z;
  double azm;
  double elev;
  double dircos;
  double vx;
  double vy;
  double vz;
  double vtot;
  double vrad;
  double vlong;
  double px;
  double py;
  double pz;
  double ptot;
  double prad;
  double plong;
  double volt;
  double e;
  double ex;
  double ey;
  double ez;
  double ke;
  double b;
  double bx;
  double by;
  double bz;
  double radius;

  // Accessory variables
  string line;
  double ionn_check = 0;
  double c = 299792.458;  // mm/us
  double amukev = 931494.095; // AMU -> keV

  // HI-PRECISION VARIABLES
  // tends variables
  mpreal ginit_x;
  mpreal ginit_y;
  mpreal ginit_z;
  mpreal ginit_azm;
  mpreal ginit_elev;
  mpreal ginit_dircos;
  mpreal ginit_vx;
  mpreal ginit_vy;
  mpreal ginit_vz;
  mpreal ginit_vtot;
  mpreal ginit_vrad;
  mpreal ginit_vlong;
  mpreal ginit_px;
  mpreal ginit_py;
  mpreal ginit_pz;
  mpreal ginit_ptot;
  mpreal ginit_prad;
  mpreal ginit_plong;
  mpreal ginit_ke;
  mpreal gsplat_tof;
  mpreal gsplat_x;
  mpreal gsplat_y;
  mpreal gsplat_z;
  mpreal gsplat_azm;
  mpreal gsplat_elev;
  mpreal gsplat_dircos;
  mpreal gsplat_vx;
  mpreal gsplat_vy;
  mpreal gsplat_vz;
  mpreal gsplat_vtot;
  mpreal gsplat_vrad;
  mpreal gsplat_vlong;
  mpreal gsplat_px;
  mpreal gsplat_py;
  mpreal gsplat_pz;
  mpreal gsplat_ptot;
  mpreal gsplat_prad;
  mpreal gsplat_plong;
  mpreal gsplat_ke;
  mpreal gsplat_radius;
//   TString *comment;

  // tfull variables
  mpreal gtof;
  mpreal gmass;
  mpreal gcharge;
  mpreal gx;
  mpreal gy;
  mpreal gz;
  mpreal gazm;
  mpreal gelev;
  mpreal gdircos;
  mpreal gvx;
  mpreal gvy;
  mpreal gvz;
  mpreal gvtot;
  mpreal gvrad;
  mpreal gvlong;
  mpreal gpx;
  mpreal gpy;
  mpreal gpz;
  mpreal gptot;
  mpreal gprad;
  mpreal gplong;
  mpreal gvolt;
  mpreal ge;
  mpreal gex;
  mpreal gey;
  mpreal gez;
  mpreal gke;
  mpreal gb;
  mpreal gbx;
  mpreal gby;
  mpreal gbz;
  mpreal gradius;

  // Accessory variables
  mpreal gc = 299792.458;  // mm/us
  mpreal gamukev = 931494.095; // AMU -> keV

  vector<particle> *particle_list = new vector<particle>;
  particle_list->push_back(xenon );
  # ifdef SINGLE_ION_SCATTER
  particle_list->push_back(xe_sng);
  # endif
  # ifdef DOUBLE_ION_SCATTER
  particle_list->push_back(xe_dbl);
  # endif
  particle_list->push_back(xray  );
  particle_list->push_back(nu    );
  particle_list->push_back(auger0);
  particle_list->push_back(auger1);

  // FILE AND TREE CREATOR
  TString file_name = TString::Format("/home/tug26830/work/fk_output/output_root_%06d-%06d/%s_%06ld.root", startrun, stoprun, filename.Data(), run_id);
  //           cout << "Writing file " << file_name.Data() << endl;
  TFile *fout = new TFile(file_name.Data(),"RECREATE");

  // TREE CREATOR
  TString tends_name = TString::Format("tends_%s", filename.Data());
  TString tends_title = TString::Format("All electrons - Init and Splat");

  TString tfull_name = TString::Format("tfull_%s", filename.Data());
  TString tfull_title = TString::Format("All electrons - Every time step");

  TTree *tends_local = new TTree("tends", tends_title.Data());
  TTree *tfull_local = new TTree("tfull", tfull_title.Data());

  // BRANCH ALLOCATION
  // tends_local
  tends_local->Branch("xenon",          &particle_list->at(0)           );
  tends_local->Branch("xe_sng",         &particle_list->at(1)           );
  tends_local->Branch("xe_dbl",         &particle_list->at(2)           );
  tends_local->Branch("xray",           &particle_list->at(1+nscat)  );
  tends_local->Branch("nu",             &particle_list->at(2+nscat)  );
  tends_local->Branch("auger0",         &particle_list->at(3+nscat)  );
  tends_local->Branch("auger1",         &particle_list->at(4+nscat)  );

  getline(data_file, line); // reads first, useless line of kinem.txt input file

  // READING NUMBER OF EVENTS
  while(in_nparticles >> nparticles){
    particle_type = 0;
    for(int i=0;i<particle_list->size();++i){
      particle_list->at(i).run_id       = 0;
      particle_list->at(i).ionn         = 0;
      particle_list->at(i).mass         = 0;
      particle_list->at(i).charge       = 0;
      particle_list->at(i).init_x       = 0;
      particle_list->at(i).init_y       = 0;
      particle_list->at(i).init_z       = 0;
      particle_list->at(i).init_azm     = 0;
      particle_list->at(i).init_elev    = 0;
      particle_list->at(i).init_dircos  = 0;
      particle_list->at(i).init_vx      = 0;
      particle_list->at(i).init_vy      = 0;
      particle_list->at(i).init_vz      = 0;
      particle_list->at(i).init_vtot    = 0;
      particle_list->at(i).init_vrad    = 0;
      particle_list->at(i).init_vlong   = 0;
      particle_list->at(i).init_px      = 0;
      particle_list->at(i).init_py      = 0;
      particle_list->at(i).init_pz      = 0;
      particle_list->at(i).init_ptot    = 0;
      particle_list->at(i).init_prad    = 0;
      particle_list->at(i).init_plong   = 0;
      particle_list->at(i).init_ke      = 0;
      particle_list->at(i).tof          = 0;
      particle_list->at(i).splat_x      = 0;
      particle_list->at(i).splat_y      = 0;
      particle_list->at(i).splat_z      = 0;
      particle_list->at(i).splat_azm    = 0;
      particle_list->at(i).splat_elev   = 0;
      particle_list->at(i).splat_dircos = 0;
      particle_list->at(i).splat_vx     = 0;
      particle_list->at(i).splat_vy     = 0;
      particle_list->at(i).splat_vz     = 0;
      particle_list->at(i).splat_vtot   = 0;
      particle_list->at(i).splat_vrad   = 0;
      particle_list->at(i).splat_vlong  = 0;
      particle_list->at(i).splat_px     = 0;
      particle_list->at(i).splat_py     = 0;
      particle_list->at(i).splat_pz     = 0;
      particle_list->at(i).splat_ptot   = 0;
      particle_list->at(i).splat_prad   = 0;
      particle_list->at(i).splat_plong  = 0;
      particle_list->at(i).splat_ke     = 0;
      particle_list->at(i).splat_radius = 0;
    }

    // DA QUI...
    for(int ite = 0; ite < nparticles; ++ite){
      // READING OF SIMION OUTPUT
      ontarget    = false;
      onmcpplane  = false;

      // INIT DATA
      data_file >> ionn >> gtof >> gmass >> gcharge >> gx >> gy >> gz >> gazm >> gelev >> gvx >> gvy >> gvz >> gvolt >> ge >> gex >> gey >> gez >> gb >> gbx >> gby >> gbz >> gke;

      gvtot    = sqrt(gvx*gvx + gvy*gvy + gvz*gvz);
      gvrad    = sqrt(gvx*gvx + gvy*gvy);
      gvlong   = gvz;
      //     dircos  = vx / vtot * TMath::RadToDeg();
      gdircos  = gvz / gvtot;
      gplong   = ((gmass * gamukev * gvlong)/gc) / (sqrt(1. - (gvtot*gvtot)/(gc*gc)));
      gprad    = ((gmass * gamukev * gvrad)/gc)  / (sqrt(1. - (gvtot*gvtot)/(gc*gc)));
      gptot    = ((gmass * gamukev * gvtot)/gc)  / (sqrt(1. - (gvtot*gvtot)/(gc*gc)));
      gpx      = ((gmass * gamukev * gvx)/gc)    / (sqrt(1. - (gvtot*gvtot)/(gc*gc)));
      gpy      = ((gmass * gamukev * gvy)/gc)    / (sqrt(1. - (gvtot*gvtot)/(gc*gc)));
      gpz      = ((gmass * gamukev * gvz)/gc)    / (sqrt(1. - (gvtot*gvtot)/(gc*gc)));
      if(particle_type == 1+nscat || particle_type == 2+nscat){
        mpreal gxtheta = acos(gvz/gvtot);
        mpreal gxphi   = atan2(gvy,gvx);
        gptot    = gke/1000.;
        gpx      = gptot * sin(gxtheta) * cos(gxphi);
        gpy      = gptot * sin(gxtheta) * sin(gxphi);
        gpz      = gptot * cos(gxtheta);
        gplong   = gpz;
        gprad    = sqrt(gpx*gpx + gpy*gpy);
      }
      gradius  = sqrt(gx*gx + gy*gy);

      // Filling of init variables
      particle_list->at(particle_type).run_id      = run_id            ;
      particle_list->at(particle_type).ionn        = ionn              ;
      particle_list->at(particle_type).mass        = gmass  .toDouble();
      particle_list->at(particle_type).charge      = gcharge.toDouble();
      particle_list->at(particle_type).init_x      = gx     .toDouble();
      particle_list->at(particle_type).init_y      = gy     .toDouble();
      particle_list->at(particle_type).init_z      = gz     .toDouble();
      particle_list->at(particle_type).init_azm    = gazm   .toDouble();
      particle_list->at(particle_type).init_elev   = gelev  .toDouble();
      particle_list->at(particle_type).init_dircos = gdircos.toDouble();
      particle_list->at(particle_type).init_vx     = gvx    .toDouble();
      particle_list->at(particle_type).init_vy     = gvy    .toDouble();
      particle_list->at(particle_type).init_vz     = gvz    .toDouble();
      particle_list->at(particle_type).init_vtot   = gvtot  .toDouble();
      particle_list->at(particle_type).init_vrad   = gvrad  .toDouble();
      particle_list->at(particle_type).init_vlong  = gvlong .toDouble();
      particle_list->at(particle_type).init_px     = gpx    .toDouble();
      particle_list->at(particle_type).init_py     = gpy    .toDouble();
      particle_list->at(particle_type).init_pz     = gpz    .toDouble();
      particle_list->at(particle_type).init_ptot   = gptot  .toDouble();
      particle_list->at(particle_type).init_prad   = gprad  .toDouble();
      particle_list->at(particle_type).init_plong  = gplong .toDouble();
      particle_list->at(particle_type).init_ke     = gke    .toDouble();

      // SPLAT DATA
      data_file >> ionn >> gtof >> gmass >> gcharge >> gx >> gy >> gz >> gazm >> gelev >> gvx >> gvy >> gvz >> gvolt >> ge >> gex >> gey >> gez >> gb >> gbx >> gby >> gbz >> gke;
      gvtot    = sqrt(gvx*gvx + gvy*gvy + gvz*gvz);
      gvrad    = sqrt(gvx*gvx + gvy*gvy);
      gvlong   = gvz;
      //     dircos  = vx / vtot * TMath::RadToDeg();
      gdircos  = gvz / gvtot;
      gplong   = ((gmass * gamukev * gvlong)/gc) / (sqrt(1. - (gvtot*gvtot)/(gc*gc)));
      gprad    = ((gmass * gamukev * gvrad)/gc)  / (sqrt(1. - (gvtot*gvtot)/(gc*gc)));
      gptot    = ((gmass * gamukev * gvtot)/gc)  / (sqrt(1. - (gvtot*gvtot)/(gc*gc)));
      gpx      = ((gmass * gamukev * gvx)/gc)    / (sqrt(1. - (gvtot*gvtot)/(gc*gc)));
      gpy      = ((gmass * gamukev * gvy)/gc)    / (sqrt(1. - (gvtot*gvtot)/(gc*gc)));
      gpz      = ((gmass * gamukev * gvz)/gc)    / (sqrt(1. - (gvtot*gvtot)/(gc*gc)));
      if(particle_type == 1+nscat || particle_type == 2+nscat){
        mpreal gxtheta = acos(gvz/gvtot);
        mpreal gxphi   = atan2(gvy,gvx);
        gptot    = gke/1000.;
        gpx      = gptot * sin(gxtheta) * cos(gxphi);
        gpy      = gptot * sin(gxtheta) * sin(gxphi);
        gpz      = gptot * cos(gxtheta);
        gplong   = gpz;
        gprad    = sqrt(gpx*gpx + gpy*gpy);
      }
      gradius  = sqrt(gx*gx + gy*gy);

      particle_list->at(particle_type).tof           = gtof   .toDouble();
      particle_list->at(particle_type).splat_x       = gx     .toDouble();
      particle_list->at(particle_type).splat_y       = gy     .toDouble();
      particle_list->at(particle_type).splat_z       = gz     .toDouble();
      particle_list->at(particle_type).splat_azm     = gazm   .toDouble();
      particle_list->at(particle_type).splat_elev    = gelev  .toDouble();
      particle_list->at(particle_type).splat_dircos  = gdircos.toDouble();
      particle_list->at(particle_type).splat_vx      = gvx    .toDouble();
      particle_list->at(particle_type).splat_vy      = gvy    .toDouble();
      particle_list->at(particle_type).splat_vz      = gvz    .toDouble();
      particle_list->at(particle_type).splat_vtot    = gvtot  .toDouble();
      particle_list->at(particle_type).splat_vrad    = gvrad  .toDouble();
      particle_list->at(particle_type).splat_vlong   = gvlong .toDouble();
      particle_list->at(particle_type).splat_px      = gpx    .toDouble();
      particle_list->at(particle_type).splat_py      = gpy    .toDouble();
      particle_list->at(particle_type).splat_pz      = gpz    .toDouble();
      particle_list->at(particle_type).splat_ptot    = gptot  .toDouble();
      particle_list->at(particle_type).splat_prad    = gprad  .toDouble();
      particle_list->at(particle_type).splat_plong   = gplong .toDouble();
      particle_list->at(particle_type).splat_ke      = gke    .toDouble();
      particle_list->at(particle_type).splat_radius  = gradius.toDouble();


      // ontarget and onmcpplane checks
      if(particle_type == 0){ // only selects the Xenon ion (and the scattered versions)
        if(particle_list->at(particle_type).splat_z >= ion_MCP_position-0.002 &&
           particle_list->at(particle_type).splat_z <= ion_MCP_position+0.001 &&
           particle_list->at(particle_type).splat_radius <= 60.){
           particle_list->at(particle_type).ontarget = true;
        }
        if(particle_list->at(particle_type).splat_z >= ion_MCP_position-1 &&
           particle_list->at(particle_type).splat_z <= ion_MCP_position+3.5){
           particle_list->at(particle_type).onmcpplane = true;
        }
      }
      # ifdef SINGLE_ION_SCATTER

      if(particle_type == 1){ // only selects the Xenon ion (and the scattered versions)
        if(particle_list->at(particle_type).splat_z >= ion_MCP_position-0.002 &&
           particle_list->at(particle_type).splat_z <= ion_MCP_position+0.001 &&
           particle_list->at(particle_type).splat_radius <= 60.){
           particle_list->at(particle_type).ontarget = true;
        }
        if(particle_list->at(particle_type).splat_z >= ion_MCP_position-1 &&
           particle_list->at(particle_type).splat_z <= ion_MCP_position+3.5){
           particle_list->at(particle_type).onmcpplane = true;
        }
      }
      # endif
      # ifdef DOUBLE_ION_SCATTER
      if(particle_type == 2){ // only selects the Xenon ion (and the scattered versions)
        if(particle_list->at(particle_type).splat_z >= ion_MCP_position-0.002 &&
           particle_list->at(particle_type).splat_z <= ion_MCP_position+0.001 &&
           particle_list->at(particle_type).splat_radius <= 60.){
           particle_list->at(particle_type).ontarget = true;
        }
        if(particle_list->at(particle_type).splat_z >= ion_MCP_position-1 &&
           particle_list->at(particle_type).splat_z <= ion_MCP_position+3.5){
           particle_list->at(particle_type).onmcpplane = true;
        }
      }
      # endif
      //       if(particle_type == 1 ){ // only selects the xray - not necessary for now, the Xrays are generated as always on target
      //         if(particle_list->at(particle_type).splat_radius <= xray_radius+0.1 &&
      //            particle_list->at(particle_type).splat_radius >= xray_radius-0.1){
      //            particle_list->at(particle_type).ontarget = true;
      //         }
      //           if(particle_list->at(particle_type).splat_radius <= xray_radius+0.1 &&
      //              particle_list->at(particle_type).splat_radius >= xray_radius-0.1){
      //              particle_list->at(particle_type).onmcpplane = true;
      //           }
      //       }
      if( particle_type == 3+nscat || particle_type == 4+nscat ){ // only selects electrons
        if(particle_list->at(particle_type).splat_z >= ele_MCP_position-0.002 &&
           particle_list->at(particle_type).splat_z <= ele_MCP_position+0.001 &&
           particle_list->at(particle_type).splat_radius <= 60.){
           particle_list->at(particle_type).ontarget = true;
        }
        if(particle_list->at(particle_type).splat_z >= ele_MCP_position-1 &&
           particle_list->at(particle_type).splat_z <= ele_MCP_position+3.5){
           particle_list->at(particle_type).onmcpplane = true;
        }
      }

      ++particle_type;
    }

    # ifdef ONLY_ION_PASSING_CUTS
    // Writing only events for which Xenon ions pass cuts
    // Check on cuts for xe_noscatter
    if( (particle_list->at(0).splat_z >= ion_MCP_position-0.1 && particle_list->at(0).splat_z <= ion_MCP_position+0.1) &&
        (particle_list->at(0).splat_radius >= 0. && particle_list->at(0).splat_radius <= CX_SR_HIGH) ){xe_nos_cuts_passed = true;}

    # ifdef SINGLE_ION_SCATTER
    // Check on cuts for xe_single
    if( (particle_list->at(1).splat_z >= ion_MCP_position-0.1 && particle_list->at(1).splat_z <= ion_MCP_position+0.1) &&
        (particle_list->at(1).splat_radius >= 0. && particle_list->at(1).splat_radius <= CX_SR_HIGH) ){xe_sng_cuts_passed = true;}
    # endif

    # ifdef DOUBLE_ION_SCATTER
    // Check on cuts for xe_double
    if( (particle_list->at(2).splat_z >= ion_MCP_position-0.1 && particle_list->at(2).splat_z <= ion_MCP_position+0.1) &&
        (particle_list->at(2).splat_radius >= 0. && particle_list->at(2).splat_radius <= CX_SR_HIGH) ){xe_dbl_cuts_passed = true;}
    # endif

    xe_all_cuts_passed = xe_nos_cuts_passed;

    # if  defined SINGLE_ION_SCATTER &&  defined DOUBLE_ION_SCATTER
    xe_all_cuts_passed = xe_nos_cuts_passed + xe_sng_cuts_passed + xe_dbl_cuts_passed;
    # endif

    # if  defined SINGLE_ION_SCATTER && !defined DOUBLE_ION_SCATTER
    xe_all_cuts_passed = xe_nos_cuts_passed + xe_sng_cuts_passed;
    # endif

    # if !defined SINGLE_ION_SCATTER &&  defined DOUBLE_ION_SCATTER
    xe_all_cuts_passed = xe_nos_cuts_passed + xe_dbl_cuts_passed;
    # endif

    # endif

    # ifdef ONLY_AUGERS_PASSING_CUTS
    // Writing only events for which electrons pass cuts
    // Check on cuts for e0
    if(  (particle_list->at(3+nscat).splat_z >= ele_MCP_position-0.1 && particle_list->at(3+nscat).splat_z <= ele_MCP_position+0.1) &&
         (particle_list->at(3+nscat).splat_radius >= CX_SR_LOW && particle_list->at(3+nscat).splat_radius <= CX_SR_HIGH) &&
         (particle_list->at(3+nscat).tof < CX_TOF) && (
        !(particle_list->at(3+nscat).tof >= nodes[0]-0.01 && particle_list->at(3+nscat).tof <= nodes[0]+0.01) &&
        !(particle_list->at(3+nscat).tof >= nodes[1]-0.01 && particle_list->at(3+nscat).tof <= nodes[1]+0.01) &&
        !(particle_list->at(3+nscat).tof >= nodes[2]-0.01 && particle_list->at(3+nscat).tof <= nodes[2]+0.01) &&
        !(particle_list->at(3+nscat).tof >= nodes[3]-0.01 && particle_list->at(3+nscat).tof <= nodes[3]+0.01) &&
        !(particle_list->at(3+nscat).tof >= nodes[4]-0.01 && particle_list->at(3+nscat).tof <= nodes[4]+0.01) &&
        !(particle_list->at(3+nscat).tof >= nodes[5]-0.01 && particle_list->at(3+nscat).tof <= nodes[5]+0.01) &&
        !(particle_list->at(3+nscat).tof >= nodes[6]-0.01 && particle_list->at(3+nscat).tof <= nodes[6]+0.01) &&
        !(particle_list->at(3+nscat).tof >= nodes[7]-0.01 && particle_list->at(3+nscat).tof <= nodes[7]+0.01) &&
        !(particle_list->at(3+nscat).tof >= nodes[8]-0.01 && particle_list->at(3+nscat).tof <= nodes[8]+0.01) &&
        !(particle_list->at(3+nscat).tof >= nodes[9]-0.01 && particle_list->at(3+nscat).tof <= nodes[9]+0.01) ) ){e0_cuts_passed = true;}

    if(particle_list->at(4+nscat).init_ke > 1.){  // Check if e1 exists (if it does, it will have init_ke >1.
      e1_exists = true;
      // Check on cuts for e1
      if(  (particle_list->at(4+nscat).splat_z >= ele_MCP_position-0.1 && particle_list->at(4+nscat).splat_z <= ele_MCP_position+0.1) &&
           (particle_list->at(4+nscat).splat_radius >= CX_SR_LOW && particle_list->at(4+nscat).splat_radius <= CX_SR_HIGH) &&
           (particle_list->at(4+nscat).tof < CX_TOF) && (
          !(particle_list->at(4+nscat).tof >= nodes[0]-0.01 && particle_list->at(4+nscat).tof <= nodes[0]+0.01) &&
          !(particle_list->at(4+nscat).tof >= nodes[1]-0.01 && particle_list->at(4+nscat).tof <= nodes[1]+0.01) &&
          !(particle_list->at(4+nscat).tof >= nodes[2]-0.01 && particle_list->at(4+nscat).tof <= nodes[2]+0.01) &&
          !(particle_list->at(4+nscat).tof >= nodes[3]-0.01 && particle_list->at(4+nscat).tof <= nodes[3]+0.01) &&
          !(particle_list->at(4+nscat).tof >= nodes[4]-0.01 && particle_list->at(4+nscat).tof <= nodes[4]+0.01) &&
          !(particle_list->at(4+nscat).tof >= nodes[5]-0.01 && particle_list->at(4+nscat).tof <= nodes[5]+0.01) &&
          !(particle_list->at(4+nscat).tof >= nodes[6]-0.01 && particle_list->at(4+nscat).tof <= nodes[6]+0.01) &&
          !(particle_list->at(4+nscat).tof >= nodes[7]-0.01 && particle_list->at(4+nscat).tof <= nodes[7]+0.01) &&
          !(particle_list->at(4+nscat).tof >= nodes[8]-0.01 && particle_list->at(4+nscat).tof <= nodes[8]+0.01) &&
          !(particle_list->at(4+nscat).tof >= nodes[9]-0.01 && particle_list->at(4+nscat).tof <= nodes[9]+0.01) ) ){e1_cuts_passed = true;}
    }
    # endif

    # ifndef ONLY_ION_PASSING_CUTS
    xe_all_cuts_passed = true;
    # endif
    # ifndef ONLY_AUGERS_PASSING_CUTS
    e0_cuts_passed = true;
    e1_exists = false;
    # endif
    e1_passer = false;
    if((e1_exists == false) ||
      (e1_exists == true && e1_cuts_passed == true)){e1_passer = true;}

    cuts_passed = xe_all_cuts_passed * e0_cuts_passed * e1_passer;

    // Filling of trees

    if(cuts_passed){
      tends_local->Fill();
    }
  }

  // ...A QUI

//   gROOT->ProcessLine(Form(".! echo 'task goes here' | cat - %s > temp && mv temp %s", ifstream_name.Data(), ifstream_name.Data()));  // adda a random text ("task goes here") to the first line of .txt input file OBSOLETE
//   gROOT->ProcessLine(Form(".! head -n -1 %s > temp.txt ; mv temp.txt %s", ifstream_name.Data(), ifstream_name.Data()));

  tends_local->Write();
  fout->Close();
  return 0;
}

#ifndef __CINT__

int main(int argc, char **argv) {
//   std::cout << "\n==========> tfull_writer <=============" << std::endl;
  int id, start, stop;
  // if not enough arguments are provided
  if ( argc < 4 ) {
//     std::cout << "Usage:" << std::endl;
//     std::cout << "./tfull_writer filename(s)_to_convert startrun stoprun [randomizes run_id]" << std::endl;
//     std::cout << "./tfull_writer filename(s)_to_convert startrun stoprun run_id" << std::endl;
    return 0;
  }

  // if enough arguments are provided
  else {
    std::vector <std::string> *source_files = new std::vector <std::string>;
    start = atoi(argv[argc-3]); // checks third to last argument
    stop  = atoi(argv[argc-2]); // checks second to last argument
    id    = atoi(argv[argc-1]); // checks last argument

    // if last argument is not a run_id
    if( id == 0 ){
      id = 666;
//       std::cout << "Univocal file run_id = " << id << std::endl << std::endl;
      for (int i = 1; i < argc-2; ++i) { source_files->push_back( argv[i] ); }  // Remember argv[0] is the path to the program, we want from argv[1] onwards
      for( int j = 0; j < source_files->size(); ++j){
        tfull_writer( source_files->at(j), start, stop, id );
//         std::cout << std::endl;
      }
    }

    // if last argument is a run_id but only the run_id has been inserted
    else if( id != 0 && argc == 4){
//       std::cout << "Usage:" << std::endl;
//       std::cout << "./tfull_writer filename(s)_to_convert startrun stoprun [randomizes run_id]" << std::endl;
//       std::cout << "./tfull_writer filename(s)_to_convert startrun stoprun run_id" << std::endl;
      return 0;
    }

    // if last argument is a run_id
    else{
//       std::cout << "Univocal file run_id = " << id << std::endl << std::endl;
      for (int i = 1; i < argc-3; ++i) { source_files->push_back( argv[i] ); }  // Remember argv[0] is the path to the program, we want from argv[1] onwards
      for( int j = 0; j < source_files->size(); ++j){
        tfull_writer( source_files->at(j), start, stop, id );
//         std::cout << std::endl;
      }
    }
  }

//   std::cout << "==> Application finished." << std::endl;

  return 0;
}

#endif /* __CINT __ */
