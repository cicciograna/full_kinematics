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
};
