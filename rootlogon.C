# include "tchain_particles_nofunc.h"
// # include "tchain_particles.h"


double pi = TMath::Pi();

void rootlogon(){
  gROOT->ProcessLine( "gErrorIgnoreLevel = 3001;");
  tree_initializer();
  cout << "Objects initialized:" << endl;
  cout << "  TChain *so (contains the information on the SO and Auger tracks)" << endl << endl;
  cout << "Available functions:" << endl;
  cout << "  FindTrack(const char* selection) - Returns a pointer to vector<int> containing the track numbers corresponding to the selection parameters" << endl;
  cout << "  DrawTrackWorld() - Draws a 3D world marking the position of source, the two MCPs and endplanes" << endl;
  cout << "  DrawFrontWorld() - Draws the splat plane on the electron side, as well as source, electron MCP and electron endplane" << endl;
//   cout << "  DrawTrack(int nentry, Color_t mcolor = 1, int point_skip = 1) - Draws a new canvas with the 'nentry' 3D track" << endl;
//   cout << "  DrawFront(int nentry, Color_t mcolor = 1, Style_t mstyle = 7) - Draws a new canvas with the 'nentry' splat positions" << endl;
//   cout << "  DrawTrackSame(int nentry, Color_t mcolor = 1, int point_skip = 1) - Draws the 'nentry' 3D track in a pre-existing canvas (works best with DrawTrackWorld)" << endl;
//   cout << "  DrawFrontSame(int nentry, Color_t mcolor = 1, Style_t mstyle = 7) - Draws the 'nentry' splat positions in a pre-existing canvas (works best with DrawFrontWorld)" << endl;
  cout << "  DrawTrackExpression(const char* selection, Color_t mcolor = 1, int point_skip = 1) - Draws the 3D track corresponding to the 'selection' parameters in a pre-existing canvas" << endl;
  cout << "  DrawFrontExpression(const char* selection, Color_t mcolor = 1, Style_t mstyle = 7) - Draws the splat positions corresponding to the 'selection' parameters in a pre-existing canvas" << endl;
  cout << endl;

  cout << "File(s) in the TChain:" << endl;
  so->GetListOfFiles()->Print();
}
