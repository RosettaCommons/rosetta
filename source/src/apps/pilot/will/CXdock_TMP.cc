// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

#define CONTACT_D2 20.25 // = 4.5*4.5
#define CONTACT_TH 10
#define NSS 8192 // sphere_15872.dat.gz  sphere_32672.dat.gz sphere_78032.dat.gz sphere_8192.dat.gz
#define C_LO 8
#define C_HI 8
#define CONTACT_DIS 3.6

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <apps/pilot/will/will_util.ihh>

using core::id::AtomID;
using core::id::DOF_ID;
using core::kinematics::FoldTree;
using basic::options::option;
using core::pose::Pose;
using core::pose::PoseAP;
using core::pose::PoseCAP;
using core::pose::PoseCOP;
using core::pose::PoseOP;
using core::Real;
using core::scoring::ScoreFunctionOP;
using core::Size;
using numeric::max;
using numeric::min;
using numeric::random::gaussian;
using numeric::random::uniform;
using numeric::rotation_matrix_degrees;
using numeric::xyzVector;
using numeric::xyzMatrix;
using numeric::x_rotation_matrix_degrees;
using numeric::y_rotation_matrix_degrees;
using numeric::z_rotation_matrix_degrees;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::format::LJ;
using ObjexxFCL::format::SS;
using ObjexxFCL::string_of;
using std::cerr;
using std::cout;
using std::string;
using std::pair;
using utility::io::izstream;
using utility::io::ozstream;
using utility::pointer::owning_ptr;
using utility::pointer::access_ptr;
using utility::pointer::ReferenceCount;
using utility::vector1;
using std::endl;
using core::import_pose::pose_from_pdb;
using ObjexxFCL::lead_zero_string_of;
using ObjexxFCL::string_of;


typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;

using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;
using core::pose::Pose;
using core::conformation::ResidueOP;

static basic::Tracer TR("CXdock");
static core::io::silent::SilentFileData sfd;

// computes the fast slide-into-contact between point counds pa & pb
// separate cb lists... could be more efficient
// 2*BIN is the min. contact distance.. same for all points
Real sicfast(
  vector1<Vec> & pa , vector1<Vec> & pb ,
  vector1<Vec> & cba, vector1<Vec> & cbb,
  const double contact_dis,
  Size & cbcount // output arg!!!
){
  const Real BIN = contact_dis / 2.0;
  // get bounds for plane hashes
  double xmx1=-9e9,xmn1=9e9,ymx1=-9e9,ymn1=9e9,xmx=-9e9,xmn=9e9,ymx=-9e9,ymn=9e9;
  for(vector1<Vec>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
    xmx1 = max(xmx1,ia->x()); xmn1 = min(xmn1,ia->x());
    ymx1 = max(ymx1,ia->y()); ymn1 = min(ymn1,ia->y());
  }
  for(vector1<Vec>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
    xmx = max(xmx,ib->x()); xmn = min(xmn,ib->x());
    ymx = max(ymx,ib->y()); ymn = min(ymn,ib->y());
  }
  xmx = min(xmx,xmx1); xmn = max(xmn,xmn1);
  ymx = min(ymx,ymx1); ymn = max(ymn,ymn1);
  int xlb = (int)floor(xmn/BIN)-2; int xub = (int)ceil(xmx/BIN)+2; // one extra on each side for correctness,
  int ylb = (int)floor(ymn/BIN)-2; int yub = (int)ceil(ymx/BIN)+2; // and one extra for outside atoms

  // insert points into hashes
  int const xsize = xub-xlb+1;
  int const ysize = yub-ylb+1;
  ObjexxFCL::FArray2D<Vec> ha(xsize,ysize,Vec(0,0,-9e9)),hb(xsize,ysize,Vec(0,0,9e9));
  for(vector1<Vec>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
    int const ix = (int)ceil(ia->x()/BIN)-xlb;
    int const iy = (int)ceil(ia->y()/BIN)-ylb;
    if( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
    if( ha(ix,iy).z() < ia->z() ) ha(ix,iy) = *ia;
  }
  for(vector1<Vec>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
    int const ix = (int)ceil(ib->x()/BIN)-xlb;
    int const iy = (int)ceil(ib->y()/BIN)-ylb;
    if( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
    if( hb(ix,iy).z() > ib->z() ) hb(ix,iy) = *ib;
  }
  // check hashes for min dis
  double mindis = 9e9;
  for(int i = 1; i <= xsize; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
    for(int j = 1; j <= ysize; ++j) {
      for(int k = -2; k <= 2; ++k) {
        if(i+k < 1 || i+k > xsize) continue;
        for(int l = -2; l <= 2; ++l) {
          if(j+l < 1 || j+l > ysize) continue;
          double const xa = ha(i  ,j  ).x();
          double const ya = ha(i  ,j  ).y();
          double const xb = hb(i+k,j+l).x();
          double const yb = hb(i+k,j+l).y();
          double const d2 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb);
          if( d2 < BIN*BIN*4.0 ) {
            double dz = hb(i+k,j+l).z() - ha(i,j).z() - sqrt(BIN*BIN*4.0-d2);
            mindis = min(dz,mindis);
          }
        }
      }
    }
  }
  // cb contact count... could be another objective func is better
  cbcount = 0;
  for(vector1<Vec>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
    for(vector1<Vec>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
      if( ib->distance_squared( (*ia) + (mindis*Vec(0,0,1)) ) < CONTACT_D2 ) {
        cbcount++;
      }
    }
  }
  return cbcount;
}

// basic info for a "hit"
struct Hit {
  int iss,irt,cbc,sym;
  core::kinematics::Stub s1,s2;
  Hit(int is, int ir, int cb, int sm) : iss(is),irt(ir),cbc(cb),sym(sm) {}
};
bool cmp(Hit i,Hit j) { return i.cbc > j.cbc; }

// finely sample every RB rotation of input structure
// dock against rotated-self, rotation based on sym
void dock(Pose const init, std::string const & fn, vector1<xyzVector<double> > const & ssamp) {

  // rotation angle samples.. in theory, should NOT be evenly spaced but this is simpler
  // probably "wastes" 1/3 of calculation with extra resolution at small rotations (somewhat
  // redundant with close axes of rot. from sphere points)
  vector1<double> asamp; for(Real i = 0; i < 180; ++i) asamp.push_back(i);

  // cache initial BB and CB coords (if no CB, use CA)
  vector1<Vec> bb0,cb0;
  for(int ir = 1; ir <= init.n_residue(); ++ir) {
    if(!init.residue(ir).is_protein()) continue;
    for(int ia = 1; ia <= ((init.residue(ir).has("CB"))?5:4); ++ia) bb0.push_back(init.xyz(AtomID(ia,ir)));
    // only count helix CBs!!!!!!!!!!
    if(init.secstruct(ir)=='H') cb0.push_back(init.xyz(AtomID(((init.residue(ir).has("CB"))?5:4),ir)));
  }

  // searches for CX, X in 2...NSYM
  TR << "searching symmetries C"<<C_LO<<" to C"<<C_HI<<"."<<std::endl;
  vector1<Mat> Rsym(C_HI,Mat::identity());
  for(int ic = C_LO; ic <= C_HI; ic++) Rsym[ic] = rotation_matrix_degrees(Vec(1,0,0),360.0/Real(ic));
  vector1<vector1<Hit> > hits(C_HI); // bit wasteful... whatever
  // loop over axes
  for(int iss = 1; iss <= ssamp.size(); ++iss) { // loop over sphere samples "ssamp"
    // status message
    if(iss%100==0) { TR << iss << " of " << NSS << " Nhits:"; for(int ic = C_LO; ic <= C_HI; ic++) TR << " " << hits[ic].size(); TR <<std::endl; }
    // axis of rot from sphere surf samples
    Vec axs = ssamp[iss];
    // loop axis of rot from sphere surf samples
    for(int irt = 1; irt <= asamp.size(); ++irt) {
      // set up point cloud A
      Mat const R = rotation_matrix_degrees( axs, asamp[irt] );
      vector1<Vec> bb1 = bb0;
      vector1<Vec> cb1 = cb0;
      for(vector1<Vec>::iterator i = bb1.begin(); i != bb1.end(); ++i) *i = R*(*i);
      for(vector1<Vec>::iterator i = cb1.begin(); i != cb1.end(); ++i) *i = R*(*i);
      for(int ic = C_LO; ic <= C_HI; ic++) {
        // set up point cloud B
        vector1<Vec> bb2 = bb1;
        vector1<Vec> cb2 = cb1;
        for(vector1<Vec>::iterator i = bb2.begin(); i != bb2.end(); ++i) *i = Rsym[ic]*(*i);
        for(vector1<Vec>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = Rsym[ic]*(*i);
        // fast slide-into-contact along Z
        Size cbc = 0;
        Real t = sicfast(bb1,bb2,cb1,cb2,CONTACT_DIS,cbc);
        // remember if good
        if(cbc >= CONTACT_TH) {
          Hit h(iss,irt,cbc,ic);
          h.s1.from_four_points(bb1[1],bb1[1],bb1[2],bb1[3]);
          h.s2.from_four_points(bb2[1],bb2[1],bb2[2],bb2[3]);
          h.s1.v += t*Vec(0,0,1);
          hits[ic].push_back(h);
        }
      }
    }
  }

  // report top 10 hits by CB count
  for(int ic = C_LO; ic <= C_HI; ic++) {
    std::sort(hits[ic].begin(),hits[ic].end(),cmp);
    for(int i = 1; i <= min((Size)10,hits[ic].size()); ++i) {
      Hit & h(hits[ic][i]);
      cout << "RESULT " << fn << " " << h.sym << " " << h.iss << " " << NSS  << " " << h.irt << " " << h.cbc << endl;
    }
  }

  // for testing, dump first hit
  Pose p(init),q(init);
  core::kinematics::Stub s(init.xyz(AtomID(1,1)),init.xyz(AtomID(2,1)),init.xyz(AtomID(3,1)));
  xform_pose_rev(p,s); xform_pose(p,hits[C_LO][1].s1);
  xform_pose_rev(q,s); xform_pose(q,hits[C_LO][1].s2);
  p.dump_pdb(utility::file_basename(fn)+"_hit1_A.pdb");
  q.dump_pdb(utility::file_basename(fn)+"_hit1_B.pdb");

}

int main(int argc, char *argv[]) {

	try {

  devel::init(argc,argv);
  using namespace basic::options::OptionKeys;

  vector1<xyzVector<double> > ssamp(NSS);
  {
    izstream is;
    basic::database::open(is,"geometry/sphere_"+str(NSS)+".dat");
    for(int i = 1; i <= NSS; ++i) {
      double x,y,z;
      is >> x >> y >> z;
      ssamp[i] = xyzVector<double>(x,y,z);
    }
    is.close();
  }

  for(Size ifn = 1; ifn <= option[in::file::s]().size(); ++ifn) {
    string fn = option[in::file::s]()[ifn];
    Pose pnat;
    TR << "searching " << fn << std::endl;
    core::import_pose::pose_from_pdb(pnat,fn);
    trans_pose(pnat,-center_of_geom(pnat,1,pnat.n_residue()));
    core::scoring::dssp::Dssp dssp(pnat);
    dssp.insert_ss_into_pose(pnat);
    if( pnat.n_residue() > 200 ) continue;
    Size cyscnt=0, nhelix=0;
    for(Size ir = 2; ir <= pnat.n_residue()-1; ++ir) {
      if(pnat.secstruct(ir) == 'H') nhelix++;
      //if(!pnat.residue(ir).is_protein()) goto cont1;
      if(pnat.residue(ir).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
      if(pnat.residue(ir).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
      if(pnat.residue(ir).name3()=="CYS") { if(++cyscnt > 3) goto cont1; }
    } goto done1; cont1: TR << "skipping " << fn << std::endl; continue; done1:
    if( nhelix < 20 ) continue;
    Pose pala(pnat);
    dock(pala,fn,ssamp);
  }

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
