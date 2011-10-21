// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

#define sym_Clo 5
#define sym_Chi 5
#define CONTACT_D2 20.25
#define CONTACT_TH 5
#define NSS 15872 // sphere_8192.dat.gz sphere_15872.dat.gz  sphere_32672.dat.gz sphere_78032.dat.gz
#define MIN_HELEX_RES 20
#define MAX_CYS_RES 3
#define MAX_NRES 200

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/init.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/Stub.hh>
#include <core/graph/Graph.hh>
#include <core/pose/annotated_sequence.hh>b
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/packer_neighbors.hh>

#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include "apps/pilot/will/will_util.hh"

typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;

using core::id::AtomID;
using basic::options::option;
using core::pose::Pose;
using core::Real;
using core::scoring::ScoreFunctionOP;
using core::Size;
using numeric::max;
using numeric::min;
using numeric::random::gaussian;
using numeric::random::uniform;
using numeric::rotation_matrix_degrees;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using ObjexxFCL::fmt::F;
using ObjexxFCL::fmt::I;
using ObjexxFCL::string_of;
using std::cerr;
using std::cout;
using std::string;
using std::pair;
using utility::io::izstream;
using utility::io::ozstream;
using utility::vector1;
using std::endl;
using core::import_pose::pose_from_pdb;
using core::kinematics::Stub;

static basic::Tracer tr("pentcb");
static core::io::silent::SilentFileData sfd;

vector1<pair<Real,Real> >
makerots(Pose const & p, Size ir, Pose const & rsd) {
  Pose pose(p);
  pose.replace_residue(ir,rsd.residue(1),true);
  core::scoring::ScoreFunction dummy_sfxn;
  dummy_sfxn.set_weight(core::scoring::fa_rep,0.44);
  dummy_sfxn( pose );
  core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task( pose );
  dummy_task->nonconst_residue_task(ir).restrict_to_repacking();
  dummy_task->nonconst_residue_task(ir).or_fix_his_tautomer(false);
  dummy_task->nonconst_residue_task(ir).and_extrachi_cutoff(1);
  dummy_task->nonconst_residue_task(ir).or_ex1(true);
  dummy_task->nonconst_residue_task(ir).or_ex2(true);
  // NO_EXtrA_CHI_SAMPLES = 0,      //0
  // EX_ONE_STDDEV,                 //1
  // EX_ONE_HALF_STEP_STDDEV,       //2
  // EX_TWO_FULL_STEP_STDDEVS,      //3
  // EX_TWO_HALF_STEP_STDDEVS,      //4
  // EX_FOUR_HALF_STEP_STDDEVS,     //5
  // EX_THREE_THIRD_STEP_STDDEVS,   //6
  // EX_SIX_QUARTER_STEP_STDDEVS,   //7
  dummy_task->nonconst_residue_task(ir).or_ex1_sample_level(core::pack::task::EX_TWO_FULL_STEP_STDDEVS);
  dummy_task->nonconst_residue_task(ir).or_ex2_sample_level(core::pack::task::EX_TWO_FULL_STEP_STDDEVS);
  core::graph::GraphOP dummy_png = core::pack::create_packer_graph( pose, dummy_sfxn, dummy_task );
  core::pack::rotamer_set::RotamerSetFactory rsf;
  core::pack::rotamer_set::RotamerSetOP rotset = rsf.create_rotamer_set( pose.residue(ir) );
  rotset->set_resid(ir);
  rotset->build_rotamers( pose, dummy_sfxn, *dummy_task, dummy_png );
  vector1<pair<Real,Real> > rots;
  for(Size irot = 1; irot <= rotset->num_rotamers(); ++irot) {
    Real r1 = rotset->rotamer(irot)->chi(1);
    Real r2 = rotset->rotamer(irot)->chi(2);
    bool seenit = false;
    for(vector1<pair<Real,Real> >::iterator i = rots.begin(); i != rots.end(); ++i) {
      if( fabs(r1-i->first) < 5.0 && fabs(r2-i->second) < 5.0 ) seenit=true;
    }
    if(!seenit) rots.push_back(pair<Real,Real>(r1,r2));
  }
  //tr << rotset->num_rotamers() << " " << rots.size() << std::endl;
  return rots;
}

bool clashcheck(Pose const & p, Vec v) {
  for(Size ir = 1; ir <= p.n_residue(); ++ir) {
    if(p.xyz(AtomID(1,ir)).distance_squared(v) < 9.0) return false;
    if(p.xyz(AtomID(2,ir)).distance_squared(v) < 9.0) return false;
    if(p.xyz(AtomID(3,ir)).distance_squared(v) < 9.0) return false;
    if(p.xyz(AtomID(4,ir)).distance_squared(v) < 9.0) return false;            
    if(p.residue(ir).name3()=="GLY") continue;
    if(p.xyz(AtomID(5,ir)).distance_squared(v) < 9.0) return false;                
  }
  return true;
}
bool clashcheckhalf(Pose const & p, Vec v) {
  for(Size ir = 1; ir <= p.n_residue(); ++ir) {
    if(p.xyz(AtomID(1,ir)).distance_squared(v) < 4.0) return false;
    if(p.xyz(AtomID(2,ir)).distance_squared(v) < 4.0) return false;
    if(p.xyz(AtomID(3,ir)).distance_squared(v) < 4.0) return false;
    if(p.xyz(AtomID(4,ir)).distance_squared(v) < 4.0) return false;            
    if(p.residue(ir).name3()=="GLY") continue;
    if(p.xyz(AtomID(5,ir)).distance_squared(v) < 4.0) return false;                
  }
  return true;
}

struct Hit {
  Vec cen,axs,ori;
  Size ir,jr,irot,jrot;
  Hit(Vec c, Vec a, Vec o, Size irs, Size jrs, Size irt, Size jrt) : cen(c),axs(a.normalized()),ori(o.normalized()),ir(irs),jr(jrs),irot(irt),jrot(jrt) {}
};

void dock(Pose const init, std::string const & fn) {
  Pose his,cys,ala;
  make_pose_from_sequence(ala,"A","fa_standard",false);
  remove_lower_terminus_type_from_pose_residue(ala,1);
  remove_upper_terminus_type_from_pose_residue(ala,1);
  make_pose_from_sequence(cys,"C[CYS_M]","fa_standard",false);
  remove_lower_terminus_type_from_pose_residue(cys,1);
  remove_upper_terminus_type_from_pose_residue(cys,1);
  make_pose_from_sequence(his,"H[HIS_DE]","fa_standard",false);
  remove_lower_terminus_type_from_pose_residue(his,1);
  remove_upper_terminus_type_from_pose_residue(his,1);
  Pose tmp(ala);

  tr << "get rotamers" << endl;
  vector1<vector1<pair<Real,Real> > > hrots(init.n_residue()),crots(init.n_residue());
  for(Size ir = 1; ir <= init.n_residue(); ++ir){
    hrots[ir] = makerots(init,ir,his);
    crots[ir] = makerots(init,ir,cys);
    //tr << ir << " " << crots[ir].size() << " " << hrots[ir].size() << endl;    
  }
  Size ihg = cys.residue(1).atom_index("ZN");
  Size ihd = his.residue(1).atom_index("ZN1");
  Size ihe = his.residue(1).atom_index("ZN2");
  Pose p(init);
  vector1<Size> nbcount(p.n_residue(),0);
  for(Size ir = 1; ir <= p.n_residue(); ++ir)
    for(Size jr = 1; jr <= p.n_residue(); ++jr)
      if(p.xyz(AtomID(2,ir)).distance_squared(p.xyz(AtomID(2,jr))) < 100.0) nbcount[ir]++;
  for(Size ir = 1; ir <= p.n_residue(); ++ir) p.replace_residue(ir,ala.residue(1),true);

  vector1<Hit> hits;
  for(Size ir = 1; ir <= p.n_residue(); ++ir) {
    if(p.residue(ir).name3()=="GLY"||p.residue(ir).name3()=="PRO") continue;
    if(nbcount[ir] < 10) continue;
    core::conformation::Residue itmp(p.residue(ir));
    p.replace_residue(ir,his.residue(1),true);
    for(Size irot = 1; irot <= hrots[ir].size(); ++irot) {
      p.set_chi(1,ir,hrots[ir][irot].first);
      p.set_chi(2,ir,hrots[ir][irot].second);
      Vec const xd = p.xyz(AtomID(ihd,ir));
      Vec const xe = p.xyz(AtomID(ihe,ir));      
      for(Size jr = 1; jr <= p.n_residue(); ++jr) {
        if(ir==jr) continue;
        if(nbcount[jr] < 10) continue;
        if(p.residue(jr).name3()=="GLY"||p.residue(jr).name3()=="PRO") continue;
        if(p.xyz(AtomID(2,jr)).distance_squared(xd) > 25.0) continue;
        core::conformation::Residue jtmp(p.residue(jr));
        p.replace_residue(jr,cys.residue(1),true);
        for(Size jrot = 1; jrot <= crots[jr].size(); ++jrot) {
          p.set_chi(1,jr,crots[jr][jrot].first);
          p.set_chi(2,jr,crots[jr][jrot].second);
          Vec const xc = p.xyz(AtomID(ihg,jr));
          bool clash = false;
          if( xc.distance_squared(xd) < 0.25 || xc.distance_squared(xe) < 0.25) {
            for(Size ia = 5; ia <= p.residue(ir).nheavyatoms(); ++ia) {
              if(ia==ihd||ia==ihe) continue;
              for(Size ja = 5; ja <= p.residue(jr).nheavyatoms(); ++ja) {
                if(ja==ihg) continue;
                if( p.xyz(AtomID(ia,ir)).distance_squared(p.xyz(AtomID(ja,jr))) < 9.0 ) clash=true;
              }
            }
            if(clash) continue;
            bool dhit = p.xyz(AtomID(ihg,jr)).distance_squared(xd) < 0.25;
            Vec xh = dhit ? xd : xe;
            Vec m = (xc+xh)/2.0;
            Vec nh = p.residue(ir).xyz(dhit?"ND1":"NE2");
            Vec sg = p.residue(jr).xyz("SG");
            if( fabs(109-angle_degrees(nh,m,sg)) > 10.0 ) continue;
            Vec a = (((m-nh).normalized()+(m-sg).normalized())/2.0).normalized();
            Vec L3 = rotation_matrix_degrees(a,90.0)*((2*m-sg)-m)+m;
            Vec L4 = rotation_matrix_degrees(a,90.0)*((2*m-nh)-m)+m;
            if(!clashcheck(p,L3)) continue;
            if(!clashcheck(p,L4)) continue;            
            tmp.replace_residue(1,p.residue(ir),false);
            tmp.set_xyz(AtomID(1,1),L3);
            tmp.dump_pdb("HIT_"+lzs(ir,3)+"_"+lzs(jr,3)+"_"+lzs(irot,3)+"_"+lzs(jrot,3)+"_H.pdb");
            tmp.replace_residue(1,p.residue(jr),false);
            tmp.set_xyz(AtomID(1,1),L4);
            tmp.dump_pdb("HIT_"+lzs(ir,3)+"_"+lzs(jr,3)+"_"+lzs(irot,3)+"_"+lzs(jrot,3)+"_C.pdb");
            tr << "HALFHIT "+lzs(ir,3)+"_"+lzs(jr,3)+"_"+lzs(irot,3)+"_"+lzs(jrot,3) << endl;
            hits.push_back(Hit(m,a,a.cross(L3-m),ir,jr,irot,jrot));
          }
        }
        p.replace_residue(jr,jtmp,false);
      }
    }
    p.replace_residue(ir,itmp,false);
  }

  Pose & pala(p);
  for(Size ih = 1; ih <= hits.size(); ++ih) {
    Hit & hi(hits[ih]);
    for(Size jh = ih+1; jh <= hits.size(); ++jh) {
      Hit & hj(hits[jh]);
      if(hi.ir==hj.ir||hi.jr==hj.ir||hi.ir==hj.jr||hi.jr==hj.jr) continue;
      if(hi.cen.distance_squared(hj.cen) < 64.0) continue;
      Vec c2cen = (hi.cen+hj.cen) / 2.0;
      if(!clashcheckhalf(pala,c2cen)) continue;
      Vec c2ori = (hi.cen-hj.cen).normalized();
      for(Size iaxs = 0; iaxs < 360; iaxs++) {
        Vec c2axs = rotation_matrix_degrees(c2ori,(Real)iaxs)*c2ori.cross(Vec(1,0,0));
        Mat c2rot = rotation_matrix_degrees(c2axs,180.0);
        // axes must be opposite ~180° 
        if(hi.axs.dot(c2rot*hj.axs) > -0.99384807753012208) continue;// cos(10°)
        // oris must be 90° rotated
        if(fabs(hi.ori.dot(c2rot*hj.ori)) > 0.17364817766693041) continue; //cos(80°)
        bool clash = false;
        // for(Size ir = 1; ir <= p.n_residue();++ir) {
        //   for(Size ia = 2; ia <= 2; ++ia) {
        //     if(!clashcheckhalf(pala,c2rot*(pala.xyz(AtomID(ia,ir))-c2cen)+c2cen)) clash=true;
        //   } if(clash) break;
        // }
        // if(clash) continue;
        tr << "FULLHIT!!!!!!" << endl;
        Pose p(pala);
        p.replace_residue(hi.ir,his.residue(1),true);
        p.replace_residue(hi.jr,cys.residue(1),true);
        p.replace_residue(hj.ir,his.residue(1),true);
        p.replace_residue(hj.jr,cys.residue(1),true);
        p.set_chi(1,hi.ir,hrots[hi.ir][hi.irot].first );
        p.set_chi(2,hi.ir,hrots[hi.ir][hi.irot].second);
        p.set_chi(1,hi.jr,crots[hi.jr][hi.jrot].first );
        p.set_chi(2,hi.jr,crots[hi.jr][hi.jrot].second);      
        p.set_chi(1,hj.ir,hrots[hj.ir][hj.irot].first );
        p.set_chi(2,hj.ir,hrots[hj.ir][hj.irot].second);      
        p.set_chi(1,hj.jr,crots[hj.jr][hj.jrot].first );
        p.set_chi(2,hj.jr,crots[hj.jr][hj.jrot].second);
        p.set_xyz(AtomID(ihg,hi.jr),hi.cen);
        p.set_xyz(AtomID(ihg,hj.jr),hj.cen);
        p.dump_pdb(utility::file_basename(fn)+"_"+lzs(ih,3)+lzs(jh,3)+lzs(iaxs,3)+"A.pdb");
        rot_pose(p,c2rot,c2cen);
        p.dump_pdb(utility::file_basename(fn)+"_"+lzs(ih,3)+lzs(jh,3)+lzs(iaxs,3)+"B.pdb");    
        //utility_exit_with_message("aosnrt");  
      }
    }
  }

}


int main(int argc, char *argv[]) {
  core::init(argc,argv);
  using namespace basic::options::OptionKeys;
  // loop over input files, do some checks, call dock
  for(Size ifn = 1; ifn <= option[in::file::s]().size(); ++ifn) {
    string fn = option[in::file::s]()[ifn];
    Pose pnat;
    tr << "checking " << fn << std::endl;
    core::import_pose::pose_from_pdb(pnat,fn);
    core::scoring::dssp::Dssp dssp(pnat);
    dssp.insert_ss_into_pose(pnat);
    Size cyscnt=0, nhelix=0;
    if( pnat.n_residue() > MAX_NRES ) goto cont1;
    for(Size ir = 2; ir <= pnat.n_residue()-1; ++ir) {
      if(pnat.secstruct(ir) == 'H') nhelix++;
      //if(!pnat.residue(ir).is_protein()) goto cont1;
      if(pnat.residue(ir).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
      if(pnat.residue(ir).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
      if(pnat.residue(ir).name3()=="CYS") { if(++cyscnt > MAX_CYS_RES) goto cont1; }
    } goto done1; cont1: tr << "skipping " << fn << std::endl; continue; done1:
    dock(pnat,fn);
  }
}

