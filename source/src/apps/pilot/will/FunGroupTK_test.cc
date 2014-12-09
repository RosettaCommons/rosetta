// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilat/will/genmatch.cc
/// @brief ???

#include <protocols/kinmatch/FunGroupTK.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
//#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/kinematics/FoldTree.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/graph/Graph.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
// AUTO-REMOVED #include <core/pack/optimizeH.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <sstream>
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>

// AUTO-REMOVED #include <time.h>

#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/will_util.ihh>


using namespace protocols::kinmatch;

using core::Real;
using core::Size;
using core::pose::Pose;
using core::kinematics::Stub;
using core::conformation::Residue;
using core::conformation::ResidueOP;
using std::string;
using utility::vector1;
using numeric::min;
using core::import_pose::pose_from_pdb;
using basic::options::option;
using numeric::min;
using numeric::max;
using numeric::constants::d::pi;
using numeric::conversions::degrees;
using utility::io::ozstream;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;

typedef utility::vector1<core::Real> Reals;
typedef utility::vector1<core::Size> Sizes;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;
typedef utility::vector1<Vec> Vecs;

static thread_local basic::Tracer TR( "FunGroupTK_test" );

inline Real const sqr(Real const r) { return r*r; }

inline Vec xyz(Pose const & p, Size const & ia, Size const & ir) {
  return p.xyz(AtomID(ia,ir));
}


core::pack::rotamer_set::RotamerSetOP get_rotset(Pose & pose, Size icys) {
  core::pack::rotamer_set::RotamerSetOP rotset;
  core::scoring::ScoreFunction dummy_sfxn;
  dummy_sfxn( pose );
  core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task( pose );
  dummy_task->initialize_from_command_line();
  dummy_task->nonconst_residue_task( icys ).and_extrachi_cutoff(10);
  dummy_task->nonconst_residue_task( icys ).restrict_to_repacking();
  dummy_task->nonconst_residue_task( icys ).or_include_current( false ); //need to do this because the residue was built from internal coords and is probably crumpled up
  dummy_task->nonconst_residue_task( icys ).or_fix_his_tautomer( true ); //since we only want rotamers for the specified restype
  core::graph::GraphOP dummy_png = core::pack::create_packer_graph( pose, dummy_sfxn, dummy_task );
  core::pack::rotamer_set::RotamerSetFactory rsf;
  rotset = rsf.create_rotamer_set( pose.residue( icys ) );
  rotset->set_resid( icys );
  rotset->build_rotamers( pose, dummy_sfxn, *dummy_task, dummy_png );
  return rotset;
}


void run_hh(std::string fname) {
  using namespace std;
  using basic::options::option;
  using namespace basic::options::OptionKeys;
  using namespace core;
  using namespace chemical;
  using namespace id;
  using namespace pose;
  using namespace scoring;

  Pose gln,his;
  make_pose_from_sequence(gln,"Q","fa_standard",false);
  make_pose_from_sequence(his,"H","fa_standard",false);
  if(gln.residue(1).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(gln,1);
  if(gln.residue(1).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(gln,1);
  if(his.residue(1).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(his,1);
  if(his.residue(1).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(his,1);

  // setup stuff
  ResidueTypeSetCAP frs=ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
  ResidueTypeSetCAP crs=ChemicalManager::get_instance()->residue_type_set( "centroid" );
  ScoreFunctionOP sfstd=get_score_function();
  ScoreFunctionOP sfcen=ScoreFunctionFactory::create_score_function("score3");

  Pose pose;
  pose_from_pdb(pose,option[in::file::s]()[1]);
		for(Size i=1; i<=pose.n_residue(); ++i) {
				core::pose::replace_pose_residue_copying_existing_coordinates(pose,i,pose.residue(i).residue_type_set().name_map("ALA"));
		}
  FunGroupTK & krs(*(new   KinFunGroupTK(pose)));
  FunGroupTK & brs(*(new BruteFunGroupTK(pose)));
  TR << "start" << std::endl;

		vector1<core::pack::rotamer_set::RotamerSetOP> rots(pose.n_residue());
		for(Size i=1; i<=pose.n_residue(); ++i) {
				core::conformation::Residue rtmp(pose.residue(i));
				pose.replace_residue(i,his.residue(1),true);
				rots[i] = get_rotset(pose,i);
				pose.replace_residue(i,rtmp,false);
		}

  time_t ktime=0,btime=0;
  Size nkh=0,nbh=0;

  for(int i=0;i<1;i++) {
    nkh=0;
    nbh=0;
    // his / cys
    for(Size ir=1;ir<=pose.n_residue();++ir) {
      core::conformation::Residue rtmp(pose.residue(ir));
      pose.replace_residue(ir,his.residue(1),true);
      for(Size irot = 1; irot <= rots[ir]->num_rotamers(); ++irot) {
								//if( ir!=22 || irot != 5 ) continue;
        pose.set_chi(1,ir,rots[ir]->rotamer(irot)->chi(1));
        pose.set_chi(2,ir,rots[ir]->rotamer(irot)->chi(2));
        bool clash = false;
        for(Size ia = 6; ia <= pose.residue(ir).nheavyatoms(); ++ia) {
          if( !krs.ifc().clash_check(pose.residue(ir).xyz(ia),ir) ) { clash=true; break; }
        }
        if(clash) continue;
        pose.set_dof(core::id::DOF_ID(AtomID(pose.residue(ir).atom_index("HE2"),ir),D ),2.1);

        Vec ne = pose.residue(ir).xyz("NE2");
        Vec he = pose.residue(ir).xyz("HE2");
        Vec axs = (ne-he).normalized();
        Vec cen = ne + 2.1*(-axs);
        Vec ori = Vec(1.911136,0,irot);
        KRSQuery q(CEN_ANG,cen,axs,ori,0.5,0.175);
        vector1<core::conformation::ResidueOP> bhits,khits;

        time_t tmp = clock();
        for(int i=0;i<1;++i) krs.place_h( q, pose.residue(ir), khits );
        ktime += clock()-tmp;
        tmp = clock();
        //brs.place_c( q, pose.residue(ir), bhits );
        btime += clock()-tmp;
        nkh += khits.size();
        nbh += bhits.size();

        if( khits.size() != bhits.size() ) {
          for(Size ih=1;ih<=khits.size();++ih) {
            core::conformation::Residue rtmp2(pose.residue(khits[ih]->seqpos()));
            pose.replace_residue(khits[ih]->seqpos(),*khits[ih],true);
            pose.dump_pdb("KIN_"+lzs(ir,3)+"_"+lzs(irot,3)+"_"+lzs(ih,2)+"_"+lzs(khits[ih]->seqpos(),3)+".pdb");
            pose.replace_residue(khits[ih]->seqpos(),rtmp2,false);
          }
          for(Size ih=1;ih<=bhits.size();++ih) {
            core::conformation::Residue rtmp2(pose.residue(bhits[ih]->seqpos()));
            pose.replace_residue(bhits[ih]->seqpos(),*bhits[ih],false);
            pose.dump_pdb("BRT_"+lzs(ir,3)+"_"+lzs(irot,3)+"_"+lzs(ih,2)+"_"+lzs(bhits[ih]->seqpos(),3)+".pdb");
            pose.replace_residue(bhits[ih]->seqpos(),rtmp2,false);
          }
          /*utility_exit_with_message*/TR << ("kin vs. brute mismatch!!! rsd: "+str(ir)+" rot: "+str(irot)+" brute: "+str(bhits.size())+" kin: "+str(khits.size())) << std::endl;
        }

      }
      pose.replace_residue(ir,rtmp,false);
    }
  }
  TR << "BRT hits: " << nbh << " BRT time: " << btime << std::endl;
  TR << "KIN hits: " << nkh << " KIN time: " << ktime << std::endl;
  TR << "K/B ratio " << (Real)btime / (Real)ktime << std::endl;

}



void run_hd(std::string fname) {
  using namespace std;
  using basic::options::option;
  using namespace basic::options::OptionKeys;
  using namespace core;
  using namespace chemical;
  using namespace id;
  using namespace pose;
  using namespace scoring;

  Pose asp,his;
  make_pose_from_sequence(asp,"D","fa_standard",false);
  make_pose_from_sequence(his,"H","fa_standard",false);
  if(asp.residue(1).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(asp,1);
  if(asp.residue(1).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(asp,1);
  if(his.residue(1).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(his,1);
  if(his.residue(1).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(his,1);

  // setup stuff
  ResidueTypeSetCAP frs=ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
  //ResidueTypeSetCAP crs=ChemicalManager::get_instance()->residue_type_set( "centroid" );
  ScoreFunctionOP sfstd=get_score_function();
  //ScoreFunctionOP sfcen=ScoreFunctionFactory::create_score_function("score3");

  Pose pose;
  pose_from_pdb(pose,option[in::file::s]()[1]);
		for(Size i=1; i<=pose.n_residue(); ++i) {
				core::pose::replace_pose_residue_copying_existing_coordinates(pose,i,pose.residue(i).residue_type_set().name_map("ALA"));
		}
  FunGroupTK & krs(*(new   KinFunGroupTK(pose)));
  FunGroupTK & brs(*(new BruteFunGroupTK(pose)));
  TR << "start" << std::endl;

		vector1<core::pack::rotamer_set::RotamerSetOP> rots(pose.n_residue());
		for(Size i=1; i<=pose.n_residue(); ++i) {
				core::conformation::Residue rtmp(pose.residue(i));
				pose.replace_residue(i,his.residue(1),true);
				rots[i] = get_rotset(pose,i);
				pose.replace_residue(i,rtmp,false);
		}

  time_t ktime=0,btime=0;
  Size nkh=0,nbh=0;

  for(int i=0;i<1;i++) {
    nkh=0;
    nbh=0;
    // his / cys
    for(Size ir=1;ir<=pose.n_residue();++ir) {
      core::conformation::Residue rtmp(pose.residue(ir));
      pose.replace_residue(ir,his.residue(1),true);
      for(Size irot = 1; irot <= rots[ir]->num_rotamers(); ++irot) {
								//if( ir!=22 || irot != 5 ) continue;
        pose.set_chi(1,ir,rots[ir]->rotamer(irot)->chi(1));
        pose.set_chi(2,ir,rots[ir]->rotamer(irot)->chi(2));
        bool clash = false;
        for(Size ia = 6; ia <= pose.residue(ir).nheavyatoms(); ++ia) {
          if( !krs.ifc().clash_check(pose.residue(ir).xyz(ia),ir) ) { clash=true; break; }
        }
        if(clash) continue;
        //pose.set_dof(core::id::DOF_ID(AtomID(pose.residue(ir).atom_index("HE2"),ir),D ),2.1);

        Vec ne = pose.residue(ir).xyz("NE2");
        Vec he = pose.residue(ir).xyz("HE2");
        Vec axs = (ne-he).normalized();
        Vec cen = he;//ne + 1.0*(-axs);
        Vec ori = Vec(1.8,0,irot);
        KRSQuery q(CEN_AXS,cen,axs,ori,0.5,0.90);
        vector1<core::conformation::ResidueOP> bhits,khits;

        time_t tmp = clock();
        //for(int i=0;i<1;++i) krs.place_d( q, pose.residue(ir), khits );
        ktime += clock()-tmp;
        tmp = clock();
        brs.place_d( q, pose.residue(ir), bhits );
        btime += clock()-tmp;
        nkh += khits.size();
        nbh += bhits.size();

        if( khits.size() != bhits.size() ) {
          for(Size ih=1;ih<=khits.size();++ih) {
            core::conformation::Residue rtmp2(pose.residue(khits[ih]->seqpos()));
            pose.replace_residue(khits[ih]->seqpos(),*khits[ih],true);
            pose.dump_pdb("KIN_"+lzs(ir,3)+"_"+lzs(irot,3)+"_"+lzs(ih,2)+"_"+lzs(khits[ih]->seqpos(),3)+".pdb");
            pose.replace_residue(khits[ih]->seqpos(),rtmp2,false);
          }
          for(Size ih=1;ih<=bhits.size();++ih) {
            core::conformation::Residue rtmp2(pose.residue(bhits[ih]->seqpos()));
            pose.replace_residue(bhits[ih]->seqpos(),*bhits[ih],false);
            pose.dump_pdb("BRT_"+lzs(ir,3)+"_"+lzs(irot,3)+"_"+lzs(ih,2)+"_"+lzs(bhits[ih]->seqpos(),3)+".pdb");
            pose.replace_residue(bhits[ih]->seqpos(),rtmp2,false);
          }
          /*utility_exit_with_message*/TR << ("kin vs. brute mismatch!!! rsd: "+str(ir)+" rot: "+str(irot)+" brute: "+str(bhits.size())+" kin: "+str(khits.size())) << std::endl;
        }

      }
      pose.replace_residue(ir,rtmp,false);
    }
  }
  TR << "BRT hits: " << nbh << " BRT time: " << btime << std::endl;
  TR << "KIN hits: " << nkh << " KIN time: " << ktime << std::endl;
  TR << "K/B ratio " << (Real)btime / (Real)ktime << std::endl;

}


int main (int argc, char *argv[]) {

	try {



  devel::init(argc,argv);

  using basic::options::option;
  using namespace basic::options::OptionKeys;

  run_hd(option[in::file::s]()[1]);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}




