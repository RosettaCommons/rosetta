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

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
#include <utility/io/izstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/mynamespaces.ihh>
#include <apps/pilot/will/will_util.ihh>


using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;
using core::pose::Pose;
using core::conformation::ResidueOP;

static thread_local basic::Tracer TR( "genI213" );
static core::io::silent::SilentFileData sfd;


inline Real sqr(Real const r) { return r*r; }
inline Real sigmoidish_neighbor( Real const & sqdist ) {
  if( sqdist > 9.*9. ) {
    return 0.0;
  } else if( sqdist < 6.*6. ) {
    return 1.0;
  } else {
    Real dist = sqrt( sqdist );
    return sqr(1.0  - sqr( (dist - 6.) / (9. - 6.) ) );
  }
}

vector1<Size> get_scanres(Pose const & pose) {
  vector1<Size> scanres;
  if(basic::options::option[basic::options::OptionKeys::willmatch::residues].user()) {
    TR << "input scanres!!!!!!" << std::endl;
    scanres = basic::options::option[basic::options::OptionKeys::willmatch::residues]();
  } else {
    for(Size i = 1; i <= pose.n_residue(); ++i) {
      if(!pose.residue(i).has("N" )) { continue; }
      if(!pose.residue(i).has("CA")) { continue; }
      if(!pose.residue(i).has("C" )) { continue; }
      if(!pose.residue(i).has("O" )) { continue; }
      if(!pose.residue(i).has("CB")) { continue; }
      if(pose.residue(i).name3()=="PRO") { continue; }
      scanres.push_back(i);
    }
  }
  return scanres;
}

void run() {

  using namespace basic::options::OptionKeys;

  core::chemical::ResidueTypeSetCAP  rs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
  Pose init,cys;
  core::import_pose::pose_from_pdb(init,*rs,option[in::file::s]()[1]);
  make_pose_from_sequence(cys,"C","fa_standard",false);
  remove_lower_terminus_type_from_pose_residue(cys,1);
  remove_upper_terminus_type_from_pose_residue(cys,1);
	add_variant_type_to_pose_residue(cys,"DISULF_PARTNER",1);

	for(Size ir = 1; ir <= init.n_residue(); ++ir) {
		init.replace_residue(ir,cys.residue(1),true);
		replace_pose_residue_copying_existing_coordinates(init,ir,init.residue(ir).residue_type_set().name_map("ALA"));
	}

	protocols::scoring::ImplicitFastClashCheck ifc(init,3.0);

  // Size nres = init.n_residue();
  // ScoreFunctionOP sf = core::scoring::get_score_function();
  ScoreFunctionOP sf = new core::scoring::symmetry::SymmetricScoreFunction(core::scoring::get_score_function());

  Pose pose = init;

  Mat R1 = rotation_matrix_degrees(Vec(0,0,1), 120.0);
  Mat R2 = rotation_matrix_degrees(Vec(0,0,1),-120.0);

  vector1<Size> scanres = get_scanres(pose);

  for(vector1<Size>::const_iterator iiter = scanres.begin(); iiter != scanres.end(); ++iiter) {
    Size irsd = *iiter;
    //if( pose.residue(irsd).name3()=="GLY" || pose.residue(irsd).name3()=="PRO" ) continue;
    ResidueOP rprev = pose.residue(irsd).clone();
    pose.replace_residue(irsd,cys.residue(1),true);
		for(Real ch1 = 0.0; ch1 < 360.0; ch1+=5.0) {
			pose.set_chi(1,irsd,ch1);
			for(Real ch2 = 0.0; ch2 < 360.0; ch2+=10.0) {
				pose.set_chi(2,irsd,ch2);
				Vec const  SG  = pose.residue(irsd).xyz( "SG");
				Vec const PSG  = pose.residue(irsd).xyz("PSG");
				Vec const PCB1 = pose.residue(irsd).xyz("PCB1");
				Vec const PCB2 = pose.residue(irsd).xyz("PCB2");
				if(!ifc.clash_check( SG,irsd)) continue;
				if(!ifc.clash_check(PSG,irsd)) continue;
				for(vector1<Size>::const_iterator jiter = scanres.begin(); jiter != scanres.end(); ++jiter) {
					Size jrsd = *jiter;
					Real dtgt = pose.residue(jrsd).xyz("CB").distance(PSG);
					if(dtgt < 10.0) continue;
					ResidueOP rprev2 = pose.residue(jrsd).clone();
					pose.replace_residue(jrsd,cys.residue(1),true);
					Pose tmp(pose);
					alignaxis(tmp,(PSG-PCB1).normalized(),(tmp.residue(jrsd).xyz("SG")-tmp.residue(jrsd).xyz("CB")).normalized());
					trans_pose(tmp, PSG - tmp.residue(jrsd).xyz("SG") );

					for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
						for(Size ia = 1; ia <= 5; ++ia) {
							if(!ifc.clash_check( tmp.xyz(AtomID(ia,ir)) )) goto cont1;
						}
					}
					goto done1; cont1:
					pose.replace_residue(jrsd,*rprev2,false);
					continue; done1:

					pose.dump_pdb("test1.pdb");
					tmp .dump_pdb("test2.pdb");
					utility_exit_with_message("arinestorae");
				}

			}
		}
    pose.replace_residue(irsd,*rprev,false);
  }

}

int main (int argc, char *argv[]) {

	try {

  devel::init(argc,argv);
  run();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

