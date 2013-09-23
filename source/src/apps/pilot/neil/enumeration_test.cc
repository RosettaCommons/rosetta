// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


#include <basic/database/open.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/options/keys/holes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/packing/HolesParams.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <fstream>
#include <iostream>
#include <numeric/random/random.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <protocols/docking/util.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/moves/Mover.hh>
#include <sstream>
#include <string>
#include <utility/io/izstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
// Added 100922
// AUTO-REMOVED #include <basic/prof.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>
// AUTO-REMOVED #include <numeric/model_quality/rms.hh>
// AUTO-REMOVED #include <core/scoring/dssp/Dssp.hh>
// AUTO-REMOVED #include <protocols/simple_moves/MinMover.hh>
// AUTO-REMOVED #include <pstream.h>
// AUTO-REMOVED #include <time.h>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
//Auto Headers


static basic::Tracer TR("enumeration_test");



using std::string;
using ObjexxFCL::string_of;
using namespace core;
using utility::vector1;
using core::id::AtomID;
using namespace core::scoring::packing;
using namespace ObjexxFCL::format;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;

int
count_int_CBs_clashes (pose::Pose const &pose, Size i, Size j, Real contact_dist, Real clash_dist) {

  Real const contact_dist_sq = contact_dist * contact_dist;
  Real const clash_dist_sq = clash_dist * clash_dist;

  int int_cb_count = 0;
  Size itype = 5;
  Size jtype = 5;

  conformation::symmetry::SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
  Size nsub = sym_info->subunits();
  
  if (i > nsub || j > nsub) {
    utility_exit_with_message("There are not that many subunits!");
  }

  for (Size ir=1; ir<=pose.n_residue(); ir++) { // Iterate over all residues in the subunit i
    Size isub = sym_info->subunit_index(ir);
    if (isub != i) continue;
    if (chemical::name_from_aa(pose.aa(ir)) == "GLY") { // Use CAs instead of CBs for GLYs
      itype = 4;
    } else {
      itype = 5;
    }
    for (Size jr=1; jr<=pose.n_residue(); jr++) { // Then iterate over all residues in subunit j
      Size jsub = sym_info->subunit_index(jr);
      if (jsub != j) continue;
      if (chemical::name_from_aa(pose.aa(jr)) == "GLY") {
        jtype = 4;
      } else {
        jtype = 5;
      }
      if (pose.xyz(AtomID(itype,ir)).distance_squared(pose.xyz(AtomID(jtype,jr))) <= contact_dist_sq) { // Check if the CBs/CAs are within contact dist
        int_cb_count++;
        for (Size ia=1; ia<=itype; ia++) {
          for (Size ja=1; ja<=jtype; ja++) {
//            if (pose.xyz(AtomID(ia,ir)).distance_squared(pose.xyz(AtomID(ja,jr))) <= clash_dist_sq) return -1;
            if (pose.xyz(AtomID(ia,ir)).distance_squared(pose.xyz(AtomID(ja,jr))) <= clash_dist_sq) { // Test for clashes
              if ( (((ia == 1) && (ja == 4)) || ((ia == 4) && (ja == 1))) && (pose.xyz(AtomID(ia,ir)).distance_squared(pose.xyz(AtomID(ja,jr))) >= 6.76) ) { // But don't count bb-bb h-bonds as clashes
                continue;
              } else {
                return -1;
              }
            }
/*//////////
            
            if (pose.xyz(AtomID(ia,ir)).distance_squared(pose.xyz(AtomID(ja,jr))) <= clash_dist_sq) {
              TR<< "CLASH: Chain " << i << ", Atom " << ia << ", Resi " << ir << ", " << chemical::name_from_aa(pose.aa(ir)) << " is sqrt(" << pose.xyz(AtomID(ia,ir)).distance_squared(pose.xyz(AtomID(ja,jr))) << ") from Chain " << j << ", Atom " << ja << ", Resi " << jr << ", " << chemical::name_from_aa(pose.aa(jr)) << ", jtype=" << jtype << std::endl;
              return -1;
            }
*//////////
          }
        }
      }
    }
  }

  return int_cb_count;

}

void
move_pose(pose::Pose & pose, Size sym_jump, Vec dt, Real dr = 0.0) {

  kinematics::Jump j = pose.jump(sym_jump);

  Vec t = j.get_translation();
  Mat r = j.get_rotation();

  j.set_translation(t+dt);
  j.set_rotation(numeric::x_rotation_matrix_degrees(dr)*r);
  pose.set_jump(sym_jump,j);

}

int
main (int argc, char *argv[])
{
	try {
  using namespace core;
	using namespace basic;
	using namespace options;
  using namespace OptionKeys;
  using namespace pose;
  using namespace core::conformation::symmetry;

  devel::init(argc,argv);
  chemical::ResidueTypeSetCAP resi_set = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");

	// Get some of the command line options for docking
	Real dt = 1;
	Real offset = 1;
	if (options::option[in::olig_search::neg_r]() == 1) {
		dt = -1;
		offset = -1;
	}
	//TR << string_of(dt) << "\t" << string_of(offset) << std::endl;
	Size nsub_bb = options::option[in::olig_search::num_subs_building_block]();
	Size nsub_total = options::option[in::olig_search::num_subs_total]();

  // Read in pose
  Pose pose;
  import_pose::pose_from_pdb(pose, option[in::file::s]()[1], resi_set);

  // If you are docking, change all non-glycine residues to ala to cut down on the size of the pose.
  if (!options::option[in::olig_search::dump_pdb]()) {
    for (Size i=1; i<=pose.n_residue(); i++) {
      if (chemical::name_from_aa(pose.aa(i)) == "GLY") {
        core::pose::replace_pose_residue_copying_existing_coordinates(pose, i, resi_set->name_map("GLY"));
      } else {
        core::pose::replace_pose_residue_copying_existing_coordinates(pose, i, resi_set->name_map("ALA"));
      }
    }
  }

  // Handle all of the symmetry stuff
//  pose.dump_pdb("pre_symm.pdb");
  core::pose::symmetry::make_symmetric_pose(pose);
//  pose.dump_pdb("pst_symm.pdb");

  SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
  std::map<Size,SymDof> dofs = sym_info->get_dofs();
  int sym_jump = 0;

  for(std::map<Size,SymDof>::iterator i = dofs.begin(); i != dofs.end(); i++) {
    Size jump_num = i->first;
    if (sym_jump == 0) {
      sym_jump = jump_num;
    } else {
      utility_exit_with_message("ERROR: Can only handle one subunit!");
    }
  }

  if (sym_jump == 0) {
    utility_exit_with_message("ERROR: No jump defined!");
  }

  // If in::olig_search::dump_pdb is set, dump the pdb and get the hell outta here.
  if (options::option[in::olig_search::dump_pdb]()) {
		utility::vector1<Real> radial_disps = options::option[in::olig_search::radial_disp]();
		utility::vector1<Real> angles = options::option[in::olig_search::angle]();
		if (radial_disps.size() != angles.size()) {
			utility_exit_with_message("ERROR: radial_disp and angle input vectors are not the same length!");
		} else {
			for (Size i=1; i<=radial_disps.size(); i++) {
		    move_pose(pose, sym_jump, Vec(radial_disps[i],0,0), angles[i]+1);
				// If desired, only spit out chain A
				if (options::option[in::olig_search::dump_chainA_only]()) {
					std::string fn = options::option[in::olig_search::prefix]()+options::option[in::olig_search::pdbID]()+"_"+string_of(radial_disps[i])+"_"+string_of(angles[i])+"_A.pdb";
  	      utility::io::ozstream out( fn );
					Size atom_counter = 1;
					for (Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
				    Size isub = sym_info->subunit_index(ir);
				    if (isub != 1) break;
		        core::io::pdb::dump_pdb_residue(pose.residue(ir),atom_counter,out);
					}
	        out.close();
				// Otherwise, spit out the entire assembly
				} else {
					pose.dump_pdb(options::option[in::olig_search::prefix]()+options::option[in::olig_search::pdbID]()+"_"+string_of(radial_disps[i])+"_"+string_of(angles[i])+".pdb");
				}
				move_pose(pose, sym_jump, Vec(-radial_disps[i],0,0), -(angles[i]+1));
			}
    return(1);
		}
  }

  // Move the sucker, and spit out the clash/counts
//  int offset_dist = 1; // Negative sign flips orientation of building block
//  int offset_dist = -1; // Negative sign flips orientation of building block
  move_pose(pose, sym_jump, Vec(offset,0.0,0.0), 0.0);
  Real radial_disp = 0;
  for (int r=0; r<(360/nsub_bb); r++) { // Will need to move away from hard-coding eventually
    move_pose(pose, sym_jump, Vec(-radial_disp,0.0,0.0), 1.0);
    radial_disp = 0;
    int int_cb_count = -1;
    while (int_cb_count != 0) {
      int tmp_int_cb_count = 0;
      int_cb_count = 0;
//      int dt = 1; // distance moved; negative sign flips orientation
//      int dt = -1; // distance moved; negative sign flips orientation
      move_pose(pose, sym_jump, Vec(dt,0.0,0.0));
      radial_disp = radial_disp + dt;
      for (Size j=(nsub_bb+1); j<=(nsub_total); j++) { // This is architecture-dependent. Need to generalize.
        tmp_int_cb_count = count_int_CBs_clashes(pose, 1, j, 10.0, 3.0);
        if (tmp_int_cb_count == -1) {
          int_cb_count = tmp_int_cb_count;
          break;
        } else {
          int_cb_count += tmp_int_cb_count;
        }
      }
      Real total_radial_disp = radial_disp + offset;
      TR<< string_of(total_radial_disp) << "\t" << r << "\t" << int_cb_count << std::endl;
      // If you don't dump the pdb, it goes much faster.
/*
      if ((int_cb_count != -1) && (int_cb_count != 0)) {
      if (int_cb_count != 0) {
        std::stringstream ss;
        ss << "32_3_tst_" << total_radial_disp << "_" << r << ".pdb";
        pose.dump_pdb(ss.str());
      }
      }
*/
    }
  }

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
}
