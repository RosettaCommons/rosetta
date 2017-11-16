// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk

#include <iostream>
#include <iomanip>

#include <devel/init.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/util.hh>

#include <core/io/pdb/build_pose_as_is.hh>

#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicWrapper.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>

#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/after_opts.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>


using namespace core;
using namespace basic::options;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;

static basic::Tracer TR( "apps.pilot.constellation_list.main" );


utility::vector1<char> list_allowable_mutations( char const starting_aa ) {

	utility::vector1<char> allowed_list;
	// TR << "Found residue " << chemical::oneletter_code_from_aa( aa ) << std::endl;

	if ( starting_aa == 'G' ) return allowed_list;

	// general cases
	// anything other than Gly can become Gly, anything other than Gly/Ala can become Ala
	allowed_list.push_back('G');
	if ( starting_aa != 'A' ) {
		allowed_list.push_back('A');
	}

	// special cases
	// Thr can become Ser
	if ( starting_aa == 'T' ) {
		allowed_list.push_back('S');
	}
	// Ile can become Val
	if ( starting_aa == 'I' ) {
		allowed_list.push_back('V');
	}
	// Tyr can become Phe
	if ( starting_aa == 'Y' ) {
		allowed_list.push_back('F');
		allowed_list.push_back('L');
	}
	// Phe can become Leu
	if ( starting_aa == 'F' ) {
		allowed_list.push_back('L');
	}
	// Trp can become Leu
	if ( starting_aa == 'W' ) {
		allowed_list.push_back('L');
	}

	// JIMMY ADD MORE HERE, eg. F->L, W->L

	return allowed_list;
}


void zero_occ_for_deleted_atoms(
	pose::Pose & pose,
	core::Size seqpos,
	char const target_aa
) {

	// set the occupancy to zero for any target residue atoms that don't need to be printed
	// leave the occupancy at one for any target residue atoms that *do* need to be printed

	char const starting_aa = chemical::oneletter_code_from_aa( pose.aa(seqpos) );

	conformation::Residue const & rsd( pose.residue(seqpos) );

	// we never need to print backbone atoms or hydrogens
	for ( Size i=1; i<= rsd.natoms(); ++i ) {
		if ( rsd.atom_is_hydrogen(i) || rsd.atom_is_backbone(i) ) {
			pose.pdb_info()->occupancy( seqpos, i, 0. );
		}
	}

	// if it's a mutation to Gly, we're done
	if ( target_aa == 'G' ) return;

	// for anything other than a mutation to Gly, suppress the C-beta
	pose.pdb_info()->occupancy( seqpos, rsd.atom_index("CB"), 0. );
	// if it's a mutation to Ala, we're done
	if ( target_aa == 'A' ) return;

	if ( ( starting_aa == 'I' ) && ( target_aa == 'V' ) ) {
		// suppress everything other than the CD1
		Size atom_inx_to_keep = rsd.atom_index("CD1");
		for ( Size i=1; i<= rsd.natoms(); ++i ) {
			if ( i != atom_inx_to_keep ) {
				pose.pdb_info()->occupancy( seqpos, i, 0. );
			}
		}
		return;
	}

	if ( ( starting_aa == 'T' ) && ( target_aa == 'S' ) ) {
		// suppress everything other than the CG2
		Size atom_inx_to_keep = rsd.atom_index("CG2");
		for ( Size i=1; i<= rsd.natoms(); ++i ) {
			if ( i != atom_inx_to_keep ) {
				pose.pdb_info()->occupancy( seqpos, i, 0. );
			}
		}
		return;
	}

	if ( ( starting_aa == 'Y' ) && ( target_aa == 'F' ) ) {
		// suppress everything other than the OH
		Size atom_inx_to_keep = rsd.atom_index("OH");
		for ( Size i=1; i<= rsd.natoms(); ++i ) {
			if ( i != atom_inx_to_keep ) {
				pose.pdb_info()->occupancy( seqpos, i, 0. );
			}
		}
		return;
	}

	if ( ( starting_aa == 'Y' ) && ( target_aa == 'L' ) ) {
		// suppress everything other than the OH, CE1, CE2, CZ
		Size inx1 = rsd.atom_index("OH");
		Size inx2 = rsd.atom_index("CE1");
		Size inx3 = rsd.atom_index("CE2");
		Size inx4 = rsd.atom_index("CZ");
		for ( Size i=1; i<= rsd.natoms(); ++i ) {
			if ( (i != inx1) && (i != inx2) && (i != inx3) && (i !=inx4) ) {
				pose.pdb_info()->occupancy( seqpos, i, 0. );
			}
		}
		return;
	}

	if ( ( starting_aa == 'F' ) && ( target_aa == 'L' ) ) {
		// suppress everything other than the CE1, CE2, CZ
		Size inx1 = rsd.atom_index("CE1");
		Size inx2 = rsd.atom_index("CE2");
		Size inx3 = rsd.atom_index("CZ");
		for ( Size i=1; i<= rsd.natoms(); ++i ) {
			if ( (i != inx1) && (i != inx2) && (i != inx3) ) {
				pose.pdb_info()->occupancy( seqpos, i, 0. );
			}
		}
		return;
	}

	if ( ( starting_aa == 'W' ) && ( target_aa == 'L' ) ) {
		// suppress everything other than the NE1, CE2, CE3, CZ2, CZ3, CH2
		Size inx1 = rsd.atom_index("NE1");
		Size inx2 = rsd.atom_index("CE2");
		Size inx3 = rsd.atom_index("CE3");
		Size inx4 = rsd.atom_index("CZ2");
		Size inx5 = rsd.atom_index("CZ3");
		Size inx6 = rsd.atom_index("CH2");
		for ( Size i=1; i<= rsd.natoms(); ++i ) {
			if ( (i != inx1) && (i != inx2) && (i != inx3) && (i != inx4) && (i != inx5) && (i != inx6) ) {
				pose.pdb_info()->occupancy( seqpos, i, 0. );
			}
		}
		return;
	}

	// JIMMY ADD MORE HERE, STARTING WITH T->S AND Y->F

	TR << "DANGER DANGER - COULD NOT MAKE A REQUESTED MUTATION!!" << std::endl;
	TR << "REQUESTED " << starting_aa << " TO " << target_aa << " BUT NO INFO ON WHAT ATOMS TO SUPPRESS....." << std::endl;

	return;

}

OPT_KEY( Integer, target_resnum )
OPT_KEY( String, target_chain )


/// General testing code
int
main( int argc, char * argv [] )
{

	NEW_OPT( target_resnum, "this residue must be included, presumably because it impacts function", -1 );
	NEW_OPT( target_chain, "which chain the target residue is on", "A" );

	devel::init(argc, argv);

	TR << "Starting constellation_list" << std::endl;

	std::string const tmp_chain = option[ target_chain ];
	if ( tmp_chain.length() != 1 ) {
		TR << "ERROR!! Chain ID should be one character" << std::endl;
		exit(1);
	}
	char const target_pdb_chain = tmp_chain[0];
	int const target_pdb_number = option[ target_resnum ];

	TR << "Target residue is " << target_pdb_number << " on chain " << target_pdb_chain << std::endl;

	// create pose for native pose from pdb
	pose::Pose pose_init;
	std::string const input_pdb_name ( basic::options::start_file() );
	core::import_pose::pose_from_file( pose_init, input_pdb_name , core::import_pose::PDB_file);

	// set target_rosetta_resnum to Rosetta internal resid for the residue to be mutated
	core::Size target_rosetta_resnum = 0;
	for ( core::Size j = 1; j <= pose_init.size(); ++j ) {
		if ( ( pose_init.pdb_info()->chain(j) == target_pdb_chain ) &&
				( pose_init.pdb_info()->number(j) == target_pdb_number ) ) {
			target_rosetta_resnum = j;
		}
	}
	if ( target_rosetta_resnum == 0 ) {
		TR << "ERROR!! Could not find residue/chain" << std::endl;
		exit(1);
	}

	// scoring function - for speed use a custom scorefxn containing only fa_atr, no need to evaluate other terms
	scoring::ScoreFunctionOP scorefxn( scoring::get_score_function() );
	scorefxn->reset();
	core::Real const fa_atr_weight = 0.8;
	scorefxn->set_weight(fa_atr,fa_atr_weight);
	(*scorefxn)( pose_init );

	// make a list of residues with vdw atr contacts (to the sidechain only, if possible?)
	// note that the structure is already scored (above)
	utility::vector1<bool> interacting_residue( pose_init.size(), false );

	core::scoring::EnergyMap unweighted_emap;
	core::conformation::Residue target_residue = pose_init.residue( target_rosetta_resnum );
	core::Real const interaction_score_threshold = -0.3;
	for ( Size i = 1; i <= pose_init.size(); ++i ) {
		if ( i == target_rosetta_resnum ) continue;
		unweighted_emap.clear();
		Vector const tgt_centroid = core::scoring::compute_sc_centroid(target_residue);
		Vector const pi_centroid = core::scoring::compute_sc_centroid(pose_init.residue(i));
		Real const tgt_rad = core::scoring::compute_sc_radius( target_residue, tgt_centroid );
		Real const pi_rad = core::scoring::compute_sc_radius( pose_init.residue(i), pi_centroid );
		core::scoring::eval_scsc_sr2b_energies( target_residue, pose_init.residue(i), tgt_centroid, pi_centroid,
			tgt_rad, pi_rad, pose_init, *scorefxn, unweighted_emap );
		Real const weighted_fa_atr = fa_atr_weight * unweighted_emap[ fa_atr ];
		if ( weighted_fa_atr < interaction_score_threshold ) interacting_residue.at(i)=true;
	}

	core::Size constellation_number = 0;

	// for each two-residue combination, loop over all allowed mutations at each position
	for ( Size j = 1; j <= pose_init.size(); ++j ) {
		if ( interacting_residue.at(j) ) {
			TR << "Building constellations with target and rosetta resnum:  " << j << std::endl;

			utility::vector1<char> allowable_target_mutations =
				list_allowable_mutations( chemical::oneletter_code_from_aa( pose_init.aa(target_rosetta_resnum) ) );
			utility::vector1<char> allowable_secondary_mutations =
				list_allowable_mutations( chemical::oneletter_code_from_aa( pose_init.aa(j) ) );

			for ( Size tmut=1; tmut <= allowable_target_mutations.size(); ++tmut ) {
				pose::Pose target_mut_pose = pose_init;
				zero_occ_for_deleted_atoms( target_mut_pose, target_rosetta_resnum, allowable_target_mutations.at(tmut) );

				for ( Size jmut=1; jmut <= allowable_secondary_mutations.size(); ++jmut ) {
					pose::Pose secondary_mut_pose = target_mut_pose;
					zero_occ_for_deleted_atoms( secondary_mut_pose, j, allowable_secondary_mutations.at(jmut) );

					++constellation_number;
					std::ostringstream outConstel_name;
					outConstel_name << "constel_" << target_pdb_number << "_" << constellation_number << ".pdb";
					utility::io::ozstream outConstel_stream;
					outConstel_stream.open(outConstel_name.str(), std::ios::out);
					outConstel_stream << "HEADER   CONST NUM " << constellation_number << " TARGET MUTATION: " << chemical::oneletter_code_from_aa( pose_init.aa(target_rosetta_resnum) ) << target_pdb_number << allowable_target_mutations.at(tmut) << "  SECONDARY_MUTATION: " << chemical::oneletter_code_from_aa( pose_init.aa(j) ) << j << allowable_secondary_mutations.at(jmut) << std::endl;

					// print the atoms that would be removed by this mutation to a pdb name containing constellation_number
					core::io::pdb::FileData fd;
					std::string data;
					utility::vector1< core::Size > residues_to_print;
					residues_to_print.push_back(target_rosetta_resnum);
					residues_to_print.push_back(j);
					fd.init_from_pose( secondary_mut_pose, residues_to_print );
					data = core::io::pdb::PDB_DReader::createPDBData(fd);
					outConstel_stream.write( data.c_str(), data.size() );
					outConstel_stream.close();
					outConstel_stream.clear();

					// make this set of mutations to the protein, print this out
					// JK note: currently we make the mutation by calling pack_rotamers,
					// but this doesn't necessarily put the new sidechain where the old atoms used to be (eg. for F->L mutation).
					// Is this good or bad? Bad in that a matching small-mol won't fit,
					// but good in that the new sidechain may not adopt that conformation so we wouldn't want to test that one...?
					pose::Pose mut_pose_init = pose_init;
					// setup a packer task
					pack::task::PackerTaskOP mut_task( pack::task::TaskFactory::create_packer_task( mut_pose_init ));
					mut_task->set_bump_check( false );
					mut_task->initialize_from_command_line();
					mut_task->or_include_current( true );
					// restrict packer task to sequence positions of interest
					utility::vector1<bool> allow_redesign( mut_pose_init.size(), false );
					allow_redesign.at(target_rosetta_resnum) = true;
					allow_redesign.at(j) = true;
					mut_task->restrict_to_residues( allow_redesign );
					// set the two mutations to make
					utility::vector1< bool > target_mut_site( core::chemical::num_canonical_aas, false );
					target_mut_site.at( chemical::aa_from_oneletter_code( allowable_target_mutations.at(tmut) ) ) = true;
					mut_task->nonconst_residue_task(target_rosetta_resnum).restrict_absent_canonical_aas( target_mut_site );
					utility::vector1< bool > secondary_mut_site( core::chemical::num_canonical_aas, false );
					secondary_mut_site.at( chemical::aa_from_oneletter_code( allowable_secondary_mutations.at(jmut) ) ) = true;
					mut_task->nonconst_residue_task(j).restrict_absent_canonical_aas( secondary_mut_site );
					// call pack_rotamers
					scoring::ScoreFunctionOP repack_scorefxn(get_score_function());
					pack::pack_rotamers( mut_pose_init, *repack_scorefxn, mut_task );
					// write the output pdb
					std::ostringstream outProtein_name;
					outProtein_name << "protein_" << target_pdb_number << "_" << constellation_number << ".pdb";
					mut_pose_init.dump_pdb( outProtein_name.str() );

				}

			}

		}
	}

	TR << "Remember to set the out::file::suppress_zero_occ_pdb_output flag !!!" << std::endl;
	TR << "Done, found " << constellation_number << " total constellations." << std::endl;

	return 0;


}
