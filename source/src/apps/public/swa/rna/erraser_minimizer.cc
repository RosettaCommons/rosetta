// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief

#include <protocols/rna/RNA_ProtocolUtil.hh>
#include <protocols/rna/RNA_IdealCoord.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

#include <core/init/init.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>

#include <protocols/viewer/viewers.hh>

#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>

//////////////////////////////////////////////////
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
///////////////////////////////////////////////////

#include <core/pose/PDBInfo.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/types.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
//////////////////////////////////////////////////////////


// C++ headers
#include <iostream>
#include <string>

#include <utility/excn/Exceptions.hh>

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;

OPT_KEY ( String, out_pdb )
OPT_KEY ( Boolean, vary_geometry )
OPT_KEY ( Boolean, constrain_P )
OPT_KEY ( Boolean, ready_set_only )
OPT_KEY ( Boolean, skip_minimize )
OPT_KEY ( Boolean, attempt_pyrimidine_flip )
OPT_KEY ( IntegerVector, fixed_res )
OPT_KEY ( IntegerVector, cutpoint_open )

////////////////////////////////////////////////////////////////////////
bool
check_num_in_vector ( int input_num, utility::vector1< int > const & input_vector  ) {
	for (Size i = 1; i <= input_vector.size(); ++i) {
		if (input_num == input_vector[i]) return true;
	}
	return false;
}
////////////////////////////////////////////////////////////////////////
void
translate_residue ( conformation::Residue & rsd,
                    Vector const & nbr_atom_xyz ) {
	Vector const translate ( nbr_atom_xyz - rsd.nbr_atom_xyz() );

	for ( Size i = 1; i <= rsd.natoms(); ++i ) {
		rsd.set_xyz ( i, rsd.xyz ( i ) + translate );
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////
bool
check_in_bonded_list ( core::id::AtomID const & atom_id1,
                       core::id::AtomID const & atom_id2,
                       utility::vector1< std::pair< core::id::AtomID, core::id::AtomID > > & bonded_atom_list ) {
	for ( Size n = 1; n <= bonded_atom_list.size(); n++ ) {
		if ( atom_id1 == bonded_atom_list[ n ].first && atom_id2 == bonded_atom_list[ n ].second ) return true;

		if ( atom_id2 == bonded_atom_list[ n ].first && atom_id1 == bonded_atom_list[ n ].second ) return true;
	}

	return false;
}



//////////////////////////////////////////////////////////////////////////////////////////////
bool
check_in_bond_angle_list ( core::id::AtomID const & atom_id1,
                           core::id::AtomID const & atom_id2,
                           core::id::AtomID const & atom_id3,
                           utility::vector1< std::pair< core::id::AtomID, std::pair< core::id::AtomID, core::id::AtomID > > > & bond_angle_list ) {
	for ( Size n = 1; n <= bond_angle_list.size(); n++ ) {
		if ( atom_id1 == bond_angle_list[ n ].first ) {
			if ( atom_id2 == bond_angle_list[ n ].second.first && atom_id3 == bond_angle_list[ n ].second.second ) return true;

			if ( atom_id3 == bond_angle_list[ n ].second.first && atom_id2 == bond_angle_list[ n ].second.second ) return true;
		}
	}

	return false;
}
//////////////////////////////////////////////////////////////////////////////////////////////
void
apply_ideal_coordinates ( pose::Pose const & pose, pose::Pose & pose_reference ) {
	using namespace core::chemical::rna;
	RNA_FittedTorsionInfo const rna_fitted_torsion_info;
	Real const DELTA_CUTOFF ( rna_fitted_torsion_info.delta_cutoff() );
	bool const is_use_phenix_geo = option[ basic::options::OptionKeys::rna::corrected_geo ];
	utility::vector1 <Size> pucker_conformation (pose_reference.total_residue(), 0);
		
	protocols::rna::RNA_IdealCoord ideal_coord;
	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
		if ( pose.residue ( n ).aa() == core::chemical::aa_vrt ) continue; //FCC

		Real const delta = pose.residue ( n ).mainchain_torsion ( DELTA );

		if ( delta > DELTA_CUTOFF ) { //south
			apply_ideal_c2endo_sugar_coords ( pose_reference, n );
			pucker_conformation[n] = 2; 
		} else { //north
			pucker_conformation[n] = 1;
		}
	}
	if (is_use_phenix_geo) {
		ideal_coord.apply(pose_reference, pucker_conformation, false /*donot keep torsions*/);
	}
}
////////////////////////////////////////////////////////////////
void
add_bond_constraint ( core::id::AtomID const & atom_id1,
                      core::id::AtomID const & atom_id2,
                      utility::vector1< std::pair< core::id::AtomID, core::id::AtomID > > & bonded_atom_list,
                      core::pose::Pose const & pose,
                      core::pose::Pose const & pose_reference,
                      core::scoring::constraints::ConstraintSetOP & cst_set ) {
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	std::string const & atom_name1 = pose.residue ( atom_id1.rsd() ).atom_name ( atom_id1.atomno() );
	std::string const & atom_name2 = pose.residue ( atom_id2.rsd() ).atom_name ( atom_id2.atomno() );

	if ( !pose_reference.residue ( atom_id1.rsd() ).has ( atom_name1 ) ) return;

	if ( !pose_reference.residue ( atom_id2.rsd() ).has ( atom_name2 ) ) return;

	if ( !check_in_bonded_list ( atom_id1, atom_id2,  bonded_atom_list ) ) {
		bonded_atom_list.push_back ( std::make_pair ( atom_id1, atom_id2 ) );
		Real const bond_length_sd_ ( 0.05 );
		Real const bond_length = ( pose_reference.residue ( atom_id1.rsd() ).xyz ( atom_name1 ) -
		                           pose_reference.residue ( atom_id2.rsd() ).xyz ( atom_name2 ) ).length();
		FuncOP dist_harm_func_ ( new HarmonicFunc ( bond_length, bond_length_sd_ ) );
		cst_set->add_constraint ( new AtomPairConstraint ( atom_id1 ,
		                          atom_id2,
		                          dist_harm_func_,
		                          rna_bond_geometry ) );

		if ( false ) {
			std::cout << "PUTTING CONSTRAINT ON DISTANCE: " <<
			          atom_id2.rsd() << " " << atom_name1 << "; "  <<
			          atom_id1.rsd() << " " << atom_name2 << " "  <<
			          bond_length <<
			          std::endl;
		}
	}
}
///////////////////////////////////////////////////////
void
add_bond_angle_constraint ( core::id::AtomID const & atom_id1,
                            core::id::AtomID const & atom_id2,
                            core::id::AtomID const & atom_id3,
                            utility::vector1< std::pair < core::id::AtomID, std::pair< core::id::AtomID, core::id::AtomID > > > & bond_angle_list,
                            core::pose::Pose const & pose,
                            core::pose::Pose const & pose_reference,
                            core::scoring::constraints::ConstraintSetOP & cst_set ) {
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace numeric::conversions;

	if ( atom_id2 == atom_id3 ) return;

	if ( atom_id1 == atom_id3 ) return;

	if ( atom_id1 == atom_id2 ) return;

	std::string const & atom_name1 = pose.residue ( atom_id1.rsd() ).atom_name ( atom_id1.atomno() );
	std::string const & atom_name2 = pose.residue ( atom_id2.rsd() ).atom_name ( atom_id2.atomno() );
	std::string const & atom_name3 = pose.residue ( atom_id3.rsd() ).atom_name ( atom_id3.atomno() );

	if ( !pose_reference.residue ( atom_id1.rsd() ).has ( atom_name1 ) ) return;

	if ( !pose_reference.residue ( atom_id2.rsd() ).has ( atom_name2 ) ) return;

	if ( !pose_reference.residue ( atom_id3.rsd() ).has ( atom_name3 ) ) return;

	if ( !check_in_bond_angle_list ( atom_id1, atom_id2, atom_id3, bond_angle_list ) ) {
		bond_angle_list.push_back ( std::make_pair ( atom_id1, std::make_pair ( atom_id2, atom_id3 ) ) );
		Real const bond_angle_sd_ ( radians ( 3.0 ) );
		Real const bond_angle = angle_radians (
		                          pose_reference.residue ( atom_id2.rsd() ).xyz ( atom_name2 )  ,
		                          pose_reference.residue ( atom_id1.rsd() ).xyz ( atom_name1 )  ,
		                          pose_reference.residue ( atom_id3.rsd() ).xyz ( atom_name3 )
		                        );

		if ( bond_angle < 0.001 ) std::cout << "WHAT THE HELL????????? " << std::endl;

		FuncOP angle_harm_func_ ( new HarmonicFunc ( bond_angle, bond_angle_sd_ ) );
		cst_set->add_constraint ( new AngleConstraint (
		                            atom_id2 , atom_id1, atom_id3, angle_harm_func_,	rna_bond_geometry ) );

		if ( false ) {
			std::cout << "PUTTING CONSTRAINT ON ANGLE: " <<
			          atom_id2.rsd() << " " << pose_reference.residue ( atom_id2.rsd() ).atom_name ( atom_id2.atomno() ) << "; "  <<
			          atom_id1.rsd() << " " << pose_reference.residue ( atom_id1.rsd() ).atom_name ( atom_id1.atomno() ) << "; "  <<
			          atom_id3.rsd() << " " << pose_reference.residue ( atom_id3.rsd() ).atom_name ( atom_id3.atomno() ) << " ==> "  << degrees ( bond_angle ) << " " << degrees ( bond_angle_sd_ ) <<
			          std::endl;
		}
	}
}

////////////////////////////////////////////////
bool
check_if_really_connected (
  core::pose::Pose const & pose,
  core::id::AtomID const & atom_id1,
  core::id::AtomID const & atom_id2 ) {
	if ( atom_id1.rsd() == atom_id2.rsd() ) return true;

	core::kinematics::tree::AtomCOP atom1 ( & pose.atom_tree().atom ( atom_id1 ) );
	core::kinematics::tree::AtomCOP atom2 ( & pose.atom_tree().atom ( atom_id2 ) );

	if ( atom1->parent() == atom2 ) return true;

	if ( atom2->parent() == atom1 ) return true;

	return false;
}

////////////////////////////////////////////////
bool
i_want_this_atom_to_move ( conformation::Residue const & residue2, Size const & k ) {
	if ( k > residue2.first_sidechain_atom() &&
	     k != chemical::rna::first_base_atom_index ( residue2 ) ) return false;

	if ( residue2.atom_type ( k ).name() == "VIRT" ) {
		//		std::cout << "Is this virtual? " << residue2.atom_name( k ) << std::endl;
		return false;
	}

	return true;
}

bool
i_want_this_atom_to_move ( pose::Pose const & pose, core::id::AtomID const & atom_id ) {
	return i_want_this_atom_to_move ( pose.residue ( atom_id.rsd() ) ,
	                                  atom_id.atomno() );
}

bool
is_atom_exist_in_reference ( pose::Pose const & pose, pose::Pose const & pose_reference, core::id::AtomID const & atom_id ) {
	std::string const & atom_name = pose.residue ( atom_id.rsd() ).atom_name ( atom_id.atomno() );

	if ( pose_reference.residue ( atom_id.rsd() ).has ( atom_name ) ) {
		return true;
	} else {
		std::cout << atom_name << std::endl;
		return false;
	}
}
//////////////////////////////////////////////////////////////////////////////
void
create_pose_reference (
  pose::Pose const & pose,
  pose::Pose & pose_reference ) {
	using namespace core::chemical;
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set ( RNA );
	make_pose_from_sequence ( pose_reference, pose.sequence(),	*rsd_set );
	apply_ideal_coordinates ( pose, pose_reference );
}

//////////////////////////////////////////////////////////////////////////////
// Following has not (yet) been carefully debugged.
void
vary_bond_geometry (
  core::kinematics::MoveMap & mm,
  pose::Pose & pose,
  pose::Pose const & pose_reference,
  ObjexxFCL::FArray1D< bool > & allow_insert_ ) {
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace numeric::conversions;
	ConstraintSetOP cst_set = pose.constraint_set()->clone();
	pose.constraint_set ( cst_set );
	Size const nres ( pose.total_residue() );
	std::map< AtomID, utility::vector1< AtomID > > lists_of_angle_bonded_atoms;
	std::cout << "Enter the vary_bond_geometry....." << std::endl;

	for ( Size i = 1; i <= nres; i++ )  {
		if ( pose.residue ( i ).aa() == core::chemical::aa_vrt ) continue; //FCC

		if ( !allow_insert_ ( i ) ) continue;

		conformation::Residue const & residue ( pose.residue ( i ) );

		for ( Size j = 1; j <= residue.natoms(); j++ ) {
			if ( !i_want_this_atom_to_move ( residue, j ) ) continue;

			if ( !is_atom_exist_in_reference ( pose, pose_reference , AtomID ( j, i ) ) )  continue;

			core::kinematics::tree::AtomCOP current_atom ( & pose.atom_tree().atom ( AtomID ( j, i ) ) );

			if ( current_atom->is_jump() ) continue;

			///////////////////
			core::kinematics::tree::AtomCOP input_stub_atom1 ( current_atom->input_stub_atom1() );

			if ( !input_stub_atom1 ) continue;

			if ( !i_want_this_atom_to_move ( pose, input_stub_atom1->id() ) ) continue;

			if ( !is_atom_exist_in_reference ( pose, pose_reference, input_stub_atom1->id() ) ) continue;

			mm.set ( DOF_ID ( AtomID ( j, i ), D ), true );

			if ( input_stub_atom1->is_jump() ) continue;

			core::kinematics::tree::AtomCOP input_stub_atom2 ( current_atom->input_stub_atom2() );

			///////////////////
			if ( !input_stub_atom2 ) continue;

			if ( input_stub_atom2 == current_atom ) continue;

			if ( !i_want_this_atom_to_move ( pose, input_stub_atom2->id() ) ) continue;

			if ( !is_atom_exist_in_reference ( pose, pose_reference, input_stub_atom2->id() ) ) continue;

			mm.set ( DOF_ID ( AtomID ( j, i ), THETA ), true );

			if ( input_stub_atom2->is_jump() ) continue;

			///////////////////
			core::kinematics::tree::AtomCOP input_stub_atom3 ( current_atom->input_stub_atom3() );

			if ( !input_stub_atom3 ) continue;

			if ( !i_want_this_atom_to_move ( pose, input_stub_atom3->id() ) ) continue;

			if ( !is_atom_exist_in_reference ( pose, pose_reference, input_stub_atom3->id() ) ) continue;

			if ( input_stub_atom3 == current_atom ) continue;

			mm.set ( DOF_ID ( AtomID ( j, i ), PHI ), true );
		}
	}

	utility::vector1< std::pair< AtomID, AtomID > >  bond_list;
	utility::vector1< std::pair< AtomID, std::pair< AtomID, AtomID > > > bond_angle_list;

	for ( Size i = 1; i <= nres; i++ )  {
		if ( pose.residue ( i ).aa() == core::chemical::aa_vrt ) continue; //FCC

		//Go through all bonds in pose...
		if ( !allow_insert_ ( i ) ) continue;

		conformation::Residue const & residue ( pose.residue ( i ) );

		for ( Size j = 1; j <= residue.natoms(); j++ ) {
			if ( !i_want_this_atom_to_move ( residue, j ) ) continue;

			AtomID atom_id1 ( j, i );
			utility::vector1< AtomID >  nbrs ( pose.conformation().bonded_neighbor_all_res ( atom_id1 ) );

			// Bond lengths.
			for ( Size n = 1; n <= nbrs.size(); n++ ) {
				AtomID const & atom_id2 ( nbrs[ n ] );
				conformation::Residue const & residue2 ( pose.residue ( atom_id2.rsd() ) ) ;
				Size const & k ( atom_id2.atomno() ) ;

				if ( ! check_if_really_connected ( pose, atom_id1, atom_id2 ) ) continue;

				if ( i_want_this_atom_to_move ( residue2, k ) )  {
					add_bond_constraint ( atom_id1, atom_id2,
					                      bond_list,
					                      pose, pose_reference, cst_set );
				}
			}

			// Bond angles
			for ( Size m = 1; m <= nbrs.size(); m++ ) {
				AtomID const & atom_id2 ( nbrs[ m ] );
				conformation::Residue const & residue2 ( pose.residue ( atom_id2.rsd() ) ) ;

				if ( ! check_if_really_connected ( pose, atom_id1, atom_id2 ) ) continue;

				Size const & k ( atom_id2.atomno() ) ;

				for ( Size n = 1; n <= nbrs.size(); n++ ) {
					AtomID const & atom_id3 ( nbrs[ n ] );
					conformation::Residue const & residue3 ( pose.residue ( atom_id3.rsd() ) ) ;

					if ( ! check_if_really_connected ( pose, atom_id1, atom_id3 ) ) continue;

					Size const & q ( atom_id3.atomno() ) ;

					if ( i_want_this_atom_to_move ( residue2, k ) &&
					     i_want_this_atom_to_move ( residue3, q ) )  {
						add_bond_angle_constraint ( atom_id1, atom_id2, atom_id3, bond_angle_list,
						                            pose, pose_reference, cst_set );
					}
				}

				utility::vector1< AtomID >  nbrs2 ( pose.conformation().bonded_neighbor_all_res ( atom_id2 ) );

				for ( Size n = 1; n <= nbrs2.size(); n++ ) {
					AtomID const & atom_id3 ( nbrs2[ n ] );
					conformation::Residue const & residue3 ( pose.residue ( atom_id3.rsd() ) ) ;

					if ( ! check_if_really_connected ( pose, atom_id2, atom_id3 ) ) continue;

					Size const & q ( atom_id3.atomno() ) ;

					if ( i_want_this_atom_to_move ( residue2, k ) &&
					     i_want_this_atom_to_move ( residue3, q ) )  {
						add_bond_angle_constraint ( atom_id1, atom_id2, atom_id3, bond_angle_list,
						                            pose, pose_reference, cst_set );
					}
				}
			}
		}
	}

	pose.constraint_set ( cst_set );
}
/////////////////////////////////////////////////////////////////////////////
//FCC: Adding Virtual res
int
add_virtual_res ( core::pose::Pose & pose ) {
	int nres = pose.total_residue();

	// if already rooted on virtual residue , return
	if ( pose.residue ( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) {
		std::cout << "add_virtual_res() called but pose is already rooted on a VRT residue ... continuing." << std::endl;
		return  pose.fold_tree().root();
	}

	// attach virt res there
	core::chemical::ResidueTypeSet const & residue_set = pose.residue_type ( 1 ).residue_type_set();
	core::chemical::ResidueTypeCOPs const & rsd_type_list ( residue_set.name3_map ( "VRT" ) );
	core::conformation::ResidueOP new_res ( core::conformation::ResidueFactory::create_residue ( *rsd_type_list[1] ) );
	pose.append_residue_by_jump ( *new_res , nres );
	// make the virt atom the root
	kinematics::FoldTree newF ( pose.fold_tree() );
	newF.reorder ( nres + 1 );
	pose.fold_tree ( newF );
	return ( nres + 1 );
}
///////////////////////////////////////////////////////////
void
setup_fold_tree ( pose::Pose & pose, utility::vector1< core::Size > const & cutpoint_list ) {
	using namespace chemical;
	using namespace core::scoring;
	using namespace core::chemical::rna;
	using namespace core::conformation;
	using namespace core::kinematics;
	using namespace core::id;
	Size const nres ( pose.total_residue() );
	Size const num_jumps ( cutpoint_list.size() );
	ObjexxFCL::FArray2D <int> jump_points ( 2, num_jumps );
	ObjexxFCL::FArray1D <int> cuts ( num_jumps );

	for ( Size n = 1; n <= cutpoint_list.size(); n++ ) {
		jump_points ( 1, n ) = cutpoint_list[n];
		jump_points ( 2, n ) = cutpoint_list[n] + 1;
		cuts ( n ) = cutpoint_list[n];
	}

	FoldTree f ( nres );
	f.tree_from_jumps_and_cuts ( nres, num_jumps, jump_points, cuts, 1, false );
	pose.fold_tree ( f );
}
///////////////////////////////////////////////////////////
bool
is_elem_in_list ( core::Size const elem, utility::vector1< core::Size > const & list ) {
	for ( Size i = 1; i <= list.size(); ++i ) {
		if ( elem == list[i] ) return true;
	}

	return false;
}
///////////////////////////////////////////////
void
setup_fold_tree_sample_res ( pose::Pose & pose, utility::vector1< core::Size > const & sample_res_list ) {
	using namespace chemical;
	using namespace core::scoring;
	using namespace core::conformation;
	using namespace core::kinematics;
	using namespace core::id;
	Size const nres ( pose.total_residue() );
	FoldTree f ( nres );
	bool is_segment_sample_res = is_elem_in_list ( 1, sample_res_list );

	for ( Size i = 2; i < nres; ++i ) {
		bool is_res_sample_res = is_elem_in_list ( i, sample_res_list );

		if ( is_res_sample_res != is_segment_sample_res ) {
			is_segment_sample_res = is_res_sample_res;
			f.new_jump ( i - 1, nres, i - 1 );
		}
	}

	f.new_jump ( nres - 1, nres, nres - 1 );
	f.reorder ( nres );
	pose.fold_tree ( f );
}
///////////////////////////////////////////
void
pyrimidine_flip_trial( pose::Pose & pose, 
											 utility::vector1< Size > const & fixed_res_list,
											 scoring::ScoreFunctionOP scorefxn )
{
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::chemical::rna;
	using namespace core::conformation;
	using namespace core::pose;
	using namespace core::chemical;

	Size const total_res = pose.total_residue();
	Pose screen_pose = pose;
	Real orig_score, new_score;
	orig_score = (*scorefxn) (pose);
	new_score = (*scorefxn) (screen_pose);
	std::cout << "Start pyrimidine_flip_trial. Filp residue :";		
	for (Size i = 1; i <= total_res; ++i) {
		if ( check_num_in_vector( i, fixed_res_list ) ) continue;
		Residue const & res = pose.residue(i);
		if ( res.is_RNA() && (res.aa() == na_rcy || res.aa() == na_ura)) {
			Real const orig_chi = pose.torsion( TorsionID( i, id::CHI, 1 ) );
			Real const new_chi = orig_chi + 180.0;
			screen_pose.set_torsion( TorsionID( i, id::CHI, 1 ), new_chi );
			new_score = (*scorefxn) (screen_pose);
			if (new_score < orig_score) { //Flip the chi!
				pose.set_torsion( TorsionID( i, id::CHI, 1 ), new_chi );
				orig_score = new_score;
				std::cout << ' ' << i;		
			} else { //Keep the original chi
				screen_pose.set_torsion( TorsionID( i, id::CHI, 1 ), orig_chi );
			}
		}
	}
	std::cout << std::endl;
}
///////////////////////////////////////////
void
pdb_minimizer() {
	using namespace core::pose;
	using namespace core::conformation;
	using namespace core::chemical;
	using namespace core::kinematics;
	using namespace core::scoring;
	using namespace core::chemical::rna;
	using namespace core::scoring::constraints;
	using namespace core::optimization;
	using namespace core::id;
	using namespace protocols::swa::rna;

	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->
	          									residue_type_set ( RNA );
	bool const vary_bond_geometry_ =  option[ vary_geometry ];
	bool const constrain_phosphate =  option[ constrain_P ];
	bool const ready_set_only_ =  option[ ready_set_only ];
	bool const skip_minimize_ =  option[ skip_minimize ];
	bool const attempt_pyrimidine_flip_ =  option[ attempt_pyrimidine_flip ];
	utility::vector1< core::Size > const fixed_res_list = option[ fixed_res ]();
	utility::vector1< core::Size > const cutpoint_list = option[cutpoint_open]();

	// Read in the pose from pdb.
	Pose pose;
	std::string pdb_name;
	if ( option[ in::file::native ].user() ) {
		import_pose::pose_from_pdb ( pose, *rsd_set, option[in::file::native]() );
		protocols::rna::make_phosphate_nomenclature_matches_mini(pose);
		pdb_name = option[in::file::native]();
	} else {
		utility_exit_with_message("User must specify -native option!");
	}

	std::string output_pdb_name;
	if ( option[ out_pdb ].user() ) {
		output_pdb_name = option[ out_pdb ] ();
	} else {
		output_pdb_name = pdb_name;
		size_t found = output_pdb_name.find(".pdb");
		if (found != std::string::npos) {
			if (ready_set_only_) {
				output_pdb_name.replace(found, found + 4, "_ready_set.pdb");
			} else {
				output_pdb_name.replace(found, found + 4, "_minimize.pdb");
			}			
		} else {
			if (ready_set_only_) {
				output_pdb_name.append("_ready_set.pdb");
			} else {
				output_pdb_name.append("_minimize.pdb");
			}	
		}
	}

	if (ready_set_only_) {
		pose.dump_pdb(output_pdb_name);
		return;
	}
	
	//Setup score function.
	std::string score_weight_file = "rna/rna_hires_elec_dens";
	if ( option[ basic::options::OptionKeys::score::weights ].user() ) {
		score_weight_file= option[ basic::options::OptionKeys::score::weights ]();
		std::cout << "User passed in score:weight option: " << score_weight_file << std::endl;
	}
	core::scoring::ScoreFunctionOP scorefxn =
			ScoreFunctionFactory::create_score_function ( score_weight_file );

	core::scoring::ScoreFunctionOP edens_scorefxn = new ScoreFunction;
	edens_scorefxn -> set_weight( elec_dens_atomwise, 1.0 );

	//Setup fold tree using user input or using Rhiju's function
	if ( cutpoint_list.size() == 0 ) {
		protocols::rna::figure_out_reasonable_rna_fold_tree ( pose );
	} else {
		setup_fold_tree ( pose, cutpoint_list );
	}

	//Add a virtual residue for density scoring
	Size const virtual_res_pos = add_virtual_res ( pose );
	pose::Pose const pose_full = pose;
	Size const nres ( pose.total_residue() );
	Size const nres_moving ( nres - fixed_res_list.size() );

	//Output the sequence
	std::string working_sequence = pose.sequence();
	std::cout << "Pose sequence = " << working_sequence << std::endl;
	protocols::swa::rna::Output_fold_tree_info ( pose.fold_tree(), "rna_pdb_minimizing" );

	//Try flipping the pyrimidines
	if ( attempt_pyrimidine_flip_ ) {
		pyrimidine_flip_trial( pose, fixed_res_list, scorefxn);
	}
	if ( skip_minimize_ ) {
		pose.dump_pdb(output_pdb_name);
		return;	
	}


	//Set the MoveMap, avoiding moving the virtual residue
	std::cout << "Setting up movemap ..." << std::endl;
	kinematics::MoveMap mm;
	ObjexxFCL::FArray1D< bool > allow_insert ( nres, false );
	mm.set_bb ( false );
	mm.set_chi ( false );
	mm.set_jump ( false );

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( pose.residue ( ii ).aa() != core::chemical::aa_vrt ) {
			allow_insert ( ii ) = true;
			mm.set_bb ( ii, true );
			mm.set_chi ( ii, true );
		}
	}

	kinematics::FoldTree fold_tree ( pose.fold_tree() );

	utility::vector1< core::Size > cut_upper, cut_lower;
	for ( Size i = 1; i <= fold_tree.num_jump(); ++i ) {
		Size const k = fold_tree.upstream_jump_residue ( i );
		Size const m = fold_tree.downstream_jump_residue ( i );
		cut_lower.push_back(k);
		cut_upper.push_back(m);

		if ( pose.residue ( k ).aa() != core::chemical::aa_vrt && 
		     pose.residue ( m ).aa() != core::chemical::aa_vrt ) {
			if ( fixed_res_list.size() != 0 &&
					 ! ( check_num_in_vector(k, fixed_res_list) ) &&
					 ! ( check_num_in_vector(m, fixed_res_list) ) ) {
				mm.set_jump ( i, true );
			}
		}
	}


	//Fixed res mode
	if ( fixed_res_list.size() != 0 ) {
		std::cout << "fixed res: ";
	}

	for ( Size i = 1; i <= fixed_res_list.size(); ++i ) {
		Size fixed_res_num ( fixed_res_list[i] );
		std::cout << fixed_res_num << " ";
		Real const coord_sdev ( 0.1 );
		Size const my_anchor ( virtual_res_pos ); //anchor on virtual residue
		ConstraintSetOP cst_set = pose.constraint_set()->clone();
		Residue const & rsd ( pose.residue ( fixed_res_num ) );
		Size const atm_indexP = rsd.atom_index ( "P" );
		Size const atm_indexO3 = rsd.atom_index ( "O3'" );
		Size const atm_indexOP2 = rsd.atom_index ( "OP2" );
		Size const atm_indexC6 = rsd.atom_index ( "C6" );
		Size atm_indexBase = 0;
		if ( rsd.aa() == core::chemical::na_rgu || rsd.aa() == core::chemical::na_rad ) {
			atm_indexBase = rsd.atom_index ( "N9" );
		} else if ( rsd.aa() == core::chemical::na_ura || rsd.aa() == core::chemical::na_rcy ) {
			atm_indexBase = rsd.atom_index ( "N1" );
		}	else {
			utility_exit_with_message("Fixed residue is not a RNA residue!!!!");
		}

		cst_set -> add_constraint ( new CoordinateConstraint ( AtomID ( atm_indexP, fixed_res_num ), AtomID ( 1, my_anchor ), rsd.xyz ( atm_indexP ), new HarmonicFunc ( 0.0, coord_sdev ) ) );
		cst_set -> add_constraint ( new CoordinateConstraint ( AtomID ( atm_indexO3, fixed_res_num ), AtomID ( 1, my_anchor ), rsd.xyz ( atm_indexO3 ), new HarmonicFunc ( 0.0, coord_sdev ) ) );
		cst_set -> add_constraint ( new CoordinateConstraint ( AtomID ( atm_indexBase, fixed_res_num ), AtomID ( 1, my_anchor ), rsd.xyz ( atm_indexBase ), new HarmonicFunc ( 0.0, coord_sdev ) ) );
		cst_set -> add_constraint ( new CoordinateConstraint ( AtomID ( atm_indexC6, fixed_res_num ), AtomID ( 1, my_anchor ), rsd.xyz ( atm_indexC6 ), new HarmonicFunc ( 0.0, coord_sdev ) ) );
		cst_set -> add_constraint ( new CoordinateConstraint ( AtomID ( atm_indexOP2, fixed_res_num ), AtomID ( 1, my_anchor ), rsd.xyz ( atm_indexOP2 ), new HarmonicFunc ( 0.0, coord_sdev ) ) );
		pose.constraint_set ( cst_set );
		scorefxn->set_weight ( coordinate_constraint, 10 );

		mm.set_chi ( fixed_res_num, false );
		mm.set_bb ( fixed_res_num, false );

		allow_insert(fixed_res_num) = false;

		if (fixed_res_num - 1 > 0 && 
		    ! ( check_num_in_vector( fixed_res_num - 1, fixed_res_list ) ) &&
		    ! ( check_num_in_vector( fixed_res_num, cut_lower ) ) ) {
			allow_insert(fixed_res_num) = true;
		}
		if (fixed_res_num + 1 <= nres && 
		    ! ( check_num_in_vector( fixed_res_num + 1, fixed_res_list ) ) &&
		    ! ( check_num_in_vector( fixed_res_num, cut_upper ) ) ) {
			allow_insert(fixed_res_num) = true;
		}

	}
	std::cout << std::endl;

	//constrain phosphate mode
	if (constrain_phosphate) {
		for ( Size i = 1; i <= nres; ++i ) {
			if ( pose.residue ( i ).aa() == core::chemical::aa_vrt ) continue;
			
			bool is_fixed_res = false;
			for ( Size j = 1; j <= fixed_res_list.size(); ++j ) {
				if (i == fixed_res_list[j]) {
					is_fixed_res = true;
					break;
				}
			}

			if (is_fixed_res) continue;

			Real const coord_sdev ( 0.3 );
			Size const my_anchor ( virtual_res_pos ); //anchor on virtual residue
			ConstraintSetOP cst_set = pose.constraint_set()->clone();
			Residue const & rsd ( pose.residue ( i ) );
			Size const atm_indexP = rsd.atom_index ( "P" );
			cst_set -> add_constraint ( new CoordinateConstraint ( AtomID ( atm_indexP, i ), AtomID ( 1, my_anchor ), rsd.xyz ( atm_indexP ), new HarmonicFunc ( 0.0, coord_sdev ) ) );
			pose.constraint_set ( cst_set );
			scorefxn->set_weight ( coordinate_constraint, 10 );
		}
	}

	//Vary Geometry
	if ( vary_bond_geometry_ ) {
		std::cout << "Setup vary_bond_geometry" << std::endl;
		pose::Pose pose_reference;
		create_pose_reference ( pose_full, pose_reference );
		vary_bond_geometry ( mm, pose, pose_reference, allow_insert );
	}

	Output_movemap ( mm, pose );
	scorefxn->show ( std::cout, pose );
	Real const score_before = ( (*scorefxn) (pose) );
	Real const edens_score_before = ( (*edens_scorefxn) (pose) );

	protocols::viewer::add_conformation_viewer ( pose.conformation(), "current", 400, 400 );

	//Start Minimizing the Full Structure
	Pose const start_pose = pose;
	AtomTreeMinimizer minimizer;
	float const dummy_tol ( 0.00000001 );

	std::cout << "Minimize using dfpmin with use_nb_list=true .." << std::endl;
	MinimizerOptions min_options_dfpmin ( "dfpmin", dummy_tol, true, false, false );
	min_options_dfpmin.max_iter ( std::min( 3000, std::max( 1000, int(nres_moving * 12) ) ) );
	minimizer.run ( pose, mm, *scorefxn, min_options_dfpmin );

	scorefxn -> show ( std::cout, pose );
	Real const score = ( (*scorefxn) (pose) );
	Real const edens_score = ( (*edens_scorefxn) (pose) );
	if (score > score_before + 5 || edens_score > edens_score_before * 0.9) {
		std::cout << "current_score = " << score << ", start_score = " << score_before << std::endl;
		std::cout << "current_edens_score = " << edens_score << ", start_edens_score = " << edens_score_before << std::endl;
		std::cout << "The minimization went wild!!! Try alternative minimization using dfpmin with use_nb_list=false .." << std::endl;

		pose = start_pose;

		MinimizerOptions min_options_dfpmin_no_nb ( "dfpmin", dummy_tol, false, false, false );
		min_options_dfpmin_no_nb.max_iter ( std::min( 3000, std::max( 1000, int(nres_moving * 12) ) ) );
		minimizer.run ( pose, mm, *scorefxn, min_options_dfpmin_no_nb );
		scorefxn -> show ( std::cout, pose );
		Real const score = ( (*scorefxn) (pose) );
		Real const edens_score = ( (*edens_scorefxn) (pose) );
		if (score > score_before + 5 || edens_score > edens_score_before * 0.9) {
			std::cout << "current_score = " << score << ", start_score = " << score_before << std::endl;
			std::cout << "current_edens_score = " << edens_score << ", start_edens_score = " << edens_score_before << std::endl;
			pose = start_pose;
			std::cout << "The minimization went wild again!!! Skip the minimization!!!!!" << std::endl;
		}
	}


	pose.dump_pdb ( output_pdb_name );
	std::cout << "Job completed sucessfully." << std::endl;
}
///////////////////////////////////////////////////////////////
void*
my_main ( void* ) {
	pdb_minimizer();
	exit ( 0 );
}
///////////////////////////////////////////////////////////////////////////////
int
main ( int argc, char * argv [] ) {
try {
	utility::vector1< Size > blank_size_vector;
	utility::vector1< std::string > blank_string_vector;
	NEW_OPT ( out_pdb, "name of output pdb file", "" );
	NEW_OPT ( vary_geometry, "vary geometry", false );
	NEW_OPT ( constrain_P, "constrain phosphate", false );
	NEW_OPT ( fixed_res, "optional: residues to be held fixed in minimizer", blank_size_vector );
	NEW_OPT ( cutpoint_open, "optional: chainbreak in full sequence", blank_size_vector );
	NEW_OPT ( ready_set_only, "load in and output directly for reformatting the pdb", false );
	NEW_OPT ( skip_minimize, "output the pdb without minimization", false );
	NEW_OPT ( attempt_pyrimidine_flip, "try to flip pyrimidine by 180 degree and pick the better energy conformer", false );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	core::init::init ( argc, argv );
	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////
	protocols::viewer::viewer_main ( my_main );
} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
}
}
