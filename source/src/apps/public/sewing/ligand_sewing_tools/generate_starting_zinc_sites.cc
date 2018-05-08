// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/guffysl/generate_starting_zinc_sites.cc
/// @brief Builds single helices containing two histidines at i and i+4 as the beginning of a zinc binding site
/// @author Sharon Guffy

//Devel
#include <devel/init.hh>
//Protocols
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
//Core
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/rms_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
//Basic
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
//Utility
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
//Numeric
#include <numeric/random/random.hh>
#include <numeric/constants.hh>
//tracers
using basic::Error;
using basic::Warning;
static basic::Tracer TR("apps.pilot.guffysl.generate_starting_zinc_sites");

//This app will not be in MPI--it should be fast and simple enough to run without it. Will change if this proves not to be true.

namespace local {
//Local options
basic::options::IntegerOptionKey const window_width( "window_width" );
basic::options::IntegerOptionKey const helix_length( "helix_length" );
basic::options::StringOptionKey const reference_zinc_site( "reference_zinc_site" );
basic::options::IntegerOptionKey const first_ref_his( "first_ref_his" );
basic::options::IntegerOptionKey const ref_zinc( "ref_zinc" );
basic::options::RealOptionKey const temperature( "temperature" );
basic::options::IntegerOptionKey const n_moves( "n_moves" );
basic::options::RealOptionKey const constraint_weight( "constraint_weight" );
basic::options::RealOptionKey const max_angle( "max_angle" );
basic::options::IntegerOptionKey const small_moves( "small_moves" );
basic::options::IntegerOptionKey const shear_moves( "shear_moves" );
}

void
mutate_residue( core::pose::Pose & pose, core::Size res_position, std::string target_res_name ){
	core::chemical::ResidueTypeSetCOP restype_set( pose.conformation().residue_type_set_for_conf() );
	core::conformation::ResidueOP new_res = core::conformation::ResidueFactory::create_residue( restype_set->name_map( target_res_name ), pose.residue( res_position ), pose.conformation() );
	core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( res_position ), *new_res, pose.conformation(), true );
	pose.replace_residue( res_position, *new_res, false );
	pose.conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only( res_position );
}
int main( int argc, char* argv[] ){
	try{
		//Set up options
		basic::options::option.add( local::window_width, "Minimum required number of residues before and after zinc site" ).def( 6 );
		basic::options::option.add( local::helix_length, "Length of helix in residues" ).def( 28 );
		basic::options::option.add( local::reference_zinc_site, "Source of zinc position" ).def( "2GZJ_ca_split_1.pdb" );
		basic::options::option.add( local::first_ref_his, "Residue number of first histidine in reference zinc site" ).def( 4 );
		basic::options::option.add( local::ref_zinc, "Residue number of zinc in reference zinc site" ).def( 9 );
		basic::options::option.add( local::temperature, "Temperature for Monte Carlo object" ).def( 1.0 );
		basic::options::option.add( local::n_moves, "Number of backbone move trials (i.e.sets of moves )to perform per structure" ).def( 5 );
		basic::options::option.add( local::constraint_weight, "Weight to use for zinc geometry constraints" ).def( 1 );
		basic::options::option.add( local::max_angle, "Maximum angle for small/shear moves" ).def( 7 );
		basic::options::option.add( local::small_moves, "Number of small moves to make per trial" ).def( 3 );
		basic::options::option.add( local::shear_moves, "Number of shear moves to make per trial" ).def( 5 );
		devel::init( argc, argv );
		//Initialize variables from options
		core::Size window_width = core::Size( basic::options::option[ local::window_width ].value() );
		core::Size helix_length = core::Size( basic::options::option[ local::helix_length ].value() );
		core::Real constraint_weight = basic::options::option[ local::constraint_weight ].value();
		core::Real max_angle = basic::options::option[ local::max_angle ].value();
		core::Size small_moves = core::Size( basic::options::option[ local::small_moves ].value() );
		core::Size shear_moves = core::Size( basic::options::option[ local::shear_moves ].value() );
		std::string sequence = "";
		while ( sequence.size() < helix_length ) {
			sequence += "A";
		}
		std::string reference_zinc_site_name = basic::options::option[ local::reference_zinc_site ].value();
		core::Size first_ref_his = basic::options::option[ local::first_ref_his ].value();
		core::Size zinc = basic::options::option[ local::ref_zinc ].value();
		core::Real temperature = basic::options::option[ local::temperature ].value();

		//Import zinc site that will be use for retrieving the zinc ion and the histidine chi angles
		//Default site name is "2GZJ_ca_split_1.pdb" with his at 4 and 8
		core::pose::PoseOP reference_zinc_site = core::import_pose::pose_from_file( reference_zinc_site_name );

		core::pose::Pose ref_pose;
		core::pose::make_pose_from_sequence( ref_pose, sequence, core::chemical::FA_STANDARD );
		//Now make this pose into a helix
		for ( core::Size i = 1; i <= ref_pose.total_residue(); ++i ) {
			TR << "Setting phi and psi of residue " << i << std::endl;
			ref_pose.set_phi( i, -60.0);
			ref_pose.set_psi( i, -40.0 );
			ref_pose.set_omega( i, 180 );
		}
		//This initial helix will be used as the starting point for all of our poses.

		//We're going to minimize the zinc position and chi angles given strict distance/angle constraints
		core::conformation::Residue sample_zn_res( reference_zinc_site->residue( zinc ) ); //We'll need this later
		//Create score function
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		//Set weights for distance, angle, and dihedral constraints
		scorefxn->set_weight( core::scoring::atom_pair_constraint, constraint_weight );
		scorefxn->set_weight( core::scoring::angle_constraint, constraint_weight );
		scorefxn->set_weight( core::scoring::dihedral_constraint, constraint_weight );

		//Create functions and constants for constraints
		core::Real const ZN_BOND_LENGTH = 2.2;
		core::Real const ZN_BOND_ANGLE = 109.5 * numeric::constants::r::deg2rad;
		core::Real const HIS_BOND_ANGLE = 126.5 * numeric::constants::r::deg2rad;
		core::Real const PLANAR_DIHEDRAL = 180 * numeric::constants::r::deg2rad; //dihedral about histidine ring
		//Make the functions I'll use for these constraints
		core::scoring::func::HarmonicFuncOP distance_func( new core::scoring::func::HarmonicFunc( ZN_BOND_LENGTH, 0.05 ) );
		core::scoring::func::CircularHarmonicFuncOP zn_angle_func( new core::scoring::func::CircularHarmonicFunc( ZN_BOND_ANGLE, 5 * numeric::constants::r::deg2rad ) );
		core::scoring::func::CircularHarmonicFuncOP his_angle_func( new core::scoring::func::CircularHarmonicFunc( HIS_BOND_ANGLE, 5 * numeric::constants::r::deg2rad ) );
		core::scoring::func::CircularHarmonicFuncOP dihedral_func( new core::scoring::func::CircularHarmonicFunc( PLANAR_DIHEDRAL, 5 * numeric::constants::r::deg2rad ) );

		//This minmover will optimize the zinc position before we refine the helix
		core::kinematics::MoveMapOP zn_mm( new core::kinematics::MoveMap );
		zn_mm->set_bb( false );
		zn_mm->set_chi( true );
		zn_mm->set_jump( true );
		protocols::minimization_packing::MinMover zn_min( zn_mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.01, true );

		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
		movemap->set_bb( true );
		movemap->set_chi( false );
		movemap->set_jump( false );

		protocols::minimization_packing::MinMover initial_min( movemap, scorefxn, "lbfgs_armijo_nonmonotone" ,0.001, true );
		//Create movers that we'll use to sample the pose
		//This is shamelessly copied from the PyRosetta tutorial refinement protocol
		core::pose::Pose pose = ref_pose;
		//Make sure the constraints are in place like they should be--this should be unnecessary, but I'm leaving it here to try if I get weird things
		//pose.add_constraints( ref_pose->constraint_set()->get_all_constraints() );
		//Create movemap
		//Small mover
		protocols::simple_moves::SmallMoverOP small_mover( new protocols::simple_moves::SmallMover( movemap, temperature, small_moves ) );
		//This is copied, but it's also the maximum angle (approximately) by which helices vary from ideal phi/psi
		small_mover->angle_max( max_angle );
		//Shear mover
		protocols::simple_moves::ShearMoverOP shear_mover( new protocols::simple_moves::ShearMover( movemap, temperature, shear_moves ) );
		shear_mover->angle_max( max_angle );
		//MinMover
		protocols::minimization_packing::MinMoverOP min_mover( new protocols::minimization_packing::MinMover( movemap, scorefxn, "lbfgs_armijo_nonmonotone", 0.01, true ) );
		//I think I might actually do a RandomMover instead of a sequence mover to get higher acceptance rates--I'll change to that if I don't like the results I'm seeing
		protocols::moves::SequenceMoverOP sequence_mover( new protocols::moves::SequenceMover( small_mover, shear_mover, min_mover ) );
		//MonteCarlo object
		protocols::moves::MonteCarloOP monte_carlo( new protocols::moves::MonteCarlo( pose, *scorefxn, temperature ) );
		//TrialMover
		protocols::moves::TrialMoverOP trial_mover( new protocols::moves::TrialMover( sequence_mover, monte_carlo ) );
		//RepeatMover
		protocols::moves::RepeatMover repeat_mover( trial_mover, core::Size( basic::options::option[ local::n_moves ].value() ) );
		//Iterate over nstruct
		for ( core::Size i = 1; i <= core::Size( basic::options::option[ basic::options::OptionKeys::out::nstruct ].value() ); ++i ) {
			pose = ref_pose;
			//Mutate two residues at i and i+4 (where (1 + window ) < i < ( length of helix - window ) ) to histidines with proper chi angles (as measured in native zinc sites)
			//Randomly choose the first histidine's position
			core::Size first_his = numeric::random::random_range( window_width + 1, helix_length - window_width - 4 );
			mutate_residue( pose, first_his, "HIS_D" );
			mutate_residue( pose, first_his + 4, "HIS_D" );
			//Set chi angles
			//His has 2 chi values
			//First his
			for ( core::Size i = 1; i <= pose.residue( first_his ).nchi(); ++ i ) {
				pose.set_chi( i, first_his, reference_zinc_site->residue( first_ref_his ).chi( i ) );
			}
			//Second his
			for ( core::Size i = 1; i <= pose.residue( first_his + 4 ).nchi(); ++ i ) {
				pose.set_chi( i, first_his + 4, reference_zinc_site->residue( first_ref_his + 4 ).chi( i ) );
			}

			//Align reference_zinc_site to pose at the two histidines
			//Generate atom map: we'll want to align on all atoms of those two histidines
			core::id::AtomID_Map< core::id::AtomID > atom_map;
			core::pose::initialize_atomid_map( atom_map, *reference_zinc_site, core::id::GLOBAL_BOGUS_ATOM_ID );
			//We'll only set values for those two residues
			for ( core::Size his1_atom = 1; his1_atom <= reference_zinc_site->residue( first_ref_his ).nheavyatoms(); ++his1_atom ) {
				core::id::AtomID const id1( his1_atom, first_ref_his );
				core::id::AtomID const id2( his1_atom, first_his );
				atom_map.set( id1, id2 );
			}
			for ( core::Size his2_atom = 1; his2_atom <= reference_zinc_site->residue( first_ref_his + 4 ).nheavyatoms(); ++his2_atom ) {
				core::id::AtomID const id1( his2_atom, first_ref_his + 4);
				core::id::AtomID const id2( his2_atom, first_his + 4);
				atom_map.set( id1, id2 );
			}
			core::scoring::superimpose_pose( *reference_zinc_site, pose, atom_map );

			//Append zinc ion: create a new one at the proper XYZ, and then append_residue_by_jump to the first histidine.
			pose.append_residue_by_jump( sample_zn_res, first_his, "NE2", "ZN", false );

			//Now we'll add constraints to the pose
			//Distance constraints
			core::id::AtomID his1_ne2( pose.residue( first_his ).atom_index( "NE2" ), first_his );
			core::id::AtomID his1_cd2( pose.residue( first_his ).atom_index( "CD2" ), first_his );
			core::id::AtomID his1_ce1( pose.residue( first_his ).atom_index( "CE1" ), first_his );
			core::id::AtomID his1_cg( pose.residue( first_his ).atom_index( "CG" ), first_his );
			core::id::AtomID his1_nd1( pose.residue( first_his ).atom_index( "ND1" ), first_his );
			core::id::AtomID zn_atom( 1, helix_length + 1 );
			core::id::AtomID his2_ne2( pose.residue( first_his + 4 ).atom_index( "NE2" ), first_his + 4 );
			core::id::AtomID his2_cd2( pose.residue( first_his + 4 ).atom_index( "CD2" ), first_his + 4 );
			core::id::AtomID his2_ce1( pose.residue( first_his + 4 ).atom_index( "CE1" ), first_his + 4 );
			core::id::AtomID his2_cg( pose.residue( first_his + 4 ).atom_index( "CG" ), first_his + 4);
			core::id::AtomID his2_nd1( pose.residue( first_his + 4 ).atom_index( "ND1" ), first_his + 4 );
			core::scoring::constraints::AtomPairConstraintOP his1_distance_cst( new core::scoring::constraints::AtomPairConstraint( his1_ne2, zn_atom, distance_func ) );
			pose.add_constraint( his1_distance_cst );
			core::scoring::constraints::AtomPairConstraintOP his2_distance_cst( new core::scoring::constraints::AtomPairConstraint( his2_ne2, zn_atom, distance_func ) );
			pose.add_constraint( his2_distance_cst );
			//Angle constraints
			core::scoring::constraints::AngleConstraintOP zn_angle_cst( new core::scoring::constraints::AngleConstraint( his1_ne2, zn_atom, his2_ne2, zn_angle_func ) );
			pose.add_constraint( zn_angle_cst );
			core::scoring::constraints::AngleConstraintOP his1_angle1_cst( new core::scoring::constraints::AngleConstraint(  his1_cd2, his1_ne2, zn_atom, his_angle_func ) );
			pose.add_constraint( his1_angle1_cst );
			core::scoring::constraints::AngleConstraintOP his1_angle2_cst( new core::scoring::constraints::AngleConstraint(  his1_ce1, his1_ne2, zn_atom, his_angle_func ) );
			pose.add_constraint( his1_angle2_cst );
			core::scoring::constraints::AngleConstraintOP his2_angle1_cst( new core::scoring::constraints::AngleConstraint(  his2_cd2, his2_ne2, zn_atom, his_angle_func ) );
			pose.add_constraint( his2_angle1_cst );
			core::scoring::constraints::AngleConstraintOP his2_angle2_cst( new core::scoring::constraints::AngleConstraint(  his2_ce1, his2_ne2, zn_atom, his_angle_func ) );
			pose.add_constraint( his2_angle2_cst );
			//Dihedral constraints
			core::scoring::constraints::DihedralConstraintOP his1_dihedral1_cst( new core::scoring::constraints::DihedralConstraint( his1_cg, his1_cd2, his1_ne2, zn_atom, dihedral_func ) );
			pose.add_constraint( his1_dihedral1_cst );
			core::scoring::constraints::DihedralConstraintOP his1_dihedral2_cst( new core::scoring::constraints::DihedralConstraint( his1_nd1, his1_ce1, his1_ne2, zn_atom, dihedral_func ) );
			pose.add_constraint( his1_dihedral2_cst );
			core::scoring::constraints::DihedralConstraintOP his2_dihedral1_cst( new core::scoring::constraints::DihedralConstraint( his2_cg, his2_cd2, his2_ne2, zn_atom, dihedral_func ) );
			pose.add_constraint( his2_dihedral1_cst );
			core::scoring::constraints::DihedralConstraintOP his2_dihedral2_cst( new core::scoring::constraints::DihedralConstraint( his2_nd1, his2_ce1, his2_ne2, zn_atom, dihedral_func ) );
			pose.add_constraint( his2_dihedral2_cst );

			//NOW, that's over with
			//Minimize the pose
			//We're going to minimize ONLY the chi angles and jumps
			//ref_pose.dump_pdb( "before_min.pdb" );
			zn_min.apply( pose );
			initial_min.apply( pose );
			//Hopefully that didn't just break everything
			//ref_pose.dump_pdb( "starting_zinc_site.pdb" );



			//Now we're ready to sample
			//Reset monte_carlo
			monte_carlo->reset( pose );
			//Apply movers
			repeat_mover.apply( pose );
			monte_carlo->recover_low( pose );
			//Dump pdb
			pose.dump_pdb( "zinc_site_" + utility::to_string( first_his ) + "_" + utility::to_string( i ) + ".pdb" );
		}
		return 0;
	}
catch ( utility::excn::Exception const &e ){
	std::cout << "Caught exception " << e.msg() << std::endl;
	return -1;
}
}
