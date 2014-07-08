// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/CyclizationMover.hh
/// @brief Implimentation file for CyclizationMover
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

// unit headers
#include <protocols/simple_moves/CyclizationMover.hh>

// protocols headers
#include <protocols/simple_moves/MinMover.hh>

// core headers
#include <core/types.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>

#include <core/scoring/func/HarmonicFunc.hh>

#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>

#include <core/kinematics/MoveMap.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// numeric
#include <numeric/constants.hh>

// c++ headers
#include <string>
#include <algorithm>

namespace protocols {
namespace simple_moves {

static basic::Tracer TR("protocols.simple_moves.CyclizationMover");

/// @brief Default constructor
CyclizationMover::CyclizationMover( core::Size chain_to_cyclize, bool add_constraints = true, bool minimize = true, core::Size minimization_rebuild_rounds = 3 ) :
  protocols::moves::Mover("CyclizationMover"),
  chain_to_cyclize_( chain_to_cyclize ),
	nterm_rsd_num_( 0 ),
	cterm_rsd_num_( 0 ),
	add_constraints_( add_constraints ),
	minimize_( minimize ),
	minimization_rebuild_rounds_( minimization_rebuild_rounds ),
	score_fxn_( 0 ),
	move_map_( 0 )
{}

CyclizationMover::CyclizationMover( core::Size chain_to_cyclize, bool add_constraints, bool minimize, core::Size minimization_rebuild_rounds, core::scoring::ScoreFunctionOP score_fxn, core::kinematics::MoveMapOP move_map ) :
chain_to_cyclize_( chain_to_cyclize ),
	nterm_rsd_num_( 0 ),
	cterm_rsd_num_( 0 ),
	add_constraints_( add_constraints ),
	minimize_( minimize ),
	minimization_rebuild_rounds_( minimization_rebuild_rounds ),
	score_fxn_( score_fxn ),
  move_map_( move_map )
{}

/// @brief Sets up inter residue cyclic connections and potentially adds constraints, and minimizes the pose
void
CyclizationMover::apply( core::pose::Pose & pose )
{
	using namespace core;

	// check to see if specified chain exists
	runtime_assert( chain_to_cyclize_ <= pose.conformation().num_chains() );

	// get residue numbers of chain N-terminus and C-terminus
	nterm_rsd_num_ = pose.conformation().chain_begin( chain_to_cyclize_ );
	cterm_rsd_num_ = pose.conformation().chain_end( chain_to_cyclize_ );

	// setup residue connections and constraints
	setup_connections( pose );

	// add constraints to maintain the cyclization to the pose
	if( add_constraints_ ) {
		setup_constraints ( pose );
	}

	// minimize the pose to bring
	if ( minimize_ && minimization_rebuild_rounds_ > 0 ) {
		setup_scorefunction();
		setup_minimizer( pose );
		for ( Size i( 1 ); i <= minimization_rebuild_rounds_; ++i ) {
			minimize_rebuild( pose );
		}
	}
}

/// @brief Modifes terminal ResidueTypes to cyclized variants and then uses connformation::detect_bonds() to create a residue connection
void
CyclizationMover::setup_connections( core::pose::Pose & pose )
{
	using namespace core;
	using namespace chemical;

  TR << "Setting up residue connections..." << std::endl;

	// make sure the seqence positions have been initialized to meaningful values
	runtime_assert( nterm_rsd_num_ != 0 && cterm_rsd_num_ != 0 );

	// make sure the two residuetypes are peptide or peptoid as that is all the NtermConnect and CtermConnect patches apply to
	runtime_assert( pose.residue( nterm_rsd_num_ ).type().is_protein() || pose.residue( nterm_rsd_num_ ).type().is_peptoid() );
	runtime_assert( pose.residue( cterm_rsd_num_ ).type().is_protein() || pose.residue( cterm_rsd_num_ ).type().is_peptoid() );

	// get types name for N-terminus and C-terminus (manipulating strings like this is a little hacky)
	std::string nterm_connect_type_name( pose.residue( nterm_rsd_num_ ).type().name3() + "_p:NtermConnect" );
	std::string cterm_connect_type_name( pose.residue( cterm_rsd_num_ ).type().name3() + "_p:CtermConnect" );

	// remove spaces if the name3 really only has 2 letters, Damn it Tim!
	nterm_connect_type_name.erase( std::remove( nterm_connect_type_name.begin(), nterm_connect_type_name.end(), ' ' ), nterm_connect_type_name.end() );
	cterm_connect_type_name.erase( std::remove( cterm_connect_type_name.begin(), cterm_connect_type_name.end(), ' ' ), cterm_connect_type_name.end() );

	// get CtermConnect and NtermConnect variant types
	ResidueTypeSetCAP rsd_type_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	ResidueType const & nterm_connect_type( rsd_type_set->name_map( nterm_connect_type_name ) );
	ResidueType const & cterm_connect_type( rsd_type_set->name_map( cterm_connect_type_name ) );

	// replace N-terminus and C-terminus with connect varients aligning new residue to old ones
	pose::replace_pose_residue_copying_existing_coordinates( pose, nterm_rsd_num_, nterm_connect_type );
	pose::replace_pose_residue_copying_existing_coordinates( pose, cterm_rsd_num_, cterm_connect_type );

	// tell rosetta to detect the new bonds
	pose.conformation().detect_bonds();
	pose.conformation().rebuild_residue_connection_dependent_atoms( nterm_rsd_num_, 2 );
	pose.conformation().rebuild_residue_connection_dependent_atoms( cterm_rsd_num_, 2 );

	// TODO: replace these calls with something more specfic
	//pose.conformation().show_residue_connections();
  //pose.fold_tree();
}

/// @brief Creates constraints to maintain the proper conformation and adds it to the pose.
void
CyclizationMover::setup_constraints( core::pose::Pose & pose )
{
	// TODO: make the constraints funcs based on mm_strech, mm_bend, mm_twist

	using namespace core;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring::constraints;
	using namespace scoring::func;

  TR << "Setting up constraints to maintain cycle based on polymeric base types..." << std::endl;

	// get connect variants of ResidueTypes of N-terminus and C-terminus
	ResidueType const & nterm_connect_type(	pose.residue( nterm_rsd_num_ ).type() );
	ResidueType const & cterm_connect_type(	pose.residue( cterm_rsd_num_ ).type() );

	// get base variants of ResidueTypes of N-terminus and C-terminus
	ResidueTypeSetCAP rsd_type_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	// remove spaces if the name3 really only has 2 letters, Damn it Tim!
	std::string nterm_base_type_name( nterm_connect_type.name3() );
	std::string cterm_base_type_name( cterm_connect_type.name3() );
	nterm_base_type_name.erase( std::remove( nterm_base_type_name.begin(), nterm_base_type_name.end(), ' ' ), nterm_base_type_name.end() );
  cterm_base_type_name.erase( std::remove( cterm_base_type_name.begin(), cterm_base_type_name.end(), ' ' ), cterm_base_type_name.end() );

	ResidueType const & nterm_base_type( rsd_type_set->name_map( nterm_base_type_name ) );
	ResidueType const & cterm_base_type( rsd_type_set->name_map( cterm_base_type_name ) );

	// make sure the two residuetypes are peptide or peptoid as that is all the NtermConnect and CtermConnect patches apply to
	runtime_assert( nterm_connect_type.is_protein() || nterm_connect_type.is_peptoid() );
	runtime_assert( cterm_connect_type.is_protein() || cterm_connect_type.is_peptoid() );

	// get ideal internal coordinates of nterm lower connect and cterm upper connect
	AtomICoor nterm_base_lc_icoor( nterm_base_type.lower_connect().icoor() );
	AtomICoor cterm_base_uc_icoor( cterm_base_type.upper_connect().icoor() );

	// TR << "DEBUG NTERM LOWER ICOOR:\t" << nterm_base_type.lower_connect_id() << "\t"
	// << nterm_base_lc_icoor.d()     << "\t" << nterm_base_lc_icoor.stub_atom1().atomno()  << "\t"
	// << nterm_base_lc_icoor.theta() << "\t" << nterm_base_lc_icoor.stub_atom2().atomno()  << "\t"
	// << nterm_base_lc_icoor.phi()   << "\t" << nterm_base_lc_icoor.stub_atom3().atomno()  << "\t"
	// << std::endl;

	// TR << "DEBUG CTERM UPPER ICOOR:\t" << cterm_base_type.upper_connect_id() << "\t"
	// << cterm_base_uc_icoor.d()     << "\t" << cterm_base_uc_icoor.stub_atom1().atomno()  << "\t"
	// << cterm_base_uc_icoor.theta() << "\t" << cterm_base_uc_icoor.stub_atom2().atomno()  << "\t"
	// << cterm_base_uc_icoor.phi()   << "\t" << cterm_base_uc_icoor.stub_atom3().atomno()  << "\t"
	// << std::endl;

	// create AtomIDs for the two central atoms and the atoms one bond away from them
	id::AtomID nterm_n(  nterm_base_lc_icoor.stub_atom1().atomno(), nterm_rsd_num_ );
	id::AtomID nterm_ca( nterm_base_lc_icoor.stub_atom2().atomno(), nterm_rsd_num_ );
	id::AtomID cterm_c(  cterm_base_uc_icoor.stub_atom1().atomno(), cterm_rsd_num_ );
	id::AtomID cterm_ca( cterm_base_uc_icoor.stub_atom2().atomno(), cterm_rsd_num_ );

	// create an AtomPairConstraint between the two central atoms
	TR << "Nterm atom " << nterm_base_lc_icoor.stub_atom1().atomno() << " is " << nterm_base_lc_icoor.d() << " from its lower connect" << std::endl;
	TR << "Cterm atom " << cterm_base_uc_icoor.stub_atom1().atomno() << " is " << cterm_base_uc_icoor.d() << " from its upper connect" << std::endl;
	Real ap_cst_length( ( nterm_base_lc_icoor.d() + cterm_base_uc_icoor.d() ) / 2 );
	TR << "Adding AtomPairConstraint of length:" << ap_cst_length << std::endl;
	ConstraintOP b1( new AtomPairConstraint( nterm_n, cterm_c, new HarmonicFunc( ap_cst_length, 0.01 ) ) );
	pose.add_constraint( b1 );

	// create two AngleConstraints
	TR << "Nterm atom " << nterm_base_lc_icoor.stub_atom2().atomno() << " and " << nterm_base_lc_icoor.stub_atom1().atomno() << " make and angle of " << nterm_base_lc_icoor.theta() << "with its lower connect" << std::endl;
	Real a1_cst_radian( numeric::constants::r::pi - nterm_base_lc_icoor.theta() );
	TR << "Adding AngleConstraint of radian: " << a1_cst_radian << std::endl;
	ConstraintOP a1( new AngleConstraint( nterm_ca, nterm_n, cterm_c, new HarmonicFunc( a1_cst_radian, 0.1 ) ) );
	pose.add_constraint( a1 );

	TR << "Cterm atom " << cterm_base_uc_icoor.stub_atom2().atomno() << " and " << cterm_base_uc_icoor.stub_atom1().atomno() << " make and angle of " << cterm_base_uc_icoor.theta() << "with its upper connect" << std::endl;
	Real a2_cst_radian( numeric::constants::r::pi - cterm_base_uc_icoor.theta() );
	TR << "Adding AngleConstraint of radian: " << a2_cst_radian << std::endl;
	ConstraintOP a2( new AngleConstraint( nterm_n, cterm_c, cterm_ca, new HarmonicFunc( a2_cst_radian, 0.1 ) ) );
	pose.add_constraint( a2 );

}

/// @brief Setup the score function with the appropriate weights on the constraints.
/// If user has provided a score function the and set the weights to a specific value,
/// that value is kept. If weight is not set, prints a warning.
void
CyclizationMover::setup_scorefunction()
{
	using namespace core;
	using namespace scoring;

	if ( score_fxn_ == 0 ) {
		TR << "Creating score function and setting geometric constraint weights to 1" << std::endl;
		score_fxn_ = get_score_function();
		score_fxn_->set_weight(	atom_pair_constraint, 10 );
		score_fxn_->set_weight(	angle_constraint, 10 );
		score_fxn_->set_weight(	dihedral_constraint, 10 );
	} else if ( score_fxn_->get_weight( atom_pair_constraint ) == 0 ) {
		TR << "WARNING: atom_pair_constraint weight set to zero. Cyclization constraints will not work properly." << std::endl;
	} else if ( score_fxn_->get_weight( angle_constraint ) == 0 ) {
		TR << "WARNING: angle_constraint weight set to zero. Cyclization constraints will not work properly" << std::endl;
	} else if ( score_fxn_->get_weight( dihedral_constraint ) == 0 ) {
		TR << "WARNING: dihedral_constraint weight set to zero. Cyclization constraints will not work properly" << std::endl;
	}
}


/// @brief
void
CyclizationMover::setup_minimizer( core::pose::Pose & pose )
{
	using namespace core;
	using namespace kinematics;

	if ( move_map_ == 0 ) {
		TR << "Creating move map and setting all backbone and sidechain DOFs to movable" << std::endl;
		move_map_ = new MoveMap();
		move_map_->set_bb_true_range( pose.conformation().chain_begin( chain_to_cyclize_ ), pose.conformation().chain_end( chain_to_cyclize_ ) );
		move_map_->set_chi_true_range( pose.conformation().chain_begin( chain_to_cyclize_ ), pose.conformation().chain_end( chain_to_cyclize_ ) );
	}
}

/// @brief Minimize the pose using the score function, and constraints we setup
void
CyclizationMover::minimize_rebuild( core::pose::Pose & pose )
{
	using namespace core;
	using namespace protocols;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	simple_moves::MinMoverOP min_mover( new simple_moves::MinMover( move_map_, score_fxn_, option[ run::min_type ].value(), 0.01, true ) );

	min_mover->apply( pose );

	// rebuild the conection dependant atoms
	pose.conformation().rebuild_residue_connection_dependent_atoms( nterm_rsd_num_, 2 );
	pose.conformation().rebuild_residue_connection_dependent_atoms( cterm_rsd_num_, 2 );
}

}//simple_moves
}//protocols
