// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/PeptideStapleMover.cc
/// @brief Peptide Staples are covalent links between i/i+[4,7] residues in a helix. Applying a PeptideStapleMover creates one of these links.
/// @author Jacob Corn and Andrew Leaver-Fay

// Unit Headers
#include <protocols/simple_moves/PeptideStapleMover.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <core/id/AtomID.hh>

#include <core/kinematics/MoveMap.hh>

//#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
//#include <basic/options/after_opts.hh>
//#include <basic/options/util.hh>

#include <core/pose/Pose.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/optimization/MinimizerOptions.hh>

//#include <core/pack/task/operation/TaskOperation.hh>
//#include <core/pack/task/TaskFactory.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

#include <protocols/simple_moves/MinMover.hh>
//#include <protocols/simple_moves/PackRotamersMover.hh>

#include <numeric/constants.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
using basic::T;
using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "protocols.simple_moves.PeptideStapleMover" );

// C++ Headers
#include <sstream>

#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>


// ObjexxFCL Headers

namespace protocols {
namespace simple_moves {

PeptideStapleMover::PeptideStapleMover( core::Size const staple_start, core::Size const staple_gap=4 ) :
	protocols::moves::Mover(),
	seqpos_(staple_start),
	staple_gap_(staple_gap)
{
	protocols::moves::Mover::type( "PeptideStapleMover" );
}

/// @brief Places and minimizes a PeptideStaple (i+4 or i+7, depending on staple_gap) at seqpos in pose
void PeptideStapleMover::apply( core::pose::Pose & pose )
{

	// BUNCH OF SANITY CHECKS
	// test whether we're running off the pose
	if ( (seqpos_ + staple_gap_) > pose.total_residue()  )
	{
		TR << "A staple gap of " << staple_gap_ << " runs off the end of the pose. Aborting staple insertion at residue " << seqpos_ << std::endl;
		return;
	}
	// test connection points for terminii
	if ( pose.residue( seqpos_ ).is_terminus() || pose.residue( seqpos_+staple_gap_ ).is_terminus() )
	{
		TR << "Staple insertion is not supported at chain terminii. Aborting staple insertion at residue " << seqpos_ << std::endl;
		return;
	}
	// test residues between connections for secondary structure, jumps, and terminii
	for( Size i = seqpos_; i <= seqpos_ + staple_gap_; ++i ) {
		if( pose.secstruct( i ) != 'H' ) {
			TR << "Secondary structure at residue " << i << " is " << pose.secstruct(i) << ", but stapling along non-helix residues is untested!" << std::endl;
		}
		if( pose.fold_tree().is_jump_point( i ) ) {
			TR << "Peptide stapling across jumps is untested!" << std::endl;
		}
	}

	core::chemical::ResidueTypeSet const & residue_set ( pose.residue(1).residue_type_set() );
	core::chemical::ResidueType const * stapleA_type = NULL;
	core::chemical::ResidueType const * stapleB_type = NULL;

	if ( staple_gap_ == 4 ) {
		stapleA_type = &residue_set.name_map("STAPLE08A") ;
		stapleB_type = &residue_set.name_map("STAPLE08B") ;
	}

	// create staple residues on top of
	core::conformation::Residue stapleA_res ( *stapleA_type, pose.residue(seqpos_), pose.conformation() );
	core::conformation::Residue stapleB_res ( *stapleB_type, pose.residue(seqpos_ + staple_gap_), pose.conformation() );

	pose.replace_residue( seqpos_, stapleA_res, true );
	pose.replace_residue( seqpos_ + staple_gap_, stapleB_res, true);
	pose.conformation().detect_bonds();

	derive_staple_constraints_( pose );
	minimize_( pose );
}

std::string
PeptideStapleMover::get_name() const {
	return "PeptideStapleMover";
}

/// @brief Derives the constraints to be used for minimizing the stapled residue
/// @author Andrew Leaver-Fay
void PeptideStapleMover::derive_staple_constraints_( core::pose::Pose & pose )
{
	typedef core::Size Size;
	typedef core::Real Real;
	using namespace core::scoring::constraints;
	using namespace core::chemical;
	using namespace core::id;
	//ConstraintSetOP cst_set( new ConstraintSet() );

	TR << "residue " << seqpos_ << " " << pose.residue( seqpos_ ).name() << std::endl;

	//	if ( pose.residue( ii ).name() == "STAPLE08A" ) {
	Size const seqpos_conn_atom = pose.residue( seqpos_ ).type().residue_connection( 3 ).atomno();
	Size const seqpos_vc_atom = pose.residue( seqpos_ ).atom_index( "VC" );

	Size const jj = pose.residue( seqpos_ ).connect_map( 3 ).resid();
	Size const jj_conn_atom = pose.residue( jj ).type().residue_connection( 3 ).atomno();
	Size const jj_vc_atom = pose.residue( jj ).atom_index( "VC" );

	TR << "seqpos conn: " << pose.residue( seqpos_ ).atom_name( seqpos_conn_atom );
	TR << " seqpos vc: " << pose.residue( seqpos_ ).atom_name( seqpos_vc_atom );
	TR << " jj conn: " << pose.residue( jj ).atom_name( jj_conn_atom );
	TR << " jj vc: " << pose.residue( jj ).atom_name( jj_vc_atom );

	{
	core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, 0.001 ) );
	ConstraintOP apc1( new AtomPairConstraint(
		AtomID( seqpos_conn_atom, seqpos_ ),
		AtomID( jj_vc_atom, jj ),
		fx ) );
	pose.add_constraint( apc1 );
	}

	{
	core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, 0.001 ) );
	ConstraintOP apc2( new AtomPairConstraint(
		AtomID( jj_conn_atom, jj ),
		AtomID( seqpos_vc_atom, seqpos_ ),
		fx ) );
	pose.add_constraint( apc2 );
	}

	Size const seqpos_conn_atom_base = pose.residue( seqpos_ ).type().icoor( seqpos_conn_atom ).stub_atom( 1 ).atomno();
	Size const jj_conn_atom_base = pose.residue( jj ).type().icoor( jj_conn_atom ).stub_atom( 1 ).atomno();

	Real const seqpos_ideal_angle = numeric::constants::d::pi - pose.residue( seqpos_ ).type().icoor( seqpos_vc_atom ).theta();
	Real const jj_ideal_angle = numeric::constants::d::pi - pose.residue( jj ).type().icoor( jj_vc_atom ).theta();

	TR << "seqpos_conn_atom_base: " << seqpos_conn_atom_base << " " << pose.residue( seqpos_ ).atom_name( seqpos_conn_atom_base ) << " " << seqpos_ideal_angle << std::endl;
	TR << "jj_conn_atom_base: " << jj_conn_atom_base << " " << pose.residue( jj ).atom_name( jj_conn_atom_base ) << " " << jj_ideal_angle << std::endl;

	{
	core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( jj_ideal_angle, 0.01 ) );
	ConstraintOP bac1( new AngleConstraint(
		AtomID( seqpos_conn_atom, seqpos_),
		AtomID( seqpos_vc_atom, seqpos_ ),
		AtomID( jj_conn_atom_base, jj ),
		fx ) );
	pose.add_constraint( bac1 );
	}

	{
	core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( seqpos_ideal_angle, 0.01 ) );
	ConstraintOP bac2( new AngleConstraint(
		AtomID( seqpos_conn_atom_base, seqpos_),
		AtomID( jj_vc_atom, jj ),
		AtomID( jj_conn_atom, jj ),
		fx ) );
	pose.add_constraint( bac2 );
	}

	Real cross_connection_dihedral_val (0);
/*	if ( basic::options::option[ cross_connection_dihedral ].user() ) {
		cross_connection_dihedral_val = basic::options::option[ cross_connection_dihedral ]();
	}
*/

	{
	core::scoring::func::FuncOP fx( new core::scoring::func::CircularHarmonicFunc( cross_connection_dihedral_val * numeric::constants::d::pi / 180.0, 0.1 ) );
	ConstraintOP dcst( new DihedralConstraint(
		AtomID( seqpos_conn_atom_base, seqpos_ ),
		AtomID( seqpos_conn_atom, seqpos_ ),
		AtomID( jj_conn_atom, jj ),
		AtomID( jj_conn_atom_base, jj ),
		fx ) );
	pose.add_constraint( dcst );
	}

	utility::vector1< AtomIndices > const & seqpos_chi( pose.residue( seqpos_ ).chi_atoms() );
	for ( Size kk = 1; kk <= seqpos_chi.size(); ++kk ) {
		Real ideal_dihedral = pose.residue( seqpos_ ).type().icoor( seqpos_chi[ kk ][ 4 ] ).phi();
		TR << "seqpos ideal dihedral: " << ideal_dihedral << " " << ideal_dihedral * 180 / numeric::constants::d::pi << std::endl;
		core::scoring::func::FuncOP fx( new core::scoring::func::CircularHarmonicFunc( ideal_dihedral, 0.1 ) );
		ConstraintOP dcst( new DihedralConstraint(
			AtomID( seqpos_chi[ kk ][ 1 ], seqpos_ ),
			AtomID( seqpos_chi[ kk ][ 2 ], seqpos_ ),
			AtomID( seqpos_chi[ kk ][ 3 ], seqpos_ ),
			AtomID( seqpos_chi[ kk ][ 4 ], seqpos_ ),
			fx ) );
		pose.add_constraint( dcst );
	}

	utility::vector1< AtomIndices > const & jj_chi( pose.residue( jj ).chi_atoms() );
	for ( Size kk = 1; kk <= jj_chi.size(); ++kk ) {
		Real ideal_dihedral = pose.residue( jj ).type().icoor( jj_chi[ kk ][ 4 ] ).phi();
		TR << "jj ideal dihedral: " << ideal_dihedral << " " << ideal_dihedral * 180 / numeric::constants::d::pi << std::endl;
		core::scoring::func::FuncOP fx( new core::scoring::func::CircularHarmonicFunc( ideal_dihedral, 0.1 ) );
		ConstraintOP dcst( new DihedralConstraint(
			AtomID( jj_chi[ kk ][ 1 ], jj ),
			AtomID( jj_chi[ kk ][ 2 ], jj ),
			AtomID( jj_chi[ kk ][ 3 ], jj ),
			AtomID( jj_chi[ kk ][ 4 ], jj ),
			fx ) );
		pose.add_constraint( dcst );
	}
	//pose.add_constraint( cst_set );
}


/// @brief Minimizes a placed PeptideStaple using the constraints derived with derive_stape_constraints_
/// @author Andrew Leaver-Fay
void PeptideStapleMover::minimize_( core::pose::Pose & pose )
{
	using namespace core::scoring;
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	scorefxn->set_weight( atom_pair_constraint, 1 );
	scorefxn->set_weight( angle_constraint, 1 );
	scorefxn->set_weight( dihedral_constraint, 1 );
	scorefxn->set_weight( coordinate_constraint, 1 );

	TR << "Total score before staple minimization: " << (*scorefxn)( pose ) << std::endl;

	// the movable dof's
	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_chi( seqpos_, true );
	mm->set_chi( seqpos_ + staple_gap_, true);
	protocols::simple_moves::MinMover min_mover( mm, scorefxn, "dfpmin_armijo_nonmonotone_atol", 10.0 /*tolerance*/,
		true /*use_nblist*/, false /*deriv_check*/, false /* non verbose-deriv-check, default*/ );
	min_mover.min_options()->nblist_auto_update( true );
	min_mover.apply( pose );

	TR << "Total score after staple minimization: " << (*scorefxn)( pose ) << std::endl;
}

} // moves
} // protocols
