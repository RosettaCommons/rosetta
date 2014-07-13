// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief protocols that are specific to RNA_Minimizer
/// @detailed
/// @author Rhiju Das (rhiju@stanford.edu), Parin Sripakdeevong (sripakpa@stanford.edu)


#include <protocols/farna/RNA_Minimizer.hh>
#include <protocols/toolbox/AllowInsert.hh>
#include <protocols/stepwise/sampling/rna/util.hh> //Parin Sripakdeevong
#include <protocols/stepwise/sampling/util.hh> // for figuring out moving chainbreaks
#include <protocols/farna/util.hh>
#include <protocols/farna/RNA_LoopCloser.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <protocols/simple_moves/ConstrainToIdealMover.hh>
#include <basic/options/option.hh>

//Minimizer stuff
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

//Packer stuff
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

//C++ headers
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#ifdef WIN32
#include <ctime>
#endif

// option key includes
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>

using namespace core;
using basic::T;

static basic::Tracer TR( "protocols.rna.RNA_Minimizer" ) ;

namespace protocols {
namespace farna {

RNA_Minimizer::RNA_Minimizer():
    Mover(),
    deriv_check_( false ),
    use_coordinate_constraints_( true ),
    coord_sdev_( 10.0 * std::sqrt(10.0) ), // awkward, but matches an old setting.
    coord_cst_weight_( 1.0 ),
    rounds_( basic::options::option[ basic::options::OptionKeys::rna::minimize_rounds ] ),
    skip_o2prime_trials_( false ),
    perform_minimizer_run_( true ),
    vary_bond_geometry_( false ),
    include_default_linear_chainbreak_( true  ),
    do_dump_pdb_( false ),
    move_first_rigid_body_( false ),
		close_loops_( true ),
    min_type_( "dfpmin" ) //Parin S. Jan 12, 2012
{
    Mover::type("RNA_Minimizer");
    if ( basic::options::option[ basic::options::OptionKeys::score::weights ].user() ) {
        scorefxn_ = core::scoring::get_score_function();
    } else {
        scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::RNA_HIRES_WTS );
    }
}

/// @details  Apply the RNA full atom minimizer.
///
void RNA_Minimizer::apply( core::pose::Pose & pose	)
{

	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	time_t pdb_start_time = time(NULL);
	scoring::constraints::ConstraintSetOP save_pose_constraints = pose.constraint_set()->clone();
	if (pose.constraint_set()->has_constraints() ){
		if ( !scorefxn_->has_nonzero_weight( atom_pair_constraint ) )  scorefxn_->set_weight( atom_pair_constraint, 1.0 );
		if ( !scorefxn_->has_nonzero_weight( coordinate_constraint ) ) scorefxn_->set_weight( coordinate_constraint, 1.0 );
	}
	if ( vary_bond_geometry_ ) scorefxn_->set_weight( rna_bond_geometry, 1.0 );
	if( include_default_linear_chainbreak_ && !scorefxn_->has_nonzero_weight( linear_chainbreak ) )	scorefxn_->set_weight( linear_chainbreak, 5.0 );

	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions options( min_type_, dummy_tol, use_nblist, deriv_check_, deriv_check_ );
	options.nblist_auto_update( true );

	/////////////////////////////////////////////////////
	kinematics::MoveMap mm;
	if (!allow_insert_) allow_insert_ = new toolbox::AllowInsert( pose ); // initialized to let all dofs move.
	update_allow_insert_with_extra_minimize_res( pose );
 	setup_movemap( mm, pose );

	scoring::constraints::ConstraintSetOP pose_constraints_without_coordinate_tethers = pose.constraint_set()->clone();
	if (use_coordinate_constraints_) core::scoring::constraints::add_coordinate_constraints( pose, coord_sdev_ );
	scoring::constraints::ConstraintSetOP pose_constraints_with_coordinate_tethers = pose.constraint_set()->clone();

	utility::vector1< Size > moving_chainbreaks = stepwise::sampling::figure_out_moving_chain_break_res( pose, mm );
	RNA_LoopCloser rna_loop_closer;

	Real const fa_rep_final( scorefxn_->get_weight( fa_rep ) );

	for (Size r = 1; r <= rounds_; r++ ) {

		ScoreFunctionOP minimize_scorefxn_ = scorefxn_->clone();

		Real const suppress = static_cast<Real>(r)/rounds_;
		minimize_scorefxn_->set_weight( fa_rep, fa_rep_final * suppress  );

		if (!skip_o2prime_trials_) o2prime_trials( pose, minimize_scorefxn_ );

		//Prevent explosions on first minimize.
		// this is silly. in first round, just make sure coordinate_constraint is one, and use constraints 'supplemented' with coordinate constraints.
		if ( use_coordinate_constraints_) {
			if ( r == 1 ) {
				if ( !minimize_scorefxn_->has_nonzero_weight( coordinate_constraint ) ) minimize_scorefxn_->set_weight( coordinate_constraint, coord_cst_weight_ );
				pose.constraint_set( pose_constraints_with_coordinate_tethers );
			} else {
				pose.constraint_set( pose_constraints_without_coordinate_tethers );
			}
		}
		TR << "Minimizing...round= " << r << std::endl;

		if (perform_minimizer_run_) minimizer.run( pose, mm, *minimize_scorefxn_, options );

		if (close_loops_)		rna_loop_closer.apply( pose, moving_chainbreaks );

	}

	time_t pdb_end_time = time(NULL);

	pose.constraint_set( save_pose_constraints );

	(*scorefxn_)( pose );
	scorefxn_->show( std::cout, pose );

	TR << "RNA minimizer finished in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;

}


///////////////////////////////////////////////////////////////////////////////
	////////////////////////
std::string
RNA_Minimizer::get_name() const {
	return "RNA_Minimizer";
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_Minimizer::show(std::ostream & output) const
{
	Mover::show(output);
	output <<   "Deriv check:              " << (deriv_check_ ? "True" : "False")  <<
				"\nSkip o2prime trials:       " << (skip_o2prime_trials_ ? "True" : "False") <<
				"\nPerform minimizer run:    " << (perform_minimizer_run_ ? "True" : "False") <<
				"\nVary bond geometry:       " << (vary_bond_geometry_ ? "True" : "False") <<
				"\nDump pdb:                 " << (do_dump_pdb_ ? "True" : "False") <<
				"\nMove first rigid body:    " << (move_first_rigid_body_ ? "True" : "False") <<
				"\nMin type:                 " << min_type_ <<
				"\nScore function:           " << scorefxn_->get_name() <<
				"\nUse coordinate constraints:        " << (use_coordinate_constraints_ ? "True" : "False")  <<
				"\nInclude default linear chainbreak: " << (include_default_linear_chainbreak_ ? "True" : "False");

}


///////////////////////////////////////////////////////////////////////////////
// Make this its own Mover?
void
RNA_Minimizer::o2prime_trials(
  core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP const & packer_scorefxn_ ) const
{

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	//task->initialize_from_command_line(); //Jan 20, 2012 Testing.

	for (Size i = 1; i <= pose.total_residue(); i++) {
		if ( !pose.residue(i).is_RNA() ) continue;

		task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
		task->nonconst_residue_task(i).or_include_current( true );

	}

	TR << "Orienting 2' hydroxyls..." << std::endl;

	pack::rotamer_trials( pose, *packer_scorefxn_, task);

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// should unify with StepWiseMinimizer -- nice function to get movemap is available inside allow_insert now.
void
RNA_Minimizer::setup_movemap( kinematics::MoveMap & mm, pose::Pose & pose ) {

	using namespace core::id;
	using namespace core::chemical::rna;

	Size const nres( pose.total_residue() );

	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( false );

	// torsions
	for  (Size i = 1; i <= nres; i++ )  {

		if ( !pose.residue(i).is_RNA() ) continue;

		for (Size j = 1; j <= ( NUM_RNA_MAINCHAIN_TORSIONS + pose.residue(i).type().nchi() ); j++) {

			id::TorsionID rna_torsion_id( i, id::BB, j );
			if ( j > NUM_RNA_MAINCHAIN_TORSIONS) rna_torsion_id = id::TorsionID( i, id::CHI, j - NUM_RNA_MAINCHAIN_TORSIONS );
			if ( !allow_insert_->get( rna_torsion_id, pose.conformation() ) ) continue;

			// this is not general. Sigh:
			if ( pose.residue(i).has_variant_type("VIRTUAL_PHOSPHATE") && ( j == 1 || j == 2 || j == 3 ) ) continue;
			mm.set( rna_torsion_id, true );

		}
	}

	// jumps
	for (Size n = 1; n <= pose.fold_tree().num_jump(); n++ ){
		Size const jump_pos1( pose.fold_tree().upstream_jump_residue( n ) );
		std::string const jump_atom1( pose.fold_tree().upstream_atom( n ) );
		AtomID jump_atom_id1( 1, jump_pos1 );
		if (jump_atom1.size() > 0) jump_atom_id1 = named_atom_id_to_atom_id( NamedAtomID( jump_atom1, jump_pos1 ), pose );

		Size const jump_pos2( pose.fold_tree().downstream_jump_residue( n ) );
		std::string const jump_atom2( pose.fold_tree().downstream_atom( n ) );
		AtomID jump_atom_id2( 1, jump_pos2 );
		if (jump_atom2.size() > 0) jump_atom_id2 = named_atom_id_to_atom_id( NamedAtomID( jump_atom2, jump_pos2 ), pose );

		if ( moveable_jump( jump_atom_id1, jump_atom_id2, *allow_insert_ ) ) 	mm.set_jump( n, true );
	}

	// allow rigid body movements... check for virtual residue at end and at least two chunks with jumps to it.
	protocols::farna::let_rigid_body_jumps_move( mm, pose, move_first_rigid_body_ );

	for ( Size n = 1; n <= extra_minimize_chi_res_.size(); n++ )	mm.set( id::TorsionID( extra_minimize_chi_res_[n], id::CHI, 1 ), true );

	// vary bond geometry
	if ( vary_bond_geometry_ ) {
		// Let additional degrees of freedom vary -- but apply constraints to stay near
		// ideal bond lengths and angles!
		protocols::simple_moves::ConstrainToIdealMover CTIMover;
		core::kinematics::MoveMapOP mmop(mm.clone());
		CTIMover.set_movemap(mmop);
		CTIMover.set_AllowInsert(allow_insert_);
		CTIMover.apply(pose);
		mm = (*CTIMover.get_movemap());
	}

}

/////////////////////////////////////////////////////////////////////////////////
void
RNA_Minimizer::update_allow_insert_with_extra_minimize_res( pose::Pose const & pose ){

	if ( extra_minimize_res_.size() == 0 ) return;

	utility::vector1< id::AtomID > atom_ids_to_move;

	for ( Size n = 1; n <= extra_minimize_res_.size(); n++ ){

		Size const i = extra_minimize_res_[n];
		runtime_assert( pose.residue( i ).is_RNA() );

		for ( Size j = 1; j <= pose.residue(i).natoms(); j++ ){
			if ( pose.residue(i).is_virtual( j ) ) continue;
			atom_ids_to_move.push_back( id::AtomID( j, i ) );
		}

		if ( pose.fold_tree().is_cutpoint( i ) ) continue;
		if ( i > pose.total_residue() ) continue;
		if ( !pose.residue( i+1 ).is_RNA() ) continue;

		// go ahead and minimize backbone torsions up to next pucker.
		atom_ids_to_move.push_back( named_atom_id_to_atom_id( id::NamedAtomID( " OP2", i+1 ), pose ) );
		atom_ids_to_move.push_back( named_atom_id_to_atom_id( id::NamedAtomID( " OP1", i+1 ), pose ) );
		atom_ids_to_move.push_back( named_atom_id_to_atom_id( id::NamedAtomID( " P  ", i+1 ), pose ) );
		atom_ids_to_move.push_back( named_atom_id_to_atom_id( id::NamedAtomID( " O5'", i+1 ), pose ) );

	}

	for ( Size n = 1; n <= atom_ids_to_move.size(); n++ ){
		if ( allow_insert_->has_domain( atom_ids_to_move[n] ) ) allow_insert_->set( atom_ids_to_move[n],  true );
	}
	// We should do a double check that we're not introducing movement that would mess up a domain?

}

/////////////////////////////////////////////////////////////////////////////////
void
RNA_Minimizer::set_score_function( core::scoring::ScoreFunctionOP const & scorefxn ){
	scorefxn_ = scorefxn->clone();
}

/////////////////////////////////////////////////////////////////////////////////
void
RNA_Minimizer::set_allow_insert(toolbox::AllowInsertOP allow_insert ){
	allow_insert_ = allow_insert->clone();
}

/////////////////////////////////////////////////////////////////////////////////
void
RNA_Minimizer::set_extra_minimize_res( utility::vector1< core::Size > setting ){
	extra_minimize_res_ = setting;
}

/////////////////////////////////////////////////////////////////////////////////
void
RNA_Minimizer::set_extra_minimize_chi_res( utility::vector1< core::Size > setting ){
	extra_minimize_chi_res_ = setting;
}

/////////////////////////////////////////////////////////////////////////////////
std::ostream &operator<< ( std::ostream &os, RNA_Minimizer const &mover )
{
	mover.show(os);
	return os;
}

} //farna
} //protocols

