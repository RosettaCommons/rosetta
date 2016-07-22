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
/// @details
/// @author Rhiju Das (rhiju@stanford.edu), Parin Sripakdeevong (sripakpa@stanford.edu)


#include <protocols/farna/movers/RNA_Minimizer.hh>
#include <protocols/farna/libraries/RNA_ChunkLibrary.hh> // for ROSETTA_LIBRARY_DOMAIN
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <protocols/stepwise/modeler/rna/util.hh> //Parin Sripakdeevong
#include <protocols/stepwise/modeler/util.hh> // for figuring out moving chainbreaks
#include <protocols/farna/util.hh>
#include <protocols/farna/movers/RNA_LoopCloser.hh>
#include <protocols/farna/options/RNA_MinimizerOptions.hh>
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

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>

using namespace core;
using basic::T;

static THREAD_LOCAL basic::Tracer TR( "protocols.farna.movers.RNA_Minimizer" );

namespace protocols {
namespace farna {
namespace movers {

RNA_Minimizer::RNA_Minimizer( options::RNA_MinimizerOptionsCOP options /* = 0 */ ):
	Mover(),
	options_( options ),
	coord_sdev_( 10.0 * std::sqrt(10.0) ), // awkward, but matches an old setting.
	coord_cst_weight_( 1.0 ),
	perform_minimizer_run_( true ),
	include_default_linear_chainbreak_( true  ),
	close_loops_( true ),
	min_type_( "lbfgs_armijo_nonmonotone" ) //Parin S. Jan 12, 2012
{
	Mover::type("RNA_Minimizer");
}

/// @details  Apply the RNA full atom minimizer.
///
void RNA_Minimizer::apply( core::pose::Pose & pose )
{

	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( scorefxn_ == 0 ) scorefxn_ = get_rna_hires_scorefxn(); //->clone();

	time_t pdb_start_time = time(NULL);
	scoring::constraints::ConstraintSetOP save_pose_constraints = pose.constraint_set()->clone();
	if ( pose.constraint_set()->has_constraints() ) {
		if ( !scorefxn_->has_nonzero_weight( atom_pair_constraint ) )  scorefxn_->set_weight( atom_pair_constraint, 1.0 );
		if ( !scorefxn_->has_nonzero_weight( base_pair_constraint ) )  scorefxn_->set_weight( base_pair_constraint, 1.0 );
		if ( !scorefxn_->has_nonzero_weight( coordinate_constraint ) ) scorefxn_->set_weight( coordinate_constraint, 1.0 );
	}
	if ( options_ == 0 ) options_ = options::RNA_MinimizerOptionsOP( new options::RNA_MinimizerOptions );
	if ( options_->vary_bond_geometry() ) scorefxn_->set_weight( rna_bond_geometry, 1.0 );
	if ( include_default_linear_chainbreak_ && !scorefxn_->has_nonzero_weight( linear_chainbreak ) ) scorefxn_->set_weight( linear_chainbreak, 5.0 );

	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions minimizer_options( min_type_, dummy_tol, use_nblist, options_->deriv_check(), options_->deriv_check() );
	minimizer_options.nblist_auto_update( true );

	/////////////////////////////////////////////////////
	kinematics::MoveMap mm;
	if ( atom_level_domain_map_input_ != 0 ) {
		atom_level_domain_map_ = atom_level_domain_map_input_->clone();
	} else {
		atom_level_domain_map_ = toolbox::AtomLevelDomainMapOP( new toolbox::AtomLevelDomainMap( pose ) ); // initialized to let all dofs move.
	}
	if ( options_->minimize_bps() ) update_atom_level_domain_map_to_move_rosetta_library_chunks();
	update_atom_level_domain_map_with_extra_minimize_res( pose );
	setup_movemap( mm, pose );

	scoring::constraints::ConstraintSetOP pose_constraints_without_coordinate_tethers = pose.constraint_set()->clone();
	if ( options_->minimizer_use_coordinate_constraints() ) core::scoring::constraints::add_coordinate_constraints( pose, coord_sdev_ );
	scoring::constraints::ConstraintSetOP pose_constraints_with_coordinate_tethers = pose.constraint_set()->clone();

	utility::vector1< Size > moving_chainbreaks = stepwise::modeler::figure_out_moving_chain_break_res( pose, mm );
	RNA_LoopCloser rna_loop_closer;

	Real const fa_rep_final( scorefxn_->get_weight( fa_rep ) );

	for ( Size r = 1; r <= options_->minimize_rounds(); r++ ) {

		ScoreFunctionOP minimize_scorefxn_ = scorefxn_->clone();

		Real const suppress = static_cast<Real>(r)/options_->minimize_rounds();
		minimize_scorefxn_->set_weight( fa_rep, fa_rep_final * suppress  );

		if ( !options_->skip_o2prime_trials() ) o2prime_trials( pose, minimize_scorefxn_ );

		//Prevent explosions on first minimize.
		// this is silly. in first round, just make sure coordinate_constraint is one, and use constraints 'supplemented' with coordinate constraints.
		if ( options_->minimizer_use_coordinate_constraints() ) {
			if ( r == 1 ) {
				if ( !minimize_scorefxn_->has_nonzero_weight( coordinate_constraint ) ) minimize_scorefxn_->set_weight( coordinate_constraint, coord_cst_weight_ );
				pose.constraint_set( pose_constraints_with_coordinate_tethers );
			} else {
				pose.constraint_set( pose_constraints_without_coordinate_tethers );
			}
		}
		TR << "Minimizing...round= " << r << std::endl;

		if ( perform_minimizer_run_ ) minimizer.run( pose, mm, *minimize_scorefxn_, minimizer_options );

		if ( close_loops_ )  rna_loop_closer.apply( pose, moving_chainbreaks );

	}

	time_t pdb_end_time = time(NULL);

	pose.constraint_set( save_pose_constraints );

	(*scorefxn_)( pose );
	scorefxn_->show( TR, pose );

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
	output <<
		"Deriv check:              " << (options_->deriv_check() ? "True" : "False")  <<
		"\nSkip o2prime trials:       " << (options_->skip_o2prime_trials() ? "True" : "False") <<
		"\nPerform minimizer run:    " << (perform_minimizer_run_ ? "True" : "False") <<
		"\nVary bond geometry:       " << (options_->vary_bond_geometry() ? "True" : "False") <<
		"\nDump pdb:                 " << (options_->dump_pdb() ? "True" : "False") <<
		"\nMove first rigid body:    " << (options_->move_first_rigid_body() ? "True" : "False") <<
		"\nMin type:                 " << min_type_ <<
		"\nScore function:           " << scorefxn_->get_name() <<
		"\nUse coordinate constraints:        " << (options_->minimizer_use_coordinate_constraints() ? "True" : "False")  <<
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

	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( !pose.residue(i).is_RNA() ) continue;

		task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
		task->nonconst_residue_task(i).or_include_current( true );

	}

	TR << "Orienting 2' hydroxyls..." << std::endl;

	pack::rotamer_trials( pose, *packer_scorefxn_, task);

}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// should unify with StepWiseMinimizer -- nice function to get movemap is available inside atom_level_domain_map now.
void
RNA_Minimizer::setup_movemap( kinematics::MoveMap & mm, pose::Pose & pose ) {

	using namespace core::id;
	using namespace core::chemical::rna;

	Size const nres( pose.total_residue() );

	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( false );

	// torsions
	for  ( Size i = 1; i <= nres; i++ )  {

		if ( !pose.residue(i).is_RNA() ) continue;

		for ( Size j = 1; j <= ( NUM_RNA_MAINCHAIN_TORSIONS + pose.residue(i).type().nchi() ); j++ ) {

			id::TorsionID rna_torsion_id( i, id::BB, j );
			if ( j > NUM_RNA_MAINCHAIN_TORSIONS ) rna_torsion_id = id::TorsionID( i, id::CHI, j - NUM_RNA_MAINCHAIN_TORSIONS );
			if ( !atom_level_domain_map_->get( rna_torsion_id, pose.conformation() ) ) continue;

			// this is not general. Sigh:
			if ( pose.residue(i).has_variant_type( core::chemical::VIRTUAL_PHOSPHATE ) &&
					( j == 1 || j == 2 || j == 3 ) ) {
				continue;
			}
			mm.set( rna_torsion_id, true );

		}
	}

	// jumps
	for ( Size n = 1; n <= pose.fold_tree().num_jump(); n++ ) {
		Size const jump_pos1( pose.fold_tree().upstream_jump_residue( n ) );
		std::string const jump_atom1( pose.fold_tree().upstream_atom( n ) );
		AtomID jump_atom_id1( 1, jump_pos1 );
		if ( jump_atom1.size() > 0 ) jump_atom_id1 = named_atom_id_to_atom_id( NamedAtomID( jump_atom1, jump_pos1 ), pose );

		Size const jump_pos2( pose.fold_tree().downstream_jump_residue( n ) );
		std::string const jump_atom2( pose.fold_tree().downstream_atom( n ) );
		AtomID jump_atom_id2( 1, jump_pos2 );
		if ( jump_atom2.size() > 0 ) jump_atom_id2 = named_atom_id_to_atom_id( NamedAtomID( jump_atom2, jump_pos2 ), pose );

		if ( moveable_jump( jump_atom_id1, jump_atom_id2, *atom_level_domain_map_ ) )  mm.set_jump( n, true );
	}

	// allow rigid body movements... check for virtual residue at end and at least two chunks with jumps to it.
	protocols::farna::let_rigid_body_jumps_move( mm, pose, options_->move_first_rigid_body() );

	for ( Size n = 1; n <= options_->extra_minimize_chi_res().size(); n++ ) mm.set( id::TorsionID( options_->extra_minimize_chi_res()[n], id::CHI, 1 ), true );

	// vary bond geometry
	if ( options_->vary_bond_geometry() ) simple_moves::setup_vary_rna_bond_geometry( mm, pose, atom_level_domain_map_ );

}

/////////////////////////////////////////////////////////////////////////////////
void
RNA_Minimizer::update_atom_level_domain_map_to_move_rosetta_library_chunks()
{
	atom_level_domain_map_->update_to_move_chunks_with_domain( libraries::ROSETTA_LIBRARY_DOMAIN );
}


/////////////////////////////////////////////////////////////////////////////////
void
RNA_Minimizer::update_atom_level_domain_map_with_extra_minimize_res( pose::Pose const & pose ){

	if ( options_->extra_minimize_res().size() == 0 ) return;

	utility::vector1< id::AtomID > atom_ids_to_move;

	for ( Size n = 1; n <= options_->extra_minimize_res().size(); n++ ) {

		Size const i = options_->extra_minimize_res()[n];
		runtime_assert( pose.residue( i ).is_RNA() );

		for ( Size j = 1; j <= pose.residue(i).natoms(); j++ ) {
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

	for ( Size n = 1; n <= atom_ids_to_move.size(); n++ ) {
		if ( atom_level_domain_map_->has_domain( atom_ids_to_move[n] ) ) atom_level_domain_map_->set( atom_ids_to_move[n],  true );
	}
	// We should do a double check that we're not introducing movement that would mess up a domain?

}

/////////////////////////////////////////////////////////////////////////////////
void
RNA_Minimizer::set_score_function( core::scoring::ScoreFunctionCOP scorefxn ){
	scorefxn_ = scorefxn->clone();
}

/////////////////////////////////////////////////////////////////////////////////
std::ostream &operator<< ( std::ostream &os, RNA_Minimizer const &mover )
{
	mover.show(os);
	return os;
}

} //movers
} //farna
} //protocols

