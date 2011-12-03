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
/// @author Rhiju Das


#include <protocols/rna/RNA_Minimizer.hh>
#include <protocols/toolbox/AllowInsert.hh>
// AUTO-REMOVED #include <protocols/rna/RNA_ProtocolUtil.hh>
// AUTO-REMOVED #include <protocols/moves/MonteCarlo.fwd.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
// AUTO-REMOVED #include <core/kinematics/AtomTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/pose/annotated_sequence.hh>
#include <core/scoring/rna/RNA_Util.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_FittedTorsionInfo.hh>

// AUTO-REMOVED #include <protocols/viewer/viewers.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/CoordinateConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>

#include <protocols/simple_moves/ConstrainToIdealMover.hh>

#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>

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

// AUTO-REMOVED #include <numeric/conversions.hh>

// External library headers


//C++ headers
#include <vector>
#include <string>
#include <sstream>
// AUTO-REMOVED #include <fstream>
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

static basic::Tracer TR( "protocols.rna.rna_minimizer" ) ;

namespace protocols {
namespace rna {

RNA_Minimizer::RNA_Minimizer():
  Mover(),
	deriv_check_( false ),
  use_coordinate_constraints_( true ),
	skip_o2star_trials_( false ),
	vary_bond_geometry_( false )
{
	Mover::type("RNA_Minimizer");
	if ( basic::options::option[ basic::options::OptionKeys::score::weights ].user() ) {
		scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( basic::options::option[ basic::options::OptionKeys::score::weights ]() );
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

	TR << "Check it! SEQUENCE " << pose.sequence() << std::endl;

	time_t pdb_start_time = time(NULL);

	//	protocols::viewer::add_conformation_viewer( pose.conformation(), "minimizer", 400, 400 );

	if (pose.constraint_set()->has_constraints() ){
		scorefxn_->set_weight( atom_pair_constraint, 1.0 );
		scorefxn_->set_weight( angle_constraint, 1.0 );
	}

	if ( vary_bond_geometry_ ){
		scorefxn_->set_weight( rna_bond_geometry, 1.0 );
	}

	// New as of 2011 (based on results with minimizing in SWA work
	scorefxn_->set_weight( linear_chainbreak, 5.0 );

	(*scorefxn_)(pose);
	//scorefxn_->show( std::cout, pose );

	/////////////////////////////////////////////////////
	//Need to be careful with coordinate constraints -- remove them later???
	/////////////////////////////////////////////////////
	scoring::constraints::ConstraintSetOP save_pose_constraints = pose.constraint_set()->clone();
	if (use_coordinate_constraints_) core::scoring::constraints::add_coordinate_constraints( pose );

	/////////////////////////////////////////////////////
	// minimzer
	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, deriv_check_, deriv_check_ );
	//	MinimizerOptions options( "dfpmin_armijo_nonmonotone", dummy_tol, use_nblist, deriv_check_, deriv_check_ );
	options.nblist_auto_update( true );

	kinematics::MoveMap mm;

	if (!allow_insert_) allow_insert_ = new toolbox::AllowInsert( pose ); // initialized to let all dofs move.
	setup_movemap( mm, pose );

// 	mm.set_bb( false );
// 	mm.set_chi( false );
// 	mm.set_jump( false );
// 	//id::TorsionID torsion_id( 1, id::CHI, 1 );
// 	id::TorsionID torsion_id( 1, id::BB, 3 );
// 	id::DOF_ID dof_id( pose.conformation().dof_id_from_torsion_id( torsion_id) );
// 	mm.set( dof_id, true );

	//Could/should be private data.
	Real const coord_cst_weight = 0.1;
	Size const rounds( option[ basic::options::OptionKeys::rna::minimize_rounds ] );
	Real const fa_rep_final( scorefxn_->get_weight( fa_rep ) );

	(*scorefxn_)( pose );
	scorefxn_->show( std::cout, pose );

	for (Size r = 1; r <= rounds; r++ ) {

		ScoreFunctionOP minimize_scorefxn_ = scorefxn_->clone();

		Real const suppress = static_cast<Real>(r)/rounds;
		minimize_scorefxn_->set_weight( fa_rep, fa_rep_final * suppress  );

		//		pose.dump_pdb( "before_o2star_trials.pdb" );
		if (!skip_o2star_trials_) o2star_trials( pose, minimize_scorefxn_ );
		//		pose.dump_pdb( "after_o2star_trials.pdb" );

		//scorefxn_->show( std::cout, pose );

		//Prevent explosions on first minimize.
		if (r == 1 && use_coordinate_constraints_) minimize_scorefxn_->set_weight( coordinate_constraint, coord_cst_weight );
		TR << "Minimizing..." << std::endl;

		minimizer.run( pose, mm, *minimize_scorefxn_, options );

		//		io::pdb::dump_pdb( pose, "minimize_round"+string_of( r)+".pdb" );
		//minimize_scorefxn_->show( std::cout, pose );

	}

	time_t pdb_end_time = time(NULL);

	//scorefxn_->show( std::cout, pose );

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


///////////////////////////////////////////////////////////////////////////////
// Make this its own Mover?
void
RNA_Minimizer::o2star_trials(
  core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP const & packer_scorefxn_,
	bool const do_pack_instead_of_rotamer_trials /* = false */ ) const
{

	//TR << "Repacking 2'-OH ... " << std::endl;

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();

	for (Size i = 1; i <= pose.total_residue(); i++) {
		if ( !pose.residue(i).is_RNA() ) continue;
		task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
		//		task->nonconst_residue_task(i).or_ex4( true );
		task->nonconst_residue_task(i).or_include_current( true );
		// How about bump check?
	}

	TR << "Orienting 2' hydroxyls..." << std::endl;

	if (do_pack_instead_of_rotamer_trials ){
		pack::pack_rotamers( pose, *packer_scorefxn_, task);
	} else {
		pack::rotamer_trials( pose, *packer_scorefxn_, task);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_Minimizer::setup_movemap( kinematics::MoveMap & mm, pose::Pose & pose ) {

	using namespace core::id;
	using namespace core::scoring::rna;

	Size const nres( pose.total_residue() );

	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( false );

	for  (Size i = 1; i <= nres; i++ )  {

			for (Size j = 1; j <= NUM_RNA_TORSIONS; j++) {
				id::TorsionID rna_torsion_id( i, id::BB, j );
				if ( j > NUM_RNA_MAINCHAIN_TORSIONS) rna_torsion_id = id::TorsionID( i, id::CHI, j - NUM_RNA_MAINCHAIN_TORSIONS );

				// CALL ALLOW INSERT TO SEE IF TORSION IS MOVEABLE [don't copy code! create allow_insert->get( TorsionID ). ]
				if ( !allow_insert_->get( rna_torsion_id, pose.conformation() ) ) continue;

				DOF_ID dof_id( pose.conformation().dof_id_from_torsion_id( rna_torsion_id ) );
				mm.set( dof_id, true );

			}

			// Following is deprecated, now that allow_insert can report on individual atoms or torsions.
			//		if ( !allow_insert_->get( i ) ) continue;
			//		mm.set_bb(  i, true );
			//		mm.set_chi( i, true );

			//Chi (glycosidic bond to base)
			//mm.set( pose.conformation().dof_id_from_torsion_id( id::TorsionID( i, id::CHI, 1 ) ), true );
			//2'-hydroxyl
			//		mm.set( pose.conformation().dof_id_from_torsion_id( id::TorsionID( i, id::CHI, 4 ) ), true );

	}

	for (Size n = 1; n <= pose.fold_tree().num_jump(); n++ ){
		Size const jump_pos1( pose.fold_tree().upstream_jump_residue( n ) );
		std::string const jump_atom1( pose.fold_tree().upstream_atom( n ) );
		AtomID const jump_atom_id1( named_atom_id_to_atom_id( NamedAtomID( jump_atom1, jump_pos1 ), pose ) );

		Size const jump_pos2( pose.fold_tree().downstream_jump_residue( n ) );
		std::string const jump_atom2( pose.fold_tree().downstream_atom( n ) );
		AtomID const jump_atom_id2( named_atom_id_to_atom_id(NamedAtomID( jump_atom2, jump_pos2 ), pose ) );

		// No -- should look at actual atoms that are connected by jump, not residues.
		//if (allow_insert_->get( jump_pos1 ) || allow_insert_->get( jump_pos2 ) ) 	 mm.set_jump( n, true );

		if ( allow_insert_->get( jump_atom_id1 ) || allow_insert_->get( jump_atom_id2 ) ) 	 mm.set_jump( n, true );

	}

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


void
RNA_Minimizer::set_score_function( core::scoring::ScoreFunctionOP const & scorefxn ){
	scorefxn_ = scorefxn->clone();
}

void
RNA_Minimizer::set_allow_insert(toolbox::AllowInsertOP allow_insert ){
	allow_insert_ = allow_insert;
}

} // namespace rna
} // namespace protocols
