// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseProteinMinimizer
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @details
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/legacy/modeler/protein/StepWiseProteinMinimizer.hh>
#include <protocols/stepwise/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/output_util.hh>
#include <protocols/stepwise/modeler/movemap/util.hh>
#include <protocols/stepwise/modeler/util.hh>

//////////////////////////////////
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

#include <string>

//Auto Headers
#include <utility/vector1.hh>
using namespace core;
using core::Real;
using ObjexxFCL::format::F;

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.legacy.modeler.protein.StepWiseProteinMinimizer" );

using namespace protocols::stepwise::modeler;
using namespace protocols::stepwise::modeler::protein;
using namespace protocols::stepwise::modeler::movemap;

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace protein {


//////////////////////////////////////////////////////////////////////////
StepWiseProteinMinimizer::StepWiseProteinMinimizer( utility::vector1< pose::PoseOP > const & pose_list,
	utility::vector1< Size > const & moving_residues ):
	Mover(),
	moving_residues_( moving_residues ),
	pose_list_( pose_list )
{
	initialize_parameters();
}

//////////////////////////////////////////////////////////////////////////
//destructor
StepWiseProteinMinimizer::~StepWiseProteinMinimizer()
{}

/////////////////////
std::string
StepWiseProteinMinimizer::get_name() const {
	return "StepWiseProteinMinimizer";
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinMinimizer::initialize_parameters(){
	Mover::type( "StepWiseProteinMinimizer" );
	move_takeoff_torsions_ = true;
	rescore_only_ = false;
	move_jumps_between_chains_ = false;
	//fa_scorefxn_ = core::scoring::get_score_function();
	min_type_ = "lbfgs_armijo_nonmonotone"; // used to be dfpmin
	cartesian_ = true;
	min_tolerance_ = 0.000025 ; // used to be 0.00000025
	use_coordinate_constraints_ = true;
	num_pose_minimize_ = 0; // signal to minimize all
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinMinimizer::apply( core::pose::Pose & pose )
{
	using namespace core::optimization;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::pose;

	clock_t const time_start( clock() );

	ConstraintSetOP cst_set;
	if ( use_coordinate_constraints_ ) cst_set = pose.constraint_set()->clone();

	utility::vector1< std::pair<core::Size,core::Size> > disulfides;
	core::conformation::disulfide_bonds(pose.conformation(), disulfides);
	runtime_assert( fa_scorefxn_ != 0 );

	CartesianMinimizer cart_minimizer;
	AtomTreeMinimizer minimizer;
	bool const use_nblist( true );
	MinimizerOptions options( min_type_, min_tolerance_, use_nblist, false, false );
	options.nblist_auto_update( true );

	kinematics::MoveMap mm_start, mm;
	//  figure_out_moving_residues( mm_start, pose, fixed_res_, move_takeoff_torsions_, move_jumps_between_chains_ );
	//  output_movemap( mm_start, pose, TR );

	// testing unification...
	utility::vector1< Size > working_minimize_res;
	for ( Size n = 1; n <= pose.total_residue(); n++ ) { if ( !fixed_res_.has_value( n ) ) working_minimize_res.push_back( n );}
	figure_out_stepwise_movemap( mm_start, pose, working_minimize_res, move_takeoff_torsions_ );
	//  output_movemap( mm_start, pose, TR );

	mm = mm_start;
	Real const original_coordinate_cst_weight = fa_scorefxn_->get_weight( coordinate_constraint );

	utility::vector1< PoseOP > output_pose_list;
	for ( Size n = 1; n <= pose_list_.size(); n++ ) {

		if ( num_pose_minimize_ > 0 &&  n > num_pose_minimize_ ) break;

		pose = *pose_list_[ n ];

		// Following are necessary because poses from clustering went thorugh silent struct and lost their constraints & disulfide information.
		if ( cst_set ) pose.constraint_set( cst_set );
		pose.conformation().fix_disulfides( disulfides );

		Real const score_original = (*fa_scorefxn_)( pose );

		// The movemap has all dofs for "non-fixed residues" free to move.
		// We can also let sidechains minimize in fixed-residues -- for
		// speed only look at neighbors of moving residues.
		mm = mm_start;
		let_neighboring_chis_minimize( mm, pose );

		if ( !rescore_only_ ) {

			if ( use_coordinate_constraints_ ) {
				// One minimize with loose coordinate tethers to make sure the pose doesn't blow up.
				core::scoring::constraints::add_coordinate_constraints( pose );
				if ( fa_scorefxn_->has_zero_weight( coordinate_constraint) ) fa_scorefxn_->set_weight( coordinate_constraint, 1.0 );
				minimizer.run( pose, mm, *fa_scorefxn_, options );
				// Now a regular minimize.
				pose.constraint_set( cst_set ); // return original constraints (no added coordinate constraints)
				fa_scorefxn_->set_weight( coordinate_constraint, original_coordinate_cst_weight );
			}

			// for poses with chainbreaks, do an initial minimization with a weak linear_chainbreak term. (anneal it in.)
			if ( pose_has_chainbreak( pose ) ) {

				Real const linear_chainbreak_weight_original = fa_scorefxn_->get_weight( linear_chainbreak );
				if ( linear_chainbreak_weight_original < 20.0 ) std::cout << "WARNING!! Your linear_chainbreak weight is " << F(8,3,linear_chainbreak_weight_original ) << ", which is less than recommended (20.0) " << std::endl;

				fa_scorefxn_->set_weight( linear_chainbreak, linear_chainbreak_weight_original * 0.25 );
				if ( !cartesian_ ) {
					minimizer.run( pose, mm, *fa_scorefxn_, options );
				}
				fa_scorefxn_->set_weight( linear_chainbreak, linear_chainbreak_weight_original );
			}

			if ( cartesian_ ) {
				cart_minimizer.run( pose, mm, *fa_scorefxn_, options );
			} else {
				minimizer.run( pose, mm, *fa_scorefxn_, options );
			}

		}

		output_pose_list.push_back( pose.clone() );
		TR.Debug << "Score minimized from " << F(8,3, score_original) << " to " << F(8,3,(*fa_scorefxn_)( pose )) << std::endl;
	}

	pose_list_ = output_pose_list;

	TR.Debug << "Total time in StepWiseProteinMinimizer: " <<
		static_cast<Real>(clock() - time_start) / CLOCKS_PER_SEC << std::endl;

}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinMinimizer::let_neighboring_chis_minimize(
	core::kinematics::MoveMap & mm,
	core::pose::Pose & pose ){

	using namespace core::scoring;

	(*fa_scorefxn_)( pose );
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	for ( Size n = 1; n <= moving_residues_.size(); n++ ) {
		Size const i = moving_residues_[ n ];
		if ( pose.residue(i).is_protein() ) { // these should be activated, but make sure . VIRTUAL_SIDE_CHAIN issue!
			mm.set_chi( i, true );
		}
	}

	for ( Size n = 1; n <= moving_residues_.size(); n++ ) {

		Size const i = moving_residues_[ n ];

		for ( graph::Graph::EdgeListConstIter
				iter = energy_graph.get_node( i )->const_edge_list_begin();
				iter != energy_graph.get_node( i )->const_edge_list_end();
				++iter ) {

			Size j( (*iter)->get_other_ind( i ) );
			if ( pose.residue(j).has_variant_type( core::chemical::VIRTUAL_RESIDUE_VARIANT ) ) continue;

			if ( pose.residue(j).is_protein() ) {
				mm.set_chi( j, true );
			} else if ( pose.residue(j).is_RNA() ) {
				mm.set( id::TorsionID( j, id::CHI, 4), true ); // 2'-OH.
			}

		}
	}

}

//////////////////////////////////////////////////////////////////////////
bool
StepWiseProteinMinimizer::pose_has_chainbreak( pose::Pose const & pose ){
	// this is pretty conservative -- actually  there might be
	// cases where the pose has a chainbreak but the minimized dofs would
	// not affect the relative positions of the chainbreak residues.
	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( pose.residue_type(i).has_variant_type( "CUTPOINT_UPPER" ) ) return true;
		if ( pose.residue_type(i).has_variant_type( "CUTPOINT_LOWER" ) ) return true;
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinMinimizer::set_min_type( std::string const & min_type ){
	min_type_ = min_type;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinMinimizer::set_min_tolerance( Real const & min_tolerance ){
	min_tolerance_ = min_tolerance;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinMinimizer::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
	fa_scorefxn_ = scorefxn;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinMinimizer::set_fixed_res( utility::vector1< core::Size > const & fixed_res ){
	fixed_res_ = fixed_res;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinMinimizer::set_calc_rms_res( utility::vector1< core::Size > const & calc_rms_res ){
	calc_rms_res_ = calc_rms_res;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinMinimizer::set_cartesian( bool const setting ){
	cartesian_ = setting;
	if ( cartesian_ ) min_type_ = "lbfgs_armijo_nonmonotone";
}


} //protein
} //modeler
} //legacy
} //stepwise
} //protocols
