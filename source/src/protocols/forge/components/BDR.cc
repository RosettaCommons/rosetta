// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/components/BDR.cc
/// @brief
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/components/BDR.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/SegmentInsert.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/methods/chainbreak_eval.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/forge/methods/util.hh>

// project headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.hh>

// C++ headers
#include <utility>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace forge {
namespace components {


static basic::Tracer TR( "protocols.forge.components.BDR" );


/// @brief default constructor
BDR::BDR() :
	Super( "BDR" ),
	use_fullmer_( false ),
	use_sequence_bias_( false ),
	max_linear_chainbreak_( 0.07 ),
	centroid_loop_mover_str_( "RemodelLoopMover" ),
	redesign_loop_neighborhood_( true ),
	dr_cycles_( 3 ),
	centroid_sfx_( core::scoring::ScoreFunctionFactory::create_score_function( "remodel_cen" ) ),
	fullatom_sfx_( core::scoring::get_score_function() )
{}


/// @brief copy constructor
BDR::BDR( BDR const & rval ) :
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	manager_( rval.manager_ ),
	design_info_( rval.design_info_ ),
	use_fullmer_( rval.use_fullmer_ ),
	use_sequence_bias_( rval.use_sequence_bias_ ),
	max_linear_chainbreak_( rval.max_linear_chainbreak_ ),
	centroid_loop_mover_str_( rval.centroid_loop_mover_str_ ),
	redesign_loop_neighborhood_( rval.redesign_loop_neighborhood_ ),
	resfile_( rval.resfile_ ),
	dr_cycles_( rval.dr_cycles_ ),
	centroid_sfx_( rval.centroid_sfx_->clone() ),
	fullatom_sfx_( rval.fullatom_sfx_->clone() )
{
	if ( rval.vlb_.get() ) {
		vlb_ = VarLengthBuildOP( new VarLengthBuild( *rval.vlb_ ) );
	}
}


/// @brief default destructor
BDR::~BDR() = default;


/// @brief clone this object
BDR::MoverOP BDR::clone() const {
	return BDR::MoverOP( new BDR( *this ) );
}


/// @brief create this type of object
BDR::MoverOP BDR::fresh_instance() const {
	return BDR::MoverOP( new BDR() );
}


/// @brief the centroid level score function, default "remodel_cen"
BDR::ScoreFunction const & BDR::centroid_scorefunction() const {
	return *centroid_sfx_;
}


/// @brief the full-atom level score function
BDR::ScoreFunction const & BDR::fullatom_scorefunction() const {
	return *fullatom_sfx_;
}


/// @brief add instruction to the manager of this BDR (no copy)
/// @param[in] bi BuildInstruction
/// @param[in] aa_during_design_refine The allowed amino acid sequence
///  during design.  Only applicable to BuildInstructions like
///  SegmentRebuild and SegmentInsert.  Make sure the length of this
///  string matches up properly.  Default empty string.
void BDR::add_instruction(
	BuildInstructionOP bi,
	String const & aa_during_design_refine
)
{
	manager_.add( bi );
	if ( !aa_during_design_refine.empty() ) {
		design_info_.push_back( std::make_pair( bi->original_interval(), aa_during_design_refine ) );
	}

	// additional instruction means we'll need a new re-init the VLB, so
	// go ahead and drop the existing one
	vlb_.reset();
}


/// @brief create directed dependency between two instructions
void BDR::create_directed_dependency(
	BuildInstructionOP u,
	BuildInstructionOP v
)
{
	manager_.create_directed_dependency( u, v );
}


/// @brief set the centroid level score function
void BDR::centroid_scorefunction( ScoreFunction const & sfx ) {
	centroid_sfx_ = sfx.clone();
}


/// @brief set the centroid level score function
void BDR::centroid_scorefunction( ScoreFunctionOP sfx ) {
	centroid_sfx_ = sfx->clone();
}


/// @brief set the full-atom level score function
void BDR::fullatom_scorefunction( ScoreFunction const & sfx ) {
	fullatom_sfx_ = sfx.clone();
}


/// @brief set the full-atom level score function
void BDR::fullatom_scorefunction( ScoreFunctionOP sfx ) {
	fullatom_sfx_ = sfx->clone();
}


/// @brief apply defined moves to given Pose
void BDR::apply( Pose & pose ) {
	using core::pose::metrics::CalculatorFactory;
	using core::pose::metrics::PoseMetricCalculatorOP;
	using basic::MetricValue;
	using core::scoring::dssp::Dssp;
	using protocols::moves::MS_SUCCESS;
	using protocols::moves::FAIL_DO_NOT_RETRY;
	using protocols::moves::FAIL_RETRY;
	using protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator;
	using protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator;

	// assign secondary structure
	Dssp dssp( pose );
	dssp.insert_ss_into_pose( pose );

	// do centroid build
	if ( !centroid_build( pose ) ) { // build failed
		set_last_move_status( FAIL_RETRY );
		return;
	}

	// setup calculators
	CalculatorFactory::Instance().remove_calculator( neighborhood_calc_name() );
	CalculatorFactory::Instance().register_calculator(
		neighborhood_calc_name(),
		PoseMetricCalculatorOP( new NeighborhoodByDistanceCalculator( manager_.union_of_intervals_containing_undefined_positions() ) )
	);

	// do design-refine iteration
	if ( dr_cycles_ > 0 ) {

		if ( !design_refine( pose ) ) { // design-refine failed
			set_last_move_status( FAIL_RETRY );
			return;
		}

	}

	// if we've gotten to this point, then the structure has been
	// built properly
	set_last_move_status( MS_SUCCESS );

	// setup the PoseMetricCalculators and add them to the evaluators in the
	// JobOutputter
	CalculatorFactory::Instance().remove_calculator( loops_buns_polar_calc_name() );
	CalculatorFactory::Instance().remove_calculator( neighborhood_buns_polar_calc_name() );

	CalculatorFactory::Instance().register_calculator(
		loops_buns_polar_calc_name(),
		PoseMetricCalculatorOP( new BuriedUnsatisfiedPolarsCalculator(
		"default",
		"default",
		manager_.union_of_intervals_containing_undefined_positions()
		) )
	);

	MetricValue< std::set< Size > > loops_neighborhood;
	pose.metric( neighborhood_calc_name(), "neighbors", loops_neighborhood );
	CalculatorFactory::Instance().register_calculator(
		neighborhood_buns_polar_calc_name(),
		PoseMetricCalculatorOP( new BuriedUnsatisfiedPolarsCalculator(
		"default",
		"default",
		loops_neighborhood.value()
		) )
	);
}


std::string
BDR::get_name() const {
	return "BDR";
}

/// @brief run the centroid level build stage
/// @return true if loop closed, false otherwise
bool BDR::centroid_build(
	Pose & pose
) {
	using core::scoring::ScoreFunctionOP;
	using core::scoring::ScoreFunctionFactory;
	using protocols::moves::MS_SUCCESS;

	using core::util::switch_to_residue_type_set;
	using protocols::forge::methods::restore_residues;
	using protocols::toolbox::pose_manipulation::construct_poly_ala_pose;

	// safety, clear the energies object
	pose.energies().clear();

	// make backup Pose for transferring sidechains
	Pose archive_pose = pose;
	Pose modified_archive_pose = archive_pose;
	manager_.modify( modified_archive_pose );

	// ensure modified_archive_pose is completely full-atom, otherwise mismatch
	// will occur when restoring sidechains at the end of the procedure
	bool mod_ap_is_full_atom = true;
	for ( Size i = 1, ie = modified_archive_pose.size(); mod_ap_is_full_atom && i != ie; ++i ) {
		mod_ap_is_full_atom &= ( modified_archive_pose.residue( i ).type().mode() == core::chemical::FULL_ATOM_t );
	}

	if ( !mod_ap_is_full_atom ) {
		core::util::switch_to_residue_type_set( modified_archive_pose, core::chemical::FULL_ATOM_t );
	}

	// flip to poly-ala-gly-pro-disulf pose
	utility::vector1< Size > protein_residues;
	for ( Size i = 1, ie = pose.size(); i <= ie; ++i ) {
		if ( pose.residue( i ).is_protein() ) {
			protein_residues.push_back( i );
		}
	}

	construct_poly_ala_pose( pose, protein_residues, true, true, true );

	// run VLB to build the new section, if no segments have been added/deleted
	// we use the same VLB so that fragment caching works properly
	if ( !vlb_.get() ) {
		vlb_ = VarLengthBuildOP( new VarLengthBuild( manager_ ) );
	}

	vlb_->scorefunction( centroid_sfx_ );
	vlb_->vall_memory_usage( VLB_VallMemoryUsage::CLEAR_IF_CACHING_FRAGMENTS );
	vlb_->use_fullmer( use_fullmer_ );
	vlb_->max_linear_chainbreak( max_linear_chainbreak_ );
	vlb_->loop_mover_str( centroid_loop_mover_str_ );

	if ( use_sequence_bias_ ) {
		vlb_->original_sequence( archive_pose.sequence() );
	}

	vlb_->apply( pose );

	if ( vlb_->get_last_move_status() == MS_SUCCESS ) {

		// record the used manager w/ all mapping info
		manager_ = vlb_->manager();

		// safety, clear all the energies before restoring full-atom residues and
		// scoring
		pose.energies().clear();

		// Swap back original sidechains.  At the moment this is a two step process
		// in case any sidechains from SegmentInsert and the like that aren't in the
		// original archive pose need to be transferred.
		restore_residues( modified_archive_pose, pose );
		restore_residues( manager_.original2modified(), archive_pose, pose );

		// go ahead and score w/ full-atom here; we do this in case there are no
		// design-refine cycles -- it's useful to have e.g. rama in the output
		(*fullatom_sfx_)( pose );

		return true; // loop closed

	} else {

		pose = archive_pose;

	}

	return false; // false if loop not closed
}


/// @brief run the design-refine stage
/// @return currently always true
bool BDR::design_refine(
	Pose & pose
)
{
	using core::kinematics::FoldTree;
	using core::pack::task::operation::RestrictResidueToRepacking;
	using core::pack::task::operation::RestrictResidueToRepackingOP;
	using core::pack::task::operation::RestrictToRepacking;
	using core::pack::task::operation::TaskOperationCOP;
	using core::scoring::ScoreFunctionOP;
	using core::scoring::ScoreFunctionFactory;
	using protocols::forge::build::SegmentInsert;
	using protocols::loops::Loops;
	using protocols::loops::loop_mover::refine::LoopMover_Refine_CCD;
	using protocols::minimization_packing::PackRotamersMover;
	using protocols::toolbox::task_operations::RestrictToNeighborhoodOperation;

	using core::pose::annotated_to_oneletter_sequence;
	using protocols::forge::methods::intervals_to_loops;
	using protocols::forge::methods::linear_chainbreak;
	using protocols::loops::remove_cutpoint_variants;

	using Positions = protocols::forge::build::BuildManager::Positions;

	// collect new regions/positions
	std::set< Interval > loop_intervals = manager_.intervals_containing_undefined_positions();
	Original2Modified original2modified_interval_endpoints = manager_.original2modified_interval_endpoints();

	// collect loops
	loops::LoopsOP loops( new Loops( intervals_to_loops( loop_intervals.begin(), loop_intervals.end() ) ) );

	// refine Mover used doesn't setup a fold tree, so do it here
	FoldTree loop_ft = protocols::forge::methods::fold_tree_from_loops( pose, *loops );

	// save original fold tree
	FoldTree original_ft = pose.fold_tree();

	// define the score function
	ScoreFunctionOP sfx = fullatom_sfx_->clone();

	// setup the design TaskFactory
	TaskFactoryOP design_tf = generic_taskfactory();
	design_tf->push_back( TaskOperationCOP( new RestrictToNeighborhoodOperation( neighborhood_calc_name() ) ) );

	if ( !redesign_loop_neighborhood_ ) {
		// set repack only for non-loop positions
		RestrictResidueToRepackingOP repack_op( new RestrictResidueToRepacking() );

		Positions new_positions = manager_.new_positions();
		for ( Size i = 1, ie = pose.size(); i != ie; ++i ) {
			if ( new_positions.find( i ) == new_positions.end() ) {
				repack_op->include_residue( i );
			}
		}

		design_tf->push_back( repack_op );
	}

	// add explicit residue types to design task factory if requested
	for ( DesignInfo::const_iterator di = design_info_.begin(), die = design_info_.end(); di != die; ++di ) {
		if ( !di->second.empty() ) {
			String const aa = core::pose::annotated_to_oneletter_sequence( di->second );

			debug_assert( original2modified_interval_endpoints.find( di->first.left ) != original2modified_interval_endpoints.end() );

			if ( aa.find( SegmentInsert::insertion_char() ) != String::npos ) { // SegmentInsert style
				process_insert_design_string( di->first, aa, original2modified_interval_endpoints, design_tf );
			} else { // SegmentRebuild style
				process_continuous_design_string( di->first, aa, original2modified_interval_endpoints, design_tf );
			}
		}
	}

	// setup the refine TaskFactory
	TaskFactoryOP refine_tf = generic_taskfactory();
	refine_tf->push_back( TaskOperationCOP( new RestrictToNeighborhoodOperation( neighborhood_calc_name() ) ) );
	refine_tf->push_back( TaskOperationCOP( new RestrictToRepacking() ) );

	// safety, clear the energies object
	pose.energies().clear();

	// run design-refine cycle
	for ( Size i = 0; i < dr_cycles_; ++i ) {

		// design the new section
		PackRotamersMover design( sfx );
		design.task_factory( design_tf );
		design.apply( pose );

		// set loop topology
		pose.fold_tree( loop_ft );

		// refine the new section
		LoopMover_Refine_CCD refine( loops, sfx );
		refine.false_movemap( manager_.movemap_as_OP() );
		refine.set_task_factory( refine_tf );
		refine.apply( pose );

		// remove cutpoint variants -- shouldn't this happen at the end
		// of the refine Mover?
		remove_cutpoint_variants( pose );

		// set original topology
		pose.fold_tree( original_ft );
	}

	// must score one last time since we've removed variants and set
	// new topology, otherwise component energies not correct for
	// e.g. structure output
	(*sfx)( pose );

	// evaluate all chainbreaks using linear chainbreak
	bool cbreaks_pass = true;
	for ( auto l = loops->begin(), le = loops->end(); l != le && cbreaks_pass; ++l ) {
		if ( l->cut() > 0 ) {
			Real const c = linear_chainbreak( pose, l->cut() );
			TR << "design_refine: final chainbreak = " << c << std::endl;
			cbreaks_pass = c <= max_linear_chainbreak_;
		}
	}

	return cbreaks_pass;
}


/// @brief return a TaskFactory useable as a starting point for either
///  design or refinement
BDR::TaskFactoryOP BDR::generic_taskfactory() {
	using core::pack::task::operation::IncludeCurrent;
	using core::pack::task::operation::InitializeFromCommandline;
	using core::pack::task::operation::ReadResfile;
	using core::pack::task::operation::ReadResfileOP;
	using core::pack::task::operation::TaskOperationCOP;
	using core::pack::task::TaskFactory;
	using core::pack::task::operation::NoRepackDisulfides;

	TaskFactoryOP tf( new TaskFactory() );

	tf->push_back( TaskOperationCOP( new InitializeFromCommandline() ) ); // also inits -ex options
	tf->push_back( TaskOperationCOP( new IncludeCurrent() ) ); // enforce keeping of input sidechains
	tf->push_back( TaskOperationCOP( new NoRepackDisulfides() ) );

	// load resfile op only if requested
	if ( !resfile_.empty() ) {
		ReadResfileOP rrf( new ReadResfile() );
		rrf->filename( resfile_ );
		tf->push_back( rrf );
	}

	return tf;
}


/// @brief process a continuous design string, adding appropriate operations
///  to the TaskFactory
void BDR::process_continuous_design_string(
	Interval const & original_interval,
	String const & design_str,
	Original2Modified const & original2modified_interval_endpoints,
	TaskFactoryOP design_tf
)
{
	using core::pack::task::operation::RestrictAbsentCanonicalAAS;
	using core::pack::task::operation::TaskOperationCOP;
	using core::chemical::aa_from_oneletter_code;

	Size const offset = original2modified_interval_endpoints.find( original_interval.left )->second;
	for ( Size i = 0, ie = design_str.length(); i < ie; ++i ) {
		utility::vector1< bool > allowed_aa_types( 20, false );

		switch ( design_str.at( i ) ) {
		case 's' : // surface case, no CFWY
			allowed_aa_types = allowed_surface_aa();
			break;
		case '.' : // protocol default design
			continue;
		default : // regular case, single aa type
			allowed_aa_types[ aa_from_oneletter_code( design_str.at( i ) ) ] = true;
			break;
		}

		design_tf->push_back( TaskOperationCOP( new RestrictAbsentCanonicalAAS( i + offset, allowed_aa_types ) ) );
	}
}


/// @brief process a design string containing an insert, adding appropriate
///  operations to the TaskFactory
void BDR::process_insert_design_string(
	Interval const & original_interval,
	String const & design_str,
	Original2Modified const & original2modified_interval_endpoints,
	TaskFactoryOP design_tf
)
{
	using core::pack::task::operation::RestrictAbsentCanonicalAAS;
	using core::pack::task::operation::RestrictResidueToRepacking;
	using core::pack::task::operation::RestrictResidueToRepackingOP;
	using core::pack::task::operation::TaskOperationCOP;
	using protocols::forge::build::Interval;
	using protocols::forge::build::SegmentInsert;

	using core::chemical::aa_from_oneletter_code;

	char const insert_char = SegmentInsert::insertion_char();

	// Figure out the number of residues in each section.
	Interval const interval(
		original2modified_interval_endpoints.find( original_interval.left )->second,
		original2modified_interval_endpoints.find( original_interval.right )->second
	);

	Size const insert_char_idx = design_str.find( insert_char );
	Size const left_nres = insert_char_idx;
	Size const right_nres = design_str.size() - left_nres - 1;
	Size const insert_nres = interval.length() - left_nres - right_nres;

	// Make setup easy by building a new design string to expand the
	// insertion character into a series of the insertion character
	// the size of the insert.
	String aa = design_str;
	aa.replace( insert_char_idx, 1, insert_nres, insert_char );

	// setup TaskOperations
	RestrictResidueToRepackingOP repack_op( new RestrictResidueToRepacking() );

	Size const left_offset = interval.left;
	for ( Size i = 0, ie = aa.size(); i < ie; ++i ) {
		utility::vector1< bool > allowed_aa_types( 20, false );

		if ( aa.at( i ) == insert_char ) { // repack only
			repack_op->include_residue( i + left_offset );
			continue;
		} else if ( aa.at( i ) == 's' ) { // surface case, no CFWY
			allowed_aa_types = allowed_surface_aa();
		} else if ( aa.at( i ) == '.' ) { // protocol default design
			continue;
		} else { // regular case, single aa type
			allowed_aa_types[ aa_from_oneletter_code( aa.at( i ) ) ] = true;
		}

		design_tf->push_back( TaskOperationCOP( new RestrictAbsentCanonicalAAS( i + left_offset, allowed_aa_types ) ) );
	}

	design_tf->push_back( repack_op );
}


/// @brief return a boolean vector specifying allowed a.a. when designing
///  on the surface
utility::vector1< bool > const & BDR::allowed_surface_aa() {
	using core::chemical::aa_from_oneletter_code;

	static String surface_aa = "ADEGHIKLMNPQRSTV";
	static utility::vector1< bool > v( 20, false );

	for ( char i : surface_aa ) {
		v[ aa_from_oneletter_code( i ) ] = true;
	}

	return v;
}


} // namespace components
} // namespace forge
} // namespace protocols
