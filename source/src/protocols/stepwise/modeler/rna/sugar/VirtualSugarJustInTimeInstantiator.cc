// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/sugar/VirtualSugarJustInTimeInstantiator.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/rna/sugar/VirtualSugarJustInTimeInstantiator.hh>
#include <protocols/stepwise/modeler/rna/sugar/VirtualSugarSampler.hh>
#include <protocols/stepwise/modeler/rna/sugar/util.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <core/chemical/rna/util.hh> // for SYN
#include <ObjexxFCL/string.functions.hh>
#include <utility/tools/make_vector1.hh>
#include <basic/Tracer.hh>

///////////////////////////////////////////////////////////////////////////////////
//////////////////////Build previously virtualized sugar/////////////////////
// A virtualized sugar occurs when a previous move was a 'floating base'
// step, which only samples euler angles of a base, but virtualized the
// attached sugar and any residues connecting that nucleotide to the
//'instantiated' body of the RNA.
//
// There are potentially four different virtualized sugar positions,
//  and here is a diagram of the most extreme case, merging of
//  two helices whose edge base pairs were all somehow created through
//  bulge skipping moves, leaving 4 virtual riboses (marked *):
//
//             A1 -- U16                 |
//             C2 -- G15                 | POSE 1
// bulge -> (U3)       (U14) <- bulge    |
//            *G4 -- C13*                |
//    moving-> |      x <- chain break to close
//            *G5 -- C12*                |
// bulge -> (U6)       (U11) <- bulge    | POSE2
//             C7 -- G10                 |
//             C8 -- G9                  |
//
// (1) Anchor_sugar (most common case) :
//      This is the sugar of the nucleotide that was built immediately
//      before the current/moving nucleotide along the chain. Corresponds
//      to sugar of residue (moving_res - 1) if current step is built in
//      the forward (3') direction. Likewise, corresponds to sugar of
//      residue (moving_res + 1) if current step built in the backward
//      (5') direction.
//
// (2) Current_sugar (rare) :
//      The sugar of the current/moving nucleotide. This sugar can
//      be virtual in the situation where the current step is a step to
//      combine two moving_elements that were previously built with SWA. It is
//      possible that the current nucleotide (in the moving moving_element) was
//      previously built with a 'floating base' step and therefore will
//      contain a virtual sugar.
//
// (3) Five_prime_chain_break_sugar :
//      Occurs only in the context where the current step is a
//      chain-closure (closing loop) step. This is the sugar of the
//      nucleotide immediately 5' of the chain-closure phosphate group.
//      Five_prime_chain_break_sugar can be identical to current_sugar
//      in certain situations. The code below contains check statements to
//      prevent double-counting in these situations.
//
// (4) Three_prime_chain_break_sugar :
//      Occurs only in the context where the current step is a
//      chain-closure (closing loop) step. This is the sugar of the
//      nucleotide immediately 3' of the chain-closure phosphate group.
//      Three_prime_chain_break_sugar can be identical to current_sugar
//      in certain situations. The code below contains check statements to
//      prevent double-counting in these situations.

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.rna.sugar.VirtualSugarJustInTimeInstantiator" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace sugar {

//Constructor
VirtualSugarJustInTimeInstantiator::VirtualSugarJustInTimeInstantiator( working_parameters::StepWiseWorkingParametersCOP & working_parameters ):
	working_parameters_( working_parameters ),
	moving_res_(  working_parameters_->working_moving_res() ),
	rebuild_bulge_mode_( working_parameters_->rebuild_bulge_mode() ),
	keep_base_fixed_( false ),
	moving_res_legacy_( false )
{}

//Destructor
VirtualSugarJustInTimeInstantiator::~VirtualSugarJustInTimeInstantiator()
{}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
VirtualSugarJustInTimeInstantiator::apply( core::pose::Pose & pose ){
	success_ = do_the_modeler( pose );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
VirtualSugarJustInTimeInstantiator::do_the_modeler( core::pose::Pose & pose ){

	output_title_text( "Build previously virtualize sugar", TR.Debug );

	// special: moving_res and its anchor res are treated specially. They might be involved in a closed cutpoint,
	// but even if not, we find possible sugar conformations which will be screened for stereochemistry in ConnectionSampler.
	utility::vector1< Size > check_res = utility::tools::make_vector1( moving_res_ );
	if ( pose.residue_type( moving_res_ ).is_RNA() && working_parameters_->working_reference_res() > 0 ) {
		check_res.push_back( working_parameters_->working_reference_res() );
	}

	// definitely need to instantiate & sample sugars if they are about to get involved in chain closure...
	utility::vector1< Size > cutpoints_closed = figure_out_moving_cutpoints_closed( pose, working_parameters_->working_moving_partition_res() );
	for ( Size n = 1; n <= cutpoints_closed.size(); n++ ) {
		check_res.push_back( cutpoints_closed[ n ]     );
		check_res.push_back( cutpoints_closed[ n ] + 1 );
	}

	if ( options_->allow_rebuild_bulge_mode() && rebuild_bulge_mode_ ) {
		// this was used in old SWA RNA workflow -- probably can be deprecated, along with VirtualResidueRNA variant types.
		TR << TR.Red << "In: REBUILD_BULGE_MODE " << " for residue " << moving_res_ << TR.Reset << std::endl;
		core::pose::rna::apply_virtual_rna_residue_variant_type( pose,  moving_res_ );
	}

	bool success( true );
	for ( Size n = 1; n <= check_res.size(); n++ ) {
		Size const i = check_res[ n ];
		success = get_sugar_modeling_set( pose, i );
		if ( !success ) break;
	}

	return success;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
VirtualSugarJustInTimeInstantiator::sampled_sugar_index( Size const i ) {
	for ( Size n = 1; n <= sugar_modeling_sets_.size(); n++ ) {
		if ( sugar_modeling_sets_[ n ]->moving_res == i ) return n;
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
VirtualSugarJustInTimeInstantiator::get_sugar_modeling_set( pose::Pose & viewer_pose, Size const i ){
	if ( sampled_sugar_index( i ) > 0 ) return true;

	SugarModelingOP sugar_modeling( new SugarModeling() );
	if ( moving_res_legacy_ && i == moving_res_ && working_parameters_->floating_base() ) return true;

	bool did_setup = setup_sugar_modeling( viewer_pose, i, *sugar_modeling);
	if ( !did_setup ) return true;

	PoseOP pose_save = viewer_pose.clone();
	do_sugar_sampling( viewer_pose, *sugar_modeling, ObjexxFCL::string_of( i ) );

	if ( sugar_modeling->pose_list.size() == 0 ) return false;
	std::sort( sugar_modeling->pose_list.begin(), sugar_modeling->pose_list.end(), sort_pose_by_score );
	sugar_modeling_sets_.push_back( sugar_modeling );
	residue_alternative_sets_.push_back( convert_sugar_modeling_to_residue_alternative_set( *sugar_modeling ) );
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
VirtualSugarJustInTimeInstantiator::do_sugar_sampling( pose::Pose & viewer_pose, SugarModeling & sugar_modeling, std::string const name ){

	VirtualSugarSampler virtual_sugar_sampler( working_parameters_, sugar_modeling );
	virtual_sugar_sampler.set_tag( name );
	virtual_sugar_sampler.set_scorefxn( scorefxn_ );
	virtual_sugar_sampler.set_integration_test_mode( options_->integration_test_mode() );
	virtual_sugar_sampler.set_do_minimize( options_->virtual_sugar_do_minimize() );
	if ( options_->integration_test_mode() )  virtual_sugar_sampler.set_do_minimize( false ); // override.
	virtual_sugar_sampler.set_use_phenix_geo( options_->use_phenix_geo() );
	virtual_sugar_sampler.set_legacy_mode( options_->virtual_sugar_legacy_mode() );
	virtual_sugar_sampler.set_keep_base_fixed( options_->virtual_sugar_keep_base_fixed() );
	virtual_sugar_sampler.set_do_screens( options_->virtual_sugar_do_screens() );
	bool const pick_one_virtual_sugar_ = options_->choose_random() && options_->virtual_sugar_do_minimize() /* takes a long time */;
	virtual_sugar_sampler.set_choose_random( pick_one_virtual_sugar_ );

	virtual_sugar_sampler.apply( viewer_pose );

	if ( sugar_modeling.pose_list.size() == 0 ) {
		TR.Debug << name + " sugar modeling unsuccessful!" << std::endl;
		return false;
	}
	return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
is_sugar_virtual( pose::Pose const & pose, Size const & n ){
	return ( pose.residue( n ).has_variant_type( core::chemical::VIRTUAL_RIBOSE ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
VirtualSugarJustInTimeInstantiator::setup_sugar_modeling( pose::Pose const & pose, Size const moving_res, SugarModeling & sugar_modeling ){
	using namespace core::chemical::rna;

	if ( !is_sugar_virtual( pose, moving_res ) ) return false;

	Size reference_res( 0 );
	if ( options_->virtual_sugar_do_screens() ) reference_res = get_reference_res_for_virtual_sugar_based_on_fold_tree( pose, moving_res );

	sugar_modeling = SugarModeling( moving_res, reference_res );
	sugar_modeling.set_base_and_pucker_state( pose, working_parameters_ );

	// model bulge?
	// this is assumed to be the residue immediately adjacent to the moving_residue,
	// in the same direction as moving_res. This would be the case in Parin's
	// usual dinucleotide move.  But we want to be able to totally leave out that
	// filler base, in which case there's no bulge to rebuild. -- rhiju
	Size const bulge_res_ = sugar_modeling.bulge_res;
	// also no point in doing chain closure if split across partitions...
	utility::vector1< Size > const & moving_partition_res = working_parameters_->working_moving_partition_res();
	if ( bulge_res_ == sugar_modeling.reference_res ||
			bulge_res_ == 0 ||
			!pose.residue_type( bulge_res_ ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ||
			moving_partition_res.has_value( sugar_modeling.moving_res ) != moving_partition_res.has_value( sugar_modeling.reference_res ) ) {
		TR.Debug <<  "CHECKING BULGE AT " << bulge_res_ << " is it bulge? " <<  pose.residue_type( bulge_res_ ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) << "  moving_res: " << sugar_modeling.moving_res << "  reference_res " << sugar_modeling.reference_res << std::endl;
		sugar_modeling.bulge_res   = 0;
		sugar_modeling.bulge_suite = 0;
	}

	return true;
}

////////////////////////////////////////////////////////////////////////////////////
void
VirtualSugarJustInTimeInstantiator::instantiate_sugars_recursively(  pose::Pose const & pose,
	utility::vector1< pose::PoseOP > & pose_data_list,
	utility::vector1< SugarModelingOP > const & sugar_modeling_sets,
	utility::vector1< Size > const & sugar_modeling_set_indices ) const {

	Size const which_set = sugar_modeling_set_indices.size() + 1;
	Size const num_sets = sugar_modeling_sets.size();
	if ( which_set <= num_sets ) {
		for ( Size n = 1; n <= sugar_modeling_sets[ which_set ]->pose_list.size(); n++ ) {
			utility::vector1< Size > indices_new = sugar_modeling_set_indices;
			indices_new.push_back( n );
			instantiate_sugars_recursively( pose, pose_data_list, sugar_modeling_sets, indices_new );
		}
	} else {
		runtime_assert( sugar_modeling_sets.size() == sugar_modeling_set_indices.size() );
		PoseOP start_pose = pose.clone();
		tag_into_pose( *start_pose, "" );
		for ( Size q = 1; q <= num_sets; q++ ) instantiate_sugar( *start_pose, *sugar_modeling_sets[ q ], sugar_modeling_set_indices[ q ] );
		pose_data_list.push_back( start_pose );
	}
}

////////////////////////////////////////////////////////////////////////////////////
void
VirtualSugarJustInTimeInstantiator::prepare_from_prior_sampled_sugar_jobs( pose::Pose const & pose,
	utility::vector1< PoseOP > & pose_data_list,
	bool const pose_explosion_legacy /* = false */) const {
	if ( pose_explosion_legacy ) {
		runtime_assert( modeler_sugar() );
		utility::vector1< Size > sugar_modeling_set_indices;
		instantiate_sugars_recursively( pose, pose_data_list, sugar_modeling_sets_, sugar_modeling_set_indices );
		if ( options_->virtual_sugar_legacy_mode() ) minimize_sugar_sets_legacy( pose, pose_data_list  );
	} else {
		PoseOP start_pose = pose.clone();
		instantiate_sugars_at_cutpoint_closed( *start_pose );
		pose_data_list.push_back( start_pose );
	}
}

//////////////////////////////////////////////////
void
VirtualSugarJustInTimeInstantiator::instantiate_sugars_at_cutpoint_closed( pose::Pose & pose ) const {
	for ( Size n = 1; n <= sugar_modeling_sets_.size(); n++ ) {
		Size const & sugar_res = sugar_modeling_sets_[n]->moving_res;
		if ( pose.residue( sugar_res ).has_variant_type( chemical::CUTPOINT_UPPER ) ||
				pose.residue( sugar_res ).has_variant_type( chemical::CUTPOINT_LOWER ) ||
				( sugar_res < pose.total_residue() && !pose.fold_tree().is_cutpoint( sugar_res ) ) ||
				( sugar_res > 1 && !pose.fold_tree().is_cutpoint( sugar_res - 1 ) ) ) {
			instantiate_sugar( pose, *sugar_modeling_sets_[ n ], 1 );
		}
	}
}


///////Ok, finally have to remove clashes that may arise due to the fact that the floating base sugar modeler and minimization were done individually of each other///
// NO! Floating bases are bulged. No reason for a virtual sugar sampler to move them. That can be another Sampler's job. --rd2013.
void
VirtualSugarJustInTimeInstantiator::minimize_sugar_sets_legacy( pose::Pose const & pose,
	utility::vector1< pose::PoseOP > & pose_data_list ) const {
	Pose pose_copy = pose;
	utility::vector1< SugarModeling > sugar_modeling_sets;
	for ( Size n = 1; n <= sugar_modeling_sets_.size(); n++ ) sugar_modeling_sets.push_back( *sugar_modeling_sets_[n] );
	minimize_all_sampled_floating_bases( pose_copy, sugar_modeling_sets, pose_data_list,
		scorefxn_, working_parameters_, true /*virtual_sugar_is_from_prior_step*/ );
}

////////////////////////////////////////////////////////////////////////////////////
// for backwards compatibility
void
VirtualSugarJustInTimeInstantiator::prepare_from_prior_sampled_sugar_jobs_for_chain_break( pose::Pose const & pose,
	utility::vector1< pose::PoseOP > & pose_data_list ) const {
	runtime_assert( modeler_sugar_at_chain_break() );
	utility::vector1< SugarModelingOP > sugar_modeling_sets_for_chainbreak = get_sugar_modeling_sets_for_chainbreak();
	utility::vector1< Size > sugar_modeling_set_indices;
	pose.dump_pdb( "BEFORE_INST.pdb" );
	instantiate_sugars_recursively( pose, pose_data_list, sugar_modeling_sets_for_chainbreak, sugar_modeling_set_indices );
	pose_data_list[1]->dump_pdb( "INSTANTIATED.pdb" );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< SugarModelingOP >
VirtualSugarJustInTimeInstantiator::get_sugar_modeling_sets_for_chainbreak() const {
	utility::vector1< SugarModelingOP > sets;
	for ( Size n = 1; n <= sugar_modeling_sets_.size(); n++ ) {
		Size const & sugar_res = sugar_modeling_sets_[n]->moving_res;
		if ( sugar_res == working_parameters_->working_moving_res()    ) continue;
		if ( sugar_res == working_parameters_->working_reference_res() ) continue;
		sets.push_back( sugar_modeling_sets_[ n ] );
	}
	return sets;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
VirtualSugarJustInTimeInstantiator::instantiate_sugar( pose::Pose & pose,
	SugarModeling const & sugar_modeling,
	Size const sugar_ID ) const {
	if ( sugar_modeling.pose_list.size() > 0 ) {
		tag_into_pose( pose,   tag_from_pose(pose) + tag_from_pose( *sugar_modeling.pose_list[ sugar_ID ] ) );
		copy_bulge_res_and_sugar_torsion( sugar_modeling, pose, ( *sugar_modeling.pose_list[ sugar_ID ] ), true /*instantiate_sugar*/ );
	} else {
		tag_into_pose( pose, tag_from_pose( pose ) + "_null" );
	}
}

/////////////////////
SugarModeling const &
VirtualSugarJustInTimeInstantiator::anchor_sugar_modeling(){
	Size const anchor_set_idx = sampled_sugar_index( working_parameters_->working_reference_res() );
	// kind of a hack, for backwards compatibility.
	if ( anchor_set_idx == 0 ) {
		SugarModelingOP blank_sugar_modeling( new SugarModeling() );
		return *blank_sugar_modeling;
	}
	return *sugar_modeling_sets_[ sampled_sugar_index( working_parameters_->working_reference_res() ) ];
}

/////////////////////
bool
VirtualSugarJustInTimeInstantiator::modeler_sugar() const{
	return ( num_sets() > 0 );
}

////////////////////////////////////////////////////////////////////////////////////
sampler::copy_dofs::ResidueAlternativeSet const &
VirtualSugarJustInTimeInstantiator::residue_alternative_set( Size const n ){
	runtime_assert( n <= residue_alternative_sets_.size() );
	return *residue_alternative_sets_[ n ];
}

////////////////////////////////////////////////////////////////////////////////////
bool
VirtualSugarJustInTimeInstantiator::modeler_sugar_at_chain_break() const{

	return ( get_sugar_modeling_sets_for_chainbreak().size() > 0 );
}

/////////////////////
void
VirtualSugarJustInTimeInstantiator::set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn ){ scorefxn_ = scorefxn; }

/////////////////////
std::string
VirtualSugarJustInTimeInstantiator::get_name() const {
	return "VirtualSugarJustInTimeInstantiator";
}

//////////////////////////////////////////////////////////////////
void
VirtualSugarJustInTimeInstantiator::set_options( options::StepWiseModelerOptionsCOP options ){
	options_ = options;
}

} //sugar
} //rna
} //modeler
} //stepwise
} //protocols
