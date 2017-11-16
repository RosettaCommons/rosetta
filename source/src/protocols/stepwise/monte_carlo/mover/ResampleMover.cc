// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/mover/ResampleMover.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/monte_carlo/mover/ResampleMover.hh>
#include <protocols/stepwise/monte_carlo/mover/ResampleMoverCreator.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMasterMover.hh> // just for setup_unified_stepwise_modeler...
#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMoveSelector.hh>
#include <protocols/stepwise/monte_carlo/mover/TransientCutpointHandler.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/modeler/StepWiseModeler.hh>
#include <protocols/stepwise/modeler/packer/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/chemical/rna/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>
#include <numeric/random/random.hh>

//Req'd on WIN32
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>

using namespace core;
using namespace protocols::stepwise::modeler;

static basic::Tracer TR( "protocols.stepwise.monte_carlo.mover.ResampleMover" );


//////////////////////////////////////////////////////////////////////////
// Makes a choice, based on current pose, and information in full_model_info
//  as to where to resample nucleotides or chunks.
//
// Does not add or delete  [see AddOrDeleteMover for those functionalities].
//
//////////////////////////////////////////////////////////////////////////

using namespace core::pose::full_model_info;


namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

//Constructor
ResampleMover::ResampleMover( protocols::stepwise::modeler::StepWiseModelerOP stepwise_modeler ):
	stepwise_modeler_( stepwise_modeler ),
	swa_move_selector_( StepWiseMoveSelectorOP( new StepWiseMoveSelector ) ),
	options_( protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptionsCOP( new protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptions ) ),
	minimize_single_res_( false ),
	slide_docking_jumps_( false )
{}

ResampleMover::ResampleMover() :
	stepwise_modeler_( nullptr ),
	swa_move_selector_( StepWiseMoveSelectorOP( new StepWiseMoveSelector ) ),
	options_( protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptionsCOP( new protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptions ) ),
	minimize_single_res_( false ),
	slide_docking_jumps_( false )
{
	stepwise_modeler_->set_native_pose( get_native_pose() );
}

//Destructor
ResampleMover::~ResampleMover()
{}

protocols::moves::MoverOP
ResampleMover::clone() const {
	return ResampleMoverOP( new ResampleMover( *this ) );
}

void
ResampleMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & /*pose*/ ) {

	// WARNING: don't use this with RNA scorefunctions; something goes wrong?
	// EnergyMethodOptions, suite bonus directories...
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data)->clone();
	stepwise_modeler_ = setup_unified_stepwise_modeler( options_, scorefxn_ );
	minimize_single_res_ = tag->getOption< bool >( "minimize_single_res", minimize_single_res_ );
	slide_docking_jumps_ = tag->getOption< bool >( "slide_docking_jumps", slide_docking_jumps_ );
	full_move_description_ = tag->getOption< std::string >( "swa_move" );
}

void
ResampleMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	//rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "minimize_single_res", xsct_rosetta_bool, "Run in single residue minimization mode", "false" );
	attlist + XMLSchemaAttribute::attribute_w_default( "slide_docking_jumps", xsct_rosetta_bool, "Slide docking jumps", "false" );
	attlist + XMLSchemaAttribute::required_attribute( "swa_move", xs_string, "String describing the SWA move to perform" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

////////////////////////////////////////////////////////////////////
std::string
ResampleMover::get_name() const {
	return mover_name();
}

std::string ResampleMover::mover_name() {
	return "ResampleMover";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResampleMover::apply( pose::Pose & pose ) {
	if ( full_move_description_ == "" ) {
		std::string move_type;
		apply( pose, move_type );
	} else {
		StepWiseMove swa_move = StepWiseMove( utility::string_split( full_move_description_ ), const_full_model_info( pose ).full_model_parameters() );
		apply( pose, swa_move );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This version randomly chooses the move element.
// This may be deprecated soon since we are moving move selection out
// to StepWiseMasterMover.
//////////////////////////////////////////////////////////////////////////
bool
ResampleMover::apply( pose::Pose & pose,
	std::string & move_type ){

	utility::vector1< StepWiseMove > swa_moves;
	swa_move_selector_->set_options( options_ );
	swa_move_selector_->get_resample_move_elements( pose, swa_moves );

	if ( swa_moves.size() == 0 ) return false;
	StepWiseMove const & swa_move = numeric::random::rg().random_element( swa_moves );

	return apply( pose, swa_move, move_type );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
ResampleMover::apply( pose::Pose & pose,
	StepWiseMove const & swa_move ){
	std::string dummy_move_type;
	return apply( pose, swa_move, dummy_move_type );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
ResampleMover::apply( pose::Pose & pose,
	StepWiseMove const & swa_move,
	std::string & move_type ){

	using namespace protocols::stepwise;
	using namespace protocols::stepwise::modeler::rna;
	using namespace protocols::stepwise::monte_carlo;
	using namespace core::pose::full_model_info;

	TR << "About to remodel move_element " << swa_move << std::endl;
	move_type = to_string( swa_move.move_type() );
	std::transform(move_type.begin(), move_type.end(), move_type.begin(), ::tolower); // this is why we love C

	// What needs to be set in StepwiseModeler:
	Size remodel_res( 0 ), remodel_suite( 0 ), cutpoint_suite( 0 );

	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
	Size const num_attachments = swa_move.attachments().size();
	bool const is_single_attachment = ( num_attachments == 1 );
	Size const move_element_size = swa_move.move_element().size();

	if ( is_single_attachment ) {
		remodel_res = get_remodel_res( swa_move, pose );
		if ( slide_docking_jumps_ && swa_move.attachment_type() == JUMP_DOCK ) slide_jump_randomly( pose, remodel_res );
	} else { // an internal residue or move_element, with two attachments.
		runtime_assert( num_attachments == 2 );
		runtime_assert( swa_move.attachments()[1].attachment_type() == BOND_TO_PREVIOUS );
		runtime_assert( swa_move.attachments()[2].attachment_type() == BOND_TO_NEXT );
		if ( move_element_size == 1 ) { // single residue
			remodel_res   = res_list.index( swa_move.moving_res() );
			remodel_suite = remodel_res - 1;
			cutpoint_suite  = res_list.index( swa_move.moving_res() );
		} else { // fixed move_element
			Size const attached_res_prev = swa_move.attachments()[1].attached_res();
			Size const attached_res_next = swa_move.attachments()[2].attached_res();
			remodel_res   = res_list.index( attached_res_prev ); // this is now the *suite*
			remodel_suite = remodel_res;
			cutpoint_suite  = res_list.index( attached_res_next ) - 1;
			runtime_assert( static_cast<int>(cutpoint_suite) == res_list.index( attached_res_next - 1 ) );
		}
	}
	runtime_assert( remodel_res > 0 );

	bool did_mutation( false );
	// based on 'n' in full_model_info.full_sequence
	if ( move_element_size == 1 ) did_mutation = mutate_res_if_allowed( pose, full_to_sub( swa_move.moving_res(), pose ) );
	bool just_min_after_mutation = ( did_mutation && ( numeric::random::rg().uniform() < options_->just_min_after_mutation_frequency() ) );

	if ( just_min_after_mutation ) {
		stepwise_modeler_->set_moving_res_and_reset( 0 );
		stepwise_modeler_->set_working_prepack_res( packer::figure_out_working_interface_res( pose, remodel_res ) );
	} else {
		stepwise_modeler_->set_moving_res_and_reset( remodel_res );
		stepwise_modeler_->set_figure_out_prepack_res( true );
	}
	utility::vector1< Size > const & moving_res = get_moving_res_from_full_model_info( pose );
	if ( ! minimize_single_res_ ) stepwise_modeler_->set_working_minimize_res( moving_res );

	if ( is_single_attachment ) {
		stepwise_modeler_->apply( pose );
	} else {
		runtime_assert( options_->allow_internal_local_moves() );
		TR << "Going to set up TRANSIENT_CUTPOINT_HANDLER with remodel_suite " << remodel_suite << " cutpoint_suite " << cutpoint_suite << std::endl;
		TransientCutpointHandler cutpoint_handler( remodel_suite, cutpoint_suite );
		if ( ! minimize_single_res_ ) cutpoint_handler.set_minimize_res( moving_res );

		cutpoint_handler.put_in_cutpoints( pose );
		stepwise_modeler_->apply( pose );
		cutpoint_handler.take_out_cutpoints( pose );
	}

	if ( did_mutation ) move_type += "-mut";
	if ( just_min_after_mutation ) move_type = "mut";

	return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
ResampleMover::get_remodel_res( StepWiseMove const & swa_move, pose::Pose const & pose ) const {
	using namespace core::pose::full_model_info;
	runtime_assert( swa_move.attachments().size() == 1 );

	// remodel res will be the first residue in the moving element that is immediatly downstream of the attachment residue.
	Size remodel_res( 0 );
	MoveElement const & move_element = swa_move.move_element();
	for ( Size n = 1; n <= move_element.size(); n++ ) {
		Size const & moving_res = full_to_sub( move_element[ n ], pose );
		Size const parent_res = pose.fold_tree().get_parent_residue( moving_res );
		if ( parent_res > 0 && sub_to_full( parent_res, pose ) == swa_move.attached_res() ) {
			remodel_res = moving_res; break;
		}
	}

	// we may have to reroot pose -- the attachment point might be 'downstream' of the moving element.
	// the rerooting will actually occur later in the Modeler.
	if ( remodel_res == 0 ) {
		Size const & moving_res = full_to_sub( swa_move.attached_res(), pose );
		if ( move_element.has_value( sub_to_full( pose.fold_tree().get_parent_residue( moving_res ), pose ) ) ) {
			remodel_res = moving_res;
		}
	}

	runtime_assert( remodel_res > 0 );
	return remodel_res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
filter_for_proximity( pose::Pose const & pose,
	utility::vector1< Size > & partition_res,
	Size const center_res ) {
	using namespace core::chemical::rna;
	runtime_assert( partition_res.has_value( center_res ) );
	static Distance const proximity_cutoff( 8.0 );
	utility::vector1< Size > filtered_partition_res;
	Vector const & center_xyz = pose.residue( center_res ).xyz( default_jump_atom( pose.residue_type( center_res ) ) );
	for ( Size n = 1; n <= partition_res.size(); n++ ) {
		Size const new_res = partition_res[ n ];
		Vector const & new_xyz =  pose.residue( new_res ).xyz( default_jump_atom( pose.residue_type( new_res ) ) );
		if ( ( new_xyz - center_xyz ).length() < proximity_cutoff ) filtered_partition_res.push_back( new_res );
	}
	partition_res = filtered_partition_res;
	runtime_assert( partition_res.has_value( center_res ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResampleMover::slide_jump_randomly( pose::Pose & pose, Size & remodel_res ) const {

	using namespace core::kinematics;
	using namespace core::chemical::rna;

	FoldTree f = pose.fold_tree();
	Size const jump_nr = f.get_jump_that_builds_residue( remodel_res );
	Size const reference_res = f.upstream_jump_residue( jump_nr );

	// check if user-inputted jump. Then can't move.
	utility::vector1< Size > const & jump_res_map = const_full_model_info( pose ).jump_res_map();
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
	if ( jump_res_map[ res_list[ remodel_res ]  ] > 0 &&
			jump_res_map[ res_list[ reference_res ] ] > 0 ) {
		runtime_assert( jump_res_map[ res_list[ remodel_res ] ] == jump_res_map[ res_list[ reference_res ] ] );
		return;
	}
	utility::vector1< Size > root_partition_res, moving_partition_res;
	figure_out_root_and_moving_partition_res( pose, remodel_res, root_partition_res, moving_partition_res );

	if ( options_->local_redock_only() ) {
		filter_for_proximity( pose, root_partition_res, reference_res );
		filter_for_proximity( pose, moving_partition_res, remodel_res );
	}

	// need to make sure JUMP_DOCK remain in different chains!
	utility::vector1< Size > dock_domain_map = figure_out_dock_domain_map_from_full_model_info_const( pose );
	utility::vector1< std::pair< Size, Size > > possible_jump_pairs;
	for ( Size i = 1; i <= root_partition_res.size(); i++ ) {
		Size const & root_res = root_partition_res[ i ];
		for ( Size j = 1; j <= moving_partition_res.size(); j++ ) {
			Size const & move_res = moving_partition_res[ j ];
			if ( dock_domain_map[ root_res ] != dock_domain_map[ move_res ] ) possible_jump_pairs.push_back( std::make_pair( root_res, move_res ) );
		}
	}
	std::pair< Size, Size > const new_jump_pair = numeric::random::rg().random_element( possible_jump_pairs );
	Size const new_reference_res = new_jump_pair.first;
	Size const new_remodel_res   = new_jump_pair.second;

	f.slide_jump( jump_nr, new_reference_res, new_remodel_res );
	f.set_jump_atoms( jump_nr, default_jump_atom( pose.residue_type( new_reference_res ) ),
		default_jump_atom( pose.residue_type( new_remodel_res ) ) );

	pose.fold_tree( f );
	TR << "Slid jump from: " << remodel_res << " -- " << reference_res <<
		" to: " << new_remodel_res << " -- " << new_reference_res << std::endl;

	remodel_res = new_remodel_res;
}

std::string ResampleMoverCreator::keyname() const {
	return ResampleMover::mover_name();
}

protocols::moves::MoverOP
ResampleMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ResampleMover );
}

void ResampleMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResampleMover::provide_xml_schema( xsd );
}

} //mover
} //monte_carlo
} //stepwise
} //protocols
