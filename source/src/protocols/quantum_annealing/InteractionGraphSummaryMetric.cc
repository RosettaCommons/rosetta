// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/quantum_annealing/InteractionGraphSummaryMetric.cc
/// @brief A simple metric that allows an InteractionGraph to be written out in a format that external annealers (including quantum annealers) can read.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Unit headers
#include <protocols/quantum_annealing/InteractionGraphSummaryMetric.hh>
#include <protocols/quantum_annealing/InteractionGraphSummaryMetricCreator.hh>

// Core headers
#include <core/simple_metrics/StringMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSetsFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/interaction_graph/util.hh>

// Protocols includes
#include <protocols/rosetta_scripts/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "protocols.quantum_annealing.InteractionGraphSummaryMetric" );


namespace protocols {
namespace quantum_annealing {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
InteractionGraphSummaryMetric::InteractionGraphSummaryMetric():
	core::simple_metrics::StringMetric(),
	task_factory_( nullptr ),
	scorefxn_( nullptr ),
	short_version_(false)
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
InteractionGraphSummaryMetric::~InteractionGraphSummaryMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
InteractionGraphSummaryMetric::InteractionGraphSummaryMetric( InteractionGraphSummaryMetric const & ) = default;

core::simple_metrics::SimpleMetricOP
InteractionGraphSummaryMetric::clone() const {
	return utility::pointer::dynamic_pointer_cast< core::simple_metrics::SimpleMetric >( utility::pointer::make_shared<InteractionGraphSummaryMetric>( *this ) );
}

std::string
InteractionGraphSummaryMetric::name() const {
	return name_static();
}

std::string
InteractionGraphSummaryMetric::name_static() {
	return "InteractionGraphSummaryMetric";

}
std::string
InteractionGraphSummaryMetric::metric() const {

	return "IG Summary";
}

/// @brief Set the task factory to be used for interaxction graph setup.
/// @details Doesn't clone input.
void
InteractionGraphSummaryMetric::set_task_factory(
	core::pack::task::TaskFactoryCOP task_factory_in
) {
	task_factory_ = task_factory_in;
}

/// @brief Set the scorefunction used for scoring.
/// @details Copied, not cloned.
void
InteractionGraphSummaryMetric::set_scorefunction(
	core::scoring::ScoreFunctionCOP sfxn_in
) {
	runtime_assert_string_msg( sfxn_in != nullptr, "Error in protocols::quantum_annealing::InteractionGraphSummaryMetric::set_scorefunction(): The scorefunction pointer provided to this function cannot be null!" );
	scorefxn_ = sfxn_in;
}


void
InteractionGraphSummaryMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	SimpleMetric::parse_base_tag( tag );
	core::pack::task::TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory != nullptr ) set_task_factory( new_task_factory );

	core::scoring::ScoreFunctionOP new_sfxn( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );
	runtime_assert_string_msg( new_sfxn != nullptr, "Error in protocols::quantum_annealing::InteractionGraphSummaryMetric::parse_my_tag(): Could not parse scorefunction." );
	set_scorefunction( new_sfxn );

	set_short_version( tag->getOption<bool>("skip_pose_reconstruction_info", false) );
}

void
InteractionGraphSummaryMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	rosetta_scripts::attributes_for_parse_score_function( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default ( "skip_pose_reconstruction_info", xsct_rosetta_bool, "If true, then this metric only stores a summary of the interaction graph.  If false (the default), then it stores both the interaction graph summary and full information for reconstructing the pose.  False by default.", "false" );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for writing out a description of the Rosetta-calculated interaction graph, which can be read in by external annealers or optimizers.", attlist);
}

std::string
InteractionGraphSummaryMetric::calculate(const core::pose::Pose & pose ) const {
	std::string const errmsg( "Error in protocols::quantum_annealing::InteractionGraphSummaryMetric::calculate(): " );

	//Check for scorefunction:
	runtime_assert_string_msg( scorefxn_ != nullptr, errmsg + "Scoring function was not specified before this function was called." );

	//Unfortunately, we have to clone the pose, since we need a nonconst pose.  This is not the most expensive thing
	//that we'll be doing, however.
	core::pose::PoseOP pose_copy( pose.clone() );
	pose_copy->update_residue_neighbors();

	core::pack::task::PackerTaskCOP my_task;
	if ( task_factory_ == nullptr ) {
		TR << "No task factory or task operations specified,  Creating default PackerTask." << std::endl;
		my_task = core::pack::task::TaskFactory::create_packer_task( *pose_copy );
	} else {
		my_task = task_factory_->create_task_and_apply_taskoperations( *pose_copy );
	}

	if ( core::pose::symmetry::is_symmetric( *pose_copy ) ) {
		my_task = core::pack::make_symmetric_task( *pose_copy, my_task );
	}

	core::pack::rotamer_set::RotamerSetsOP rot_sets( core::pack::rotamer_set::RotamerSetsFactory::create_rotamer_sets( *pose_copy ) );

	core::pack::interaction_graph::AnnealableGraphBaseOP intgraph;

	core::pack::pack_rotamers_setup( *pose_copy, *scorefxn_, my_task, rot_sets, intgraph );

	std::stringstream igstream;

	core::pack::interaction_graph::get_annealable_graph_summary( igstream, *pose_copy, *rot_sets, intgraph, short_version() );

	return igstream.str();
}

void
InteractionGraphSummaryMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	InteractionGraphSummaryMetric::provide_xml_schema( xsd );
}

std::string
InteractionGraphSummaryMetricCreator::keyname() const {
	return InteractionGraphSummaryMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
InteractionGraphSummaryMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new InteractionGraphSummaryMetric );

}

} //protocols
} //quantum_annealing

#ifdef    SERIALIZATION



template< class Archive >
void
protocols::quantum_annealing::InteractionGraphSummaryMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::StringMetric>( this ) );
	arc( CEREAL_NVP( task_factory_ ) );
	arc( CEREAL_NVP( scorefxn_ ) );
	arc( CEREAL_NVP( short_version_ ) ); //bool
}

template< class Archive >
void
protocols::quantum_annealing::InteractionGraphSummaryMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::StringMetric >( this ) );

	core::pack::task::TaskFactoryOP local_task;
	arc( local_task );
	task_factory_ = local_task; // copy the non-const pointer(s) into the const pointer(s)

	core::scoring::ScoreFunctionOP local_scorefxn;
	arc( local_scorefxn );
	scorefxn_ = local_scorefxn; //Again, non-const to const.

	arc( short_version_ ); //bool
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::quantum_annealing::InteractionGraphSummaryMetric );
CEREAL_REGISTER_TYPE( protocols::quantum_annealing::InteractionGraphSummaryMetric )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_quantum_annealing_InteractionGraphSummaryMetric )
#endif // SERIALIZATION





