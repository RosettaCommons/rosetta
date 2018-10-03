// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/burial_metrics/PolarGroupBurialPyMolStringMetric.cc
/// @brief A string metric that generates a string of PyMol commands to colour a structure's polar
/// groups based on burial.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Unit headers
#include <protocols/analysis/burial_metrics/PolarGroupBurialPyMolStringMetric.hh>
#include <protocols/analysis/burial_metrics/PolarGroupBurialPyMolStringMetricCreator.hh>

// Core headers
#include <core/simple_metrics/StringMetric.hh>
#include <core/simple_metrics/util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/BuriedUnsatPenalty.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/ResidueLevelTask_.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

//Protocols headers
#include <protocols/rosetta_scripts/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/pointer/memory.hh>


static basic::Tracer TR( "protocols.analysis.burial_metrics.PolarGroupBurialPyMolStringMetric" );


namespace protocols {
namespace analysis {
namespace burial_metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
/// @details Initializes ScoreFunction to nullptr.
PolarGroupBurialPyMolStringMetric::PolarGroupBurialPyMolStringMetric():
	core::simple_metrics::StringMetric(),
	sfxn_(nullptr),
	verbose_(false)
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PolarGroupBurialPyMolStringMetric::~PolarGroupBurialPyMolStringMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
/// @details Copies, but does not clone, input scorefunction.
PolarGroupBurialPyMolStringMetric::PolarGroupBurialPyMolStringMetric( PolarGroupBurialPyMolStringMetric const & src ):
	core::simple_metrics::StringMetric( src ),
	sfxn_(src.sfxn_),
	verbose_(src.verbose_)
{

}

core::simple_metrics::SimpleMetricOP
PolarGroupBurialPyMolStringMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new PolarGroupBurialPyMolStringMetric( *this ) );

}

/// @brief Name of the class
/// @details Calls name_static().
std::string
PolarGroupBurialPyMolStringMetric::name() const {
	return name_static();
}

///@brief Name of the class for creator.
std::string
PolarGroupBurialPyMolStringMetric::name_static() {
	return "PolarGroupBurialPyMolStringMetric";

}

///@brief Name (descriptive) of the metric.
std::string
PolarGroupBurialPyMolStringMetric::metric() const {
	return "polar_group_burial_pymol";
}

/// @brief Given an XML tag, configure an instance of this class.
void
PolarGroupBurialPyMolStringMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap  )
{
	SimpleMetric::parse_base_tag( tag );
	set_verbose( tag->getOption<bool>("verbose", false) );
	set_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, "scorefxn", datamap ) );
}

void
PolarGroupBurialPyMolStringMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(  "verbose", xsct_rosetta_bool, "If true, this metric outputs to tracer as well as to the pose.  False by default",  "false"  );

	protocols::rosetta_scripts::attributes_for_parse_score_function_w_description_when_required( attlist, "scorefxn", "A scorefunction, defined previously in the RosettaScripts XML, which will be used to extract the buried_unsaturated_penalty's definition of burial.  Required option." );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A string metric that generates a string of PyMol commands to colour a structure's polar groups based on burial.", attlist);
}

std::string
PolarGroupBurialPyMolStringMetric::calculate(
	const core::pose::Pose &original_pose
) const {
	runtime_assert_string_msg(sfxn_ != nullptr, "Error in PolarGroupBurialPyMolStringMetric::calculate(): A scorefunction must be passed to this instance of this class before it can be used to calculate metrics!");

	core::pose::Pose pose( original_pose ); //Copy the pose.
	bool const is_symmetric( core::pose::symmetry::is_symmetric(pose) );

	if ( pose.energies().data().has( core::scoring::EnergiesCacheableDataType::BURIED_UNSAT_HBOND_GRAPH ) ) {
		pose.energies().data().clear( core::scoring::EnergiesCacheableDataType::BURIED_UNSAT_HBOND_GRAPH ); //To prevent reuse.
	}

	core::pack::guidance_scoreterms::buried_unsat_penalty::BuriedUnsatPenalty scoreterm( sfxn_->energy_method_options() );
	core::scoring::ScoreFunctionOP sfxn( is_symmetric ? utility::pointer::make_shared< core::scoring::symmetry::SymmetricScoreFunction >() : utility::pointer::make_shared< core::scoring::ScoreFunction >() );
	sfxn->set_energy_method_options( sfxn_->energy_method_options() );
	sfxn->set_weight( core::scoring::buried_unsatisfied_penalty, 1.0 );
	core::pack::task::PackerTaskOP task( new core::pack::task::PackerTask_(pose) );
	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(i))).prevent_repacking();
	}

	core::pack::task::PackerTaskOP possibly_symmetric_task(nullptr);
	if ( is_symmetric ) {
		task->request_symmetrize_by_intersection();
		possibly_symmetric_task = core::pack::make_new_symmetric_PackerTask_by_requested_method( pose, task );
	} else {
		possibly_symmetric_task = task;
	}
	debug_assert( possibly_symmetric_task != nullptr );

	core::pack::rotamer_set::RotamerSetsOP rotsets( nullptr );
	if ( is_symmetric ) {
		rotsets = utility::pointer::make_shared< core::pack::rotamer_set::symmetry::SymmetricRotamerSets >();
	} else {
		rotsets = utility::pointer::make_shared< core::pack::rotamer_set::RotamerSets >();
	}
	debug_assert( rotsets != nullptr );
	rotsets->set_task( possibly_symmetric_task );
	rotsets->initialize_pose_for_rotsets_creation(pose);

	utility::graph::GraphOP packer_neighbor_graph( core::pack::create_packer_graph( pose, *sfxn, possibly_symmetric_task ) );
	rotsets->build_rotamers(pose, *sfxn, packer_neighbor_graph);
	rotsets->prepare_sets_for_packing(pose, *sfxn);
	scoreterm.set_prevent_pruning(true);
	scoreterm.set_up_residuearrayannealableenergy_for_packing( pose, *rotsets, *sfxn );

	std::stringstream outstream;
	scoreterm.provide_pymol_commands_to_show_groups( outstream, pose );

	if ( verbose_ ) TR << "\n" << outstream.str() << std::endl;

	return outstream.str();
}

/// @brief Set the scorefunction.
/// @details Used for definition of burial.  Input pointer is copied without cloning.
void
PolarGroupBurialPyMolStringMetric::set_scorefxn(
	core::scoring::ScoreFunctionCOP sfxn_in
) {
	runtime_assert_string_msg( sfxn_in != nullptr, "Error in PolarGroupBurialPyMolStringMetric::set_scorefxn(): A null pointer was passed to this function!" );
	sfxn_ = sfxn_in;
}

/// @brief Set whether this metric should also output to tracer when evaluated.
void
PolarGroupBurialPyMolStringMetric::set_verbose(
	bool const setting
) {
	verbose_  = setting;
}

void
PolarGroupBurialPyMolStringMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PolarGroupBurialPyMolStringMetric::provide_xml_schema( xsd );
}

std::string
PolarGroupBurialPyMolStringMetricCreator::keyname() const {
	return PolarGroupBurialPyMolStringMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
PolarGroupBurialPyMolStringMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new PolarGroupBurialPyMolStringMetric );

}

} //protocols
} //analysis
} //burial_metrics






