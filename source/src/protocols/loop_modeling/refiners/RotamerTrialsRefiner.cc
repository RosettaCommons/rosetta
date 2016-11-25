// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Unit headers
#include <protocols/loop_modeling/refiners/RotamerTrialsRefiner.hh>
#include <protocols/loop_modeling/refiners/RotamerTrialsRefinerCreator.hh>
#include <protocols/loop_modeling/refiners/packing_helper.hh>
#include <protocols/loop_modeling/utilities/rosetta_scripts.hh>

// Core headers
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// RosettaScripts headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace loop_modeling {
namespace refiners {

using namespace std;
using core::pose::Pose;
using core::pack::task::TaskFactoryOP;
using core::scoring::ScoreFunctionOP;

// XRW TEMP protocols::moves::MoverOP RotamerTrialsRefinerCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new RotamerTrialsRefiner );
// XRW TEMP }

// XRW TEMP std::string RotamerTrialsRefinerCreator::keyname() const {
// XRW TEMP  return "RotamerTrialsRefiner";
// XRW TEMP }

void RotamerTrialsRefiner::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose) {

	LoopMover::parse_my_tag(tag, data, filters, movers, pose);
	utilities::set_scorefxn_from_tag(*this, tag, data);
	utilities::set_task_factory_from_tag(*this, tag, data);
}

bool RotamerTrialsRefiner::do_apply(Pose & pose) {
	return packing_helper(pose, this, core::pack::rotamer_trials);
}

ScoreFunctionOP RotamerTrialsRefiner::get_score_function() {
	return get_tool<ScoreFunctionOP>(ToolboxKeys::SCOREFXN);
}

void RotamerTrialsRefiner::set_score_function(ScoreFunctionOP score_function) {
	set_tool(ToolboxKeys::SCOREFXN, score_function);
}

TaskFactoryOP RotamerTrialsRefiner::get_task_factory() {
	return get_tool<TaskFactoryOP>(ToolboxKeys::TASK_FACTORY);
}

TaskFactoryOP RotamerTrialsRefiner::get_task_factory(TaskFactoryOP fallback) {
	return get_tool<TaskFactoryOP>(ToolboxKeys::TASK_FACTORY, fallback);
}

void RotamerTrialsRefiner::set_task_factory(TaskFactoryOP task_factory) {
	set_tool(ToolboxKeys::TASK_FACTORY, task_factory);
}

std::string RotamerTrialsRefiner::get_name() const {
	return mover_name();
}

std::string RotamerTrialsRefiner::mover_name() {
	return "RotamerTrialsRefiner";
}

void RotamerTrialsRefiner::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	// add score function tags
	utilities::attributes_for_set_scorefxn_from_tag( attlist );
	// add task operation tags
	utilities::attributes_for_set_task_factory_from_tag( attlist );
	// add LoopMover base class tags
	XMLSchemaSimpleSubelementList subelement_list;
	XMLSchemaComplexTypeGenerator ct_gen;
	LoopMover::define_composition_schema( xsd, ct_gen, subelement_list );
	ct_gen.element_name( mover_name() )
		.description(
		"Use rotamer trials to quickly optimize the sidechains in and around the loop being sampled."
		"Rotamer trials is a greedy algorithm for packing sidechains. Each residue is considered only"
		" once, and is placed in its optimal rotamer" )
		.add_attributes( attlist  )
		.write_complex_type_to_schema( xsd );
}

std::string RotamerTrialsRefinerCreator::keyname() const {
	return RotamerTrialsRefiner::mover_name();
}

protocols::moves::MoverOP
RotamerTrialsRefinerCreator::create_mover() const {
	return protocols::moves::MoverOP( new RotamerTrialsRefiner );
}

void RotamerTrialsRefinerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RotamerTrialsRefiner::provide_xml_schema( xsd );
}


}
}
}
