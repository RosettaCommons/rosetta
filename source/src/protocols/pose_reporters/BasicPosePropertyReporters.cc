// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/pose_reporters/BasicPoseReporters.hh
/// @brief  Collection of simple pose property reporters
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_pose_reporters_BasicPoseReporters_CC
#define INCLUDED_protocols_pose_reporters_BasicPoseReporters_CC

// Unit Headers
#include <protocols/rosetta_scripts/PosePropertyReporter.hh>
#include <protocols/pose_reporters/BasicPosePropertyReporters.hh>
#include <protocols/pose_reporters/BasicPosePropertyReporterCreators.hh>
#include <protocols/rosetta_scripts/PosePropertyReporterFactory.hh>
// Package headers
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/filters/FilterFactory.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/rms_util.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
// C++ Headers
#include <string>

static THREAD_LOCAL basic::Tracer TR( "protocols.pose_reporters.BasicPosePropertyReporters" );

namespace protocols {
namespace pose_reporters {

////////////////////////////////////////////////////////////////////////
// EnergyReporter

// Creator
protocols::rosetta_scripts::PosePropertyReporterOP EnergyReporterCreator::create_reporter() const {
	return protocols::rosetta_scripts::PosePropertyReporterOP( new EnergyReporter() );
}
void
EnergyReporterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	EnergyReporter::provide_xml_schema( xsd );
}

// Reporter
EnergyReporter::EnergyReporter() :
	scorefxn_(/* NULL */),
	scoretype_(core::scoring::dummy_score_type)
{
}

void
EnergyReporter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "scorefunction", xs_string, "Name of score function weights to use" )
		+ XMLSchemaAttribute( "term", xs_string, "Score term to evaluate" );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( name() )
		.description( "XRW TO DO" )
		.add_attributes( attlist )
		.add_optional_name_attribute()
		.complex_type_naming_func( & rosetta_scripts::PosePropertyReporterFactory::complex_type_name_for_pose_reporter )
		.write_complex_type_to_schema( xsd );
}

core::Real EnergyReporter::report_property( core::pose::Pose & p ) const
{
	core::Real r = 0;

	if ( scoretype_ < core::scoring::dummy_score_type ) {
		r = scorefxn_->score_by_scoretype(p, scoretype_);
	} else {
		r = scorefxn_->score(p);
	}

	return r;
}

void EnergyReporter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
)
{
	scorefxn_ = tag->hasOption("scorefunction") ?
		core::scoring::ScoreFunctionFactory::create_score_function( tag->getOption<std::string>("scorefunction") ) :
		core::scoring::get_score_function();

	if ( tag->hasOption("term") ) {
		scoretype_ = core::scoring::ScoreTypeManager::score_type_from_name(tag->getOption<std::string>("term"));
	}
}


////////////////////////////////////////////////////////////////////////
// FilterReporter

// Creator
protocols::rosetta_scripts::PosePropertyReporterOP FilterReporterCreator::create_reporter() const {
	return protocols::rosetta_scripts::PosePropertyReporterOP( new FilterReporter() );
}

void
FilterReporterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	FilterReporter::provide_xml_schema( xsd );
}

// Reporter
FilterReporter::FilterReporter() :
	filter_(/* NULL */)
{
}

core::Real FilterReporter::report_property( core::pose::Pose & pose ) const
{
	core::Real r = 0;

	if ( filter_ ) {
		r = filter_->report_sm(pose);
		setPoseExtraScore(pose, filter_->get_user_defined_name(), (float)r);
	} else {
		TR << "No filter instance; cannot score pose!" << std::endl;
	}

	return r;
}


void
FilterReporter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	protocols::filters::FilterFactory::get_instance()->define_filter_xml_schema( xsd );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "filter", xs_string, "Name attribute of filter defined previously in the RosettaScript" );

	XMLSchemaSimpleSubelementList subelements;
	subelements.add_group_subelement( & filters::FilterFactory::filter_xml_schema_group_name );


	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( name() )
		.description( "XRW TO DO" )
		.set_subelements_repeatable( subelements, 0, 1 )
		.add_attributes( attlist )
		.add_optional_name_attribute()
		.complex_type_naming_func( & rosetta_scripts::PosePropertyReporterFactory::complex_type_name_for_pose_reporter )
		.write_complex_type_to_schema( xsd );
}


void FilterReporter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
)
{
	using namespace utility::tag;
	using namespace protocols::rosetta_scripts;

	TagCOP filter_tag(nullptr);

	if ( tag->hasOption("filter") ) {
		// Find a filter by name defined somewhere upstream in the script
		std::string filter_name( tag->getOption<std::string>("filter") );
		RosettaScriptsParser parser;
		filter_tag = parser.find_rosettascript_tag(
			tag,
			"FILTERS",
			"name",
			filter_name
		);

		if ( !filter_tag ) {
			throw utility::excn::EXCN_RosettaScriptsOption("Cannot find filter named \"" + filter_name + "\"");
		}

	} else {
		// Filter is defined inline (first child tag)
		utility::vector0< TagCOP > tags( tag->getTags() );
		for ( utility::vector0< TagCOP >::const_iterator it = tags.begin(); it != tags.end(); ++it ) {
			filter_tag = *it;
			break;
		}
	}

	if ( filter_tag ) {
		filter_  = protocols::filters::FilterFactory::get_instance()->newFilter( filter_tag, data, filters, movers, pose );
	}

	if ( !filter_ ) {
		std::ostringstream s;
		s << "Cannot create filter from script tag: " << tag;
		throw utility::excn::EXCN_RosettaScriptsOption(s.str());
	}
}

////////////////////////////////////////////////////////////////////////
// RMSDReporter

// Creator
protocols::rosetta_scripts::PosePropertyReporterOP RMSDReporterCreator::create_reporter() const {
	return protocols::rosetta_scripts::PosePropertyReporterOP( new RMSDReporter() );
}

void
RMSDReporterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	RMSDReporter::provide_xml_schema( xsd );
}


// Reporter
RMSDReporter::RMSDReporter() :
	mode_(MODE_NONE)
{
}

/// @brief Report the RMSD between poses p1 and p2
/// Depending on the mode, this can be either all atom or CA RMSD, optionally with residue selection
core::Real RMSDReporter::report_property( core::pose::Pose & p1, core::pose::Pose & p2 ) const
{
	core::Real r = 0;

	switch(mode_) {
	case MODE_CA_rmsd :
		if ( residues_.size() > 0 ) {
			r = core::scoring::CA_rmsd( p1, p2, residues_ );
		} else {
			r = core::scoring::CA_rmsd( p1, p2 );
		}
		break;

	case MODE_all_atom_rmsd :
		if ( residues_.size() > 0 ) {
			r = core::scoring::all_atom_rmsd( p1, p2, residues_ );
		} else {
			r = core::scoring::all_atom_rmsd( p1, p2 );
		}
		break;

	default :
		TR.Warning << "Unknown RMSD report mode! This should not happen." << std::endl;
	}

	TR.Debug << p1.sequence() << " :: " << p2.sequence() << " RMSD: " << r << std::endl;

	return r;
}


void
RMSDReporter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	XMLSchemaRestriction rmsd_mode;
	rmsd_mode.name( "rmsd_mode" );
	rmsd_mode.base_type( xs_string );
	rmsd_mode.add_restriction( xsr_enumeration, "CA" );
	rmsd_mode.add_restriction( xsr_enumeration, "all_atom" );
	xsd.add_top_level_element( rmsd_mode );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "mode", "rmsd_mode", "Calculate CA or all_atom RMSD?" );

	core::pose::attributes_for_get_resnum_list( attlist, xsd, "residues" );
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( name() )
		.description( "Reports the Calpha or full-atom RMSD of a pose (optionally within a specified residue range)" )
		.add_attributes( attlist )
		.add_optional_name_attribute()
		.complex_type_naming_func( & rosetta_scripts::PosePropertyReporterFactory::complex_type_name_for_pose_reporter )
		.write_complex_type_to_schema( xsd );
}




void RMSDReporter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /* data */,
	protocols::filters::Filters_map & /* filters */,
	protocols::moves::Movers_map const & /* movers */,
	core::pose::Pose const & pose
)
{
	using namespace utility::tag;

	if ( tag->hasOption("mode") ) {
		std::string mode = tag->getOption<std::string>("mode");
		if ( mode == "CA" ) {
			mode_ = MODE_CA_rmsd;
		} else if ( mode == "all_atom" ) {
			mode_ = MODE_all_atom_rmsd;
		} else {
			throw utility::excn::EXCN_RosettaScriptsOption("Unknown RMSD reporter mode: " + mode);
		}
	}

	if ( tag->hasOption("residues") ) {
		std::string residues = tag->getOption<std::string>("residues");
		// Why do we have std::set in some places and std::list elsewhere?
		std::set< core::Size > residues_list ( core::pose::get_resnum_list(residues, pose) );
		residues_.assign( residues_list.begin(), residues_list.end() );
	}

	if ( mode_ == MODE_NONE ) {
		throw utility::excn::EXCN_RosettaScriptsOption("RMSD reporter mode not specified. Choose form: CA, all_atom");
	}
}

////////////////////////////////////////////////////////////////////////


} // pose_reporters
} // protocols

#endif //INCLUDED_protocols_pose_reporters_BasicPoseReporters_CC
