// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/drug_design/RDKitMetricsMover.cc
/// @brief Add RDKit Metrics to a pose
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit Headers
#include <protocols/drug_design/RDKitMetricsMover.hh>
#include <protocols/drug_design/RDKitMetricsMoverCreator.hh>

// Package Headers
#include <protocols/moves/Mover.hh>
#include <core/chemical/rdkit/util.hh>
#include <core/chemical/rdkit/RestypeToRDMol.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <protocols/rosetta_scripts/util.hh>

// Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Debugging output
#include <utility/io/ozstream.hh>

// External headers

// C/C++ headers
#include <string>

namespace protocols {
namespace drug_design {

static basic::Tracer TR("protocols.drug_design.RDKitMetricsMover");

std::string
RDKitMetricsMoverCreator::keyname() const
{
	return RDKitMetricsMover::mover_name();
}

protocols::moves::MoverOP
RDKitMetricsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RDKitMetricsMover );
}

void RDKitMetricsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RDKitMetricsMover::provide_xml_schema( xsd );
}

/// @brief default constructor
RDKitMetricsMover::RDKitMetricsMover() {}

/// @brief destructor
RDKitMetricsMover::~RDKitMetricsMover() {}

/// @brief clone this object
protocols::moves::MoverOP
RDKitMetricsMover::clone() const
{
	return protocols::moves::MoverOP( new RDKitMetricsMover( *this ) );
}

/// @brief create this type of object
protocols::moves::MoverOP
RDKitMetricsMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new RDKitMetricsMover() );
}

void
RDKitMetricsMover::apply( Pose & pose ) {
	if ( chains_.size() == 0 ) {
		utility_exit_with_message("Must set chains for RDKitMetricsMover!");
	}

	for ( core::Size ii(1); ii <= chains_.size(); ++ii ) {
		if ( chains_[ii].size() != 1 || ! core::pose::has_chain(chains_[ii], pose) ) {
			utility_exit_with_message("In RDKitMetricsMover, the Pose doesn't have a chain '"+chains_[ii]+"'!");
		}
		add_scores_for_chain( pose, chains_[ii] );
	}
}// apply

void
RDKitMetricsMover::add_scores_for_chain( Pose & pose, std::string const & chain ) const {
	if ( chain.size() != 1 ) {
		utility_exit_with_message("In RDKitMetricsMover, chain designation '"+chain+"' shoud be a single letter.");
	}
	for ( core::Size ii: core::pose::get_resnums_for_chain( pose, chain[0] ) ) {
		add_scores_for_residue( pose, ii);
	}
}

void
RDKitMetricsMover::add_scores_for_residue( Pose & pose, core::Size resid ) const {
	std::string tag( get_tag( pose, resid ) );

	// Map of name:description pairs
	std::map< std::string, std::string > const metrics( core::chemical::rdkit::get_metric_names() );

	// Use a neutral molecule without hydrogens for calculating the metric values - should be the default
	core::chemical::MutableResidueTypeOP restype( utility::pointer::make_shared< core::chemical::MutableResidueType >( pose.residue_type(resid) ) );
	core::chemical::rdkit::RestypeToRDMol converter(*restype);
	::RDKit::RWMolOP rdmol(converter.Mol() );

	for ( std::map< std::string, std::string >::const_iterator itr( metrics.begin() ), itr_end( metrics.end() );
			itr != itr_end; ++itr ) {
		std::string const & metric( itr->first );
		core::Real metric_value( core::chemical::rdkit::rdkit_metric( *rdmol, metric ) );
		// "lm" here is for ligand metric (in parallel to the "if" interface tag.
		core::pose::setPoseExtraScore(pose, "lm_"+tag+"_"+metric, metric_value );
	}

}

std::string
RDKitMetricsMover::get_tag( Pose & pose, core::Size resid ) const {
	core::pose::PDBInfoCOP pdb_info( pose.pdb_info() );
	std::stringstream tag;
	tag << pdb_info->chain(resid) << pdb_info->number(resid);
	if ( pdb_info->icode(resid) && pdb_info->icode(resid) != ' ' ) {
		tag << pdb_info->icode(resid);
	}
	return tag.str();
}

std::string
RDKitMetricsMover::get_name() const {
	return mover_name();
}

std::string
RDKitMetricsMover::mover_name() {
	return "RDKitMetricsMover";
}

void RDKitMetricsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "chains", xsct_chain_cslist, "The chain of the ligand to compute the metrics for" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Add RDKit metrics for a given ligand to the Pose's extra scores.",
		attlist );

}

/// @brief parse xml file
void
RDKitMetricsMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & )
{
	std::string const chains_str = tag->getOption<std::string>("chains");
	chains( utility::string_split(chains_str, ',') );
}

} // ns drug_design
} // ns protocols
