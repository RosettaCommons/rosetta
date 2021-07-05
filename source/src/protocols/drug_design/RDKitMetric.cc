// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/drug_design/RDKitMetric.cc
/// @brief A SimpleMetric which measures properties calcualted by RDKit on a ligand.
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <protocols/drug_design/RDKitMetric.hh>
#include <protocols/drug_design/RDKitMetricCreator.hh>

#include <core/chemical/rdkit/RestypeToRDMol.hh>
#include <core/chemical/rdkit/util.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "protocols.drug_design.RDKitMetric" );


namespace protocols {
namespace drug_design {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
RDKitMetric::RDKitMetric():
	core::simple_metrics::RealMetric()
{}

RDKitMetric::RDKitMetric( core::select::residue_selector::ResidueSelectorCOP residue, std::string const & metric_name ):
	core::simple_metrics::RealMetric(),
	residue_( residue )
{
	rdkit_metric( metric_name );
}

core::simple_metrics::SimpleMetricOP
RDKitMetric::clone() const {
	return utility::pointer::make_shared< RDKitMetric >( *this );
}

void
RDKitMetric::rdkit_metric( std::string const & setting ) {
	if ( setting.size() != 0 ) {
		std::map< std::string, std::string > const metrics( core::chemical::rdkit::get_metric_names() );
		if ( metrics.count( setting ) == 0 ) {
			TR.Error << "Metric " << setting << " not recognized by Rosetta as a valid RDKit metric. Valid metrics are:\n\n";
			for ( std::map< std::string, std::string >::const_iterator itr( metrics.begin() ), itr_end( metrics.end());
					itr != itr_end; ++itr ) {
				TR.Error << itr->first << "\t-\t" << itr->second << "\n";
			}
			TR << std::endl;
			utility_exit_with_message("Rosetta doesn't understand '"+setting+"' as an RDKit metric.");
		}
	}
	rdkit_metric_ = setting;
}

std::string
RDKitMetric::name() const {
	return name_static();
}

std::string
RDKitMetric::name_static() {
	return "RDKitMetric";

}
std::string
RDKitMetric::metric() const {
	return rdkit_metric_;
}

void
RDKitMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	SimpleMetric::parse_base_tag( tag );

	residue_selector( core::select::residue_selector::parse_residue_selector( tag, datamap ) );
	// Residue setting
	rdkit_metric( tag->getOption< std::string >("metric_name") );
}

void
RDKitMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute("metric_name", xs_string, "The RDKit metric to calculate.");

	attributes_for_parse_residue_selector( attlist, "residue_selector",
		"Selector to pick which residues to use." );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"Use RDKit to calculate metric values for a particular ResidueType.", attlist);
}

core::Real
RDKitMetric::calculate(const core::pose::Pose & pose ) const {
	if ( rdkit_metric_.size() == 0 ) {
		utility_exit_with_message("Must set metric in RDKitMetric.");
	}
	if ( residue_ == nullptr ) {
		utility_exit_with_message("Must set residue for RDKitMetric.");
	}
	utility::vector1< core::Size > resnums = core::select::get_residues_from_subset( residue_->apply( pose ) );
	if ( resnums.size() != 1 ) {
		TR.Error << "Can only apply RDKitMetric to a single residue. Got " << resnums.size() << " residues instead." << std::endl;
		utility_exit_with_message("Can only apply RDKitMetric to one residue.");
	}
	core::Size resnum( resnums[1] );

	core::chemical::MutableResidueTypeOP restype( utility::pointer::make_shared< core::chemical::MutableResidueType >(pose.residue(resnum).type()) );
	TR << "Calculating RDKit metric '" << rdkit_metric_ << "' value for residue " << resnum << ", of type " << restype->name() << std::endl;

	core::Real retval(0);

	// Use a neutral, hydrogen free molecule for calculating the metric values - should be the default
	core::chemical::rdkit::RestypeToRDMol converter(*restype);
	::RDKit::RWMolOP rdmol(converter.Mol() );
	retval = core::chemical::rdkit::rdkit_metric( *rdmol, rdkit_metric_ );

	return retval;
}


void
RDKitMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	RDKitMetric::provide_xml_schema( xsd );
}

std::string
RDKitMetricCreator::keyname() const {
	return RDKitMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
RDKitMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< RDKitMetric >();
}

} //drug_design
} //protocols


#ifdef    SERIALIZATION



template< class Archive >
void
protocols::drug_design::RDKitMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::RealMetric>( this ) );
	arc( CEREAL_NVP( residue_ ) );
	arc( CEREAL_NVP( rdkit_metric_ ) );
}

template< class Archive >
void
protocols::drug_design::RDKitMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::RealMetric >( this ) );
	arc( residue_ );
	arc( rdkit_metric_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::drug_design::RDKitMetric );
CEREAL_REGISTER_TYPE( protocols::drug_design::RDKitMetric )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_drug_design_RDKitMetric )
#endif // SERIALIZATION




