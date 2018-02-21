// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/metrics/SequenceMetric.cc
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Unit headers
#include <protocols/multistage_rosetta_scripts/cluster/metrics/SequenceMetric.hh>
#include <protocols/multistage_rosetta_scripts/cluster/metrics/SequenceMetricCreator.hh>
#include <protocols/multistage_rosetta_scripts/cluster/util.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/XMLSchemaValidation.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>
#include <utility/pointer/memory.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace multistage_rosetta_scripts {
namespace cluster {
namespace metrics {

SequenceMetric::SequenceMetric() :
	sequence_( "" )
{}

SequenceMetric::SequenceMetric( std::string sequence ) :
	sequence_( std::move( sequence) )
{}

SequenceMetric::~SequenceMetric(){}

platform::Real
SequenceMetric::distance( SequenceMetric const & other ) const {
	platform::Size const min_size =
		( other.sequence().length() < sequence_.length() ? other.sequence().length() : sequence_.length() );

	platform::Size num_mutations = 0;

	for ( platform::Size i = 0; i < min_size; ++i ) {
		if ( sequence_[ i ] != other.sequence_[ i ] ) {
			++num_mutations;
		}
	}

	return platform::Real( num_mutations );
}

void
SequenceMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "selector", xs_string,
		"Name of the optional residue selector to use. If used, this will only use the selected residues to create the sequence." );

	xsd_type_definition_w_attributes( xsd, "Sequence",
		"Stores the amino acid sequence of the pose or subset of the pose. Distance is the number of mutations.", attlist );
}

void SequenceMetric::parse_my_tag(
	core::pose::Pose const & pose,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) {
	if ( tag->hasOption("selector") ) {
		std::string const selector_str = tag->getOption< std::string >( "selector" );

		core::select::residue_selector::ResidueSelectorCOP selector =
			datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_str )->clone();

		analyze( pose, selector );
	} else {
		analyze( pose, 0 );
	}
}


void
SequenceMetric::analyze(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSelectorCOP selector
){
	if ( sequence_.length() > 0 ) {
		return;
	} else if ( selector ) {
		std::string full_sequence = " " + pose.sequence();//adding dummy char for 1-indexing
		utility::vector1< bool > const count_this_residue = selector->apply( pose );
		sequence_ = "";

		for ( core::Size resid = 1; resid <= count_this_residue.size(); ++resid ) {
			if ( count_this_residue[ resid ] ) {
				sequence_ += full_sequence[ resid ];
			}
		}
	} else {
		sequence_ = pose.sequence();
	}
}

ClusterMetricOP
SequenceMetricCreator::create_metric( ) const {
	return utility::pointer::make_shared< SequenceMetric >();
}

void
SequenceMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	return SequenceMetric::provide_xml_schema( xsd );
}

} // metrics
} // cluster
} // multistage_...
} // protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::multistage_rosetta_scripts::cluster::metrics::SequenceMetric::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::multistage_rosetta_scripts::cluster::ClusterMetric >( this ) );
	arc( CEREAL_NVP( sequence_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::multistage_rosetta_scripts::cluster::metrics::SequenceMetric::load( Archive & arc ) {
	arc( cereal::base_class< protocols::multistage_rosetta_scripts::cluster::ClusterMetric >( this ) );
	arc( sequence_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::multistage_rosetta_scripts::cluster::metrics::SequenceMetric );
CEREAL_REGISTER_TYPE( protocols::multistage_rosetta_scripts::cluster::metrics::SequenceMetric )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_multistage_rosetta_scripts_cluster_metrics_SequenceMetric )
#endif // SERIALIZATION
