// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/metrics/SequenceMetric.hh
/// @brief  Stores the amino acid sequence of the pose or subset of the pose. Distance is the number of mutations.
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_multistage_rosetta_scripts_cluster_metric_SequenceMetric_hh
#define INCLUDED_protocols_multistage_rosetta_scripts_cluster_metric_SequenceMetric_hh

// unit headers
#include <protocols/multistage_rosetta_scripts/cluster/metrics/SequenceMetric.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// package headers
#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetric.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <string>

#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace multistage_rosetta_scripts {
namespace cluster {
namespace metrics {

class SequenceMetric : public ClusterMetric {

public:
	SequenceMetric();

	SequenceMetric( std::string sequence );

	~SequenceMetric() override;

	///@brief Just counts number of mutations. Does not align.
	/// Might someday use BLOSUM62 to create different weights
	platform::Real distance( SequenceMetric const & other ) const;

	///@brief Just counts number of mutations. Does not align.
	/// Might someday use BLOSUM62 to create different weights
	platform::Real distance( ClusterMetric const & other ) const override {
		return distance( static_cast< SequenceMetric const & > ( other ) );
	}

	///@brief Store sequence from the pose
	void analyze(
		core::pose::Pose const &,
		core::select::residue_selector::ResidueSelectorCOP selector
	);

public://XML stuff
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void parse_my_tag (
		core::pose::Pose const & pose,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datacache
	) override;

public://getters/setters
	std::string const & sequence() const {
		return sequence_;
	}

	void set_sequence( std::string setting ) {
		sequence_ = setting;
	}

private:
	std::string sequence_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace metrics
} // namespace cluster
} // namespace multistage_rosetta_scripts
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_multistage_rosetta_scripts_cluster_metrics_SequenceMetric )
#endif // SERIALIZATION


#endif
