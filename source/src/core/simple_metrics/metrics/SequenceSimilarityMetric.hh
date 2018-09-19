// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/SequenceSimilarityMetric.hh
/// @brief compares the sequences of the native structure (using native flag) to the sequence of a given pose using BLOSUM62.
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_core_simple_metrics_metrics_SequenceSimilarityMetric_hh
#define INCLUDED_core_simple_metrics_metrics_SequenceSimilarityMetric_hh

#include <core/simple_metrics/metrics/SequenceSimilarityMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <core/types.hh>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace core {
namespace simple_metrics {
namespace metrics {

/// @brief compares the sequences of the native structure (using native flag) to the sequence of a given pose using BLOSUM 62.
class SequenceSimilarityMetric : public core::simple_metrics::RealMetric{

public:

	SequenceSimilarityMetric();

	SequenceSimilarityMetric( select::residue_selector::ResidueSelectorCOP selector );

	SequenceSimilarityMetric(SequenceSimilarityMetric const & src);

	~SequenceSimilarityMetric() override;

	core::Real
	calculate( pose::Pose const & pose ) const override;

	static core::Real
	score( std::string const & seq1, std::string const & seq2, bool normalize );

public:

	///@brief the selector will decide which positions will be scored
	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector ){
		selector_ = selector;
	}

	void
	set_residue_selector(
		core::select::residue_selector::ResidueSelectorCOP selector,
		bool apply_selector_to_native
	){
		selector_ = selector;
		set_apply_selector_to_native( apply_selector_to_native );
	}

	void
	set_native_pose( core::pose::PoseCOP setting ){
		native_pose_ = setting;
	}

	void
	set_apply_selector_to_native( bool setting ){
		apply_selector_to_native_ = setting;
	}

	void
	set_native_pose( core::pose::PoseCOP setting, bool apply_selector_to_native ){
		set_native_pose( setting );
		set_apply_selector_to_native( apply_selector_to_native );
	}

	void
	set_normalize( bool setting ){
		normalize_ = setting;
	}

public:

	///@brief Name of the class
	std::string
	name() const override;

	///@brief Name of the class for creator.
	static
	std::string
	name_static();

	///@brief Name of the metric
	std::string
	metric() const override;

public:
	/// @brief called by parse_my_tag -- should not be used directly
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	SimpleMetricOP
	clone() const override;

private:
	core::select::residue_selector::ResidueSelectorCOP selector_;
	mutable core::pose::PoseCOP native_pose_;//mutable so that we can load the native flag at the last minute if this has not yet been assigned
	bool apply_selector_to_native_;
	bool normalize_;
};

} //core
} //simple_metrics
} //metrics



#endif //INCLUDED_core_metrics_SequenceSimilarityMetric_hh





