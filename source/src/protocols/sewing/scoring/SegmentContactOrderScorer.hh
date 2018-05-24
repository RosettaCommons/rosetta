// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SegmentContactOrderScorer.hh
///
/// @brief Favors assemblies whose segments form contacts with segments distant in the assembly
/// @author Sharon Guffy (guffy@email.unc.edu)

#ifndef INCLUDED_protocols_sewing_scoring_SegmentContactOrderScorer_hh
#define INCLUDED_protocols_sewing_scoring_SegmentContactOrderScorer_hh

//Unit headers
#include <protocols/sewing/scoring/SegmentContactOrderScorer.fwd.hh>
#include <protocols/sewing/scoring/MotifScorer.hh>

//Package headers
#include <protocols/sewing/data_storage/SmartAssembly.hh>

//Core headers
#include <core/types.hh>
//Utility headers

namespace protocols {
namespace sewing  {
namespace scoring {

class SegmentContactOrderScorer : public MotifScorer {

public:

	///@brief default construct
	SegmentContactOrderScorer();

	virtual ~SegmentContactOrderScorer()=default;

	SegmentContactOrderScorer( SegmentContactOrderScorer const & );

	static std::string type_name();

	core::Real
	score(
		data_storage::SmartAssemblyCOP assembly
	) override;

	core::Real
	contact_order_score(
		data_storage::SmartAssemblyCOP assembly
	);

	core::Real
	get_weight() const override;

	void
	set_weight( core::Real ) override;

	std::string
	get_name() const override{
		return "SegmentContactOrderScorer";
	}

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & );

	void
	set_options_from_tag(
		utility::tag::TagCOP scorer_tag,
		basic::datacache::DataMap& datamap,
		protocols::filters::Filters_map const & filtermap,
		protocols::moves::Movers_map const & movermap,
		core::pose::Pose const & pose) override;

private:
	core::Real weight_=1.0;
};


} //scoring namespace
} //sewing namespace
} //protocols namespace

#endif
