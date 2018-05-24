// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TopNMotifScorer.hh
///
/// @brief Returns the normalized motif score for only the best N segments
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_sewing_scoring_TopNMotifScorer_hh
#define INCLUDED_protocols_sewing_scoring_TopNMotifScorer_hh

//Unit headers
#include <protocols/sewing/scoring/TopNMotifScorer.fwd.hh>
#include <protocols/sewing/scoring/MotifScorer.hh>

//Package headers
//Core headers

//Utility headers

namespace protocols {
namespace sewing  {
namespace scoring {

class TopNMotifScorer : public MotifScorer {

public:

	///@brief default construct
	TopNMotifScorer();

	virtual ~TopNMotifScorer()=default;

	TopNMotifScorer( TopNMotifScorer const & );


	static std::string type_name();

	//Virtual because derived classes override this as well
	virtual
	core::Real
	score(
		data_storage::SmartAssemblyCOP assembly
	) override;

	core::Real
	norm_motif_score(
		data_storage::SmartAssemblyCOP assembly
	);

	core::Real
	full_motif_score(
		data_storage::SmartAssemblyCOP assembly
	);


	virtual std::string
	get_name() const override;

	virtual core::Real
	get_weight() const override;

	virtual void
	set_weight( core::Real ) override;

	virtual core::Size
	get_scores_to_keep() const;

	virtual void
	set_scores_to_keep(core::Size);


	virtual void
	set_options_from_tag(
		utility::tag::TagCOP scorer_tag,
		basic::datacache::DataMap& datamap,
		protocols::filters::Filters_map const & filtermap,
		protocols::moves::Movers_map const & movermap,
		core::pose::Pose const & pose) override;

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & );

private:
	core::Real weight_=1.0;
	utility::vector1<core::Real> scores_;
	core::Size scores_to_keep_=1.0;
};


} //scoring namespace
} //sewing namespace
} //protocols namespace

#endif
