// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file InterModelMotifScorer.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_sewing_scoring_InterModelMotifScorer_hh
#define INCLUDED_protocols_sewing_scoring_InterModelMotifScorer_hh

//Unit headers
#include <protocols/sewing/scoring/InterModelMotifScorer.fwd.hh>
#include <protocols/sewing/scoring/MotifScorer.hh>

//Package headers
#include <protocols/sewing/data_storage/SmartAssembly.fwd.hh>

//Core headers
#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/chemical/ResidueTypeSet.hh>

//Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace sewing  {
namespace scoring {

class InterModelMotifScorer : public MotifScorer {

public:

	///@brief default construct
	InterModelMotifScorer();

	virtual ~InterModelMotifScorer()=default;

	InterModelMotifScorer( InterModelMotifScorer const & );

	core::Real
	score(
		data_storage::SmartAssemblyCOP assembly
	) override;

	core::Real
	full_motif_score(
		data_storage::SmartAssemblyCOP assembly
	);

	//core::Real
	//norm_motif_score(
	// AssemblyCOP assembly
	//);
	std::string
	get_name() const override;

	static std::string type_name();

	virtual core::Real
	get_weight() const override;

	virtual void
	set_weight( core::Real ) override;


	void
	set_options_from_tag(
		utility::tag::TagCOP scorer_tag,
		basic::datacache::DataMap& datamap,
		protocols::filters::Filters_map const & filtermap,
		protocols::moves::Movers_map const & movermap,
		core::pose::Pose const & pose) override;

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & );


private:
	core::Real weight_=10.0;
	//private:
	//core::Real last_score_;
};


} //scoring namespace
} //sewing namespace
} //protocols namespace

#endif
