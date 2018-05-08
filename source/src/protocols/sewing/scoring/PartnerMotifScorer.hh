// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file PartnerMotifScorer.hh
///
/// @brief
/// @author Frank Teets

#ifndef INCLUDED_protocols_sewing_scoring_PartnerMotifScorer_hh
#define INCLUDED_protocols_sewing_scoring_PartnerMotifScorer_hh

//Unit headers
#include <protocols/sewing/scoring/PartnerMotifScorer.fwd.hh>
#include <protocols/sewing/scoring/MotifScorer.hh>

//Package headers
#include <protocols/sewing/data_storage/SmartAssembly.hh>

//Core headers
#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/chemical/ResidueTypeSet.hh>

//Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace sewing  {
namespace scoring {

class PartnerMotifScorer : public MotifScorer {

public:

	///@brief default construct
	PartnerMotifScorer();

	virtual ~PartnerMotifScorer()=default;

	PartnerMotifScorer( PartnerMotifScorer const & );



	static std::string type_name();

	core::Real
	score(
		data_storage::SmartAssemblyCOP assembly
	) override;

	core::Real
	interface_motif_score(
		data_storage::SmartAssemblyCOP assembly
	);

	virtual core::Real
	get_weight() const override;

	virtual void
	set_weight( core::Real ) override;

	std::string
	get_name() const override{
		return "PartnerMotifScorer";
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
	//private:
	//core::Real last_score_;
};


} //scoring namespace
} //sewing namespace
} //protocols namespace

#endif
