// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SubsetPartnerMotifScorer.hh
///
/// @brief
/// @author Frank Teets

#ifndef INCLUDED_protocols_sewing_scoring_SubsetPartnerMotifScorer_hh
#define INCLUDED_protocols_sewing_scoring_SubsetPartnerMotifScorer_hh

//Unit headers
#include <protocols/sewing/scoring/SubsetPartnerMotifScorer.fwd.hh>
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

class SubsetPartnerMotifScorer : public MotifScorer {

public:

	///@brief default construct
	SubsetPartnerMotifScorer();

	virtual ~SubsetPartnerMotifScorer()=default;

	SubsetPartnerMotifScorer( SubsetPartnerMotifScorer const & );



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
	get_starts_and_ends() const;

	void
	set_starts_and_ends(std::string);



	core::Size
	get_region_start() const;

	core::Size
	get_region_end() const;

	void
	set_region_start( core::Size );

	void
	set_region_end( core::Size );


	std::string
	get_name() const override{
		return "SubsetPartnerMotifScorer";
	};

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
	core::Real weight_=1;
	std::string starts_and_ends_=""; //I don't think this variable actually does anything anymore
	core::Size region_start_=1;
	core::Size region_end_=1;
};


} //scoring namespace
} //sewing namespace
} //protocols namespace

#endif
