// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file IntraDesignTerminusMotifScorer.hh
///
/// @brief
/// @author Frank Teets

#ifndef INCLUDED_protocols_sewing_scoring_IntraDesignTerminusMotifScorer_hh
#define INCLUDED_protocols_sewing_scoring_IntraDesignTerminusMotifScorer_hh

//Unit headers
#include <protocols/sewing/scoring/IntraDesignTerminusMotifScorer.fwd.hh>
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

class IntraDesignTerminusMotifScorer : public MotifScorer {

public:

	///@brief default construct
	IntraDesignTerminusMotifScorer();

	virtual ~IntraDesignTerminusMotifScorer(){}

	IntraDesignTerminusMotifScorer( IntraDesignTerminusMotifScorer const & );

	static std::string type_name();

	core::Real
	score(
		data_storage::SmartAssemblyCOP assembly
	) override;

	core::Real
	terminus_motif_score(
		data_storage::SmartAssemblyCOP assembly
	);

	virtual core::Real
	get_weight() const override;

	virtual void
	set_weight( core::Real ) override;

	//Getters and setters
	core::Real
	get_optimum_distance() const;

	core::Real
	get_maximum_unpenalized_variance() const;

	void
	set_optimum_distance( core::Real );

	void
	set_maximum_unpenalized_variance( core::Real );




	std::string
	get_name() const override{
		return "IntraDesignTerminusMotifScorer";
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
	core::Real weight_=1;
	//core::Size partner_residue_;
	core::Real optimum_distance_=0;
	core::Real maximum_unpenalized_variance_=0;
	//char terminus_;
	//private:
	//core::Real last_score_;
};


} //scoring namespace
} //sewing namespace
} //protocols namespace

#endif
