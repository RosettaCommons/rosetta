// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TerminusMotifScorer.hh
///
/// @brief
/// @author Frank Teets

#ifndef INCLUDED_protocols_sewing_scoring_TerminusMotifScorer_hh
#define INCLUDED_protocols_sewing_scoring_TerminusMotifScorer_hh

//Unit headers
#include <protocols/sewing/scoring/TerminusMotifScorer.fwd.hh>
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

class TerminusMotifScorer : public MotifScorer {

public:

	///@brief default construct
	TerminusMotifScorer();

	virtual ~TerminusMotifScorer()=default;

	TerminusMotifScorer( TerminusMotifScorer const & );

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

	//Additional getters and setters
	core::Size
	get_partner_residue() const;

	core::Real
	get_optimum_distance() const;

	core::Real
	get_maximum_unpenalized_variance() const;

	char
	get_terminus() const;

	void
	set_partner_residue( core::Size );

	void
	set_optimum_distance( core::Real );

	void
	set_maximum_unpenalized_variance( core::Real );

	void
	set_terminus( char );

	std::string
	get_name() const override{
		return "TerminusMotifScorer";
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
	core::Real weight_=0;
	core::Size partner_residue_=1;
	core::Real optimum_distance_=0;
	core::Real maximum_unpenalized_variance_=0;
	char terminus_='X';
};


} //scoring namespace
} //sewing namespace
} //protocols namespace

#endif
