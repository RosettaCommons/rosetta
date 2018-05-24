// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LigandScorer.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_sewing_scoring_LigandScorer_hh
#define INCLUDED_protocols_sewing_scoring_LigandScorer_hh

//Unit headers
#include <protocols/sewing/scoring/LigandScorer.fwd.hh>
#include <protocols/sewing/scoring/LigandAssemblyScorer.hh>

//Package headers
#include <protocols/sewing/data_storage/SmartAssembly.fwd.hh>

//Core headers
#include <core/scoring/func/Func.fwd.hh>
#include <core/chemical/AtomTypeSet.hh>

//Utility headers
#include <utility/vector1.hh>
#include <numeric/constants.hh>

namespace protocols {
namespace sewing  {
namespace scoring {

class LigandScorer : public LigandAssemblyScorer {

public:

	///@brief default construct
	LigandScorer();

	virtual ~LigandScorer()=default;
	LigandScorer( LigandScorer const & );


	core::Real
	score(
		data_storage::SmartAssemblyCOP assembly
	) override;

	//core::Real
	//norm_motif_score(
	// AssemblyCOP assembly
	//);
	std::string
	get_name() const override;

	core::Real
	get_weight() const override;

	void
	set_weight( core::Real ) override;

	core::Real
	get_angle_multiplier( core::Real theta_positive );

	core::Real
	get_last_score() const override;

	void
	set_last_score( core::Real ) override;

	core::Real
	get_old_last_score() const override;

	void
	set_old_last_score( core::Real ) override;


	//Additional getters and setters
	core::scoring::func::FuncCOP
	get_func() const;

	core::Real
	get_cutoff_distance() const;

	core::Real
	get_cutoff_angle() const;

	void
	set_func( core::scoring::func::FuncOP );

	void
	set_cutoff_distance( core::Real );

	void
	set_cutoff_angle( core::Real );

	static std::string type_name();

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
	core::Real last_score_;
	core::Real old_last_score_;

	core::Real weight_=1.0;
	core::scoring::func::FuncOP func_;
	core::Real cutoff_distance_=5.0;
	core::Real cutoff_angle_=numeric::constants::r::pi_over_2/2.0;

	core::chemical::AtomTypeSetCAP atom_types_;
	//private:
	//core::Real last_score_;
};


} //scoring namespace
} //sewing namespace
} //protocols namespace

#endif
