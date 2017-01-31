// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/SimpleHbondsToAtomFilter.hh
/// @brief Simple filter for detercting Hbonds to atom with energy < energy cutoff
/// @author Benjamin Basanta (basantab@uw.edu)

#ifndef INCLUDED_protocols_simple_filters_SimpleHbondsToAtomFilter_hh
#define INCLUDED_protocols_simple_filters_SimpleHbondsToAtomFilter_hh

// Unit headers
#include <protocols/simple_filters/SimpleHbondsToAtomFilter.fwd.hh>
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/Tag.fwd.hh> //transcluded from Filter.hh
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Filter.hh

namespace protocols {
namespace simple_filters {

///@brief Simple filter for detercting Hbonds to atom with energy < energy cutoff
class SimpleHbondsToAtomFilter : public protocols::filters::Filter {

public:
	SimpleHbondsToAtomFilter();

	// destructor (important for properly forward-declaring smart-pointer members)
	~SimpleHbondsToAtomFilter() override;

	/// @brief returns true if the structure passes the filter, false otherwise
	bool
	apply( core::pose::Pose const & pose ) const override;

	/// @brief required for reporting score values
	core::Real
	report_sm( core::pose::Pose const & pose ) const override;

	/// @brief allows printing data to a stream
	void
	report( std::ostream & os, core::pose::Pose const & pose ) const override;

public:
	core::Size get_target_residue() const;
	void set_target_residue( core::Size target_residue );
	core::Size get_n_partners() const;
	void set_n_partners( core::Size n_partners );
	std::string get_target_atom_name() const;
	void set_target_atom_name( std::string atom_name );
	core::scoring::ScoreFunctionOP get_scorefxn() const;
	void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	core::Real get_hb_e_cutoff() const;
	void set_hb_e_cutoff( core::Real hb_e_cutoff );
	core::Size compute( core::pose::Pose const & pose ) const;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief parse XML tag (to use this Filter in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	clone() const override;

private:
	std::string target_atom_name_;
	core::Size target_residue_;
	core::Size n_partners_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Real hb_e_cutoff_;

};

} //protocols
} //simple_filters

#endif //INCLUDED_protocols_simple_filters_SimpleHbondsToAtomFilter_hh
