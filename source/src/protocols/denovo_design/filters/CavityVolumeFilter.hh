// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/filters/CavityVolumeFilter.hh
/// @brief Tom's Denovo Protocol. This is freely mutable and used for playing around with stuff
/// @details
/// @author Tom Linsky (tlinsky@gmail.com)
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitue.org) -- Moved from devel to protocols.


#ifndef INCLUDED_protocols_denovo_design_filters_CavityVolumeFilter_hh
#define INCLUDED_protocols_denovo_design_filters_CavityVolumeFilter_hh

// Unit headers
#include <protocols/denovo_design/filters/CavityVolumeFilter.fwd.hh>

// Basic headers

// Project headers
#include <protocols/filters/Filter.hh>


#include <core/pose/Pose.fwd.hh>

#include <core/select/residue_selector/ResidueSelector.fwd.hh>

//// C++ headers
#include <string>



namespace protocols {
namespace denovo_design {
namespace filters {

class CavityVolumeFilter : public protocols::filters::Filter {
public:

	/// @brief Initialize CavityVolumeFilter
	CavityVolumeFilter();

	/// @brief virtual constructor to allow derivation
	~CavityVolumeFilter() override;

	/// @brief Parses the CavityVolumeFilter tags
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	/// @brief Return the name of this mover.
	virtual std::string get_name() const;

	/// @brief return a fresh instance of this class in an owning pointer
	protocols::filters::FilterOP clone() const override;

	/// @brief Apply the CavityVolumeFilter. Overloaded apply function from filter base class.
	protocols::filters::FilterOP fresh_instance() const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	bool apply( core::pose::Pose const & pose ) const override;

	core::Real compute( core::pose::Pose const & pose ) const;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public: //CitationManager functions:

	/// @brief Provide the citation.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

private:   // private functions

private:   // options

private:   // other data
	/// @brief residue selector to choose residues near which to scan
	core::select::residue_selector::ResidueSelectorCOP selector_;
};


} // filters
} // denovo_design
} // protocols

#endif
