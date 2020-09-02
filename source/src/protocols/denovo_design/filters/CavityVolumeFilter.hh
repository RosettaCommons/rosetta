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
#include <basic/citation_manager/UnpublishedModuleInfo.fwd.hh>

// Project headers
#include <protocols/filters/Filter.hh>
#include <protocols/fldsgn/topology/HelixPairing.fwd.hh>
#include <protocols/fldsgn/topology/HSSTriplet.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>
#include <protocols/parser/BluePrint.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/sc/MolecularSurfaceCalculator.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.fwd.hh>

//// C++ headers
#include <string>

#include <core/io/silent/silent.fwd.hh>
#include <utility/vector1.hh>


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

	/// @brief Does this filter indicate that it is unpublished (and, by extension, that the author should be
	/// included in publications resulting from it)?  Yes, it does.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	bool filter_is_unpublished() const override;

	/// @brief Provide a list of authors and their e-mail addresses, as strings.
	/// @returns A pairs of ("Tom Linsky", "tlinsky@gmail.com").
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP > provide_authorship_info_for_unpublished() const override;

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
