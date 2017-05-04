// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/BuriedSurfaceAreaFilter.hh
/// @brief Headers for the BuriedSurfaceAreaFilter.
/// @details Calculates buried surface area (exposed surface area minus total surface area, on a per-residue basis).  Accepts
/// a residue selector to allow buried subsets to be considered.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_BuriedSurfaceAreaFilter_hh
#define INCLUDED_protocols_simple_filters_BuriedSurfaceAreaFilter_hh

// Unit headers
#include <protocols/simple_filters/BuriedSurfaceAreaFilter.fwd.hh>
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/Tag.fwd.hh> //transcluded from Filter.hh
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Filter.hh

namespace protocols {
namespace simple_filters {

enum BSAF_mode {
	BSAF_all_atoms = 1, //keep first
	BSAF_hydrophobic_atoms,
	BSAF_polar_atoms, //keep second-to-last
	BSAF_end_of_list = BSAF_polar_atoms //keep last
};

///@brief Calculates buried surface area (exposed surface area minus total surface area, on a per-residue basis).  Accepts a residue selector to allow buried subsets to be considered.
class BuriedSurfaceAreaFilter : public protocols::filters::Filter {

public:
	BuriedSurfaceAreaFilter();

	// destructor (important for properly forward-declaring smart-pointer members)
	~BuriedSurfaceAreaFilter() override;

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

public: //Setters:

	/// @brief Sets the residue selector to use to select a subset of residues for which to calculate
	/// buried surface area.
	/// @details Copies the input owning pointer; does not clone.  This means that residue selectors could be
	/// shared with other Rosetta modules.
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

	/// @brief Set whether only hydrophobic residues and alanine (FAMILYVW) are considered.  False by default.
	/// @details The selection FAMILYVW is combined with the residue selector, if specified, using AND logic.
	void set_select_only_FAMILYVW( bool const setting );

	/// @brief Set the cutoff buried surface area below which (or above which, if filter_out_low_ is false) structures
	/// are discarded.
	void set_cutoff_buried_surface_area( core::Real const &setting );

	/// @brief Set whether structures with less than the cutoff buried area or greater than the cutoff
	/// buried area are discarded.
	/// @details If true, structures with less than the cutoff buried area are discrded.
	void set_filter_out_low( bool const setting );

	/// @brief Set the atom mode (the subset of atoms to use for the calculation) by string.
	///
	void set_atom_mode( std::string const &setting );

	/// @brief Set the atom mode (the subset of atoms to use for the calculation).
	///
	void set_atom_mode( BSAF_mode const setting );

public: //Getters:

	/// @brief Gets the residue selector to use to select a subset of residues for which to calculate
	/// buried surface area.
	inline core::select::residue_selector::ResidueSelectorCOP residue_selector() const { return residue_selector_; }

	/// @brief Get whether only hydrophobic residues and alanine (FAMILYVW) are considered.  False by default.
	/// @details The selection FAMILYVW is combined with the residue selector, if specified, using AND logic.
	inline bool select_only_FAMILYVW() const { return select_only_FAMILYVW_; }

	/// @brief Get the cutoff buried surface area below which (or above which, if filter_out_low_ is false) structures
	/// are discarded.
	inline core::Real const & cutoff_buried_surface_area() const { return cutoff_area_; }

	/// @brief Get whether structures with less than the cutoff buried area or greater than the cutoff
	/// buried area are discarded.
	/// @details If true, structures with less than the cutoff buried area are discrded.
	inline bool filter_out_low() const { return filter_out_low_; }

	/// @brief Get the atom mode (the subset of atoms to use).
	///
	inline BSAF_mode atom_mode() const { return atom_mode_; }

private: //Methods:

	/// @brief Common function called by apply(), report(), and report_sm().
	/// @details Does the actual computation.  Returns true if the filter passes and false if it fails.
	bool compute( std::ostream &os, core::pose::Pose const &pose, core::Real &buried_surf_area ) const;

	/// @brief Is this a residue type for which we can calculate a total SASA value?
	///
	bool is_allowed_type( char const oneletter_code ) const;

	/// @brief Is this a hydrophobic residue type (FAMILYVW)?
	///
	bool is_FAMILYVW( char const oneletter_code ) const;

	/// @brief Given a pose (input) and a ResidueSubset (output), compute the residues that should be
	/// operated on.
	void compute_residue_selection( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset &selection ) const;

private: //Data:

	/// @brief A residue selector to use to select a subset of residues for which to calculate
	/// buried surface area.  Unused if not specified (i.e. all residues are used by default).
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;

	/// @brief If true, only hydrophobic residues and alanine (FAMILYVW) are considered.  False by default.
	/// @details The selection FAMILYVW is combined with the residue selector, if specified, using AND logic.
	bool select_only_FAMILYVW_;

	/// @brief The subset of atoms to consider.  Default is all atoms, but only hydrophobic atoms or only polar atoms
	/// are both possible, too.
	BSAF_mode atom_mode_;

	/// @brief The cutoff buried surface area below which (or above which, if filter_out_low_ is false) structures
	/// are discarded.
	/// @details Defaults to 500 square Angstroms (arbitrarily chosen).
	core::Real cutoff_area_;

	/// @brief If true, structures with less than the cutoff buried area are discrded.
	/// @details True by default.
	bool filter_out_low_;

};

} //protocols
} //simple_filters

#endif //INCLUDED_protocols_simple_filters_BuriedSurfaceAreaFilter_hh
