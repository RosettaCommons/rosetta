// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/CycpepSymmetryFilter.hh
/// @brief A filter that examines a cyclic peptide's structure and returns TRUE if and only if it has a desired backbone symmetry.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_CycpepSymmetryFilter_hh
#define INCLUDED_protocols_cyclic_peptide_CycpepSymmetryFilter_hh

// Unit headers
#include <protocols/cyclic_peptide/CycpepSymmetryFilter.fwd.hh>
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace cyclic_peptide {

///@brief A filter that examines a cyclic peptide's structure and returns TRUE if and only if it has a desired backbone symmetry.
class CycpepSymmetryFilter : public protocols::filters::Filter {

public:

	/// @brief Constructor.
	CycpepSymmetryFilter();

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~CycpepSymmetryFilter() override;

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

	/// @brief Sets the repeats in the symmetry that we're looking for
	/// (e.g. 2 for c2 or c2/m symmetry, 3 for c3, etc.).
	void set_symm_repeats( core::Size const repeats_in );

	/// @brief Sets whether we're considering mirror (e.g. c4/m, c6/m) or non-mirror
	/// (e.g. c4, c6) symmetry.
	inline void set_mirror_symm( bool const symm_in ) { mirror_symm_ = symm_in; }

	/// @brief Sets the residue selector.
	///
	void set_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

	/// @brief Set the cutoff, in degrees, that two mainchain dihedral values must lie within in order for two residues to
	/// be considered to have the "same" value for that mainchain degree of freedom.
	void set_angle_threshold( core::Real const &setting );

	/// @brief Gets the repeats in the symmetry that we're looking for
	/// (e.g. 2 for c2 or c2/m symmetry, 3 for c3, etc.).
	inline core::Size symm_repeats() const { return symm_repeats_; }

	/// @brief Gets whether we're considering mirror (e.g. c4/m, c6/m) or non-mirror
	/// (e.g. c4, c6) symmetry.
	inline bool mirror_symm() const { return mirror_symm_; }

	/// @brief The cutoff, in degrees, that two mainchain dihedral values must lie within in order for two residues to
	/// be considered to have the "same" value for that mainchain degree of freedom.
	inline core::Real const &angle_threshold() const { return angle_threshold_; }

	/// @brief Gets the residue selector.
	///
	inline core::select::residue_selector::ResidueSelectorCOP selector() const { return selector_; }

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	clone() const override;

private:

	/****************************************
	PRIVATE FUNCTIONS
	*****************************************/

	/// @brief Given a pose and an optional ResidueSelector, get a list of residues.
	void get_residues( core::pose::Pose const &pose, core::select::residue_selector::ResidueSelectorCOP selector, utility::vector1<core::Size> &reslist_out ) const;

	/// @brief Checks that the geometry to which we're applying this filter is actually a
	/// cyclic peptide.
	bool is_cyclic_peptide( core::pose::Pose const &pose, utility::vector1< core::Size > const &residues ) const;

	/// @brief Returns TRUE if the number of mainchain torsions differs, or if the values differ; FALSE otherwise.
	/// @brief If flip is true, then torsion values are inverted prior to comparison.
	bool mainchain_torsions_differ( core::pose::Pose const &pose, core::Size const res1, core::Size const res2, bool const flip ) const;

private:

	/****************************************
	PRIVATE VARIABLES
	*****************************************/

	/// @brief The number of symmetry repeats (e.g. 2 for c2 or c2/m symmetry, 3 for c3, etc.).
	/// @details Defaults to 2.
	core::Size symm_repeats_;

	/// @brief Are we considing mirror symmetry?
	/// @details False by default; can only be true if symm_repeats_ is even.
	bool mirror_symm_;

	/// @brief A ResidueSelector to select a cyclic peptide part of a pose.
	/// @details Unused if left null.  (In this case, the filter is applied to the whole pose).
	core::select::residue_selector::ResidueSelectorCOP selector_;

	/// @brief The cutoff, in degrees, that two mainchain dihedral values must lie within in order for two residues to
	/// be considered to have the "same" value for that mainchain degree of freedom.
	/// @details Defaults to 10 degrees.
	core::Real angle_threshold_;

};

} //protocols
} //cyclic_peptide

#endif //INCLUDED_protocols_cyclic_peptide_CycpepSymmetryFilter_hh
