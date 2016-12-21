// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/OversaturatedHbondAcceptorFilter.hh
/// @brief Headers for the OversaturatedHbondAcceptorFilter.  This filter flags poses containing more than two
/// hydrogen bonds to an oxygen atom, a common pathology that results from Rosetta's pairwise-decomposible
/// scorefunction, which can't penalize excessive hydrogen bonds.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_OversaturatedHbondAcceptorFilter_hh
#define INCLUDED_protocols_cyclic_peptide_OversaturatedHbondAcceptorFilter_hh

// Unit headers
#include <protocols/cyclic_peptide/OversaturatedHbondAcceptorFilter.fwd.hh>
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace cyclic_peptide {

/// @brief This filter flags poses containing more than two hydrogen bonds to an oxygen atom, a common pathology
/// that results from Rosetta's pairwise-decomposible scorefunction, which can't penalize excessive hydrogen bonds.
class OversaturatedHbondAcceptorFilter : public protocols::filters::Filter {

public:

	/// @brief Default constructor.
	///
	OversaturatedHbondAcceptorFilter();

	/// @brief Copy constructor.
	///
	OversaturatedHbondAcceptorFilter( OversaturatedHbondAcceptorFilter const &src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	///
	~OversaturatedHbondAcceptorFilter() override;

	/// @brief Required in the context of the parser/scripting scheme.
	/// @details Make a copy of this object and return an owning pointer to the copy.
	protocols::filters::FilterOP
	clone() const override;

	/// @brief Required in the context of the parser/scripting scheme.
	///
	protocols::filters::FilterOP
	fresh_instance() const override;

public:

	// ---------- PUBLIC FUNCTIONS ------------------

	virtual std::string
	get_name() const;

	/// @brief Returns true if the structure passes the filter, false otherwise.
	///
	bool
	apply( core::pose::Pose const & pose ) const override;

	/// @brief Required for reporting score values.
	///
	core::Real
	report_sm( core::pose::Pose const & pose ) const override;

	/// @brief Allows printing a summary of this filter to a stream.
	///
	void
	report( std::ostream & os, core::pose::Pose const & pose ) const override;

	/// @brief Parse XML tag (to use this Mover in RosettaScripts).
	///
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

public:

	// ---------- SETTERS AND GETTERS ---------------

	/// @brief Set the maximum allowed number of instances of an oversaturated
	/// hydrogen bond acceptor.
	void set_max_allowed_oversaturated( core::Size const max_allowed );

	/// @brief Get the maximum allowed number of instances of an oversaturated
	/// hydrogen bond acceptor.
	core::Size max_allowed_oversaturated() const;

	/// @brief Set the threshold for considering something to be a
	/// hydrogen bond.
	void set_hbond_energy_cutoff( core::Real const &input_value );

	/// @brief Get the threshold for considering something to be a
	/// hydrogen bond.
	core::Real const & hbond_energy_cutoff() const;

	/// @brief Set whether we only consider mainchain hydrogen bond
	/// donors and acceptors.
	void set_consider_mainchain_only ( bool const input_setting );

	/// @brief Get whether we only consider mainchain hydrogen bond
	/// donors and acceptors.
	bool consider_mainchain_only() const;

	/// @brief Set the residue selector for donor residues.
	/// @details Clones the input.
	void set_donor_selector( core::select::residue_selector::ResidueSelectorCOP donor_selector_in );

	/// @brief Get the residue selector for donor residues.
	///
	core::select::residue_selector::ResidueSelectorCOP donor_selector() const;

	/// @brief Set the residue selector for acceptor residues.
	/// @details Clones the input.
	void set_acceptor_selector( core::select::residue_selector::ResidueSelectorCOP acceptor_selector_in );

	/// @brief Get the residue selector for acceptor residues.
	///
	core::select::residue_selector::ResidueSelectorCOP acceptor_selector() const;

	/// @brief Set the scorefunction.
	/// @details Clones the input.
	void set_scorefxn( core::scoring::ScoreFunctionCOP sfxn_in);

	/// @brief Get the scorefunction.
	///
	core::scoring::ScoreFunctionCOP scorefxn() const;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	// ---------- PRIVATE FUNCTIONS -----------------

	/// @brief Given an atom, return the maximum number of hydrogen bonds that an atom of this type
	/// is allowed to receive.
	/// @details For now, this function returns 0 if this is not an acceptor, 2 if it is.  This is crude,
	/// of course.  Ultimately, it should figure out how many hbonds an atom can actually receive.
	core::Size max_allowed_hbonds( core::pose::Pose const &pose, core::Size const res_index, core::Size const atom_index) const;

	/// @brief The function that actually calculates the value that this filter returns, called by the apply(),
	/// report(), and report_sm() functions.
	/// @details Returns the number of atoms receiving more than the allowed number of hydrogen bonds.
	core::Size compute( core::pose::Pose const &pose ) const;

	// ---------- PRIVATE MEMBER VARIABLES ----------

	/// @brief What is the maximum allowed number of instances of an oversaturated
	/// hydrogen bond acceptor?
	/// @details Defaults to 0.
	core::Size max_allowed_oversaturated_;

	/// @brief The energy cutoff for considering something to be a hydrogen bond.
	/// @details Defaults to -0.1.
	core::Real hbond_energy_cutoff_;

	/// @brief Should we only consider mainchain hydrogen bond donors and acceptors?
	/// @details Defaults to true.
	bool consider_mainchain_only_;

	/// @brief Optional residue selector for the residues that could possibly be donating
	/// hydrogen bonds.
	core::select::residue_selector::ResidueSelectorCOP donor_selector_;

	/// @brief Optional residue selector for the residues that could possibly be accepting
	/// hydrogen bonds.
	core::select::residue_selector::ResidueSelectorCOP acceptor_selector_;

	/// @brief The scorefunction to use for hydrogen bond scoring.
	/// @details If no scorefunction is provided, then the default scorefunction is used.
	core::scoring::ScoreFunctionOP scorefxn_;


};

} //protocols
} //cyclic_peptide

#endif //INCLUDED_protocols_cyclic_peptide_OversaturatedHbondAcceptorFilter_hh
