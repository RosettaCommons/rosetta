// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/cyclic_peptide/RamaMutationSelector.hh
/// @brief  Selects positions that would have a rama_prepro score below a given threshold IF mutated to a given residue type.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_RamaMutationSelector_HH
#define INCLUDED_protocols_cyclic_peptide_RamaMutationSelector_HH

// Unit headers
#include <protocols/cyclic_peptide/RamaMutationSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace cyclic_peptide {

/// @brief Selects positions that would have a rama_prepro score below a given threshold IF mutated to a given residue type.
class RamaMutationSelector : public core::select::residue_selector::ResidueSelector {
public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSubset ResidueSubset;

public:

	/// @brief Constructor.
	RamaMutationSelector();

	/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
	RamaMutationSelector(RamaMutationSelector const & src);

public:

	/// @brief Destructor.
	~RamaMutationSelector() override;

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.

	ResidueSelectorOP clone() const override;

	/// @brief "Apply" function.
	/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
	/// indicating whether each residue is selected ("true") or not ("false").
	ResidueSubset apply( core::pose::Pose const & pose ) const override;

	/// @brief XML parse.
	/// @details Parse RosettaScripts tags and set up this mover.
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) override;

	/// @brief Get the mover class name.
	std::string
	get_name() const override;

	/// @brief Get the mover class name.
	static std::string
	class_name();

	/// @brief Provide XSD information, enabling mechanical validation of input XML.
	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public: //Setters

	/// @brief Set the residue type (full name) to which we are considering
	/// mutations.  If set to an empty string, selection is based on rama_prepro scoring
	/// of the residue type at each position.
	void set_target_type( std::string const &type_in );

	/// @brief Set the score threshold.
	/// @details Positions which, when mutated to the target type, have a rama_prepro
	/// score above this threshold are not selected.
	void set_score_threshold( core::Real const &threshold_in );

	/// @brief Set the weight multiplier of the rama_prepro term.
	/// @details Defaults to 0.45 to match beta_nov15.
	void set_rama_prepro_multiplier( core::Real const &multiplier_in );

public: //Getters

	/// @brief Get the residue type (full name) to which we are considering
	/// mutations.  If set to an empty string, selection is based on rama_prepro scoring
	/// of the residue type at each position.
	inline std::string const & target_type() const { return target_type_; }

	/// @brief Get the score threshold.
	/// @details Positions which, when mutated to the target type, have a rama_prepro
	/// score above this threshold are not selected.
	inline core::Real const &score_threshold() const { return score_threshold_; }

	/// @brief Get the weight multiplier of the rama_prepro term.
	/// @details Defaults to 0.45 to match beta_nov15.
	inline core::Real const &rama_prepro_multiplier() const { return rama_prepro_multiplier_; }

private: // Private methods

	/// @brief Is a position in a pose a terminal position?
	/// @details Returns true if (a) it has a terminal type, (b) its lower_connect is missing,
	/// (c) its lower_connect is unconnected, (d) its upper_connect is missing, or (e) its
	/// upper_connect is unconnected.
	bool is_terminus( core::conformation::Residue const &rsd ) const;

private: // Private data

	/// @brief The residue type (full name) to which we are considering
	/// mutations.  If left blank, selection is based on rama_prepro scoring
	/// of the residue type at each position.
	std::string target_type_;

	/// @brief The rama_prepro score threshold above which residues are not
	/// selected.
	core::Real score_threshold_;

	/// @brief The weight in the scorefunction of the rama_prepro term.
	/// @details Defaults to 0.45 to match beta_nov15.
	core::Real rama_prepro_multiplier_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION


};


} //protocols
} //cyclic_peptide

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_cyclic_peptide_RamaMutationSelector )
#endif // SERIALIZATION

#endif //INCLUDEDprotocols_cyclic_peptide_RamaMutationSelector_HH
