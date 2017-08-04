// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/residue_selectors/LigandMetalContactSelector.hh
/// @brief  This residue selector takes a selector or residue number of a ligand and returns any residues in contact with metal atoms in the ligand.
/// @author Allison Watwood (acw538@msstate.edu)

#ifndef INCLUDED_protocols_residue_selectors_LigandMetalContactSelector_HH
#define INCLUDED_protocols_residue_selectors_LigandMetalContactSelector_HH

// Unit headers
#include <protocols/residue_selectors/LigandMetalContactSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>

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
namespace residue_selectors {

/// @brief This residue selector takes a selector or residue number of a ligand and returns any residues in contact with metal atoms in the ligand.
class LigandMetalContactSelector : public core::select::residue_selector::ResidueSelector {
public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSubset ResidueSubset;

public:

	/// @brief Constructor.
	LigandMetalContactSelector();

	/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
	LigandMetalContactSelector(LigandMetalContactSelector const & src);

public:

	/// @brief Destructor.
	~LigandMetalContactSelector() override;

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.

	core::select::residue_selector::ResidueSelectorOP clone() const override;

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

	std::string get_resnum_string() const;
	void set_resnum_string(std::string const &);

	core::select::residue_selector::ResidueSelectorCOP get_input_set_selector() const;
	void set_input_set_selector(core::select::residue_selector::ResidueSelectorCOP);

	bool get_using_a_resselect() const;
	void set_using_a_resselect(bool);

	core::Real get_dist_cutoff_multiplier() const;
	void set_dist_cutoff_multiplier(core::Real const);

private:
	void calculate_ligand_resnums( core::pose::Pose const &, std::set< core::Size > &) const;

private:

	std::string resnum_string_;
	core::select::residue_selector::ResidueSelectorCOP input_set_selector_;
	bool using_a_resselect_;
	core::Real dist_cutoff_multiplier_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION


};


} //protocols
} //residue_selectors

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_residue_selectors_LigandMetalContactSelector )
#endif // SERIALIZATION

#endif //INCLUDED protocols_residue_selectors_LigandMetalContactSelector_HH
