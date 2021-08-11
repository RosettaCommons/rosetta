// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ConfluentSelector.hh
/// @brief  selector that selects all the residues between those selected by another selector that are not interspersed with those selected by a third selector
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_ConfluentSelector_HH
#define INCLUDED_core_select_residue_selector_ConfluentSelector_HH

// Unit headers
#include <core/select/residue_selector/ConfluentSelector.fwd.hh>

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

namespace core {
namespace select {
namespace residue_selector {

/// @brief selector that selects all the residues between those selected by another selector that are not interspersed with those selected by a third selector
class ConfluentSelector : public core::select::residue_selector::ResidueSelector {
public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSubset ResidueSubset;

public:

	/// @brief Constructor.
	ConfluentSelector();

	/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
	//ConfluentSelector(ConfluentSelector const & src);

public:

	/// @brief Destructor.
	~ConfluentSelector() override;

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

	void set_terminus_selector(core::select::residue_selector::ResidueSelectorCOP terminus_selector);

	void set_breaking_selector(core::select::residue_selector::ResidueSelectorCOP breaking_selector);

	/// @brief Get the mover class name.
	std::string
	get_name() const override;

	/// @brief Get the mover class name.
	static std::string
	class_name();

	/// @brief Provide XSD information, enabling mechanical validation of input XML.
	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::select::residue_selector::ResidueSelectorCOP terminus_selector_;
	core::select::residue_selector::ResidueSelectorCOP breaking_selector_;



};


} //residue_selectors
} //pose_sewing
} //protocols


#endif //INCLUDEDcore_select_residue_selector_ConfluentSelector_HH
