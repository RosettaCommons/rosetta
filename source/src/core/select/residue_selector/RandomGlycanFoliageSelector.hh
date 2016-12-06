// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/RandomGlycanFoliageSelector.hh
/// @brief  Selects a random carbohydrate residue from a subset or selector, then selects the rest of the glycan foliage.  Used for sampling.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_RandomGlycanFoliageSelector_HH
#define INCLUDED_core_select_residue_selector_RandomGlycanFoliageSelector_HH

// Unit headers
#include <core/select/residue_selector/RandomGlycanFoliageSelector.fwd.hh>
#include <core/select/residue_selector/GlycanResidueSelector.fwd.hh>
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

/// @brief Selects a random carbohydrate residue from a subset or selector, then selects the rest of the glycan foliage.  Used for sampling.
class RandomGlycanFoliageSelector : public core::select::residue_selector::ResidueSelector {
public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSubset ResidueSubset;

public:

	/// @brief Constructor.
	RandomGlycanFoliageSelector();

	/// @brief Constructor passing a subset from which to choose from
	RandomGlycanFoliageSelector( ResidueSubset const & subset );

	/// @brief Constructor passing a selector, from which to generate a subset on apply and from which to choose the roots from.
	RandomGlycanFoliageSelector( ResidueSelectorOP selector );

	/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
	RandomGlycanFoliageSelector(RandomGlycanFoliageSelector const & src);

public:

	/// @brief Set a subset to select the glycan root and subsequent foliage on.
	void
	set_subset(ResidueSubset const & subset);

	void
	/// @brief Set a selector to set the glycan root and subsequent foliage on.
	set_selector( ResidueSelectorCOP selector);

public:

	/// @brief Destructor.
	virtual
	~RandomGlycanFoliageSelector();

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.
	virtual
	ResidueSelectorOP clone() const;

	/// @brief "Apply" function.
	/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
	/// indicating whether each residue is selected ("true") or not ("false").
	virtual
	ResidueSubset apply( core::pose::Pose const & pose ) const;

	/// @brief XML parse.
	/// @details Parse RosettaScripts tags and set up this mover.
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	/// @brief Get the mover class name.
	virtual
	std::string
	get_name() const;

	/// @brief Get the mover class name.
	static std::string
	class_name();

	/// @brief Provide XSD information, enabling mechanical validation of input XML.
	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	/// @brief Setup anyting nessessary for this class.
	void
	setup();

private:

	ResidueSelectorCOP selector_;
	ResidueSubset subset_;

};


} //core
} //select
} //residue_selector


#endif //INCLUDEDcore/select/residue_selector_RandomGlycanFoliageSelector_hh
