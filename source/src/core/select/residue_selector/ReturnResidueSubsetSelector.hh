// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ReturnResidueSubsetSelector.hh
/// @brief  A simple selector that returns the set subset.  This to enable simplification of code-based interfaces to residue selectors, so that one may accept only selectors, but using this selector, we can set subsets.  This greatly reduces the interface complexity and code-complexity arising from accepting BOTH ResidueSubsets and ResidueSelectors (Which I'm terribly sick of doing at this point).
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_ReturnResidueSubsetSelector_HH
#define INCLUDED_core_select_residue_selector_ReturnResidueSubsetSelector_HH

// Unit headers
#include <core/select/residue_selector/ReturnResidueSubsetSelector.fwd.hh>

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

namespace core {
namespace select {
namespace residue_selector {

/// @brief
/// A simple selector that returns the set subset.
///  This is to enable simplification of code-based interfaces to residue selectors,
///  so that one may accept only selectors, but using this selector, we can set subsets.
///
/// This greatly reduces the c++ interface complexity and
///  private variable - complexity arising from accepting BOTH ResidueSubsets and ResidueSelectors
///  (Which I'm terribly sick of doing at this point).
///
class ReturnResidueSubsetSelector : public core::select::residue_selector::ResidueSelector {
public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSubset ResidueSubset;

public:

	/// @brief Constructor.
	ReturnResidueSubsetSelector();

	ReturnResidueSubsetSelector( ResidueSubset const & subset );

	/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
	ReturnResidueSubsetSelector(ReturnResidueSubsetSelector const & src);

	/// @brief Destructor.
	~ReturnResidueSubsetSelector() override;

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.
	ResidueSelectorOP clone() const override;

public:

	///@brief Set the ResidueSubset, which will be returned at apply-time.
	void
	set_residue_subset(ResidueSubset const & subset );

	/// @brief "Apply" function.
	/// @details Return the set subset.
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

public:

	///@brief Return an editable reference to the stored residue subset.
	///  Please be careful with this.  Use with knowledge.
	///    It is here to speed up some parts instead of creating a new subset at each point in some protocol.
	///
	ResidueSubset &
	subset();

private:

	ResidueSubset subset_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION


};


} //core
} //select
} //residue_selector


#endif //INCLUDEDcore/select/residue_selector_ReturnResidueSubsetSelector_hh
