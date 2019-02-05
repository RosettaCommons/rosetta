// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResiduePropertySelector.hh
/// @brief  A residue selector that selects based on set residue properties.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_ResiduePropertySelector_HH
#define INCLUDED_core_select_residue_selector_ResiduePropertySelector_HH

// Unit headers
#include <core/select/residue_selector/ResiduePropertySelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/chemical/ResidueProperty.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

enum basic_selection_logic{
	and_logic = 1,
	or_logic,
};

static const std::map< std::string, basic_selection_logic >
logic_map( { { "or_logic", or_logic }, { "and_logic", and_logic } } );

/// @brief A residue selector that selects based on set residue properties.
///
/// @details
///  Default is to use AND logic for multiple properties.  This can be changed via set_selection_logic.
///
class ResiduePropertySelector : public core::select::residue_selector::ResidueSelector {
public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSubset ResidueSubset;

public:

	/// @brief Defualt Constructor.
	ResiduePropertySelector();

	/// @brief Constructor setting a single property for selection
	ResiduePropertySelector( core::chemical::ResidueProperty property );

	/// @brief Constructor setting a list of properties for selection
	ResiduePropertySelector( utility::vector1< chemical::ResidueProperty > properties );

	/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
	//ResiduePropertySelector(ResiduePropertySelector const & src);

public:

	/// @brief Destructor.
	~ResiduePropertySelector() override;

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.

	ResidueSelectorOP clone() const override;

	/// @brief "Apply" function.
	/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
	/// indicating whether each residue is selected ("true") or not ("false").

	ResidueSubset apply( core::pose::Pose const & pose ) const override;

public:

	///@brief Set a single property for selection
	void
	set_property( core::chemical::ResidueProperty property );

	///@brief Add a property to the list for selection
	void
	add_property( core::chemical::ResidueProperty property );

	///@brief Set a list of properties for selection
	void
	set_properties( utility::vector1< core::chemical::ResidueProperty > properties  );

	///@brief Set the logic for multiple sets of properties
	/// Default is AND logic.
	///
	void
	set_selection_logic( basic_selection_logic logic);

public:

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


private:

	utility::vector1< core::chemical::ResidueProperty > properties_;
	basic_selection_logic logic_ = and_logic;



#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION


};


} //core
} //select
} //residue_selector

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_select_residue_selector_ResiduePropertySelector )
#endif // SERIALIZATION

#endif //INCLUDEDcore_select_residue_selector_ResiduePropertySelector_HH
