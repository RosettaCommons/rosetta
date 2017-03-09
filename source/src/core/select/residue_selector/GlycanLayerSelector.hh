// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/GlycanLayerSelector.hh
/// @brief  A selector for choosing glycan residues based on their layer - as measured by the residue distance to the start of the glycan tree.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_GlycanLayerSelector_HH
#define INCLUDED_core_select_residue_selector_GlycanLayerSelector_HH

// Unit headers
#include <core/select/residue_selector/GlycanLayerSelector.fwd.hh>

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
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

/// @brief A selector for choosing glycan residues based on their layer - as measured by the residue distance to the start of the glycan tree.
///  If no layer is set, will select all glycan residues.
class GlycanLayerSelector : public core::select::residue_selector::ResidueSelector {
public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSubset ResidueSubset;

public:

	/// @brief Constructor.
	GlycanLayerSelector();

	/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
	//GlycanLayerSelector(GlycanLayerSelector const & src);

public:

	/// @brief Destructor.
	~GlycanLayerSelector() override;

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.

	ResidueSelectorOP clone() const override;

	/// @brief "Apply" function.
	/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
	/// indicating whether each residue is selected ("true") or not ("false").

	ResidueSubset apply( core::pose::Pose const & pose ) const override;

public:
	
	///@brief Set the layer we will be returning.
	void
	set_layer( core::Size start, core::Size end);
	
	///@brief Set the layer as all residues greater than or equal to this number (such as the end of the tree)
	void
	set_layer_as_greater_than_or_equal_to( core::Size start );
	
	///@brief Set the layer as all residue less or equal to this number (the beginning of the tree).
	void
	set_layer_as_less_than_or_equal_to( core::Size end );
	
	
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

	core::Size start_ = 0;
	core::Size end_ = 0;
	
	core::Size start_from_as_layer_ = 0;
	core::Size end_for_layer_ = 0;
	
	bool range_set_ = false;
	
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
CEREAL_FORCE_DYNAMIC_INIT( core_select_residue_selector_GlycanLayerSelector )
#endif // SERIALIZATION

#endif //INCLUDEDcore_select_residue_selector_GlycanLayerSelector_HH
