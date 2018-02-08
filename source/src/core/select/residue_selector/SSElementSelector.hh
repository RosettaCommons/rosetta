// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/SSElementSelector.hh
/// @brief  The SSElementSelector selects based on the secondary element using DSSP.
/// @author TJ Brunette (tjbrunette@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_SSElementSelector_HH
#define INCLUDED_core_select_residue_selector_SSElementSelector_HH

// Unit headers
#include <core/select/residue_selector/SSElementSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

/// @brief The SSElementSelector allows you to select objects based on their secondary structure defintion.
///  Language:
///  "A,B,C" or "A,B", "c_term","n_term","middle"
///     A = numerical position of THAT type of secondary structure. So if the secondary structure was Helix-Sheet-Helix. 2,H would refer to the second helix
///         A could be positive or negative depending if you're counting from the n_term(positive) or c_term negative
///     B = secondary structure type as assigned by DSSP , "H","E","L"
///     C = S=Start, E=End, M=Middle
///
///  you can also select with n_term or c_term.
///  Examples:
///  selection="3,H,S" to_selection="4,H,E"  Would select from the start of Helix 3 to the end of helix 4.
///     selection="-1,H,S" to_selection="4,L,E"  Would select from the start of the last helix and go to the end of the fourth loop. (order is not important)
///     selection="n_term" to_selection="3,L,M" Would select from the n_term to the middle of the third loop
class SSElementSelector : public ResidueSelector {
public:
	struct SSElement{
		Size start_res;
		Size end_res;
		std::string type;
		SSElement(Size start_res_i, Size end_res_i, Size type_res_i){
			start_res = start_res_i;
			end_res = end_res_i;
			type = type_res_i;
		}
	};

public:
	SSElementSelector();

	/// @brief Copy constructor
	///
	SSElementSelector( SSElementSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	virtual ResidueSelectorOP clone() const;

	/// @brief Destructor
	virtual ~SSElementSelector();

	/// @brief parses the secondary structure elements
	utility::vector1<SSElement> parse_ss(core::pose::Pose const & pose) const;

	// @brief This gets the SS_element when the user asks for things like the second helix
	SSElementSelector::SSElement get_SSElement(utility::vector1<SSElementSelector::SSElement> ss_elements, int goal_position, std::string type,std::string description) const;

	// @brief converts the strings shown below to protein residues
	utility::vector1<Size> convert_string_to_residues(utility::vector1<SSElementSelector::SSElement> ss_elements, std::string description) const;

	// @brief creates the subset residue which is returned. Uses the start_selection_residues and the end_selection_residues
	ResidueSubset combine_residue_selections(utility::vector1<Size> start_selection_residues, utility::vector1<Size>  end_selection_residues,core::pose::Pose const & pose) const;

	/// @brief apply.
	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;

	/// @brief tag parsing
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	);

	virtual
	std::string
	get_name() const;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private: // data members
	std::string start_;
	std::string end_;
	std::string chain_;
	core::Size reassign_short_terminal_loop_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //namespace residue_selector
} //namespace select
} //namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_SSElementSelector )
#endif // SERIALIZATION


#endif
