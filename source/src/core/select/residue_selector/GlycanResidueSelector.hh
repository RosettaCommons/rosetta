// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/GlycanResidueSelector.hh
/// @brief  A ResidueSelector for carbohydrates and individual carbohydrate trees.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_GlycanResidueSelector_HH
#define INCLUDED_core_select_residue_selector_GlycanResidueSelector_HH

// Unit headers
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

namespace core {
namespace select {
namespace residue_selector {

/// @brief A ResidueSelector for carbohydrates and individual carbohydrate trees.
///  Selects all Glycan residues if no option is given or the branch going out from the root residue.
///  Selecting from root residues allows you to choose the whole glycan branch or only tips, etc.
///
class GlycanResidueSelector : public core::select::residue_selector::ResidueSelector {
public:

	/// @brief Constructor.
	GlycanResidueSelector();

	/// @brief Constructor to select tree residues from branch roots.
	///   See set_branch_residues for more
	GlycanResidueSelector( utility::vector1< bool > root_residues, bool include_root = false );

	GlycanResidueSelector( core::Size root_residue, bool include_root = false );

	/// @brief Set the residue(s) to select from.  These can be the branch points of the glycans or
	///  carbohydrate residues from which to select the downstream branch from.
	///
	///  Like the rest of a tree from a particular position.  That position could be the trunk or individual branches, which keep branching out.
	///
	///  Note that the Subset will not include the Root residue by default.
	void
	set_select_from_branch_residues( utility::vector1< bool > root_residues );

	/// @brief Set the residue to select from.  These can be the branch points of the glycans or
	///  carbohydrate residues from which to select the downstream branch from.
	///
	///  Like the rest of a tree from a particular position.  That position could be the trunk or individual branches, which keep branching out.
	///
	///  Note that the Subset will not include the Root residue by default.
	void
	set_select_from_branch_residue( core::Size root_residue );

	///@brief Option to include the root(s) we are selecting from.
	/// Default FALSE.
	void
	set_include_root(bool include_root);

	/// @brief Destructor.
	virtual ~GlycanResidueSelector();

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.
	virtual ResidueSelectorOP clone() const;

	/// @brief "Apply" function.
	/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
	/// indicating whether each residue is selected ("true") or not ("false").
	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;

	/// @brief XML parse.
	/// @details Parse RosettaScripts tags and set up this mover.
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	/// @brief Get the mover class name.
	virtual
	std::string
	get_name() const;

	/// @brief Get the mover class name.
	static std::string class_name();

	/// @brief Provide XSD information, enabling mechanical validation of input XML.
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	utility::vector1< bool > root_residues_;
	core::Size root_residue_;
	utility::vector1< std::string > parsed_positions_;
	std::string ref_pose_name_;
	bool include_root_;

};


} //core
} //select
} //residue_selector


#endif //INCLUDEDcore/select/residue_selector_GlycanResidueSelector_hh
