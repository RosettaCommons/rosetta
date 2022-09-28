// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/drug_design/ApplyChemistryMover.hh
/// @brief Apply a given Chemistry modifier to the ResidueType at a given position, then replace the ResidueType at that position.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_ApplyChemistryMover_HH
#define INCLUDED_protocols_drug_design_ApplyChemistryMover_HH

// Unit headers
#include <protocols/drug_design/ApplyChemistryMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/chemistries/Chemistry.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace drug_design {

///@brief Apply a given Chemistry modifier to the ResidueType at a given position, then replace the ResidueType at that position.
class ApplyChemistryMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	ApplyChemistryMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	ApplyChemistryMover( ApplyChemistryMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~ApplyChemistryMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:

	core::select::residue_selector::ResidueSelectorCOP
	residue_selector() const { return residue_selector_; }

	std::string const & residues() const { return residues_; }

	core::Size residue_id() const { return residue_id_; }

	/// @details returns by value and not by reference so you get a COP and can't modify the chemistries
	utility::vector1< protocols::chemistries::ChemistryCOP >
	chemistries() const { return chemistries_; }

	std::string const & new_name() const { return new_name_; }

	std::string const & tag() const { return tag_; }

	/// @brief Set the residues to act on -- Residue selector
	/// The set value will take precedence over previously set residue specifications
	void residue_selector( core::select::residue_selector::ResidueSelectorCOP setting ) {
		residue_selector_ = setting;
	}

	/// @brief Set the residues to act on -- comma-separated list of (string) residues
	/// The set value will take precedence over previously set residue specifications
	/// (and will clear the residue selector)
	void residues( std::string const & setting ) {
		residue_selector_ = nullptr; // If we explicitly set residues_, disable the others
		residues_ = setting;
	}

	/// @brief Set the residues to act on -- Pose number
	/// The set value will take precedence over previously set residue specifications
	/// (and will clear the residue selector & residue string settings)
	void residue_id( core::Size setting ) {
		residue_selector_ = nullptr; // If we explicitly set residue_id, disable the others
		residues_ = "";
		residue_id_ = setting;
	}

	void add_chemistry( protocols::chemistries::ChemistryOP setting );

	void new_name( std::string const & setting ) { new_name_ = setting; }

	void tag( std::string const & setting ) { tag_ = setting; }

	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	//ApplyChemistryMover & operator=( ApplyChemistryMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private: // methods

private: // data

	/// @brief Which Residue selector to use.
	/// If set, will be used in preference to residue_ and residue_id_
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;

	/// @brief Which residue in the pose to substitute.
	/// If not empty residues_ will be used in preference to residue_id_
	std::string residues_;
	core::Real residue_id_ = 0;

	/// @brief Which chemistries to use?
	utility::vector1< protocols::chemistries::ChemistryOP > chemistries_;

	/// @brief If set, this is the new name for the Residue type.
	/// (Only if there's only one type.)
	std::string new_name_;

	/// @brief What tag to use on the residue type to distinguish it from other residues
	/// May be added multiple times, if it results in a name conflict.
	std::string tag_ = "mod";
};

std::ostream &
operator<<( std::ostream & os, ApplyChemistryMover const & mover );

} //protocols
} //drug_design

#endif //protocols_drug_design_ApplyChemistryMover_HH
