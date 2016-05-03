// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/PhiSelector.hh
/// @brief  A ResidueSelector that selects alpha-amino acids that are either in the positive phi or negative
/// phi region of Ramachandran space (depending on user preferences).
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_core_select_residue_selector_PhiSelector_HH
#define INCLUDED_core_select_residue_selector_PhiSelector_HH

// Unit headers
#include <core/select/residue_selector/PhiSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

namespace core {
namespace select {
namespace residue_selector {

/// @brief A ResidueSelector that selects alpha-amino acids that are either in the positive phi or negative phi region of Ramachandran space (depending on user preferences).
class PhiSelector : public core::select::residue_selector::ResidueSelector {

public:

	/// @brief Constructor.
	///
	PhiSelector();

	/// @brief Destructor.
	///
	virtual ~PhiSelector();

	/// @brief Clone function.
	/// @details Copy this object and return owning pointer to the copy (created on the heap).
	virtual ResidueSelectorOP clone() const;

	/// @brief "Apply" function.
	/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
	/// indicating whether each residue is selected ("true") or not ("false").
	virtual core::select::residue_selector::ResidueSubset apply( core::pose::Pose const & pose ) const;

	/// @brief XML parse.
	/// @details Parse RosettaScripts tags and set up this mover.
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	/// @brief Get the mover class name.
	///
	virtual std::string get_name() const;

	/// @brief Get the mover class name.
	///
	static std::string class_name();

	/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
	///
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Are we selecting the residues in the positive phi region?
	/// @details Default true.
	inline bool select_positive_phi() const { return select_positive_phi_; }

	/// @brief Set whether we're selecting residues in the positive phi region.
	///
	inline void set_select_positive_phi( bool const setting ) { select_positive_phi_ = setting; }

	/// @brief Are we ignoring residues with no upper connection?
	/// @details Default true.
	inline bool ignore_unconnected_upper() const { return ignore_unconnected_upper_; }

	/// @brief Set whether we're ignoring residues with no upper connection.
	///
	inline void set_ignore_unconnected_upper( bool const setting ) { ignore_unconnected_upper_ = setting; }

private: //Private member variables:


	/// @brief Are we selecting the residues in the positive phi region?
	/// @details Default true.
	bool select_positive_phi_;

	/// @brief Are we ignoring residues with no upper connection?
	/// @details Default true.
	bool ignore_unconnected_upper_;

};


} //core
} //select
} //residue_selector


#endif //INCLUDEDcore/select/residue_selector_PhiSelector_hh
