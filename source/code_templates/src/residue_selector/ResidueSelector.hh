// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   --path--/--class--.hh
/// @brief  --brief--
/// @author --name-- (--email--)

#ifndef INCLUDED_--path_underscore--_--class--_HH
#define INCLUDED_--path_underscore--_--class--_HH

// Unit headers
#include <--path--/--class--.fwd.hh>

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

--namespace--

/// @brief --brief--
class --class-- : public core::select::residue_selector::ResidueSelector {
public:
	// derived from base class

	/// @brief Constructor.
	///
	--class--();

	/// @brief Destructor.
	///
	virtual ~--class--();

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.
	virtual core::select::residue_selector::ResidueSelectorOP --class--::clone() const;

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
	virtual
	std::string
	get_name() const;

	/// @brief Get the mover class name.
	///
	static std::string class_name();

	/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
	///
	static void provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd );

private:


};


--end_namespace--


#endif //INCLUDED--path--_--class--_hh
