// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/ResidueInMembraneSelector.hh
/// @brief  The ResidueInMembraneSelector selects residues either in or out of the membrane
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)

#ifndef INCLUDED_core_select_residue_selector_ResidueInMembraneSelector_HH
#define INCLUDED_core_select_residue_selector_ResidueInMembraneSelector_HH

// Unit headers
#include <core/select/residue_selector/ResidueInMembraneSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>

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

/// @brief The ResidueInMembraneSelector returns a ResidueSubset, i.e. a utility::vector1< bool > containing
/// 'true' for residue positions which are either in or out of the membrane. default is in.
class ResidueInMembraneSelector : public ResidueSelector {
public:
	// derived from base class
	ResidueInMembraneSelector();
	//ResidueInMembraneSelector( std::string const & res_name_str );
	virtual ~ResidueInMembraneSelector();

	/// @brief Copy constructor
	///
	ResidueInMembraneSelector( ResidueInMembraneSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	virtual ResidueSelectorOP clone() const;

	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	);

	virtual
	std::string
	get_name() const;

	static std::string class_name();


	// //unit-specific
	// /// @brief sets the comma-separated string of residue names to be selected
	// void set_residue_names( std::string const & res_name_str );
	//
	// /// @brief sets the comma-separated string of 3-character residue names to be selected
	// void set_residue_name3( std::string const & res_name3_str );

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	void set_select_in_membrane( bool const & select_in_membrane );

private: // data members
	bool select_in_membrane_;

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
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_ResidueInMembraneSelector )
#endif // SERIALIZATION


#endif
