// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/ChainSelector.hh
/// @brief  The ChainSelector selects all the residues from a given chain, given either as a number
///         or a character.  Clears the passed ResidueSubset.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_ChainSelector_HH
#define INCLUDED_core_select_residue_selector_ChainSelector_HH

// Unit headers
#include <core/select/residue_selector/ChainSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <list>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

class ChainSelector : public ResidueSelector {
public:
	// derived from base class
	ChainSelector();

	/// @brief Copy constructor
	///
	ChainSelector( ChainSelector const &src);

	// takes a comma-separated list of chains (e.g. "A,B")
	ChainSelector( std::string chains );
	virtual ~ChainSelector();

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	virtual ResidueSelectorOP clone() const;

	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	virtual
	std::string
	get_name() const;

	static std::string class_name();
	static void provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd );

	utility::vector1< std::string > const &
	chain_strings() const;

	void set_chain_strings( utility::vector1< std::string > const & );

private:

	void select_chain_by_index( core::pose::Pose const &, ResidueSubset &, core::Size ) const;
	void select_chain_by_pdb_chain_char( core::pose::Pose const &, ResidueSubset &, char ) const;

private: // data members

	utility::vector1< std::string > chain_strings_;

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
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_ChainSelector )
#endif // SERIALIZATION


#endif
