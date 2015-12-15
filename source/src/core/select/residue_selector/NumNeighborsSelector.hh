// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/NumNeighborsSelector.hh
/// @brief  The NumNeighborsSelector identifies residues that have at least X neighbors within a distance Y.
///  Clears the given ResidueSubset.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_NumNeighborsSelector_HH
#define INCLUDED_core_select_residue_selector_NumNeighborsSelector_HH

// Unit headers
#include <core/select/residue_selector/NumNeighborsSelector.fwd.hh>

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

class NumNeighborsSelector : public ResidueSelector {
public:
	// derived from base class
	NumNeighborsSelector();

	/// @brief Copy constructor
	///
	NumNeighborsSelector( NumNeighborsSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	virtual ResidueSelectorOP clone() const;

	NumNeighborsSelector( Size threshold, core::Real distance_cutoff );
	//// Undefined, commenting out to fix PyRosetta build  NumNeighborsSelector( bool count_water, Size threshold, core::Real distance_cutoff );
	virtual ~NumNeighborsSelector();

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

	bool count_water() const;
	Size threshold() const;
	Real distance_cutoff() const;
	void count_water( bool setting );
	void threshold( Size setting );
	void distance_cutoff( Size setting );

private:
	bool count_water_;
	Size threshold_;
	Real distance_cutoff_;

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
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_NumNeighborsSelector )
#endif // SERIALIZATION


#endif
