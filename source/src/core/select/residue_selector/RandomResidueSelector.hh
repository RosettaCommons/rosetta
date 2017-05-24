// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/RandomResidueSelector.hh
/// @brief  The RandomResidueSelector selects residues using a string containing residue names
/// @author Tom Linsky (tlinsky@uw.edu))

#ifndef INCLUDED_core_select_residue_selector_RandomResidueSelector_HH
#define INCLUDED_core_select_residue_selector_RandomResidueSelector_HH

// Unit headers
#include <core/select/residue_selector/RandomResidueSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueVector.fwd.hh>
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

/// @brief The RandomResidueSelector returns a ResidueSubset, i.e. a utility::vector1< bool > containing
/// 'true' for N randomly selected residue positions, where N is a user-specified integer.
class RandomResidueSelector : public ResidueSelector {
public:
	RandomResidueSelector();
	RandomResidueSelector( ResidueSelectorCOP selector, Size const num_residues );
	RandomResidueSelector( ResidueSelectorCOP selector, Size const num_residues, bool const select_res_cluster, Real const distance_cutoff );
	virtual ~RandomResidueSelector();

	virtual ResidueSelectorOP
	clone() const;

	virtual ResidueSubset
	apply( core::pose::Pose const & pose ) const;

	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data );

	virtual std::string
	get_name() const;

	static std::string
	class_name();

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	ResidueSubset
	subset_from_randomized_vector( core::pose::Pose const & pose, ResidueVector const & random_order_residue_set ) const;

private: // data members
	ResidueSelectorCOP selector_;
	Size num_residues_;
	bool select_res_cluster_;
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

#endif
