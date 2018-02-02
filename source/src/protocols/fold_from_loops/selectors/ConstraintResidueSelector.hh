// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fold_from_loops/ConstraintResidueSelector.hh
/// @brief  Selects those residues that have distance/angle/dihedral constraints
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL

#ifndef INCLUDED_protocols_fold_from_loops_ConstraintResidueSelector_HH
#define INCLUDED_protocols_fold_from_loops_ConstraintResidueSelector_HH

// Unit headers
#include <protocols/fold_from_loops/selectors/ConstraintResidueSelector.fwd.hh>

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

namespace protocols {
namespace fold_from_loops {
namespace selectors {

class ConstraintResidueSelector : public core::select::residue_selector::ResidueSelector {
public:
	// derived from base class
	ConstraintResidueSelector();

	/// @brief Copy constructor
	///
	ConstraintResidueSelector( ConstraintResidueSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	core::select::residue_selector::ResidueSelectorOP clone() const override;

	~ConstraintResidueSelector() override;

	core::select::residue_selector::ResidueSubset apply( core::pose::Pose const & pose ) const override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) override;

	std::string
	get_name() const override;

	bool distance() const { return distance_;}
	void distance( bool pick ){ distance_ = pick; };
	bool angle() const { return angle_;}
	void angle( bool pick ){ angle_ = pick; };
	bool dihedral() const { return dihedral_;}
	void dihedral( bool pick ){ dihedral_ = pick; };
	static bool default_distance() { return true; };
	static bool default_angle() { return true; };
	static bool default_dihedral() { return true; };
	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	bool distance_;
	bool angle_;
	bool dihedral_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //namespace selectors
} //namespace fold_from_loops
} //namespace protocols


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_fold_from_loops_ConstraintResidueSelector )
#endif // SERIALIZATION


#endif
