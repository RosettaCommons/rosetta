// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.hh
/// @brief  rotamer set implementation class for symmetric packing
/// @author Ingemar Andre


#ifndef INCLUDED_core_pack_rotamer_set_symmetry_SymmetricRotamerSet__hh
#define INCLUDED_core_pack_rotamer_set_symmetry_SymmetricRotamerSet__hh

//Unit headers
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.fwd.hh>

//Package headers
#include <core/pack/rotamer_set/RotamerSet_.hh>

//Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>

#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace rotamer_set {
namespace symmetry {

/// @brief Container for a set of rotamers for use in packing.
/// Rotamers are sorted into groups of the same residue type.
/// Offsets into these rotamer groups are maintained by this class, as is
/// information concerning the "original rotamer" -- the rotamer
/// present on the input pose before packing began.
/// symmetrical version of RotamerSet_
class SymmetricRotamerSet_ : public RotamerSet_
{
public:
	typedef conformation::symmetry::SymmetricConformation SymmetricConformation;
	typedef conformation::symmetry::SymmetryInfoCOP SymmetryInfoCOP;

public:
	SymmetricRotamerSet_();
	virtual ~SymmetricRotamerSet_();

	/// @brief Computes the packers "one body energies" for the set of rotamers.
	virtual
	void
	compute_one_body_energies(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,
		utility::vector1< core::PackerEnergy > & energies ) const;

	void
	PackerEnergyMultiply(
		utility::vector1< core::PackerEnergy > & energies,
		core::Size factor ) const;

	void
	PackerEnergyAdd(
		utility::vector1< core::PackerEnergy > & energies,
		utility::vector1< core::PackerEnergy > const & add ) const;

	void
	PackerEnergySubtract(
		utility::vector1< core::PackerEnergy > & energies,
		utility::vector1< core::PackerEnergy > const & subtract ) const;

	RotamerSetOP
	orient_rotamer_set_to_symmetric_partner(
		pose::Pose const & pose,
		conformation::Residue const & residue_in,
		int const & sympos,
		RotamerSet const & rotset_in
	) const;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace symmetry
} // namespace rotamer_set
} // namespace pack
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_rotamer_set_symmetry_SymmetricRotamerSet_ )
#endif // SERIALIZATION


#endif // INCLUDED_core_pack_RotamerSet_RotamerSet__HH

