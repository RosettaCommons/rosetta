// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ShortRangeTwoBodyEnergy.hh
/// @brief  Short Range Two Body Energy Method base class declaration
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_scoring_methods_ShortRangeTwoBodyEnergy_hh
#define INCLUDED_core_scoring_methods_ShortRangeTwoBodyEnergy_hh

#include <core/scoring/methods/TwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

class ShortRangeTwoBodyEnergy : public TwoBodyEnergy
{
public:
	typedef TwoBodyEnergy parent;

public:

	/// @brief Constructor with EnergyMethodCreator to provide to the EnergyMethod
	/// grandparent the list of the ScoreTypes this EnergyMethod is responsible for
	/// computing
	ShortRangeTwoBodyEnergy( EnergyMethodCreatorOP );

	virtual ~ShortRangeTwoBodyEnergy();


	/// @brief how far apart must two heavy atoms be to have a zero interaction energy?
	///
	/// @details If hydrogen atoms interact at the same range as heavy atoms, then
	/// this distance should build-in a 2 * max-bound-h-distance-cutoff buffer.
	/// There is an improper mixing here between run-time aquired chemical knowledge
	/// (max-bound-h-distance-cutoff) and compile time aquired scoring knowledge
	/// (max atom cutoff); this could be resolved by adding a boolean
	/// uses_hydrogen_interaction_distance() to the SRTBEnergy class along with
	/// a method of the ChemicalManager max_bound_h_distance_cutoff().
	virtual
	Distance
	atomic_interaction_cutoff() const = 0;


	/// @brief A derived class should return true for this function if it implements its own
	/// versions of the backbone_backbone_energy, backbone_sidechain_energy and
	/// sidechain_sidechain_energy functions.  The default sidechain_sidechain_energy implemented
	/// by the TwoBodyEnergy base class calls residue_pair_energy.  If the derived class implements its own
	/// versions of these functions, then calling code may avoid calling it on pairs of residues
	/// that are "provably distant" based on a pair of bounding spheres for a sidechains and
	/// backbones and this method's atomic_interaction_cutoff energy method.
	virtual
	bool
	divides_backbone_and_sidechain_energetics() const;


	/// @brief Batch computation of rotamer pair energies.  Need not be overriden in
	/// derived class -- by default, iterates over all pairs of rotamers,
	/// and calls derived class's residue_pair_energy method.  Since short range rotamer pairs
	/// may not need calculation, the default method looks at blocks of residue type pairs
	/// and only calls the residue_pair_energy method if the rotamer pairs are within range
	virtual
	void
	evaluate_rotamer_pair_energies(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
	) const;


	/// @brief Batch computation of rotamer/background energies.  Need not be overriden
	/// in derived class -- by default, iterates over all rotamers in the set, and calls
	/// derived class's residue_pair_energy method for each one against the background rotamer
	/// Since short range rotamer pairs may not need calculation, the default method
	/// looks at blocks of residue type pairs and only calls the residue_pair_energy method
	/// if the rotamer pairs are within range
	virtual
	void
	evaluate_rotamer_background_energies(
		conformation::RotamerSetBase const & set,
		conformation::Residue const & residue,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		utility::vector1< core::PackerEnergy > & energy_vector
	) const;

	/// @brief Batch computation of rotamer/background energies.  Need not be overriden
	/// in derived class -- by default, iterates over all rotamers in the set, and calls
	/// derived class's residue_pair_energy method for each one against the background rotamer
	/// Since short range rotamer pairs may not need calculation, the default method
	/// looks at blocks of residue type pairs and only calls the residue_pair_energy method
	/// if the rotamer pairs are within range
	virtual
	void
	evaluate_rotamer_background_energy_maps(
		conformation::RotamerSetBase const & set,
		conformation::Residue const & residue,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		utility::vector1< EnergyMap > & emaps
	) const;


};


}
}
}

#endif
