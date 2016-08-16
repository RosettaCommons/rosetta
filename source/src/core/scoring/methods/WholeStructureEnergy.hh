// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/methods/WholeStructureEnergy.hh
/// @brief EnergyMethod that applies to a whole structure and not on a residue-by-residue basis
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_scoring_methods_WholeStructureEnergy_hh
#define INCLUDED_core_scoring_methods_WholeStructureEnergy_hh

// Unit headers
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/EnergyMethod.hh>

// Project headers
#include <core/types.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

/// @brief Base class for EnergyMethods which are meaningful only on entire structures,
/// for example, the Radius of Gyration.  These EnergyMethods do all of their work in
/// the "finalize_total_energy" section of score function evaluation.
class WholeStructureEnergy : public EnergyMethod {
public:
	typedef EnergyMethod parent;

public:

	/// @brief Constructor with EnergyMethodCreator to list the ScoreTypes
	/// computed by this WholeStructureEnergy.
	WholeStructureEnergy( EnergyMethodCreatorOP );

	virtual ~WholeStructureEnergy() {}


	EnergyMethodType
	method_type() const;

	/// @brief how far apart must two heavy atoms be to have a zero interaction energy?
	///
	/// @details If hydrogen atoms interact at the same range as heavy atoms, then
	/// this distance should build-in a 2 * max-bound-h-distance-cutoff buffer.
	/// There is an improper mixing here between run-time aquired chemical knowledge
	/// (max-bound-h-distance-cutoff) and compile time aquired scoring knowledge
	/// (max atom cutoff); this could be resolved by adding a boolean
	/// uses_hydrogen_interaction_distance() to the SRTBEnergy class along with
	/// a method of the ChemicalManager max_bound_h_distance_cutoff().
	///
	/// This method allows the WholeStructureEnergy class to define which edges
	/// should be included in the EnergyGraph so that during the finalize() method
	/// the Energy class can iterate across the EnergyGraph.  This iteration occurrs
	/// in the SecondaryStructureEnergy class, where the edges must span 12 angstroms
	/// between the centroids.  Arguably, the SecondaryStructureEnergy class could use
	/// the TwelveANeighborGraph (a context graph) and not require that the EnergyGraph
	/// span such long distances.
	virtual
	Distance
	atomic_interaction_cutoff() const { return 0.0; }

};

} // methods
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
