// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/annealing/ResidueArrayAnnealableEnergy.hh
/// @brief  Annealable method interrface for score types evaluated over explicit list of residues.
/// @author Alex Ford (fordas@uw.edu)

#ifndef INCLUDED_core_scoring_annealing_ResidueArrayAnnealableEnergy_hh
#define INCLUDED_core_scoring_annealing_ResidueArrayAnnealableEnergy_hh

#include <core/scoring/annealing/ResidueArrayAnnealableEnergy.fwd.hh>
#include <core/scoring/annealing/RotamerSets.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/types.hh>

namespace core {
namespace scoring {
namespace annealing {

class ResidueArrayAnnealableEnergy
{
public:

	/// @brief Constructor.
	///
	ResidueArrayAnnealableEnergy();

	/// @brief Copy constructor.
	///
	ResidueArrayAnnealableEnergy( ResidueArrayAnnealableEnergy const &src );

	/// @brief Destructor.
	///
	virtual ~ResidueArrayAnnealableEnergy();

	/// @brief Calculate the energy given a vector of const-owning pointers to Residue objects.
	/// @details Must be implemented by derived classes.
	virtual core::Real calculate_energy( utility::vector1< core::conformation::ResidueCOP > const &resvect, core::Size const substitution_position = 0 ) const = 0;

	/// @brief ResidueArrayAnnealableEnergy objects may optionally cache data within the EnergyMethod prior to a packer run.
	/// This function is defined as doing nothing by default, but can be redefined on a per-EnergyMethod basis to cache whatever
	/// data are necessary.
	/// @details Note that this is generally intended so that data can be cached FROM the Pose or ScoreFunction, not TO the Pose or ScoreFunction.
	/// There are exceptions, though: sometimes it is necessary to cache a large, reusable object within the Pose itself (e.g. in the datacache
	/// of the Energies object in the Pose), and for that reason the pose is nonconst, here.
	/// @note The method is nonconst to permit caching within the ResidueArrayAnnealableEnergy-derived EnergyMethod.
	virtual void set_up_residuearrayannealableenergy_for_packing ( core::pose::Pose &pose, core::pack::rotamer_set::RotamerSets const &rotamersets, core::scoring::ScoreFunction const &sfxn);

	/// @brief Allows the ResidueArrayAnnealableEnergy to clean up cached data, either within the EnergyMethod or in the pose, after
	/// a packer run.
	/// @details Base class version does nothing; may be overridden by derived classes.
	virtual void clean_up_residuearrayannealableenergy_after_packing( core::pose::Pose &pose );

	/// @brief What to do when a substitution that was considered is accepted.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	virtual void commit_considered_substitution();

};

}
}
}

#endif
