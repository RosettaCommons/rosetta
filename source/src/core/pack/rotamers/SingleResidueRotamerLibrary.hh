// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamers/SingleResidueRotamerLibrary.hh
/// @brief  SingleResidueRotamerLibrary virtual base class
/// @author Andrew Leaver-Fay
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDE_core_pack_rotamers_SingleResidueRotamerLibrary_hh
#define INCLUDE_core_pack_rotamers_SingleResidueRotamerLibrary_hh

// Unit Headers
#include <core/pack/rotamers/SingleResidueRotamerLibrary.fwd.hh>

// Package Headers
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh> // ChiVector

//Project Headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/io/izstream.fwd.hh>

// Numeric headers
#include <numeric/random/random.fwd.hh>

namespace core {
namespace pack {
namespace rotamers {

/// @brief  SingleResidueRotamerLibrary pure virtual base class
class SingleResidueRotamerLibrary : public utility::pointer::ReferenceCount
{
public:
	virtual
	~SingleResidueRotamerLibrary();

	virtual
	Real
	rotamer_energy_deriv(
		conformation::Residue const & rsd,
		dunbrack::RotamerLibraryScratchSpace & scratch
	) const = 0;

	virtual
	Real
	rotamer_energy(
		conformation::Residue const & rsd,
		dunbrack::RotamerLibraryScratchSpace & scratch
	) const = 0;

	/// @brief Returns the energy of the lowest-energy rotamer accessible to the given residue
	/// (based on e.g. its current phi and psi values).
	/// If curr_rotamer_only is true, then consider only the idealized version of the
	/// residue's current rotamer (local optimum); otherwise, consider all rotamers (global optimum).
	virtual
	Real
	best_rotamer_energy(
		conformation::Residue const & rsd,
		bool curr_rotamer_only,
		dunbrack::RotamerLibraryScratchSpace & scratch
	) const = 0;

	/// @brief Pick a rotamer for the input residue according to the rotamer probability
	/// distribution and assign chi angles to the input rsd.  If perturb_from_rotamer_center
	/// is true, then push the rotamer off from the center; for chi angles with a normal
	/// distribution, the perturbation is taken from a Gaussian random number with a standard
	/// deviation matching the chi angle's standard deviation.  For chi angles that are not
	/// normally distributed, the behavior is open to the derived classe's interpretation.
	virtual
	void
	assign_random_rotamer_with_bias(
		conformation::Residue const & rsd,
		pose::Pose const & pose, // DOUG MAY CAUSE PROBLEMS BUILDONG
		dunbrack::RotamerLibraryScratchSpace & scratch,
		numeric::random::RandomGenerator & RG,
		dunbrack::ChiVector & new_chi_angles,
		bool perturb_from_rotamer_center
	) const = 0;

	virtual
	void
	fill_rotamer_vector(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		pack::task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const& existing_residue,
		utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
		bool buried,
		RotamerVector & rotamers
	) const = 0;

	/// @brief Filter a RotamerVector by "bump energy" of a rotamer:
	/// All rotamers with bump energies over a certain threshold will be discarded
	/// Exception: if all rotamers are over the threshold, one rotamer (with the lowest
	/// bump energy) will be reserved.
	/// The vector "rotamers" will be modified "in-place"
	virtual
	void
	bump_filter(
		RotamerVector & rotamers,
		core::Size resid,
		scoring::ScoreFunction const & sf,
		pose::Pose const & pose,
		task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph
	) const;

	/// @brief Computes the "bump energy" of a rotamer: the bump energy is the
	/// sum of rotamer's interactions with 1) the backbone-and-side chains of
	/// neighboring residues that are held fixed during this repacking optimization
	/// and 2) the backbones of neighboring residues that are changable during this
	/// repacking optimization.
	virtual
	core::PackerEnergy
	bump_check(
		core::conformation::ResidueCOP rotamer,
		core::Size resid,
		scoring::ScoreFunction const & sf,
		pose::Pose const & pose,
		task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph
	) const;

	/// @brief Adds the current rotamer to rotamer vector, if the Rotlib supports it
	/// @details This is in this class mainly because of historical
	/// behavior of certain rotamer libraries not supporting current rotamers
	virtual
	core::Size
	current_rotamer(
		RotamerVector & rotamers,
		core::Size resid,
		task::PackerTask const & task,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const & existing_residue
	) const;

	/// @brief Generate an "emergency rotamer" if we don't have any
	/// @details This is in this class mainly because of historical
	/// behavior of certain rotamer libraries not supporting current rotamers
	virtual
	void
	emergency_rotamer(
		RotamerVector & rotamers,
		core::Size resid,
		pose::Pose const & pose,
		task::PackerTask const & task,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const & existing_residue
	) const;

	/// @brief Add a virtualized sidechain to the rotamer vector if
	/// settings call for it.
	RotamerVector
	virtual_sidechain(
		RotamerVector const & rotamers,
		core::Size resid,
		pose::Pose const & pose,
		task::PackerTask const & task,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const & existing_residue
	) const;

	//XRW_B_T1
	/*
	virtual
	SingleResidueRotamerLibraryOP
	coarsify(coarse::Translator const &map) const = 0;
	*/
	//XRW_E_T1

	virtual
	void
	write_to_file( utility::io::ozstream &out ) const = 0;

	/// @brief Equality test for equivalence.
	/// Two SingleResidueRotamerLibraries test equal if and only if they represent the exact same behavior
	virtual
	bool
	operator ==( SingleResidueRotamerLibrary const & ) const;

};

} // namespace rotamers
} // namespace pack
} // namespace core

#endif
