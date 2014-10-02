// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/InterchainEnergy.cc
/// @brief  Statistically derived rotamer pair potentials
/// @detailed For docking (or between chains) only those residues at the interface
///						and between the two interfaces need to be evaluated
/// @author Monica Berrondo

#ifndef INCLUDED_protocols_scoring_InterchainPotential_hh
#define INCLUDED_protocols_scoring_InterchainPotential_hh

#include <core/types.hh>

// Unit headers
#include <core/scoring/AtomVDW.fwd.hh>
#include <protocols/scoring/InterchainPotential.fwd.hh>
// AUTO-REMOVED #include <protocols/scoring/InterfaceInfo.hh>

// Package headers
// AUTO-REMOVED #include <protocols/scoring/Interface.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <basic/datacache/CacheableData.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

#include <protocols/scoring/InterfaceInfo.fwd.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace scoring {

using namespace core::conformation;

/// @brief Singleton class to hold the interface-derived statistics for residue-pair
/// scores at protein/protein interfaces.
/// @details This previously derived from the EnvPairPotential, which is in no way
/// necessary because the two classes have nothing in common; rather, the
/// InterchainPairEnergy and InterchainEnvEnergy classes can hold a pointer to both
/// the InterchainPotential and the EnvPairPotential.
class InterchainPotential : public utility::SingletonBase< InterchainPotential > {
public:
	friend class utility::SingletonBase< InterchainPotential >;

public:
	void
	compute_interface( core::pose::Pose & pose ) const;

	void
	finalize( core::pose::Pose & pose ) const;

	///
	void
	evaluate_env_score(
		core::pose::Pose const & pose,
		core::conformation::Residue const & rsd,
		core::Real & env_score
	) const;

	///
	void
	evaluate_contact_score(
		core::pose::Pose const & pose,
		core::Real & contact_score
	) const;

	///
	void
	evaluate_pair_and_vdw_score(
		core::pose::Pose const & pose,
		core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2,
		core::Real & pair_score,
		core::Real & vdw_score
	) const;

	// Commention out to make PyRosetta compile (undefined in .cc file)
	//core::Size interface_residues( core::pose::Pose const & pose ) const;

	InterfaceInfo const & interface_from_pose( core::pose::Pose const & ) const;
	InterfaceInfo & nonconst_interface_from_pose( core::pose::Pose & ) const;

private:

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static InterchainPotential * create_singleton_instance();

	InterchainPotential();
	InterchainPotential( InterchainPotential const & src );
	InterchainPotential & operator = ( InterchainPotential const & rhs );

	// const-ref to scoring database
	core::scoring::AtomVDW const & atom_vdw_;

	ObjexxFCL::FArray2D< core::Real > interchain_env_log_;
	ObjexxFCL::FArray2D< core::Real > interchain_pair_log_;
};

} // ns scoring
} // ns core

#endif
