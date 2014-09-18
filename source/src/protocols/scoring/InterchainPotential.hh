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
#include <core/scoring/EnvPairPotential.hh>

// Package headers
// AUTO-REMOVED #include <protocols/scoring/Interface.hh>
#include <core/conformation/Residue.fwd.hh>

#include <basic/datacache/CacheableData.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

#include <protocols/scoring/InterfaceInfo.fwd.hh>

// Utility headers
#include <utility/vector1.hh>


#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <atomic>
#include <mutex>
#endif
#endif

namespace protocols {
namespace scoring {

using namespace core::conformation;

class InterchainPotential : public core::scoring::EnvPairPotential {

public:

	static InterchainPotential * get_instance();

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

#ifdef MULTI_THREADED
#ifdef CXX11
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:
	static std::mutex singleton_mutex_;
#endif
#endif

private:

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static InterchainPotential * create_singleton_instance();

	InterchainPotential();
	InterchainPotential( InterchainPotential const & src );
	InterchainPotential & operator = ( InterchainPotential const & rhs );

#if defined MULTI_THREADED && defined CXX11
	static std::atomic< InterchainPotential * > instance_;
#else
	static InterchainPotential * instance_;
#endif

	// const-ref to scoring database
	core::scoring::AtomVDW const & atom_vdw_;

	ObjexxFCL::FArray2D< core::Real > interchain_env_log_;
	ObjexxFCL::FArray2D< core::Real > interchain_pair_log_;
};

} // ns scoring
} // ns core

#endif
