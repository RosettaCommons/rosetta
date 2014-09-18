// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/scoring/geometric_solvation/ExactOccludedHbondSolEnergy.hh
/// @brief  Solvation model based on penalizing potential for Hbonding to solvent
/// @author John Karanicolas


#ifndef INCLUDED_core_scoring_geometric_solvation_ExactOccludedHbondSolEnergy_hh
#define INCLUDED_core_scoring_geometric_solvation_ExactOccludedHbondSolEnergy_hh

#include <core/types.hh>

// Project headers
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBEvalTuple.fwd.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
// AUTO-REMOVED #include <core/chemical/AtomTypeSet.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>

// C++ headers
#include <map>

#include <core/chemical/AtomTypeSet.fwd.hh>
#include <utility/vector1.hh>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <atomic>
#include <mutex>
#endif
#endif

namespace core {
namespace scoring {
namespace geometric_solvation {


// singleton class
class GridInfo {

public:
	static GridInfo * get_instance();

	// accessors
	core::Size xnum_points() const { return xnum_points_; };
	core::Size ynum_points() const { return ynum_points_; };
	core::Size znum_points() const { return znum_points_; };
	core::Real xstep() const { return xstep_; };
	core::Real ystep() const { return ystep_; };
	core::Real zstep() const { return zstep_; };
	core::Real xorigin() const { return xorigin_; };
	core::Real yorigin() const { return yorigin_; };
	core::Real zorigin() const { return zorigin_; };

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
	//private constructor
	GridInfo();
	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static GridInfo * create_singleton_instance();

private:
	/// @brief static data member holding pointer to the singleton class itself
#if defined MULTI_THREADED && defined CXX11
	static std::atomic< GridInfo * > instance_;
#else
	static GridInfo * instance_;
#endif

	// private member data
	core::Size xnum_points_, ynum_points_, znum_points_;
	core::Real xstep_, ystep_, zstep_;
	core::Real xorigin_, yorigin_, zorigin_;

};



// singleton class
class WaterWeightGridSet {

public:
	static WaterWeightGridSet * get_instance();

	std::vector < std::vector < std::vector <core::Real> > > const &
	get_water_weight_grid( hbonds::HBEvalType const & hbond_eval_type ) const;

	core::Real
	get_sum_water_weight_grid( hbonds::HBEvalType const & hbond_eval_type ) const;

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
	//private constructor
	WaterWeightGridSet();

	core::Real fill_water_grid( std::vector < std::vector < std::vector <core::Real> > > & water_weights,
		hbonds::HBEvalTuple const & hbond_eval_type, GridInfo const & grid_info, bool const water_is_donor);

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static WaterWeightGridSet * create_singleton_instance();

private:
	/// @brief static data member holding pointer to the singleton class itself
#if defined MULTI_THREADED && defined CXX11
	static std::atomic< WaterWeightGridSet * > instance_;
#else
	static WaterWeightGridSet * instance_;
#endif


	// private member data
	std::map< hbonds::HBEvalType, std::vector < std::vector < std::vector <core::Real> > > > all_water_weights_;
	std::map< hbonds::HBEvalType, core::Real> sum_all_water_weights_;

	hbonds::HBondOptionsOP   hbondoptions_;
	hbonds::HBondDatabaseCOP hb_database_;
};

typedef std::map< hbonds::HBEvalType, std::vector < std::vector < std::vector <core::Real> > > >::const_iterator all_water_weights_iterator;
typedef std::map< hbonds::HBEvalType, core::Real>::const_iterator sum_water_weights_iterator;


class ExactOccludedHbondSolEnergy : public methods::ContextDependentOneBodyEnergy  {
public:
	typedef methods::ContextDependentOneBodyEnergy  parent;

public:

	ExactOccludedHbondSolEnergy(
		bool const exact_occ_skip_Hbonders = false,
		bool const exact_occ_include_Hbond_contribution = false,
		bool const exact_occ_pairwise = false,
		bool const exact_occ_pairwise_by_res = false,
		bool const exact_occ_split_between_res = false,
		bool const exact_occ_self_res_occ = true,
		core::Real const occ_radius_scaling = 1.,
		bool const verbose = false
	);

	ExactOccludedHbondSolEnergy( ExactOccludedHbondSolEnergy const & src );

	~ExactOccludedHbondSolEnergy();

	virtual methods::EnergyMethodOP clone() const;

	virtual void setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual void setup_for_packing(pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const;

	virtual void setup_for_derivatives( pose::Pose &pose, ScoreFunction const &  ) const;

	virtual void setup_for_minimizing(pose::Pose & pose, ScoreFunction const & , kinematics::MinimizerMapBase const &) const;

	virtual void residue_energy(
		conformation::Residue const & polar_rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	virtual void indicate_required_context_graphs( utility::vector1< bool > &  ) const {};

	virtual bool defines_intrares_energy( EnergyMap const & ) const { return false; }

	virtual Distance atomic_interaction_cutoff() const;

	core::Real
	compute_grid_constant(
	  hbonds::HBEvalTuple const & hbond_eval_type
	) const;

	core::Real compute_polar_group_sol_energy(
		pose::Pose const & pose,
		conformation::Residue const & polar_rsd,
	  Size const polar_atom,
		bool const restrict_to_single_occluding_residue = false,
		Size const single_occluding_resinx = 0,
		bool const restrict_to_single_occluding_atom = false,
		Size const single_occluding_atominx = 0
  ) const;

	core::Real compute_polar_group_sol_energy(
    pose::Pose const & pose,
		conformation::Residue const & polar_rsd,
		core::Size const polar_atomno,
		GridInfo const & grid_info,
		core::Real const & grid_constant,
		std::vector < std::vector < std::vector <core::Real> > > const & water_weights,
		bool const restrict_to_single_occluding_residue = false,
		core::Size const single_occluding_resinx = 0,
		bool const restrict_to_single_occluding_atom = false,
		core::Size const single_occluding_atominx = 0
	) const;


private:

	bool const exact_occ_skip_Hbonders_;
	bool const exact_occ_include_Hbond_contribution_;
	bool const exact_occ_pairwise_;
	bool const exact_occ_pairwise_by_res_;
	bool const exact_occ_split_between_res_;
	bool const exact_occ_self_res_occ_;
	core::Real const occ_radius_scaling_;
	hbonds::HBondOptionsOP   hbondoptions_;
	hbonds::HBondDatabaseCOP hb_database_;
	bool const verbose_;
	chemical::AtomTypeSetCAP atom_type_set_ptr_;


	// note: this would really be a local variable, it's just that we don't want to pay
	// the price to reallocate memory every time. Making it member data keeps it persistent,
	// but then we can't alter it in const member functions. For this reason, it's mutable...
	mutable std::vector < std::vector < std::vector <bool> > > occluded_sites_;

	virtual
	core::Size version() const;

};

} // geometric_solvation
} // scoring
} // core

#endif

