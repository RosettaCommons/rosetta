// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/buns/BuriedUnsatPolarsFinder2.hh
/// @brief
/// @details
/// @author Kevin Houlihan (khouli@unc.edu)
/// @author Bryan Der

#ifndef INCLUDED_devel_buns_BuriedUnsatisfiedPolarsCalculator2_hh
#define INCLUDED_devel_buns_BuriedUnsatisfiedPolarsCalculator2_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>

#include <basic/options/option.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <set>

#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace devel {
namespace buns {

class BuriedUnsatisfiedPolarsCalculator2 : public core::pose::metrics::EnergyDependentCalculator {

public:

	BuriedUnsatisfiedPolarsCalculator2(
		std::string weak_bunsat_calc_
	);

	BuriedUnsatisfiedPolarsCalculator2(
		std::string weak_bunsat_calc,
		std::set< Size > special_region
	);

	core::pose::metrics::PoseMetricCalculatorOP clone()
	const override
	{
		return core::pose::metrics::PoseMetricCalculatorOP( new BuriedUnsatisfiedPolarsCalculator2( name_of_weak_bunsat_calc_ ) );
	};

	std::string const & name_of_weak_bunsat_calc()
	const
	{
		return name_of_weak_bunsat_calc_;
	}

	void
	set_layered_sasa(bool val)
	{
		layered_sasa_ = val;
		this->notify_energy_change();
	}

	void
	set_generous_hbonds(bool val)
	{
		generous_hbonds_ = val;
		this->notify_energy_change();
	}

	void
	set_sasa_burial_cutoff(core::Real const & val)
	{
		sasa_burial_cutoff_ = val;
		this->notify_energy_change();
	}

	void
	set_AHD_cutoff(core::Real const & val)
	{
		AHD_cutoff_ = val;
		this->notify_energy_change();
	}

	void
	set_dist_cutoff(core::Real const & val)
	{
		dist_cutoff_ = val;
		this->notify_energy_change();
	}

	void
	set_hxl_dist_cutoff(core::Real const & val)
	{
		hxl_dist_cutoff_ = val;
		this->notify_energy_change();
	}

	void
	set_sulph_dist_cutoff(core::Real const & val)
	{
		sulph_dist_cutoff_ = val;
		this->notify_energy_change();
	}

	void
	set_metal_dist_cutoff(core::Real const & val)
	{
		metal_dist_cutoff_ = val;
		this->notify_energy_change();
	}

protected:
	void lookup( std::string const & key, basic::MetricValueBase * valptr ) const override;
	std::string print( std::string const & key ) const override;
	void recompute( core::pose::Pose const & this_pose ) override;

private:
	virtual void generous_hbond() const;

	virtual void bunsats_thorough_check(
		core::pose::Pose const & pose,
		core::id::AtomID_Map< bool > & bunsat_thorough_atomid_map);

	virtual
	bool
	single_bunsat_thorough_check(
		core::pose::Pose const & pose,
		core::id::AtomID const & bunsat_candidate_atom_id
	);

	virtual
	void
	bunsat_donor_nbr_residue_check(
		core::pose::Pose const & pose,
		core::id::AtomID const & bunsat_candidate_atom_id,
		core::conformation::Residue const & bunsat_rsd,
		numeric::xyzVector< core::Real > const & bunsat_xyz,
		Size const test_resi,
		Size & num_hbonds
	);

	virtual
	void
	bunsat_acc_nbr_residue_check(
		core::pose::Pose const & pose,
		core::id::AtomID const & bunsat_candidate_atom_id,
		core::conformation::Residue const & bunsat_rsd,
		numeric::xyzVector< core::Real > const & bunsat_xyz,
		Size const & test_resi,
		Size & num_hbonds
	);

	bool
	metal_check(
		core::conformation::Residue const & test_rsd,
		numeric::xyzVector< core::Real > const & bunsat_xyz,
		numeric::xyzVector< core::Real > const & test_xyz
	) const;

	bool
	adjacent_bbbb_check(
		Size const & bunsat_resi,
		std::string const & bunsat_atom_name,
		Size const & test_resi,
		std::string const & test_atom_name
	) const;

	bool
	self_scsc(
		core::conformation::Residue const & bunsat_rsd,
		core::Size const & bunsat_resi,
		Size const & bunsat_atom_num,
		core::conformation::Residue const & test_rsd,
		Size const & test_resi,
		Size const & test_atom_num
	) const;

	bool
	sulphur_bond_check(
		core::conformation::Residue const & test_rsd,
		Size const & test_atom_num,
		numeric::xyzVector< core::Real > const & bunsat_xyz,
		numeric::xyzVector< core::Real > const & test_xyz
	) const;

	bool
	don_geom_check(
		core::pose::Pose const & pose,
		core::conformation::Residue const & bunsat_rsd,
		Size const & bunsat_atom_num,
		numeric::xyzVector< core::Real > const & bunsat_xyz,
		numeric::xyzVector< core::Real > const & test_xyz
	) const;

	bool
	acc_geom_check(
		core::pose::Pose const & pose,
		numeric::xyzVector< core::Real > const & bunsat_xyz,
		core::conformation::Residue const & test_rsd,
		Size const & test_atom_num,
		numeric::xyzVector< core::Real > const & test_xyz
	) const;


private:
	void assert_calculators();

	void
	show();

	Size all_bur_unsat_polars_;
	Size special_region_bur_unsat_polars_;
	core::id::AtomID_Map< bool > atom_bur_unsat_;
	utility::vector1< Size > residue_bur_unsat_polars_;

private:
	static
	core::Size satisfaction_cutoff( std::string atom_type );
	std::string name_of_weak_bunsat_calc_;
	std::set< Size > special_region_;
	bool layered_sasa_;
	bool generous_hbonds_;
	core::Real sasa_burial_cutoff_;
	core::Real AHD_cutoff_;
	core::Real dist_cutoff_;
	core::Real hxl_dist_cutoff_;
	core::Real sulph_dist_cutoff_;
	core::Real metal_dist_cutoff_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	BuriedUnsatisfiedPolarsCalculator2();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // BuriedUnsatisfiedPolarsCalculator2

} // namespace buns
} // namespace devel

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( devel_buns_BuriedUnsatisfiedPolarsCalculator2 )
#endif // SERIALIZATION


#endif
