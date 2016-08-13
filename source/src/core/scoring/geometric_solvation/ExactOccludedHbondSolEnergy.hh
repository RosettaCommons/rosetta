// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/geometric_solvation/ExactOccludedHbondSolEnergy.hh
/// @brief  Solvation model based on penalizing potential for Hbonding to solvent
/// @author John Karanicolas


#ifndef INCLUDED_core_scoring_geometric_solvation_ExactOccludedHbondSolEnergy_hh
#define INCLUDED_core_scoring_geometric_solvation_ExactOccludedHbondSolEnergy_hh

#include <core/scoring/geometric_solvation/ExactOccludedHbondSolEnergy.fwd.hh>
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/etable/EtableEnergy.fwd.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBEvalTuple.fwd.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/types.hh>
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>
#include <map>

namespace core {
namespace scoring {
namespace geometric_solvation {

// singleton class
class GridInfo : public utility::SingletonBase< GridInfo >
{
public:
	friend class utility::SingletonBase< GridInfo >;

public:

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

private:

	//private constructor
	GridInfo();

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static GridInfo * create_singleton_instance();

private:

	// private member data
	core::Size xnum_points_, ynum_points_, znum_points_;
	core::Real xstep_, ystep_, zstep_;
	core::Real xorigin_, yorigin_, zorigin_;
};


// singleton class
class WaterWeightGridSet : public utility::SingletonBase< WaterWeightGridSet >
{

public:

	typedef std::vector < std::vector < std::vector <core::Real> > > Grid;

	friend class utility::SingletonBase< WaterWeightGridSet >;

	Grid const &
	get_water_weight_grid( hbonds::HBEvalType const & hbond_eval_type ) const;

	core::Real
	get_sum_water_weight_grid( hbonds::HBEvalType const & hbond_eval_type ) const;

	/// @brief prints a given xz-plane of a water grid
	void print_water_weight_grid_xz_plane(hbonds::HBEvalType const & hbond_eval_type, int const y) const;

private:

	//private constructor
	WaterWeightGridSet();

	core::Real fill_water_grid( Grid & water_weights,
		hbonds::HBEvalTuple const & hbond_eval_type, GridInfo const & grid_info, bool const water_is_donor);

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static WaterWeightGridSet * create_singleton_instance();

private:

	// private member data
	std::map< hbonds::HBEvalType, Grid > all_water_weights_;
	std::map< hbonds::HBEvalType, core::Real> sum_all_water_weights_;

	hbonds::HBondOptionsOP   hbondoptions_;
	hbonds::HBondDatabaseCOP hb_database_;
};

typedef std::map< hbonds::HBEvalType, WaterWeightGridSet::Grid >::const_iterator all_water_weights_iterator;
typedef std::map< hbonds::HBEvalType, core::Real >::const_iterator sum_water_weights_iterator;


class ExactOccludedHbondSolEnergy : public methods::ContextDependentOneBodyEnergy  {

public:

	typedef methods::ContextDependentOneBodyEnergy  parent;

	ExactOccludedHbondSolEnergy(
		etable::Etable const & etable_in,
		bool const analytic_etable_evaluation,
		bool const exact_occ_skip_Hbonders = false,
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

	void init_hbond_data( pose::Pose const& pose) const;

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

	/// @brief computes the desolvation energy (i.e., either SHO or median LK) of a donor atom
	core::Real compute_donor_atom_energy(
		conformation::Residue const & polar_rsd,
		core::Size polar_resnum,
		core::Size const polar_atom,
		pose::Pose const & pose
	) const;

	/// @brief computes the desolvation energy (i.e., either SHO or median LK) of an acceptor atom
	core::Real compute_acceptor_atom_energy(
		conformation::Residue const & polar_rsd,
		core::Size polar_resnum,
		core::Size const polar_atom,
		pose::Pose const & pose
	) const;

	/// @brief computes the SHO energy of a donor atom
	core::Real compute_sho_donor_atom_energy(
		conformation::Residue const & polar_rsd,
		core::Size polar_resnum,
		core::Size const polar_atom,
		pose::Pose const & pose) const;

	/// @brief computes the SHO energy of an acceptor atom
	core::Real compute_sho_acceptor_atom_energy(
		conformation::Residue const & polar_rsd,
		core::Size polar_resnum,
		core::Size const polar_atom,
		pose::Pose const & pose) const;

	/// @brief computes the SHO energy of a polar group
	core::Real compute_polar_group_sol_energy(
		pose::Pose const & pose,
		conformation::Residue const & polar_rsd,
		Size const polar_atom,
		bool const restrict_to_single_occluding_residue = false,
		Size const single_occluding_resinx = 0,
		bool const restrict_to_single_occluding_atom = false,
		Size const single_occluding_atominx = 0
	) const;

	/// @brief returns the LK energy of a given atom due to all its neighboring residues
	core::Real get_atom_lk_energy(
		core::Size const atom_idx,
		core::conformation::Residue const& res,
		core::pose::Pose const& ps ) const;

	hbonds::HBondSetCOP hbond_set() const {return hbond_set_;}

	virtual
	core::Size version() const;

private:

	void allocate_grid_of_occluded_sites();

	/// @brief computes energy of fully buried polar group
	core::Real compute_fully_buried_ene() const;

	/// @brief computes the grid constant for a given polar group (i.e., the denominator in the solvation energy equation)
	core::Real compute_grid_constant(
		core::scoring::hbonds::HBEvalTuple const & hbond_eval_type,
		core::Real fully_buried_ene ) const;

	/// @brief computes the SHO energy of a polar group
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

	/// @brief returns the LK energy of a given atom due to a given residue
	core::Real get_atom_lk_energy(
		core::Size atom_idx,
		conformation::Residue const& res,
		conformation::Residue const& occ_res) const;

	bool const exact_occ_skip_Hbonders_;
	bool const exact_occ_pairwise_;
	bool const exact_occ_pairwise_by_res_;
	bool const exact_occ_split_between_res_;
	bool const exact_occ_self_res_occ_;
	core::Real const occ_radius_scaling_;
	hbonds::HBondOptionsOP   hbondoptions_;
	hbonds::HBondDatabaseCOP hb_database_;
	chemical::AtomTypeSetCAP atom_type_set_ptr_;

	// note: this would really be a local variable, it's just that we don't want to pay
	// the price to reallocate memory every time. Making it member data keeps it persistent,
	// but then we can't alter it in const member functions. For this reason, it's mutable...
	mutable WaterWeightGridSet::Grid occluded_sites_;

	/// @brief stores H-bond info for the entire pose.
	core::scoring::hbonds::HBondSetOP hbond_set_;

	/// @brief computes LK energies for H-bonded atoms
	etable::EtableEvaluatorOP etable_evaluator_;

	/// @brief interatomic squared distance below which LK energies are computed
	core::Real lk_safe_max_dis2_;

	bool const verbose_;
};

/// @brief creates an ExactOccludedHbondSolEnergy object according to command-line options
ExactOccludedHbondSolEnergyOP create_ExactSHOEnergy_from_cmdline(methods::EnergyMethodOptions const & options);

} // geometric_solvation
} // scoring
} // core

#endif
