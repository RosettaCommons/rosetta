// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/voids_penalty_energy/VoidsPenaltyEnergy.hh
/// @brief Headers for an EnergyMethod intended for packing, which penalizes solutions in which the total volume to fill differs greatly
/// from the total volume of the current set of rotamers.
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

#ifndef INCLUDED_core_pack_guidance_scoreterms_voids_penalty_energy_VoidsPenaltyEnergy_hh
#define INCLUDED_core_pack_guidance_scoreterms_voids_penalty_energy_VoidsPenaltyEnergy_hh

// Unit headers
#include <core/pack/guidance_scoreterms/voids_penalty_energy/VoidsPenaltyEnergy.fwd.hh>
#include <core/pack/guidance_scoreterms/voids_penalty_energy/VoidsPenaltyVoxelGrid.fwd.hh>

// Package headers
#include <core/scoring/annealing/ResidueArrayAnnealableEnergy.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueProperty.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

// Project headers
#include <core/types.hh>
#include <map>
#include <string>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace voids_penalty_energy {

/// @brief VoidsPenaltyEnergy ("voids_penalty" score term), an EnergyMethod intended for packing, which penalizes solutions in which the total volume to fill differs greatly
/// from the total volume of the current set of rotamers.
class VoidsPenaltyEnergy : public core::scoring::methods::WholeStructureEnergy, public core::scoring::annealing::ResidueArrayAnnealableEnergy {
public:
	typedef core::scoring::methods::WholeStructureEnergy parent1;
	typedef core::scoring::annealing::ResidueArrayAnnealableEnergy parent2;

public:

	/// @brief Options constructor.
	///
	VoidsPenaltyEnergy( core::scoring::methods::EnergyMethodOptions const &options );

	/// @brief Copy constructor.
	///
	VoidsPenaltyEnergy( VoidsPenaltyEnergy const &src );

	/// @brief Default destructor.
	///
	virtual ~VoidsPenaltyEnergy();

	/// @brief Clone: create a copy of this object, and return an owning pointer
	/// to the copy.
	core::scoring::methods::EnergyMethodOP clone() const override;

	/// @brief VoidsPenaltyEnergy is context-independent and thus indicates that no context graphs need to be maintained by
	/// class Energies.
	void indicate_required_context_graphs( utility::vector1< bool > &context_graphs_required ) const override;

	/// @brief VoidsPenaltyEnergy is version 1.0 right now.
	///
	core::Size version() const override;

	/// @brief Actually calculate the total energy
	/// @details Called by the scoring machinery.
	/// @note VoidsPenaltyEnergy::finalize_total_energy() can return a slightly different energy than was computed
	/// during packing.  This is because reachable volume cannot be computed (since we don't have a rotamer set), so
	/// total buried volume is used in the calculation.  This will also be a whole-pose calculation, and won't just
	/// focus on the designable region (since there's no "designable region" when this function is called).
	void finalize_total_energy( core::pose::Pose & pose, core::scoring::ScoreFunction const &, core::scoring::EnergyMap & totals ) const override;

	/// @brief Calculate the total energy given a vector of const owning pointers to residues.
	/// @details Called directly by the ResidueArrayAnnealingEvaluator during packer runs.
	core::Real calculate_energy( utility::vector1< core::conformation::ResidueCOP > const &resvect, core::Size const substitution_position = 0 ) const override;

	/// @brief Get a summary of all loaded data.
	///
	void report() const;

	/// @brief Cache data from the pose in this EnergyMethod in anticipation of scoring.
	///
	void set_up_residuearrayannealableenergy_for_packing ( core::pose::Pose &pose, core::pack::rotamer_set::RotamerSets const &rotamersets, core::scoring::ScoreFunction const &sfxn) override;

	/// @brief Called at the beginning of atom tree minimization, this method
	/// allows the derived class the opportunity to initialize pertinent data
	/// that will be used during minimization.  During minimzation, the chemical
	/// structure of the pose is constant, so assumptions on the number of atoms
	/// per residue and their identities are safe so long as the pose's Energies
	/// object's "use_nblist()" method returns true.
	/// @details This just disables this score term during minimization, in the case of the VoidsPenaltyEnergy.
	void setup_for_minimizing( pose::Pose & pose, core::scoring::ScoreFunction const & sfxn, kinematics::MinimizerMapBase const &minmap ) const override;

	/// @brief Called after minimization.
	/// @details Re-enables the score term after minimization.
	void finalize_after_minimizing( pose::Pose & pose ) const override;

	/// @brief Set whether this term is disabled except during packing.
	inline void set_disabled_except_during_packing( bool const setting ) { disabled_except_during_packing_ = setting; }

	/// @brief Get whether this term is disabled except during packing.
	inline bool disabled_except_during_packing( ) const { return disabled_except_during_packing_; }

private:

	/******************
	Private functions:
	******************/

	/// @brief given a voxel grid object, set its parameters using stored values (originally
	/// from user settings or the options system).
	void configure_voxel_grid( VoidsPenaltyVoxelGrid &voxel_grid ) const;

	/******************
	Private variables:
	******************/

	/// @brief Cone dot product cutoff, used by the voxel grid calculator.
	/// @details The cutoff value for the dot product of a cone vector and a cone base-test point vector below which we declare the test point not to be
	/// within the cone.  Effectively, this is the cone width.  Lower values make broader cones.  Default 0.1.  Can range from 1.0 (infinitely thin
	/// cone) to -1.0 (full spherical volume), with 0.0 represeting all points on one side of the plane perpendicular to the cone vector.
	core::Real cone_dotproduct_cutoff_;

	/// @brief Cone distance cutoff, used by the voxel grid calculator.
	/// @details The cutoff value for the distance from the cone base at which we are considered no longer to be within the cone.  Defaults
	/// to 8.0 Angstroms.
	core::Real cone_distance_cutoff_;

	/// @brief The minimum number of cones in which a voxel must lie in order for that voxel to be considered "buried".
	/// @details Defaults to 6 cones.
	core::Size containing_cones_cutoff_;

	/// @brief Voxel grid voxel size, used by the voxel grid calculator.
	/// @details In Angstroms.  Length of one side of a cube.
	core::Real voxel_size_;

	/// @brief Voxel grid lateral padding, used by the voxel grid calculator.
	/// @details This is added to the bounding box of the pose when setting up the voxel grid.
	core::Real voxel_grid_padding_;

	/// @brief Is the voids_penalty energy term disabled except during packing?
	/// @details If true (the default), the term is only evaluated during packing.  If false, it is evaluated during
	/// packing or scoring (but not minimizing).
	bool disabled_except_during_packing_;

	/***********************
	Cached data for packing:
	************************/

	/// @brief The volumes of rotamers, cached from a precomputation for packing.
	mutable utility::vector1< std::map < core::conformation::ResidueCOP, core::Real > > rotamer_volumes_;

	/// @brief The total volume of the space to fill during packing, cached from a precomputation.
	mutable core::Real volume_to_fill_;

	/// @brief Are we currently minimizing?
	/// @details The VoidsPenaltyEnergy is not evaluated during minimization.
	mutable bool minimizing_;

};

} // voids_penalty_energy
} // guidance_scoreterms
} // pack
} // core


#endif // INCLUDED_core_pack_EtableEnergy_HH
