// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/mhc_epitope_energy/MHCEpitopeEnergy.hh
/// @brief An energy encapsulating prediction of MHC-peptide binding, to enable deimmunization by mutagenic epitope deletion
/// @details This energy evaluates contiguous peptides in a chain and thus is not pairwise decomposable and instead leverages the ResidueaArrayAnnealableEnergy machinery. For efficiency, when a substitution is proposed, only peptides containing that position are re-evaluated.
/// The code is largely based on (via copying and modifying) NetChargeEnergy (helpers) and HBNetEnergy (only updating around substitution position)
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

#ifndef INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopeEnergy_hh
#define INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopeEnergy_hh

// Unit headers
#include <core/scoring/mhc_epitope_energy/MHCEpitopeEnergy.fwd.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopeEnergySetup.fwd.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopeEnergySetup.hh>

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

// Project headers
#include <core/types.hh>
#include <map>
#include <string>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace mhc_epitope_energy {

/// @brief A cache of epitope scores (considered and commited) at each position,
/// to enable efficient updating of only those scores affected by a substitution during annealing.
/// @details This is separated out into a class to parallel the MHCEpitopeEnergySetup helpers,
/// as each such instance will affect the score separately.
class ScoreCache : public utility::pointer::ReferenceCount {
public:
	/// @brief Construct a cache for the given number of residues in a pose
	/// (residues will be maintained and passed by the ResidueArrayAnnealableEnergy machinery)
	ScoreCache(int nres);
	ScoreCache(ScoreCache const &);
	~ScoreCache();

	ScoreCacheOP clone() const;

	/// @brief Commited scores
	utility::vector1< core::Real > scores_;
	/// @brief The range of residues affected by the currently considered substitution
	core::Size considered_start_, considered_stop_;
	/// @brief Scores affected by the currently considered substitution
	utility::vector1< core::Real > considered_scores_;
	/// @brief Scores for the native pose, in case the score is relative rather than absolute
	utility::vector1< core::Real > native_scores_;
};

class MHCEpitopeEnergy : public core::scoring::methods::WholeStructureEnergy, public core::scoring::annealing::ResidueArrayAnnealableEnergy {
public:
	typedef core::scoring::methods::WholeStructureEnergy parent1;
	typedef core::scoring::annealing::ResidueArrayAnnealableEnergy parent2;

public:

	/// @brief Options constructor.
	///
	MHCEpitopeEnergy( core::scoring::methods::EnergyMethodOptions const &options );

	/// @brief Copy constructor.
	///
	MHCEpitopeEnergy( MHCEpitopeEnergy const &src );

	/// @brief Default destructor.
	///
	virtual ~MHCEpitopeEnergy();

	/// @brief Clone: create a copy of this object, and return an owning pointer
	/// to the copy.
	core::scoring::methods::EnergyMethodOP clone() const override;

	/// @brief MHCEpitopeEnergy is context-independent and thus indicates that no context graphs need to be maintained by
	/// class Energies.
	void indicate_required_context_graphs( utility::vector1< bool > &context_graphs_required ) const override;

	/// @brief MHCEpitopeEnergy is version 1.0 right now.
	///
	core::Size version() const override;

	/// @brief Actually calculate the total energy
	/// @details Called by the scoring machinery.
	void finalize_total_energy( core::pose::Pose & pose, ScoreFunction const &, EnergyMap & totals ) const override;

	/// @brief Calculate the total energy given a vector of const owning pointers to residues.
	/// @details Called directly by the ResidueArrayAnnealingEvaluator during packer runs.
	core::Real calculate_energy(
		utility::vector1< core::conformation::ResidueCOP > const &resvect,
		utility::vector1< core::Size > const & rotamer_ids,
		core::Size const substitution_position = 0
	) const override;

	/// @brief What to do when a substitution that was considered is accepted.
	void commit_considered_substitution() override;

	/// @brief Get a summary of all loaded data.
	///
	void report() const;

	/// @brief Cache data from the pose in this EnergyMethod in anticipation of scoring.
	///
	void set_up_residuearrayannealableenergy_for_packing ( core::pose::Pose &pose, core::pack::rotamer_set::RotamerSets const &rotamersets, core::scoring::ScoreFunction const &sfxn) override;

	/// @brief Disable this energy during minimization.
	void setup_for_minimizing( pose::Pose & pose, ScoreFunction const & sfxn, kinematics::MinimizerMapBase const & minmap ) const override;

	/// @brief Re-enable this energy after minimization.
	void finalize_after_minimizing( pose::Pose & pose ) const override;

private:

	/******************
	Private functions:
	******************/

	/// @brief Get a const-access pointer to a setup helper object.
	///
	inline MHCEpitopeEnergySetupCOP setup_helper_cop( core::Size const index ) const {
		runtime_assert_string_msg( index > 0 && index <= setup_helpers_.size(), "Error in core::scoring::mhc_epitope_energy::MHCEpitopeEnergy::setup_helper_cop(): Index of setup helper is out of range." );
		return utility::pointer::dynamic_pointer_cast<MHCEpitopeEnergySetup const>(setup_helpers_[index]);
	}

	/// @brief Get the number of setup helper objects:
	///
	inline core::Size setup_helper_count() const { return setup_helpers_.size(); }

	/// @brief Given a pose, pull out the MHCEpitopeEnergySetup objects stored in SequenceConstraints in the pose and
	/// append them to the setup_helpers_ vector, returning a new vector.  This also generates a vector of masks simultaneously.
	/// @param [in] pose The pose from which the MHCEpitopeEnergySetupCOPs will be extracted.
	/// @param [out] setup_helpers The output vector of MHCEpitopeEnergySetupCOPs that is the concatenation of those stored in setup_helpers_ and those from the pose.
	/// @param [out] masks The output vector of ResidueSubsets, which will be equal in size to the helpers vector.
	/// @details The output vectors are first cleared by this operation.
	void get_helpers_from_pose(
		core::pose::Pose const &pose,
		utility::vector1< MHCEpitopeEnergySetupCOP > &setup_helpers,
		utility::vector1< core::select::residue_selector::ResidueSubset > &masks,
		utility::vector1< core::Real > &cst_weights
	) const;

	/// @brief Check for symmetry and set up the symm_multipliers_ vector so that scoring is done just on the asymmetric part and scaled up accordingly
	void setup_symmetry(core::pose::Pose const &pose) const;

	/// @brief Debugging function that will check total_ against a full_rescore total score, and returns the difference.
	/// @details This function essentially checks the integrity of the cache.  It will add significant
	/// time to the run if used during packing, as every substitution will re-calculate the score for the
	/// entire pose.  Should be used for debugging only.
	core::Real difference_btw_cached_and_full_rescore(
		utility::vector1< core::conformation::ResidueCOP > const & resvect,
		utility::vector1< MHCEpitopeEnergySetupCOP > const &setup_helpers,
		utility::vector1< core::select::residue_selector::ResidueSubset > const &masks,
		utility::vector1< core::Real > const &cst_weights
	) const;

	/// @brief Compute and return the total score for the whole sequence (resvect),
	/// according to the given energy setups (helpers) and their focused positions (masks).
	/// @details Cache indicates whether to store the score in the score cache, while native indicates whether this is the native / reference pose for relative scoring
	/// Note that method can be called either by finalize_score (setup_helpers_ and masks_ from pose)
	/// or by calculate_energy within packing simulated annealing (set_helpers_/masks_for_packing_ initialized at start of packing)
	core::Real full_rescore(
		utility::vector1< core::conformation::ResidueCOP > const & resvect,
		utility::vector1< MHCEpitopeEnergySetupCOP > const &setup_helpers,
		utility::vector1< core::select::residue_selector::ResidueSubset > const &masks,
		utility::vector1< core::Real > const &cst_weights,
		bool cache, bool native
	) const;

	/// @brief Update and return the total _considered_ score for the sequence (resvect)
	/// by computing deltas for those position potentially affected by a considered
	/// mutation (at substitution_position).
	/// @details Updates the cached per-position _considered_ scores
	core::Real update_score(
		utility::vector1< core::conformation::ResidueCOP > const & resvect,
		core::Size const substitution_position,
		utility::vector1< core::Real > const &cst_weights
	) const;

	/******************
	Private variables:
	******************/

	/// @brief Is this energy disabled (e.g. for minimization)?
	mutable bool disabled_;

	/// @brief The vector of helper objects that store all of the data for setting up this scoring function.
	/// @details Initialized on scoreterm initialization.
	utility::vector1<MHCEpitopeEnergySetupOP> setup_helpers_;

	/// @brief A cache of the MHCEpitopeEnergySetup objects to be used for packing.
	/// @details This is the list in setup_helpers_ plus all of the MHCEpitopeEnergySetup objects
	/// stored in MHCEpitopeConstraint objects in the pose.  Initialized by the ResidueArrayAnnealingEvaluator
	/// in its initialize() function.
	utility::vector1<MHCEpitopeEnergySetupCOP> setup_helpers_for_packing_;

	/// @brief A vector of masks corresponding to the MHCEpitopeEnergySetup objects in setup_helpers_for_packing_.
	/// @detalis These are generated from ResidueSelectors.
	utility::vector1< core::select::residue_selector::ResidueSubset > setup_helper_masks_for_packing_;

	/// @brief A vector of cst weights corresponding to the MHCEpitopeEnergySetup objects in setup_helpers_for_packing_.
	utility::vector1< core::Real > setup_helper_weights_for_packing_;

	/// @brif A cache for each of the MHCEpitopeEnergySetup objects in setup_helpers_for_packing_.
	utility::vector1< ScoreCacheOP > score_caches_;

	/// @brief The current commited total energy
	mutable core::Real total_;
	/// @brief The current considered total energy (when considering a substitution, before commiting it)
	mutable core::Real considered_total_;

	/// @brief Is the protein symmetric?
	mutable bool symm_;
	/// @brief If it is symmetric, only maintain scores (considered and accepted) for the asymmetric portions,
	/// and multiply them by these factors (a constant for each position).
	/// @details The value 0 indicates that the position is in a symmetric clone.
	mutable utility::vector1<core::Size> symm_multipliers_;
};

} // mhc_epitope_energy
} // scoring
} // core


#endif // INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopeEnergy_hh
