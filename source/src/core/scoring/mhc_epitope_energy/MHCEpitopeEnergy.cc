// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/mhc_epitope_energy/MHCEpitopeEnergy.cc
/// @brief An energy encapsulating prediction of MHC-peptide binding, to enable deimmunization by mutagenic epitope deletion
/// The code is largely based on (via copying and modifying) NetChargeEnergy (helpers) and HBNetEnergy (only updating around substitution position)
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

// Unit headers
#include <core/scoring/mhc_epitope_energy/MHCEpitopeEnergy.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopeEnergyCreator.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopeEnergySetup.hh>

// Package headers
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <utility/numbers.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/aa_composition_energy/SequenceConstraint.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopeConstraint.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

// Options system
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// File I/O
#include <basic/database/open.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

// Other Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/memory.hh>

#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictor.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorMatrix.hh>

namespace core {
namespace scoring {
namespace mhc_epitope_energy {

static basic::Tracer TR("core.scoring.mhc_epitope_energy.MHCEpitopeEnergy");

ScoreCache::ScoreCache(int nres)
: utility::pointer::ReferenceCount(),
	scores_(nres, 0),
	considered_start_(0),
	considered_stop_(0),
	considered_scores_(nres, 0),
	native_scores_(nres, 0)
{}

ScoreCache::ScoreCache(ScoreCache const &src)
: utility::pointer::ReferenceCount(),
	scores_(src.scores_),
	considered_start_(src.considered_start_),
	considered_stop_(src.considered_stop_),
	considered_scores_(src.considered_scores_),
	native_scores_(src.native_scores_)
{}

ScoreCache::~ScoreCache() = default;

ScoreCacheOP ScoreCache::clone() const
{
	return ScoreCacheOP( utility::pointer::make_shared<ScoreCache>(*this));
}


/// @brief This must return a fresh instance of the MHCEpitopeEnergy class, never an instance already in use.
///
core::scoring::methods::EnergyMethodOP
MHCEpitopeEnergyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const &options ) const
{
	return utility::pointer::make_shared< MHCEpitopeEnergy >( options );
}

/// @brief Defines the score types that this energy method calculates.
///
ScoreTypes
MHCEpitopeEnergyCreator::score_types_for_method() const
{
	ScoreTypes sts;
	sts.push_back( mhc_epitope );
	return sts;
}

/// @brief Options constructor.
///
MHCEpitopeEnergy::MHCEpitopeEnergy ( core::scoring::methods::EnergyMethodOptions const &options ) :
	parent1( utility::pointer::make_shared< MHCEpitopeEnergyCreator >() ),
	parent2( ),
	disabled_(false),
	setup_helpers_(),
	setup_helpers_for_packing_(),
	setup_helper_masks_for_packing_(),
	setup_helper_weights_for_packing_(),
	score_caches_()
{
	//The following reads from disk the first time only, and caches the data in memory:
	setup_helpers_ = core::scoring::ScoringManager::get_instance()->get_cloned_mhc_epitope_setup_helpers( options );
	if ( TR.Debug.visible() ) report();
}

/// @brief Copy constructor.
///
MHCEpitopeEnergy::MHCEpitopeEnergy( MHCEpitopeEnergy const &src ) :
	parent1( utility::pointer::make_shared< MHCEpitopeEnergyCreator >() ),
	parent2( src ),
	disabled_( src.disabled_ ),
	setup_helpers_(), //CLONE the helper data below; don't copy them.
	setup_helpers_for_packing_(), //CLONE these below, too -- don't copy them.
	setup_helper_masks_for_packing_(  src.setup_helper_masks_for_packing_ ),
	setup_helper_weights_for_packing_( src.setup_helper_weights_for_packing_ ),
	score_caches_() // below
{
	for ( core::Size i=1, imax=src.setup_helpers_.size(); i<=imax; ++i ) {
		setup_helpers_.push_back( src.setup_helpers_[i]->clone() );
	}
	for ( core::Size i=1, imax=src.setup_helpers_for_packing_.size(); i<=imax; ++i ) {
		setup_helpers_for_packing_.push_back( src.setup_helpers_for_packing_[i]->clone() );
	}
	for ( core::Size i=1, imax=src.score_caches_.size(); i<=imax; ++i ) {
		score_caches_.push_back( src.score_caches_[i]->clone() );
	}
}

/// @brief Default destructor.
///
MHCEpitopeEnergy::~MHCEpitopeEnergy() = default;

/// @brief Clone: create a copy of this object, and return an owning pointer
/// to the copy.
core::scoring::methods::EnergyMethodOP MHCEpitopeEnergy::clone() const {
	return utility::pointer::make_shared< MHCEpitopeEnergy >(*this);
}

/// @brief MHCEpitopeEnergy is context-independent and thus indicates that no context graphs need to be maintained by
/// class Energies.
void MHCEpitopeEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const
{
	//Do nothing.
	return;
}

/// @brief MHCEpitopeEnergy is version 1.0 right now.
///
core::Size MHCEpitopeEnergy::version() const
{
	return 1; // Initial versioning
}

/// @brief Actually calculate the total energy
/// @details Called by the scoring machinery.
void MHCEpitopeEnergy::finalize_total_energy( core::pose::Pose & pose, ScoreFunction const &, EnergyMap & totals ) const
{
	if ( disabled_ ) return; //Do nothing when this energy is disabled.

	TR.Debug << "finalizing score" << std::endl;

	// Need to make sure we have the symmetry information stored, in order to multiply accordingly
	if ( symm_multipliers_.size() == 0 ) setup_symmetry(pose);

	//Number of residues:
	core::Size const nres( pose.size() );

	//Vector of residue owning pointers:
	utility::vector1< core::conformation::ResidueCOP > reslist;
	reslist.reserve(nres);

	//Populate the vector with const owning pointers to the residues:
	for ( core::Size ir=1; ir<=nres; ++ir ) {
		reslist.push_back( pose.residue(ir).get_self_ptr() );
	}

	// For diagnostic purposes, make a note of where epitopes are broken due to
	// defined chain breaks and residues that aren't the 20 AAs
	for ( core::Size ir(1); ir<=nres; ++ir ) {
		// If it's a symmetric clone, don't directly add to the score -- its contribution comes by way of a multiplicative factor on the original in the asymmetric unit
		if ( symm_multipliers_[ir] == 0 ) continue;
		// If the chain ID number changes between ir-1 and ir (and we aren't on the first residue), report it
		if ( ir>1 && reslist[ir]->chain() != reslist[ir-1]->chain() ) {
			TR.Debug << "chain break betwen " << ir-1 << " and " << ir << std::endl;
		}
		// If the restype is a ligand or non-canonical AA, report it.
		if ( reslist[ir]->type().aa() > core::chemical::AA::num_canonical_aas ) {
			TR.Debug << "Ligand or non-canonical AA " << reslist[ir]->type().aa() << " at " << ir << std::endl;
		}
	}

	//Get the MHCEpitopeEnergySetup objects from the pose and append them to the setup_helpers_ list, making a new setup_helpers list:
	utility::vector1< MHCEpitopeEnergySetupCOP > setup_helpers;
	utility::vector1< core::select::residue_selector::ResidueSubset > masks;
	utility::vector1< core::Real > cst_weights;
	get_helpers_from_pose( pose, setup_helpers, masks, cst_weights ); //Pulls MHCEpitopeEnergySetupCOPs from pose; generates masks from ResidueSelectors simultaneously.
	runtime_assert( masks.size() == setup_helpers.size() ); //Should be guaranteed to be true.
	runtime_assert( cst_weights.size() == setup_helpers.size() ); //Should be guaranteed to be true.

	totals[ mhc_epitope ] += full_rescore( reslist, setup_helpers, masks, cst_weights, false, true ); //Using the vector of residue owning pointers, calculate the energy (unweighted) and set the mhc_epitope to this value.
}

/// @brief Calculate the total energy given a vector of const owning pointers to residues.
/// @details Called directly by the ResidueArrayAnnealingEvaluator during packer runs.  Requires
/// that set_up_residuearrayannealablenergy_for_packing() be called first.
core::Real
MHCEpitopeEnergy::calculate_energy(
	utility::vector1< core::conformation::ResidueCOP > const &reslist,
	utility::vector1< core::Size > const &,
	core::Size const substitution_position /* = 0 */
) const {
	if ( disabled_ ) return 0.0; //Do nothing when disabled.
	if ( substitution_position == 0 ) {
		// No substitution => do a full rescore
		return full_rescore(reslist, setup_helpers_for_packing_, setup_helper_masks_for_packing_, setup_helper_weights_for_packing_, true, false);
	} else {
		// Consider impact of substitutions on scores of peptides overlapping it
		return update_score(reslist, substitution_position, setup_helper_weights_for_packing_);
	}
}

/// @brief What to do when a substitution that was considered is accepted.
void
MHCEpitopeEnergy::commit_considered_substitution() {
	// Commit the considered total
	total_ = considered_total_;

	// Commit the considered values in the caches
	for ( core::Size h=1, hmax=score_caches_.size(); h<=hmax; ++h ) {
		ScoreCacheOP cache(score_caches_[h]);
		for ( core::Size i = cache->considered_start_; i <= cache->considered_stop_; i++ ) {
			cache->scores_[i] = cache->considered_scores_[i];
		}
	}
}

/// @brief Get a summary of all loaded data.
///
void MHCEpitopeEnergy::report() const {
	if ( !TR.Debug.visible() ) return; //Do nothing if I don't have a tracer.

	TR.Debug << std::endl << "Summary of data loaded by MHCEpitopeEnergy object:" << std::endl;

	for ( core::Size i=1, imax=setup_helper_count(); i<=imax; ++i ) {
		TR.Debug << "MHCEpitopeEnergySetup #" << i << ":" << std::endl;
		TR.Debug << setup_helper_cop(i)->report();
	}

	TR.Debug << std::endl;

	TR.Debug.flush();
}

/// @brief Cache data from the pose in this EnergyMethod in anticipation of scoring.
///
void
MHCEpitopeEnergy::set_up_residuearrayannealableenergy_for_packing (
	core::pose::Pose & pose,
	core::pack::rotamer_set::RotamerSets const & /*rotamer_sets*/,
	core::scoring::ScoreFunction const & /*sfxn*/
) {
	disabled_ = false; //Enabled for packing.

	core::Size const nres( pose.size() );

	// Set up helpers
	get_helpers_from_pose( pose, setup_helpers_for_packing_, setup_helper_masks_for_packing_, setup_helper_weights_for_packing_ );

	runtime_assert( setup_helper_masks_for_packing_.size() == setup_helpers_for_packing_.size() ); //Should be guaranteed to be true.
	runtime_assert( setup_helper_weights_for_packing_.size() == setup_helpers_for_packing_.size() ); //Should be guaranteed to be true.

	// Summaries of the various helpers
	if ( TR.Debug.visible() ) {
		for ( core::Size i=1, imax=setup_helper_masks_for_packing_.size(); i<=imax; ++i ) {
			TR.Debug << "Setup helper " << i << ":" << std::endl;
			TR.Debug << "Mask: ";
			for ( core::Size j=1, jmax=setup_helper_masks_for_packing_[i].size(); j<=jmax; ++j ) {
				TR.Debug << (setup_helper_masks_for_packing_[i][j] ? "1" : "0");
			}
			TR.Debug << std::endl;
			TR.Debug << setup_helpers_for_packing_[i]->report();
		}
		TR.Debug.flush();
	}

	// Set up score caches
	score_caches_.clear();
	for ( core::Size h=1, hmax=setup_helpers_for_packing_.size(); h<=hmax; ++h ) {
		ScoreCacheOP cache(utility::pointer::make_shared<ScoreCache>(nres));
		score_caches_.push_back(cache);
	}

	// Need to make sure we have the symmetry information stored, in order to multiply accordingly
	if ( symm_multipliers_.size() == 0 ) setup_symmetry(pose);

	// Finally populate the caches by full rescore
	// Create a vector of size nres
	utility::vector1< core::conformation::ResidueCOP > reslist;
	reslist.reserve(nres);
	// Populate with pointers to each residue
	for ( core::Size ir=1; ir<=nres; ++ir ) {
		reslist.push_back( pose.residue(ir).get_self_ptr() );
	}
	// Rescore over entire reslist, caching both native and xformed scores. Masks are passed through.
	full_rescore(reslist, setup_helpers_for_packing_, setup_helper_masks_for_packing_, setup_helper_weights_for_packing_, true, true);

	// Report the native scores
	if ( TR.Debug.visible() ) {
		for ( core::Size h=1, hmax=setup_helpers_for_packing_.size(); h<=hmax; ++h ) {
			TR.Debug << "native scores #" << h << ": ";
			for ( core::Size ir=1; ir<=nres; ++ir ) {
				TR.Debug << score_caches_[h]->scores_[ir] << " ";
			}
			TR.Debug << std::endl;
		}
	}
}

/// @brief Disable this energy during minimization.
void
MHCEpitopeEnergy::setup_for_minimizing( pose::Pose & /*pose*/, ScoreFunction const & /*sfxn*/, kinematics::MinimizerMapBase const & /*minmap*/ ) const {
	TR.Debug << "Disabling MHCEpitopeEnergy during minimization." << std::endl;
	disabled_ = true;
}

/// @brief Re-enable this energy after minimization.
void
MHCEpitopeEnergy::finalize_after_minimizing( pose::Pose & /*pose*/ ) const {
	TR.Debug << "Re-enabling MHCEpitopeEnergy following minimization." << std::endl;
	disabled_ = false;
}

/// @brief Given a pose, pull out the MHCEpitopeEnergySetup objects stored in SequenceConstraints in the pose and
/// append them to the setup_helpers_ vector, returning a new vector.  This also generates a vector of masks simultaneously.
/// @param [in] pose The pose from which the MHCEpitopeEnergySetupCOPs will be extracted.
/// @param [out] setup_helpers The output vector of MHCEpitopeEnergySetupCOPs that is the concatenation of those stored in setup_helpers_ and those from the pose.
/// @param [out] masks The output vector of ResidueSubsets, which will be equal in size to the helpers vector.
/// @details The output vectors are first cleared by this operation.
void
MHCEpitopeEnergy::get_helpers_from_pose(
	core::pose::Pose const &pose,
	utility::vector1< MHCEpitopeEnergySetupCOP > &setup_helpers,
	utility::vector1< core::select::residue_selector::ResidueSubset > &masks,
	utility::vector1< core::Real > &cst_weights
) const {
	setup_helpers.clear();
	masks.clear();
	cst_weights.clear();
	if ( setup_helpers_.size() > 0 ) {
		setup_helpers = setup_helpers_; //Copy the setup_helpers_ list.
		// TODO: drop any without a specified predictor? (as below with cur_cst)
		masks.resize( setup_helpers_.size(), core::select::residue_selector::ResidueSubset( pose.size(), true ) ); //All of the helpers in the setup_helpers_ list should be applied globally.
		cst_weights.resize( setup_helpers_.size(), 1.0 );
	}

	core::Size const n_sequence_constraints( pose.constraint_set()->n_sequence_constraints() );
	if ( n_sequence_constraints > 0 ) {
		// Loop over each set of sequence constraints in the pose.
		for ( core::Size i=1; i<=n_sequence_constraints; ++i ) {
			MHCEpitopeConstraintCOP cur_cst( utility::pointer::dynamic_pointer_cast<MHCEpitopeConstraint const>( pose.constraint_set()->sequence_constraint(i) ) );
			if ( !cur_cst ) continue; //Continue if this isn't an MHCEpitopeConstraint.
			if ( cur_cst->mhc_epitope_energy_setup()->is_default() ) {
				TR << "got a constraint with no specified predictor; ignoring" << std::endl;
				continue;
			}
			setup_helpers.push_back( cur_cst->mhc_epitope_energy_setup() ); //Append the MHCEpitopeEnergySetup object stored in the current sequence constraint to the list to be used.
			core::select::residue_selector::ResidueSelectorCOP selector( cur_cst->selector() ); //Get the ResidueSelector in the current sequence constraint object, if there is one.  (May be NULL).
			if ( selector ) { //If we have a ResidueSelector, generate a mask from the pose and store it in the masks list.
				masks.push_back( selector->apply( pose ) );
			} else { //If not, add an all-true mask
				masks.push_back( core::select::residue_selector::ResidueSubset( pose.size(), true ) );
			}

			//Get the cst weight for the current contraint.
			cst_weights.push_back( cur_cst->get_cst_weight() );
		}
	}

	runtime_assert( setup_helpers.size() == masks.size() ); //Should be guaranteed true.
	runtime_assert( setup_helpers.size() == cst_weights.size() ); //Should be guaranteed true.

	TR.Debug << "End of get_helpers_from_pose_" << std::endl;
	TR.Debug << "Number of helpers: " << setup_helpers.size() << std::endl;
	TR.Debug << "Masks: " << masks << std::endl;
	TR.Debug << "Weights: " << cst_weights << std::endl;

	return;
}

// ----------------------------------------
// Details of MHC epitope scoring

/// @brief Helper for testing -- returns the AA sequence, with X for noncanonical and space between chains
std::string get_seq(utility::vector1< core::conformation::ResidueCOP > const & reslist,
	utility::vector1<core::Size> & symm_multipliers)
{
	std::stringstream output("");

	core::Size const stop = reslist.size();
	bool broken = true;
	for ( core::Size start=1; start <= stop; start++ ) {
		core::conformation::ResidueCOP res = reslist[start];
		if ( symm_multipliers[start] == 0 || (start > 1 && res->chain() != reslist[start-1]->chain()) ) {
			if ( !broken ) {
				broken = true;
				output << " ";
			}
		} else {
			broken = false;
			if ( res->type().aa() > core::chemical::AA::num_canonical_aas ) {
				output << "X";
			} else {
				output << res->name1();
			}
		}
	}

	return output.str();
}

void MHCEpitopeEnergy::setup_symmetry(core::pose::Pose const &pose) const
{
	symm_ = core::pose::symmetry::is_symmetric( pose );
	if ( !symm_ ) {
		// Single subunit; every residue gets a multiplier of 1
		for ( core::Size ir=1; ir<=pose.size(); ++ir ) symm_multipliers_.push_back(1);
	} else {
		// Get the symmetry info
		core::conformation::symmetry::SymmetricConformationCOP symmconf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( pose.conformation_ptr() ) );
		core::conformation::symmetry::SymmetryInfoCOP symm_info = symmconf->Symmetry_Info();
		for ( core::Size ir=1; ir<=pose.size(); ++ir ) {
			// If the residue is independent and not virtual, set the multiplier to 1 + the number of symmetric clones (i.e. the number of copies of that residue).
			// If the residue is not independent, we don't count it, as it will be taken care of by the multiplier on the independent copy. So set to 0.
			// The result of this should be that every residue in the ASU should be have its multiplier set to the number of copies it has in the symmetric pose.
			// The residues in the symmetric copies should be set to 0 so MHCEpitopeEnergy doesn't look at them during scoring.
			core::Size multiplier = (symm_info->independent_residues()[ir] && !symm_info->is_virtual(ir)) ? (1+symm_info->bb_clones(ir).size()) : 0;
			symm_multipliers_.push_back(multiplier);
		}
		if ( TR.Debug.visible() ) {
			// Report the multipliers
			TR.Debug << "symmetry multipliers:";
			for ( core::Size ir=1; ir<=pose.size(); ++ir ) TR.Debug << " " << symm_multipliers_[ir];
			TR.Debug << std::endl;
		}
	}
}

core::Real MHCEpitopeEnergy::difference_btw_cached_and_full_rescore(
	utility::vector1< core::conformation::ResidueCOP > const & resvect,
	utility::vector1< MHCEpitopeEnergySetupCOP > const &setup_helpers,
	utility::vector1< core::select::residue_selector::ResidueSubset > const &masks,
	utility::vector1< core::Real > const &cst_weights
) const {
	//Store the current total_ in tmp
	core::Real tmp = total_;

	//Do a full re-score, passing through the relevant parameters.
	//The last two variables (cache and native) should always be false.
	core::Real full = full_rescore(resvect, setup_helpers, masks, cst_weights, false, false);

	//full_rescore will update total_ to be the full_rescore value.
	//We don't want to lose gradual accumulation of differences between the cached score and the full_rescore, so reset total_ to its value before running full_rescore.
	total_ = tmp;

	//Return the different between considered_total_ and full.  Should be 0 if the cache is working.
	core::Real diff = std::abs(considered_total_ - full);
	TR.Debug << "Difference between cached considered_total_ and full rescore is " << diff << std::endl;

	if ( diff > 0.1 ) {
		TR.Error << "The difference between the currently cached mhc_energy score and the fully rescored pose is " << diff << std::endl;
		TR.Error << "This shouldn't happen, and likely indicates a bug with how mhc_energy scores are being cached with your pose!" << std::endl;
	}

	return(diff);
}

core::Real MHCEpitopeEnergy::full_rescore(
	utility::vector1< core::conformation::ResidueCOP > const & reslist,
	utility::vector1< MHCEpitopeEnergySetupCOP > const &setup_helpers,
	utility::vector1< core::select::residue_selector::ResidueSubset > const &masks,
	utility::vector1< core::Real > const &cst_weights,
	bool cache, bool native
) const {
	if ( setup_helpers.size() == 0 ) return 0;

	// Running total score over all peptides, all helpers
	total_ = 0;
	// For testing, the raw total score
	core::Real raw_total = 0;

	// TODO: could flip the order of loops, so only extract peptide once
	// but since predictors don't necessarily use peptides of the same length, will get more complicated

	// Loop over helpers
	core::Size ihelpermax=setup_helpers.size();
	for ( core::Size ihelper=1; ihelper<=ihelpermax; ++ihelper ) {
		MHCEpitopeEnergySetupCOP helper( setup_helpers[ihelper] );
		if ( helper->is_default() ) continue; // default contributes 0 to score, so just bypass

		// Loop over peptide starting positions
		// Generate a "peptide" of A's of the appropriate length
		core::Size len = helper->get_peptide_length();
		std::string peptide(len, 'A');
		core::Size const stop = reslist.size()-len+1; // Set stop to the last residue where an epitope could start
		// Loop over all peptide starting positions
		for ( core::Size start=1; start <= stop; start++ ) {
			//BJY: It is tempting to throw this code into its own function and call it from full_rescore and update_score, to avoid duplication.
			//BJY: The treatment of the masks and cache are a bit different between the two.  TODO.
			//CBK: yeah, that's why it's duplicated -- the amount of duplication is probably less than the amount of ugliness to handle that
			if ( symm_multipliers_[start] == 0 ) continue; // not in asymmetric unit, so doesn't contribute
			if ( !masks[ihelper][start] ) continue; // doesn't contribute to this because it's masked
			bool chain_break=false, noncanon=false;
			core::Size i=0; // need outside of loop, to jump if not contiguous
			// Check for chainbreaks, noncannonicals, and symmetry breaks.
			for ( ; i < len; i++ ) {
				if ( symm_multipliers_[start+i] == 0 ) {
					// Went into symmetric clone; treat like chain break
					chain_break = true;
					break;
				}
				core::conformation::ResidueCOP res = reslist[start+i];
				// If the chain id is different from the previous residue's chain id, flag it as a chainbreak.
				if ( i>0 && res->chain() != reslist[start+i-1]->chain() ) {
					chain_break = true;
					break;
				}
				// If the restype is not a canonical AA, flag it as noncanonical.
				// This will catch ligands too.
				if ( res->type().aa() > core::chemical::AA::num_canonical_aas ) {
					noncanon = true;
					break;
				}
				// If the residue is not flagged, change peptide to include it at the appropriate position.
				peptide[i] = res->name1();
			}
			if ( chain_break || noncanon ) {
				if ( chain_break ) i--; // will pick up with this AA (just different chain) vs. will skip over it (noncanonical)
				if ( cache ) {
					// Zero out the contributions from peptides spanning the position
					for ( core::Size j=0; j<=i && start+j<=stop; j++ ) {
						score_caches_[ihelper]->scores_[start+j] = 0;
					}
				}
				start += i; // jump to this point, since everything else before it will also be broken (note that i will also be incremented as normal in for-loop)
			} else {
				// Deal with raw and xformed scores in separate steps so can cache and xform appropriately
				// Start with raw score from epitope predictor
				core::Real score = helper->raw_score(peptide);
				raw_total += score;
				if ( native ) {
					if ( cache ) {
						// Cache the native score
						score_caches_[ihelper]->native_scores_[start] = score;
					}
					// Xform it (it is the native score)
					score = helper->xform(score, score);
				} else {
					// Xform it (using the cached native score)
					score = helper->xform(score, score_caches_[ihelper]->native_scores_[start]);
				}
				if ( cache ) {
					// Cache the xformed score
					score_caches_[ihelper]->scores_[start] = score;
				}
				//total_ should be the score * any symmetry multipliers * any cst weights (if applicble)
				total_ += score * symm_multipliers_[start] * cst_weights[ihelper];
			}
		}
	}

	TR.Debug << "full total: " << total_ << "; raw_total: " << raw_total << " " << get_seq(reslist, symm_multipliers_)<< std::endl;

	return total_;
}

core::Real MHCEpitopeEnergy::update_score(
	utility::vector1< core::conformation::ResidueCOP > const & reslist,
	core::Size const subst_resid,
	utility::vector1< core::Real > const &cst_weights
) const {
	// Same basic thing as the total, but only around subst_resid, and caching considered scores and summing considered total

	// Start from current total and update it around subst_resid
	considered_total_ = total_;

	// TODO: see note in full_rescore about loop order

	// Loop over helpers
	core::Size ihelpermax=setup_helpers_for_packing_.size();
	for ( core::Size ihelper=1; ihelper<=ihelpermax; ++ihelper ) {
		MHCEpitopeEnergySetupCOP helper( setup_helpers_for_packing_[ihelper] );
		if ( helper->is_default() ) continue; // default contributes 0 to score, so just bypass
		ScoreCacheOP cache(score_caches_[ihelper]);

		// Loop over peptide starting positions
		core::Size len = helper->get_peptide_length();
		std::string peptide(len, 'A');

		// Starting peptide: upstream by len, as long as in the same chain as the substitution
		if ( subst_resid <= len ) cache->considered_start_ = 1;
		else cache->considered_start_ = subst_resid-len+1;
		core::Size chain(reslist[subst_resid]->chain());
		while ( reslist[cache->considered_start_]->chain() != chain ) cache->considered_start_++;

		// Ending peptide: the substitution position itself, as long as there's len worth of protein after it
		cache->considered_stop_ = reslist.size()-len+1;
		if ( subst_resid < cache->considered_stop_ ) cache->considered_stop_ = subst_resid;

		for ( core::Size start = cache->considered_start_; start <= cache->considered_stop_; start++ ) {
			if ( symm_multipliers_[start] == 0 ) continue; // not in asymmetric unit
			if ( !setup_helper_masks_for_packing_[ihelper][start] ) continue; // doesn't contribute to this one
			bool chain_break=false, noncanon=false;
			core::Size i=0; // need outside of loop, to jump if not contiguous
			for ( ; i < len; i++ ) {
				if ( symm_multipliers_[start+i] == 0 ) {
					// Went into symmetric clone; treat like chain break
					chain_break = true;
					break;
				}
				core::conformation::ResidueCOP res = reslist[start+i];
				if ( i>0 && res->chain() != reslist[start+i-1]->chain() ) {
					chain_break = true;
					break;
				}
				if ( res->type().aa() > core::chemical::AA::num_canonical_aas ) {
					noncanon = true;
					break;
				}
				peptide[i] = res->name1();
			}
			if ( chain_break || noncanon ) {
				if ( chain_break ) i--; // will pick up with this AA (just different chain) vs. will skip over it (noncanonical)
				for ( core::Size j=0; j<=i && start+j<=cache->considered_stop_; j++ ) {
					// Zero out the contributions from peptides spanning the position
					cache->considered_scores_[start+j] = 0;
					considered_total_ -= symm_multipliers_[start+j] * cst_weights[ihelper] * cache->scores_[start+j];
				}
				start += i; // jump to this point, since everything else before it will also be broken (note that i will also be incremented as normal in for-loop)
			} else {
				// Compute the epitope (xforming as needed)
				cache->considered_scores_[start] = helper->score(peptide, score_caches_[ihelper]->native_scores_[start]);
				// Update the considered total by subtracting out the contribution from what was there
				// and adding in the contribution from what is being considered to be there now
				considered_total_ += symm_multipliers_[start] * cst_weights[ihelper] * (cache->considered_scores_[start] - cache->scores_[start]);
			}
		}
	}

	// Periodically give the sequence and its score, and if in debug mode, check the difference between cached score and full_rescore.
	static int reps = 0;
	if ( reps < 100 ) {
		++reps;
	} else {
		reps = 0;
		TR.Debug << "considered total " << considered_total_ << " " << get_seq(reslist, symm_multipliers_) << std::endl;

		// This does a full re-score after each 100 updates to check the accuracy of the score cache.
		// Enabled for debugging, though it should be turned off if run times are an issue.
		// Doesn't execute in release mode.
		if ( true ) {
			// For debugging, calculate the difference between cached and full rescore.
			// Allowing for some slop in the match. If there's too much drift, we should full rescore periodically....
			debug_assert( 0.1 > difference_btw_cached_and_full_rescore(reslist, setup_helpers_for_packing_, setup_helper_masks_for_packing_, setup_helper_weights_for_packing_) );
		}
	}

	// This significantly adds to the length of a run, as it will do a full re-score with each update.
	// If cache accuracy issues are becoming a problem, you can check the cache every update to help debug.
	// Turn the following if statement to true, and turn the above if statement to false (to avoid checking the score twice).
	if ( false ) {
		// For debugging, calculate the difference between cached and full rescore.
		// Allowing for some slop in the match. If there's too much drift, we should full rescore periodically....
		debug_assert( 0.1 > difference_btw_cached_and_full_rescore(reslist, setup_helpers_for_packing_, setup_helper_masks_for_packing_, setup_helper_weights_for_packing_) );
	}

	return considered_total_;
}

} // mhc_epitope_energy
} // scoring
} // core
