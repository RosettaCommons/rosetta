// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file   FlexPepDockingPoseMetrics.hh
///
/// @brief metrics calculations specific for FlexPepDock (at least for now)
/// @date March 29th, 2009
/// @author Barak Raveh / Nir London

#ifndef INCLUDED_protocols_flexpep_docking_FlexPepDockingPoseMetrics_hh
#define INCLUDED_protocols_flexpep_docking_FlexPepDockingPoseMetrics_hh

#include <protocols/flexpep_docking/FlexPepDockingFlags.hh>
#include <protocols/flexpep_docking/FlexPepDockingPoseMetrics.fwd.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.fwd.hh>

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace flexpep_docking {

class FlexPepDockingPoseMetrics
{

public:

	typedef bool (*t_predicate_func)(core::pose::Pose const &, core::pose::Pose const &, core::Size, core::Size);


	FlexPepDockingPoseMetrics(FlexPepDockingFlagsCOP flags)
	: flags_(flags)
	{} // TODO: move to .cc

	core::Real calc_frac_native_contacts ( core::pose::Pose const& native,
		core::pose::Pose const& final,
		core::Real threashold ) const;

	/////////////////////////////////////////////////////////////////////////////
	///
	/// @brief calculate fractions of atoms that are at the same location as native
	///
	//  @details
	//  Calculate fraction of atoms that are less than k Angstroms between native and
	//  final pose, using only residues marked in res_subset, and only atoms qualifying to
	//  predicate function (see core/scoring/rms_utils.hh for definitions of predicates)
	//
	//  @param
	//  pose1 - the first structure to be assessed
	//  pose2 - the second structure to be assessed
	//  predicate - a predicate on atoms included in the assesment (e.g., CA atoms only)
	//  k - threshold in angstroms of atoms to be included in the assesment
	//  ngood [output] - absolute number of atoms that are less than k Angstroms
	//
	// @return
	// fraction of atoms that are less than k Angstroms between poses
	////////////////////////////////////////////////////////////////////////////

	core::Real calc_frac_atoms_kA_to_native(
		core::pose::Pose const& pose1, core::pose::Pose const& pose2,
		ObjexxFCL::FArray1D_bool const & res_subset,
		t_predicate_func predicate,
		double k,
		core::Size& ngood
	) const;


	// check all sequential Kmers in the peptide, and output the best RMS Kmer
	core::Real best_Kmer_rms(
		core::pose::Pose const& pose1,
		core::pose::Pose const& pose2,
		t_predicate_func predicate,
		core::Size k
	) const;

	/////////////////////////////////////////////////////////////////////////////
	///
	//  @details
	//  Calculate the RMSD in phi.psi angle over peptide between pose1 and pose2
	//
	//  @param
	//  pose1 - the first structure to be assessed
	//  pose2 - the second structure to be assessed
	//
	// @return
	// phi/psi torsion-RMSD between peptide backbones
	////////////////////////////////////////////////////////////////////////////
	core::Real calc_phipsi_RMSD(
		core::pose::Pose const& pose1,
		core::pose::Pose const& pose2,
		ObjexxFCL::FArray1D_bool const & res_subset
	) const;


	//calculaer different metrics for the interface of the given pose;
	std::map < std::string, core::Real >
	calc_interface_metrics
	(core::pose::Pose & pose, core::Size rb_jump, core::scoring::ScoreFunctionOP scorefxn);

	/////////////////////////////////////////
	// Calculate peptide score with and w/o
	// fa_ref sequence reference energy
	/////////////////////////////////////////
	void
	calc_pep_scores
	( core::pose::Pose const & pose, core::Real& pepScore, core::Real& pepScore_noref ) const;

	// ======== Accessor Methods ========

	void set_flags(FlexPepDockingFlagsCOP flags){
		flags_ = flags;
	}


private:
	FlexPepDockingFlagsCOP flags_;

	// helper method to check if two given residues are within "threashold" of each other
	bool isInContact ( core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real threashold ) const;


}; // end class FlexPepDockingPoseMetrics

} // end namespace flexPepDocking
} // end namespace protocols

#endif
