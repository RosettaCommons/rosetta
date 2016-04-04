// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief take out the "compute_frag_bias()" from CartesianSampler
/// @author Ray Wang wangyr@uw.edu
//
#ifndef INCLUDED_protocols_hybridization_FragmentBiasAssigner_hh
#define INCLUDED_protocols_hybridization_FragmentBiasAssigner_hh

#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragSet.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>

#include <utility/vector1.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace hybridization {

using namespace core;

class FragmentBiasAssigner : public utility::pointer::ReferenceCount
{
public:
	// constructor
	//FragmentBiasAssigner();
	FragmentBiasAssigner( pose::Pose &pose );

	void init( pose::Pose &pose );

	// everytime when you call these function the probability is going to be additive to the frag_bias_ container
	void compute_frag_bias( utility::vector1<numeric::random::WeightedSampler> &frag_bias, // output
		pose::Pose &pose,
		utility::vector1<core::fragment::FragSetOP> fragments );

	void include_residues( std::set< core::Size > residues_to_include );

	void exclude_residues( std::set< core::Size > residues_to_exclude );

	/////////////////////////////////////////////////////////////
	void uniform();

	////////////////////////////////////
	// ray's protocol to select residues to refine
	void automode( pose::Pose &pose,
		Real score_cut );

	/////////////////////////////////////////////////////////////
	void user( std::set<core::Size> user_pos,
		protocols::loops::LoopsOP loops );

	////////////////////////////////////
	// individual method

	void density_nbr( pose::Pose &pose );

	void rama( pose::Pose &pose,
		core::Real weight=0.2 );

	void geometry( pose::Pose &pose,
		core::Real weight=1.0 );

	void density( pose::Pose &pose );

	////////////////////////////////
	void old_rama( pose::Pose &pose );
	void chainbreak( pose::Pose &pose );
	void bfactors( pose::Pose &pose );
	void fragbias_reporter( pose::Pose &pose );
	void cumulate_probability(){ cumulative_=true; }

	void set_rsd_wdw_to_assign_prob( int wdw=0 ){ rsd_wdw_size_=wdw; }
	void set_wdw_to_freeze( int wdw=0 ){ wdw_to_freeze_=wdw; }
	void set_score_threshold( Real threshold=123456789 ){ score_threshold_=threshold; }

private:
	// functions
	void cal_perrsd_score( pose::Pose &pose,
						   core::scoring::ScoreType const &score_type,
		utility::vector1<core::Real> &perrsd_score,
		Real weight );

	void cal_zscore( utility::vector1<core::Real> const &input_v,
		utility::vector1<core::Real> &zscore_v,
		bool negating=false);

	// This function calls assign_prob_with_rsd_wdw(rsn) to assign probability to the residue with a window controlled by "rsd_wdw_size_".
	void assign_fragprobs( utility::vector1<core::Real> const &perrsd_score,
		Real threshold );

	void assign_prob_with_rsd_wdw( int rsn );


	// variables
	core::Size nres_, n_symm_subunit_;
	int wdw_to_freeze_;
	int rsd_wdw_size_; // to assign prob
	Real score_threshold_; // to assign prob
	bool cumulative_;
	bool fragProbs_assigned_;
	core::conformation::symmetry::SymmetryInfoCOP symminfo_;

	////////////////////////////////////
	// containers
	// the per-residue based container for all the methods to add values into it
	utility::vector1<core::Real> fragmentProbs_;

	// the per residue container for score terms
	utility::vector1<core::Real> perrsd_dens_;
	utility::vector1<core::Real> perrsd_nbrdens_;
	utility::vector1<core::Real> perrsd_rama_;
	utility::vector1<core::Real> perrsd_geometry_;

	//utility::vector1<numeric::random::WeightedSampler> frag_bias_;
}; // class FragmentBiasAssigner

} // hybridization
} // protocol
#endif
