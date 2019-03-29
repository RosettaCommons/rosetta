// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/EntropyEstimator.cc
///
/// @brief
/// @author Hahnbeom Park and Frank DiMaio

#ifndef INCLUDED_protocols_ligand_docking_GALigandDock_EntropyEstimator_hh
#define INCLUDED_protocols_ligand_docking_GALigandDock_EntropyEstimator_hh

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <map>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

struct ChiInfo
{
	core::Real E;
	utility::vector1< core::Real > ligchis;
	std::map< core::Size, utility::vector1< core::Real > > recchis;
};

/// @brief
/// Estimates entropy change by ligand binding using short MC simulation
/// @details
/// This class takes a full complex pose and ligand seqpos index to calculate
/// entropy change upon ligand binding. MC simulation perturbs chi angles defined for
/// the ligand at free state, and estimates effective entropy loss by binding
/// by processing torsion probability during simulation assuming that ligand gets
/// completely after binding.

class EntropyEstimator //this is not a mover
{
public:

	EntropyEstimator( core::scoring::ScoreFunctionOP sfxn,
		core::pose::Pose const & pose,
		core::Size const ligid
	);

	~EntropyEstimator(){}

	core::Real apply( core::pose::Pose const& pose ) const;

	void set_niter( core::Size setting ){ niter_ = setting; }

private:

	/// @brief Estimates per-chiangle weights used for per-chiangle-weighted entropy estimation mode
	void
	get_chi_weight( core::pose::Pose const &pose_ref );

	/// @brief torsion entropy calculation function inside estimate_Stors
	core::Real
	analyze_trajectory( utility::vector1< ChiInfo > const &chitrj,
		core::pose::Pose const &pose,
		utility::vector1< std::pair< core::Size, core::Size > > const &chidefs,
		core::Size const ligid,
		core::Real const Emin,
		core::Real const RT,
		bool const run_on_ligand,
		bool const run_on_receptor
	) const;

	/// @brief Runs MC and returns torsion entropy change
	core::Real
	estimate_Stors( core::pose::Pose pose,
		core::Size const ligid,
		utility::vector1< ChiInfo > &chitrj,
		utility::vector1< core::Size > const &flexscs,
		utility::vector1< std::pair< core::Size, core::Size > > const &chidefs,
		bool const run_on_ligand,
		bool const run_on_receptor
	) const;

	/// @brief main perturb function in MC
	void
	perturb( core::pose::Pose &pose,
		core::Size const ligid,
		core::Size const nligchi,
		utility::vector1< core::Size > const &flexscs,
		bool &pert_ligand
	) const;

	/// @brief initialize MC trj datastructure
	utility::vector1< ChiInfo >
	setup_trj( core::pose::Pose const &pose,
		utility::vector1< core::Size > &flexscs,
		utility::vector1< std::pair< core::Size, core::Size > > &chidefs
	) const;

	/// @brief get list of sidechains contacting to ligand
	utility::vector1< core::Size >
	get_contacting_reslist( core::pose::Pose const &pose,
		utility::vector1< std::pair< core::Size, core::Size > > &chidefs
	) const;

	core::Real
	get_temperature( core::Size const it ) const;

	/// @brief perturb ligand chis
	void
	update_chis( core::conformation::Residue const &rsd,
		utility::vector1< core::Real > &chis ) const;

	/// @brief perturb receptor sidechains
	void
	update_flexscs( core::pose::Pose const &pose,
		utility::vector1< std::pair< core::Size, core::Size > > const &chidefs,
		std::map< core::Size, utility::vector1< core::Real > > &chis
	) const;

	core::Size
	chis2rotid( utility::vector1< core::Real > const &chis ) const;

private:
	// Basic
	core::Size ligid_;
	core::scoring::ScoreFunctionOP sfxn_;
	core::Real run_apostate_;
	core::Real run_holostate_;

	// MC mover params
	core::Real P_randomize_;
	core::Real maxpert_;
	bool minimize_;

	// MC schedule parameters
	//core::Size report_every_iter_;
	core::Size sample_every_iter_;
	core::Real temp_i_;
	core::Real temp_f_;
	core::Size niter_;
	core::Size iter_collect_begin_;

	// dG weight params
	core::Real wRG_;
	core::Real wtors_;

	// weight on chis
	bool weighted_;
	utility::vector1< core::Size > chiweights_;


}; //class EntropyEstimator

}
}
}

#endif // INCLUDED_protocols_ligand_docking_GALigandDock_EntropyEstimator_hh
