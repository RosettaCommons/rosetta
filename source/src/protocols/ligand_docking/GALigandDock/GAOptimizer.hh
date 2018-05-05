// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock/GAOptimizer.hh
///
/// @brief
/// @author Hahnbeom Park and Frank DiMaio

#ifndef INCLUDED_protocols_ligand_docking_GALigandDock_GAOptimizer_hh
#define INCLUDED_protocols_ligand_docking_GALigandDock_GAOptimizer_hh

#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

// represent options for a single "phase" of ligand docking
struct GADockStageParams {
	GADockStageParams() :
		repeats(0), pool(0), rmsthreshold(0), pmut(0), maxiter(0), packcycles(0), smoothing(0)
	{}

	GADockStageParams(
		core::Size repeats_in, core::Size pool_in,
		core::Real rmsthreshold_in, core::Real pmut_in, core::Real maxiter_in, core::Real packcycles_in, core::Real smoothing_in,
		utility::vector1<core::Real> const &ramp_schedule_in
	) :
		repeats(repeats_in),
		pool(pool_in),
		rmsthreshold(rmsthreshold_in),
		pmut(pmut_in),
		maxiter(maxiter_in),
		packcycles(packcycles_in),
		smoothing(smoothing_in),
		ramp_schedule(ramp_schedule_in)
	{}

	core::Size repeats;  // repeat these parameters for 'n' generations
	core::Size pool;     // size of the pool

	core::Real rmsthreshold; // enforce dissimilarity using this RMS
	core::Real pmut;         // mutation probability for this stage
	core::Real maxiter;      // maximum minimize iterations
	core::Real packcycles;   // number of packer cycles (x nSCs)
	core::Real smoothing;    // grid smoothing
	utility::vector1<core::Real> ramp_schedule;  // fa_rep ramping schedule in optimization
};


class GAOptimizer  {
public:
	GAOptimizer( GridScorerOP grid );
	~GAOptimizer();

	void run( LigandConformers & genes );

	void set_protocol( utility::vector1< GADockStageParams > const &protocol_in ) { protocol_ = protocol_in; }

	void set_native( LigandConformer const native ) { nativegene_ = native; }

	// optimize one generation
	void optimize_generation(
		LigandConformers & genes,
		utility::vector1<core::Real> const &ramping );

	void show_status(
		LigandConformers & genes, std::string comment="",
		bool calculate_native_rmsd = true,
		bool verbose = false );

	void set_max_rot_cumulative_prob( core::Real newval ) { max_rot_cumulative_prob_ = newval; }
	void set_rot_energy_cutoff( core::Real newval ) { rot_energy_cutoff_ = newval; }
	void set_favor_native( core::Real newval ) { favor_native_ = newval; }
	void set_altcrossover( bool newval ) { altcrossover_ = newval; }

private:
	//// HELPER FUNCTIONS
	// set up rotamer set
	void
	initialize_rotamer_set_and_scores(
		LigandConformer lig
	);

	// reset tags for a generation
	void update_tags( LigandConformers &genes ) const;

	// generate putative next generation
	void next_generation( LigandConformers const & genes, LigandConformers & genes_new, core::Size, core::Real );

	// update our pool
	void update_pool( LigandConformers & genes, LigandConformers & genes_new, core::Size, core::Real );

private:
	//// DATA
	LigandConformer nativegene_; // -in:file:native (for reporting)

	// protocol
	GridScorerOP scorefxn_;
	utility::vector1< GADockStageParams > protocol_;
	bool altcrossover_;

	// rotamer data
	core::Real max_rot_cumulative_prob_;
	core::Real rot_energy_cutoff_;
	core::Real favor_native_;
	utility::vector1< PlaceableRotamers > rotamer_data_;
	RotamerPairEnergies rotamer_energies_;
};

typedef utility::pointer::shared_ptr< GAOptimizer > GAOptimizerOP;
typedef utility::pointer::shared_ptr< GAOptimizer const > GAOptimizerCOP;

}
}
}

#endif // INCLUDED_protocols_ligand_docking_GAOptimizer_HH
