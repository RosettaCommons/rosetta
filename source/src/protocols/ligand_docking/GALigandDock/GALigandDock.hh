// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock.cc
///
/// @brief
/// @author Hahnbeom Park and Frank DiMaio

#ifndef INCLUDED_protocols_ligand_docking_GALigandDock_GALigandDock_hh
#define INCLUDED_protocols_ligand_docking_GALigandDock_GALigandDock_hh

#include <protocols/ligand_docking/GALigandDock/LigandConformer.hh>
#include <protocols/ligand_docking/GALigandDock/GAOptimizer.hh>
#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>

#include <protocols/ligand_docking/GALigandDock/GALigandDock.fwd.hh>

#include <protocols/moves/Mover.hh>

#include <queue>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

// docking protocol
class GALigandDock : public protocols::moves::Mover {
public:
	GALigandDock();

	// main apply of GA ligand docking
	void apply( Pose & pose ) override;

	// as this is a one->many mover, use this function to return multiple outputs from a single call to apply
	core::pose::PoseOP get_additional_output() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		filters::Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	// generate an initial set of perturbed structures
	LigandConformers
	generate_perturbed_structures(
		LigandConformer const &gene_initial,
		GridScorerOP gridscorer,
		core::Size
	) const;

	// create and initialize the grid optimizer
	GAOptimizerOP
	get_optimizer(
		LigandConformer const &gene_initial,
		GridScorerOP gridscorer
	) const;

	// final optimization cycle
	void
	final_exact_cartmin(
		core::Size nneigh,
		LigandConformers & genes,
		utility::vector1< core::pose::PoseOP > &poses
	);

	// alternate final optimization cycle
	void
	final_exact_scmin(
		LigandConformers & genes,
		utility::vector1< core::pose::PoseOP > &poses
	);

	//
	void
	final_exact_rtmin(
		LigandConformers & genes,
		utility::vector1< core::pose::PoseOP > &poses
	);

	//
	void
	final_exact_nirtmin(
		LigandConformers & genes,
		utility::vector1< core::pose::PoseOP > &poses
	);

private:
	core::scoring::ScoreFunctionOP scfxn_; // scorefunction to be used in docking

	core::pose::PoseOP pose_native_; // native pose (for reporting purposes only)

	// grid-building parameters
	core::Real grid_, padding_, hashsize_;
	core::Size subhash_;
	bool exact_, debug_;   // debugging options

	std::string runmode_;

	//  ... define the "movable" parts
	std::string ligid_, sidechains_;
	core::Real sc_edge_buffer_;

	// protocol options
	bool altcrossover_;
	core::Real max_rot_cumulative_prob_, rot_energy_cutoff_;
	core::Real init_oversample_;
	std::string initial_pool_; // pdbs to include in initial pool
	std::string final_exact_minimize_; // do the last iteration exactly?
	std::string fast_relax_script_file_; // script file for fast relax (EXACT ONLY!)
	core::Real favor_native_; // give a bonus to input rotamer

	// dump score or pose for all the samples generated
	std::string report_all_samples_;
	bool rtmin_nonideal_;

	// per-cycle defaults
	core::Size ngen_, npool_;
	core::Real rmsdthreshold_; // enforce dissimilarity using this RMS
	core::Real pmut_;         // mutation probability for this stage
	core::Real maxiter_;      // maximum minimize iterations
	core::Real packer_cycles_;   // number of packer cycles (x nSCs)
	core::Real smoothing_;    // grid smoothing
	utility::vector1<core::Real> ramp_schedule_;  // fa_rep ramping schedule in optimization

	// detailed protocol
	utility::vector1< GADockStageParams > protocol_;

	// handle multiple outputs
	std::queue<LigandConformer> remaining_compact_outputs_;
	std::queue<core::pose::PoseOP> remaining_expanded_outputs_;
};

}
}
}

#endif // INCLUDED_protocols_ligand_docking_GALigandDock_HH
