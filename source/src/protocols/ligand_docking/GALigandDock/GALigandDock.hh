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
#include <protocols/ligand_docking/GALigandDock/LigandAligner.hh>

#include <protocols/ligand_docking/GALigandDock/GALigandDock.fwd.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <protocols/moves/Mover.hh>

#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <queue>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

// helper class to manage multiple outputs
struct StructInfo {
	core::io::silent::SilentStructOP str;
	core::scoring::constraints::ConstraintSetOP cst;
	core::Real rms, E, ligandscore, recscore;
	core::Size ranking_prerelax;
	std::string ligandname;
};

class StructInfoComp {
public:
	bool operator() ( StructInfo &a, StructInfo &b ) {  return (a.E > b.E); }
};

/// @brief helper class to manage multiple outputs
class OutputStructureStore {
public:
	OutputStructureStore() {}

	void
	push( core::pose::Pose const &pose, core::Real E,
		core::Real rms=0.0,
		core::Real ligandscore=0.0,
		core::Real recscore=0.0,
		core::Size ranking_prerelax=0,
		std::string ligandname=""
	) {
		StructInfo newstruct;

		core::io::silent::SilentFileOptions opts;
		newstruct.str = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary", opts);
		newstruct.str->fill_struct( pose );
		newstruct.rms = rms;
		newstruct.ligandscore = ligandscore;
		newstruct.cst = pose.constraint_set()->clone();
		newstruct.E = E;
		newstruct.recscore = recscore;
		newstruct.ranking_prerelax = ranking_prerelax;
		newstruct.ligandname = ligandname;

		struct_store_.push( newstruct );
	}

	void
	pop( core::pose::Pose &pose, core::Real &E,
		core::Real &rms, core::Real &ligandscore,
		core::Real &recscore,
		core::Size &ranking_prerelax,
		std::string &ligandname
	)
	{
		struct_store_.top().str->fill_pose( pose );
		pose.constraint_set( struct_store_.top().cst );
		rms = struct_store_.top().rms;
		E = struct_store_.top().E;
		ligandscore = struct_store_.top().ligandscore;
		recscore = struct_store_.top().recscore;
		ranking_prerelax = struct_store_.top().ranking_prerelax;
		ligandname = struct_store_.top().ligandname;

		struct_store_.pop();
	}

	core::pose::PoseOP
	pop()
	{
		if ( !has_data() ) return nullptr;
		core::pose::PoseOP retval (new core::pose::Pose);

		core::Real rms, E, ligscore, recscore;
		core::Size ranking_prerelax;
		std::string ligandname;
		pop(*retval, E, rms, ligscore, recscore, ranking_prerelax, ligandname );

		core::pose::setPoseExtraScore( *retval, "ligandname", ligandname);
		core::pose::setPoseExtraScore( *retval, "lig_rms", rms);
		core::pose::setPoseExtraScore( *retval, "ligscore", ligscore );
		core::pose::setPoseExtraScore( *retval, "recscore", recscore );
		core::pose::setPoseExtraScore( *retval, "ranking_prerelax", ranking_prerelax );

		return retval;
	}

	void
	clear(){
		struct_store_ = std::priority_queue< StructInfo, std::vector<StructInfo> , StructInfoComp >();
	}

	bool
	has_data( ) {
		return ( struct_store_.size() > 0 );
	}

	core::Size
	size( ){ return struct_store_.size(); }

private:
	std::priority_queue< StructInfo, std::vector<StructInfo> , StructInfoComp > struct_store_;
};


/// @brief
/// Ligand Docking protocol using Genetic Algorithm
/// @details
/// Runs ligand docking using pre-computed full-atom grid score, binding-motif search,
/// and genetic algorithm search on Ligand & flexible sidechains at receptor conformations.
/// This docking method supports full on-the-fly search of ligand conformation hence
/// ligand "conformer" generation is not required.

class GALigandDock : public protocols::moves::Mover {
public:
	GALigandDock();

	/// @brief main apply of GA ligand docking
	void apply( Pose & pose ) override;

	/// @details as this is a one->many mover, use this function to return multiple outputs from a single call to apply
	core::pose::PoseOP get_additional_output() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief parse options based on umbrella runmode option
	void setup_params_for_runmode( std::string runmode ); // access from pyrosetta etc...

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

	/// @brief main function running on single ligand type
	core::pose::Pose
	run_docking( LigandConformer const &gene_initial,
		GridScorerOP gridscore,
		LigandAligner &aligner,
		OutputStructureStore &outputs );

	utility::vector1< core::Size >
	get_movable_scs( core::pose::Pose const &pose,
		GridScorerCOP gridscore,
		core::Size const lig_resno ) const;

	void
	idealize_and_repack_pose( core::pose::Pose &pose,
		utility::vector1< core::Size > const &movable_scs,
		core::Size const lig_resno ) const;

	/// @brief generate object containing binding-motif search results
	LigandAligner
	setup_ligand_aligner( core::pose::Pose const & pose,
		core::Size const lig_resno,
		utility::vector1< core::Size > movable_scs_in_ref // call by value
	) const;

	/// @brief generate an initial set of perturbed structures
	LigandConformers
	generate_perturbed_structures(
		LigandConformer const &gene_initial,
		GridScorerOP gridscorer,
		core::Size npool,
		LigandAligner aligner //call by value
	) const;

	/// @brief replace ligand to a new ideal res named by ligand_name
	core::pose::PoseOP
	make_starting_pose_for_virtual_screening( core::pose::Pose const &pose_apo,
		core::Size const &lig_resno,
		std::string const ligand_name
	) const;

	/// @brief load the initial pool
	void
	load_initial_pool(
		LigandConformer const &gene_initial,
		LigandConformers &genes_sel
	) const;

	// load the initial pool
	void
	load_reference_pool(
		LigandConformer const &gene_initial,
		utility::vector1< ConstraintInfo > & ref_ligs
	) const;

	/// @brief create and initialize the grid optimizer
	GAOptimizerOP
	get_optimizer(
		LigandConformer const &gene_initial,
		GridScorerOP gridscorer
	) const;

	core::Real
	calculate_free_receptor_score( core::pose::Pose pose, // call by value
		core::Size const lig_resno,
		utility::vector1< core::Size > const& moving_scs,
		bool simple=true
	) const;

	core::Real
	calculate_free_ligand_score( core::conformation::Residue const ligand ) const;

	void
	apply_coord_cst_to_sctip( core::pose::PoseOP pose,
		utility::vector1< core::Size > const& moving_scs
	) const;

	/// @brief pre cart-min ligand before docking
	void
	premin_ligand( core::pose::Pose &pose, core::Size const lig_resno ) const;

	/// @brief  final optimization cycle
	void
	final_exact_cartmin(
		core::Size nneigh,
		LigandConformer &gene,
		core::pose::Pose &pose
		//bool dualrelax = false
	);

	/// @brief  alternate final optimization cycle
	void
	final_exact_scmin(
		LigandConformer const &gene,
		core::pose::Pose &pose
	);

	/// @brief  alternate final optimization cycle
	void
	final_cartligmin(
		LigandConformer const &gene,
		core::pose::Pose &pose
	);

	/// @brief  solvate pose before final relax
	void
	final_solvate(
		LigandConformer &gene,
		core::pose::Pose &pose
	);

private:
	core::scoring::ScoreFunctionOP scfxn_; // scorefunction to be used in docking
	core::scoring::ScoreFunctionOP scfxn_relax_; // scorefunction to be used in relax
	std::string runmode_; // preset modes; see setup_params_for_runmode

	core::pose::PoseOP pose_native_; // native pose (for reporting purposes only)

	// grid-building parameters
	core::Real grid_, padding_, hashsize_;
	core::Size subhash_;
	bool exact_, debug_;   // debugging options
	core::Real fa_rep_grid_;
	core::Real grid_bound_penalty_;

	//  ... define the "movable" parts
	std::string ligid_, sidechains_;
	core::Real sc_edge_buffer_;
	bool move_water_;

	// protocol options
	bool use_pharmacophore_;
	bool altcrossover_;
	core::Real max_rot_cumulative_prob_, rot_energy_cutoff_;
	core::Real random_oversample_, reference_oversample_, reference_frac_;
	bool reference_frac_auto_;
	std::string initial_pool_, reference_pool_; // pdbs to include in initial pool
	bool premin_ligand_;
	bool sample_ring_conformers_;

	std::string final_exact_minimize_; // do the last iteration exactly?
	bool cartmin_lig_, min_neighbor_;  // more final relax properties

	bool final_solvate_;
	bool redefine_flexscs_at_relax_;
	std::string fast_relax_script_file_; // script file for fast relax (EXACT ONLY!)
	std::vector< std::string > fast_relax_lines_; // in explicit texts
	core::Real favor_native_; // give a bonus to input rotamer
	bool optimize_input_H_; // optimize_h_mode_ at the beginning; goes to grid construction
	bool full_repack_before_finalmin_;
	core::Size nrelax_;
	core::Size nreport_;
	bool estimate_dG_;

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

	// handle multiple ligand types for virtual screening
	utility::vector1< std::string > multiple_ligands_;

	// handle multiple outputs
	OutputStructureStore remaining_outputs_;
};

}
}
}

#endif // INCLUDED_protocols_ligand_docking_GALigandDock_HH
