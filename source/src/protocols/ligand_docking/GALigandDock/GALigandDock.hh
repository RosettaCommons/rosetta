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
#include <protocols/ligand_docking/GALigandDock/TorsionSampler.hh>

#include <protocols/ligand_docking/GALigandDock/GALigandDock.fwd.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <protocols/moves/Mover.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>

#include <core/kinematics/FoldTree.hh>

#include <queue>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

typedef std::map< std::pair < core::Size, core::Size >, std::vector < core::Size > > HbondMap;

// Enum for opt H mode
enum OptHMode {
	OPTH_NONE,
	OPTH_FULL,
	OPTH_NO_FLIP_HNQ,
	OPTH_FLEXIBLE_SIDECHAINS,
	OPTH_REDEFINE_SIDECHAINS
};

// helper class to manage multiple outputs
struct StructInfo {
	core::io::silent::SilentStructOP str;
	core::scoring::constraints::ConstraintSetOP cst;
	core::Real rms, E, dH, ligandscore, recscore, complexscore;
	core::Size ranking_prerelax;
	std::string ligandname;
};

class StructInfoComp {
public:
	bool operator() ( StructInfo &a, StructInfo &b ) {  return (a.E > b.E); }
};

class StructInfoCompdH {
public:
	bool operator() ( StructInfo &a, StructInfo &b ) {  return (a.dH > b.dH); }
};

/// @brief helper class to manage multiple outputs
class OutputStructureStore {
public:
	OutputStructureStore() {}

	void
	push( core::pose::Pose const &pose, core::Real E,
		core::Real rms=0.0,
		core::Real complexscore=0.0,
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
		newstruct.dH = complexscore - recscore - ligandscore;
		newstruct.recscore = recscore;
		newstruct.complexscore = complexscore;
		newstruct.ranking_prerelax = ranking_prerelax;
		newstruct.ligandname = ligandname;

		struct_store_.push( newstruct );
		struct_store_dH_.push( newstruct );
	}

	void
	pop( core::pose::Pose &pose, core::Real &E,
		core::Real &rms,
		core::Real &complexscore,
		core::Real &ligandscore,
		core::Real &recscore,
		core::Size &ranking_prerelax,
		std::string &ligandname
	)
	{
		struct_store_.top().str->fill_pose( pose );
		pose.constraint_set( struct_store_.top().cst );
		rms = struct_store_.top().rms;
		E = struct_store_.top().E;
		complexscore = struct_store_.top().complexscore;
		ligandscore = struct_store_.top().ligandscore;
		recscore = struct_store_.top().recscore;
		ranking_prerelax = struct_store_.top().ranking_prerelax;
		ligandname = struct_store_.top().ligandname;

		struct_store_.pop();
	}

	void
	dH_pop( core::pose::Pose &pose, core::Real &E,
		core::Real &rms,
		core::Real &complexscore,
		core::Real &ligandscore,
		core::Real &recscore,
		core::Size &ranking_prerelax,
		std::string &ligandname
	)
	{
		struct_store_dH_.top().str->fill_pose( pose );
		pose.constraint_set( struct_store_dH_.top().cst );
		rms = struct_store_dH_.top().rms;
		E = struct_store_dH_.top().E;
		complexscore = struct_store_dH_.top().complexscore;
		ligandscore = struct_store_dH_.top().ligandscore;
		recscore = struct_store_dH_.top().recscore;
		ranking_prerelax = struct_store_dH_.top().ranking_prerelax;
		ligandname = struct_store_dH_.top().ligandname;

		struct_store_dH_.pop();
	}

	core::pose::PoseOP
	pop( std::string metric="score" )
	{
		if ( !has_data() ) return nullptr;
		core::pose::PoseOP retval (new core::pose::Pose);

		core::Real rms, E, ligscore, recscore, complexscore;
		core::Size ranking_prerelax;
		std::string ligandname;
		if ( metric == "dH" ) {
			dH_pop(*retval, E, rms, complexscore, ligscore, recscore, ranking_prerelax, ligandname);
		} else {
			pop(*retval, E, rms, complexscore, ligscore, recscore, ranking_prerelax, ligandname );
		}

		core::pose::setPoseExtraScore( *retval, "ligandname", ligandname);
		core::pose::setPoseExtraScore( *retval, "lig_rms", rms);
		core::pose::setPoseExtraScore( *retval, "ligscore", ligscore );
		core::pose::setPoseExtraScore( *retval, "recscore", recscore );
		core::pose::setPoseExtraScore( *retval, "complexscore", complexscore );
		core::pose::setPoseExtraScore( *retval, "ranking_prerelax", ranking_prerelax );
		core::pose::setPoseExtraScore( *retval, "dH", complexscore-recscore-ligscore );

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

	void
	resize( core::Size n ) {
		std::priority_queue< StructInfo, std::vector<StructInfo> , StructInfoComp > struct_store_temp_;
		std::priority_queue< StructInfo, std::vector<StructInfo> , StructInfoCompdH > struct_store_dH_temp_;
		while ( struct_store_temp_.size()<n && struct_store_.size()>0 ) {
			struct_store_temp_.push( struct_store_.top() );
			struct_store_dH_temp_.push( struct_store_dH_.top() );
			struct_store_.pop();
			struct_store_dH_.pop();
		}
		struct_store_ = struct_store_temp_;
		struct_store_dH_ = struct_store_dH_temp_;
	}

private:
	std::priority_queue< StructInfo, std::vector<StructInfo> , StructInfoComp > struct_store_;
	std::priority_queue< StructInfo, std::vector<StructInfo> , StructInfoCompdH > struct_store_dH_;
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
		basic::datacache::DataMap & data
	) override;

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

	void
	eval_docked_pose_helper( core::pose::Pose &pose,
		utility::vector1< core::Size > const& lig_ids,
		utility::vector1< core::Size > &movable_scs );

	void
	eval_docked_pose( core::pose::Pose &pose,
		utility::vector1< core::Size > const& lig_ids );

	utility::vector1< core::Size >
	get_movable_scs( core::pose::Pose const &pose,
		GridScorerCOP gridscore,
		utility::vector1 < core::Size > const &lig_resnos ) const;

	void
	idealize_and_repack_pose( core::pose::Pose &pose,
		utility::vector1< core::Size > const &movable_scs,
		utility::vector1< core::Size > const &lig_resnos ) const;

	/// @brief generate object containing binding-motif search results
	LigandAligner
	setup_ligand_aligner( core::pose::Pose const & pose,
		utility::vector1< core::Size > const &lig_resnos,
		utility::vector1< core::Size > movable_scs_in_ref, // call by value
		utility::vector1< core::conformation::Residue > const &rsds_to_build_grids,
		utility::vector1< bool > const &use_sc_only_in_grid
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

	/// @brief load the initial pool
	void
	load_template_pool(
		LigandConformer const &gene_initial,
		LigandConformers &genes_sel,
		core::Size nsel,
		utility::vector1< ConstraintInfo > & template_cst_infos
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
		utility::vector1< core::Size > const &lig_resnos,
		utility::vector1< core::Size > const& moving_scs,
		bool simple=true
	) const;

	//calculates the fraction of matching atom types for two residues
	//useful for similarity of pseudo-symmetiric lipids
	core::Real
	match_by_atomtype(
		core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2
	);

	std::pair < core::Real, core::Real >
	compare_hbonds_to_native( HbondMap const& native_hbond_map,
		HbondMap const& lig_hbond_map
	) const;

	core::Real
	calculate_free_ligand_score( core::pose::Pose pose, // call by value
		utility::vector1< core::Size > const &lig_resnos ) const;

	void
	apply_coord_cst_to_sctip( core::pose::PoseOP pose,
		utility::vector1< core::Size > const& moving_scs
	) const;

	/// @brief pre cart-min ligand before docking
	void
	premin_ligand(
		core::pose::Pose &pose,
		utility::vector1 < core::Size > const &lig_resnos ) const;

	/// @brief  final optimization cycle
	void
	final_exact_cartmin(
		core::Size nneigh,
		LigandConformer &gene,
		core::pose::Pose &pose,
		core::pack::task::PackerTaskOP task
		//bool dualrelax = false
	);

	/// @brief  alternate final optimization cycle
	void
	final_exact_scmin(
		LigandConformer const &gene,
		core::pose::Pose &pose,
		core::pack::task::PackerTaskOP task
	);

	void
	final_exact_ligmin_helper(
		LigandConformer const & gene,
		core::pose::Pose &pose,
		core::kinematics::MoveMapOP movemap,
		core::scoring::ScoreFunctionOP scfxn_local,
		core::Real const& fa_rep_weight,
		core::Real const& coordinate_cst_weight,
		core::Real const& torlerance,
		core::Size const& maxiter
	);

	void
	final_exact_ligmin(
		LigandConformer const &gene,
		core::pose::Pose &pose,
		core::pack::task::PackerTaskOP task
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

	void
	final_optH(
		core::Size const &optH_mode,
		LigandConformer const &gene,
		core::pose::Pose &pose,
		core::pack::task::PackerTaskOP task,
		utility::vector1< core::Size > &contact_scs
	);

	/// @brief  get ligand residue ids from pose
	void
	get_ligand_resids(core::pose::Pose const &pose,
		utility::vector1 < core::Size >& lig_resids) const;

	core::pose::PoseOP
	replace_ligand(core::pose::Pose pose_complex, core::pose::Pose& pose_ligand, bool align=true) const;

	void
	add_macrocycle_constraints(
		core::pose::Pose &pose,
		utility::vector1< core::Size > const &ligids
	) const;

	core::Size
	auto_determine_optH_mode(
		core::pose::Pose const &pose,
		utility::vector1< core::Size > const & movable_scs
	) const;

private:
	core::scoring::ScoreFunctionOP scfxn_; // scorefunction to be used in docking
	core::scoring::ScoreFunctionOP scfxn_relax_; // scorefunction to be used in relax
	std::string runmode_; // preset modes; see setup_params_for_runmode
	std::string top_pose_metric_;
	bool debug_report_; // print out debug info

	core::pose::PoseOP pose_native_; // native pose (for reporting purposes only)

	// grid-building parameters
	core::Real grid_, padding_, hashsize_, grid_radius_;
	core::Size subhash_;
	bool exact_, debug_;   // debugging options
	core::Real fa_rep_grid_;
	core::Real grid_bound_penalty_;

	//  ... define the "movable" parts
	std::string ligid_, sidechains_;
	core::select::residue_selector::ResidueSelectorCOP frozen_residues_;
	core::Real sc_edge_buffer_;
	bool move_water_;

	// protocol options
	bool use_pharmacophore_, aligner_fastmode_;
	core::Real max_rot_cumulative_prob_, rot_energy_cutoff_;
	core::Real random_oversample_, reference_oversample_, reference_frac_, init_dens_weight_;
	bool reference_frac_auto_;
	std::string initial_pool_, reference_pool_, template_pool_; // pdbs to include in initial pool
	core::Size n_template_;
	bool premin_ligand_;
	bool sample_ring_conformers_;
	core::Real torsion_sampler_percentage_;
	TorsionSamplerCOP torsion_sampler_;
	core::Real contact_distance_;
	bool freeze_ligand_backbone_;
	bool freeze_ligand_;
	bool macrocycle_ligand_;
	core::Size shuffle_ligands_;

	core::Real skeleton_threshold_const_; //constant used for determining skeleton threshold
	core::Size neighborhood_size_; //size of neighborhood to search during erosion
	core::Real skeleton_radius_;
	std::string method_for_radius_;
	bool advanced_map_erosion_;
	bool print_initial_pool_;
	bool altcrossover_;
	bool single_mutation_;
	core::Real local_res_;

	bool calculate_native_density_;

	std::string final_exact_minimize_; // do the last iteration exactly?
	bool cartmin_lig_, min_neighbor_;  // more final relax properties

	bool final_solvate_;
	bool turnon_flexscs_at_relax_;
	bool redefine_flexscs_at_relax_;
	std::string fast_relax_script_file_; // script file for fast relax (EXACT ONLY!)
	std::vector< std::string > fast_relax_lines_; // in explicit texts
	core::Real favor_native_; // give a bonus to input rotamer
	bool optimize_input_H_; // optimize_h_mode_ at the beginning; goes to grid construction
	bool pre_optH_relax_;
	bool auto_final_optH_;
	core::Size final_optH_mode_;
	bool full_repack_before_finalmin_;
	core::Size nrelax_;
	core::Size nreport_;
	bool estimate_dG_;
	std::string entropy_method_;
	bool estimate_buns_;
	bool use_dalphaball_;
	utility::vector1< core::Size > hb_resids_;
	bool hb_resids_include_bb_;
	std::string hb_resids_metric_;

	bool use_mean_maxRad_;
	core::Real stdev_multiplier_;
	core::Real multi_ligands_maxRad_;


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
	utility::vector1< std::string > ligand_file_list_;

	// atom ids to align back to reference after updates
	utility::vector1< core::id::AtomID > align_reference_atom_ids_;

	// handle multiple outputs
	OutputStructureStore remaining_outputs_;

	//used for rooting initial_pool_ pdbs to same as input
	bool is_virtual_root_;
	core::kinematics::FoldTree input_fold_tree_;

	bool has_density_map_;
	bool output_ligand_only_;
};


}
}
}

#endif // INCLUDED_protocols_ligand_docking_GALigandDock_HH
