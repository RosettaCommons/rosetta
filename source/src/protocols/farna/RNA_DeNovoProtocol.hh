// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_DeNovo_Protocol.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_RNA_DeNovoProtocol_HH
#define INCLUDED_protocols_rna_RNA_DeNovoProtocol_HH

#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/rna/RNA_DataReader.fwd.hh>
#include <protocols/farna/RNA_Fragments.fwd.hh>
#include <protocols/farna/RNA_StructureParameters.hh>
#include <protocols/farna/RNA_LoopCloser.fwd.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/farna/RNA_FragmentMover.fwd.hh>
#include <protocols/farna/RNA_Minimizer.fwd.hh>
#include <protocols/farna/RNA_Relaxer.fwd.hh>
#include <protocols/farna/RNA_ChunkLibrary.fwd.hh>

//Oooh.
#include <ObjexxFCL/FArray1D.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <list>

#include <core/io/silent/silent.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace farna {

/// @brief The RNA de novo structure modeling protocol
class RNA_DeNovoProtocol: public protocols::moves::Mover {
public:

	/// @brief Construct the protocol object given
	/// the RNA fragment library to use.
	RNA_DeNovoProtocol(
										 Size const nstruct,
										 std::string const silent_file,
										 bool const heat_structure = true,
										 bool const minimize_structure = false,
										 bool const relax_structure = false,
										 bool const allow_bulge = false );

	~RNA_DeNovoProtocol();

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const;

	/// @brief Apply the RNA denovo modeling protocol to the input pose
	void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	virtual void show(std::ostream & output=std::cout) const;

	void
	set_dump_pdb( bool const setting ){ dump_pdb_ = setting; };

	void
	set_temperature( core::Real const setting ){ m_Temperature_ = setting; };

	void
	set_jump_library_file( std::string const jump_library_file ) {
		jump_library_file_ = jump_library_file;
	}

	void
	set_all_rna_fragments_file( std::string const file ) {
		all_rna_fragments_file_ = file;
	}

	void
	set_vall_torsions_file( std::string const file ) {
		all_rna_fragments_file_ = file;
	}

	void
	set_rna_params_file( std::string const file ) {
		rna_params_file_ = file;
	}

	void
	set_rna_data_file( std::string const file ) {
		rna_data_file_ = file;
	}

	void
	set_chunk_pdb_files( utility::vector1< std::string > const & chunk_pdb_files ) {
		chunk_pdb_files_ = chunk_pdb_files;
	}

	void
	set_chunk_silent_files( utility::vector1< std::string > const & chunk_silent_files ) {
		chunk_silent_files_ = chunk_silent_files;
	}

	void
	set_input_res( utility::vector1< Size > const & input_res ) {
		input_res_ = input_res;
	}

	void
	ignore_secstruct( bool const setting ) { ignore_secstruct_ = setting; }

	// No longer works -- need to specify allow_insert from a "params file"
	//	void
	//	set_allow_insert( FArray1D <bool> const & allow_insert  ){ allow_insert_ = allow_insert; }

    void
    jump_change_frequency( core::Real const value ){ jump_change_frequency_ = value; }

    void
    set_close_loops( bool const setting ){
        close_loops_at_end_ = setting;
        if ( close_loops_at_end_ ) binary_rna_output_ = true;
    }

    void
    set_close_loops_after_each_move( bool const setting ){ close_loops_after_each_move_ = setting; }

    void
    simple_rmsd_cutoff_relax( bool const setting ){ simple_rmsd_cutoff_relax_ = setting; }

    void
    output_lores_silent_file( bool const setting ){ output_lores_silent_file_ = setting; }

	void
	set_filter_lores_base_pairs( bool const setting ){ filter_lores_base_pairs_ = setting; }

	void
	set_filter_lores_base_pairs_early( bool const setting ){
		filter_lores_base_pairs_early_ = setting;
		if (filter_lores_base_pairs_early_) filter_lores_base_pairs_ = true;
	}

	void
	set_filter_chain_closure( bool const setting ){ filter_chain_closure_ = setting; }

	void
	set_filter_chain_closure_distance( core::Real const setting ){ filter_chain_closure_distance_ = setting; }

	void
	set_filter_chain_closure_halfway( bool const setting ){ filter_chain_closure_halfway_ = setting; }

	void
	set_binary_rna_output( bool const setting ){ binary_rna_output_ = setting; }

	void
	set_suppress_bp_constraint( core::Real const setting ){ suppress_bp_constraint_ = setting; }

	void
	set_lores_scorefxn( std::string const & lores_scorefxn ){ lores_scorefxn_ = lores_scorefxn; }

	void
	set_vary_bond_geometry( bool const setting ){
		vary_bond_geometry_ = setting;
		if ( vary_bond_geometry_ ) binary_rna_output_ = true;
	}

	void
	output_to_silent_file( core::pose::Pose & pose, std::string const & silent_file, std::string const & out_file_tag, bool const score_only = false ) const;

	void
	align_and_output_to_silent_file( core::pose::Pose & pose, std::string const & silent_file, std::string const & out_file_tag ) const;

	void
	set_staged_constraints( bool const setting ){ staged_constraints_ = setting; }

	void
	set_allow_consecutive_bulges( bool const setting ){ allow_consecutive_bulges_ = setting; };

	void
	set_chainbreak_weight( core::Real setting ){ chainbreak_weight_ = setting; }

	void
	set_linear_chainbreak_weight( core::Real setting ){ linear_chainbreak_weight_ = setting; }

	void
	set_allowed_bulge_res( utility::vector1< core::Size > const & setting ){ allowed_bulge_res_ = setting; };

	void
	set_move_first_rigid_body( bool const setting ){ move_first_rigid_body_ = setting; }

	void
	set_root_at_first_rigid_body( bool const setting ){ root_at_first_rigid_body_ = setting; }

	void
	set_output_filters( bool const setting ){ output_filters_ = setting; }

	void
	set_autofilter( bool const setting ){ autofilter_ = setting; }

	void
	set_monte_carlo_cycles( Size const setting ){ monte_carlo_cycles_ = setting;  user_defined_cycles_ = true; }

	void
	set_rounds( Size const setting ){ rounds_ = setting; }

	void
	set_extra_minimize_res( utility::vector1< core::Size > setting );

	void
	set_extra_minimize_chi_res( utility::vector1< core::Size > setting );

	void
	set_refine_pose_list( utility::vector1<core::pose::PoseOP> const & setting );

	void
	set_refine_pose( Size const setting ){ refine_pose_ = setting; }

	void
	set_bps_moves( Size const setting ){ bps_moves_ = setting; }

	void
	set_minimizer_use_coordinate_constraints( Size const setting ){ minimizer_use_coordinate_constraints_ = setting; }

private:

	void
	initialize_movers( core::pose::Pose & pose );

	void
	setup_monte_carlo_cycles( core::pose::Pose const & pose );

	utility::vector1< Size >
	get_moving_res( core::pose::Pose const & pose ) const;

	void
	initialize_lores_silent_file();

	void
	initialize_constraints( core::pose::Pose & pose );

	void
	initialize_scorefxn( core::pose::Pose & pose );

	void
	initialize_tag_is_done();

	void
	setup_rigid_body_mover( core::pose::Pose const & pose, core::Size const r );

	void
	final_score( core::pose::Pose & pose );

	void
	output_silent_struct( core::io::silent::SilentStruct & s,
												core::io::silent::SilentFileData & silent_file_data,
												std::string const & silent_file,
												core::pose::Pose & pose,
												std::string const out_file_tag,
												bool const score_only = false ) const;

	void
	do_random_moves( core::pose::Pose & pose );

	void
	randomize_rigid_body_orientations( core::pose::Pose & pose );

	void
	update_denovo_scorefxn_weights( Size const r );

	Size
	figure_out_constraint_separation_cutoff( Size const r, Size const  max_dist );

	void
	update_pose_constraints( Size const r, core::pose::Pose & pose );

	void
	update_frag_size( Size const r );

	void
	random_fragment_trial( core::pose::Pose & pose );

	bool
	random_chunk_trial( core::pose::Pose & pose );

	void
	random_jump_trial( core::pose::Pose & pose );

	void
	RNA_move_trial( core::pose::Pose & pose );

	void
	add_number_base_pairs( core::pose::Pose const & pose_input, core::io::silent::SilentStruct & s ) const;

	void
	add_number_native_base_pairs( core::pose::Pose & pose, core::io::silent::SilentStruct & s ) const;

	void
	calc_rmsds( core::io::silent::SilentStruct & s, core::pose::Pose & pose, std::string const & out_file_tag ) const;

	bool
	check_score_filter( core::Real const lores_score_, std::list< core::Real > & all_lores_score_ );

	void
	apply_chem_shift_data(core::pose::Pose & pose, std::string const out_file_tag);

	void
	add_chem_shift_info(core::io::silent::SilentStruct & silent_struct, core::pose::Pose const & const_pose) const;

private:

	// protocol-specific data ... need to be specified as input.
	Size const nstruct_;
	Size rounds_;
	Size monte_carlo_cycles_;
	Size const monte_carlo_cycles_max_default_;
	bool user_defined_cycles_;
	std::string all_rna_fragments_file_;
	std::string const silent_file_;
	std::string lores_silent_file_;
	bool const heat_structure_;
	bool dump_pdb_;
	bool const minimize_structure_;
	bool const relax_structure_;
	bool ignore_secstruct_;

	bool do_close_loops_;
	bool close_loops_at_end_;
	bool close_loops_in_last_round_;
	bool close_loops_after_each_move_;

	bool simple_rmsd_cutoff_relax_;
	bool allow_bulge_, allow_consecutive_bulges_;

  bool const use_chem_shift_data_;

	// parameters
	core::Real m_Temperature_; // default temperature for monte carlo
	Size frag_size_;

	protocols::moves::MonteCarloOP monte_carlo_;

	protocols::farna::RNA_FragmentsOP all_rna_fragments_;
	protocols::farna::RNA_FragmentMoverOP rna_fragment_mover_;
	protocols::farna::RNA_LoopCloserOP rna_loop_closer_;
	protocols::farna::RNA_ChunkLibraryOP rna_chunk_library_;
	protocols::farna::RNA_MinimizerOP rna_minimizer_;
	protocols::farna::RNA_RelaxerOP rna_relaxer_;
	protocols::rigid::RigidBodyPerturbMoverOP rigid_body_mover_;

	std::string rna_params_file_;
	std::string rna_data_file_;
	std::string jump_library_file_;
	// This object will actually hold the jump library...
	RNA_StructureParametersOP rna_structure_parameters_;
	core::io::rna::RNA_DataReaderOP rna_data_reader_;

	bool output_lores_silent_file_;

	bool filter_lores_base_pairs_;
	bool filter_lores_base_pairs_early_;

	bool filter_chain_closure_;
	core::Real filter_chain_closure_distance_;
	bool filter_chain_closure_halfway_;

	bool vary_bond_geometry_;
	bool binary_rna_output_;

	core::scoring::constraints::ConstraintSetOP constraint_set_;

	core::Real jump_change_frequency_;

	std::map< std::string, bool > tag_is_done_;

	std::string lores_scorefxn_;
	core::scoring::ScoreFunctionOP denovo_scorefxn_;
	core::scoring::ScoreFunctionOP hires_scorefxn_;
	core::scoring::ScoreFunctionOP chem_shift_scorefxn_;
	core::scoring::ScoreFunctionOP initial_denovo_scorefxn_;
	core::scoring::ScoreFunctionOP final_scorefxn_;
	core::scoring::rna::RNA_LowResolutionPotential local_rna_low_resolution_potential_;

	utility::vector1< std::string > chunk_pdb_files_;
	utility::vector1< std::string > chunk_silent_files_;
	utility::vector1< core::Size > input_res_;
	utility::vector1< core::Size > allowed_bulge_res_;
	utility::vector1< core::Size > extra_minimize_res_;
	utility::vector1< core::Size > extra_minimize_chi_res_;
	core::Real chunk_coverage_;

	bool staged_constraints_;

	core::Real chainbreak_weight_;
	core::Real linear_chainbreak_weight_;

	bool titrate_stack_bonus_;
	bool move_first_rigid_body_;
	bool root_at_first_rigid_body_;
	core::Real suppress_bp_constraint_;

	bool output_filters_;
	core::Real lores_score_early_;
	core::Real lores_score_final_;
	bool autofilter_;
	core::Real autofilter_score_quantile_;
	std::list< core::Real > all_lores_score_final_;
	bool refine_from_silent_;
	utility::vector1<core::pose::PoseOP> refine_pose_list_;
	bool refine_pose_;
	bool bps_moves_;
	bool minimizer_use_coordinate_constraints_;

}; // class RNA_DeNovoProtocol

std::ostream &operator<< ( std::ostream &os, RNA_DeNovoProtocol const &mover );

} //farna
} //protocols

#endif
