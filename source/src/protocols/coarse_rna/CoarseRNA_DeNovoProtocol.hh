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
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_coarse_rna_CoarseRNA_DeNovoProtocol_HH
#define INCLUDED_protocols_coarse_rna_CoarseRNA_DeNovoProtocol_HH

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/farna/movers/RNA_FragmentMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/coarse_rna/MultipleDomainMover.fwd.hh>
#include <protocols/farna/libraries/RNA_ChunkLibrary.fwd.hh>
#include <protocols/farna/setup/RNA_DeNovoPoseSetup.fwd.hh>
#include <core/io/rna/RNA_DataReader.fwd.hh>
#include <protocols/coarse_rna/CoarseRNA_LoopCloser.fwd.hh>

#include <utility/vector1.hh>

//Oooh.

//// C++ headers
#include <string>

//Auto Headers
namespace protocols {
namespace coarse_rna {

/// @brief The RNA de novo structure modeling protocol
class CoarseRNA_DeNovoProtocol: public protocols::moves::Mover {
public:

	/// @brief Construct the protocol object given
	/// the RNA fragment library to use.
	CoarseRNA_DeNovoProtocol(
		Size const nstruct,
		Size const monte_carlo_cycles,
		std::string const silent_file );

	~CoarseRNA_DeNovoProtocol();

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const;

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	void
	set_dump_pdb( bool const setting ){ dump_pdb_ = setting; };

	void
	set_force_ideal_chainbreak( bool const & setting ){ force_ideal_chainbreak_ = setting; }

	void
	set_check_pairing_dists( bool const & setting ){ check_pairing_dists_ = setting; }

	void
	set_add_base_pair_constraints( bool const & setting ){ add_base_pair_constraints_ = setting; }

	void
	set_jump_library_file( std::string const jump_library_file ) {
		jump_library_file_ = jump_library_file;
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
	set_chunk_silent_files( utility::vector1< std::string > const & chunk_silent_files ) {
		chunk_silent_files_ = chunk_silent_files;
	}

	void
	set_input_res(  utility::vector1< core::Size > const & setting ){ input_res_ = setting; }


	void
	set_lores_scorefxn( std::string const & lores_scorefxn ){ lores_scorefxn_ = lores_scorefxn; }

	void
	output_to_silent_file( core::pose::Pose & pose, std::string const & silent_file, std::string const & out_file_tag, bool const score_only = false ) const;

	void
	set_temperature( core::Real const & setting){ m_Temperature_ = setting; }

	void
	set_staged_constraints( bool const & setting){ staged_constraints_ = setting; }

	void
	set_sim_anneal( bool const & setting){ sim_anneal_ = setting; }

	void
	set_close_loops( bool const & setting){ close_loops_ = setting; }

	void
	set_choose_best_solution( bool const & setting){ choose_best_solution_ = setting; }

	void
	set_freeze_domains( bool const & setting){ freeze_domains_ = setting; }

private:

	void
	initialize_movers( core::pose::Pose & pose );

	void
	initialize_tag_is_done();

	core::Real
	get_temperature( Size const & r, Size const & rounds ) const;

	void
	do_random_fragment_insertions( core::pose::Pose & pose );

	void
	RNA_move_trial( core::pose::Pose & pose );

	void
	random_fragment_trial( core::pose::Pose & pose );

	void
	random_domain_move_trial( core::pose::Pose & pose );

	void
	initialize_constraints( core::pose::Pose & pose );

	Size
	figure_out_constraint_separation_cutoff( Size const & r, Size const & rounds, Size const & max_dist );

	void
	update_pose_constraints( Size const & r, Size const & rounds, core::pose::Pose & pose );

	void
	update_domain_rot_trans_mag( Size const & r, Size const & rounds );

	void
	fill_pairing_dists( core::pose::Pose & pose );

	void
	check_new_pairing_dists( core::pose::Pose & pose, Size const & frag_pos );

private:

	// protocol-specific data ... need to be specified as input.
	Size const nstruct_;
	Size const monte_carlo_cycles_;
	Size const rounds_;
	std::string const silent_file_;
	bool freeze_domains_;
	bool dump_pdb_;
	core::Real domain_move_frequency_;

	// parameters
	core::Real m_Temperature_; // default temperature for monte carlo
	bool sim_anneal_;

	core::scoring::constraints::ConstraintSetOP constraint_set_;
	bool staged_constraints_;

	Size frag_size_;

	protocols::moves::MonteCarloOP monte_carlo_;

	std::string rna_params_file_;
	std::string rna_data_file_;
	std::string all_rna_fragments_file_;
	std::string jump_library_file_;

	std::map< std::string, bool > tag_is_done_;

	std::string lores_scorefxn_;
	core::scoring::ScoreFunctionOP denovo_scorefxn_;

	protocols::farna::RNA_FragmentMoverOP frag_mover_;
	protocols::farna::RNA_DeNovoPoseSetupOP rna_structure_parameters_;
	core::io::rna::RNA_DataReaderOP rna_data_reader_;
	protocols::farna::RNA_ChunkLibraryOP rna_chunk_library_;
	protocols::coarse_rna::CoarseRNA_LoopCloserOP rna_loop_closer_;
	protocols::farna::MultipleDomainMoverOP multiple_domain_mover_;

	utility::vector1< std::string > chunk_silent_files_;
	utility::vector1< core::Size > input_res_;

	bool close_loops_;
	bool choose_best_solution_;
	bool force_ideal_chainbreak_;
	bool add_base_pair_constraints_;
	bool check_pairing_dists_;
	bool view_monte_carlo_;

	utility::vector1< core::Real > pairing_dists_;

}; // class CoarseRNA_DeNovoProtocol


}
} // protocols

#endif
