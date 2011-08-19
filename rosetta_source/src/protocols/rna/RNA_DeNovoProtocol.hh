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


#ifndef INCLUDED_protocols_rna_RNA_DeNovoProtocol_hh
#define INCLUDED_protocols_rna_RNA_DeNovoProtocol_hh

#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/rna/RNA_DataReader.fwd.hh>
#include <protocols/rna/RNA_FragmentsClasses.hh>
#include <protocols/rna/RNA_StructureParameters.fwd.hh>
#include <protocols/rna/RNA_StructureParameters.hh>
#include <protocols/rna/RNA_LoopCloser.fwd.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/rna/RNA_FragmentMover.fwd.hh>
#include <protocols/rna/RNA_Minimizer.fwd.hh>
#include <protocols/rna/RNA_Relaxer.fwd.hh>
#include <protocols/rna/RNA_ChunkLibrary.fwd.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/scoring/ScoreFunction.hh>

//Oooh.
#include <ObjexxFCL/FArray1D.hh>

//// C++ headers
#include <cstdlib>
#include <string>

//Auto Headers
#include <core/io/silent/silent.fwd.hh>


namespace protocols {
namespace rna {

/// @brief The RNA de novo structure modeling protocol
class RNA_DeNovoProtocol: public protocols::moves::Mover {
public:

	/// @brief Construct the protocol object given
	/// the RNA fragment library to use.
	RNA_DeNovoProtocol(
		Size const nstruct,
		Size const monte_carlo_cycles,
		std::string const silent_file,
		bool const heat_structure = true,
		bool const minimize_structure = false,
		bool const relax_structure = false );

	~RNA_DeNovoProtocol();

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const;

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	void
	set_dump_pdb( bool const setting ){ dump_pdb_ = setting; };

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
	set_chunk_silent_files( utility::vector1< std::string > const & chunk_silent_files ) {
		chunk_silent_files_ = chunk_silent_files;
	}

	void
	ignore_secstruct( bool const setting ) { ignore_secstruct_ = setting; }

	// No longer works -- need to specify allow_insert from a "params file"
	//	void
	//	set_allow_insert( FArray1D <bool> const & allow_insert  ){ allow_insert_ = allow_insert; }

	void
	jump_change_frequency( core::Real const value ){ jump_change_frequency_ = value; }

	void
	close_loops( bool const setting ){
		close_loops_ = setting;
		if ( close_loops_ ) binary_rna_output_ = true;
	}

	void
	simple_rmsd_cutoff_relax( bool const setting ){ simple_rmsd_cutoff_relax_ = setting; }

	void
	output_lores_silent_file( bool const setting ){ output_lores_silent_file_ = setting; }

	void
	set_filter_lores_base_pairs( bool const setting ){ filter_lores_base_pairs_ = setting; }

	void
	set_binary_rna_output( bool const setting ){ binary_rna_output_ = setting; }

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
	set_staged_constraints( bool const setting ){ staged_constraints_ = setting; }

private:

	void
	initialize_movers( core::pose::Pose & pose );

	void
	initialize_lores_silent_file();

	void
	initialize_constraints( core::pose::Pose & pose );

	void
	initialize_tag_is_done();


	void
	output_silent_struct(
		core::io::silent::SilentStruct & s,
		core::io::silent::SilentFileData & silent_file_data,
		std::string const & silent_file,
		core::pose::Pose & pose,
		std::string const out_file_tag,
		bool const score_only = false
	) const;

	void
	do_random_fragment_insertions( core::pose::Pose & pose );

	void
	update_denovo_scorefxn_weights( Size const & r, Size const & rounds );

	Size
	figure_out_constraint_separation_cutoff( Size const & r, Size const & rounds, Size const & max_dist );

	void
	update_pose_constraints( Size const & r, Size const & rounds, core::pose::Pose & pose );

	void
	update_frag_size( Size const & r, Size const & rounds );

	void
	random_fragment_trial( core::pose::Pose & pose );

	void
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

private:

	// protocol-specific data ... need to be specified as input.
	Size const nstruct_;
	Size const monte_carlo_cycles_;
	std::string all_rna_fragments_file_;
	std::string const silent_file_;
	std::string lores_silent_file_;
	bool const heat_structure_;
	bool dump_pdb_;
	bool const minimize_structure_;
	bool const relax_structure_;
	bool ignore_secstruct_;
	bool close_loops_;
	bool close_loops_after_each_move_;
	bool simple_rmsd_cutoff_relax_;

	// parameters
	core::Real m_Temperature_; // default temperature for monte carlo
	Size frag_size_;

	protocols::moves::MonteCarloOP monte_carlo_;

	protocols::rna::RNA_FragmentsOP all_rna_fragments_;
	protocols::rna::RNA_FragmentMoverOP rna_fragment_mover_;
	protocols::rna::RNA_LoopCloserOP rna_loop_closer_;
	protocols::rna::RNA_ChunkLibraryOP rna_chunk_library_;
	protocols::rna::RNA_MinimizerOP rna_minimizer_;
	protocols::rna::RNA_RelaxerOP rna_relaxer_;

	std::string rna_params_file_;
	std::string rna_data_file_;
	std::string jump_library_file_;
	// This object will actually hold the jump library...
	RNA_StructureParametersOP rna_structure_parameters_;
	RNA_DataReaderOP rna_data_reader_;

	bool output_lores_silent_file_;
	bool filter_lores_base_pairs_;
	bool vary_bond_geometry_;
	bool binary_rna_output_;

	core::scoring::constraints::ConstraintSetOP constraint_set_;

	core::Real jump_change_frequency_;

	std::map< std::string, bool > tag_is_done_;

	std::string lores_scorefxn_;
	core::scoring::ScoreFunctionOP denovo_scorefxn_;
	core::scoring::ScoreFunctionOP initial_denovo_scorefxn_;
	core::scoring::rna::RNA_LowResolutionPotential local_rna_low_resolution_potential_;

	utility::vector1< std::string > chunk_silent_files_;
	core::Real chunk_coverage_;

	bool staged_constraints_;

}; // class RNA_DeNovoProtocol



}
} // protocols

#endif
