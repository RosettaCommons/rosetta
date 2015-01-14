// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/mpi_refinement/MPI_Refinement.hh
/// @brief
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_mpi_refinement_MPI_Refinement_hh
#define INCLUDED_protocols_mpi_refinement_MPI_Refinement_hh

#include <protocols/wum/SilentStructStore.hh>
#include <protocols/wum/MPI_WorkUnitManager.hh>
#include <protocols/mpi_refinement/MultiObjective.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/io/silent/ProteinSilentStruct.hh>
#include <string>
#include <vector>

// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <core/kinematics/Jump.hh>

namespace protocols {
namespace mpi_refinement {

class MPI_Refinement: public protocols::wum::MPI_WorkUnitManager {
  public:
    MPI_Refinement( char machine_letter );

		void set_defaults();

    virtual ~MPI_Refinement(){}

	protected: // overloaded functions

	  // --none-- this is a pure virtual class

	protected: // added functions

		void load_structures_from_cmdline_into_library(
			 protocols::wum::SilentStructStore &library );

		void save_state(std::string prefix = "default" );

		void save_state_auto();

		void load_state(std::string prefix = "default" );

		void print_stats();

	void print_summary( std::string const prefix = "SUMM " );

	void print_library(protocols::wum::SilentStructStore &library,
										 std::string const prefix = "LIB " );

		// adding arriving structures to library
		virtual bool add_structures_to_library( protocols::wum::SilentStructStore &new_structs, std::string add_algorithm = "" );

		virtual bool add_structure_to_library( core::io::silent::SilentStructOP pss, std::string add_algorithm = "" );

		bool add_structure_to_library_direct( core::io::silent::SilentStruct &pss );

		bool add_structure_to_library_add_n_replace( core::io::silent::SilentStruct &pss );

		bool add_structure_to_library_single_replace( core::io::silent::SilentStruct &pss );

	  void update_library_NSGAII( protocols::wum::SilentStructStore &new_structs );

	  void setup_multi_objective();

	  void send_sortedpick_library_structs( core::Size dest_rank, 
																					core::Size nsend,
																					std::string const scorename,
																					bool const weighted,
																					bool const inverse = false,
																					core::Real const kT = 0.2 );

	  void send_random_library_structs( core::Size dest_rank, core::Size nsend );

		void send_random_library_struct( core::Size dest_rank, core::Size ssid ) const ;

		void limit_library();

  	void shave_library( protocols::wum::SilentStructStore &new_structs,
												std::string const scorename,
												core::Real const frac ) const;

		void dump_structures( const protocols::wum::SilentStructStore &new_structs, 
													bool score_only = true,
													std::string const prefix = "" ) const;

		void set_ident_string( std::string new_ident ){ ident_string_ = new_ident; }

	  void report_time( ) const;

	protected: // accesors
		const std::string &mpi_resume(){ return mpi_resume_; }

		core::Size & totaltime_loophash(){ return totaltime_loophash_; }

		protocols::wum::SilentStructStore  &library_central(){ return library_central_;}
		protocols::wum::SilentStructStore  &library_ref(){ return library_ref_;}


		const std::string & mpi_feedback( ){ return mpi_feedback_; }

		void set_mpi_feedback( const std::string &mpi_feedback){ mpi_feedback_ = mpi_feedback; }


		core::Size  max_lib_size(){ return max_lib_size_; }
		core::Size  max_ref_lib_size(){ return max_ref_lib_size_; }

		void set_max_lib_size( core::Size max_lib_size){ max_lib_size_ = max_lib_size; }
		void set_max_ref_lib_size( core::Size max_ref_lib_size){ max_ref_lib_size_ = max_ref_lib_size; }

		core::Real score( const core::io::silent::SilentStructOP &ss ) const;

	  void retag_library( protocols::wum::SilentStructStore &store,
												std::string const prefix ) const;

		std::string format_silent_struct( const core::io::silent::SilentStructOP &ss ) const;

		core::Real score( const core::io::silent::SilentStruct &ss ) const;

		std::string format_silent_struct( const core::io::silent::SilentStruct &ss ) const;

protected:
	// Multiobjective
	MultiObjectiveOP fobj_;

private:

	// parameters
	core::Size  max_lib_size_;
	core::Size  max_ref_lib_size_;
	core::Size  save_state_interval_;
	std::string mpi_feedback_;
	core::Real  mpi_metropolis_temp_;
	core::Real  rms_limit_;
	std::string objective_function_;

	std::string mpi_resume_;
	std::string jobname_;

	// Live data
	protocols::wum::SilentStructStore library_central_;
	// ref_structures used for generating structures; seperate it from curretly generated guys
	protocols::wum::SilentStructStore library_ref_;
	core::Size last_save_state_;

protected: // statistics can be protected such that derived classes can modify their values.
		// Statistics
	  core::Size starttime_;
		core::Size totaltime_loophash_;
		core::Size n_loophash_;
		core::Size totaltime_batchrelax_;
		core::Size n_batchrelax_;
		core::Size total_structures_;
		core::Size total_structures_relax_;
		core::Size total_metropolis_;
		core::Size total_metropolis_accepts_;

	// Native 
	bool native_given_;
	core::pose::Pose native_pose_;

private:
		std::string ident_string_; // Unique identified for this job.
  	std::string sim_replace_obj_;
};




} // namespace loops
} // namespace protocols



#endif
