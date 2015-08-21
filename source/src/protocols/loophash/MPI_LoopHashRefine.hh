// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loophash/MPI_LoopHashRefine.hh
/// @brief
/// @author Mike Tyka


#ifndef INCLUDED_protocols_loophash_MPI_LoopHashRefine_hh
#define INCLUDED_protocols_loophash_MPI_LoopHashRefine_hh

#include <protocols/wum/SilentStructStore.hh>
#include <protocols/wum/MPI_WorkUnitManager.hh>

#include <core/types.hh>
#include <utility/vector1.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <string>
#include <vector>

#include <core/kinematics/Jump.hh>


namespace protocols {
namespace loophash {


class MPI_LoopHashRefine: public protocols::wum::MPI_WorkUnitManager {
public:
	MPI_LoopHashRefine( char machine_letter );

	void set_defaults();

	virtual ~MPI_LoopHashRefine(){}

protected: // overloaded functions

	// --none-- this is a pure virtual class

protected: // added functions

	void load_structures_from_cmdline_into_library( core::Size structure_read_offset );

	void save_state(std::string prefix = "default" );

	void save_state_auto();

	void load_state(std::string prefix = "default" );

	void print_stats();

	void print_library();

	// adding arriving structures to library
	virtual bool add_structures_to_library( protocols::wum::SilentStructStore &new_structs, std::string add_algorithm = "" );

	virtual bool add_structure_to_library( core::io::silent::SilentStruct &pss, std::string add_algorithm = "" );

	bool add_structure_to_library_direct( core::io::silent::SilentStruct &pss );

	bool add_structure_to_library_add_n_replace( core::io::silent::SilentStruct &pss );

	bool add_structure_to_library_single_replace( core::io::silent::SilentStruct &pss );

	void send_random_library_struct( core::Size dest_rank, core::Size ssid ) const ;

	void limit_library();

	void dump_structures( const protocols::wum::SilentStructStore &new_structs, bool score_only = true ) const;

	void set_ident_string( std::string new_ident ){ ident_string_ = new_ident; }

protected: // accesors
	const std::string &mpi_resume(){ return mpi_resume_; }

	core::Size & totaltime_loophash(){ return totaltime_loophash_; }

	protocols::wum::SilentStructStore  &library_central(){ return library_central_;}


	const std::string & mpi_feedback( ){ return mpi_feedback_; }

	void set_mpi_feedback( const std::string &mpi_feedback){ mpi_feedback_ = mpi_feedback; }


	core::Size  max_lib_size(){ return max_lib_size_; }

	void set_max_lib_size( core::Size max_lib_size){ max_lib_size_ = max_lib_size; }

	core::Real objective_function( const core::io::silent::SilentStructOP &ss ) const;

	core::Real score( const core::io::silent::SilentStructOP &ss ) const;

	std::string format_silent_struct( const core::io::silent::SilentStructOP &ss ) const;


	core::Real objective_function( const core::io::silent::SilentStruct &ss ) const;

	core::Real score( const core::io::silent::SilentStruct &ss ) const;

	std::string format_silent_struct( const core::io::silent::SilentStruct &ss ) const;

private:
	// parameters
	core::Size  max_lib_size_;
	core::Size  save_state_interval_;
	std::string mpi_feedback_;
	core::Real  mpi_metropolis_temp_;
	core::Real  rms_limit_;
	std::string objective_function_;


	std::string mpi_resume_;
	std::string jobname_;

	// Live data
	protocols::wum::SilentStructStore  library_central_;
	core::Size last_save_state_;

protected: // statistics can be protected such that derived classes can modify their values.
	// Statistics
	core::Size totaltime_loophash_;
	core::Size n_loophash_;
	core::Size totaltime_batchrelax_;
	core::Size n_batchrelax_;
	core::Size total_structures_;
	core::Size total_structures_relax_;
	core::Size total_metropolis_;
	core::Size total_metropolis_accepts_;

private:
	std::string ident_string_; // Unique identified for this job.
};


} // namespace loops
} // namespace protocols


#endif


