// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loophash/MPI_LoopHashRefine_Master.hh
/// @brief
/// @author Mike Tyka

#ifndef INCLUDED_protocols_loophash_MPI_LoopHashRefine_Master_hh
#define INCLUDED_protocols_loophash_MPI_LoopHashRefine_Master_hh

#include <protocols/wum/SilentStructStore.hh>
#include <protocols/wum/MPI_WorkUnitManager.hh>
#include <protocols/loophash/MPI_LoopHashRefine.hh>

#include <core/types.hh>
#include <utility/vector1.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <string>
#include <vector>

namespace protocols {
namespace loophash {


class MPI_LoopHashRefine_Master: public MPI_LoopHashRefine {
public:
	MPI_LoopHashRefine_Master( core::Size my_emperor,
	                           core::Size master_rank ):
		MPI_LoopHashRefine( 'M' ),
		my_emperor_( my_emperor ),
		master_rank_( master_rank )
		{
			set_defaults();
		}

	virtual ~MPI_LoopHashRefine_Master(){};

	void set_defaults();

public:
	virtual void go();

protected: // overloaded functions

	virtual void init();

	virtual void process_inbound_wus();

	virtual void process_outbound_wus();

protected: // Added functions

	void create_loophash_WUs( const core::io::silent::SilentStructOP &start_struct );

	void add_relax_batch( protocols::wum::SilentStructStore &start_decoys );

	void check_library_expiry_dates();

	void load_sample_weight();

	virtual bool add_structure_to_library( core::io::silent::SilentStruct &pss, std::string add_algorithm = "" );

	void report_structure_to_emperor(  core::io::silent::SilentStructOP &ss ) ;

	void report_structure_to_emperor(  core::io::silent::SilentStruct &pss ) ;

	core::Real ev_objective_function( core::io::silent::SilentStructOP &ss );

protected: // Accesors

	core::Size master_rank(){ return master_rank_; }

	core::Size my_emperor(){ return my_emperor_; }
private:

	// parameters
	core::Size max_loophash_per_structure_;
	core::Size batch_relax_chunks_;
	core::Size batch_relax_absolute_max_;
	core::Size outbound_wu_buffer_size_;
	core::Size loophash_split_size_;
	core::Size library_expiry_time_;
	core::Size expire_after_rounds_;
	bool       mpi_master_save_score_only_;
	// local store
	protocols::wum::SilentStructStore  to_be_relaxed_;

	// sampling weight (is not constant!)
	std::string sample_weight_str_;

	// static settings
	const core::Size my_emperor_;
	const core::Size master_rank_;

	std::string objective_function;
};


} // namespace loops
} // namespace protocols


#endif


