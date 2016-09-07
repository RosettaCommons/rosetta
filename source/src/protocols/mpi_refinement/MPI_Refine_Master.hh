// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/mpi_refinement/MPI_Refine_Master.hh
/// @brief
/// @author Mike Tyka
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_mpi_refinement_MPI_Refine_Master_hh
#define INCLUDED_protocols_mpi_refinement_MPI_Refine_Master_hh

#include <protocols/wum/SilentStructStore.hh>
#include <protocols/wum/MPI_WorkUnitManager.hh>
#include <protocols/mpi_refinement/MPI_Refinement.hh>
#include <protocols/mpi_refinement/Scheduler.hh>

#include <core/types.hh>
#include <utility/vector1.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <string>
#include <vector>

namespace protocols {
namespace mpi_refinement {

class MPI_Refine_Master: public MPI_Refinement {
public:
	MPI_Refine_Master( core::Size my_emperor,
		core::Size master_rank ):
		MPI_Refinement( 'M' ),
		my_emperor_( my_emperor ),
		master_rank_( master_rank )
	{
		set_defaults();
	}

	~MPI_Refine_Master() override= default;

	void set_defaults();

public:
	void go() override;

protected: // overloaded functions

	void init() override;

	void process_inbound_wus() override;

	void process_outbound_wus() override;

	// not an overloaded, but for clarity - same level with process_inbound/outbound
	bool process_round();

	bool process_termination();

protected: // Added functions

	void create_WUs( const core::io::silent::SilentStructOP &start_struct,
		core::Size const i_ss );

	void assign_loop_info( core::io::silent::SilentStructOP ss ) const;

	// Simpler re-relax
	void add_relax_simple( protocols::wum::SilentStructStore &start_decoys,
		core::Size const rerelax_type );

	void check_library_expiry_dates();

	void load_sample_weight();

	bool add_structure_to_library( core::io::silent::SilentStructOP ss, std::string add_algorithm = "" ) override;

	void feedback_structures_to_emperor( bool get_feedback,
		std::string const pick_strategy,
		std::string const objfunction );

	void feedback_structure_to_emperor( core::io::silent::SilentStructOP &ss ) ;

	void feedback_structure_to_emperor(  core::io::silent::SilentStruct &pss ) ;

	core::Size n_to_gen() const { return scheduler_.n_to_gen(); }

protected: // Accesors

	core::Size master_rank(){ return master_rank_; }

	core::Size my_emperor(){ return my_emperor_; }
private:

	core::pose::Pose
	get_average_structure( protocols::wum::SilentStructStore &decoys,
		utility::vector1< core::Size > const touse,
		std::string const columnname,
		bool const minimize,
		bool const calcdev = false) const;

private:
	// parameters
	core::Size max_sample_per_structure_;
	core::Size batch_relax_chunks_;
	core::Size batch_relax_absolute_max_;
	core::Size outbound_wu_buffer_size_;
	core::Size loophash_split_size_;
	core::Size library_expiry_time_;
	core::Size expire_after_rounds_;
	bool       mpi_master_save_score_only_;
	std::string loophash_scan_type_;
	//core::Size n_outbound_seek_next_;
	//core::Size n_outbound_round_;
	core::Size ngen_serial_per_structure_; // default ngen when not given by scheduler
	core::Real prob_terminus_ramapert_;

	//core::Size n_rerelaxed_; //
	//core::Size n_to_rerelax_; //

	// Slave control
	utility::vector1< int > myslaves_;
	std::map< int, bool > slaves_finished_;

	// Scheduler
	Scheduler scheduler_;
	bool sent_termination_;
	bool got_termination_signal_;
	bool asked_for_feedback_;
	core::Size sch_stage_;

	// sampling weight (is not constant!)
	std::string sample_weight_str_;

	// static settings
	const core::Size my_emperor_;
	const core::Size master_rank_;

	// communication with Emperor
	core::Size n_lib_modified_;
	core::Size nchange_emperor_call_;

	// structure id based on run name
	std::map< std::string, core::Size > ssids_by_name_;

};

} // namespace loops
} // namespace protocols


#endif
