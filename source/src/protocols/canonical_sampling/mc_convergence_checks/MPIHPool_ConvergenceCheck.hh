// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_canonical_sampling_mc_convergence_checks_MPIHPool_ConvergenceCheck_hh
#define INCLUDED_protocols_canonical_sampling_mc_convergence_checks_MPIHPool_ConvergenceCheck_hh

#include <protocols/canonical_sampling/mc_convergence_checks/MPIHPool_ConvergenceCheck.fwd.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/HierarchicalLevel.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/HierarchicalLevel.fwd.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/Pool_ConvergenceCheck.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/pool_util.hh>
#include <core/pose/Pose.hh>


#include <core/io/silent/SilentStruct.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <utility/vector1.hh>
#include <core/types.hh>

#ifdef USEMPI
#include <mpi.h>
#endif


namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {


class MPIHPool_RMSD : public Pool_RMSD {

public:

	typedef utility::vector1< core::Size > Address;

	MPIHPool_RMSD( std::string const& silent_file, core::Size levels ); //

	//~MPIHPool_RMSD();

	using Pool_RMSD::evaluate_and_add;

	core::Size evaluate_and_add(
		core::pose::Pose const& pose,
		std::string& best_decoy,
		core::Real& best_rmsd);


	void set_discovered_out( std::string const& newout );

	std::string const& get_discovered_out();

	void set_transition_threshold( core::Real threshold );

	void set_nresidues( core::Size nres );

	core::Size get_nresidues();

	void max_cache_size( core::Size max_cache );

	void finalize(); //need to put in

	bool is_in_neighborhood( Address & q_address, Address & ref_address );

	core::Real resolved_level_best_rmsd( Address& addr, utility::vector1< core::Real > & rmsd );

	bool is_new_structure( Address & address,
		utility::vector1< core::Real > & radii,
		utility::vector1< core::Real > & rmsds );

	bool is_new_structure( Address & address,
		utility::vector1< core::Real > & radii,
		core::Real & rmsds );

	core::Size find_address( Address & query_addr, utility::vector1< Address > & address_database );

	void
	print_address( Address & addr );

private:

	void send_receive_and_write_structures( bool winning_rank, core::pose::Pose const& pose );

	void write_decoys_to_hierarchy( core::io::silent::SilentFileData& sfd, core::io::silent::SilentStructOP& ss, Address& ss_addr, core::Size new_level_begins);

	void write_headers_to_hierarchy( core::io::silent::SilentStructOP& ss );


	void
	buf_to_address( Address & addr, int* addr_buf, core::Size index );

	void
	address_to_buf( Address & addr, int* addr_buf, core::Size index );

	bool
	is_my_structure();

	core::Size
	any_node_finished();

	void
	update_comm( core::Size newsize );

	void
	create_comm( int* ranks_to_include, int new_size );

	void
	initialize();

	void
	resolve_address_and_assign_tag( Address& new_addr, core::Size& new_level_start, std::string& new_candidate_tag );

	void
	receive_and_output_structures( core::io::silent::SilentFileData&, core::Size num_structures_to_write );

	void
	assign_tag( std::string const& address_tag, core::Size assigned_id_num, std::string & newtag );

	void
	assign_tag( Address& address_tag, core::Size assigned_id_num, std::string & newtag );

	//void
	//increment_pool_size( core::Size num_to_add );

	bool
	get_next_candidate(); //returns false if there are no more structures to process

	void
	receive_silent_struct_any_source( core::io::silent::SilentFileData&, core::io::silent::SilentStructOP & ss, Address& ss_addr, core::Size& new_level );

	void
	send_silent_struct_to_rank( core::io::silent::SilentFileData&, core::io::silent::SilentStructOP & ss, Address& ss_addr, core::Size& new_level_begins );

	void
	send_silent_struct_to_rank( core::io::silent::SilentFileData&, core::io::silent::SilentStructOP & ss, Address& ss_addr, core::Size& new_level_begins, core::Size rank);

	void
	prepare_send_new_coords( bool send_coords );

	void scan_output_and_setup_to_receive();

	void address_to_string( Address & address_buf, core::Size index, std::string & address_tag );

	void string_to_address( Address & address_buf, core::Size index, std::string & address_tag );

	protocols::canonical_sampling::mc_convergence_checks::HierarchicalLevel hlevel_;
	core::Size pool_size_;
	core::Size num_structures_added_;
	core::Size npes_;
	core::Size rank_;
	core::Size pool_npes_;
	core::Size pool_rank_;
	std::string new_decoys_out_;
	core::Size nresidues_;
	core::Size const nlevels_; //number of levels cannot change during the simulation
	bool tracer_visible_;
	bool first_time_writing_;

	utility::vector1< core::Real > level_radii_;
	Address current_address_;
	utility::vector1< core::Real > current_best_rmsds_;
	Address best_address_;
	std::string current_address_str_;

	DataBuffer buf_;
	core::Size current_trajectory_state_;

#ifdef USEMPI
	static MPI_Comm MPI_COMM_POOL;
#endif
};

} //mc_convergence_checks
} //moves
} //protocols

#endif
