// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_canonical_sampling_mc_convergence_checks_HierarchicalLevel_hh
#define INCLUDED_protocols_canonical_sampling_mc_convergence_checks_HierarchicalLevel_hh

#include <core/pose/Pose.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/HierarchicalLevel.fwd.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/Pool_ConvergenceCheck.hh>
#include <core/types.hh>
#include <core/io/silent/SilentStruct.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <utility/vector1.hh>
#include <list>

//Auto Headers
namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {

typedef ObjexxFCL::FArray2D<double> FArray2D_double;
typedef utility::vector1< core::Size > Address;

struct PoolData {
public:
	PoolData();
	PoolData( Pool_RMSD_OP pool, utility::vector1<core::Size> & address, core::Size nlevels );
	PoolData( Pool_RMSD_OP pool, utility::vector1<core::Size> & address);

	void
	setup( core::Size nlevels );

	Pool_RMSD_OP pool_;
	utility::vector1<core::Size> address_;
	core::Size nlevels_;

};


class HierarchicalLevel : public Pool_RMSD {

public:

	HierarchicalLevel(core::Size nlevels, std::string const & filename);

	HierarchicalLevel( core::Size nlevels );

	HierarchicalLevel(core::Size level, core::Size max_levels );

	std::string&
	filename();

	//  void
	//  union_of_addresses( Address& addr1, Address& addr2, Address& union_addr );

	void
	fill_top_level( Pool_RMSD_OP & top_level_pool );

	void
	get_addresses_at_level(  utility::vector1< Address > & );

	void
	debug_print_hierarchy();

	void
	debug_print_size_per_level();

	void
	print_addresses_at_level();

	HierarchicalLevelOP
	next_level();

	core::Size
	size();

	core::Size
	level_size( core::Size level );

	core::Size
	pool_size( Address & address, core::Size level );

	core::Size
	top_level_pool_size();

	Address &
	address( core::Size index ); //returns address of pool at index

	void
	radius( core::Real radius );

	core::Real
	radius();

	core::Real
	level_radius( core::Size level );

	core::Size
	level();

	core::Size
	nlevels();

	void
	max_cache_size( core::Size max_size );

	core::Size
	max_cache_size();

	bool
	has_next_level();

	void
	clear();

	void
	add_new( core::pose::Pose const& pose,
		std::string & tag,
		Address & address,
		bool write_to_file,
		core::Size new_level );

	void
	add_new( core::io::silent::SilentStruct &,
		std::string & tag,
		Address & address,
		bool write_to_file,
		core::Size new_level );

	bool
	add_new( FArray2D_double&,
		std::string & tag,
		Address & address);

	bool
	add_new( FArray2D_double&,
		std::string & tag,
		Address & address,
		bool resolve_remaining_levels );

	void
	add_to_top_level( FArray2D_double&, std::string & tag );


	bool
	equal_addresses(Address & address1,
		Address & address2);

	core::Size
	evaluate( core::pose::Pose const& pose,
		std::string & best_decoy,
		utility::vector1< core::Real > & best_rmsd,
		Address & address);

	core::Size
	evaluate( core::io::silent::SilentStruct const& ,
		std::string & best_decoy,
		utility::vector1< core::Real > & best_rmsd,
		Address & address);

	core::Size
	evaluate( FArray2D_double&,
		std::string & best_decoy,
		utility::vector1< core::Real > & best_rmsd,
		Address & address);

	core::Size
	evaluate( FArray2D_double&,
		std::string & best_decoy,
		utility::vector1< core::Real > & best_rmsd,
		Address & address,
		bool reset_all_levels,
		bool load_if_missing );


	/**
	core::Size
	evaluate( FArray2D_double&,
	std::string & best_decoy,
	core::Real & best_rmsd,
	Address & best_indices,
	utility::vector1< Address > & prev_address,
	utility::vector1< core::Size > & prev_indices );
	**/
	/**
	std::list<PoolData> &
	get_cache(); //for debugging
	**/
	// void
	//  write_headers_to_all_files_in_hierarchy();

	std::string
	lib_full_path(Address & address ); //made public for debugging

	bool
	address_exists_in_cache( Address & addr );

	bool
	level_find( Address & address, core::Size const& level, std::list<PoolData>::iterator& );

	bool
	pool_exists( Address & address ); //simply checks if appropriate filename exists

	void
	round( core::Real& rmsd );

	void
	test_round();

	// Undefined, commenting out to fix PyRosetta build
	//core::Size
	//next_free_address( Address& address, core::Size level );

	void
	write_headers_to_hierarchy( core::io::silent::SilentStructOP& ss );

	void
	debug_print_addresses();


	core::Size
	first_zero_pos( Address & address );


protected:

private:

	void
	sort_pool( Pool_RMSD_OP& pool_ptr );

	core::Size
	find_address( Address& query, utility::vector1< Address >& addresses );

	void
	string_to_address( std::string&, Address& );

	void
	address_to_string( std::string&, Address& );

	void
	sort_addresses( utility::vector1< Address >& addresses );

	bool
	less_than( Address& addr1, Address& addr2 );

	void
	max_level_index( Address&);


	void
	reset( Address & addr,
		utility::vector1< core::Real > & rms );

	void
	print_address( Address & address );


	core::Size
	num_matching_levels( Address & address1,
		Address & address2 );

	void
	add_elem_to_cache( FArray2D_double & coords, std::string & tag, Address & input_addr );

	std::list<PoolData>::iterator
	find(Address & address);

	PoolData
	load_pool( Address & address );


	/**
	std::string
	lib_full_path(std::string tag_origin );
	**/

	std::list<PoolData> pool_cache_;
	core::Size max_cache_size_;
	core::Size const level_;
	core::Size const max_levels_;
	core::Real radius_;
	core::Size num_clusters_in_level_;
	std::string filename_;
	core::Size next_free_address_index_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool first_time_writing_;

	HierarchicalLevelOP next_level_;
};

}
}
}

#endif //INCLUDED_protocols_canonical_sampling_mc_convergence_check_HierarchicalLevel_HH
