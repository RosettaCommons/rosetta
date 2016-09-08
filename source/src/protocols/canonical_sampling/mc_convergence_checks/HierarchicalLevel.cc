// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <protocols/canonical_sampling/mc_convergence_checks/HierarchicalLevel.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/HierarchicalLevel.fwd.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/Pool_ConvergenceCheck.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>


#include <ObjexxFCL/FArray2D.hh>
#include <protocols/toolbox/DecoySetEvaluation.hh>
#include <protocols/toolbox/superimpose.hh>

#include <utility/file/FileName.hh>
#include <utility/file/PathName.hh>
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/prof.hh>

#include <set>
#include <list>
#include <sstream>
#include <fstream>
#include <iostream>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mc.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/Jump.hh>

//this is actually only for debugging
#ifdef USEMPI
#include <mpi.h>
#endif

namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {

static THREAD_LOCAL basic::Tracer TR( "HierarchicalLevel" );
//  core::Real const MAX_RADIUS = 9999;

using namespace basic;

PoolData::PoolData():
	pool_(),
	address_(),
	nlevels_(0)
{}

PoolData::PoolData( Pool_RMSD_OP pool, utility::vector1<core::Size> & address, core::Size nlevels ):
	pool_(pool),
	address_(address)
{
	setup(nlevels);
}

PoolData::PoolData( Pool_RMSD_OP pool, Address & address):
	pool_(pool),
	address_(address),
	nlevels_(0)
{}

void
PoolData::setup( core::Size nlevels ){
	nlevels_ = nlevels;
	address_.resize(nlevels,0);
}


void
HierarchicalLevel::get_addresses_at_level( utility::vector1< Address >& addresses_at_level ) {
	core::Size num_addresses_at_level = 0;
	for ( std::list< PoolData >::iterator itr = pool_cache_.begin(), end = pool_cache_.end();
			itr != end; ++itr ) {
		Address addr = (*itr).address_;
		addresses_at_level.push_back( addr );
		num_addresses_at_level++;
	}
}

HierarchicalLevel::HierarchicalLevel( core::Size maxlevels ):
	pool_cache_(), //cache is filled as needed
	max_cache_size_( 100 ),
	level_(1),
	max_levels_( maxlevels ),
	radius_(),
	num_clusters_in_level_(0),
	next_free_address_index_(0)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1< core::Size >  universal_address( max_levels_ , 0 );
	universal_address[ 1 ] = 1;
	Pool_RMSD_OP poolrms( new Pool_RMSD() );
	pool_cache_.push_back( PoolData( poolrms, universal_address ) );
	num_clusters_in_level_++;
	radius_ = option[ cluster::K_radius ]()[ level_ ];

	if ( TR.visible() ) { TR.Debug << "new hierarchical-level created, level=" << level_ << " radius=" << radius_ << std::endl;}
	if ( ( level_ + 1 ) <= max_levels_ ) {
		next_level_ = HierarchicalLevelOP( new HierarchicalLevel( (core::Size)(level_ + 1), max_levels_ ) );
	} else {
		next_level_.reset();
	}
}

//filename assumed to specify structures at top-most level
HierarchicalLevel::HierarchicalLevel( core::Size maxlevels, std::string filename ):
	pool_cache_(), //cache is filled as needed
	max_cache_size_( 100 ),
	level_(1),
	max_levels_( maxlevels ),
	radius_(),
	num_clusters_in_level_(0),
	filename_(filename),
	next_free_address_index_(0)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1< core::Size >  universal_address( max_levels_ , 0 );
	universal_address[ 1 ] = 1;
	Pool_RMSD_OP poolrms( new Pool_RMSD() );
	pool_cache_.push_back( PoolData( poolrms, universal_address ) );
	num_clusters_in_level_++;
	radius_ = option[ cluster::K_radius ]()[ level_ ];

	if ( TR.visible() ) { TR.Debug << "new hierarchical-level created, level=" << level_ << " radius=" << radius_ << std::endl; }
	if ( ( level_ + 1 ) <= max_levels_ ) {
		next_level_ = HierarchicalLevelOP( new HierarchicalLevel( (core::Size)(level_ + 1), max_levels_ ) );
	} else {
		next_level_.reset();
	}
}

HierarchicalLevel::HierarchicalLevel( core::Size level, core::Size max_levels ):
	pool_cache_(), //cache is filled as needed
	max_cache_size_( 100 ),
	level_(level),
	max_levels_(max_levels),
	radius_(0),
	num_clusters_in_level_(0),
	next_free_address_index_(0)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	radius_ = option[ cluster::K_radius ]()[ level_ ];
	if ( TR.visible() ) { TR.Debug << "new hierarchical-level created level=" << level_ << " radius=" << radius_ << std::endl; }
	if ( level_ < max_levels_ ) {
		next_level_ = HierarchicalLevelOP( new HierarchicalLevel( level_ + 1, max_levels ) );
	} else {
		next_level_.reset();
	}
}

std::string &
HierarchicalLevel::filename() {
	return filename_;
}

HierarchicalLevelOP
HierarchicalLevel::next_level() {
	return next_level_;
}

void
HierarchicalLevel::radius( core::Real radius ) {
	radius_ = radius;
}

core::Real
HierarchicalLevel::radius(){
	return radius_;
}

core::Real
HierarchicalLevel::level_radius( core::Size level ) {
	core::Real level_radius = 0;
	if ( level_ != level && has_next_level() ) {
		level_radius = next_level_->level_radius( level );
	} else {
		runtime_assert( level_ == level );
		level_radius = radius_;
	}
	return level_radius;
}

core::Size
HierarchicalLevel::level(){
	return level_;
}

core::Size
HierarchicalLevel::nlevels() {
	return max_levels_;
}

void
HierarchicalLevel::max_cache_size( core::Size max_size ) {
	max_cache_size_ = max_size;
	if ( has_next_level() ) {
		next_level_->max_cache_size( max_size );
	}
}

core::Size
HierarchicalLevel::max_cache_size() {
	return max_cache_size_;
}


bool
HierarchicalLevel::has_next_level(){
	return ( next_level_ != 0 );
}

//if pool exists on file but not in memory, will not load pool from file
std::list<PoolData>::iterator
HierarchicalLevel::find( utility::vector1< core::Size > & address ) {
	PROF_START( basic::HIERARCHICAL_FIND );
	if ( pool_cache_.size() == 0 ) {
		return pool_cache_.end();
	}
	std::list<PoolData>::iterator itr = pool_cache_.begin();
	while ( itr != pool_cache_.end() ) {
		if ( equal_addresses( address, (*itr).address_ ) ) {
			break;
		}
		++itr;
	}
	// This PROF_STOP used to be following the return!
	PROF_STOP( basic::HIERARCHICAL_FIND );
	return itr;
}

bool
HierarchicalLevel::level_find( Address & addr, core::Size const& level, std::list<PoolData>::iterator& itr ) { //find a particular address at a particular level
	PROF_START( basic::HIERARCHICAL_FIND );
	if ( pool_cache_.size() == 0  && level == level_ ) {
		itr = pool_cache_.end();
		return false;
	}

	if ( level_ != level  && has_next_level() ) {
		bool found_addr = next_level_->level_find( addr, level, itr );
		return found_addr;
	} else if ( level_ == level ) {
		if ( TR.visible() ) {
			TR.Debug << "searching for address: ";
			for ( core::Size ii = 1; ii <= addr.size(); ii++ ) {
				TR.Debug  << addr[ ii ] << " ";
			}
			TR.Debug << " at level: " << level << std::endl;
		}
		runtime_assert( level_ == level );
		runtime_assert( first_zero_pos( addr ) > level ); //not looking for an address that is unresolved
		itr = find( addr );
		//ek added 7/29/10
		bool load_if_missing = false; //never load if missing!
		if ( load_if_missing && itr == pool_cache_.end() && pool_exists( addr ) ) {  //pool exists on file, but not in memory
			PoolData pd = load_pool( addr ); //load pool into cache and return itr to pool
			pool_cache_.push_front( pd );
			itr = pool_cache_.begin();
			Address tmp = (*itr).address_;
			if ( TR.visible() ) {
				TR.Debug << "found address at level: " << level << " : ";
				for ( core::Size ii = 1; ii <= tmp.size(); ii++ ) {
					TR.Debug << tmp[ ii] << " ";
				}
				TR.Debug << " size of pool: " << (*itr).pool_->size() << std::endl;
			}
			return true;
			//ek added 7/29/10
		} else if ( itr == pool_cache_.end() ) { //if pool doesn't exists on file
			return false;
		} else {
			Address tmp = (*itr).address_;
			if ( TR.visible() ) {
				TR.Debug << "found address at level: " << level << " : ";
				for ( core::Size ii = 1; ii <= tmp.size(); ii++ ) {
					TR.Debug << tmp[ ii] << " ";
				}
				TR.Debug << " size of pool: " << (*itr).pool_->size() << std::endl;
			}
			return true;
		}
	}
	TR.Error << "you got to a point in level_find that you shouldn't have gotten to" << std::endl;
	PROF_STOP( basic::HIERARCHICAL_FIND );
	return false;
}


core::Size
HierarchicalLevel::num_matching_levels( utility::vector1< core::Size > & address1,
	utility::vector1< core::Size > & address2 ){
	runtime_assert( address1.size() == address2.size() );
	core::Size num_matching_levels = 0;
	for ( core::Size ii = 1; ii < address1.size(); ii++ ) {
		if ( address1[ ii ] == address2[ ii ] ) {
			num_matching_levels++;
		}
	}
	return num_matching_levels;
}


bool
HierarchicalLevel::equal_addresses( utility::vector1<core::Size> & address1,
	utility::vector1<core::Size> & address2) {
	//runtime_assert( address1.size() == address2.size() );
	if ( address1.size() != address2.size() ) return false;
	bool is_same = true;
	for ( core::Size ii = 1; ii <= address1.size(); ii++ ) {
		if ( address1[ ii ] != address2[ ii ] && address1[ ii ] != 0 && address2[ ii ] != 0 ) {
			is_same = false;
		}
	}
	return is_same;
}


void
HierarchicalLevel::add_new( core::pose::Pose const& pose,
	std::string & tag,
	utility::vector1< core::Size > & address,
	bool write_to_file,
	core::Size new_level) {
	ObjexxFCL::FArray2D_double coord( 3, pose.size(), 0.0 );
	protocols::toolbox::fill_CA_coords( pose, coord );
	utility::vector1< core::Size > tmp_address = address;
	if ( new_level == 0 ) {
		new_level = tmp_address.size() + 1; //no new levels
	} else {
		for ( core::Size ii = new_level; ii <= tmp_address.size(); ii++ ) {
			tmp_address[ ii ] = 0;
		} //reset so we can write to all new-levels
	}
	add_new( coord, tag, address );
	if ( write_to_file ) {
		if ( TR.visible() ) { TR.Debug << "does this file exist? " << lib_full_path(tmp_address) << " " << utility::file::file_exists(lib_full_path( tmp_address )) << std::endl; }
		core::io::silent::SilentFileData sfd;
		sfd.strict_column_mode( true );
		core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
		ss->fill_struct(pose, tag);

		for ( core::Size ii = ( new_level - 1 ); ii <= tmp_address.size(); ii++ ) { //need to add to very last level as well
			std::string silent_outfile = lib_full_path( tmp_address );
			utility::file::FileName file( silent_outfile );
			if ( !utility::file::file_exists(silent_outfile) ) {
				if ( TR.visible() ) { TR.Debug << "this file " << silent_outfile << " does not exist. attempting to create" << std::endl; }
				utility::file::create_directory_recursive( file.path() );
				utility::file::create_blank_file( file.name() );
				std::ofstream os;
				os.open( (file.name()).c_str() );
				ss->print_header( os ); //temporary fix
				os.close();
			} else {
				std::ofstream os;
				os.open( (file.name()).c_str(), std::ios::app );
				ss->print_header( os ); //temporary fix
				os.close();
			}
			//runtime_assert( utility::file::file_exists(silent_outfile) ) ;
			if ( TR.visible() ) { TR.Debug << "writing " << ss->decoy_tag() << " to " << silent_outfile << std::endl; }
			std::ofstream os;
			os.open( (silent_outfile).c_str(), std::ios::app );
			sfd._write_silent_struct( (*ss), os, false );
			os.close();
			if ( ii < tmp_address.size() ) {
				tmp_address[ ii + 1 ] = address[ ii + 1 ];
			}
		}
	}
}

void
HierarchicalLevel::add_new( core::io::silent::SilentStruct & ss,
	std::string & tag,
	utility::vector1< core::Size > & address,
	bool write_to_file,
	core::Size new_level ) {
	ObjexxFCL::FArray2D_double coord( ss.get_CA_xyz() );
	utility::vector1< core::Size > tmp_address = address;
	if ( new_level == 0 ) {
		new_level = tmp_address.size() + 1;
	} else {
		for ( core::Size ii = new_level; ii <= tmp_address.size(); ii++ ) {
			tmp_address[ ii ] = 0;
		}
	}
	add_new( coord, tag, address );
	if ( write_to_file ) {
		//actually requires some book-keeping to keep track of which ones to write
		for ( core::Size ii = ( new_level - 1 ); ii <= address.size(); ii++ ) {
			core::io::silent::SilentFileData sfd;
			sfd.strict_column_mode( true );
			std::string silent_outfile = lib_full_path( tmp_address );
			utility::file::FileName file( silent_outfile ); //temp fix
			if ( !utility::file::file_exists(silent_outfile) ) {
				utility::file::create_directory_recursive( file.path() );
				utility::file::create_blank_file( file.name() );
				std::ofstream os;
				os.open( (file.name()).c_str() );//temp fix
				ss.print_header( os ); //temporary fix
				os.close();
			} else {
				std::ofstream os;
				os.open( (file.name()).c_str(), std::ios::app );
				ss.print_header( os ); //temporary fix
				os.close();
			}
			runtime_assert( utility::file::file_exists(silent_outfile));
			std::ofstream os;
			os.open( (silent_outfile).c_str(), std::ios::app );
			sfd._write_silent_struct( ss, os, false );
			os.close();
			if ( ii < tmp_address.size() ) {
				tmp_address[ ii + 1 ] = address[ ii + 1 ];
			}
		}
	}

}


core::Size
HierarchicalLevel::evaluate(
	core::pose::Pose const& pose,
	std::string & best_decoy,
	utility::vector1< core::Real > & best_rmsd,
	utility::vector1< core::Size > & best_address
) {
	ObjexxFCL::FArray2D_double coord( 3, pose.size(), 0.0 );
	protocols::toolbox::fill_CA_coords( pose, coord );
	core::Size last_level_address = evaluate( coord, best_decoy, best_rmsd, best_address );
	return last_level_address;
}

core::Size
HierarchicalLevel::evaluate(
	core::io::silent::SilentStruct const& ss,
	std::string & best_decoy,
	utility::vector1< core::Real > & best_rmsd,
	utility::vector1< core::Size > & best_address
) {
	ObjexxFCL::FArray2D_double coord( ss.get_CA_xyz() );
	core::Size last_level_address = evaluate( coord, best_decoy, best_rmsd, best_address );
	return last_level_address;
}

core::Size
HierarchicalLevel::first_zero_pos( Address & address ) {
	core::Size first_zero_pos = 0;
	for ( first_zero_pos = 1; first_zero_pos <= address.size(); first_zero_pos++ ) {
		if ( address[ first_zero_pos ] == 0 ) break;
	}
	return first_zero_pos;
}

void
HierarchicalLevel::max_level_index( Address & addr) {
	addr.resize( max_levels_, 0 );
	addr[ 1 ] = 1;
	for ( core::Size ii = 2; ii <= addr.size(); ii++ ) {
		addr[ ii ] = 0;
	}
	if ( level_ == 1 && has_next_level() ) {
		next_level_->max_level_index( addr );
	} else {
		if ( next_free_address_index_ > 0 ) {
			addr[ level_ ] = next_free_address_index_;
		} else {
			std::string test_file = lib_full_path( addr );
			while ( utility::file::file_exists( test_file ) ) {
				addr[ level_ ] = addr[ level_ ] + 1;
				test_file = lib_full_path( addr );
			}
			addr[ level_ ] = addr[ level_ ] - 1;
		}
	}
}

void
HierarchicalLevel::write_headers_to_hierarchy( core::io::silent::SilentStructOP& ss ) {
	Address max_levels_seen( max_levels_, 0 );
	max_level_index( max_levels_seen );
	Address path_address( max_levels_, 0 );
	for ( core::Size ii = 1; ii <= max_levels_seen.size(); ii++ ) {
		for ( core::Size jj = 1; jj <= max_levels_seen[ ii ]; jj++ ) {
			path_address[ ii ] = jj;
			std::string filename = lib_full_path( path_address );
			std::ofstream os;
			os.open( filename.c_str(), std::ios::app );
			ss->print_header( os );
			os.close();
		}
	}
}

void
HierarchicalLevel::sort_pool( Pool_RMSD_OP& pool_ptr ) {
	PROF_START( basic::SORT_POOL );
	if ( TR.visible() ) TR.Debug << "now going to sort pool" << std::endl;
	Pool_RMSD_OP sorted_pool( new Pool_RMSD() );
	utility::vector1< Address > tags;
	utility::vector1< Address > initial_order;
	tags.reserve( pool_ptr->size() );
	if ( TR.visible() ) TR.Debug << "there are " << pool_ptr->size() << " structures in the starting pool-rmsd" << std::endl;
	initial_order.reserve( pool_ptr->size() );

	for ( core::Size ii = 1; ii <= pool_ptr->size(); ii++ ) {
		Address addr;
		string_to_address( pool_ptr->get_tag( ii ), addr );
		tags.push_back( addr );
		initial_order.push_back( addr );
	}

	sort_addresses( tags );
	FArray2D_double coords;
	for ( core::Size ii = 1; ii <= pool_ptr->size(); ii++ ) {
		core::Size index = find_address( tags[ ii ], initial_order );
		runtime_assert( index <= pool_ptr->size());
		pool_ptr->get( index, coords );
		std::string newtag = pool_ptr->get_tag( index );
		if ( TR.visible() ) TR.Debug << "adding structures with tag: " << newtag << " and " << coords.u2() << " residues  index: " << index << " out of " << pool_ptr->size() << std::endl;
		sorted_pool->add( coords, coords.u2(), newtag );
	}
	pool_ptr = sorted_pool;
	PROF_STOP( basic::SORT_POOL );
}

core::Size
HierarchicalLevel::find_address( Address& query, utility::vector1< Address >& addresses ) {
	PROF_START( basic::HIERARCHICAL_FIND_ADDRESS );
	core::Size ii;
	for ( ii = 1; ii <= addresses.size(); ii++ ) {
		if ( equal_addresses( addresses[ ii ], query ) ) {
			break;
		}
	}
	// AMW: this return used to be before PROF_STOP
	PROF_STOP( basic::HIERARCHICAL_FIND_ADDRESS );
	return ii;
}

void
HierarchicalLevel::sort_addresses( utility::vector1< Address >& addresses ) { //horribly inefficient, but just for starters
	Address tmp;
	PROF_START( basic::HIERARCHICAL_SORT_ADDRESSES );
	for ( core::Size ii = 1; ii <= addresses.size(); ii++ ) {
		core::Size min_index = ii;
		Address min_index_addr = addresses[ ii ];
		for ( core::Size jj = ii + 1; jj <= addresses.size(); jj++ ) {
			if ( less_than( addresses[ jj ], min_index_addr  ) ) { // jj less than min_index
				min_index = jj;
				min_index_addr = addresses[ jj ];
			}
		}
		if ( ii != min_index ) {
			tmp = addresses[ ii ];
			addresses[ ii ] = min_index_addr;
			addresses[ min_index ] = tmp;
		}
	}
	PROF_STOP(  basic::HIERARCHICAL_SORT_ADDRESSES );
}

void
HierarchicalLevel::string_to_address( std::string& addr, Address& result_addr ) {
	core::Size pos= addr.find(".",0);
	core::Size length = addr.length();
	while ( pos < length ) {
		core::Size newpos = addr.find( ".", pos+1 );
		result_addr.push_back( atoi(addr.substr( pos+1, newpos - pos - 1).c_str() ) );
		pos = newpos;
	}
}

void
HierarchicalLevel::address_to_string(std::string& newtag, Address& addr ) {
	std::ostringstream q;
	for ( core::Size ii = 1; ii < addr.size(); ii++ ) {
		q << addr[ ii ] << ".";
	}
	q << addr[ addr.size() ];
	newtag = "new." + q.str();
}

bool
HierarchicalLevel::less_than( Address& addr1, Address& addr2 ) { //addr1 less than addr2?
	if ( addr1.size() < addr2.size() ) {
		return true; //probably doesn't have correct format, just put in front
	}
	for ( core::Size ii = 1; ii <= addr1.size(); ii++ ) {
		if ( addr1[ ii ] != addr2[ ii ] ) { //ignore if numbers are equal
			if ( addr1[ ii ] < addr2[ ii ] ) {
				return true;
			} else {
				return false;
			}
		}
	}
	return false;
}


/**
//this function not used
core::Size
HierarchicalLevel::next_free_address( Address & query_address, core::Size level ) {
runtime_assert( level != 1 );
core::Size next_index = 0;
if( level_ != level && has_next_level() ) {
next_index = next_level_->next_free_address( query_address, level );
}
else {
runtime_assert( level_ == level );
if ( next_free_address_index_ > 0 ) {
next_free_address_index_++;
next_index = next_free_address_index_;
} else {
runtime_assert( first_zero_pos( query_address ) == level ); //this level is unresolved
Address address = query_address;
address[ level_ ] = pool_cache_.size() + 1;
for( core::Size ii = level_ + 1; ii <= address.size(); ii++ ) {
address[ ii ] = 0; //just to be safe
}
std::string test_file = lib_full_path( address );
while( utility::file::file_exists( test_file ) ) {
address[ level_ ] = address[ level_ ] + 1;
test_file = lib_full_path( address );
}
next_free_address_index_ = address[ level_ ];
next_index = next_free_address_index_;
}
}
return next_index;
}
**/

void
HierarchicalLevel::reset( Address & addr,
	utility::vector1< core::Real > & rms ) {
	for ( core::Size ii =1; ii <= addr.size(); ii++ ) {
		addr[ ii ] = 0;
		rms[ ii ] = 0.0;
	}
}

core::Size
HierarchicalLevel::evaluate(
	FArray2D_double& coords,
	std::string & best_decoy,
	utility::vector1< core::Real > & best_rmsd,
	Address & best_indices
) {
	return evaluate( coords, best_decoy, best_rmsd, best_indices, true, true );//when you evaluate, start from top-level
}

core::Size
HierarchicalLevel::evaluate(
	FArray2D_double& coords,
	std::string & best_decoy,
	utility::vector1< core::Real > & best_rmsd,
	Address & best_indices,
	bool reset_all_levels, //default true,
	bool load_if_missing_from_cache //default true
) {
	using namespace basic;
	PROF_START( basic::HIERARCHICAL_EVALUATE );
	if ( level_ == 1  ) {
		best_rmsd.resize( max_levels_, 0.0 );
		best_indices.resize( max_levels_, 0 );
		if ( reset_all_levels ) {
			reset( best_indices, best_rmsd );
			best_indices[ 1 ] = 1;
		}
	}

	if ( !reset_all_levels && level_ < max_levels_ && best_indices[ level_ + 1 ] != 0 ) {
		return next_level_->evaluate( coords, best_decoy, best_rmsd, best_indices, reset_all_levels, load_if_missing_from_cache );
	}
	//std::string current_best_decoy="";
	std::list<PoolData>::iterator addr_itr  = find( best_indices );
	PoolData matching_pool;
	Pool_RMSD_OP matching_pool_ptr;
	bool element_exists_in_cache = true;

	if ( TR.visible() ) { TR.Debug << "address at time of evaluation for level " << level_ << " is: "; print_address( best_indices ); }
	if ( addr_itr == pool_cache_.end() ) { //address not in cache, but exists as file
		//      if( pool_exists( best_indices ) ) { //ek elminated un-necessary file-checks 10/8/2010
		if ( TR.visible() ) TR.Debug << "address doesn't exist in cache, but exists on file... loading..." << std::endl;
		element_exists_in_cache = false;
		if ( load_if_missing_from_cache ) {
			while ( pool_cache_.size() >= max_cache_size_ ) { //make some room
				pool_cache_.pop_back();
			}
			if ( TR.visible() ) {
				TR.Debug << "EVAL: address "; print_address(best_indices);
				TR << "loading " << lib_full_path( best_indices ) << " from file" << std::endl;
			}
			matching_pool = load_pool( best_indices );
			pool_cache_.push_front( matching_pool );
		} else {
			if ( TR.visible() ) TR.Debug << "element doesn't exist in cache and we're not allowed to load from file, exiting eval" << std::endl;
			return 0; //element doesn't exist in cache, and we're not allowed to load from file
		}
		//} //ek elminated un-necessary file-checks 10/8/2010
		/** //ek elminated un-necessary file-checks 10/8/2010
		else {
		if( TR.visible() ) TR.Debug << "pool desired from address doesn't exist! exiting eval "; print_address( best_indices );
		return 0;
		}
		**/
	} else {
		matching_pool = (*addr_itr);
	}
	matching_pool_ptr = matching_pool.pool_;
	if ( element_exists_in_cache && pool_cache_.size() > 1 ) {
		if ( TR.visible() ) { TR.Debug << "EVAL: element exists in cache, shuffling to front" << std::endl;}
		pool_cache_.erase( addr_itr ); //take PoolData and put to front
		pool_cache_.push_front( matching_pool );
	}

	//evaluation begins
	core::Real best_level_rmsd = 0.0;
	core::Size level_address;
	runtime_assert( pool_cache_.begin() == find( best_indices ));
	level_address = ((*pool_cache_.begin()).pool_)->evaluate( coords, best_decoy, best_level_rmsd );
	// AMW: cppcheck will flag this because they think it's cmath round
	round( best_level_rmsd ); //ek 8/3/2010

	core::Size last_level_address = 0;
	//tolerance = 0.05;
	if ( has_next_level() ) {
		bool outside_known_structures = false;
		if ( TR.visible() ) TR.Debug << "going down to next level" << std::endl;
		if ( best_level_rmsd > radius_  ) {
			// if( best_level_rmsd > radius_ + tolerance ) { // 10/12/10 ek added in tolerance for very close cases
			// if( best_level_rmsd > radius_ + tolerance ) { // 10/12/10 ek added in tolerance for very close cases
			if ( TR.visible() ) { TR.Debug << "best level rmsd: " << best_level_rmsd << " greater than radius: " << radius_ << std::endl; }
			outside_known_structures = true;
			best_indices[ level_ + 1 ] = 0;
			best_rmsd[ level_ ] = best_level_rmsd;
			if ( TR.visible() ) TR.Debug << "assigned values 0 and " << best_level_rmsd << " to best_indices and best_rmsd " << std::endl;
		} else {
			best_indices[ level_  + 1 ] = level_address;
			best_rmsd[ level_ ] = best_level_rmsd;
			if ( TR.visible() ) TR.Debug << "assigned values " << level_address << " and " << best_level_rmsd << " to best_indices and best_rmsd " << std::endl;
		}
		if ( !outside_known_structures ) {
			if ( TR.visible() ) { TR.Debug << "EVAL: now going one level down to evaluate recursively " << std::endl; }
			last_level_address = next_level_->evaluate( coords, best_decoy, best_rmsd, best_indices, reset_all_levels, load_if_missing_from_cache );
		}
	} else { //last-level
		if ( TR.visible() ) TR.Debug << "reached last level: assigning a value of " << level_ << " and rms: " << best_level_rmsd << std::endl;
		best_rmsd[ level_ ] = best_level_rmsd;
		return level_address;
	}
	PROF_STOP( basic::HIERARCHICAL_EVALUATE );
	return last_level_address;
}

void
HierarchicalLevel::round( core::Real& best_rmsd ) {
	PROF_START( basic::HIERARCHICAL_ROUND );
	//round to nearest 0.01
	if ( TR.visible() ) TR.Debug << "input to round is: " << best_rmsd << " ";
	core::Real t = best_rmsd * 100;
	core::Real remainder = t - floor(t);
	//TR.Debug << "remainder: " << remainder << std::endl;
	if ( remainder >= 0.5 ) {
		remainder *= 10;
		remainder = ceil(remainder);
		remainder /= 10;
		//TR.Debug << " ceil-remainder: " << remainder << std::endl;
	} else {
		remainder *= 10;
		remainder = floor(remainder);
		remainder /= 10;
		//TR.Debug << "floor-remainder: " << remainder << std::endl;
	}

	//TR.Debug << "remainder: " << remainder << " floor(t): " << floor(t) << std::endl;
	best_rmsd =  (remainder + floor(t))/100;
	//hack to cut off last digit when best_rmsd is x.0y such that 2.01 --> 2.0
	if ( (int)(floor( best_rmsd * 10.0 )) % 10 == 0 ) {
		best_rmsd = (floor(best_rmsd * 10)) / 10;
	}
	if ( TR.visible() ) TR.Debug << " output of round: " << best_rmsd << std::endl;
	PROF_STOP(  basic::HIERARCHICAL_ROUND  );
}

void
HierarchicalLevel::test_round() {
	TR.Debug << "testing ROUND:\n";
	for ( core::Real t = 0; t <= 10.0; t += 0.0001 ) {
		core::Real rounded = t;
		// AMW: cppcheck will flag this because it will think this is cmath's round
		round(rounded);
		TR.Debug << t << " " << rounded << std::endl;
	}
	TR.Debug << "end testing ROUND\n";
	exit(1);

}

core::Size
HierarchicalLevel::size(){
	return pool_cache_.size();
}


core::Size
HierarchicalLevel::level_size( core::Size level ) {
	core::Size level_size = 0;
	if ( level_ != level ) {
		runtime_assert( level <= max_levels_ );
		if ( has_next_level() ) {
			level_size = next_level_->level_size( level );
		}
	} else {
		level_size =  size();
	}
	return level_size;
}

core::Size
HierarchicalLevel::pool_size( Address & addr, core::Size level ) {
	PROF_START( basic::HIERARCHICAL_POOL_SIZE );
	//    core::Size first_zero = first_zero_pos( addr );
	//    core::Size pool_size = 0;
	std::list< PoolData >::iterator itr;
	if ( level_find( addr, level, itr ) ) {
		if ( TR.visible() ) TR.Debug << "pool size from pool_size function is:  " << (*itr).pool_->size() << std::endl;
		return (*itr).pool_->size();
	} else {
		return 0;
	}
	PROF_STOP( basic::HIERARCHICAL_POOL_SIZE );
}


core::Size
HierarchicalLevel::top_level_pool_size() {
	runtime_assert( level_ == 1 );
	runtime_assert( pool_cache_.size() == 1 );
	return (*pool_cache_.begin()).pool_->size();
}

void
HierarchicalLevel::add_elem_to_cache( FArray2D_double & coords, std::string & tag, Address & input_addr ) {
	PROF_START( basic::HIERARCHICAL_ADD_ELEM_TO_CACHE );
	runtime_assert( find(input_addr) == pool_cache_.end() );
	PoolData newpooldat = PoolData( Pool_RMSD_OP( new Pool_RMSD() ), input_addr );
	newpooldat.pool_->add( coords, coords.u2(), tag );
	while ( pool_cache_.size() >= max_cache_size_ ) {
		pool_cache_.pop_back();
	}
	pool_cache_.push_front( newpooldat );
	PROF_STOP( basic::HIERARCHICAL_ADD_ELEM_TO_CACHE );
}

bool HierarchicalLevel::add_new(
	FArray2D_double& coord,
	std::string & tag,
	Address & address){
	return add_new( coord, tag, address, false );
}


bool HierarchicalLevel::add_new(
	FArray2D_double& coord,
	std::string & tag,
	Address & address,
	bool resolve_remaining_levels //default is false
) {
	using namespace basic;
	PROF_START( basic::HIERARCHICAL_ADD );
	//what if we're adding to a pool that exists but we haven't loaded yet?
	//check for existence, load if exists, and then add
	if ( !address_exists_in_cache( address ) && pool_exists( address ) ) {
		TR.Error << " CANNOT ADD IF ADDRESS DOESN'T EXIST IN CACHE. exiting without adding!" << std::endl;
		utility_exit();
		return false;
	}

	bool added_level = false;
	core::Size address_at_level = address[ level_ ];
	if ( ( address_at_level > 0  && find( address ) == pool_cache_.end() ) ||
			address_at_level == 0 ) {
		Address new_level_addr = address;
		for ( core::Size ii = ( level_ + 1 ); ii <= new_level_addr.size(); ii++ ) {
			new_level_addr[ ii ] = 0;
		}
		add_elem_to_cache( coord, tag, new_level_addr );
		added_level = true;
	}
	if ( has_next_level() ) {
		if ( address[ level_ + 1 ] == 0  ) {
			std::list<PoolData>::iterator itr;
			if ( level_find( address, level_, itr ) ) {
				if ( !added_level ) {
					address[ level_ + 1 ] = (*itr).pool_->size() + 1;
				} else {
					runtime_assert( (*itr).pool_->size() > 0 );
					address[ level_ + 1 ] = (*itr).pool_->size(); //to avoid double-counting
				}
			} else {
				address[ level_ + 1 ] = 1;
			}
		}
		bool added_levels_to_child = next_level_->add_new( coord, tag, address, resolve_remaining_levels );
		if ( added_levels_to_child  && !added_level ) { //second condition avoids duplicate addition
			( *find( address ) ).pool_->add( coord, coord.u2(), tag );
		}
	} else {
		if ( address_at_level > 0 && find( address ) != pool_cache_.end() && !added_level ) { //no new levels, adding to lowest pool
			( *find( address ) ).pool_->add( coord, coord.u2(), tag );
		}
	}
	PROF_STOP( basic::HIERARCHICAL_ADD );
	return added_level;
}


void
HierarchicalLevel::add_to_top_level(FArray2D_double& coords, std::string & tag) {
	//if the structure is new in the top-level and not in our branch,
	// we need to be able to evaluate this structure, but we don't want to clutter our
	// hierarchy by creating a tree-structure... so we avoid it by only adding to the top-level
	// if we ever hit this structure later in the trajectory, we'll load the branch from file
	runtime_assert(level_ == 1 && pool_cache_.size() == 1);
	(*pool_cache_.begin()).pool_->add( coords, coords.u2(), tag);
}

void
HierarchicalLevel::print_addresses_at_level(){
	for ( std::list< PoolData >::iterator it = pool_cache_.begin(), end = pool_cache_.end(); it != end; ++it ) {
		print_address( (*it).address_ );
	}
}

void
HierarchicalLevel::debug_print_addresses() {
	if ( TR.visible() ) TR.Debug << "level: " << level_ << " ";
	for ( std::list<PoolData>::iterator itr = pool_cache_.begin(), end = pool_cache_.end(); itr != end; ++itr ) {
		Address addr = (*itr).address_;
		if ( TR.visible() ) {
			for ( core::Size ii = 1; ii <= addr.size(); ii++ ) {
				TR.Debug << addr[ ii ] << " ";
			}
			TR.Debug << " , ";
		}
	}
	if ( TR.visible() ) TR.Debug << std::endl;
	if ( has_next_level() ) {
		next_level_->debug_print_addresses();
	}
}

void
HierarchicalLevel::debug_print_hierarchy(){
	// runtime_assert( level_ == 1 );
	if ( TR.visible() ) { TR.Debug << "level: " << level_ << " size of level: " << size() << std::endl; }
	for ( std::list< PoolData >::iterator itr = pool_cache_.begin(), end = pool_cache_.end(); itr != end; ++itr ) {
		Pool_RMSD_OP poolop = (*itr).pool_;
		if ( TR.visible() ) {
			TR.Debug << "pool at address: "; print_address( (*itr).address_ );
			TR.Debug << " expecting: " << poolop->size() << std::endl;
			for ( core::Size ii = 1; ii <= poolop->size(); ii++ ) {
				TR.Debug << poolop->get_tag( ii ) << std::endl;
			}
			TR.Debug << "end pool at this level" << std::endl;
		}
	}
	if ( has_next_level() ) {
		next_level_->debug_print_hierarchy();
	}
}

void
HierarchicalLevel::debug_print_size_per_level() {
	if ( TR.visible() ) TR.Debug << "size per level: " << level_ << " " << this->size() << " ";
	for ( std::list< PoolData >::iterator itr = pool_cache_.begin(), end = pool_cache_.end(); itr != end; ++itr ) {
		if ( TR.visible() ) {
			TR.Debug << "address: ";
			for ( core::Size ii = 1; ii <= (*itr).address_.size(); ii++ ) {
				TR.Debug << (*itr).address_[ ii ] << " ";
			}
			TR.Debug << ", size: " << (*itr).pool_->size() << " | ";
			TR.Debug << std::endl;
		}
	}
	if ( has_next_level() ) {
		next_level_->debug_print_size_per_level();
	}

}

bool
HierarchicalLevel::address_exists_in_cache( Address & address ) {
	bool has_address = (find(address) != pool_cache_.end()) || ( address[ level_ ] == 0 );
	if ( this->has_next_level() && next_level_->level() <= address.size() && address[ level_ + 1 ] != 0 ) {
		//TR.Debug << "going one level down from " << level_ << " to " << next_level_->level() << std::endl;
		has_address = has_address && next_level_->address_exists_in_cache( address );
	}
	//TR.Debug << "address exists in cache: " << has_address << ": "; print_address(address);
	return has_address;

}

void
HierarchicalLevel::fill_top_level( Pool_RMSD_OP & top_level_pool ) {
	if ( level_ == 1 ) {
		sort_pool( top_level_pool );
		(*pool_cache_.begin()).pool_ = top_level_pool;
	}
}


PoolData
HierarchicalLevel::load_pool( utility::vector1< core::Size > & address ) {
	PROF_START( basic::LOAD_HIERARCHY );
	std::string file = lib_full_path(address);
	/**
	utility::io::izstream testin( file.c_str() );
	TR.Debug << "stream is " << (testin.good() ? "good" : "bad" ) << std::endl;
	while( !testin.good() ) {
	TR.Debug << "stream is " << (testin.good() ? "good" : "bad" ) << std::endl;
	testin.clear();
	testin.close();
	sleep( 5 );
	testin.open( file.c_str() );
	}
	**/
	Pool_RMSD_OP pool_ptr( new Pool_RMSD(file) );
	core::Size ntries = 5;
	while ( pool_ptr->size() == 0 && ntries > 0 ) { //we never try to open empty files
#ifdef _WIN32
#ifndef WIN32
		Sleep(5000);  // REQUIRED FOR WINDOWS
#endif
#else
		sleep(5);
#endif
		pool_ptr = Pool_RMSD_OP( new Pool_RMSD( file ) );
		ntries--;
	}
	runtime_assert( pool_ptr->size() > 0 );
	/**
	if( pool_ptr->size() == 0 ) {
	//runtime_assert( utility::file::file_exists( file ) );
	while( pool_ptr->size() == 0 ) {
	sleep( 5 );
	pool_ptr = new Pool_RMSD( file );
	}
	}
	**/
	sort_pool( pool_ptr );
	PROF_STOP( basic::LOAD_HIERARCHY );
	return PoolData( pool_ptr, address );
}

bool
HierarchicalLevel::pool_exists( utility::vector1< core::Size > & address) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::string filename = lib_full_path( address );
	return utility::file::file_exists( filename ); //ek old way... replaced 10/8/2010
	/**
	utility::io::izstream testin( filename.c_str() );
	TR.Debug << "testing stream. stream is " << ( testin.good() ? "good" : "bad") << std::endl;
	core::Size ntries = 3;
	while( !testin.good() && ntries > 0 ) {
	sleep( 5 );
	ntries--;
	}
	if( !testin.good() ) {
	std::cerr << "expected file: " << filename << " was not found in file-system.. fatal error!" << std::endl;
	TR.Debug << "expected file: " << filename << " was not found in file-system... fatal error!" << std::endl;
	utility_exit();
	}
	**/

}

void
HierarchicalLevel::clear(){
	pool_cache_.clear();
}

utility::vector1< core::Size > &
HierarchicalLevel::address( core::Size index ) {
	core::Size count = 0;
	std::list< PoolData >::iterator itr = pool_cache_.begin();
	while ( count < index ) {
		++itr;
	}
	return (*itr).address_;
}


//for debugging
void
HierarchicalLevel::print_address( utility::vector1< core::Size > & address){
	for ( core::Size ii = 1; ii <= address.size(); ii++ ) {
		TR.Debug << address[ ii ] << " ";
	}
	TR.Debug << std::endl;
}

//get desired lib file name from tag
std::string HierarchicalLevel::lib_full_path( utility::vector1< core::Size > & address ){

	PROF_START( basic::LIB_FULL_PATH );
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace std;
	//TR.Debug << "getting lib_full_path of "; print_address( address );
	runtime_assert(address[ 1 ] > 0 ); //top-level must have something in it

	std::string subdir_prefix("sub_");
	std::string libdir_prefix("c_");
	std::string libdir_suffix("_lib");
	std::string rootpath_suffix("_dir");
	std::string s("/");
	int nofsubdir = option[cluster::K_n_sub];
	//std::string rootpath = utility::file::FileName(option[mc::hierarchical_pool]()).path();
	std::string defaultpath = utility::file::FileName(option[mc::hierarchical_pool]()).base() + rootpath_suffix + s;

	//Size n = address.size();


	std::ostringstream fnstream;
	//fnstream << "c_" << setfill ('0') << setw (5) << address[n-1] << ".out";

	core::Size last_num_not_zero;
	if ( address[ 1 ] == 0 ) {
		last_num_not_zero = 1;
	} else {
		for ( last_num_not_zero = 1; last_num_not_zero < address.size(); last_num_not_zero++ ) { //last number is not a level address
			if ( address[ last_num_not_zero + 1 ] == 0 ) {
				break;
			}
		}
	}
	runtime_assert( last_num_not_zero > 0 ); //expect no elements to be sdasdzero

	fnstream << "c_" << setfill ('0') << setw (5) << address[last_num_not_zero] << ".out";
	std::string fn(fnstream.str());
	//TR.Debug << "silent file we're lookin for: " << fnstream.str() << std::endl;
	for ( Size i=last_num_not_zero; i >= 1; i-- ) {
		ostringstream pathstream1;
		ostringstream pathstream2;

		//main dir
		if ( i>1 ) pathstream1 << libdir_prefix << setfill ('0') << setw (5) << address[i-1] << libdir_suffix << s ;
		//sub dir
		pathstream2 << subdir_prefix << setfill ('0') << setw (3) << int((address[i]-1)/nofsubdir) << s;
		fn = pathstream1.str()+pathstream2.str()+fn;
		//TR.Debug << "fn is now " << fn << std::endl;
	}
	if ( TR.visible() ) {
		TR.Debug << "LIB_FULL_PATH: from address: ";
		for ( core::Size ii = 1; ii <= address.size(); ii++ ) {
			TR.Debug << address[ ii ] << " ";
		}
		TR.Debug << " got file-name: " << defaultpath+fn;
	}
	PROF_STOP( basic::LIB_FULL_PATH );
	return defaultpath+fn;

}

/**
std::string HierarchicalLevel::lib_full_path(std::string tag_origin)
{
using namespace std;
using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

std::string subdir_prefix("sub_");
std::string libdir_prefix("c_");
std::string libdir_suffix("_lib");
std::string s("/");
//int nofsubdir = option[cluster::K_n_sub];
std::string rootpath = utility::file::FileName(option[mc::known_structures]()).path();
std::string defaultpath = utility::file::FileName(option[mc::known_structures]()).base() + libdir_suffix + s;

std::string tag = tag_origin + ".0000";

Size pos=1;
Size len=tag.length();
utility::vector1<Size> address;
while (pos<len)
{
Size newpos = tag.find(".", pos+1);
address.push_back(atoi(tag.substr(pos+1,newpos-pos-1).c_str()));
pos = newpos;
}
return lib_full_path( address );
}
**/


}
}
}
