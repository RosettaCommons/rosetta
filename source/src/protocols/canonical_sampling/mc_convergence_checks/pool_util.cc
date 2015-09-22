// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#include <basic/prof.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/pool_util.hh>
#include <basic/Tracer.hh>


static THREAD_LOCAL basic::Tracer TR( "HPool_util" );

#ifdef USEMPI
#include <mpi.h>
#endif

namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {

typedef ObjexxFCL::FArray2D<double> FArray2D_double;
typedef utility::vector1< core::Size > Address;
using namespace core;

DataBuffer::DataBuffer():
	neighbor_addresses_(0),
	coords_transfer_buffer_(0),
	coords_receiving_buffer_(0),
	coords_(),
	temp_coords_(),
	num_new_neighbors_(0),
	memory_offset_(0),
	int_buf1_(0),
	winning_ranks_(0),
	candidate_nbr_index_(0),
	candidate_best_rmsd_(0.0),
	candidate_best_rmsds_(0),
	candidate_address_(),
	winning_tag_(""),
	winning_address_(),
	new_level_begins_(0),
	best_candidate_rmsds_(),
	is_a_neighbor_(),
	finished_(0)
{}

void
DataBuffer::setup( int num_slave_nodes, int nresidues, int nlevels ){
	TR.Debug << "setting up pool-util: num-slave-nodes: " << num_slave_nodes <<
		" num-residues: " << nresidues <<
		" nlevels: " << nlevels << std::endl;

	neighbor_addresses_ = new int[ num_slave_nodes * nlevels ];
	coords_transfer_buffer_ = new double[ (nresidues * 3) ];
	coords_receiving_buffer_ = new double[ (nresidues * 3 * num_slave_nodes ) ];
	coords_ = FArray2D_double( 3, nresidues, 0.0 );
	temp_coords_ = FArray2D_double( 3, nresidues, 0.0 );
	num_new_neighbors_ = 0;
	memory_offset_ = new int[ num_slave_nodes ];
	int_buf1_ = new int[ num_slave_nodes * nlevels ];
	winning_ranks_ = new int[ num_slave_nodes * ( nlevels + 1 ) ];
	candidate_best_rmsds_ = new double[ num_slave_nodes ];
	candidate_address_.resize( nlevels, 0 );
	winning_address_.resize( nlevels, 0 );
	candidate_coords_ = FArray2D_double( 3, nresidues, 0.0);
	is_a_neighbor_ = utility::vector1<bool>( num_slave_nodes, false );
	finished_ = new int[ num_slave_nodes ];
}

DataBuffer::~DataBuffer() {
	delete [] neighbor_addresses_;
	delete [] coords_transfer_buffer_;
	delete [] coords_receiving_buffer_;
	delete [] memory_offset_;
	delete [] int_buf1_;
	delete [] winning_ranks_;
	delete [] candidate_best_rmsds_;
	delete [] finished_;
}

void
DataBuffer::address_to_buf( Address & address, int* buf, core::Size start_index  ) {
	for ( core::Size ii = 1; ii <= address.size(); ii++ ) {
		buf[ ii - 1 + start_index ] = address[ ii ];
	}
}

void
DataBuffer::buf_to_address( Address & address, int* buf, core::Size start_index ){
	for ( core::Size ii = 1; ii <= address.size(); ii++ ) {
		address[ ii ] = buf[ ii - 1 + start_index ];
	}
}


void
DataBuffer::farray_to_array( core::Size index,
	core::Size /*num_to_add*/,
	FArray2D_double const&  coords,
	double* coord_buf ){
	PROF_START( basic::FARRAY_MANIPULATION );
	for ( int ii = 1; ii <= coords.u1(); ii++ ) {
		for ( int jj = 1; jj <= coords.u2(); jj++ ) {
			coord_buf[index++] = coords( ii, jj );
		}
	}
	PROF_STOP( basic::FARRAY_MANIPULATION );
}

void
DataBuffer::farray_to_array( core::Size index,
	FArray2D_double const&  coords,
	double* coord_buf ){
	farray_to_array( index, 1, coords, coord_buf );
}

void
DataBuffer::array_to_farray( core::Size index,
	core::Size /*num_to_add*/,
	FArray2D_double & coords,
	double* coord_buf){
	PROF_START( basic::FARRAY_MANIPULATION );
	for ( int ii = 1; ii <= coords.u1(); ii++ ) {
		for ( int jj = 1; jj <= coords.u2(); jj++ ) {
			coords( ii, jj ) =  coord_buf[index++];
		}
	}
	PROF_STOP( basic::FARRAY_MANIPULATION );
}

void
DataBuffer::array_to_farray( core::Size index,
	FArray2D_double & coords,
	double* coord_buf){
	array_to_farray( index, 1, coords, coord_buf );
}


}
}
}


