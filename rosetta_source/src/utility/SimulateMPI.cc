// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/SimulateMPI.cc
/// @brief
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#include <utility/SimulateMPI.hh>
#include <ObjexxFCL/StaticIndexRange.hh>

namespace utility {

// initialize private static data
SimulateMPIData * SimulateMPI::simulation_;

SimulateMPIData::SimulateMPIData() :
	mpi_rank_(0),
	mpi_nprocs_(0),
	proc_string_buf_(),
	proc_char_buf_(),
	proc_integer_buf_(),
	proc_integers_buf_(),
	proc_double_buf_(),
	proc_doubles_buf_()
{}

SimulateMPIData::SimulateMPIData( const SimulateMPIData & src ) :
	mpi_rank_(src.mpi_rank_),
	mpi_nprocs_(src.mpi_nprocs_),
	proc_string_buf_(src.proc_string_buf_),
	proc_char_buf_(src.proc_char_buf_),
	proc_integer_buf_(src.proc_integer_buf_),
	proc_integers_buf_(src.proc_integers_buf_),
	proc_double_buf_(src.proc_double_buf_),
	proc_doubles_buf_(src.proc_doubles_buf_)
{}

SimulateMPI::SimulateMPI() {}

SimulateMPI::SimulateMPI( const SimulateMPI & ) {}

void
SimulateMPI::initialize_simulation( int nprocs ) {
	simulation_ = new SimulateMPIData();
	SimulateMPI::set_mpi_nprocs(nprocs);
}

bool
SimulateMPI::simulate_mpi() {
	return simulation_;
}


void
SimulateMPI::set_mpi_rank(int value) {
	assert(simulation_);
	assert(0 <= value);
	assert(value < mpi_nprocs());
	simulation_->mpi_rank_ = value;
}


int
SimulateMPI::mpi_rank() {
	assert(simulation_);
	return simulation_->mpi_rank_;
}

void
SimulateMPI::set_mpi_nprocs(int value) {
	assert(simulation_);
	simulation_->mpi_nprocs_ = value;

	ObjexxFCL::StaticIndexRange dim(0,value-1);

	simulation_->proc_string_buf_.dimension(dim,dim);
	simulation_->proc_char_buf_.dimension(dim,dim);
	simulation_->proc_integer_buf_.dimension(dim,dim);
	simulation_->proc_integers_buf_.dimension(dim,dim);
	simulation_->proc_double_buf_.dimension(dim,dim);
	simulation_->proc_doubles_buf_.dimension(dim,dim);
}

int
SimulateMPI::mpi_nprocs() {
	return simulation_->mpi_nprocs_;
}

std::string
SimulateMPI::receive_string_from_node(
	int source
) {
	assert(simulation_);
	assert(0 <= source);
	assert(source < mpi_nprocs());
	std::string message( simulation_->proc_string_buf_(source,mpi_rank()) );
	simulation_->proc_string_buf_(source,mpi_rank()) = "";
	return message;
}

void
SimulateMPI::send_string_to_node(
	int destination,
	std::string const & message
) {
	assert(simulation_);
	assert(0 <= destination);
	assert(destination < mpi_nprocs());
	simulation_->proc_string_buf_(mpi_rank(),destination) = message;
}


char
SimulateMPI::receive_char_from_node(
	int source
) {
	assert(simulation_);
	assert(0 <= source);
	assert(source < mpi_nprocs());

	char message( simulation_->proc_char_buf_(source,mpi_rank()) );
	simulation_->proc_char_buf_(source,mpi_rank()) = 0;
	return message;
}

void
SimulateMPI::send_char_to_node(
	int destination,
	char message
) {
	assert(simulation_);
	assert(0 <= destination);
	assert(destination < mpi_nprocs());

	simulation_->proc_char_buf_(mpi_rank(),destination) = message;
}


int
SimulateMPI::receive_integer_from_node(
	int source
) {
	assert(simulation_);
	assert(0 <= source);
	assert(source < mpi_nprocs());

	int message( simulation_->proc_integer_buf_(source,mpi_rank()) );
	simulation_->proc_integer_buf_(source,mpi_rank()) = 0;
	return message;
}

void
SimulateMPI::send_integer_to_node(
	int destination,
	int message
) {
	assert(simulation_);
	assert(0 <= destination);
	assert(destination < mpi_nprocs());

	simulation_->proc_integer_buf_(mpi_rank(),destination) = message;
}

vector1< int >
SimulateMPI::receive_integers_from_node(
	int source
) {
	assert(simulation_);
	assert(0 <= source);
	assert(source < mpi_nprocs());

	vector1< int > message( simulation_->proc_integers_buf_(source,mpi_rank()) );
	simulation_->proc_integers_buf_(source,mpi_rank()) = vector1< int >();
	return message;
}

void
SimulateMPI::send_integers_to_node(
	int destination,
	vector1< int > const & message
) {
	assert(simulation_);
	assert(0 <= destination);
	assert(destination < mpi_nprocs());

	simulation_->proc_integers_buf_(mpi_rank(),destination) = message;
}


double
SimulateMPI::receive_double_from_node(
	int source
) {
	assert(simulation_);
	assert(0 <= source);
	assert(source < mpi_nprocs());

	double message( simulation_->proc_double_buf_(source,mpi_rank()) );
	simulation_->proc_double_buf_(source,mpi_rank()) = 0;
	return message;
}

void
SimulateMPI::send_double_to_node(
	int destination,
	double message
) {
	assert(simulation_);
	assert(0 <= destination);
	assert(destination < mpi_nprocs());

	simulation_->proc_double_buf_(mpi_rank(),destination) = message;
}

vector1< double >
SimulateMPI::receive_doubles_from_node(
	int source
) {
	assert(simulation_);
	assert(0 <= source);
	assert(source < mpi_nprocs());

	vector1< double > message( simulation_->proc_doubles_buf_(source,mpi_rank()) );
	simulation_->proc_doubles_buf_(source,mpi_rank()) = vector1<double>();
	return message;
}

void
SimulateMPI::send_doubles_to_node(
	int destination,
	vector1< double > const & message
) {
	assert(simulation_);
	assert(0 <= destination);
	assert(destination < mpi_nprocs());

	simulation_->proc_doubles_buf_(mpi_rank(),destination) = message;
}


} // utility
