// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/canonical_sampling/BiasEnergyMover.cc
/// @brief BiasEnergy methods implemented
/// @author

#ifndef INCLUDED_protocols_canonical_sampling_AsyncMPITemperingBase_hh
#define INCLUDED_protocols_canonical_sampling_AsyncMPITemperingBase_hh


#include <protocols/canonical_sampling/BiasEnergy.hh>

#include <ObjexxFCL/string.functions.hh>
#include <utility/exit.hh>
#include <math.h>
#include <iomanip>


namespace protocols {
namespace canonical_sampling {

template <typename FLOAT>
void _toggle( FLOAT *& ptr1, FLOAT *& ptr2 ) {
	FLOAT *ptrX = ptr1;
	ptr1 = ptr2;
	ptr2 = ptrX;
}

template< typename T>
void BiasEnergy::Histogram<T>::toggle() {
	_toggle( data_, recv_buf_ );
}

template< typename T>
BiasEnergy::Histogram<T>::Histogram( core::Real grid_min, core::Real grid_max, core::Size ngrid_cells /* = 100 */ ) :
data_( NULL ),
recv_buf_( NULL ),
grid_min_( grid_min ),
grid_max_( grid_max ),
delta_grid_( (grid_max-grid_min)/ngrid_cells )
{
	ngrid_cells_ = ngrid_cells;
	runtime_assert( ngrid_cells_ > 0 );
	data_ = new ValueType[ ngrid_cells_ ];
	recv_buf_ = new ValueType[ ngrid_cells_ ];
	reset();
}

template< typename T>
BiasEnergy::Histogram<T>::Histogram( BiasEnergy::Histogram<T> const& other ) :
Parent( other ),
grid_min_( other.grid_min_ ),
grid_max_( other.grid_max_ ),
delta_grid_( other.delta_grid_ ),
ngrid_cells_( other.ngrid_cells_ )
{
	data_ = new ValueType[ other.ngrid_cells_ ];
	recv_buf_ = new ValueType[ other.ngrid_cells_ ];
	copy_data( other );
}

template< typename T>
void BiasEnergy::Histogram<T>::reset() {
	for ( core::Size i = 0; i < ngrid_cells_; i++ ) {
		data_[ i ] = 0.0;
		recv_buf_[ i ] = 0.0;
	}
}

template< typename T>
void BiasEnergy::Histogram<T>::copy_data( BiasEnergy::Histogram<T> const& other ) {
	for ( core::Size i = 0; i < ngrid_cells_; i++ ) {
		data_[ i ] = other.data_[ i ];
		recv_buf_[ i ] = other.recv_buf_[ i ];
	}
}

template< typename T>
BiasEnergy::Histogram<T>& BiasEnergy::Histogram<T>::operator=( Histogram<T> const& other ) {
	grid_min_ = other.grid_min_;
	grid_max_ = other.grid_max_;
	delta_grid_ = other.delta_grid_;
	ngrid_cells_ = other.ngrid_cells_;
	if ( other.ngrid_cells_ != ngrid_cells_ ) {
		delete [] data_;
		delete [] recv_buf_;
		data_ = new ValueType[ other.ngrid_cells_ ];
		recv_buf_ = new ValueType[ other.ngrid_cells_ ];
	}
	copy_data( other );
	return *this;
}

template< typename T>
BiasEnergy::Histogram<T>::~Histogram() {
	delete [] data_;
	delete [] recv_buf_;
}

#ifdef USEMPI

template< typename T>
MPI_Datatype mpi_data_type() {
	utility_exit_with_message("not implemented");
	return MPI_FLOAT;
}

template<>
MPI_Datatype mpi_data_type<float>() {
	return MPI_FLOAT;
}

template<>
MPI_Datatype mpi_data_type<int>() {
	return MPI_INT;
}


template< typename T>
void BiasEnergy::Histogram<T>::mpi_exchange( int exchange_partner, MPI_Comm const& mpi_comm ) {
	MPI_Request reqs[2];
	MPI_Status stats[2];
	int const mpi_GRID_XCHANGE = 30;
	MPI_Isend( data_, ngrid_cells_, mpi_data_type<T>(), exchange_partner, mpi_GRID_XCHANGE, mpi_comm, &reqs[0]);
	MPI_Irecv( recv_buf_, ngrid_cells_, mpi_data_type<T>(), exchange_partner, mpi_GRID_XCHANGE, mpi_comm, &reqs[1]);
	MPI_Waitall(2, reqs, stats);
}
#endif

template< typename T>
bool BiasEnergy::Histogram<T>::check_range( core::Real val ) const {
	return val > grid_min_ && val <= grid_max_ ;
}

template< typename T>
core::Size BiasEnergy::Histogram<T>::cell_index( core::Real val ) const {
	//max(1,min(ebin_num,ceil((Vtrial_unbiased-e_min)/grid_dE)));
	int index( ceil(( val - grid_min_ ) / delta_grid_ ) );
	runtime_assert( ngrid_cells_ > 1 );
	runtime_assert( index >= 1 );
	runtime_assert( index >= 1 && index <= (int) ngrid_cells_ );
	return index;
}

template< typename T>
void BiasEnergy::Histogram<T>::write_to_stream( std::ostream& os ) const {
	os << "HISTOGRAM " << ngrid_cells_ << " " << grid_min_ << " " << grid_max_ << std::endl;
	os << "GRID_START" << std::endl;
	core::Size line_break_ct = 20;
	for ( core::Size index=1; index <= ngrid_cells_; ++index ) {
		if ( index % line_break_ct == 0 ) os << std::endl;
		os << std::setprecision( 4 ) << at( index ) << " ";
	}
	os << "GRID_END" << std::endl;
}

template< typename T >
void BiasEnergy::Histogram<T>::write_to_string( std::string& str ) const {
	using namespace ObjexxFCL;
	str = str + "HISTOGRAM " + string_of( ngrid_cells_ ) + " " + string_of( grid_min_ ) + " " + string_of( grid_max_ ) + " GRID_START ";
	for ( core::Size index=1; index <= ngrid_cells_; ++index ) {
		str = str + string_of( at(index), 4 ) + " ";
	}
	str = str+"GRID_END ";
}

}
}

#endif // INCLUDED_protocols_canonical_sampling_AsyncMPITemperingBase_hh
