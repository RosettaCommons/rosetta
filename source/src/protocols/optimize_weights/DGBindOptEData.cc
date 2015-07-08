// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/optimize_weights/DGBindOptEData.cc
///
/// @brief
/// @author Ian W. Davis

#ifdef USEMPI
#include <mpi.h>
#endif

#include <protocols/optimize_weights/DGBindOptEData.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace optimize_weights {


DGBindOptEData::DGBindOptEData():
	deltaG_bind_(0),
	bound_(/* NULL */),
	unbound_(/* NULL */)
{
}


DGBindOptEData::~DGBindOptEData() {}


core::Real
DGBindOptEData::do_score(
	std::ostream & ostr,
	Multivec const & component_weights,
	Multivec const & vars,
	Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const ,//num_ref_dofs,
	int const ,//num_total_dofs,
	EnergyMap const & fixed_terms,
	ScoreTypes const & ,//score_list,
	ScoreTypes const & fixed_score_list,
	bool const print
) const
{
	Real boundE = 0, unboundE = 0;
	for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
		boundE   += vars[ ii ] * bound_->free_data()[ ii ];
		unboundE += vars[ ii ] * unbound_->free_data()[ ii ];
	}
	for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
		boundE   += fixed_terms[ fixed_score_list[ ii ] ] * bound_->fixed_data()[ ii ];
		unboundE += fixed_terms[ fixed_score_list[ ii ] ] * unbound_->fixed_data()[ ii ];
	}

	// PdbBind 2007 core set
	// Linear regression of components of interface_delta gave an intercept between -3 and -4;
	// including this allows a better overall fit to the real binding energy.
	Real const TdS = 3.5; // kcal/mol
	//Real const TdS = 0; // kcal/mol
	Real const pred_dG = (boundE - unboundE) - TdS;
	Real const diff_dG = pred_dG - deltaG_bind_;
	Real const sq_err = diff_dG * diff_dG;

	for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
		dE_dvars[ ii ] += component_weights[ type() ] * 2 * diff_dG * ( bound_->free_data()[ ii ] - unbound_->free_data()[ ii ] );
	}

	if( print ) {
		ostr << "DGBindOptEData bound[ " << boundE << " ] - unbound[ " << unboundE << " ] = predicted[ " << pred_dG << " ], "
			<< "predicted[ " << pred_dG << " ] - experimental[ " << deltaG_bind_ << " ] = error[ " << diff_dG << " ], "
			<< "error^2[ " << sq_err << " ], weighted error^2[ " << component_weights[ type() ]*sq_err << " ]" << std::endl;

	}

	return component_weights[ type() ] * sq_err;
}


/// @brief Return the upper and lower bound on the unweighted components at this
/// position if they are larger (or smaller) than the unweighted values already in
/// the two input EnergyMaps.
void
DGBindOptEData::range(
	ScoreTypes const & free_score_list,
	ScoreTypes const & fixed_score_list,
	EnergyMap & lower_bound,
	EnergyMap & upper_bound
) const
{
	update_range( bound_,   free_score_list, fixed_score_list, lower_bound, upper_bound );
	update_range( unbound_, free_score_list, fixed_score_list, lower_bound, upper_bound );
}


core::Size
DGBindOptEData::memory_use() const
{
	Size total = sizeof( DGBindOptEData ) +
		sizeof( SingleStructureData ) * 1 +
		sizeof( SingleStructureData ) * 1;
	total += sizeof(Real) * (bound_->free_data().size()+bound_->fixed_data().size()) * 1;
	total += sizeof(Real) * (unbound_->free_data().size()+unbound_->fixed_data().size()) * 1;
	return total;
}


#ifdef USEMPI
void
DGBindOptEData::send_to_node( int const destination_node, int const tag ) const
{
	/// 1. Experimental dG
	Real deltaG_bind = deltaG_bind_; // stupid const pointer
	MPI_Send( & deltaG_bind, 1, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	/// 2. n free
	Size n_free = bound_->free_data().size();
	//std::cout << "sending n_free to node " << destination_node << " " << n_free << std::endl;
	MPI_Send( & n_free, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	/// 3. n fixed
	Size n_fixed = bound_->fixed_data().size();
	//std::cout << "sending n_fixed to node " << destination_node  << " " << n_fixed << std::endl;
	MPI_Send( & n_fixed, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	/// Send natives, then send decoys
	Real * free_data = new Real[ n_free ];
	Real * fixed_data = new Real[ n_fixed ];
	for ( Size jj = 1; jj <= n_free; ++jj ) {
		free_data[ ( jj - 1 ) ] = bound_->free_data()[ jj ];
	}
	for ( Size jj = 1; jj <= n_fixed; ++jj ) {
		fixed_data[ ( jj - 1 ) ] = bound_->fixed_data()[ jj ];
	}

	//std::cout << "sending native free_data to node " << destination_node << " " << free_data <<  std::endl;
	/// 5. native free data
	MPI_Send( free_data, n_free, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	//std::cout << "sending native fixed_data to node " << destination_node << " " << fixed_data <<  std::endl;
	/// 6. fixed data
	MPI_Send( fixed_data, n_fixed, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	//std::cout << "Sent -- about to delete data" << std::endl;

	/// now send decoys
	Real * decoy_free_data = new Real[ n_free ];
	Real * decoy_fixed_data = new Real[ n_fixed ];
	for ( Size jj = 1; jj <= n_free; ++jj ) {
		decoy_free_data[ ( jj - 1 ) ] = unbound_->free_data()[ jj ];
	}
	for ( Size jj = 1; jj <= n_fixed; ++jj ) {
		decoy_fixed_data[ ( jj - 1 ) ] = unbound_->fixed_data()[ jj ];
	}
	/// 7. decoy free data
	//std::cout << "sending decoy free_data to node " << destination_node << std::endl;
	MPI_Send( decoy_free_data, n_free, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	/// 8. decoy fixed data
	//std::cout << "sending decoy fixed_data to node " << destination_node << std::endl;
	MPI_Send( decoy_fixed_data, n_fixed, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	delete [] free_data;
	delete [] fixed_data;

	delete [] decoy_free_data;
	delete [] decoy_fixed_data;

	OptEPositionData::send_to_node( destination_node, tag );

}


void
DGBindOptEData::receive_from_node( int const source_node, int const tag )
{
	MPI_Status stat;
	//TR << "PNatStructureOptEData::Recieving data from node... " << source_node << std::endl;

	/// 1. Experimental dG
	MPI_Recv( & deltaG_bind_, 1, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 2. n free
	Size n_free( 0 );
	MPI_Recv( & n_free, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 3. n fixed
	Size n_fixed( 0 );
	MPI_Recv( & n_fixed, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );

	/// Recieve native data first, then decoys
	Real * free_data = new Real[ n_free ];
	Real * fixed_data = new Real[ n_fixed ];

	/// 5. free data
	MPI_Recv( free_data, n_free, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 6. fixed data
	MPI_Recv( fixed_data, n_fixed, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	utility::vector1< Real > free_data_v( n_free );
	utility::vector1< Real > fixed_data_v( n_fixed );
	for ( Size jj = 1; jj <= n_free; ++jj ) {
		free_data_v[ jj ] = free_data[ ( jj - 1 ) ];
	}
	for ( Size jj = 1; jj <= n_fixed; ++jj ) {
		fixed_data_v[ jj ] = fixed_data[ ( jj - 1 ) ];
	}
	bound_ = SingleStructureDataOP( new SingleStructureData( free_data_v, fixed_data_v ) );


	delete [] free_data;// free_data = 0;
	delete [] fixed_data;// fixed_data = 0;

	//// Now receive decoy data
	free_data = new Real[ n_free ];
	fixed_data = new Real[ n_fixed ];

	/// 5. free data
	MPI_Recv( free_data, n_free, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 6. fixed data
	MPI_Recv( fixed_data, n_fixed, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	for ( Size jj = 1; jj <= n_free; ++jj ) {
		free_data_v[ jj ] = free_data[ ( jj - 1 ) ];
	}
	for ( Size jj = 1; jj <= n_fixed; ++jj ) {
		fixed_data_v[ jj ] = fixed_data[ ( jj - 1 ) ];
	}
	unbound_ = SingleStructureDataOP( new SingleStructureData( free_data_v, fixed_data_v ) );

	delete [] free_data;
	delete [] fixed_data;

	OptEPositionData::receive_from_node( source_node, tag );

}
#endif // USEMPI


} // namespace optimize_weights
} // namespace protocols
