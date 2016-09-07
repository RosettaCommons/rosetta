// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author ashworth

#include <protocols/dna/DnaInterfaceFinder.hh>

#include <protocols/dna/util.hh>

#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>

#include <utility/exit.hh> // runtime_assert

#include <utility/vector1.hh>


namespace protocols {
namespace dna {

using namespace core;
using namespace conformation;
using utility::vector1;

void
DnaInterfaceFinder::determine_interface( pose::Pose const & pose )
{
	// Does both protein and DNA
	vector1< Size > protein_positions, dna_positions;
	for ( Size index(1), end( pose.total_residue() ); index <= end; ++index ) {
		if ( pose.residue_type( index ).is_protein() ) protein_positions.push_back( index );
		else if ( pose.residue_type( index ).is_DNA() ) dna_positions.push_back( index );
	}
	determine_protein_interface( pose, protein_positions, dna_positions );
	determine_dna_interface( pose, protein_positions, dna_positions );
}

void
DnaInterfaceFinder::determine_protein_interface( pose::Pose const & pose )
{
	// all protein positions against all DNA positions
	// this is the slowest possible (intended) usage of DnaInterfaceFinder
	vector1< Size > protein_positions, dna_positions;
	for ( Size index(1), end( pose.total_residue() ); index <= end; ++index ) {
		if ( pose.residue_type( index ).is_protein() ) protein_positions.push_back( index );
		else if ( pose.residue_type( index ).is_DNA() ) dna_positions.push_back( index );
	}
	determine_protein_interface( pose, protein_positions, dna_positions );
}

void
DnaInterfaceFinder::determine_interface(
	pose::Pose const & pose,
	vector1< Size > const & protein_positions,
	vector1< Size > const & dna_positions
)
{
	determine_protein_interface( pose, protein_positions, dna_positions );
	determine_dna_interface( pose, protein_positions, dna_positions );
}

void
DnaInterfaceFinder::determine_protein_interface(
	pose::Pose const & pose,
	vector1< Size > const & protein_positions,
	vector1< Size > const & dna_positions
)
{
	for ( unsigned long protein_position : protein_positions ) {
		runtime_assert( pose.residue_type( protein_position ).is_protein() );

		Residue const & pres( pose.residue( protein_position ) );
		protein_neighbors_[ protein_position ] = DnaNeighbor();
		DnaNeighbor & neighbor( protein_neighbors_[ protein_position ] );

		Real shortest_arg_dis2(10000);
		for ( auto dna_index( dna_positions.begin() ),
				end( dna_positions.end() ); dna_index != end && !neighbor.contact(); ++dna_index ) {
			runtime_assert( pose.residue_type( *dna_index ).is_DNA() );

			Residue const & dres( pose.residue( *dna_index ) );
			neighbor.close( close_to_dna( pres, dres, close_threshold_, base_only_ ) );
			if ( neighbor.close() ) {
				if ( z_axis_dist( pres, dres ) < z_cutoff_ ) {
					Real dis2(
						argrot_dna_dis2( pose, protein_position, pres, dres, contact_threshold_, base_only_ )
					);
					if ( dis2 < shortest_arg_dis2 ) shortest_arg_dis2 = dis2;
					if ( shortest_arg_dis2 < contact_threshold_ ) neighbor.contact(true);
				}
			}
		}
	}
	initialized_ = true;
}

void
DnaInterfaceFinder::determine_dna_interface(
	pose::Pose const & pose,
	vector1< Size > const & protein_positions,
	vector1< Size > const & dna_positions
)
{
	for ( unsigned long dna_position : dna_positions ) {
		runtime_assert( pose.residue_type( dna_position ).is_DNA() );

		Residue const & dres( pose.residue( dna_position ) );
		dna_neighbors_[ dna_position ] = DnaNeighbor();
		DnaNeighbor & neighbor( dna_neighbors_[ dna_position ] );

		Real shortest_arg_dis2(10000);
		for ( auto p_index( protein_positions.begin() ),
				end( protein_positions.end() ); p_index != end && !neighbor.contact(); ++p_index ) {
			runtime_assert( pose.residue_type( *p_index ).is_protein() );

			Residue const & pres( pose.residue( *p_index ) );
			neighbor.close( close_to_dna( pres, dres, close_threshold_, base_only_ ) );
			if ( neighbor.close() ) {
				if ( z_axis_dist( pres, dres ) < z_cutoff_ ) {
					Real dis2(
						argrot_dna_dis2( pose, *p_index, pres, dres, contact_threshold_, base_only_ )
					);
					if ( dis2 < shortest_arg_dis2 ) shortest_arg_dis2 = dis2;
					if ( shortest_arg_dis2 < contact_threshold_ ) neighbor.contact(true);
				}
			}
		}
	}
	initialized_ = true;
}

DnaNeighbors const &
DnaInterfaceFinder::protein_neighbors() const
{
	runtime_assert( initialized() );
	return protein_neighbors_;
}

DnaNeighbors const &
DnaInterfaceFinder::dna_neighbors() const
{
	runtime_assert( initialized() );
	return dna_neighbors_;
}

} // namespace dna
} // namespace protocols
