// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoringManager.hh
/// @brief  Scoring manager class header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

/// Unit headers
#include <core/scoring/SS_Info.hh>

/// Package headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <ObjexxFCL/format.hh>

/// C++ Headers
#include <iostream>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

/// @details.  After a change in size, the residue types and the integer indices are all wrong.
/// Erase the old information.
void
BB_Pos::resize( int const nres )
{
	if ( N_.size() == (Size) nres ) return;

	N_.resize( nres );
	CA_.resize( nres );
	C_.resize( nres );
	O_.resize( nres );
	CB_.resize( nres );

	residue_types_.resize( nres ); std::fill( residue_types_.begin(), residue_types_.end(), static_cast< chemical::ResidueType const * > (0) );
	N_index_.resize( nres );       std::fill( N_index_.begin(), N_index_.end(), 0 );
	CA_index_.resize( nres );      std::fill( CA_index_.begin(), CA_index_.end(), 0 );
	CB_index_.resize( nres );      std::fill( CB_index_.begin(), CB_index_.end(), 0 );
	C_index_.resize( nres );       std::fill( C_index_.begin(), C_index_.end(), 0 );
	O_index_.resize( nres );       std::fill( O_index_.begin(), O_index_.end(), 0 );

}

/// @details: Optimize the common case where the sequence of the pose is not changing from
/// score function evaluation to evaluation (e.g. abinitio!)
void
BB_Pos::take_coordinates_from_pose( pose::Pose const & pose )
{
	if ( ! bbindices_up_to_date( pose ) ) {
		update_indices( pose );
	}

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		if ( rsd.is_protein() ) {
			assert( N_index_[ i ]  );
			assert( CA_index_[ i ] );
			assert( C_index_[ i ]  );
			assert( O_index_[ i ]  );

			N_ [i] = rsd.xyz( N_index_[ i ]  );
			CA_[i] = rsd.xyz( CA_index_[ i ] );
			C_ [i] = rsd.xyz( C_index_[ i ]  );
			O_ [i] = rsd.xyz( O_index_[ i ]  );
			if ( rsd.aa() == chemical::aa_gly ) {
				CB_[i] = 0.0;
			} else {
				CB_[i] = rsd.xyz( CB_index_[ i ] );
			}
		} else {
			N_ [i] = 0.0;
			CA_[i] = 0.0;
			C_ [i] = 0.0;
			O_ [i] = 0.0;
			CB_[i] = 0.0;
		}
	}

}

bool
BB_Pos::bbindices_up_to_date( pose::Pose const & pose ) const
{
	if ( N_.size() != pose.total_residue() ) return false;

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( residue_types_[ ii ] != & ( pose.residue_type( ii ) ) ) return false;
	}
	return true;
}

void
BB_Pos::update_indices( pose::Pose const & pose )
{
	resize( pose.total_residue() );

	static std::string const bbN("N");
	static std::string const bbCA("CA");
	static std::string const bbC("C");
	static std::string const bbO("O");
	static std::string const scCB("CB");

	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		residue_types_[ i ] = & ( rsd.type() );
		if ( rsd.is_protein() ) {

			N_index_[  i ] = rsd.atom_index( bbN );
			CA_index_[ i ] = rsd.atom_index( bbCA );
			C_index_[  i ] = rsd.atom_index( bbC );
			O_index_[  i ] = rsd.atom_index( bbO );

			if ( rsd.aa() == chemical::aa_gly ) {
				CB_index_[ i ] = 0;
			} else {
				CB_index_[ i ] = rsd.atom_index( scCB );
			}
		} else {
			N_index_[  i ] = 0;
			CA_index_[ i ] = 0;
			C_index_[  i ] = 0;
			O_index_[  i ] = 0;
			CB_index_[ i ] = 0;
		}
	}


}



//////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief default constructor
Strands::Strands()
{}

/// @brief total residue constructor
Strands::Strands(
	int const & total_residue
)
{
	resize( total_residue );
}

/// @brief copy constructor
Strands::Strands(
	Strands const & s
) : total_SS_dimer( s.total_SS_dimer ),
		SS_resnum( s.SS_resnum ),
		total_strands( s.total_strands ),
		SS_strand( s.SS_strand ),
		SS_dimer( s.SS_dimer ),
		SS_strand_end( s.SS_strand_end ),
		dimer_neighbor( s.dimer_neighbor )
	//strand_strand_score( s.strand_strand_score )
{}

/// @brief default destructor
Strands::~Strands()
{}

///
void
Strands::resize( Size const nres )
{
	if ( SS_resnum.size1() == nres ) return;
	SS_resnum.dimension( nres );
	SS_strand.dimension( nres );
	SS_dimer.dimension( nres );
	SS_strand_end.dimension( 2, nres );
	dimer_neighbor.dimension( 2, nres );
	clear();
}

void Strands::clear() {
	SS_resnum = 0;
	SS_strand = 0;
	SS_dimer = 0;
	SS_strand_end = 0;
	dimer_neighbor = 0;
}

/// @brief copy assignment
Strands const &
Strands::operator =( Strands const & s )
{
	if ( this != &s ) {
		total_SS_dimer = s.total_SS_dimer;
		SS_resnum = s.SS_resnum;
		total_strands = s.total_strands;
		SS_strand = s.SS_strand;
		SS_dimer = s.SS_dimer;
		SS_strand_end = s.SS_strand_end;
		dimer_neighbor = s.dimer_neighbor;
		//strand_strand_score = s.strand_strand_score;
	}
	return *this;
}


std::ostream &
operator<< ( std::ostream & out, Strands const & s )
{
	using ObjexxFCL::format::I;
	out << "Strands: " << s.total_strands << " strands, " << s.total_SS_dimer << " dimers.\n";
	for ( int i=1; i<= s.total_SS_dimer; ++i ) {
		out << "Strands: " << I(4,i) << I(4,s.SS_strand(i)) << I(4,s.SS_resnum(i)) << I(4,s.SS_dimer(s.SS_resnum(i))) <<
			I(4,s.SS_strand_end(1,i)+1) << "-" << I(4,s.SS_strand_end(2,i)-1) << '\n';
	}
	return out;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////


/// @brief default constructor
Helices::Helices()
{}

/// @brief total residue constructor
Helices::Helices(
	int const & total_residue
)
{
	resize( total_residue );
}

///
void
Helices::resize( int const nres )
{
	HH_resnum.dimension( nres );
	HH_helix_end.dimension( 2, nres );
}

Helices::Helices(
	Helices const & h
) : total_HH_dimer( h.total_HH_dimer ),
		HH_resnum( h.HH_resnum ),
		HH_helix_end( h.HH_helix_end )
{}

Helices::~Helices()
{}

Helices const &
Helices::operator =( Helices const & h )
{
	if ( this != &h ) {
		total_HH_dimer = h.total_HH_dimer;
		HH_resnum = h.HH_resnum;
		HH_helix_end = h.HH_helix_end;
	}
	return *this;
}

std::ostream &
operator<< ( std::ostream & out, Helices const & s )
{
	using ObjexxFCL::format::I;
	out << "Helices: " << s.total_HH_dimer << " dimers.\n";
	for ( int i=1; i<= s.total_HH_dimer; ++i ) {
		out << "Helices: " << I(4,i) << I(4,s.HH_resnum(i)) << I(4,s.HH_helix_end(1,i)+1) << "-" <<
			I(4,s.HH_helix_end(2,i)-1) << '\n';
	}
	return out;
}



} // ns scoring
} // ns core

