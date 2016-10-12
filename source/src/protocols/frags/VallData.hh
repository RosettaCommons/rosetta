// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_src_devel_blab_classic_frags_VallData_HH
#define INCLUDED_src_devel_blab_classic_frags_VallData_HH


// Rosetta Headers

#include <core/types.hh>
#include <protocols/frags/VallData.fwd.hh>
#include <protocols/frags/TorsionFragment.fwd.hh>

// ObjexxFCL Headers
#include <utility/vector1.hh>
#include <utility/exit.hh>

// C++ Headers
#include <cmath>
#include <map>
//#include <cstdlib>
//#include <stdio>
#include <iostream>
// #include <fstream>
// #include <sstream>

namespace protocols {
namespace frags {

inline
char
torsion2big_bin(
	core::Real const phi,
	core::Real const psi,
	core::Real const omega
)
{
	if ( std::abs( omega ) < 90 ) { // this is not quite right: should be omega BEFORE the residue...
		return 'O'; // cis-omega
	} else if ( phi >= 0.0 ) {
		if ( -100 < psi && psi <= 100 ) {
			return 'G'; // alpha-L
		} else {
			return 'E'; // E
		}
	} else {
		if ( -125 < psi && psi <= 50 ) {
			return 'A'; // helical
		} else {
			return 'B'; // extended
		}
	}
	return 'X';
}


////
class VallData {
public:
	/// default constructor
	VallData ()
	{
		core::Size const big_size( 100000 ); // rough guess
		sequence_.reserve( big_size );
		secstruct_.reserve( big_size );
		bigbin_.reserve( big_size );
		phi_.reserve( big_size );
		psi_.reserve( big_size );
		omega_.reserve( big_size );
	}

	/// constructor from input vall database file
	VallData ( std::string const & filename )
	{
		core::Size const big_size( 100000 ); // rough guess
		// prevent lots of redimensioning as we read file? does this even matter?
		sequence_.reserve( big_size );
		secstruct_.reserve( big_size );
		bigbin_.reserve( big_size );
		phi_.reserve( big_size );
		psi_.reserve( big_size );
		omega_.reserve( big_size );
		read_file( filename );
	}

	/// removes excess storage capacity to minimize memory usage
	void
	shrink()
	{
		sequence_.shrink();
		secstruct_.shrink();
		bigbin_.shrink();
		phi_.shrink();
		psi_.shrink();
		omega_.shrink();
	}

	// read from vall database file "filename"
	void
	read_file( std::string const & filename );

	/// read in one more line from Vall input file
	void
	add_line(
		const char sq,
		const char ss,
		const core::Real ph,
		const core::Real ps,
		const core::Real om
	)
	{
		runtime_assert( ss == 'H' || ss == 'E' || ss == 'L' );
		// chain info
		if ( phi_.empty() ) {
			chain_.push_back( 1 );
		} else {
			core::Real const prev_psi( psi_.back() ), prev_omega( omega_.back() );
			if ( ( std::abs( prev_psi ) + std::abs( prev_omega ) ) < 0.01 ) {
				chain_.push_back( chain_.back() + 1 );
			} else {
				chain_.push_back( chain_.back() );
			}
		}

		sequence_.push_back( sq );
		secstruct_.push_back( ss );

		phi_.push_back( ph );
		psi_.push_back( ps );
		omega_.push_back( om );

		bigbin_.push_back( torsion2big_bin( ph, ps, om ) );
	}

	utility::vector1< char > const & sequence () const {return sequence_; }
	utility::vector1< char > const & secstruct() const {return secstruct_;}
	utility::vector1< char > const & bigbin() const {return bigbin_;}

	utility::vector1< core::Real > const & phi  () const {return phi_;}
	utility::vector1< core::Real > const & psi  () const {return psi_;}
	utility::vector1< core::Real > const & omega() const {return omega_;}

	utility::vector1< core::Size > const & chain() const {return chain_;}

	//  inline
	//  bool
	//  is_lower_terminus( core::Size const pos ) const
	//  {
	//   return ( ( pos == 1 ) ||
	//        ( std::abs( psi_[ pos-1 ] ) + std::abs( omega_[ pos-1 ] ) + std::abs( phi_[ pos   ] ) < core::Real(0.01) ) );
	//  }

	//  inline
	//  bool
	//  is_upper_terminus( core::Size const pos ) const
	//  {
	//   return ( ( pos == size() ) ||
	//        ( std::abs( psi_[ pos   ] ) + std::abs( omega_[ pos   ] ) + std::abs( phi_[ pos+1 ] ) < 0.01 ) );
	//  }

	//  inline
	//  bool
	//  is_terminus( core::Size const pos ) const
	//  {
	//   return ( is_lower_terminus( pos ) || is_upper_terminus( pos ) );
	//  }


	/// number of lines in Vall database
	int size() const { return sequence_.size(); }

	/// number of chains
	core::Size
	num_chains() const
	{
		return chain_.back();
	}

	// pick fragments for a single residue position from vall database
	void
	get_frags(
		core::Size const nfrags,
		std::string const & target_seq,
		std::string const & target_ss,
		core::Real const seq_weight,
		core::Real const ss_weight,
		bool const exclude_gly,
		bool const exclude_pro,
		bool const exclude_cis_peptides, // at non-pre-Pro positions
		utility::vector1< core::Size > const & homs_to_exclude,
		SingleResidueTorsionFragmentLibrary & library,
		core::Real const bb_weight = 0.0,
		std::string const & target_bb = std::string()
	) const;

	// pick fragments for a single residue position from vall database
	void
	get_frags(
		core::Size const nfrags,
		std::string const & target_seq,
		utility::vector1< std::map< char, core::Real > > const & target_ss,
		core::Real const seq_weight,
		core::Real const ss_weight,
		bool const exclude_gly,
		bool const exclude_pro,
		bool const exclude_cis_peptides, // at non-pre-Pro positions
		utility::vector1< core::Size > const & homs_to_exclude,
		SingleResidueTorsionFragmentLibrary & library,
		core::Real const bb_weight = 0.0,
		std::string const & target_bb = std::string()
	) const;

	void
	get_cheating_frags(
		core::Size const nfrags,
		std::string const & target_seq,
		std::string const & target_ss,
		utility::vector1< core::Real > const & target_phi,
		utility::vector1< core::Real > const & target_psi,
		utility::vector1< core::Real > const & target_omega,
		core::Real const seq_weight,
		core::Real const ss_weight,
		core::Real const torsion_weight,
		core::Real const min_torsion_dev,
		core::Real const max_torsion_dev,
		utility::vector1< core::Size > const & homs_to_exclude,
		SingleResidueTorsionFragmentLibrary & library
	) const;

private:
	utility::vector1< char > sequence_;
	utility::vector1< char > secstruct_;
	utility::vector1< char > bigbin_;

	utility::vector1< core::Real > phi_;
	utility::vector1< core::Real > psi_;
	utility::vector1< core::Real > omega_;

	utility::vector1< core::Size > chain_;

};


/// handles loading the vall if necessary
void
get_frags(
	core::Size const nfrags,
	std::string const & target_seq,
	std::string const & target_ss,
	core::Real const seq_weight,
	core::Real const ss_weight,
	bool const exclude_gly,
	bool const exclude_pro,
	bool const exclude_cis_peptides,
	utility::vector1< core::Size > const & homs_to_exclude,
	SingleResidueTorsionFragmentLibrary & library,
	core::Real const bb_weight = 0.0,
	std::string const & target_bb = std::string()
);

/// handles loading the vall if necessary
void
get_frags(
	core::Size const nfrags,
	std::string const & target_seq,
	utility::vector1< std::map< char, core::Real > > const & target_ss, // HEL
	core::Real const seq_weight,
	core::Real const ss_weight,
	bool const exclude_gly,
	bool const exclude_pro,
	bool const exclude_cis_peptides,
	utility::vector1< core::Size > const & homs_to_exclude,
	SingleResidueTorsionFragmentLibrary & library,
	core::Real const bb_weight = 0.0,
	std::string const & target_bb = std::string()
);

void
get_cheating_frags(
	core::Size const nfrags,
	std::string const & target_seq,
	std::string const & target_ss,
	utility::vector1< core::Real > const & target_phi,
	utility::vector1< core::Real > const & target_psi,
	utility::vector1< core::Real > const & target_omega,
	core::Real const seq_weight,
	core::Real const ss_weight,
	core::Real const torsion_weight,
	core::Real const min_torsion_dev,
	core::Real const max_torsion_dev,
	utility::vector1< core::Size > const & homs_to_exclude,
	SingleResidueTorsionFragmentLibrary & library
);


void
dump_vall_fasta( std::string const & fasta_filename );

} // ns frags
} // ns protocols
#endif
