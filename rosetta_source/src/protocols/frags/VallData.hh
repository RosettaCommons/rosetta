// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_frags_VallData_hh
#define INCLUDED_protocols_frags_VallData_hh


// Rosetta Headers

#include <core/types.hh>
#include <protocols/frags/VallData.fwd.hh>
#include <protocols/frags/TorsionFragment.fwd.hh>

// ObjexxFCL Headers
// AUTO-REMOVED #include <utility/vector1.hh>

// C++ Headers
// #include <cmath>
// #include <cstdlib>
// #include <iostream>
// #include <fstream>
// #include <sstream>
#include <string>

//Auto Headers
#include <utility/vector1_bool.hh>


namespace protocols {
namespace frags {

using core::Real;
using core::Size;

class VallData {
public:
	/// default constructor
	VallData ()
	{
		Size const big_size( 100000 ); // rough guess
		sequence_.reserve( big_size );
		secstruct_.reserve( big_size );
		phi_.reserve( big_size );
		psi_.reserve( big_size );
		omega_.reserve( big_size );
	}

	/// constructor from input vall database file
	VallData ( std::string const & filename )
	{
		Size const big_size( 100000 ); // rough guess
		// prevent lots of redimensioning as we read file? does this even matter?
		sequence_.reserve( big_size );
		secstruct_.reserve( big_size );
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
		const Real ph,
		const Real ps,
		const Real om
	)
	{
		sequence_.push_back( sq );
		secstruct_.push_back( ss );

		phi_.push_back( ph );
		psi_.push_back( ps );
		omega_.push_back( om );
	}

	utility::vector1< char > const & sequence () const {return sequence_; }
	utility::vector1< char > const & secstruct() const {return secstruct_;}

	utility::vector1< Real > const & phi  () const {return phi_;}
	utility::vector1< Real > const & psi  () const {return psi_;}
	utility::vector1< Real > const & omega() const {return omega_;}

	/// number of lines in Vall database
	int size() const { return sequence_.size(); }

	// pick fragments for a single residue position from vall database
	void
	get_frags(
		Size const nfrags,
		std::string const & target_seq,
		std::string const & target_ss,
		Real const seq_weight,
		Real const ss_weight,
		bool const exclude_gly,
		bool const exclude_pro,
		bool const exclude_cys_peptides,
		SingleResidueTorsionFragmentLibrary & library
	) const;

private:
	utility::vector1< char > sequence_;
	utility::vector1< char > secstruct_;

	utility::vector1< Real > phi_;
	utility::vector1< Real > psi_;
	utility::vector1< Real > omega_;

};

} // ns frags
} // ns protocols

#endif
