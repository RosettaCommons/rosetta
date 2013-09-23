// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ChemicalShiftSequence.cc
/// @author James Thompson

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/ChemicalShiftSequence.hh>

#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <utility/pointer/owning_ptr.hh>

// AUTO-REMOVED #include <core/chemical/AA.hh>

#include <map>
#include <iostream>
#include <string>

#include <ObjexxFCL/format.hh>

#include <utility/vector1.hh>



namespace core {
namespace sequence {

static basic::Tracer tr( "core.sequence.ChemicalShiftSequence" );

void ChemicalShiftSequence::read_from_file(
	utility::file::FileName const & fn
) {
	using core::Size;
	using core::Real;
	using std::string;
	using utility::vector1;

	// order of amino acids in the 2nd_inCS.tab file
	// C
	// CA
	// CB
	// HA
	// HN
	// N
	static bool init_order(false);
	static std::map< string, Size > order;
	static Real const bad_shift  ( 9999.000 );
	static Real const upper_limit( 20.0 );

	if ( !init_order ) {
		order[ "C"  ] = 1;
		order[ "CA" ] = 2;
		order[ "CB" ] = 3;
		order[ "HA" ] = 4;
		order[ "HN" ] = 5;
		order[ "N"  ] = 6;
		//avg_vsigma[ "C"  ] = 0.948093982461;
		//avg_vsigma[ "CA" ] = 0.828795509906;
		//avg_vsigma[ "CB" ] = 0.969003290496;
		//avg_vsigma[ "HA" ] = 0.231641960268;
		//avg_vsigma[ "HN" ] = 0.725943902221;
		//avg_vsigma[ "N"  ] = 2.78164308475;
		alphabet_.push_back( "C"  );
		alphabet_.push_back( "CA" );
		alphabet_.push_back( "CB" );
		alphabet_.push_back( "HA" );
		alphabet_.push_back( "HN" );
		alphabet_.push_back( "N"  );
		init_order = true;
	}

	string aa_seq;
	utility::io::izstream input(fn);
	id(fn.name());

	if ( !input ) {
		string msg( "ERROR: Unable to open file " + fn.name() );
		utility_exit_with_message(msg);
	}

	//std::cout << "reading data from " << fn.name() << std::endl;

	std::string line;
	while ( line.substr(0,13) != "DATA SEQUENCE" ) {
		getline( input, line );
	}

	while ( line.substr(0,13) == "DATA SEQUENCE" ) {
		std::string dummy;
	 	std::istringstream ls( line );
		ls >> dummy >> dummy;
		while ( ls.good() ) {
			ls >> dummy;
			aa_seq += dummy;
		}
		getline(input,line);
	}
	for ( Size ii = 1; ii <= 3; ++ii ) {
		getline( input, line );
	}

	vector1< Real > prof_row( order.size(), bad_shift );
	//vector1< vector1< Real > > new_prof(
	//	aa_seq.size(), prof_row
	//	//( order.size(), bad_shift )
	//);

	vector1< vector1< Real > > new_prof;

	//std::cout << "SEQUENCE: " << aa_seq << std::endl;
	Size last_res_no(0);
	aa_seq = "";
	while( getline( input, line ) ) {
	 	std::istringstream ls( line );
		core::Size res_no;
		std::string atm_name;
		char aa;
		core::Real chemical_shift;
		ls >> res_no >> aa >> atm_name >> chemical_shift;

		if ( res_no != last_res_no ) {
			for ( Size ii = last_res_no; ii <= res_no; ++ii ) {
				new_prof.push_back(
					vector1< Real >( order.size(), bad_shift )
				);
			}
			aa_seq += aa;
		}

		last_res_no = res_no;

		if ( chemical_shift > upper_limit ) {
			chemical_shift = bad_shift;
		}

		std::cout << line << std::endl;
		std::cout << "res_no = " << res_no << std::endl;
		std::cout << "aa = " << aa << std::endl;
		std::cout << "aa_seq[" << res_no-1 << "] = " << aa_seq[res_no-1]
			<< std::endl;
		Size const idx( order[atm_name] );
		std::cout << "idx = " << idx << std::endl;
		std::cout << std::endl << std::endl;

		if ( idx != 0 ) new_prof[ res_no ][ idx ] = chemical_shift;
		//runtime_assert( aa_seq[res_no-1] == aa );
	}
	input.close();
	sequence(aa_seq);
	//std::cout << "sequence = " << sequence() << std::endl;
	//std::cout << "finished reading." << std::endl;

	sequence( aa_seq );
	profile( new_prof );

	check_internals_();
}

Real ChemicalShiftSequence::sigma(
	Size /*pos*/, Size idx
) {
	// hacky - assumed order of atoms:
	if ( idx == 1 ) {
		return 0.948093982461;
	} else if ( idx == 2 ) {
		return 0.828795509906;
	} else if ( idx == 3 ) {
		return 0.969003290496;
	} else if ( idx == 4 ) {
		return 0.231641960268;
	} else if ( idx == 5 ) {
		return 0.725943902221;
	} else if ( idx == 6 ) {
		return 2.78164308475;
	} else {
		utility_exit_with_message("Error: invalid index!");
	}
	return 0.0;
}

std::ostream & operator<<(
	std::ostream & out, const ChemicalShiftSequence & p
) {
	Size width     = 8;
	Size precision = 3;

	out << p.to_string() << std::endl;
	for ( Size i = 1; i <= p.length(); ++i ) {
		for ( Size j = 1; j <= p.width(); ++j ) {
			out << ObjexxFCL::format::F( width, precision, p.prof_row(i)[j] );
		}
		out << std::endl;
	}

	return out;
}

void ChemicalShiftSequence::check_internals_() const {
	using core::Real;
	using utility::vector1;

	runtime_assert( profile().size() == length() );

	for ( Size i = 1; i <= profile().size(); ++i ) {
		//std::cout << profile()[i].size() << " " << alphabet_.size() << std::endl;
		//runtime_assert( profile()[i].size() == alphabet_.size() );
		runtime_assert( profile()[i].size() == 6 );
	}
}


} // sequence
} // core
