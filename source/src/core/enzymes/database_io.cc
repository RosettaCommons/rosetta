// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/enzymes/database_io.hh
/// @brief   Database input/output function definitions for enzyme data.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/enzymes/database_io.hh>

// Utility header
#include <utility/io/util.hh>

// C++ header
#include <sstream>

// Basic header
#include <basic/Tracer.hh>


// Construct tracer.
static basic::Tracer TR( "core.enzymes.database_io" );


namespace core {
namespace enzymes {

EnzymeData
read_enzyme_data_from_file( std::string const & filename )
{
	using namespace std;
	using namespace utility;

	vector1< string > const lines( io::get_lines_from_file_data( filename ) );
	EnzymeData enzyme_data;

	// The first (non-comment) line contains all of the information except for the co-substrate list.
	istringstream line_word_by_word( lines[ 1 ] );

	string consensus_sequence;
	string cs_type;
	uint cs_resnum;  // the position in the consensus sequence that is potentially modified
	string atom_to_modify;
	Real efficiency;  // ratio of times this enzyme acts on a recognized site

	line_word_by_word >> consensus_sequence >> cs_type >> cs_resnum >> atom_to_modify >> efficiency;

	enzyme_data.consensus_sequence = consensus_sequence;
	if ( cs_type == "AA" ) {
		enzyme_data.cs_type = AA;
	} else if ( cs_type == "NA" ) {
		enzyme_data.cs_type = NA;
	} else if ( cs_type == "SACCHARIDE" ) {
		enzyme_data.cs_type = SACCHARIDE;
	} else {
		TR.Error << filename << " contains an invalid consensus sequence type: assuming peptide sequence." << endl;
		enzyme_data.cs_type = AA;
	}
	enzyme_data.cs_resnum = cs_resnum;
	enzyme_data.atom_to_modify = atom_to_modify;
	enzyme_data.efficiency = efficiency;

	// All the remaining lines are co-substrates or by-products, if any.
	Size const n_lines( lines.size() );
	for ( uint i( 2 ); i <= n_lines; ++i ) {  // Start at line 2.
		enzyme_data.second_substrates_or_byproducts.push_back( lines[ i ] );
	}

	return enzyme_data;
}

}  // namespace enzymes
}  // namespace core
