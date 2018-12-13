// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/relax/membrane/util.hh
/// @brief      utility functions for the MPMutateRelaxMover
/// @author     JKLeman (julia.koehler1982@gmail.com)

// Unit Headers
#include <protocols/moves/Mover.hh>

// Project Headers

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/io/util.hh>
#include <core/pose/util.hh>
#include <core/types.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "protocols.relax.membrane.util" );

namespace protocols {
namespace relax {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::moves;

////////////////////////////////////////////////////////////////////////////////
/// @brief Add mutants to private data: A163F into vectors
/// @details This is an entire line of the mutant input file
void add_mutant_to_vectors(
	Pose & pose,
	std::string mutations,
	utility::vector1< utility::vector1< char > > & wt_res,
	utility::vector1< utility::vector1< core::Size > > & resid,
	utility::vector1< utility::vector1< char > > & mut_res )
{

	using namespace utility;
	TR << "adding mutants to vectors" << std::endl;

	// initialize line vectors
	utility::vector1< char > wildtypes;
	utility::vector1< core::Size > seqids;
	utility::vector1< char > mutants;

	// initialize single point variables
	char wt, mut;
	core::Size seqid;

	// PDB numbering?
	bool pdb_numbering( false );

	// split string by whitespace, write output into vector
	utility::vector1< std::string > all_mutants = split_whitespace( mutations );

	// iterate over mutants
	for ( core::Size col = 1; col <= all_mutants.size(); ++col ) {

		// get single mutant
		std::string wt_id_mut = all_mutants[ col ];

		// PDB or pose numbering?
		// PDB numbering: A_E15L, chain, underscore, wt, resn, mutant
		// pose numbering: E15L, wt, resn, mutant
		// if string contains underscore, then PDB numbering
		if ( col == 1 ) {
			if ( wt_id_mut.find( "_" ) != std::string::npos ) {
				TR << "Mutations are in PDB numbering scheme." << std::endl;
				pdb_numbering = true;
			} else {
				TR << "Mutations are in pose numbering scheme." << std::endl;
			}
		}

		// check input file format
		if ( ! pdb_numbering && wt_id_mut.find( "_" ) != std::string::npos ) {
			utility_exit_with_message( "Mutations should either be in PDB numbering (A_E15L) where A is the chain, or in pose numbering (E15L) where the residues are renumbered without gaps, starting from 1. Quitting." );
		} else if ( pdb_numbering && wt_id_mut.find( "_" ) == std::string::npos ) {
			utility_exit_with_message( "Mutations should either be in PDB numbering (A_E15L) where A is the chain, or in pose numbering (E15L) where the residues are renumbered without gaps, starting from 1. Quitting." );
		}

		// pose numbering
		if ( pdb_numbering == false ) {

			// get wt and mutant from vector: split string by character
			wt = char( wt_id_mut[ 0 ] );
			mut = char( wt_id_mut[ wt_id_mut.size()-1 ] );

			// get residue number: get each digit
			utility::vector1< std::string > tmp;
			for ( core::Size i = 1; i <= wt_id_mut.size()-2; ++i ) {
				tmp.push_back( to_string( wt_id_mut[ i ] ) );
			}

			// join digits together to the residue number
			seqid = string2Size( join( tmp, "" ) );

		} else {
			// PDB numbering

			// get chain
			char chain = wt_id_mut[ 0 ];

			// remove the chain and the underscore from the string
			wt_id_mut.erase( 0, 2 );

			// get wt and mutant from vector: split string by character
			wt = char( wt_id_mut[ 0 ] );
			mut = char( wt_id_mut[ wt_id_mut.size()-1 ] );

			// get residue number: get each digit
			utility::vector1< std::string > tmp;
			for ( core::Size i = 1; i <= wt_id_mut.size()-2; ++i ) {
				tmp.push_back( to_string( wt_id_mut[ i ] ) );
			}

			// join digits together to the residue number
			int seqid_pdb = string2int( join( tmp, "" ) );

			// get pose resnumber from pdb numbering
			seqid = pose.pdb_info()->pdb2pose( chain, seqid_pdb );
		}

		TR << "wt " << wt << ", seqid " << seqid << ", mut " << mut << std::endl;

		// put all of the wt / seqid / mut info into the line vectors
		wildtypes.push_back( wt );
		seqids.push_back( seqid );
		mutants.push_back( mut );

	} // iterate over columns in line

	// add the line to the file (or outer vector in this case)
	wt_res.push_back( wildtypes );
	resid.push_back( seqids );
	mut_res.push_back( mutants );

} // add mutants to constructs vector

////////////////////////////////////////////////////////////////////////////////
/// @brief Check mutant file for errors
/// @details If Rosetta doesn't start crying, you're good to go
bool check_mutants_ok(
	Pose & pose,
	utility::vector1< utility::vector1< char > > wt_res,
	utility::vector1< utility::vector1< core::Size > > resid )
{
	using namespace utility;
	TR << "checking mutant file" << std::endl;

	// go through constructs
	for ( core::Size c = 1; c <= wt_res.size(); ++c ) {

		// go through mutants
		for ( core::Size m = 1; m <= wt_res[ c ].size(); ++m ) {

			// check whether wt residue has the same identity in the pose
			if ( wt_res[ c ][ m ] != pose.residue_type( resid[ c ][ m ] ).name1() ) {
				TR << "c " << c << ", m " << m << std::endl;
				TR << "residue identity in PDB " << to_string( pose.residue_type( m ).name1() ) << ", " << wt_res[ c ][ m ] << std::endl;
				TR << "Residue identity in input file doesn't match the pose!" << std::endl;
				return false;
			}

		} // inner loop
	} // outer loop
	return true;

} // check mutant file


} // membrane
} // relax
} // protocols
