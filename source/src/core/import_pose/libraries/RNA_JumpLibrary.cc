// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Library of RNA base pair jumps
/// @brief protocols that are specific to RNA_JumpLibrary
/// @details
/// @author Rhiju Das


// Rosetta Headers
#include <core/import_pose/libraries/RNA_JumpLibrary.hh>

#include <utility>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

#include <numeric/random/random.hh>
#include <core/chemical/ResidueType.hh>

#include <core/types.hh>
#include <core/chemical/ResidueType.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <fstream>
#include <iostream>

#include <utility/vector1.hh>

using namespace core::chemical::rna;

namespace core {
namespace import_pose {
namespace libraries {

/// @details Auto-generated virtual destructor
RNA_JumpLibrary::~RNA_JumpLibrary() = default;

/// @details Auto-generated virtual destructor
RNA_PairingTemplate::~RNA_PairingTemplate() = default;


static basic::Tracer tr( "protocols.rna.denovo.libraries.RNA_JumpLibrary" );

std::string
BasePairType::tag() const {
	std::string tag;
	tag += aa1;
	tag += aa2;
	tag += get_edge_from_num( edge1 );
	tag += get_edge_from_num( edge2 );
	tag += get_orientation_from_num( orientation );
	return tag;
}


//@brief constructor
RNA_PairingTemplate::RNA_PairingTemplate(
	core::kinematics::Jump const & j,
	std::string const & atom_name1,
	std::string const & atom_name2 ):
	jump_forward_ ( j ),
	jump_backward_( j ),
	atom_name1_( atom_name1 ),
	atom_name2_( atom_name2 )
{}

//@brief constructor
RNA_PairingTemplate::RNA_PairingTemplate(
	core::kinematics::Jump const & j1,
	core::kinematics::Jump const & j2,
	std::string const & atom_name1,
	std::string const & atom_name2 ):
	jump_forward_ ( j1 ),
	jump_backward_( j2 ),
	atom_name1_( atom_name1 ),
	atom_name2_( atom_name2 )
{}

/// @brief constructor
RNA_JumpLibrary::RNA_JumpLibrary( std::string const & filename ):
	jump_library_filename_( filename )
{
	//read_jumps_from_file();
}

/////////////////////////////////////////////////////////////////////////////////////
/// @details following is not actually const, but involves lazy loading to a mutable map.
void
RNA_JumpLibrary::read_jumps_from_file() const
{
	utility::io::izstream data( jump_library_filename_ );
	tr << "Reading RNA jump library: " << jump_library_filename_ << std::endl;

	if ( data.fail() )  {
		utility_exit_with_message(  "Bad jumpdata file? " + jump_library_filename_ );
	}

	std::string line;
	std::string tag;//,filename;
	Size pos1,pos2;
	char reschar1,reschar2,edgechar1,edgechar2,orientation;
	core::kinematics::Jump jump1, jump2;

	char tmpbuf[ 4 ];

	while ( getline( data,line ) ) {
		std::istringstream is( line );

		is >> tag >>
			pos1 >> edgechar1 >> pos2 >> edgechar2 >>
			orientation >> reschar1 >> reschar2;

		//  std::cout << "READ IN: " << reschar1 << " " << reschar2 << std::endl;

		is.read( tmpbuf, 1 );
		is.read( tmpbuf, 4 );
		std::string atom_name1 = "";
		for ( char n : tmpbuf ) atom_name1 += n;

		is.read( tmpbuf, 1 );
		is.read( tmpbuf, 4 );
		std::string atom_name2 = "";
		for ( char n : tmpbuf ) atom_name2 += n;

		//  std::cout << "READ IN ATOM NAMES? " << atom_name1 << " " << atom_name2 << std::endl;

		if ( is.fail() || tag != "PAIR" ) continue;

		is >> jump1;

		if ( edgechar1 == 'P' || edgechar2 == 'P' ) {
			is >> jump2;
		} else {
			jump2 = jump1;
		}

		save_in_jump_library( reschar1, reschar2,
			edgechar1, edgechar2,
			orientation,
			atom_name1, atom_name2,
			jump1, jump2 );

		jump1.reverse();
		jump2.reverse();

		save_in_jump_library( reschar2, reschar1,
			edgechar2, edgechar1,
			orientation,
			atom_name2, atom_name1,
			jump1, jump2 );

	}

	runtime_assert( !rna_pairing_template_map_.empty() );
}

/////////////////////////////////////////////////////////////////////////////////////
// This is a lot of work just to get base-phosphate jumps...
void
RNA_JumpLibrary::check_forward_backward(
	std::string & atom_name,
	bool const forward,
	core::kinematics::Jump & j,
	RNA_PairingTemplateOP const & t) const
{
	if ( atom_name == " ?  " ) {
		if ( forward ) {
			atom_name = " O5'";
		} else {
			atom_name = " O3'";
			j = t->jump_backward();
		}
	}
}

using namespace core;

char one_letter_from_rt( chemical::ResidueType const & rt ) {
	if ( rt.aa() == chemical::na_rad || rt.aa() == chemical::na_ade || rt.na_analogue() == chemical::na_rad ) return 'a';
	if ( rt.aa() == chemical::na_rcy || rt.aa() == chemical::na_cyt || rt.na_analogue() == chemical::na_rcy ) return 'c';
	if ( rt.aa() == chemical::na_rgu || rt.aa() == chemical::na_gua || rt.na_analogue() == chemical::na_rgu ) return 'g';
	if ( rt.aa() == chemical::na_ura || rt.na_analogue() == chemical::na_ura ) return 'u';
	if ( rt.name1() == 't' ) return 'u';

	utility_exit_with_message( "Sorry, but " + rt.name() + " is not yet accounted for by RNA_JumpLibrary!" );
}

/////////////////////////////////////////////////////////////////////////////////////
core::kinematics::Jump
RNA_JumpLibrary::get_random_base_pair_jump(
	core::chemical::ResidueType const & rt1,
	core::chemical::ResidueType const & rt2,
	BaseEdge const edge1,
	BaseEdge const edge2,
	BaseDoubletOrientation const orientation,
	std::string & atom_name1,
	std::string & atom_name2,
	bool & success,
	bool const forward1 /*= true*/,
	bool const forward2 /*= true*/ ) const
{

	if ( rna_pairing_template_map_.empty() ) read_jumps_from_file();

	// key for looking up the template geometry:
	char const aa1 = one_letter_from_rt( rt1 );
	char const aa2 = one_letter_from_rt( rt2 );

	BasePairType key( aa1, aa2, edge1, edge2, orientation );
	Size ntemplates = 0;

	// tr << "Looking for: " << aa1 << ' ' << aa2 << ' ' << edge1 << ' ' << edge2 << ' ' << orientation << std::endl;

	if ( rna_pairing_template_map_.find( key ) == rna_pairing_template_map_.end() ) {

		tr << "Can't seem to find a pairing inside database with aa1: " <<  aa1 << " aa2: " << aa2 << " edge1: " << edge1 << " edge2: " << edge2 << " orientation: " << orientation << std::endl;

		//  utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		success = false;
		atom_name1 = "XXXX";
		atom_name2 = "XXXX";
		return core::kinematics::Jump(); //default garbage jump.
	}

	RNA_PairingTemplateList const & templates( rna_pairing_template_map_.find( key )->second );

	ntemplates = templates.size();

	if ( ntemplates < 1 ) {

		std::cout << "Can't seem to find a pairing inside database with aa1: " <<  aa1 << " aa2: " << aa2 << " edge1: " << edge1 << " edge2: " << edge2 << " orientation: " << orientation << std::endl;

		std::cerr << "Can't seem to find a pairing inside database with aa1: " <<  aa1 << " aa2: " << aa2 << " edge1: " << edge1 << " edge2: " << edge2 << " orientation: " << orientation << std::endl;

		//  utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		success = false;
		atom_name1 = "XXXX";
		atom_name2 = "XXXX";
		return core::kinematics::Jump(); //default garbage jump.
	}

	Size const index( static_cast<Size>( numeric::random::rg().uniform() * ntemplates )  + 1 );

	core::kinematics::Jump j ( templates[ index ]->jump() );
	atom_name1 = templates[ index ]->atom_name1();
	atom_name2 = templates[ index ]->atom_name2();

	//Currently if we connect to a main chain atom -- e.g., phosphate,
	// the jump transformation depends on whether we're going forward or backward.
	// This is encoded as an atom name of "?"
	// Following seems a little inelegant right now.
	// Also, currently only allows only one of the partners to
	//  be on the RNA backbone. Again, inelegant.
	// runtime_assert( atom_name1 != " ?  " || atom_name2 != " ?  " );
	// check_forward_backward( atom_name1, forward1, j, templates[ index ] );
	// check_forward_backward( atom_name2, forward2, j, templates[ index ] );
	runtime_assert( edge1 != PHOSPHATE || edge2 != PHOSPHATE );
	if ( edge1 == PHOSPHATE && !forward1 ) j = templates[ index ]->jump_backward();
	if ( edge2 == PHOSPHATE && !forward2 ) j = templates[ index ]->jump_backward();

	success = true;

	// tr << "Found     it?" << j << ' ' << success << std::endl;

	return j;
}

/////////////////////////////////////////////////////////////////////////////////////
void
RNA_JumpLibrary::save_in_jump_library(
	Size const reschar1, Size const reschar2,
	char const edgechar1, char const edgechar2,
	char const orientationchar,
	std::string const & atom_name1,
	std::string const & atom_name2,
	core::kinematics::Jump const & jump1,
	core::kinematics::Jump const & jump2
) const {
	RNA_PairingTemplateOP p( new RNA_PairingTemplate( jump1, jump2,  atom_name1, atom_name2) );

	BaseEdge edge1( get_edge_from_char( edgechar1 ) );
	BaseEdge edge2( get_edge_from_char( edgechar2 ) );
	BaseDoubletOrientation orientation( get_orientation_from_char( orientationchar ) );

	// Save to template lists with wildcards ('X') too.
	for ( Size i = 0; i <= 1; i++ ) {

		BaseEdge const edge1_temp = i ? edge1 : ANY_BASE_EDGE;

		for ( Size j = 0; j <= 1; j++ ) {

			BaseEdge const edge2_temp = j ? edge2 : ANY_BASE_EDGE;

			for ( Size k = 0; k <= 1; k++ ) {

				BaseDoubletOrientation const orientation_temp = k ? orientation : ANY_BASE_DOUBLET_ORIENTATION;

				BasePairType base_pair_type( reschar1, reschar2, edge1_temp, edge2_temp, orientation_temp);
				rna_pairing_template_map_[ base_pair_type ].push_back( p );
			}
		}
	}
}


/////////////////////////////////////////////////////////////////////////////////////
bool
RNA_JumpLibrary::has_template( BasePairType const & base_pair_type ) const {
	if ( rna_pairing_template_map_.empty() ) read_jumps_from_file();
	return  ( rna_pairing_template_map_.find( base_pair_type ) != rna_pairing_template_map_.end() );
}

/////////////////////////////////////////////////////////////////////////////////////
RNA_PairingTemplateList const &
RNA_JumpLibrary::rna_pairing_template_map( BasePairType const & base_pair_type) const {
	runtime_assert( has_template( base_pair_type ) );
	return ( rna_pairing_template_map_.find( base_pair_type )->second );
}

} //libraries
} //denovo
} //protocols
