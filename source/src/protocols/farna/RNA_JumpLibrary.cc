// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Library of RNA base pair jumps
/// @brief protocols that are specific to RNA_JumpLibrary
/// @detailed
/// @author Rhiju Das


// Rosetta Headers
#include <protocols/farna/RNA_JumpLibrary.hh>

#include <utility/io/izstream.hh>
#include <utility/exit.hh>

#include <numeric/random/random.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <fstream>
#include <iostream>

#include <utility/vector1.hh>



namespace protocols {
namespace farna {

/// @details Auto-generated virtual destructor
RNA_JumpLibrary::~RNA_JumpLibrary() {}

/// @details Auto-generated virtual destructor
RNA_PairingTemplate::~RNA_PairingTemplate() {}


static thread_local basic::Tracer tr( "protocols.rna.RNA_JumpLibrary" );

//@brief constructor
RNA_PairingTemplate::RNA_PairingTemplate(
  core::kinematics::Jump const j,
  std::string const atom_name1,
  std::string const atom_name2 ):
  jump_forward_ ( j ),
  jump_backward_( j ),
  atom_name1_( atom_name1 ),
  atom_name2_( atom_name2 )
{}

//@brief constructor
RNA_PairingTemplate::RNA_PairingTemplate(
  core::kinematics::Jump const j1,
  core::kinematics::Jump const j2,
  std::string const atom_name1,
  std::string const atom_name2 ):
  jump_forward_ ( j1 ),
  jump_backward_( j2 ),
  atom_name1_( atom_name1 ),
  atom_name2_( atom_name2 )
{}

/////////////////////////////////////////////////////////////////////////////////////
void
RNA_JumpLibrary::read_jumps_from_file( std::string const & jump_library_filename )
{

  utility::io::izstream data( jump_library_filename );
  tr << "Reading RNA jump library: " << jump_library_filename << std::endl;

	if ( data.fail() )  {
		utility_exit_with_message(  "Bad jumpdata file? " + jump_library_filename );
	}

  std::string line;
  std::string tag,filename;
  Size pos1,pos2;
  char reschar1,reschar2,edgechar1,edgechar2,orientation;
  core::kinematics::Jump jump1, jump2;

  char tmpbuf[ 4 ];

  while ( getline( data,line ) ) {
    std::istringstream is( line );

    is >> tag >>
      pos1 >> edgechar1 >> pos2 >> edgechar2 >>
      orientation >> reschar1 >> reschar2;

    //		std::cout << "READ IN: " << reschar1 << " " << reschar2 << std::endl;

    is.read( tmpbuf, 1 );
    is.read( tmpbuf, 4 );
    std::string atom_name1 = "";
    for (int n=0; n <4;n++) atom_name1 += tmpbuf[n];

    is.read( tmpbuf, 1 );
    is.read( tmpbuf, 4 );
    std::string atom_name2 = "";
    for (int n=0; n <4;n++) atom_name2 += tmpbuf[n];

    //		std::cout << "READ IN ATOM NAMES? " << atom_name1 << " " << atom_name2 << std::endl;

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

}

/////////////////////////////////////////////////////////////////////////////////////
// This is a lot of work just to get base-phosphate jumps...
void
RNA_JumpLibrary::check_forward_backward(
   std::string & atom_name,
	 bool const forward,
	 core::kinematics::Jump & j,
	 protocols::farna::RNA_PairingTemplateOP const & t) const
{
	if (atom_name == " ?  " ) {
		if ( forward ) {
			atom_name = " O5'";
		} else {
			atom_name = " O3'";
			j = t->jump_backward();
		}
	}
}


/////////////////////////////////////////////////////////////////////////////////////
core::kinematics::Jump
RNA_JumpLibrary::get_random_base_pair_jump(
    char const aa1,
    char const aa2,
    char const edge1,
    char const edge2,
    char const orientation,
    std::string & atom_name1,
    std::string & atom_name2,
    bool & success,
		bool const forward1 /*= true*/,
		bool const forward2 /*= true*/ ) const
{
  // key for looking up the template geometry:
  BasePairType key( aa1, aa2, edge1, edge2, orientation);
  Size ntemplates = 0;

	//	tr << "Looking for: " << aa1 << ' ' << aa2 << ' ' << edge1 << ' ' << edge2 << ' ' << orientation << std::endl;

  if( rna_pairing_template_map_.find( key ) == rna_pairing_template_map_.end() ){

    std::cout << "Can't seem to find a pairing inside database with aa1: " <<  aa1 << " aa2: " << aa2 << " edge1: " << edge1 << " edge2: " << edge2 << " orientation: " << orientation << std::endl;

    std::cerr << "Can't seem to find a pairing inside database with aa1: " <<  aa1 << " aa2: " << aa2 << " edge1: " << edge1 << " edge2: " << edge2 << " orientation: " << orientation << std::endl;

    //		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
    success = false;
    atom_name1 = "XXXX";
    atom_name2 = "XXXX";
    return core::kinematics::Jump(); //default garbage jump.

  }

  RNA_PairingTemplateList const & templates( rna_pairing_template_map_.find( key )->second );

  ntemplates = templates.size();

  if( ntemplates < 1 ){

    std::cout << "Can't seem to find a pairing inside database with aa1: " <<  aa1 << " aa2: " << aa2 << " edge1: " << edge1 << " edge2: " << edge2 << " orientation: " << orientation << std::endl;

    std::cerr << "Can't seem to find a pairing inside database with aa1: " <<  aa1 << " aa2: " << aa2 << " edge1: " << edge1 << " edge2: " << edge2 << " orientation: " << orientation << std::endl;

    //		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
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
	//	runtime_assert( atom_name1 != " ?  " || atom_name2 != " ?  " );
	//	check_forward_backward( atom_name1, forward1, j, templates[ index ] );
	//	check_forward_backward( atom_name2, forward2, j, templates[ index ] );
	runtime_assert( edge1 != 'P' || edge2 != 'P' );
	if ( edge1 == 'P' && !forward1 ) j = templates[ index ]->jump_backward();
	if ( edge2 == 'P' && !forward2 ) j = templates[ index ]->jump_backward();

  success = true;

	// tr << "Found     it?" << j << ' ' << success << std::endl;

  return j;
}

/////////////////////////////////////////////////////////////////////////////////////
void
RNA_JumpLibrary::save_in_jump_library(
   Size const reschar1, Size const reschar2,
   char const edgechar1, char const edgechar2,
   char const orientation,
   std::string const & atom_name1,
   std::string const & atom_name2,
   core::kinematics::Jump const & jump1,
   core::kinematics::Jump const & jump2
 )
{

  RNA_PairingTemplateOP p( new RNA_PairingTemplate( jump1, jump2,  atom_name1, atom_name2) );

  //	std::cout << "JUMPY! " << reschar1 << " " << reschar2 << " " << edgechar1 << " " << edgechar2 << " "  << orientation << " " << jump << std::endl;

  // fill in a new base-pairing template
  //	Pairing_template_RNA t( Epos1, Epos2);

  //	Base_pair_type base_pair_type( reschar1, reschar2, edgechar1, edgechar2, orientation);
  //	pairing_template_map_RNA[ base_pair_type ].push_back( t );

  // Save to template lists with wildcards ('X') too.
  for ( Size i = 0; i <= 1; i++ ) {

    char const edgechar1_temp = i ? edgechar1:'X';

    for ( Size j = 0; j <= 1; j++ ) {

      char const edgechar2_temp = j ? edgechar2:'X';

      for ( Size k = 0; k <= 1; k++ ) {

				char const orientation_temp = k ? orientation:'X';

				BasePairType base_pair_type( reschar1, reschar2, edgechar1_temp, edgechar2_temp, orientation_temp);
				rna_pairing_template_map_[ base_pair_type ].push_back( p );

      }
    }
  }
}


} //farna
} //protocols
