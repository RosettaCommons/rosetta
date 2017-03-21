// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/fldsgn/BluePrint.cc
/// @brief the basic idea of BluePrint is on the remodel Possu wrote in rosetta++.
/// @author Nobuyasu Koga

// Unit header
#include <protocols/jd2/parser/BluePrint.hh>

// Package header

// Project header
#include <core/types.hh>
#include <core/sequence/ABEGOManager.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <boost/lexical_cast.hpp>

#include <core/kinematics/MoveMap.hh>
#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.parser.BluePrint" );

using namespace core;

namespace protocols {
namespace jd2 {
namespace parser {

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
BluePrint::BluePrint() :
	total_residue_( 0 ),
	total_residue_wolig_( 0 ),
	sequence_( "" ),
	secstruct_( "" ),
	strand_pairings_( "" ),
	helix_pairings_( "" ),
	hss_triplets_( "" )
{}


/// @brief value constructor
BluePrint::BluePrint( std::string const & filename ):
	total_residue_( 0 ),
	total_residue_wolig_( 0 ),
	sequence_( "" ),
	secstruct_( "" ),
	strand_pairings_( "" ),
	helix_pairings_( "" ),
	hss_triplets_( "" )
{
	if ( ! read_blueprint( filename ) ) {
		TR.Error << "Error in reading blueprint file " << filename << std::endl;
		exit(0);
	}
}


/// @brief destructor
BluePrint::~BluePrint()
{}


/// @brief copy constructor
BluePrint::BluePrint( BluePrint const & src ):
	ReferenceCount(),
	total_residue_( src.total_residue_ ),
	total_residue_wolig_( src.total_residue_wolig_ ),
	sequence_( src.sequence_ ),
	secstruct_( src.secstruct_ ),
	resnum_( src.resnum_ ),
	resname_( src.resname_ ),
	sstype_( src.sstype_ ),
	abego_( src.abego_ ),
	buildtype_( src.buildtype_ ),
	resnum_map_( src.resnum_map_ ),
	strand_pairings_( src.strand_pairings_ ),
	helix_pairings_( src.helix_pairings_ ),
	hss_triplets_( src.hss_triplets_ )
{}


/// @brief total residue number defined in blueprint file
BluePrint::Size
BluePrint::total_residue() const
{
	return total_residue_;
}


/// @brief total residue number defined in blueprint file
BluePrint::Size
BluePrint::total_residue_wolig() const
{
	return total_residue_wolig_;
}


/// @brief sequence defined in blueprint file
BluePrint::String
BluePrint::sequence() const
{
	return sequence_;
}


/// @brief amino acid type at a position in blueprint file
char
BluePrint::sequence( core::Size seqpos ) const
{
	debug_assert( seqpos >= 1 );
	return resname_[ seqpos ];
}


/// @brief secondary structures defined in blueprint file
BluePrint::String
BluePrint::secstruct() const
{
	return secstruct_;
}

/// @brief abego defined in blueprint file
utility::vector1< BluePrint::String >
BluePrint::abego() const
{
	return abego_;
}

/// @brief abego defined in blueprint file
BluePrint::String
BluePrint::abego( core::Size seqpos ) const
{
	runtime_assert( seqpos >= 1 && seqpos <= abego_.size() );
	return abego_[ seqpos ];
}

/// @brief secondary structure at a position in blueprint file
char
BluePrint::secstruct( core::Size seqpos ) const
{
	debug_assert( seqpos >= 1 && seqpos <= buildtype_.size() );
	return sstype_[ seqpos ];
}


/// @brief residue number of each position in blueprint file
BluePrint::Size
BluePrint::resnum( core::Size seqpos ) const
{
	debug_assert( seqpos >= 1 && seqpos <= buildtype_.size() );
	return resnum_[ seqpos ];
}


/// @brief translate residue number of pose to that of blueprint file
BluePrint::Size
BluePrint::resnum_map( Size resnum_pose ) const
{
	std::map<Size,Size>::const_iterator itr;
	itr = resnum_map_.find( resnum_pose );
	if ( itr != resnum_map_.end() ) {
		Size resnum_blueprint = itr->second;
		return resnum_blueprint;
	} else {
		return 0;
	}
}


/// @brief build type at each position in blueprint
char
BluePrint::buildtype( core::Size seqpos ) const
{
	debug_assert( seqpos >= 1 && seqpos <= buildtype_.size() );
	return buildtype_[ seqpos ];
}

/// @brief build type at each position in blueprint
BluePrint::String
BluePrint::extra( core::Size seqpos ) const
{
	if ( seqpos >= 1 && seqpos <= extra_.size() ) {
		return extra_[ seqpos ];
	} else {
		return "";
	}
}

/// @brief build type at each position in blueprint
BluePrint::String
BluePrint::insertion( core::Size i ) const
{
	debug_assert( i <= insertion_.size() && i > 0 );
	return insertion_[ i ];
}


/// @brief strand pairings defined at the line of SSPAIR in blueprint
BluePrint::String
BluePrint::strand_pairings() const
{
	return strand_pairings_;
}


/// @brief helix pairings defined at the line of HHPAIR in blueprint
BluePrint::String
BluePrint::helix_pairings() const
{
	return helix_pairings_;
}


/// @brief helix-strand-strand set at the line of HSSTRIPLET in blueprint
BluePrint::String
BluePrint::hss_triplets() const
{
	return hss_triplets_;
}


/// @brief secondary structure information
///BluePrint::SS_Info2_OP
///BluePrint::ssinfo() const
///{
/// return ss_info_;
///}


/// @brief read blueprint file
bool
BluePrint::read_blueprint( std::string const & filename )
{
	using utility::string_split;

	utility::io::izstream data( filename );
	if ( !data ) {
		TR.Error << "can not open blueprint file " << filename << std::endl;
		return false;
	}
	return read_blueprint_stream( data, filename );
}

/// @brief read blueprint file from stream
bool
BluePrint::read_blueprint_stream( std::istream & data, std::string const & filename = "" )
{
	// StranPairing info read from the line of SSPAIR
	//StrandPairings spairs;
	String spairs;

	String line;
	Size linecount( 0 ), count( 0 );
	while ( getline( data, line ) ) {
		linecount++;
		utility::vector1< String > tokens ( utility::split( line ) );

		// skip reading line that is commented out
		if ( tokens[1][0] == '#' ) continue;

		if ( tokens[1] == "FOLDINFO" || tokens[1] == "SSPAIR" || tokens[1] == "HHPAIR" || tokens[1] == "HSSTRIPLET" ) {
			// read the line of SSPAIR( FOLDINFO ), HHPAIR

			if ( tokens.size() > 2 ) {
				TR.Error << "error parsing at " << linecount << " in " << filename << std::endl;
				return false;
			}

			if ( tokens[1] == "SSPAIR" || tokens[1] == "FOLDINFO" ) {
				// set string of strand pairings
				strand_pairings_ = tokens[2];
			} else if ( tokens[1] == "HHPAIR" ) {
				// set string of helix pairings
				helix_pairings_ = tokens[2];
			} else if ( tokens[1] == "HSSTRIPLET" ) {
				// set string of hss triplets
				hss_triplets_ = tokens[2];
			}

		} else if ( tokens[1] == "INSERT" ) {
			// read the line of INSERT for reading insertion pdb file
			insertion_.push_back( tokens[2] );

		} else {

			// read the line where residue number, amino acid type, secondary structure and build type are written
			// don't need to specify build type in blueprint
			//flo aug 2012:
			//it would be nice to be able to put comments behind bluerint lines...
			for ( core::Size i(1); i <= tokens.size(); ++i ) {
				if ( tokens[i][0] == '#' ) {
					tokens.resize( i - 1 );
					break;
				}
			}
			TR.Debug << "line=" << line << " tokens=" << tokens << std::endl;
			runtime_assert( tokens.size() == 3 || tokens.size() == 4 || tokens.size() == 5 );

			count ++;
			Size ii = boost::lexical_cast<Size>( tokens[1] );
			char aa ( tokens[2][0] );
			char sec( tokens[3][0] );
			String abego("");
			if ( tokens[3].length() > 1 ) {
				core::sequence::ABEGOManager am;
				// check characters of abego is appropriate or not
				for ( Size k=2; k<=tokens[3].length(); k++ ) {
					am.symbol2index( tokens[3][k-1] );
				}
				abego = tokens[3].substr( 1, tokens[3].length() - 1 );
			} else {
				abego = "X";
			}

			// check the 3rd column where secondary structure and abego info are written
			if ( sec != 'E' && sec != 'H' && sec != 'L' && sec != 'D' ) {
				TR.Error << "unrecognized secstruct char : " << sec << " at lines " << linecount << " in " << filename << std::endl;
				return false;
			}

			resnum_.push_back( ii );
			resname_.push_back( aa );
			sstype_.push_back( sec );
			abego_.push_back( abego );
			if ( ii != 0 ) {
				resnum_map_.insert( std::map<Size, Size>::value_type( ii, count ) );
			}

			if ( tokens.size() >= 4 ) {

				char build( tokens[4][0] );
				if ( build != '.' && build != 'R' && build != 'I' && build != 'X' &&
						build != 'F' && build != 'P' && build != 'C' && build != 'W' ) {
					TR.Error << "unrecognized build char : " << build << " at lines " << linecount << " in " << filename << std::endl;
					return false;
				}
				if ( ii == 0 && build == '.' ) {
					TR.Error << "Cannot have simultaneously '0' at 1st column and '.' in 4th column in blueprint " << std::endl;
					return false;
				}
				buildtype_.push_back( build );

				if ( tokens.size() >= 5 ) {
					String extra( tokens[5] );
					extra_.push_back( extra );
				} else {
					extra_.push_back( "" );
				}

			} else {
				buildtype_.push_back( '.' );
			} // tokens.size() >=4

		}// tokens.size() == 2 || tokens.size() == 3
	} // while ( getline )

	debug_assert( resname_.size() == sstype_.size() );
	total_residue_ = sstype_.size();
	debug_assert( total_residue_ > 0 );

	total_residue_wolig_ = 0;
	for ( Size ii=1; ii<=total_residue(); ii++ ) {
		if ( extra( ii ) == "LIGAND" ) continue;
		total_residue_wolig_ ++;
	}

	for ( utility::vector1< char >::const_iterator iter = sstype_.begin(); iter != sstype_.end() ; ++iter ) {
		secstruct_ += *iter;
	}
	for ( utility::vector1< char >::const_iterator iter = resname_.begin(); iter != resname_.end() ; ++iter ) {
		sequence_ += *iter;
	}

	TR << secstruct_ << std::endl;
	TR << sequence_ << std::endl;

	//ss_info_->initialize( secstruct_ );
	//if( spairs.size() > 0 ){
	//strand_pairings_ = set_strand_pairings( ss_info_, spairs );
	//}

	return true;

} // read_blueprint


/// @brief
void
BluePrint::set_movemap( MoveMapOP & movemap )
{
	for ( Size i=1; i<=total_residue_ ; ++i ) {
		if ( buildtype( i ) == '.' ) {
			movemap->set_bb( i, true );
			movemap->set_chi( i, true );
		} else if ( buildtype( i ) == 'P' ) {
			movemap->set_bb( i, false );
			movemap->set_chi( i, true );
			TR.Debug << "At residue " << i << ", Backbone is freezed." << std::endl;
		} else if ( buildtype( i ) == 'W' ) {
			movemap->set_bb( i, true );
			movemap->set_chi( i, false );
			TR.Debug << "At residue " << i << ", Rotamer is freezed." << std::endl;
		} else if ( buildtype( i ) == 'F' ) {
			movemap->set_bb( i, false );
			movemap->set_chi( i, false );
			TR.Debug << "At residue " << i << ", both backbone and rotamer are freezed. " << std::endl;
		} else {
			movemap->set_bb( i, true );
			movemap->set_chi( i, true );
		}
	}
}

/// @brief set secondary structure into pose
void
BluePrint::insert_ss_into_pose( Pose & pose )
{
	for ( Size i=1; i<= pose.size(); i++ ) {
		Size resi = resnum_map( i );
		if ( resi == 0 ) continue;
		char ss = secstruct( resi );
		pose.set_secstruct( i, boost::lexical_cast<char>( ss ) );
	}
}

} // ns parser
} // ns jd2
} // ns protocols
