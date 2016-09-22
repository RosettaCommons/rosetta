// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#if (defined _WIN32) && (!defined WIN_PYROSETTA)
#define ZLIB_WINAPI  // REQUIRED FOR WINDOWS
#endif

#include <protocols/frags/TorsionFragment.hh>
#include <protocols/frags/VallData.hh>

// #include <devel/blab/using.hh>
// #include <devel/blab/tools.hh>
// #include <devel/blab/typedef.hh>

// Rosetta Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/types.hh>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/vector1.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <numeric/conversions.hh>

#include <basic/Tracer.hh>

//#include <devel/dna/DoubleFragment.hh>
// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>


// C++ Headers
// #include <cmath>
// #include <cstdlib>
// #include <iostream>
// #include <fstream>
// #include <sstream>

namespace protocols {
namespace frags {

using namespace core;
using namespace core::conformation;
using std::endl;
using std::string;
using utility::vector1;
static THREAD_LOCAL basic::Tracer TR( "protocols.frags.TorsionFragment" );

typedef utility::vector1< core::Size > Sizes;


/// these will go into a helper file soon
template < class T >
bool
has_element( utility::vector1< T > const & v, T const & t )
{
	return ( std::find( v.begin(), v.end(), t ) != v.end() );
}
/// never can remember the order...
template < class T >
bool
has_element( T const & t, utility::vector1< T > const & v )
{
	return ( std::find( v.begin(), v.end(), t ) != v.end() );
}


TorsionFragment::~TorsionFragment() = default;

TorsionFragment::TorsionFragment( TorsionFragment const & src ):
	ReferenceCount(),
	torsions_( src.torsions_ ),
	secstruct_( src.secstruct_ ),
	sequence_ ( src.sequence_ )
{}

TorsionFragmentOP
TorsionFragment::clone() const
{
	return TorsionFragmentOP( new TorsionFragment( *this ) );
}

///\brief insert this piece of fragment to a pose at position "begin"
///
///call pose.set_torsion which maps TorsionID to DOF_ID, so it is safe
///to use even if the folding direction is not standard as N to C.
void
TorsionFragment::insert( pose::Pose & pose, Size const begin ) const
{
	for ( Size i=1; i<= size(); ++i ) {
		utility::vector1< Real > const & bb( torsions_[i] );
		Size const seqpos( begin + i - 1 );
		for ( Size j=1; j<= bb.size(); ++j ) {
			pose.set_torsion( id::TorsionID( seqpos, id::BB, j ), bb[j] );
		}
		{ // pbhacking
			conformation::Residue const & rsd( pose.residue( seqpos ) );
			if ( rsd.is_protein() && rsd.is_upper_terminus() ) {
				// special case for psi
				Real const psi( bb[2] );
				id::AtomID const id1( rsd.atom_index("N"), seqpos ), id2( rsd.atom_index("CA"), seqpos ),
					id3( rsd.atom_index("C"), seqpos ), id4( rsd.atom_index("O"), seqpos );
				pose.conformation().set_torsion_angle( id1, id2, id3, id4, numeric::conversions::radians( psi ) );
				//     Real const actual_psi( numeric::dihedral_degrees( pose.xyz( id1 ), pose.xyz( id2 ), pose.xyz( id3 ),
				//                              pose.xyz( id4 ) ) );
				//     if ( std::abs( util::subtract_degree_angles( psi, actual_psi ) ) > 1.0 ) {
				//      std::cerr << "setting psi at terminus failed: desired= " << psi << " actual= " << actual_psi << std::endl;
				//      std::cout << "setting psi at terminus failed: desired= " << psi << " actual= " << actual_psi << std::endl;
				//     }
			}
		}
		pose.set_secstruct( seqpos, secstruct_[i] );
	}
}

void
SingleResidueTorsionFragmentLibrary::copy_fragments( SingleResidueTorsionFragmentLibrary const & src )
{
	for ( Size i=1; i<= src.size(); ++i ) {
		append_fragment( src[i].clone() );
	}
}


//////////////////////////////////////////////////////////////////////////////////
///\brief initialize fragment data from a classic Rosetta fragment library
///
///frag_size is the size of fragment in the library and nbb is the number of backbone
///torision angles stored in each fragment, for protein, this will be 3. Return false if
///there is a reading error and all the data which have been stored will be erased.
///
bool
TorsionFragmentLibrary::read_file(
	std::string const & filename,
	Size const frag_size,
	Size const nbb
)
{
	utility::io::izstream data ( filename );
	if ( !data ) {
		std::cerr << "Cannot open " + ObjexxFCL::string_of(frag_size) + "-mer fragment library file: " + filename + "\n" ;
		return false;
	}
	std::string line;
	while ( getline( data, line ) ) {
		std::istringstream line_stream( line );
		std::string tag1, tag2;
		Size position, neighbors;
		line_stream >> tag1 >> position >> tag2 >> neighbors;
		if ( line_stream.fail() || tag1 != "position:" || tag2 != "neighbors:" ) {
			std::cerr << " format errors in fragment library file: " << line << std::endl;
			resize(0);
			return false;
		}
		resize( position );
		getline(data, line); // skip blank line
		SingleResidueTorsionFragmentLibrary & current_position_library( (*this)[position] );
		for ( Size i=1; i<=neighbors; ++i ) {
			TorsionFragmentOP fragment( new TorsionFragment( frag_size, nbb ) );
			// read lines within each fragment
			std::string last_pdb; char last_chain('0'); Size last_seqpos(0);
			for ( Size j=1; j<=frag_size; ++j ) {
				getline(data, line);
				std::istringstream line_stream( line );
				std::string pdb; char chain, aa, secstruct; Size seqpos;
				line_stream >> pdb >> chain >> seqpos >> aa >> secstruct;
				if ( j == 1 ) {
					last_pdb = pdb;
					last_chain = chain;
					last_seqpos = seqpos;
				} else {
					if ( last_pdb != pdb  || last_chain != chain ) { // || last_seqpos != (seqpos-j+1) ) {
						std::cerr << "fragment reading error -- pdb, chain and seqpos mismatch\n"
							<< "position: " << position  << "; neighbor: " << i  << "; line: " << j << std::endl;
						resize(0);
						return false;
					}
					if ( last_seqpos != (seqpos-j+1) ) {
						TR.Error << "fragment reading POSSIBLE error -- seqpos mismatch: " << pdb << ' ' << chain << ' ' <<
							seqpos << " != " << last_seqpos+j-1 << " position: " << position  << "; neighbor: " << i <<
							"; line: " << j << std::endl;
					}
				}
				fragment->set_secstruct(j,secstruct);
				fragment->set_sequence (j,aa); // NEW
				for ( Size k=1; k<=nbb; ++k ) {
					Real torsion;
					line_stream >> torsion;
					fragment->set_torsion(j,k,torsion);
				} // bb torsion
			} // fragment
			// append this fragment to SingleResidueFragmentLibrary
			current_position_library.append_fragment( fragment );
			getline(data,line); // skip blank line
		} // SingleResidueFragmentLibrary
	} // TorsionFragmentLibrary
	TR.Info << "read succesfully " << frag_size << "-mer fragment library file: " << filename << std::endl;
	return true;
}


TorsionFragmentLibrary::~TorsionFragmentLibrary() = default;

/// after calling, the size will be the oldsize + current2desired_offset
///
void
TorsionFragmentLibrary::shift( int const current2desired_offset )
{
	utility::vector1< SingleResidueTorsionFragmentLibraryOP > new_fragments;

	Size const oldsize( fragments_.size() ), newsize( oldsize + current2desired_offset );
	new_fragments.resize( newsize );

	for ( int i=1; i<= int( newsize ); ++i ) {
		int j( i - current2desired_offset );
		if ( j >=1 && j <= int( oldsize ) ) {
			new_fragments[ i ] = fragments_[ j ];
		} else {
			new_fragments[ i ] = protocols::frags::SingleResidueTorsionFragmentLibraryOP( new SingleResidueTorsionFragmentLibrary() );
		}
	}
	fragments_.swap( new_fragments );
}

void
TorsionFragmentLibrary::delete_residue( Size const seqpos )
{
	Size const old_size( fragments_.size() );

	TR.Trace << "TorsionFragmentLibrary::delete_residue: " << seqpos << " old_size: " << old_size << std::endl;

	for ( Size i=1; i< seqpos; ++i ) {
		/// check to see if the fragments at this position overlap the deleted position
		/// if so, clear them
		if ( fragments_[i] && ( fragments_[i]->size() > 0 ) ) {
			Size const fragsize( (*fragments_[i])[1].size() );
			if ( i + fragsize - 1 >= seqpos ) {
				TR.Trace << "TorsionFragmentLibrary::delete_residue: erasing " << fragments_[i]->size() <<
					" fragments of size " << fragsize << " at position: " << i << std::endl;
				fragments_[i]->clear();
			}
		}
	}

	for ( Size i=seqpos; i< old_size; ++i ) {
		fragments_[i] = fragments_[i+1]; // these are just (owning) pointers... some may be null.
	}

	fragments_.resize( old_size-1 );
}


void
TorsionFragmentLibrary::copy_fragments( TorsionFragmentLibrary const & src )
{
	if ( size() < src.size() ) resize( src.size() );

	for ( Size i=1; i<= src.size(); ++i ) {
		fragments_[i]->copy_fragments( src[i] );
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////
///\brief extract a fragment library with smaller fragment size from the one with larger lize
///
/// for example, 1-mer library can ben extracted from 3-mer library. This function is set up
/// in a general way that both sizes of smaller fragment and larger fragment can be flexibly
/// specified by my_size and src_size. Also, the source library does not have to contain
/// fragment data for every residue position. For example, for loop modeling, we only need
/// fragment for loop segment and after this function is called, smaller library is created
/// only for regions in which there is data in the larger library. Lastly, smaller fragment
/// is aligned in the center of larger fragment to have data extracted.
bool
TorsionFragmentLibrary::derive_from_src_lib(
	Size my_size,
	Size src_size,
	TorsionFragmentLibrary const & src_lib
)
{
	// can only go from larger size to smaller size
	runtime_assert( my_size < src_size );
	// get continuous segment in src_lib_op
	std::map< Size, Size > seg_map;
	Size prev_nbrs(0);
	Size seg_start(0);
	for ( Size i = 1; i <= src_lib.size(); ++i ) {
		Size current_nbrs( src_lib[i].size() );
		if ( !prev_nbrs && current_nbrs ) {
			runtime_assert( !seg_start );
			seg_start = i;
		} else if ( prev_nbrs && !current_nbrs ) {
			runtime_assert( seg_start && seg_start <i );
			seg_map.insert( std::make_pair( seg_start, i-1 ) );
			seg_start = 0;
		} else if ( i == src_lib.size() && prev_nbrs && current_nbrs ) {
			runtime_assert( seg_start && seg_start < i );
			seg_map.insert( std::make_pair( seg_start, i ) );
			seg_start = 0;
		}
		prev_nbrs = current_nbrs;
	}
	runtime_assert( !seg_start );

	if ( seg_map.empty() ) { //size() ) {
		std::cerr << "Warning: source fragment library does not have any data!" << std::endl;
		return false;
	}

	// resize myself lib to correct size
	Size const size_diff(src_size - my_size);
	resize( src_lib.size() + size_diff );
	// find offset by which shorter frag is aligned to longer frag
	Size const offset( size_diff%2 == 0 ? size_diff/2 : (size_diff-1)/2 );

	// loop through each segment and derive fragments
	for ( std::map<Size, Size>::const_iterator it = seg_map.begin(),
			it_end = seg_map.end(); it != it_end; ++it ) {
		Size const seg_begin(it->first);
		Size const seg_end(it->second);
		Size lib_index, copy_start;
		for ( Size my_index = seg_begin, my_end = seg_end+size_diff; my_index <= my_end; my_index++ ) {
			// get hold of current residue frag lib
			SingleResidueTorsionFragmentLibraryOP my_residue_lib( fragments_[my_index] );
			// figure out which subset of fragment data we should copy
			if ( my_index < seg_begin + offset ) {
				lib_index = seg_begin;
				copy_start = 1 + my_index - seg_begin;
			} else if ( my_index > seg_end + offset ) {
				lib_index = seg_end;
				copy_start = 1 + my_index - seg_end;
			} else {
				lib_index = my_index - offset;
				copy_start = 1 + offset;
			}
			// do actual copying from src residue_fragment_library
			SingleResidueTorsionFragmentLibrary const & src_residue_lib( src_lib[lib_index] );
			for ( Size i = 1; i <= src_residue_lib.size(); ++i ) {
				TorsionFragment const & src_fragment( src_residue_lib[i] );
				Size const nbb( src_fragment.nbb() );
				TorsionFragmentOP my_fragment( new TorsionFragment( my_size, nbb ) );
				for ( Size j = 1, jj = copy_start; j <= my_size; ++j, ++jj ) {
					for ( Size k = 1; k <= nbb; ++k ) {
						my_fragment->set_torsion( j, k, src_fragment.get_torsion(jj, k) );
					} // each torsion
					my_fragment->set_secstruct( j, src_fragment.get_secstruct(jj) );
					my_fragment->set_sequence ( j, src_fragment.get_sequence (jj) );
				} // each frag_pos
				my_residue_lib->append_fragment( my_fragment );
			} // each neighbor
		} // each residue position
	} //seg_map

	// successful
	return true;
}
bool
TorsionFragmentLibrary::derive_from_src_lib(
	Size my_size,
	Size src_size,
	TorsionFragmentLibraryCOP src_lib_op
)
{
	return this->derive_from_src_lib( my_size, src_size, *src_lib_op );
}

void
FragLib::delete_residue( Size const seqpos )
{
	for ( auto & it : frag_map_ ) {
		it.second->delete_residue( seqpos );
	}
}

void
FragLib::copy_fragments( FragLib const & src )
{
	Sizes const fragsizes( src.frag_sizes() );
	for ( core::Size fragsize : fragsizes ) {
		library( fragsize ).copy_fragments( src.library( fragsize ) );
	}
}

void
FragLib::shift( int const current2desired_offset )
{
	for ( auto & it : frag_map_ ) {
		it.second->shift( current2desired_offset );
	}
}

TorsionFragmentLibrary const &
FragLib::library( Size const size ) const
{
	if ( frag_map_.count( size ) == 0 ) {
		utility_exit_with_message( "No frag lib for that size! "+ObjexxFCL::string_of( size ) );
	}
	return *( frag_map_.find( size )->second );
}

TorsionFragmentLibrary &
FragLib::library( Size const size )
{
	if ( frag_map_.count( size ) == 0 ) {
		TR.Info << "Creating new fragment library for frag_size " << size << endl;
		frag_map_[ size ] = protocols::frags::TorsionFragmentLibraryOP( new TorsionFragmentLibrary() );
	}
	return *( frag_map_.find( size )->second );
}

utility::vector1< Size >
FragLib::frag_sizes() const
{
	utility::vector1< Size > sizes;
	for ( auto const & it : frag_map_ ) {
		sizes.push_back( it.first );
	}
	return sizes;
}


bool
ss_length_check(
	core::Size const min_len_helix,
	core::Size const min_len_strand,
	core::pose::Pose const & pose
)
{
	std::string const ss( pose.secstruct() );

	bool passed( true );

	for ( Size i=0; i< ss.size(); ++i ) {

		if ( ss[i] == 'H' && ( i+1 == ss.size() || ss[i+1] != 'H' ) ) {
			// at the last residue of helix
			// figure out where this helix starts
			Size j(i);
			while ( j>0 && ss[j-1] == 'H' ) --j;
			runtime_assert( ss[j] == 'H' && ( j==0 || ss[j-1] != 'H' ) );
			Size const len( i-j+1 );
			if ( len < min_len_helix ) passed = false;
		}

		if ( ss[i] == 'E' && ( i+1 == ss.size() || ss[i+1] != 'E' ) ) {
			// at the last residue of strand
			// figure out where this strand starts
			Size j(i);
			while ( j>0 && ss[j-1] == 'E' ) --j;
			runtime_assert( ss[j] == 'E' && ( j==0 || ss[j-1] != 'E' ) );
			Size const len( i-j+1 );
			if ( len < min_len_strand ) passed = false;
		}
	}
	return passed;
}

bool
TorsionFragmentMover::pose_passes_ss_length_check( pose::Pose const & pose ) const
{
	return ss_length_check( min_len_helix_, min_len_strand_, pose );
}


void
TorsionFragmentMover::apply( pose::Pose & pose )
{
	last_inserted_frag_begin_ = 0; // initialize

	TorsionFragmentLibrary const & lib( fraglib_->library( frag_size_ ) );

	// figure out where we can insert
	utility::vector1< Size > windows;

	for ( Size i=1; i<= pose.size() - frag_size_ + 1; ++i ) {
		bool allowed( true );
		for ( Size j=0; j< frag_size_; ++j ) {
			// not sure if this movemap call should be allowed
			if ( !mm_->get_bb( i + j ) || !pose.residue( i+j ).is_protein() ) {
				allowed = false;
				break;
			}
		}
		if ( allowed ) windows.push_back( i );
	}

	if ( windows.empty() ) {
		TR.Warning << "no insertable positions!" << endl;
		return;
	}

	bool check_ss_lengths_this_time( check_ss_lengths_ );
	if ( check_ss_lengths_this_time && !pose_passes_ss_length_check( pose ) ) {
		check_ss_lengths_this_time = false;
		TR.Warning << "TorsionFragmentLibrary::apply: check_ss_lengths_ = true but incoming pose fails ss_check!" <<
			std::endl;
	}

	Size ntries(0);
	Pose save_pose;
	if ( check_ss_lengths_this_time ) save_pose = pose;
	while ( ntries < 1000 ) {
		++ntries;
		Size const window( numeric::random::rg().random_element( windows ) );
		Size const nfrags = lib[ window ].size();
		if ( nfrags <= 0 ) continue;
		Size const nn( static_cast< int > ( numeric::random::uniform() * nfrags ) + 1 );

		// insert fragment
		lib[ window ][ nn ].insert( pose, window );
		last_inserted_frag_begin_ = window;

		if ( check_ss_lengths_this_time && !pose_passes_ss_length_check( pose ) ) {
			TR.Trace << "TorsionFragmentMover::apply: sscheck failed! ntries=" << ntries << endl;
			pose = save_pose;
			last_inserted_frag_begin_ = 0;
			continue;
		}

		break;
	}
	if ( ntries >= 1000 ) TR.Warning << "all insertable frag windows have nfrags==0!" << endl;

} // apply

utility::vector1< Size > const empty_size_list;
//std::string const empty_string;


void
add_vall_fragments(
	utility::vector1< core::Size > const & frag_sizes,
	Size const nfrags,
	pose::Pose const & pose,
	core::kinematics::MoveMap const & mm,
	std::string const & secstruct,
	Real const seq_weight,
	Real const ss_weight,
	FragLib & frag_lib,
	utility::vector1< Size > const & homs_to_exclude, // = empty_size_list
	Real const bb_weight, // = 0
	std::string const & bigbins, // = empty string
	std::string const & inputseq // = empty string
)
{

	Size const nres( pose.size() );
	string seq( inputseq );
	if ( seq.empty() ) seq = pose.sequence();
	//string const seq( pose.sequence() );

	for ( core::Size frag_size : frag_sizes ) {
		TorsionFragmentLibrary & lib( frag_lib.library( frag_size ) );
		lib.resize( nres - frag_size + 1 );
		for ( Size i = 1; i<= nres - frag_size+1; ++i ) {
			bool allowed( true );
			for ( Size j=i; j<i+frag_size; ++j ) {
				Residue const & jrsd( pose.residue(j) );
				if ( ( j > i             && jrsd.is_lower_terminus() ) ||
						( j < i+frag_size-1 && jrsd.is_upper_terminus() ) ||
						!pose.residue(j).is_protein() ||
						!mm.get_bb( j ) ) { // this call should not be allowed!
					allowed = false;
					break;
				}
			}
			if ( !allowed ) continue;
			std::string const frag_seq( seq.substr(i-1,frag_size) );
			std::string const frag_ss( secstruct.substr(i-1,frag_size ) );
			std::string const frag_bb( bigbins.empty() ? bigbins : bigbins.substr( i-1,frag_size) );
			bool const exclude_gly( false ), exclude_pro( false ), exclude_cis_peptides( true ); // except at pre-pro
			get_frags( nfrags, frag_seq, frag_ss, seq_weight, ss_weight, exclude_gly, exclude_pro, exclude_cis_peptides,
				homs_to_exclude, lib[i], bb_weight, frag_bb );
		} // for each residue window
	}
}

void
add_vall_fragments(
	utility::vector1< core::Size > const & frag_sizes,
	Size const nfrags,
	pose::Pose const & pose,
	core::kinematics::MoveMap const & mm,
	utility::vector1< std::map< char, core::Real > > const & target_ss,
	Real const seq_weight,
	Real const ss_weight,
	FragLib & frag_lib,
	utility::vector1< Size > const & homs_to_exclude, // = empty_size_list
	Real const bb_weight, // = 0
	std::string const & bigbins, // = empty string
	std::string const & inputseq // = empty string
)
{

	Size const nres( pose.size() );
	string seq( inputseq );
	if ( seq.empty() ) seq = pose.sequence();
	//string const seq( pose.sequence() );

	for ( core::Size frag_size : frag_sizes ) {
		TorsionFragmentLibrary & lib( frag_lib.library( frag_size ) );
		lib.resize( nres - frag_size + 1 );
		for ( Size i = 1; i<= nres - frag_size+1; ++i ) {
			bool allowed( true );
			utility::vector1< std::map< char, core::Real > > frag_ss;
			for ( Size j=i; j<i+frag_size; ++j ) {
				Residue const & jrsd( pose.residue(j) );
				if ( ( j > i             && jrsd.is_lower_terminus() ) ||
						( j < i+frag_size-1 && jrsd.is_upper_terminus() ) ||
						!pose.residue(j).is_protein() ||
						j > target_ss.size() ||
						!mm.get_bb( j ) ) { // this call should not be allowed!
					allowed = false;
					break;
				}
				frag_ss.push_back( target_ss[j] );
			}
			if ( !allowed ) continue;
			std::string const frag_seq( seq.substr(i-1,frag_size) );
			//std::string const frag_ss( secstruct.substr(i-1,frag_size ) );
			std::string const frag_bb( bigbins.empty() ? bigbins : bigbins.substr( i-1,frag_size) );
			bool const exclude_gly( false ), exclude_pro( false ), exclude_cis_peptides( true ); // except at pre-pro
			get_frags( nfrags, frag_seq, frag_ss, seq_weight, ss_weight, exclude_gly, exclude_pro, exclude_cis_peptides,
				homs_to_exclude, lib[i], bb_weight, frag_bb );
		} // for each residue window
	}
}

FragLibOP
setup_vall_fragments(
	utility::vector1< core::Size > const & frag_sizes,
	Size const nfrags,
	pose::Pose const & pose,
	core::kinematics::MoveMap const & mm,
	std::string const & secstruct,
	Real const seq_weight,
	Real const ss_weight,
	utility::vector1< Size > const & homs_to_exclude // = empty_size_list
)
{
	FragLibOP frag_lib( new FragLib() );
	add_vall_fragments( frag_sizes, nfrags, pose, mm, secstruct, seq_weight, ss_weight, *frag_lib, homs_to_exclude );
	return frag_lib;
}


void
add_vall_cheating_fragments(
	utility::vector1< core::Size > const & frag_sizes,
	Size const nfrags,
	pose::Pose const & pose,
	core::kinematics::MoveMap const & mm,
	std::string const & secstruct,
	Real const seq_weight,
	Real const ss_weight,
	Real const torsion_weight,
	Real const min_torsion_dev,
	Real const max_torsion_dev,
	FragLib & frag_lib,
	utility::vector1< Size > const & homs_to_exclude // = empty_size_list
)
{
	//FragLibOP frag_lib( new FragLib() );

	Size const nres( pose.size() );
	string const seq( pose.sequence() );

	for ( core::Size frag_size : frag_sizes ) {
		TorsionFragmentLibrary & lib( frag_lib.library( frag_size ) );
		lib.resize( nres - frag_size + 1 );
		for ( Size i = 1; i<= nres - frag_size+1; ++i ) {
			bool allowed( true );
			for ( Size j=i; j<i+frag_size; ++j ) {
				Residue const & jrsd( pose.residue(j) );
				if ( ( j > i             && jrsd.is_lower_terminus() ) ||
						( j < i+frag_size-1 && jrsd.is_upper_terminus() ) ||
						!pose.residue(j).is_protein() ||
						!mm.get_bb( j ) ) { // this call should not be allowed!
					allowed = false;
					break;
				}
			}
			if ( !allowed ) continue;
			std::string const frag_seq( seq.substr(i-1,frag_size) );
			std::string const frag_ss( secstruct.substr(i-1,frag_size ) );
			utility::vector1< Real > frag_phi, frag_psi, frag_omega;
			for ( Size j=i; j<i+frag_size; ++j ) {
				frag_phi.push_back( pose.phi(j) );
				frag_psi.push_back( pose.psi(j) );
				frag_omega.push_back( pose.omega(j) );
				if ( pose.residue(j).is_lower_terminus() ) {
					runtime_assert( j == i );
					frag_phi.back() = 0.0;
				}
				if ( pose.residue(j).is_upper_terminus() ) {
					runtime_assert( j == i+frag_size-1 );
					frag_psi.back() = 0.0;
					frag_omega.back() = 0.0;
				}
			}


			get_cheating_frags( nfrags, frag_seq, frag_ss, frag_phi, frag_psi, frag_omega,
				seq_weight, ss_weight, torsion_weight,
				min_torsion_dev, max_torsion_dev,
				homs_to_exclude, lib[i] );

		} // for each residue window
	}
}


FragLibOP
setup_vall_cheating_fragments(
	utility::vector1< core::Size > const & frag_sizes,
	Size const nfrags,
	pose::Pose const & pose,
	core::kinematics::MoveMap const & mm,
	std::string const & secstruct,
	Real const seq_weight,
	Real const ss_weight,
	Real const torsion_weight,
	Real const min_torsion_dev,
	Real const max_torsion_dev,
	utility::vector1< Size > const & homs_to_exclude // = empty_size_list
)
{
	FragLibOP frag_lib( new FragLib() );
	add_vall_cheating_fragments( frag_sizes, nfrags, pose, mm, secstruct, seq_weight, ss_weight, torsion_weight,
		min_torsion_dev, max_torsion_dev, *frag_lib, homs_to_exclude );
	return frag_lib;
}


/// @details  Fill in gaps in the FragLib with fragments chosen from the vall
/// @note  Does not change the size of any of the individual fragment libraries
void
fill_in_gaps(
	Size const nfrags,
	pose::Pose const & pose,
	std::string const & secstruct,
	Real const seq_weight,
	Real const ss_weight,
	FragLib & frag_lib,
	utility::vector1< Size > const & homs_to_exclude, // = empty_size_list
	bool const allow_uninitialized_secstruct // = false
)
{
	string const seq( pose.sequence() );

	if ( secstruct == std::string( secstruct.size(), 'L' ) && !allow_uninitialized_secstruct ) {
		utility_exit_with_message("fill_in_gaps:: uninitialized secstruct?");
	}

	Sizes const frag_sizes( frag_lib.frag_sizes() );
	for ( core::Size frag_size : frag_sizes ) {
		TorsionFragmentLibrary & lib( frag_lib.library( frag_size ) );
		for ( Size i = 1; i<= lib.size(); ++i ) {
			if ( lib[i].size() > 0 ) continue; // skip windows with at least one fragment
			bool allowed( true );
			for ( Size j=i; j<i+frag_size; ++j ) {
				if ( j > pose.size() || !pose.residue(j).is_protein() ) {
					allowed = false;
					break;
				}
			}
			if ( !allowed ) continue;
			std::string const frag_seq( seq.substr(i-1,frag_size) );
			std::string const frag_ss( secstruct.substr(i-1,frag_size ) );
			TR.Trace << "filling gap in frag_lib: frag_size= " << frag_size << " window= " << i << std::endl;
			get_frags( nfrags, frag_seq, frag_ss, seq_weight, ss_weight, false, false, true, homs_to_exclude, lib[i] );
		} // for each residue window
	}
}


/// @details  Show size info about library
void
TorsionFragmentLibrary::print( std::ostream & os ) const
{
	os << "TorsionFragmentLibrary:: size= " << size() << std::endl;
	for ( Size i=1; i<= size(); ++i ) {
		os << "window: " << i << " depth: " << fragments_[i]->size() << std::endl;
	}
}

SingleResidueTorsionFragmentLibrary::~SingleResidueTorsionFragmentLibrary() = default;


///////////////////////////////////////////////////////////////////////////////
/// @details if desired_insert_pos is specified, the fragmnet will be inserted
/// into that position. Otherwise, it will be randonly inserted into a region
/// between "begin" and "end". The fragment insetion is done by changing according
/// backbone torsions in the pose and "refolding" will be handled by the Pose before
/// next time XYZ coordindates are being accessed.
/// PBMOD: allow for cases where not all windows between begin and end have fragments by
/// adding a while loop within which insert_pos is chosen and the frag depth is checked.
///
///////////////////////////////////////////////////////////////////////////////
void
insert_fragment(
	int const begin,
	int const end,
	pose::Pose & pose,
	TorsionFragmentLibrary const & lib,
	int const desired_insert_pos // = 0
)
{
	bool const paranoid( false ); // should we require exact length match between pose and fraglib?

	// number of single residue window in the fragment library
	int const frag_nres = lib.size();
	// determine frag_size from lib, cant do THIS: "int const frag_size = int( lib[begin][1].size() );"
	int frag_size(0);
	for ( int i=begin; i<= end && i <= frag_nres; ++i ) {
		if ( lib[i].size() ) {
			int const this_frag_size( lib[i][1].size() );
			runtime_assert( frag_size == 0 || frag_size == this_frag_size );
			frag_size = this_frag_size;
		}
	}
	if ( !frag_size ) utility_exit_with_message("no insertable windows in region!");
	// number of total residue in pose
	int const pose_nres = pose.size();
	// ensure fragment library and pose matches each other
	if ( paranoid && pose_nres != (frag_nres + frag_size - 1) ) {
		std::cerr << pose_nres << "  " << frag_nres << "  " << frag_size << std::endl;
		utility_exit_with_message( "Fragment files and input PDB file don't match in length!" );
	}

	// length of insertable region
	int const region_size ( end - begin + 1);
	// in case it is smaller than the size of fragment
	int const actual_frag_size
		( region_size < frag_size  ? region_size : frag_size );
	int ntries( 0 ), insert_pos( 0 ), nfrags( 0 );
	while ( ntries < 1000 && nfrags <= 0 ) {
		++ntries;
		// choose a window in the insertable region
		int const pos ( static_cast<int>( numeric::random::uniform() * ( region_size - actual_frag_size + 1) ) );
		insert_pos = ( desired_insert_pos == 0 ? begin + pos : desired_insert_pos );
		runtime_assert ( 1<= insert_pos && insert_pos <= frag_nres );
		// choose a fragment from the library at that window position
		nfrags = lib[ insert_pos ].size();
	}
	if ( ntries >= 1000 ) {
		TR.Warning << "insert_fragment: no fragments found between " << begin << " and " << end << endl;
		return;
	}
	int const nn = static_cast< int > ( numeric::random::uniform() * nfrags) + 1;
	// insert fragment
	lib[ insert_pos ][ nn ].insert( pose, insert_pos );
}

std::ostream &
operator << ( std::ostream & out, TorsionFragment const & f )
{
	out << "TorsionFragment " << f.size() << ' '<< f.nbb();
	for ( Size i=1; i<= f.size(); ++i ) {
		out << ' ' << f.get_secstruct(i);
		out << ' ' << f.get_sequence (i); // adding SEQ
		for ( Size j=1; j<= f.nbb(); ++j ) {
			out << ObjexxFCL::format::F(9,3,f.get_torsion( i, j ) );
		}
	}

	return out;
}

std::istream &
operator >> ( std::istream & data, TorsionFragment & f )
{
	string tag;
	Size size, nbb;
	data >> tag >> size >> nbb;
	if ( tag != "TorsionFragment" ) {
		data.setstate( std::ios_base::failbit );
		return data;
	}
	f.set_size_and_nbb( size, nbb );
	for ( Size i=1; i<= size; ++i ) {
		char ss,aa;
		data >> ss >> aa;
		f.set_secstruct( i, ss );
		f.set_sequence ( i, aa );
		for ( Size j=1; j<= nbb; ++j ) {
			Real torsion;
			data >> torsion;
			f.set_torsion( i, j, torsion );
		}
	}
	return data;
}

std::ostream &
operator << ( std::ostream & out, SingleResidueTorsionFragmentLibrary const & lib )
{
	out << "SingleResidueTorsionFragmentLibrary " << lib.size() << '\n';
	for ( Size i=1; i<= lib.size(); ++i ) {
		out << lib[i] << '\n';
	}
	return out;
}

std::istream &
operator >> ( std::istream & data, SingleResidueTorsionFragmentLibrary & lib )
{
	lib.clear();

	string tag;
	Size size;
	data >> tag >> size;
	if ( tag != "SingleResidueTorsionFragmentLibrary" ) {
		data.setstate( std::ios_base::failbit );
		return data;
	}

	for ( Size i=1; i<= size; ++i ) {
		TorsionFragmentOP f( new TorsionFragment() );
		data >> (*f);
		lib.append_fragment( f );
	}
	return data;
}


std::ostream &
operator << ( std::ostream & out, TorsionFragmentLibrary const & lib )
{
	out << "TorsionFragmentLibrary " << lib.size() << '\n';
	for ( Size i=1; i<= lib.size(); ++i ) {
		out << lib[i];
	}
	return out;
}

std::istream &
operator >> ( std::istream & data, TorsionFragmentLibrary & lib )
{
	lib.clear();

	string tag;
	Size size;
	data >> tag >> size;
	if ( data.fail() || tag != "TorsionFragmentLibrary" ) {
		data.setstate( std::ios_base::failbit );
		return data;
	}

	lib.resize( size );

	for ( Size i=1; i<= size; ++i ) data >> lib[i];

	return data;
}


std::ostream &
operator << ( std::ostream & out, FragLib const & lib )
{
	Sizes const fragsizes( lib.frag_sizes() );
	out << "FragLib " << fragsizes.size() << '\n';
	for ( core::Size fragsize : fragsizes ) {
		out << "frag_size " << fragsize << ' ' << lib.library( fragsize );
	}
	return out;
}

std::istream &
operator >> ( std::istream & data, FragLib & lib )
{
	lib.clear();

	string tag;
	Size size;
	data >> tag >> size;
	if ( tag != "FragLib" ) {
		data.setstate( std::ios_base::failbit );
		return data;
	}

	for ( Size i=1; i<= size; ++i ) {
		Size fragsize;
		data >> tag >> fragsize;
		if ( tag != "frag_size" ) {
			data.setstate( std::ios_base::failbit );
			return data;
		}
		data >> lib.library( fragsize );
	}
	return data;
}

/// @details  Determine length of shortest and longest contiguous segments
void
get_min_and_max_contigs( Sizes const & seg_in, core::Size & min_seg, core::Size & max_seg )
{
	Sizes seg( seg_in );
	std::sort( seg.begin(), seg.end() );
	Size last_pos( 0 ), seg_begin( 0 );
	utility::vector1< Size > sizes;
	for ( Sizes::const_iterator pos= seg.begin(); pos != seg.end(); ++pos ) {
		if ( !seg_begin ) {
			seg_begin = *pos;
		} else if ( *pos != last_pos+1 ) {
			runtime_assert( last_pos && seg_begin );
			sizes.push_back( last_pos - seg_begin + 1 );
			seg_begin = *pos;
		}
		last_pos= *pos;
	}
	sizes.push_back( last_pos - seg_begin + 1 );
	std::sort( sizes.begin(), sizes.end() );
	min_seg = sizes.front();
	max_seg = sizes.back();
}


void
insert_random_fragments_in_flexible_protein_regions(
	Sizes const & flex_protein,
	FragLib const & frag_lib,
	pose::Pose & pose
)
{
	Size const nres( pose.size() );
	vector1< Size > frag_sizes( frag_lib.frag_sizes() );
	Size min_seg, max_seg;
	get_min_and_max_contigs( flex_protein, min_seg, max_seg );
	if ( min_seg < utility::min( frag_sizes ) ) utility_exit_with_message("no frag small enough to match smallest contiguous regn");
	/// insert fragments from smallest to largest...
	std::sort( frag_sizes.begin(), frag_sizes.end() );
	Sizes all_positions; for ( Size i=1; i<= nres; ++i ) all_positions.push_back(i);
	for ( vector1< Size >::const_iterator size= frag_sizes.begin(); size!= frag_sizes.end(); ++size ) {
		Size const frag_size( *size );
		protocols::frags::TorsionFragmentLibrary const & lib( frag_lib.library( frag_size ) );
		numeric::random::random_permutation( all_positions, numeric::random::rg() );
		for ( Sizes::const_iterator pos = all_positions.begin(); pos != all_positions.end(); ++pos ) {
			Size const i( *pos );
			if ( i+frag_size-1 > nres ) continue;
			bool allowed( true );
			for ( Size k=0; k< frag_size; ++k ) {
				if ( ( k >           0 && pose.residue(i+k).is_lower_terminus() ) ||
						( k < frag_size-1 && pose.residue(i+k).is_upper_terminus() ) ||
						!has_element( flex_protein, i+k ) ) allowed = false;
			}
			if ( allowed ) {
				runtime_assert( pose.residue(i).is_protein() && pose.residue(i+frag_size-1).is_protein() &&
					pose.chain(i) == pose.chain(i+frag_size-1 ) );
				TR.Trace << "Inserting random " << frag_size << "-mer at position " << i <<
					" frag_depth= " << lib[i].size() << endl;
				// this will fail if there are no fragments in this window:
				protocols::frags::insert_fragment( i, i+frag_size-1, pose, lib, i );
			}
		}
	}
}


} // ns frags
} // ns protocols
