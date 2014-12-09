// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


// Rosetta Headers
#include <protocols/farna/RNA_FragmentsClasses.hh>
#include <protocols/farna/RNA_SecStructInfo.hh>

#include <core/chemical/rna/util.hh>
#include <core/pose/Pose.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/StaticIndexRange.hh>

#include <ObjexxFCL/format.hh>

#include <utility/io/izstream.hh>
#include <utility/exit.hh>

#include <numeric/xyzVector.hh>

#include <core/types.hh>

// C++ headers
#include <fstream>
#include <iostream>

//Auto Headers
#include <numeric/random/random.fwd.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end



namespace protocols {
namespace farna {

	using core::Size;
	using core::Real;

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	TorsionSet::TorsionSet( Size const size ){
		torsions.dimension( core::chemical::rna::NUM_RNA_TORSIONS, SRange(0, size) );
		torsion_source_name.dimension( SRange(0, size), std::string( 4, ' ' )  );
		secstruct.dimension( SRange(0, size), 'L' );
		non_main_chain_sugar_coords_defined = false;
		size_ = size;
	}

//////////////////////////////////////////////////////////////////////
	TorsionSet &
	TorsionSet::operator =(
			TorsionSet const & src
	){
		size_ = src.size_;

		for (Size offset = 0; offset < size_; offset++){
			for (Size j = 1; j <= core::chemical::rna::NUM_RNA_TORSIONS; j++ ){
				torsions( j, offset) = src.torsions( j, offset);
			}
			torsion_source_name( offset ) = src.torsion_source_name( offset );
			secstruct( offset ) = src.secstruct( offset );
		}

		non_main_chain_sugar_coords_defined = src.non_main_chain_sugar_coords_defined;

		if (non_main_chain_sugar_coords_defined) {
			non_main_chain_sugar_coords.dimension( SRange(0,size_), 3, 3 );
			for (Size offset = 0; offset < size_; offset++){
				for (Size j = 1; j <= 3; j++ ) {
					for (Size k = 1; k <= 3; k++ ) {
						non_main_chain_sugar_coords( offset, j, k ) =
							src.non_main_chain_sugar_coords( offset, j, k );
					}
				}
			}
		}

		return *this;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	Real FragmentLibrary::get_fragment_torsion( Size const num_torsion,  Size const which_frag, Size const offset ){
		return align_torsions_[ which_frag - 1 ].torsions( num_torsion, offset) ;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	TorsionSet const FragmentLibrary::get_fragment_torsion_set( Size const which_frag ){
		return align_torsions_[ which_frag - 1 ];
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	void  FragmentLibrary::add_torsion( TorsionSet const torsion_set ){
		align_torsions_.push_back( torsion_set );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	void  FragmentLibrary::add_torsion(
																		 RNA_Fragments const & vall,
																		 Size const position,
																		 Size const size
																		 )
	{
		TorsionSet torsion_set( size );

		for (Size offset = 0; offset < size; offset++){
			for (Size j = 1; j <= core::chemical::rna::NUM_RNA_TORSIONS; j++ ){
				torsion_set.torsions( j, offset) = vall.torsions( j, position+offset);
			}
			torsion_set.torsion_source_name( offset ) = vall.name( position+offset );
			torsion_set.secstruct( offset ) = vall.secstruct( position+offset );

			//Defined non-ideal geometry of sugar ring -- to keep it closed.
			if ( vall.non_main_chain_sugar_coords_defined() ){
				torsion_set.non_main_chain_sugar_coords_defined = true;
				torsion_set.non_main_chain_sugar_coords.dimension( SRange(0,size), 3, 3 );
				for (Size j = 1; j <= 3; j++ ){
					for (Size k = 1; k <= 3; k++ ){
						torsion_set.non_main_chain_sugar_coords( offset, j, k ) =
							vall.non_main_chain_sugar_coords( position+offset, j, k );
					}
				}
			} else {
				torsion_set.non_main_chain_sugar_coords_defined = false;
			}

		}

		align_torsions_.push_back( torsion_set );
	}


	/////////////////////////////////////////////////////////////////////////////////////////////////
	Size FragmentLibrary::get_align_depth() {
		return align_torsions_.size();
	}

	///////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_Fragments::pick_fragment_library( SequenceSecStructPair const & key ){

		//Don't worry, I coded a destructor in later.
		FragmentLibrary * fragment_library_p;
		fragment_library_p = new FragmentLibrary;

		std::string const RNA_string = key.first;
		std::string const RNA_secstruct_string = key.second;

		Size const size = RNA_string.length();

		runtime_assert( RNA_string.length() == RNA_secstruct_string.length() );

		// dummy initialization.
		std::string vall_current_sequence ( RNA_string );
		std::string vall_current_secstruct( RNA_secstruct_string );

		for (Size i = 1; i <= vall_size_ - size + 1; i++ ){

			bool match( true );

			for (Size offset = 0; offset < size; offset++ ){
				vall_current_sequence [offset] = vall_sequence_ ( i + offset );
				vall_current_secstruct[offset] = vall_secstruct_( i + offset );

				if ( /*vall_is_chainbreak_( i + offset ) ||*/
						 !compare_RNA_char( vall_current_sequence[offset], RNA_string[ offset ] ) ||
						 !compare_RNA_secstruct( vall_current_secstruct[offset], RNA_secstruct_string[ offset ] ) )	{
					match = false;
					break;
				}
			}

			if (match) {
				fragment_library_p->add_torsion( *this, i, size );
			}

		}


		if ( fragment_library_p->get_align_depth() == 0  ) {
			// Problem -- need to repick with less stringent requirements?
			for (Size i = 1; i <= vall_size_ - size + 1; i++ ){

				bool match( true );

				for (Size offset = 0; offset < size; offset++ ){
					vall_current_sequence [offset] = vall_sequence_ ( i + offset );

					if ( !compare_RNA_char( vall_current_sequence[offset], RNA_string[ offset ] ) ) {
						match = false;
						break;
					}
				}

				if (match) {
					fragment_library_p->add_torsion( *this, i, size );
				}

			}
		}


		std::cout << "Picked Fragment Library for sequence " << RNA_string << " " <<
			" and sec. struct " << RNA_secstruct_string << " ... found " <<
			fragment_library_p->get_align_depth() << " potential fragments" << std::endl;

		fragment_library_pointer_map[ key ] = fragment_library_p;

	}

	///////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_Fragments::pick_random_fragment(
			 TorsionSet & torsion_set,
			 std::string const RNA_string,
			 std::string const RNA_secstruct_string,
			 Size const type /* = MATCH_YR */){

		std::string const RNA_string_local = convert_based_on_match_type( RNA_string, type );

		SequenceSecStructPair const key( std::make_pair( RNA_string_local, RNA_secstruct_string ) );

		if (! fragment_library_pointer_map.count( key ) ){
			pick_fragment_library( key );
		}

		FragmentLibraryOP fragment_library_pointer = fragment_library_pointer_map[ key ];

		Size const num_frags = fragment_library_pointer->get_align_depth();

		if (num_frags == 0) { //trouble.
			std::cout << "Fragment Library: zero fragments found for " << RNA_string_local << std::endl;
			std::cerr << "Fragment Library: zero fragments found for " << RNA_string_local << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}

		Size const which_frag = static_cast <Size> ( numeric::random::uniform() * num_frags) + 1;

		torsion_set = fragment_library_pointer->get_fragment_torsion_set( which_frag );

	}

	///////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_Fragments::pick_random_fragment(
			 TorsionSet & torsion_set,
			 core::pose::Pose & pose,
			 Size const position,
			 Size const size,
			 Size const type /* = MATCH_YR */){

		std::string const & RNA_sequence( pose.sequence() );
		std::string const & RNA_string = RNA_sequence.substr( position - 1, size );

		//Desired "secondary structure".
		// TEMPORARY HACK!!!
		//		std::string const RNA_secstruct( pose.total_residue(), 'X' );
		std::string const & RNA_secstruct( protocols::farna::get_rna_secstruct( pose ) );
		std::string const & RNA_secstruct_string = RNA_secstruct.substr( position - 1, size );

		pick_random_fragment( torsion_set, RNA_string, RNA_secstruct_string, type );

	}

	///////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_Fragments::read_vall_torsions( std::string const filename ){

		//Just read in this file once.
		static bool init ( false );
		if (init) return;
		init = true;

		///////////////////////////////////////////////////////////////
		//A bunch of vectors for temporary readin.
		//At the end, transfer all the data to FArrays for faster access.
		typedef numeric::xyzVector< Real > Vector;
		utility::vector1< utility::vector1< Real > > vall_torsions;
		utility::vector1< utility::vector1< Vector > > vall_non_main_chain_sugar_coords;
		utility::vector1< char > vall_sequence;
		utility::vector1< char > vall_secstruct;
		utility::vector1< bool > vall_is_chainbreak;
		utility::vector1< utility::vector1< bool > > vall_edge_is_base_pairing;
		utility::vector1< bool > vall_makes_canonical_base_pair;
		utility::vector1< std::string > vall_name;
		vall_non_main_chain_sugar_coords_defined_ = false;


		///////////////////////////////////////////////////////////////
		std::cout << "Reading in vall_torsions file: " <<  filename << std::endl;

		//This will check in rosetta_database first.
		utility::io::izstream vall_in( filename.c_str() );
		if ( vall_in.fail() ){
			utility_exit_with_message(  "Bad vall torsions file? " + filename );
		}

		std::string line, tag;

		char dummy_char;
		bool dummy_bool;
		Real dummy_real;
		std::string dummy_string;

		Size count( 0 );
		while (  getline( vall_in, line) ){

			std::istringstream line_stream( line );

			count++;

			line_stream >> dummy_char;
			vall_sequence.push_back( dummy_char );

			utility::vector1 < Real > dummy_vec;
			for (Size i = 1; i <= core::chemical::rna::NUM_RNA_TORSIONS; i++ ) {
				line_stream >> dummy_real;
				dummy_vec.push_back( dummy_real );
			}
			vall_torsions.push_back( dummy_vec );

			line_stream >> dummy_char;

			//In the new style fragment set... keep track of C3', C2', O4' coordinates
			// explicitly, allowing for non-ideal bond lengths and bond angles.
			if ( dummy_char == 'S' ) {
				vall_non_main_chain_sugar_coords_defined_ = true;
				utility::vector1< Vector > vecs;
				Real x,y,z;
				for (Size n = 1; n <= 3; n++ ) {
					line_stream >> x >> y >> z;
					vecs.push_back( Vector( x,y,z) );
				}
				vall_non_main_chain_sugar_coords.push_back( vecs );
				line_stream >> dummy_char;
			}

			vall_secstruct.push_back( dummy_char );

			utility::vector1 < bool > dummy_vec2;
			for (Size i = 1; i <= core::chemical::rna::NUM_EDGES; i++ ) {
				line_stream >> dummy_bool;
				dummy_vec2.push_back( dummy_bool );
			}
			vall_edge_is_base_pairing.push_back( dummy_vec2 );

			//vall_is_chainbreak_( count ) = 0.0;
			line_stream >> dummy_bool;
			vall_is_chainbreak.push_back( dummy_bool );

			//In principle could look for one more string in the vall
			// torsions file as a "name", but for now just keep track
			// of line number.
			vall_name.push_back( I( 4, count ) );

		} // line_stream

		vall_size_ = count;

		vall_in.close();

		std::cout << "Lines read from vall_torsions file: " << vall_size_ << std::endl;

		///////////////////////////////////////////////////////////////
		// Permanent storage.
		vall_torsions_.dimension ( SRange(0, core::chemical::rna::NUM_RNA_TORSIONS), vall_size_ );
		vall_sequence_.dimension ( vall_size_ );
		vall_secstruct_.dimension ( vall_size_ );
		vall_is_chainbreak_.dimension ( vall_size_ );
		vall_edge_is_base_pairing_.dimension( vall_size_, core::chemical::rna::NUM_EDGES );
		vall_name_.dimension( vall_size_ );

		if ( vall_non_main_chain_sugar_coords_defined_ ) vall_non_main_chain_sugar_coords_.dimension( vall_size_, 3, 3 );

		for (Size n = 1; n <= vall_size_; n++ ) {

			for (Size i = 1; i <= core::chemical::rna::NUM_RNA_TORSIONS; i++ ) {
				vall_torsions_( i, n ) = vall_torsions[ n ][ i ];
			}

			if (vall_non_main_chain_sugar_coords_defined_) {
				for ( Size i = 1; i <= 3; i++ ) {
					vall_non_main_chain_sugar_coords_( n, i, 1 ) = vall_non_main_chain_sugar_coords[ n ][ i ].x();
					vall_non_main_chain_sugar_coords_( n, i, 2 ) = vall_non_main_chain_sugar_coords[ n ][ i ].y();
					vall_non_main_chain_sugar_coords_( n, i, 3 ) = vall_non_main_chain_sugar_coords[ n ][ i ].z();
				}
			}

			vall_sequence_( n ) = vall_sequence[ n ];
			vall_secstruct_( n ) = vall_secstruct[ n ];
			vall_is_chainbreak_( n ) = vall_is_chainbreak[ n ];
			for (Size i = 1; i <= core::chemical::rna::NUM_EDGES; i++ ) {
				vall_edge_is_base_pairing_( n , i) = vall_edge_is_base_pairing[ n ][ i ];
			}
			vall_name_( n ) = vall_name[ n ];
		}


	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	std::string const
	RNA_Fragments::convert_based_on_match_type( std::string const RNA_string, Size const type ){

		std::string RNA_string_local = RNA_string;

		Size const size = RNA_string.length();

		//Obey orders to match exactly, match pyrimidine/purine, or match all.
		if (type == MATCH_ALL){

			for (Size i = 0; i < size; i++) 	RNA_string_local[ i ] = 'n';

		} else if ( type == MATCH_YR ) {

			for (Size i = 0; i < size; i++) {
				if (RNA_string[ i ] == 'g' || RNA_string[ i ] == 'a' ){
					RNA_string_local[ i ] = 'r';
				} else {
					runtime_assert( RNA_string[ i ] == 'u' || RNA_string[ i ] == 'c' );
					RNA_string_local[ i ] = 'y';
				}
			}

		}

		return RNA_string_local;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	RNA_Fragments::compare_RNA_char( char const char1, char const char2 ) {
		//Man this is silly, there must be a more elegant way to do this.
		if (char1 == char2) return true;
		if (char1 == 'n' || char2 == 'n') return true;
		if (char1 == 'r' && (char2 == 'a' || char2 == 'g')) return true;
		if (char1 == 'y' && (char2 == 'c' || char2 == 'u')) return true;
		if (char2 == 'r' && (char1 == 'a' || char1 == 'g')) return true;
		if (char2 == 'y' && (char1 == 'c' || char1 == 'u')) return true;
		return false;
	}

	bool
	RNA_Fragments::compare_RNA_secstruct( char const char1, char const char2 ) {
		if (char1 == char2) return true;
		if (char1 == 'X' || char2 == 'X' ) return true;
		if (char1 == 'L' && ( char2 == 'N' || char2 == 'P') ) return true;
		if (char2 == 'L' && ( char1 == 'N' || char1 == 'P') ) return true;
		return false;
	}


} //farna
} //protocols

