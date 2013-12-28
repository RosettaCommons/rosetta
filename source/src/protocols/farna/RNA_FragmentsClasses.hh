// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//  CVS information:
//  $Revision: 1.1.2.1 $
//  $Date: 2005/11/07 21:05:35 $
//  $Author: rhiju $
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_rna_RNA_FragmentsClasses_hh
#define INCLUDED_protocols_rna_RNA_FragmentsClasses_hh

#include <protocols/farna/RNA_FragmentsClasses.fwd.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>
#include <map>
#include <vector>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end



/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
// Goal: to make a fragment object that can choose fragments
// "on the fly" for RNA ab inito folding.
//
// After reading in a set of torsions from, e.g., the ribosome crystal structure,
//  should be able to generate fragments of size 1, 2, or 3, with
//  exact sequence matches, partial Y/R matches, or ignoring sequence.
//
namespace protocols {
namespace farna {

    /// name MATCH_ENUM added for PyRosetta compatability
	enum MATCH_ENUM { MATCH_ALL /* 0 */, MATCH_YR /* 1 */, MATCH_EXACT /* 2 */};

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	class TorsionSet {
	public:
		TorsionSet & operator =( TorsionSet const & src );
		TorsionSet( core::Size const size );

		//make this private?
		FArray2D <core::Real>   torsions;  // dimensions: (NUM_RNA_TORSIONS, SRange(0, size) );
		FArray1D <std::string> torsion_source_name;  // dimensions: ( SRange(0, size) );
		FArray1D <char>   secstruct;

		FArray3D <core::Real>   non_main_chain_sugar_coords;
		bool  non_main_chain_sugar_coords_defined;

		inline
		core::Size get_size() const { return size_; }

	private:
		core::Size size_;

	};

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	class FragmentLibrary : public utility::pointer::ReferenceCount  {
	public:

		//constructor!
		//FragmentLibrary();

		//destructor -- necessary?
		//~FragmentLibrary();

		core::Real get_fragment_torsion(
																		core::Size const num_torsion,
																		Size const which_frag,
																		core::Size const offset );

		TorsionSet const get_fragment_torsion_set( core::Size const which_frag );

		void  add_torsion( TorsionSet const torsion_set );

		void  add_torsion(
					RNA_Fragments const & vall,
					core::Size const position,
					core::Size const size
											);

		core::Size get_align_depth();

	private:
		std::vector< TorsionSet > align_torsions_;

	};

	/////////////////////////////////////////////////////////////////////////////////////////////////
	typedef utility::pointer::owning_ptr< FragmentLibrary > FragmentLibraryOP;
	typedef std::pair< std::string, std::string > SequenceSecStructPair;
	typedef std::map< SequenceSecStructPair, FragmentLibraryOP >  FragmentLibraryPointerMap;

	/////////////////////////////////////////////////////////////////////////////////////////////////
	class RNA_Fragments : public utility::pointer::ReferenceCount {
	public:
		//Constructor -- needs vall_torsions_file to get started.
		// RNA_Fragments();
		RNA_Fragments( std::string const filename ){
			read_vall_torsions( filename );
		}

		//We need an explicit destructor, I think, because each
		// of the fragment libraries is "new". Not anymore, with
		// the owning pointer business!
		//		~RNA_Fragments(){
		//			FragmentLibraryPointerMap::iterator it;
		//			for (it = fragment_library_pointer_map.begin();
		//					 it != fragment_library_pointer_map.end();
		//					 it++ ){
		//				delete it->second;
		//			}
		//		}

		void read_vall_torsions( std::string const filename );

		//Probably the only thing that will actually get called publicly:
		void
		pick_random_fragment(
					TorsionSet & torsion_set,
          core::pose::Pose & pose,
					core::Size const position,
					core::Size const size,
					core::Size const type = MATCH_YR );

		core::Real
		torsions( core::Size const & i, core::Size const & j ) const { return vall_torsions_( i, j ); }

		std::string
		name( core::Size const & i ) const { return vall_name_( i ); }

		char
		secstruct( core::Size const & i ) const { return vall_secstruct_( i ); }

		bool
		non_main_chain_sugar_coords_defined() const { return vall_non_main_chain_sugar_coords_defined_; }

		core::Real
		non_main_chain_sugar_coords( core::Size const & i, core::Size const & j, core::Size const & k ) const{ return vall_non_main_chain_sugar_coords_( i, j, k);}

	private:

		void
		pick_random_fragment(
					TorsionSet & torsion_set,
					std::string const RNA_string,
					std::string const RNA_secstruct_string,
					core::Size const type = MATCH_YR );

		std::string const
		convert_based_on_match_type( std::string const RNA_string, core::Size const type );

		bool
		compare_RNA_char( char const char1, char const char2 );

		bool
		compare_RNA_secstruct( char const char1, char const char2 );


	private:

		// Probably should make following "vall" stuff a different object
		// and, come on, these could be a vector to save memory!
		FArray2D <core::Real> vall_torsions_;
		FArray3D <core::Real> vall_non_main_chain_sugar_coords_;
		FArray1D <char>  vall_sequence_;
		FArray1D <bool>  vall_is_chainbreak_;
		FArray2D <bool>  vall_edge_is_base_pairing_;
		FArray1D <bool>  vall_makes_canonical_base_pair_;
		FArray1D <char>  vall_secstruct_;
		FArray1D <std::string>  vall_name_;
		core::Size vall_size_;
		bool vall_non_main_chain_sugar_coords_defined_;

		// Need to hold on to some fragment libraries. These
		// will be picked on the fly when the code requires them.
		//   Indexed by sequence, e.g., AAA, AGA, GUA ... or even RYR ...  or even NNN (totally generic!)
		FragmentLibraryPointerMap fragment_library_pointer_map;

		void pick_fragment_library( SequenceSecStructPair const & key );

		void pick_random_fragment( FArray1D <core::Real> & RNA_torsions, std::string const RNA_string );

	};

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////


} //farna
} //protocols

#endif
