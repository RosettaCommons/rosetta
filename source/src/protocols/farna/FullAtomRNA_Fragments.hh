
//Auto using namespaces
//Auto using namespaces end
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
/// @author Rhiju Das

#ifndef INCLUDED_protocols_rna_FullAtomRNA_Fragments_HH
#define INCLUDED_protocols_rna_FullAtomRNA_Fragments_HH

#include <protocols/farna/FullAtomRNA_Fragments.fwd.hh>
#include <protocols/farna/RNA_Fragments.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/farna/RNA_MatchType.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>
#include <map>
#include <vector>

#include <utility/vector1.hh>


//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS


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


	class FullAtomRNA_Fragments; // defined below.

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	class FragmentLibrary : public utility::pointer::ReferenceCount  {
	public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~FragmentLibrary();

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
					FullAtomRNA_Fragments const & vall,
					core::Size const position,
					core::Size const size
											);

		core::Size get_align_depth();

	private:
		std::vector< TorsionSet > align_torsions_;

	};

	/////////////////////////////////////////////////////////////////////////////////////////////////
	typedef utility::pointer::shared_ptr< FragmentLibrary > FragmentLibraryOP;
	typedef std::pair< std::string, std::string > SequenceSecStructPair;
	typedef std::map< SequenceSecStructPair, FragmentLibraryOP >  FragmentLibraryPointerMap;

	/////////////////////////////////////////////////////////////////////////////////////////////////
	class FullAtomRNA_Fragments : public RNA_Fragments {
	public:
		//Constructor -- needs vall_torsions_file to get started.
		// RNA_Fragments();
		FullAtomRNA_Fragments( std::string const & filename );

		~FullAtomRNA_Fragments(){}

		//Probably the only thing that will actually get called publicly:
		virtual void
		apply_random_fragment(
          core::pose::Pose & pose,
					core::Size const position,
					core::Size const size,
					core::Size const type,
				toolbox::AllowInsertOP allow_insert );

		virtual bool
		is_fullatom();

		void read_vall_torsions( std::string const filename );


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

		void
		pick_random_fragment(
					TorsionSet & torsion_set,
          core::pose::Pose & pose,
					core::Size const position,
					core::Size const size,
					core::Size const type = MATCH_YR );

		void
		insert_fragment(
										core::pose::Pose & pose,
										Size const position,
										protocols::farna::TorsionSet const & torsion_set,
									toolbox::AllowInsertOP allow_insert );

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
