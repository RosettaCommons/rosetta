// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Rhiju Das

#ifndef INCLUDED_protocols_rna_FullAtomRNA_Fragments_HH
#define INCLUDED_protocols_rna_FullAtomRNA_Fragments_HH

#include <protocols/rna/denovo/fragments/FullAtomRNA_Fragments.fwd.hh>
#include <protocols/rna/denovo/fragments/RNA_Fragments.hh>
#include <protocols/rna/denovo/fragments/RNA_FragmentHomologyExclusion.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/rna/denovo/fragments/RNA_MatchType.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>
#include <map>
#include <set>
#include <vector>

#include <utility/vector1.hh>


//Auto using namespaces
// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
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
namespace rna {
namespace denovo {
namespace fragments {

class FullAtomRNA_Fragments; // defined below.

/////////////////////////////////////////////////////////////////////////////////////////////////

typedef std::tuple< std::string, std::string, RNA_FragmentHomologyExclusion, utility::vector1< SYN_ANTI_RESTRICTION > > FragmentLibraryPointerKey;
typedef std::map< FragmentLibraryPointerKey, FragmentLibraryOP >  FragmentLibraryPointerMap;

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
		RNA_FragmentHomologyExclusionCOP const & homology_exclusion,
		toolbox::AtomLevelDomainMapCOP atom_level_domain_map ) const;

	virtual bool
	is_fullatom();

	std::string
	name( core::Size const & i ) const { return vall_name_( i ); }

	char
	secstruct( core::Size const & i ) const { return vall_secstruct_( i ); }

	bool
	non_main_chain_sugar_coords_defined() const { return vall_non_main_chain_sugar_coords_defined_; }

	core::Real
	non_main_chain_sugar_coords( core::Size const & i, core::Size const & j, core::Size const & k ) const { return vall_non_main_chain_sugar_coords_( i, j, k);}

	core::Real
	torsions( core::Size const & i, core::Size const & j ) const { return vall_torsions_( i, j ); }

	FragmentLibraryOP
	get_fragment_library_pointer(
		std::string const & RNA_string,
		std::string const & RNA_secstruct_string,
		RNA_FragmentHomologyExclusionCOP const & homology_exclusion,
		utility::vector1< SYN_ANTI_RESTRICTION > const & restriction = utility::vector1< SYN_ANTI_RESTRICTION >(),
		Size const type = MATCH_YR ) const;

	void
	insert_fragment(
		core::pose::Pose & pose,
		Size const position,
		TorsionSet const & torsion_set,
		toolbox::AtomLevelDomainMapCOP atom_level_domain_map ) const;

private:

	void read_vall_torsions( std::string const & filename );

	void
	pick_random_fragment(
		TorsionSet & torsion_set,
		std::string const & RNA_string,
		std::string const & RNA_secstruct_string,
		RNA_FragmentHomologyExclusionCOP const & homology_exclusion,
		utility::vector1< SYN_ANTI_RESTRICTION > const & restriction = utility::vector1< SYN_ANTI_RESTRICTION >(),
		core::Size const type = MATCH_YR ) const;

	void
	pick_random_fragment(
		TorsionSet & torsion_set,
		core::pose::Pose & pose,
		core::Size const position,
		core::Size const size,
		RNA_FragmentHomologyExclusionCOP const & homology_exclusion,
		core::Size const type = MATCH_YR ) const;

	void pick_fragment_library( FragmentLibraryPointerKey const & key ) const;

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
	mutable FragmentLibraryPointerMap fragment_library_pointer_map;

};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


} //fragments
} //denovo
} //rna
} //protocols

#endif
