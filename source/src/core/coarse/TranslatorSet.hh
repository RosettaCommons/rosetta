// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Oliver Lange

#ifndef INCLUDED_core_coarse_TranslatorSet_hh
#define INCLUDED_core_coarse_TranslatorSet_hh

// Unit headers
// you cannot #include yourself #include <core/coarse/TranslatorSet.hh>

// Coarse headers
#include <core/coarse/Rules.fwd.hh>

// Project headers
#include <core/coarse/Translator.fwd.hh>
//#include <core/coarse/Translator.hh>
#ifdef __clang__
#include <core/coarse/Translator.hh>
#endif
#include <core/chemical/ResidueTypeSet.fwd.hh>
//#include <core/chemical/AtomTypeSet.hh>
//#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/AA.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>



// std headers
#include <ostream>
#include <map>

/* TODO:

make atoms in beads weighted
then a true centroid representation wouldn't need any change of the programming but just different rules
w=1 for CB and w=0 for all others
but all sidechain atoms in B1

*/


namespace core {
namespace coarse {

class TranslatorSet : public   utility::pointer::ReferenceCount {

	/* this is a container of Translators */
	/* for every Residue in ResidueTypeSet we find a Translator in here */
public:
	typedef std::string ResName;
	typedef std::map<ResName,TranslatorCOP> TranslatorMap;
	typedef TranslatorMap::iterator iterator;
	typedef TranslatorMap::const_iterator const_iterator;
	typedef pack::dunbrack::RotamerLibrary RotamerLibrary;
	typedef utility::pointer::access_ptr< chemical::ResidueTypeSet > ResidueTypeSetAP;
public:


	TranslatorSet( const RuleSet &rules, chemical::ResidueTypeSetCAP residue_set, ResidueTypeSetAP coarse_set);

	/// @brief prints Translator mappings
	void pretty_print( std::ostream &os ) const;

	/// @brief generate new coarseified pose from coordinates in pose_in
	void coarsify( pose::Pose &pose_out,pose::Pose const &pose_in ) const;

	/// @brief generate new coarsified residue from coordinates in fine_rsd
	conformation::ResidueOP coarsify( const conformation::Residue& fine_rsd ) const;

	// Undefined, commented out to allow pyton bindings to be builded
	//void coarsify(
	//	RotamerLibrary &rot_lib_coarse,
	//	RotamerLibrary const &rot_lib_fine
	//) const;

	//can ResName be translated ?
	bool has(ResName name) const; //silly: by construction there are no residues in residue_set but not here.

	//// @brief return the Translator for the "most generic" residue_type that is of type aa
	//// this will be the first one found now, but maybe we add the notion of "most generic"
	//// to the residue_set later on ?
	TranslatorCOP const &default_for_aa(chemical::AA aa) const;
protected: //protected since not needed apart from pretty_print....  might change
	iterator begin() { return coarse_maps_.begin(); };
	const_iterator begin() const { return coarse_maps_.begin(); };

	iterator end() { return coarse_maps_.end(); };
	const_iterator end() const { return coarse_maps_.end(); };


private:
	mutable TranslatorMap coarse_maps_;
	const chemical::ResidueTypeSetCAP residue_set_;
	const chemical::ResidueTypeSetCAP coarse_residue_set_;

};



} //namespace coarse
} // namespace core

#endif
