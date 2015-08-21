// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/jumping/PairingTemplate
/// @brief header file for ClassicAbinitio protocol
/// @details
///  from converting jumping_pairings.cc of rosetta++ into mini
///
///
///
/// @author Oliver Lange

#ifndef INCLUDED_protocols_jumping_DisulfPairingLibrary_hh
#define INCLUDED_protocols_jumping_DisulfPairingLibrary_hh

// Unit Headers
#include <protocols/jumping/DisulfPairingsList.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/MoveMap.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <map>

namespace protocols {
namespace jumping {

class DisulfTemplate  {
	typedef utility::vector1<std::string> AtomList;
public:

	//create template with identical set of three upstream and downstream atoms to define their stubs
	DisulfTemplate( std::string const& c, std::string const& s1, std::string const& s2, std::string const& s3);
	DisulfTemplate( std::string const& s1, std::string const& s2, std::string const& s3);
	core::kinematics::RT rt_;
	ObjexxFCL::FArray1D_float phi;
	ObjexxFCL::FArray1D_float psi;
	ObjexxFCL::FArray1D_float omega;
	ObjexxFCL::FArray1D_char secstruct;
	AtomList atoms_downstream_;
	AtomList atoms_upstream_;
};


class BaseDisulfPairingLibrary : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~BaseDisulfPairingLibrary();

	//Not currently using base class as copied from PairingLibrary, it's hanging
	//around because I plan to build into it

	//virtual void
	//create_jump_fragments(
	//  int const orientation,
	// int const pleating,
	// bool bWithTorsion,
	// core::fragment::FragDataOPs &
	//) const = 0;

	//virtual void
	//generate_jump_frags(
	// PairingsList const& pairings,
	// core::kinematics::MoveMap const& mm,
	// bool bWithTorsion,
	// core::fragment::FragSet& frags_accumulator
	//) = 0;
};

class DisulfPairingLibrary : public BaseDisulfPairingLibrary {
	typedef std::vector< DisulfTemplate > DisulfTemplateList;
	typedef std::map< std::pair< int, int >, DisulfTemplateList > DisulfTemplateMap;
public:
	DisulfPairingLibrary():num_of_pairings_(0) {};
	void read_from_file( std::string const& fn);
	//void read_from_file_no_filters( std::string const& fn); /*Version which does not assume the jump is a beta sheet */
	//void read_from_database();

	void
	create_jump_fragments(
		bool bWithTorsion,
		core::fragment::FragDataOPs& frags
	) const;

	core::Size
	size() const {
		return num_of_pairings_;
	}

	void
	generate_jump_frags(
		DisulfPairingsList const & pairings,
		core::kinematics::MoveMap const & mm,
		bool bWithTorsion,
		core::fragment::FragSet & frags_accumulator
	) const; // APL TESTING -- CAN THIS BE CONST?

private:
	DisulfTemplateMap pairings_;
	DisulfTemplateList all_pairings_;
	core::Size num_of_pairings_;
};

/// @brief The %StandardDisulfPairingsLibrary initializes itself in its constructor
/// from the sampling/disulfide_jump_database_wip.dat file in the database.  Users
/// should not in any circumstance invoke any of its non-const methods that are defined
/// in the subclass.  This class ought to be reworked to ensure that it is threadsafe.
class StandardDisulfPairingLibrary :
	public DisulfPairingLibrary,
	public utility::SingletonBase< StandardDisulfPairingLibrary >
{
public:
	friend class utility::SingletonBase< StandardDisulfPairingLibrary >;

private:
	StandardDisulfPairingLibrary() {};

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static StandardDisulfPairingLibrary * create_singleton_instance();
private:
	/// @brief static data member holding pointer to the singleton class itself

};


} //protocols
} //jumping

#endif

