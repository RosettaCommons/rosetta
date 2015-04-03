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


#ifndef INCLUDED_protocols_jumping_PairingLibrary_hh
#define INCLUDED_protocols_jumping_PairingLibrary_hh

// Unit Headers

// Package Headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/MoveMap.fwd.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <map>

#include <core/scoring/dssp/PairingsList.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace jumping {


class PairingTemplate  {
	typedef utility::vector1<std::string> AtomList;
public:

	//create template with identical set of three upstream and downstream atoms to define their stubs
	PairingTemplate( std::string const& c, std::string const& s1, std::string const& s2, std::string const& s3);
	PairingTemplate( std::string const& s1, std::string const& s2, std::string const& s3);
	core::kinematics::RT rt_;
	ObjexxFCL::FArray1D_float phi;
	ObjexxFCL::FArray1D_float psi;
	ObjexxFCL::FArray1D_float omega;
	ObjexxFCL::FArray1D_char secstruct;
	AtomList atoms_downstream_;
	AtomList atoms_upstream_;
};

/// @brief returns relative orientation of chains at res1 and res2
/// this is the inner product of the respective N-C vectors.
//void compute_orientation_and_pleating( core::pose::Pose const&, core::Size res1, core::Size res2, core::Size &orientation, core::Size &pleating);

class BasePairingLibrary : public utility::pointer::ReferenceCount {
public:
	virtual ~BasePairingLibrary();

	virtual void
	create_jump_fragments(
    int const orientation,
		int const pleating,
		bool bWithTorsion,
		core::fragment::FragDataOPs &
	) const = 0;

	virtual void
	generate_jump_frags(
		core::scoring::dssp::PairingsList const& pairings,
		core::kinematics::MoveMap const& mm,
		bool bWithTorsion,
		core::fragment::FragSet& frags_accumulator
	) = 0;
};

class PairingLibrary : public BasePairingLibrary {
public:
  typedef std::vector< PairingTemplate > PairingTemplateList;
  typedef std::map< std::pair< int, int >, PairingTemplateList > PairingTemplateMap;

public:
	PairingLibrary() : num_of_pairings_(0) {}
	void read_from_file( std::string const& fn );
	void read_from_file_no_filters( std::string const& fn ); /*Version which does not assume the jump is a beta sheet */
	//void read_from_database();

	/// @brief classic rosetta++ accessor
	core::kinematics::RT
	get_random_beta_sheet_jump(
		int const orientation,
		int const pleating
	) const;

	/// @brief classic rosetta++ accessor
	core::kinematics::RT
	get_random_tmh_jump(
		int const orientation,
		int const pos1,
		int const pos2
	) const;

	void
	set_tmh_jump(
		core::pose::Pose pose,
		int const jump_number,
		int const orientation,
		int const pos1,
		int const pos2
	) const;

	/// @brief puts all jump-geometries that fit the orientation and pleating into
	/// list of FragData's. Try to reuse these FragData for different Frames that have same orientation and pleating
	/// This creates Fragments with single JumpSRFD --- PairingLibrary also stores phi/psi/omega of start and end residue
	/// use bWithTorsion = true to get FragData with BBTorsionSRFD and JumpSRFD
	/// length of single FragData is
	///   noTorsion 1
	///   withTorsion 3
	/// bWithTorsion = true length of single FragData is 3   start jump end
	void
	create_jump_fragments(
    int const orientation,
		int const pleating,
		bool bWithTorsion,
		core::fragment::FragDataOPs &
	) const;

	core::Size
	size() const {
		return num_of_pairings_;
	}

	void
	generate_jump_frags(
		core::scoring::dssp::PairingsList const & pairings,
		core::kinematics::MoveMap const & mm,
		bool bWithTorsion,
		core::fragment::FragSet & frags_accumulator
	);

private:
  PairingTemplateMap pairings_;
  core::Size num_of_pairings_;
};

class SpecificGeometryLibrary : public PairingLibrary {
public:

private:
};

/// @brief This class is thread-unsafe, though, if perhaps none of its non-const
/// functions were accessible, then it wouldn't be.
class StandardPairingLibrary :
		public PairingLibrary,
		public utility::SingletonBase< StandardPairingLibrary > {
public:
	friend class utility::SingletonBase< StandardPairingLibrary >;

private:
	StandardPairingLibrary() {}

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static StandardPairingLibrary * create_singleton_instance();

};


} //protocols
} //jumping

#endif

