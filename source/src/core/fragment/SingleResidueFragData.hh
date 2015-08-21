// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragData.hh
/// @brief  A fragment as list of SingleResidue Data
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///
#ifndef INCLUDED_core_fragment_SingleResidueFragData_HH
#define INCLUDED_core_fragment_SingleResidueFragData_HH

// Unit Headers
#include <core/fragment/SingleResidueFragData.fwd.hh>

// Package Headers
#include <core/fragment/Frame.fwd.hh>

// Project Headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>


// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers

#include <utility/vector1.hh>


namespace core {
namespace fragment {


typedef utility::vector1 < Size > PositionList;

/// @brief Base class for SRFD classes
/// @detail Instances of  SRFD classes contain information on specific dofs in
/// a single residue or a jump connected to a residue
/// The classes' apply method will now how to implement the specified dofs in the give pose
/// at the given residue position


/// TODO: change SRFD interface such that apply is called like this
/// apply( pose, inframe_pos (1..length), Frame const& )
/// jumpFrags can then ask the Frame for their upstream residue and check that if they
/// are at position 1 that they do nothgin...
/// can have basic implementation that translates this apply into the old apply
// only JumpSRFD needs to overload this function
// scrap continuous apply of FragData ?!


class SingleResidueFragData : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~SingleResidueFragData();
	SingleResidueFragData( char sequence )
	: sequence_( sequence )
	{};

	SingleResidueFragData(  )
	: sequence_( 'X' )
	{};

	virtual
	SingleResidueFragDataOP
	clone() const = 0;

	/// @brief create a new instance of this object
	virtual
	SingleResidueFragDataOP create() const = 0;

	/// @brief insert fragment_data into pose at position seq_pos
	virtual bool apply( pose::Pose &, Size seq_pos ) const = 0;

	/// @brief insert fragment_data into pose at position seq_pos for dofs that are
	///  allowed to move in the movemap
	/// @return True if operation succeeds, False otherwise.
	virtual bool apply( kinematics::MoveMap const & movemap, pose::Pose & pose, Size const seqpos ) const = 0;

	/// @brief insert fragment_data sec-struct into ss-string at position seq_pos
	virtual bool apply_ss( std::string&, Size seq_pos ) const = 0;

	/// @brief insert fragment_data into pose at position seq_pos
	virtual bool steal( pose::Pose const&, Size seq_pos ) = 0;

	/// @brief insert fragment_data into pose at position given by Frame.seqpos( intra_frame_pos );
	virtual bool apply(pose::Pose &, Size const intra_frame_pos, Frame const & ) const;

	/// @brief insert fragment_data into pose at position given by Frame.seqpos( intra_frame_pos )
	///  for dofs that are allowed to move in the MoveMap
	virtual bool apply( kinematics::MoveMap const & movemap, pose::Pose & pose, Size const intra_frame_pos, Frame const & frame ) const;

	/// @brief insert fragment_data sec-struct into ss-string at position seq_pos
	virtual bool apply_ss( std::string&, Size intra_frame_pos, Frame const& ) const;

	/// @brief insert fragment_data into pose at position seq_pos
	virtual bool steal( pose::Pose const&, Size intra_frame_pos, Frame const& );

	/// @brief check weather SRFD applies to same dofs and is of same type
	virtual bool
	is_compatible( SingleResidueFragData const& ) const = 0;

	/// @brief check weather dofs can be moved
	virtual bool
	is_applicable( kinematics::MoveMap const&, Size intra_frame_pos, Frame const& ) const;

	/// @brief check whether dofs can be moved
	virtual bool
	is_applicable( kinematics::MoveMap const&, Size pos ) const = 0;


	void
	set_sequence( char const sequence )
	{
		sequence_ = sequence;
	}
	///
	char
	sequence() const
	{
		return sequence_;
	};

	virtual char secstruct() const {
		// this shouldn't be called for an SRFD that does not have secstruct
		debug_assert( 0 );
		return 'X';
	}

	virtual
	void show( std::ostream &out ) const;

	/// @brief Default implementation: noop
	virtual
	void read_data( std::istream& );


	virtual
	std::string type() const;

	static std::string _static_type_name();

protected:
	char sequence_;
};

inline std::ostream& operator<< ( std::ostream& out, SingleResidueFragData const& srfd ) {
	srfd.show( out );
	return out;
}

inline std::istream& operator>> ( std::istream& in, SingleResidueFragData& srfd ) {
	srfd.read_data( in );
	return in;
}

}
}

#endif
