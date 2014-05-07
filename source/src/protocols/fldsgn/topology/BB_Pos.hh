// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/topology/BB_Pos.hh
/// @brief bb_pos class header
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_topology_BB_Pos_hh
#define INCLUDED_protocols_fldsgn_topology_BB_Pos_hh

/// Unit headers
#include <protocols/fldsgn/topology/BB_Pos.fwd.hh>

/// Package headers
#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pose/Pose.fwd.hh>

/// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>

/// Numeric headers
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace fldsgn {
namespace topology {

//////////////////////////////////////////////////////////////////////////////////////////////////////
class BB_Pos {


	typedef core::Size Size;
	typedef core::Vector Vector; //basically a xyzVector
	typedef core::pose::Pose Pose;


public:

	void
	resize( Size const nres );

	void
	clear();

	void
	take_coordinates_from_pose( Pose const & pose );

	/// @details accessor for N's coordinate on residue i, requires take_coordinates_from_pose have been recently called.
	Vector const &
	N( Size const i ) const
	{
		return N_[i];
	}

	/// @details accessor for CA's coordinate on residue i, requires take_coordinates_from_pose have been recently called.
	Vector const &
	CA( Size const i ) const
	{
		return CA_[i];
	}

	/// @details accessor for CB's coordinate on residue i, requires take_coordinates_from_pose have been recently called.
	Vector const &
	CB( Size const i ) const
	{
		return CB_[i];
	}


	/// @details accessor for C's coordinate on residue i, requires take_coordinates_from_pose have been recently called.
	Vector const &
	C( Size const i ) const
	{
		return C_[i];
	}

	/// @details accessor for O's coordinate on residue i, requires take_coordinates_from_pose have been recently called.
	Vector const &
	O( Size const i ) const
	{
		return O_[i];
	}

	/// @
	Size size() const
	{
		return residue_types_.size();
	}

private:

	bool bbindices_up_to_date( Pose const & pose ) const;
	void update_indices( Pose const & pose );

private: // DATA

	utility::vector1< Vector > N_;
	utility::vector1< Vector > CA_;
	utility::vector1< Vector > CB_;
	utility::vector1< Vector > C_;
	utility::vector1< Vector > O_;

	/// Residue types must match those of the pose for the indices
	/// to match.
	utility::vector1< core::chemical::ResidueType const * > residue_types_;
	utility::vector1< Size > N_index_;
	utility::vector1< Size > CA_index_;
	utility::vector1< Size > CB_index_;
	utility::vector1< Size > C_index_;
	utility::vector1< Size > O_index_;

};

} // ns topology
} // ns fldsgn
} // ns devel

#endif
