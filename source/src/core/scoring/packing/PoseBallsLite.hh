// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packing/PoseBallsLite.hh
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_packing_PoseBallsLite_hh
#define INCLUDED_core_scoring_packing_PoseBallsLite_hh


// Project headers
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>

#include <core/scoring/packing/Ball.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace packing {

class PoseBallsLite : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~PoseBallsLite();

	// hmode = 0 --> no hydrogens
	// hmode = 1 --> polar H's
	// hmode = 2 --> all H's
	// in all cases, H area is added to the "parent" atom
	PoseBallsLite( core::pose::Pose const & pose, core::Size Hmode = 0, bool ignore_water = true );
	PoseBallsLite( core::pose::Pose const & pose, core::id::AtomID_Mask const & whichatoms );

	inline core::Size const & nballs() const {
		return nballs_;
	}

	inline core::Size const & id_to_index( core::id::AtomID const & aid ) const {
		return id_to_index_[aid];
	}

	inline core::id::AtomID const & index_to_id( core::Size const & index ) const {
		return index_to_id_[index];
	}

	inline Ball const & ball( core::Size const & index ) const {
		return balls_[index];
	}

	inline Ball const & ball( core::id::AtomID const & id ) const {
		return balls_[ id_to_index_[id] ];
	}

	inline Ball & ball( core::Size const & index ) { return balls_[index];	}
	inline Ball & ball( core::id::AtomID const & id ) { return balls_[ id_to_index_[id] ];	}

	inline core::Size const & atom_num( core::Size const & index ) const {
		return atom_num_[index];
	}
	inline core::Size const & res_num( core::Size const & index ) const {
		return res_num_[index];
	}


private:

	// core::pose::Pose pose_;
	core::Size nballs_;
	core::id::AtomID_Map<Size> id_to_index_;
	utility::vector1<core::id::AtomID> index_to_id_;
	utility::vector1<Ball> balls_;
	utility::vector1<core::Size>  atom_num_;
	utility::vector1<core::Size>  res_num_;

};

template< class T >
void
initialize_AtomID_Map( core::id::AtomID_Map<T> & map, PoseBallsLite const & pb ) {
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		Size res_num = pb.res_num(i);
		Size atom_num = pb.atom_num(i);
		map.resize( std::max(map.size(),res_num) );
		map.resize( res_num, std::max(map.n_atom(res_num),atom_num) );
	}
}


} // namespace packing
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_packing_PoseBallsLite_HH
