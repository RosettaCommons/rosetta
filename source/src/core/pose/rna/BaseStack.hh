// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/rna/BaseStack.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_rna_BaseStack_HH
#define INCLUDED_core_pose_rna_BaseStack_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/pose/rna/BaseStack.fwd.hh>
#include <core/chemical/rna/util.hh>
#include <core/types.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace rna {


class BaseStack : public utility::pointer::ReferenceCount {
public:

	BaseStack();

	BaseStack( core::Size const & res1, core::Size const & res2,
		core::chemical::rna::BaseDoubletOrientation const & orientation,
		core::chemical::rna::BaseStackWhichSide const & which_side  ):
		res1_( res1 ),
		res2_( res2 ),
		orientation_( orientation ),
		which_side_( which_side )
	{
	};

	~BaseStack(){};

	BaseStack
	flipped() const;

	friend
	bool operator < ( BaseStack const & lhs, BaseStack const & rhs );

	friend
	std::ostream &
	operator << ( std::ostream & out, BaseStack const & s );

	void set_res1( core::Size const & setting ){ res1_ = setting; }
	core::Size res1() const { return res1_; }

	void set_res2( core::Size const & setting ){ res2_ = setting; }
	core::Size res2() const { return res2_; }

	void set_orientation( core::chemical::rna::BaseDoubletOrientation const & setting ){ orientation_ = setting; }
	core::chemical::rna::BaseDoubletOrientation orientation() const { return orientation_; }

	void set_which_side( core::chemical::rna::BaseStackWhichSide const & setting ){ which_side_ = setting; }
	core::chemical::rna::BaseStackWhichSide which_side() const { return which_side_; }

private:

	Size res1_;
	Size res2_;
	core::chemical::rna::BaseDoubletOrientation orientation_; // 1 = antiparallel; 2 = parallel
	core::chemical::rna::BaseStackWhichSide which_side_;  // 1 = residue 2 is 3' to residue1;  2 = residue 2 is 5' to residue 1

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

typedef std::pair< Real, BaseStack > EnergyBaseStack;
typedef std::list < EnergyBaseStack > EnergyBaseStackList;
typedef utility::vector1 < BaseStack > RNA_BaseStackList;

} //rna
} //pose
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pose_rna_BaseStack )
#endif // SERIALIZATION


#endif
