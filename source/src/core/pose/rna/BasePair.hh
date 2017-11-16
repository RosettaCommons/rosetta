// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/rna/BasePair.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_rna_BasePair_HH
#define INCLUDED_core_pose_rna_BasePair_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/pose/rna/BasePair.fwd.hh>
#include <core/chemical/rna/util.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace rna {

class BasePair : public utility::pointer::ReferenceCount {

public:

	BasePair( Size const res1 = 0, Size const res2 = 0,
		core::chemical::rna::BaseEdge const edge1 = core::chemical::rna::ANY_BASE_EDGE,
		core::chemical::rna::BaseEdge const edge2 = core::chemical::rna::ANY_BASE_EDGE,
		core::chemical::rna::BaseDoubletOrientation const orientation = core::chemical::rna::ANY_BASE_DOUBLET_ORIENTATION );

	~BasePair(){}

	BasePair
	flipped() const;

	void
	print_info( std::ostream & out = std::cout ) const;

	friend
	bool operator < ( BasePair const & lhs, BasePair const & rhs );

	friend
	bool operator == ( BasePair const & lhs, BasePair const & rhs );

	friend
	std::ostream &
	operator << ( std::ostream & out, BasePair const & s );

	void set_res1( core::Size const & setting ){ res1_ = setting; }
	core::Size res1() const { return res1_; }

	void set_res2( core::Size const & setting ){ res2_ = setting; }
	core::Size res2() const { return res2_; }

	void set_edge1( core::chemical::rna::BaseEdge const & setting );
	core::chemical::rna::BaseEdge edge1() const { return edge1_; }

	void set_edge2( core::chemical::rna::BaseEdge const & setting );
	core::chemical::rna::BaseEdge edge2() const { return edge2_; }

	void set_orientation( core::chemical::rna::BaseDoubletOrientation const & setting );
	core::chemical::rna::BaseDoubletOrientation orientation() const { return orientation_; }
	core::chemical::rna::LW_BaseDoubletOrientation LW_orientation() const { return LW_orientation_; }

private:

	void
	derive_LW_orientation();

private:

	Size res1_;
	Size res2_;
	core::chemical::rna::BaseEdge edge1_;
	core::chemical::rna::BaseEdge edge2_;
	core::chemical::rna::BaseDoubletOrientation orientation_; // 1 = antiparallel; 2 = parallel

	// following is derived based on BasePairOrientation
	core::chemical::rna::LW_BaseDoubletOrientation LW_orientation_; // 1 = cis; 2 = trans


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

typedef std::pair< Real, BasePair > EnergyBasePair;
typedef std::list < EnergyBasePair > EnergyBasePairList;

} //rna
} //pose
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pose_rna_BasePair )
#endif // SERIALIZATION


#endif
