// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Alex Ford (fordas@u.washington.edu)
//
#ifndef INCLUDED_protocols_hotspot_hashing_StubGenerator_hh
#define INCLUDED_protocols_hotspot_hashing_StubGenerator_hh

// Unit headers
#include <protocols/hotspot_hashing/StubGenerator.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace hotspot_hashing {

class StubGenerator
{
public:
	typedef  numeric::xyzMatrix< core::Real > Matrix;
	typedef  numeric::xyzVector< core::Real > Vector;

public:
	static core::conformation::ResidueOP getStubByName( std::string name );

	static void placeResidueAtTransform( core::pose::Pose & pose, core::conformation::ResidueCOP sourceResidue, core::kinematics::Stub transform, core::Size & residuejumpindex, core::Size & residueindex);

	static void placeResidueOnPose(core::pose::Pose & pose, core::conformation::ResidueCOP residue);

	//@brief Moves residue into the transform's reference frame via local2global.
	static void moveIntoStubFrame(core::conformation::ResidueOP residue, core::kinematics::Stub transform);

	//@brief Moves residue from transform's reference frame via global2local.
	static void moveFromStubFrame(core::conformation::ResidueOP residue, core::kinematics::Stub transform);

	//@brief Returns a stub generated from the residue's orient atoms.
	//
	//Used in recapitulating residue orientations when the residue type is known.
	static core::kinematics::Stub residueStubOrientFrame(core::conformation::ResidueCOP const residue);

	//@brief Returns a stub generated from the residue's backbone and centroid location.
	//
	// Used in orienting residues of varying type in a generation direction of "interest".
	//
	// Transform aligns CA atom to <0, 0, 0>
	// CA->SC heavyatom centroid vector along <1,0,0>
	// CA->C vector on the XY plane (<CA->C> * <0,0,1> == 0)
	static core::kinematics::Stub residueStubCentroidFrame(core::conformation::ResidueCOP const residue);

	static Vector residueStubCentroid(core::conformation::ResidueCOP const residue);
private:
	StubGenerator();
};

} // namespace hotspot_hashing
} // namespace protocols

#endif
