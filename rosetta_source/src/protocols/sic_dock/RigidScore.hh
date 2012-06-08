// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_sic_dock_RigidScore_hh
#define INCLUDED_protocols_sic_dock_RigidScore_hh

#include <protocols/sic_dock/RigidScore.fwd.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID_Map.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/sic_dock/xyzStripeHashPose.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace sic_dock {

class RigidScore : public utility::pointer::ReferenceCount {
protected:
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef core::kinematics::Stub Stub;
	typedef numeric::xyzVector<Real> Vec;
	typedef numeric::xyzMatrix<Real> Mat;
	typedef utility::vector1<Vec> Vecs;
	typedef utility::vector1<Real> Reals;
public:

	virtual ~RigidScore();

	virtual
	void
	setup_for_scoring(
		Pose const & pose1
	) = 0;

	virtual
	void
	setup_for_scoring(
		Pose const & pose1,
		Pose const & pose2
	) = 0;

	virtual
	core::Real
	score(
		Stub const & x1,
		Stub const & x2
	) const = 0;

};


class CBScore : public RigidScore {
public:
	CBScore();
	virtual ~CBScore(){}
	void setup_for_scoring(	Pose const & pose1 );
	void setup_for_scoring( Pose const & pose1,	Pose const & pose2 );
	core::Real score( Stub const & x1, Stub const & x2 );
private:
	xyzStripeHashPoseCOP hash1,hash2;
};

class Linkerscore : public RigidScore {
public:
	Linkerscore();
	virtual ~Linkerscore(){}
	void setup_for_scoring(	Pose const & pose1 );
	void setup_for_scoring( Pose const & pose1,	Pose const & pose2 );
	core::Real score( Stub const & x1, Stub const & x2 );
private:
	// loophash crap
};

class EdgeStandScore : public RigidScore {
public:
	EdgeStandScore();
	virtual ~EdgeStandScore(){}
	void setup_for_scoring(	Pose const & pose1 );
	void setup_for_scoring( Pose const & pose1,	Pose const & pose2 );
	core::Real score( Stub const & x1, Stub const & x2 );
private:
	Vecs donors,acceptors;
};

class HelixScore : public RigidScore {
public:
	HelixScore();
	virtual ~HelixScore(){}
	void setup_for_scoring(	Pose const & pose1 );
	void setup_for_scoring( Pose const & pose1,	Pose const & pose2 );
	core::Real score( Stub const & x1, Stub const & x2 );
private:
};

class BuriedPolarScore : public RigidScore {
public:
	BuriedPolarScore();
	virtual ~BuriedPolarScore(){}
	void setup_for_scoring(	Pose const & pose1 );
	void setup_for_scoring( Pose const & pose1,	Pose const & pose2 );
	core::Real score( Stub const & x1, Stub const & x2 );
private:
	Vecs polars;
};

} // namespace sic_dock
} // namespace protocols

#endif
