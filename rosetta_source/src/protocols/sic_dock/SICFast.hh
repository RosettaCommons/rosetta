// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_sic_dock_SICFast_hh
#define INCLUDED_protocols_sic_dock_SICFast_hh

#include <protocols/sic_dock/SICFast.fwd.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID_Map.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/sic_dock/xyzStripeHashPose.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace sic_dock {

class SICFast : public utility::pointer::ReferenceCount {
public:
	typedef numeric::xyzVector<core::Real> Vec;

	SICFast();

	virtual ~SICFast(){}

	void init(
		core::pose::Pose const & pose1
	);

	void init(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2
	);

	void init(
		core::pose::Pose const & pose1,
		core::id::AtomID_Map<core::Real> const & clash_atoms1,
		core::id::AtomID_Map<core::Real> const & score_atoms1
	);

	void init(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::id::AtomID_Map<core::Real> const & clash_atoms1, // currently >0 means include in clash check
		core::id::AtomID_Map<core::Real> const & clash_atoms2, // could be nice if was clash radius
		core::id::AtomID_Map<core::Real> const & score_atoms1, // weight for a "contact"
		core::id::AtomID_Map<core::Real> const & score_atoms2  // 
	);

	// return distace xmob*pose1 must move along ori to contact xfix*pose2
	// moves xmob 
	double slide_into_contact(
		core::kinematics::Stub       & xmob,
		core::kinematics::Stub const & xfix,
		Vec                            ori,
		double                       & score
	);

private:
	double CTD,CLD,CTD2,CLD2,BIN;
	xyzStripeHashPoseOP xh1c_,xh1s_;
	xyzStripeHashPoseOP xh2c_,xh2s_;
	utility::vector1<double> w1_,w2_;
};


int flood_fill3D(int i, int j, int k, ObjexxFCL::FArray3D<double> & grid, double t);

void termini_exposed(core::pose::Pose const & pose, bool & ntgood, bool & ctgood );


} // namespace sic_dock
} // namespace protocols

#endif
