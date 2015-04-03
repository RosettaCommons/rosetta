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
#include <protocols/sic_dock/types.hh>

#include <utility/vector1.hh>
#include <core/id/AtomID_Map.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/pose/xyzStripeHashPose.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/kinematics/Stub.fwd.hh>


namespace protocols {
namespace sic_dock {

class SICFast : public utility::pointer::ReferenceCount {
public:
	typedef numeric::xyzVector<platform::Real> Vec;

	SICFast();
	SICFast(core::Real clash_dis);

	virtual ~SICFast();

	void init(
		core::pose::Pose const & pose1
	);

	void init(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2
	);

	void init(
		core::pose::Pose const & pose1,
		core::id::AtomID_Map<platform::Real> const & clash_atoms1
	);

	void init(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::id::AtomID_Map<platform::Real> const & clash_atoms1, // currently >0 means include in clash check
		core::id::AtomID_Map<platform::Real> const & clash_atoms2  // could be nice if was clash radius
	);

	// return distace xmob*pose1 must move along ori to contact xfix*pose2
	double slide_into_contact(
		Xform const & xmob,
		Xform const & xfix,
		Vec            ori
	) const;

	double slide_into_contact(
		Xforms const & xmob,
		Xforms const & xfix,
		Vec             ori
	) const;

	// convenience adaptor
	double slide_into_contact_DEPRICATED(
		core::kinematics::Stub const & xmob,
		core::kinematics::Stub const & xfix,
		Vec                            ori
	) const;

private:
	double CLD, CLD2, BIN;
	// KAB - below variables commented out (-Wunused-private-field) on 2014-09-11
	// double CTD, CTD2;

	core::pose::xyzStripeHashPose *h1_,*h2_;
};


} // namespace sic_dock
} // namespace protocols

#endif
