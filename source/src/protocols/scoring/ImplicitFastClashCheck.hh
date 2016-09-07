// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/scoring/ImplicitFastClashCheck.hh
/// @brief  does implicit fast clash checking WRT the provided pose
/// @author Will Sheffler (will@sheffler.me)

#ifndef INCLUDED_protocols_scoring_ImplicitFastClashCheck_hh
#define INCLUDED_protocols_scoring_ImplicitFastClashCheck_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzTriple.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <platform/types.hh>

#include <core/kinematics/Stub.fwd.hh>
#include <utility/io/ozstream.fwd.hh>


namespace protocols {
namespace scoring {

class ImplicitFastClashCheck : public utility::pointer::ReferenceCount {

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~ImplicitFastClashCheck() override;

	ImplicitFastClashCheck(
		core::pose::Pose const & pose_in,
		core::Real clash_dis
	);

	ImplicitFastClashCheck(
		core::pose::Pose const & pose_in,
		core::Real clash_dis,
		utility::vector1<core::Size> ignore
	);

	ImplicitFastClashCheck(
		utility::vector1<core::pose::Pose> const & poses_in,
		core::Real clash_dis,
		utility::vector1<core::Size> ignore
	);

	void
	init_clash_check(
		utility::vector1<core::pose::Pose> const & poses,
		core::Real neighbor_cutoff,
		utility::vector1<core::Size> ignore
	);

	bool
	clash_check(
		numeric::xyzVector<core::Real> const & pp
	) const;

	platform::uint
	clash_count(
		numeric::xyzVector<core::Real> const & pp
	) const;

	bool
	clash_check(
		numeric::xyzVector<core::Real> const & pp,
		core::Size resno
	) const;

	// bool
	// clash_check(
	//  core::pose::Pose const & pose,
	//  core::Size refrsd
	// ) const;
	//
	// bool
	// clash_check(
	//  core::kinematics::Stub const & stub,
	//  numeric::xyzVector<core::Real> pos
	// ) const;
	//
	bool
	clash_check_trimer(
		core::pose::Pose const & pose,
		Size refrsd
	) const;

	void
	dump_debug_pdb(
		utility::io::ozstream & out,
		core::kinematics::Stub const & stub,
		char chain = 'Z'
	) const;

	void
	dump_debug_pdb(
		std::string const & fname,
		core::kinematics::Stub const & stub,
		char chain = 'Z'
	) const;

	bool
	clash_check_test( numeric::xyzVector<core::Real> const & pp ) const;

	core::Size size() {
		return points_.size();
	}

private:

	core::pose::PoseCOP pose_;

	utility::vector1<numeric::xyzVector<core::Real> > points_;

	utility::vector1<core::Size> resno_;

	utility::vector1<core::Size> atomno_;

	ObjexxFCL::FArray3D< utility::vector1<unsigned int> > cubes_;

	numeric::xyzVector<core::Real> bbl_;

	numeric::xyzTriple< core::Size > cube_dim_;

	core::Real side_inv_, neighbor_cutoff_, neighbor_cutoff_sq_;

	ObjexxFCL::FArray3D< utility::vector1<numeric::xyzVector<core::Real> > > cubes_ca_;
	numeric::xyzVector<core::Real> bbl_ca_;
	numeric::xyzTriple< core::Size > cube_dim_ca_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::Real side_inv_ca_, neighbor_cutoff_ca_, neighbor_cutoff_sq_ca_;

};

} // namespace scoring {
} // namespace protocols {

#endif // INCLUDED_protocols_scoring_methods_ImplicitFastClashCheck_hh
