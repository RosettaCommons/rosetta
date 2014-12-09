// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// includes
#ifndef INCLUDED_core_scoring_motif_util_hh
#define INCLUDED_core_scoring_motif_util_hh

#include <core/scoring/motif/motif_hash_stuff.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>
#include <core/pose/xyzStripeHashPose.fwd.hh>
#include <numeric/xyzVector.hh>
#include <utility/fixedsizearray1.hh>
#include <numeric/geometry/hashing/SixDHasher.hh>
#include <numeric/xyzTransform.hh>
#include <numeric/HomogeneousTransform.hh>
#include <boost/unordered_map.hpp>

namespace core {
namespace scoring {
namespace motif {

//types

	using core::id::AtomID;
	using core::pose::Pose;
	using core::pose::PoseCOP;
	using core::pose::PoseCAP;
	using core::Real;
	using core::scoring::ScoreFunctionOP;
	using core::Size;
	using std::string;
	using utility::vector1;
	using numeric::geometry::hashing::Real3;
	using numeric::geometry::hashing::Real6;
	using core::pose::xyzStripeHashPose;
	using core::pose::xyzStripeHashPoseCAP;

	typedef utility::vector1<Real> Reals;
	typedef utility::vector1<Size> Sizes;
	typedef utility::vector1<int>  Ints;
	typedef utility::vector1<float> Floats;
	typedef utility::vector1<bool> Bools;
	typedef numeric::xyzVector<core::Real> Vec;
	typedef numeric::xyzMatrix<core::Real> Mat;
	typedef numeric::xyzTransform<core::Real> Xform;

	std::ostream & operator<<(std::ostream & out, Real6 const & r6);

	RM_Type rpm_type1(RPM_Type const & type);
	RM_Type rpm_type2(RPM_Type const & type);

	Reals get_sasa(Pose const & pose, Real const & probesize);
	Reals get_nbrs(Pose const & pose);

	// Real get_sc_sasa(Pose const & pose, Size const & ir);
	// Real get_sc_sasa(Pose const & pose, Size const & ir, Size const & jr);
	// Real get_sc_sasa_2rsd(Pose const & pose, Size ir, Size jr);

	void xform_pose( core::pose::Pose & pose, Xform const & s, Size sres=1, Size eres=0 );
	Mat random_rotation();

	std::string tag_from_pdb_fname(string const & fname0);

	// Real rt6_rt6_bb_dis2_explicit_stupid(
	// 	Real6 const & x1,
	// 	Real6 const & x2
	// );

	// Real rt6_rt6_dis2(Real6 const & x1, Real6 const & x2, Real const & lever);
	// Real rt6_rt6_bb_dis2(Real6 const & x1, Real6 const & x2);
	Real6 inverse_rt6(Real6 const & rt);
	// Real6 rt_to_real6(core::kinematics::RT const & rt);
	// core::kinematics::RT real6_to_rt(Real6 const & rt6);

	Xform get_residue_pair_xform(Pose const & pose1, Size ir, Pose const & pose2, Size jr, RPM_Type const & type=BB_BB);
	Real6 get_residue_pair_rt6  (Pose const & pose1, Size ir, Pose const & pose2, Size jr, RPM_Type const & type=BB_BB);
	Xform get_residue_pair_xform(Pose const & pose , Size ir, Size jr, RPM_Type const & type=BB_BB);
	Real6 get_residue_pair_rt6  (Pose const & pose , Size ir, Size jr, RPM_Type const & type=BB_BB);

	void set_residue_pair_xform(Xform const & x, Pose & pose , Size ir, Size jr, RPM_Type const & type=BB_BB);

	Real align_motif_pose              ( Pose & pose, Pose const & paln1, Size const & ir, Pose const & paln2, Size const & jr, RPM_Type const & type );
	// Real align_motif_pose_NCAC_super   ( Pose & pose, Pose const & paln1, Size const & ir, Pose const & paln2, Size const & jr, RPM_Type const & type );
	Real align_motif_pose_break        ( Pose & pose, Pose const & paln1, Size const & ir, Pose const & paln2, Size const & jr, RPM_Type const & type );
	Real align_motif_pose_by_one_frame ( Pose & pose, Pose const & paln1, Size const & ir, Pose const & paln2, Size const & jr, RPM_Type const & type );
	Real align_motif_pose_super        ( Pose & pose, Pose const & paln1, Size const & ir, Pose const & paln2, Size const & jr, RPM_Type const & type );

	Real6 get_bins(Real c, Real a);

	core::id::AtomID_Mask get_motif_atom_mask( Pose const & motif_pose, RPM_Type const & type, bool with_Hpol=false );

	void HACK_dump_helix(Pose const & pose, string fname, int beg, int end);
	int HACK_dump_helices(Pose const & pose, string tag, int nres, int minlen=10);

	inline Real angle_distance(Real const & a, Real const & b){
		return numeric::min( fabs(a-b), fabs(fabs(a-b)-360.0) );
	}


	core::Real aa_trustworthiness(char aa);


}
}
}

#endif


