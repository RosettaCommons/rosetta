// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_sic_dock_sic_dock_hh
#define INCLUDED_protocols_sic_dock_sic_dock_hh


#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID_Map.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/sic_dock/xyzStripeHashPose.fwd.hh>


namespace protocols {
namespace sic_dock {

	class SICFast {
	public:
		typedef numeric::xyzVector<core::Real> Vec;

		SICFast();

		virtual ~SICFast();

		void init(
							core::pose::Pose const & pose1,
							core::pose::Pose const & pose2,
							core::id::AtomID_Map<core::Real> const & clash_atoms1,
							core::id::AtomID_Map<core::Real> const & clash_atoms2,
							core::id::AtomID_Map<core::Real> const & score_atoms1,
							core::id::AtomID_Map<core::Real> const & score_atoms2
							);

		// void init(
		// 					core::pose::Pose         const & cmp1in,
		// 					utility::vector1<Vec>            cmp1cbs,
		// 					utility::vector1<double> const & cmp1wts,
		// 					core::pose::Pose         const & cmp2in,
		// 					utility::vector1<Vec>            cmp2cbs,
		// 					utility::vector1<double> const & cmp2wts
		// 					);

		double slide_into_contact(core::kinematics::Stub const & xa,
															core::kinematics::Stub const & xb,
															utility::vector1<Vec>          pa,
															utility::vector1<Vec>          pb,
															utility::vector1<Vec>  const & cba,
															utility::vector1<Vec>  const & cbb,
															Vec                            ori,
															double                       & score
															);


	private:
		double refine_mindis_with_xyzHash(core::kinematics::Stub const & xa, core::kinematics::Stub const & xb,
																			utility::vector1<Vec> const & pa, utility::vector1<Vec> const & pb, double mindis, Vec ori);

		double get_score(core::kinematics::Stub const & xa, core::kinematics::Stub const & xb,
										 utility::vector1<Vec> const & cba, utility::vector1<Vec> const & cbb, Vec ori, double mindis);

		void rotate_points(utility::vector1<Vec> & pa, utility::vector1<Vec> & pb, Vec ori);
		void get_bounds(utility::vector1<Vec> & pa, utility::vector1<Vec> & pb);
		void fill_plane_hash(utility::vector1<Vec> & pa, utility::vector1<Vec> & pb);
		double get_mindis_with_plane_hashes();

	private:
		double xmx1,xmn1,ymx1,ymn1,xmx,xmn,ymx,ymn;
		double CTD,CLD,CTD2,CLD2,BIN;
		int xlb,ylb,xub,yub;
		xyzStripeHashPose *xh2_bb_,*xh2_cb_;
		utility::vector1<Vec>    clash1_,clash2_;
		utility::vector1<Vec>    score1_,score2_;
		utility::vector1<double> w1_,w2_;
		ObjexxFCL::FArray2D<Vec> ha,hb;

	};




	int flood_fill3D(int i, int j, int k, ObjexxFCL::FArray3D<double> & grid, double t);




	// core::pose::Pose ntrim(core::pose::Pose const & pose, int Nsym) {
	// 	core::pose::Pose trimmed(pose);
	// 	for(int i=1; i <= Nsym; ++i) {
	// 		trimmed.conformation().delete_residue_slow(1); // This is probably not quite right.
	// 	}
	// 	return trimmed;
	// }

	// void trim_tails(core::pose::Pose & pose, int Nsym) {
	// 	// for(int i = 1; i < )
	// }

	void termini_exposed(core::pose::Pose const & pose, bool & ntgood, bool & ctgood );


} // namespace sic_dock
} // namespace protocols

#endif
