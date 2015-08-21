// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#ifndef INCLUDED_protocols_sic_dock_scores_TrisBpyScore_hh
#define INCLUDED_protocols_sic_dock_scores_TrisBpyScore_hh

#include <protocols/sic_dock/scores/TrisBpyScore.fwd.hh>
#include <protocols/sic_dock/RigidScore.hh>

#include <core/scoring/motif/motif_hash_stuff.fwd.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/sic_dock/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>
#include <core/pose/xyzStripeHashPose.fwd.hh>


namespace protocols {
namespace sic_dock {
namespace scores {

using core::pose::Pose;

class TrisBpyScore : public RigidScore {
public:
	TrisBpyScore(Pose const & pose, Real const & tolerance=0.3, Real const & min_bpy_contacts_=15);
	virtual ~TrisBpyScore();
	core::Real score      ( Xforms const & x1s, Xforms const & x2s ) const;
	core::Real score_extra( Xforms const & x1s, Xforms const & x2s, Real&cbc,Real&err,int&chirl ,Xform&xbpy) const;
	void      dump_trisbpy( Xforms const & x1s, Xforms const & x2s, std::string const & fname ) const;
	void show(std::ostream & out                                      , int width=10) const;
	void show(std::ostream & out, Xforms const & x1s, Xforms const & x2s, int width=10) const;
	std::string type() const { return "TrisBpy"; }
private:
	Pose const & pose_;
	numeric::Xforms bbx_;
	core::pose::xyzStripeHashPoseCOP cc_;
	Vecs cb_;
	Real const tolerance_;
	Real const min_bpy_contacts_;
};


} // namespace scores
} // namespace sic_dock
} // namespace protocols

#endif
