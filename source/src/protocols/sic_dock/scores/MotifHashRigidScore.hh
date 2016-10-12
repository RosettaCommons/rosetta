// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#ifndef INCLUDED_protocols_sic_dock_scores_MotifHashRigidScore_hh
#define INCLUDED_protocols_sic_dock_scores_MotifHashRigidScore_hh

#include <protocols/sic_dock/scores/MotifHashRigidScore.fwd.hh>
#include <protocols/sic_dock/RigidScore.hh>

#include <core/scoring/motif/motif_hash_stuff.fwd.hh>
#include <core/pose/xyzStripeHashPose.fwd.hh>

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

class MotifHashRigidScore : public RigidScore {
	typedef std::map<std::string,Real> Stats;
	typedef std::pair<int,int> intint;
	typedef std::pair<numeric::xyzVector<float>,int> VecIR;
public:
	MotifHashRigidScore(core::pose::Pose const & pose1, core::pose::Pose const & pose2);
	virtual ~MotifHashRigidScore();
	core::Real score     ( Xforms const & x1, Xforms const & x2 ) const;
	core::Real score_meta( Xforms const & x1, Xforms const & x2, int & nsheet, Real & rawscore, Real & sselem_score , Real & coverage , Real & res_score, Real & sheetsc, int & nres, int &Nhh, int &Nee, int &Neh, int &Nh, int &Ne, int &Nl ) const;

	int dump_matching_motifs( core::pose::Pose   const & pose1, core::pose::Pose   const & pose2, std::ostream & out, int & count, core::pose::xyzStripeHashPoseCOP clash_check=NULL, bool print=false ) const;
	int dump_matching_motifs( Xforms const & x1s  , Xforms const & x2s  , std::ostream & out, core::pose::xyzStripeHashPoseCOP clash_check=NULL, bool print=false ) const;

	std::string type() const { return "MotifHash"; }

	void show(std::ostream & out                                      , int width=10) const;
	void show(std::ostream & out, Xforms const & x1, Xforms const & x2, int width=10) const;
	void show(std::ostream & out, Xforms const & x1, Xforms const & x2, core::pose::xyzStripeHashPoseCOP clash_check, int width=10) const;

	core::scoring::motif::MotifHashCOP motif_hash() const { return mh_; }

	core::Size nhashlookups() const { return nhashlookups_; }

private:
	MotifHashRigidScore(){}

	core::pose::Pose pose1_,pose2_;
	numeric::Xforms bbx1_,bbx2_;
	mutable core::scoring::motif::MotifHashCOP mh_;
	core::scoring::motif::XformScoreCOP xs_, xsee_, xseh_, xshe_, xshh_, xspp_;
	protocols::fldsgn::topology::SS_Info2 *ssinfo1_, *ssinfo2_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::Size nss1_,nss2_;
	core::pose::xyzStripeHashPose* reshash_;
	utility::vector1< std::pair<numeric::xyzVector<float>,int> > reslist_;
	bool hash_pose1_;
	std::string ss1_,ss2_;
	mutable core::Size nhashlookups_;
};

} // namespace scores
} // namespace sic_dock
} // namespace protocols

#endif
