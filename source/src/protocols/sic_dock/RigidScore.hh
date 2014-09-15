
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spac, int width=10es:t -*-
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
#include <protocols/sic_dock/types.hh>
#include <protocols/sic_dock/xyzStripeHashPoseWithMeta.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

namespace protocols {
namespace sic_dock {


class RigidScore : public utility::pointer::ReferenceCount {
protected:
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef numeric::xyzVector<Real> Vec;
	typedef numeric::xyzMatrix<Real> Mat;
	typedef utility::vector1<Vec> Vecs;
	typedef utility::vector1<Real> Reals;
	typedef utility::vector1<Size> Sizes;
	typedef numeric::Xforms Xforms;
	typedef utility::vector1<RigidScoreCOP> Scores;
public:

	virtual ~RigidScore() {}

	virtual core::Real score( Xforms const & x1, Xforms const & x2 ) const = 0;

	virtual std::string type() const = 0;
	virtual void show(std::ostream & out                                      , int width=10) const;
	virtual void show(std::ostream & out, Xforms const & x1, Xforms const & x2, int width=10) const;
};


class CBScore : public RigidScore {
public:
	// Undefined, commenting out to fix PyRosetta build  CBScore(Real clash_dis, Real contact_dis);
	virtual ~CBScore(){}
	CBScore(
		Pose const & pose1,
		Pose const & pose2,
		Real clash_dis,
		Real contact_dis
	);
	CBScore(
		Pose const & pose1,
		Pose const & pose2,
		Real clash_dis,
		Real contact_dis,
		core::id::AtomID_Map<core::Real> const & weights1,
		core::id::AtomID_Map<core::Real> const & weights2
	);
	core::Real score( Xforms const & x1, Xforms const & x2 ) const;
	std::string type() const { return "CBScore"; }
	// Undefinded, commenting out to fix PyRosetta build  Real dist2_score( Real const & sqdist ) const;
//private:
	bool const hash_pose1_;
	core::Real const clash_dis_, contact_dis_;
	Reals const weights_;
	Vecs const points_;
	xyzStripeHashPoseWithMeta const xyzhash_;
	// Pose const & pose1_,pose2_;
};

template<typename T>
inline T CBScore_dist_score( T const & sqdist, T const & start, T const & stop ) {
	if( sqdist > stop*stop ) {
		return (T)0.0;
	} else if( sqdist < start*start ) {
		return (T)1.0;
	} else {
		core::Real dist = sqrt( sqdist );
		return (stop-dist)/(stop-start);
		//return sqr(1.0	- sqr( (dist - start) / (stop - start) ) );
	}
}

class LinkerScore : public RigidScore {
public:
	LinkerScore(
		Pose const & pose1,
		Pose const & pose2,
		Size max_loop_len,
		Size lookup_radius,
		std::string const & outtag
	);
	virtual ~LinkerScore(){}
	core::Real  score( Xforms const & x1, Xforms const & x2 ) const;
	// Undefined, commenting out to fix PyRosetta build  void dump_linkers( Xform const & x1, Xform const & x2 ) const;
	bool dump_linkers( Xform const & x1, Xform const & x2, std::string const & out_perfix ) const;
	std::string type() const { return "LinkerScore"; }
private:
	static protocols::loophash::LoopHashLibraryOP loop_hash_library_;
	Sizes const loopsizes_;
	core::Size lookup_radius_;
	Pose const & pose1_,pose2_;
	TermInfo lowers1_,uppers1_,lowers2_,uppers2_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// Real max_dis2_;
	std::string outtag_;
};



struct AtomIDHashFunction { std::size_t operator ()(core::id::AtomID const aid) const {
	return aid.atomno()+65536*aid.rsd();
} };

// NOTE!!!!
// my convention here is that AtomIDs for pose2 are residue number 1000000 + actual
//
class ConstraintSetScore : public RigidScore {
public:
	ConstraintSetScore(
		Pose const & pose1,
		Pose const & pose2,
		core::scoring::constraints::ConstraintSet const & cstset
	);
	virtual ~ConstraintSetScore(){}
	core::Real score( Xforms const & x1, Xforms const & x2 ) const;
	std::string type() const { return "ConstraintSetScore"; }
private:
        // KAB - below variables commented out (-Wunused-private-field) on 2014-09-11
        // Pose const & pose1_;
	// Pose const & pose2_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::scoring::constraints::ConstraintSet const & cstset_;
	core::scoring::constraints::ConstraintCOPs const csts_;
	boost::unordered_map<core::id::AtomID,Vec,AtomIDHashFunction> start_coords_;
};





class EdgeStandScore : public RigidScore {
public:
	EdgeStandScore();
	virtual ~EdgeStandScore(){}
	// Undefined, commenting out to fix PyRosetta build  core::Real score( Xforms const & x1, Xforms const & x2 ) const;
private:
	Vecs donors,acceptors;
};

class HelixScore : public RigidScore {
public:
	HelixScore();
	virtual ~HelixScore(){}
	// Undefined, commenting out to fix PyRosetta build  core::Real score( Xforms const & x1, Xforms const & x2 ) const;
private:
};

class BuriedPolarScore : public RigidScore {
public:
	BuriedPolarScore(); // c'tor should store the buriend unsat polar coords
	virtual ~BuriedPolarScore(){}
	// Undefined, commenting out to fix PyRosetta build
	// core::Real score( Xforms const & x1, Xforms const & x2 ) const;
private:
	Vecs polars;
};

////// composite scores

class JointScore : public RigidScore {
public:
	JointScore();
	JointScore(
		Scores scores,
		Reals weights
	);
	void add_score(RigidScoreCOP score, Real weight);

	virtual ~JointScore(){}
	std::string type() const { return "JointScore"; }
	void show(std::ostream & out                                      , int width=10) const;
	void show(std::ostream & out, Xforms const & x1, Xforms const & x2, int width=10) const;
	core::Real score( Xforms const & x1, Xforms const & x2 ) const;
private:
	Scores scores_;
	Reals weights_;
  // KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
  // Real minscore_hack_,maxscore_hack_;

};


// class CachedScore : public RigidScore {
// public:
// 	CachedScore(RigidScoreCOP score);
// 	virtual ~CachedScore(){}
// 	core::Real score( Xforms const & x1, Xforms const & x2 ) const;
// private:
// 	RigidScoreCOP score_;
// 	// some kind of 6 dof hash
// };






} // namespace sic_dock
} // namespace protocols

#endif
