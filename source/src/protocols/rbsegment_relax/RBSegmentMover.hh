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
/// @author Frank DiMaio
/// @author Srivatsan Raman

#ifndef INCLUDED_protocols_rbsegment_relax_RBSegmentMover_hh
#define INCLUDED_protocols_rbsegment_relax_RBSegmentMover_hh

#include <protocols/moves/Mover.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

#include <protocols/rbsegment_relax/RBSegment.hh>
#include <protocols/hybridization/util.hh>

// C++ Headers
#include <map>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>


namespace protocols {
namespace rbsegment_relax {

//////////////////////////////////////////////////////////
///@brief Performs a rigid-body movement on a segment of a protein
///       Derived classes must implement 'getCoordinateTransform' and 'apply'
/////////////////////////////////////////////////////////
class RBSegmentMover : public protocols::moves::Mover {
protected:
	RBSegment segment_;

	/// @brief Helper function to get a coordinate transformation from 3 points:
	///   the origin, a point specifying the +z axis, and a point specifying the x-z plane
	static numeric::xyzMatrix< core::Real > coordTransformFromThreePoints(
	                                       numeric::xyzVector< core::Real > ori,
	                                       numeric::xyzVector< core::Real > zAxis,
	                                       numeric::xyzVector< core::Real > xzPlane );

	/// @brief Helper function computes center of mass of the segment in the pose
	inline numeric::xyzVector< core::Real > getCoM( core::pose::Pose const & pose );

public:
	/// @brief constructor
	RBSegmentMover() {}

	/// @brief constructor
	RBSegmentMover( RBSegment const & seg ) :
		segment_(seg) {	}

	/// @brief Apply the rigid-body fragment mover to a pose.  Must be defined by derived classes.
	virtual void apply( core::pose::Pose & pose ) = 0;

	/// @brief Set the segment this mover is working on
	virtual void set_segment( RBSegment const & seg) { segment_ = seg; }

	virtual std::string get_name() const;

	/// @brief Returns: (a) the matrix that rotates global coordinates into local fragment coordinates and
	/// (b) the center of rotation of the fragment.
	/// Default implementation sets them to global coords
	virtual void getCoordinateTransformation(
		core::pose::Pose const & /*pose*/ ,
		numeric::xyzVector< core::Real > &rotationCenter,
		numeric::xyzMatrix< core::Real > &coordinateTransform
	)
	{
		rotationCenter = numeric::xyzVector< core::Real >(0,0,0);
		coordinateTransform.clear();
		coordinateTransform.xx(1); coordinateTransform.yy(1); coordinateTransform.zz(1);
	}

	/// @brief Set mover-specific movement parameters.  Do nothing by default.
	virtual void set_movement( core::Real p1=0.0, core::Real p2=0.0, core::Real p3=0.0, core::Real p4=0.0);

	/// @brief steal movement params from a segment definition
	virtual void set_movement( RBSegment const & rb) {
		core::Real p1,p2,p3,p4;
		rb.get_movement( p1,p2,p3,p4 );
		set_movement( p1,p2,p3,p4 );
	}

	/// @brief Apply an arbitrary rotation specified by Euler angles (in degrees!)
	inline
	void applyRotation( core::pose::Pose & pose, core::Real alpha, core::Real beta, core::Real gamma ) {
		double sa = std::sin( numeric::conversions::radians( alpha ) );
		double ca = std::cos( numeric::conversions::radians( alpha ) );
		double sb = std::sin( numeric::conversions::radians( beta  ) );
		double cb = std::cos( numeric::conversions::radians( beta  ) );
		double sg = std::sin( numeric::conversions::radians( gamma ) );
		double cg = std::cos( numeric::conversions::radians( gamma ) );

		numeric::xyzMatrix< core::Real > rotation;
		rotation.xx( -sa*cb*sg + ca*cg ); rotation.xy(  ca*cb*sg + sa*cg ); rotation.xz(  sb*sg );
		rotation.yx( -sa*cb*cg - ca*sg ); rotation.yy(  ca*cb*cg - sa*sg ); rotation.yz(  sb*cg );
		rotation.zx(  sa*sb );            rotation.zy( -ca*sb );            rotation.zz(  cb );

		applyRotation( pose, rotation );
	}

	/// @brief Apply an arbitrary rotation specified by a rotation matrix
	void applyRotation( core::pose::Pose & pose, numeric::xyzMatrix< core::Real > rotation );

	/// @brief Apply an arbitrary translation
	void applyTranslation( core::pose::Pose & pose, numeric::xyzVector< core::Real > translation );

	/// @brief Apply a rotation followed by a translation (does not recompute coordinate transformation between the two!)
	void applyTransformation( core::pose::Pose & pose, numeric::xyzMatrix< core::Real > rotation , numeric::xyzVector< core::Real > translation );

	/// @brief Apply a spin of the specified angle (in degrees) about arbitrary axis
	void applySpin( core::pose::Pose & pose, numeric::xyzVector< core::Real > rotationAxis, core::Real degrees );

	/// @brief (re)set the starting and ending residues of this transform
	void setResidueRange( RBSegment const & seg );

	/// @brief Get the the starting and ending residues of transform
	RBSegment const & getResidueRange();

	/// @brief Print debugging info
	void print() {
		int nsegs = segment_.nContinuousSegments();
		std::cerr << nsegs << " segment object\n";
		for (int i=1; i<=nsegs; ++i) {
			std::cerr << "    " << i << ":  " << segment_[i].start() << " to " << segment_[i].end() << std::endl;
		}
	}
};


///////////////////////////////////////////
///
///  Random movements wrt the helical axis
///
///////////////////////////////////////////
class HelicalGaussianMover : public RBSegmentMover {
private:
	core::Real sigAxisT_;
	core::Real sigAxisR_;
	core::Real sigOffAxisT_;
	core::Real sigOffAxisR_;

public:
	/// @brief constructor:
	/// @li sigAxisR: the stdev of rotation along the helical axis
	/// @li sigAxisT: the stdev of movement along the helical axis
	/// @li sigOffAxisR: the stdev of rotation normal to the helical axis
	/// @li sigOffAxisT: the stdev of movement normal to the helical axis
	HelicalGaussianMover( RBSegment const & seg,
	                      core::Real sigAxisR=0.0,
	                      core::Real sigAxisT=0.0,
	                      core::Real sigOffAxisR=0.0,
	                      core::Real sigOffAxisT=0.0 ) :
		RBSegmentMover(seg),
		sigAxisT_( sigAxisT ),
		sigAxisR_( sigAxisR ),
		sigOffAxisT_( sigOffAxisT ),
		sigOffAxisR_( sigOffAxisR )
	{ }

	/// @brief constructor:
	/// @li sigAxisR: the stdev of rotation along the helical axis
	/// @li sigAxisT: the stdev of movement along the helical axis
	/// @li sigOffAxisR: the stdev of rotation normal to the helical axis
	/// @li sigOffAxisT: the stdev of movement normal to the helical axis
	HelicalGaussianMover( core::Real sigAxisR=0.0, core::Real sigAxisT=0.0, core::Real sigOffAxisR=0.0, core::Real sigOffAxisT=0.0 ) :
		sigAxisT_( sigAxisT ),
		sigAxisR_( sigAxisR ),
		sigOffAxisT_( sigOffAxisT ),
		sigOffAxisR_( sigOffAxisR )
	{ }

	/// @brief set movement parameters
	/// @li sigAxisR: the stdev of rotation along the helical axis
	/// @li sigAxisT: the stdev of movement along the helical axis
	/// @li sigOffAxisR: the stdev of rotation normal to the helical axis
	/// @li sigOffAxisT: the stdev of movement normal to the helical axis
	virtual void set_movement( core::Real sigAxisR=0.0, core::Real sigAxisT=0.0, core::Real sigOffAxisR=0.0, core::Real sigOffAxisT=0.0)
	{
		if (sigAxisT != 0.0) sigAxisT_ = sigAxisT;
		if (sigAxisR != 0.0) sigAxisR_ = sigAxisR;
		if (sigOffAxisT != 0.0) sigOffAxisT_ = sigOffAxisT;
		if (sigOffAxisR != 0.0) sigOffAxisR_ = sigOffAxisR;

		std::cerr << "HelicalGaussianMover : Setting params (" << sigAxisR_ << ", " << sigAxisT_ << ", " << sigOffAxisR_ << ", " << sigOffAxisT_ << ")\n";
	}

	virtual protocols::moves::MoverOP clone() const { return RBSegmentMoverOP(new HelicalGaussianMover(*this)); }

	/// @brief Apply a +1 or -1 residue "shift" to this helix
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	/// @brief Get the fragment center / matrix that rotates global coordinates into local fragment coordinates.
	/// Defined such that +z points to the C-terminal end of the helix axis,
	/// +x from the helix axis to the N-terminal residue
	void getCoordinateTransformation(  core::pose::Pose const & pose,
																		 numeric::xyzVector< core::Real > &rotationCenter,
																		 numeric::xyzMatrix< core::Real > &coordinateTransform );
};

///////////////////////////////////////////
///
///  "Register shift" a segment by one amino-acid
///  Works in both centroid and all-atom
///
///////////////////////////////////////////
class SequenceShiftMover : public RBSegmentMover {
public:
	/// @brief constructor
	SequenceShiftMover(RBSegment const & seg, core::Size magnitude=4) :
		RBSegmentMover(seg),
		magnitude_(magnitude),
		last_move_(0)
	{}

	/// @brief constructor
	SequenceShiftMover(core::pose::Pose const & pose_in, RBSegment const & seg, core::Size magnitude=4) :
		RBSegmentMover(seg),
		magnitude_(magnitude),
		last_move_(0),
		ref_pose_(pose_in)
	{
		core::Size nres = hybridization::get_num_residues_nonvirt( pose_in );
		offsets_.resize(nres,0);
	}

	/// @brief constructor
	SequenceShiftMover(RBResidueRange const &rb, core::Size magnitude=4) :
		RBSegmentMover(RBSegment(rb)),
		magnitude_(magnitude),
		last_move_(0)
	{}

	SequenceShiftMover() :
		RBSegmentMover(),
		magnitude_(4),
		last_move_(0)
	{}

	/// @brief clone this object
	virtual protocols::moves::MoverOP
	clone() const { return RBSegmentMoverOP(new SequenceShiftMover(*this)); }

	/// @brief set movement parameters.  ignore all input args
	virtual void
	set_movement( core::Real /*sigAxisR=0.0*/, core::Real /*sigAxisT=0.0*/, core::Real /*sigOffAxisR=0.0*/, core::Real /*sigOffAxisT=0.0*/)	{ }

	/// @brief Apply a + or - residue "shift" to this helix
	void
	apply( core::pose::Pose & pose );

	/// @brief Last move accepted; update offsets_
	void
	trigger_accept();


	/// @brief Return a score: the minimum number of block shifts accumulated
	int
	score();


	/// @brief IF ref_pose is given, get a list of residues to rebuild via fragment insertion
	loops::LoopsOP
	get_residues_to_rebuild();

private:
	core::Size magnitude_;
	int last_move_;

	core::pose::Pose ref_pose_;
	utility::vector1< int > offsets_;

	virtual std::string get_name() const;
};

///////////////////////////////////////////
///
///  Generic random segment mover
///
///////////////////////////////////////////
class GaussianRBSegmentMover : public RBSegmentMover {
private:
	core::Real sigma_rot, sigma_trans;

public:
	/// @brief constructor
	/// @li sigT: the stdev of movement along the helical axis
	/// @li sigR: the stdev of rotation along the helical axis
	GaussianRBSegmentMover(RBSegment const & seg, core::Real sigR = 0.0, core::Real sigT = 0.0) :
		RBSegmentMover(seg),
		sigma_rot(sigR),
		sigma_trans(sigT)
	{}

	GaussianRBSegmentMover(core::Real sigR = 0.0, core::Real sigT = 0.0) :
		sigma_rot(sigR),
		sigma_trans(sigT)
	{}

	/// @brief set movement parameters
	/// @li sigR: the stdev of rotation
	/// @li sigT: the stdev of movement
	virtual void set_movement( core::Real sigAxisR=0.0, core::Real sigAxisT=0.0, core::Real sigOffAxisR=0.0, core::Real sigOffAxisT=0.0)
	{
		if (sigAxisT != 0.0) sigma_trans = sigAxisT;
		if (sigAxisR != 0.0) sigma_rot = sigAxisR;

		std::cerr << "GaussianRBSegmentMover : Setting params (" << sigma_rot << ", " << sigma_trans << ")\n";
		if (sigOffAxisR != 0.0 || sigOffAxisT != 0.0 )
			std::cerr << "GaussianRBSegmentMover : Ignore params (" << sigOffAxisR << ", " << sigOffAxisT << ")\n";
	}

	/// @brief clone this object
	virtual protocols::moves::MoverOP clone() const { return RBSegmentMoverOP(new GaussianRBSegmentMover(*this)); }

	/// @brief Randomly perturb the segment
	void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	/// @brief Get the fragment center / matrix that rotates global coordinates into local fragment coordinates.
	/// Defined such that +z points to the C-terminal end of the helix axis,
	/// +x from the helix axis to the N-terminal residue
	void getCoordinateTransformation( core::pose::Pose const & pose,
																		numeric::xyzVector< core::Real > &rotationCenter,
																		numeric::xyzMatrix< core::Real > &coordinateTransform );
};

///////////////////////////////////////////
///
///  Strand Twisting
///
///////////////////////////////////////////
class StrandTwistingMover : public RBSegmentMover {
private:
	core::Real sigAxisT_;
	core::Real sigAxisR_;
	core::Real sigOffAxisT_;
	core::Real sigOffAxisR_;

public:
	/// @brief constructor:
	/// @li sigAxisR: the stdev of rotation along the helical axis
	/// @li sigAxisT: the stdev of movement along the helical axis
	/// @li sigOffAxisR: the stdev of rotation normal to the helical axis
	/// @li sigOffAxisT: the stdev of movement normal to the helical axis
	StrandTwistingMover ( RBSegment const & seg,
	                      core::Real sigAxisR=0.0,
	                      core::Real sigAxisT=0.0,
	                      core::Real sigOffAxisR=0.0,
	                      core::Real sigOffAxisT=0.0 ) :
		RBSegmentMover(seg),
		sigAxisT_( sigAxisT ),
		sigAxisR_( sigAxisR ),
		sigOffAxisT_( sigOffAxisT ),
		sigOffAxisR_( sigOffAxisR )
	{ }

	/// @brief constructor:
	/// @li sigAxisR: the stdev of rotation along the helical axis
	/// @li sigAxisT: the stdev of movement along the helical axis
	/// @li sigOffAxisR: the stdev of rotation normal to the helical axis
	/// @li sigOffAxisT: the stdev of movement normal to the helical axis
	StrandTwistingMover( core::Real sigAxisR=0.0, core::Real sigAxisT=0.0, core::Real sigOffAxisR=0.0, core::Real sigOffAxisT=0.0 ) :
		sigAxisT_( sigAxisT ),
		sigAxisR_( sigAxisR ),
		sigOffAxisT_( sigOffAxisT ),
		sigOffAxisR_( sigOffAxisR )
	{ }

	/// @brief set movement parameters
	/// @li sigAxisR: the stdev of rotation along the helical axis
	/// @li sigAxisT: the stdev of movement along the helical axis
	/// @li sigOffAxisR: the stdev of rotation normal to the helical axis
	/// @li sigOffAxisT: the stdev of movement normal to the helical axis
	virtual void set_movement( core::Real sigAxisR=0.0, core::Real sigAxisT=0.0, core::Real sigOffAxisR=0.0, core::Real sigOffAxisT=0.0)
	{
		if (sigAxisT != 0.0) sigAxisT_ = sigAxisT;
		if (sigAxisR != 0.0) sigAxisR_ = sigAxisR;
		if (sigOffAxisT != 0.0) sigOffAxisT_ = sigOffAxisT;
		if (sigOffAxisR != 0.0) sigOffAxisR_ = sigOffAxisR;

		std::cerr << "StrandTwistingMover : Setting params (" << sigAxisR_ << ", " << sigAxisT_ << ", " << sigOffAxisR_ << ", " << sigOffAxisT_ << ")\n";
	}

	virtual protocols::moves::MoverOP clone() const { return RBSegmentMoverOP(new StrandTwistingMover(*this)); }

	/// @brief Apply a +1 or -1 residue "shift" to this helix
	void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	/// @brief Get the fragment center / matrix that rotates global coordinates into local fragment coordinates.
	/// Defined such that +z points to the C-terminal end of the helix axis,
	/// +x from the helix axis to the N-terminal residue
	void getCoordinateTransformation(
		core::pose::Pose const & pose,
		numeric::xyzVector< core::Real > &rotationCenter,
		numeric::xyzMatrix< core::Real > &coordinateTransform
	);
};



}
}

#endif
