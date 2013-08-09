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
#include <protocols/rbsegment_relax/RBSegmentMover.hh>
// AUTO-REMOVED #include <protocols/rbsegment_relax/util.hh>

// Rosetta Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>

// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>

// AUTO-REMOVED #include <core/scoring/constraints/Constraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/BoundConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/CoordinateConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AmbiguousConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <basic/basic.hh>
#include <basic/Tracer.hh>

// Random number generator
// AUTO-REMOVED #include <numeric/xyzVector.io.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>

//
#include <string>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace rbsegment_relax {

using namespace core;

static numeric::random::RandomGenerator rbseg_RG(18632);
static basic::Tracer TR("protocols::moves::RBSegmentMover");

std::string
RBSegmentMover::get_name() const {
	return "RBSegmentMover";
}

void RBSegmentMover::set_movement( core::Real , core::Real , core::Real , core::Real ) { }


/// @brief Helper function to get a coordinate transformation from 3 points:
///   the origin, a point specifying the +z axis, and a point specifying the x-z plane
numeric::xyzMatrix< Real > RBSegmentMover::coordTransformFromThreePoints(
                                     numeric::xyzVector< Real > ori,
                                     numeric::xyzVector< Real > zAxis,
                                     numeric::xyzVector< Real > xzPlane ) {
	numeric::xyzVector< Real > xDir, yDir, zDir;
	numeric::xyzMatrix< Real > trans;

	zDir = zAxis - ori;
	zDir.normalize();

	xDir = xzPlane - ori;
	xDir.project_normal(zDir);
	xDir.normalize();

	yDir = zDir.cross(xDir);

	trans.col_x(xDir);
	trans.col_y(yDir);
	trans.col_z(zDir);

	return trans;
}

numeric::xyzVector< Real > RBSegmentMover::getCoM( core::pose::Pose const & pose ) {
	numeric::xyzVector< Real > CoM( 0, 0, 0 );
	int N = 0;

	for (Size i=1; i<=segment_.nContinuousSegments(); ++i) {
		// apply to transformation to every atom in this segment
		Size i_start = std::max(segment_[i].start(), (Size)1);
		Size i_end   = std::min(segment_[i].end(),  pose.total_residue());

		for ( Size j = i_start; j <= i_end; ++j ) {
			CoM += pose.residue(j).xyz("CA");
			N++;
		}
	}
	CoM = CoM * (1.0/N);

	return CoM;
}



/// @brief Apply an arbitrary rotation specified by a rotation matrix
void  RBSegmentMover::applyRotation( core::pose::Pose & pose, numeric::xyzMatrix< Real > rotation ) {
	//	get transformation
	numeric::xyzVector< Real > origin;
	numeric::xyzMatrix< Real > local2global, global2local;

	getCoordinateTransformation( pose, origin, local2global);
	global2local = numeric::inverse( local2global );

	utility::vector1< id::AtomID > atm_ids;
	utility::vector1< numeric::xyzVector< core::Real> > atm_xyzs;

	numeric::xyzVector< Real > localX, localRX, globalRX;
	for (Size i=1; i<=segment_.nContinuousSegments(); ++i) {
		// apply to transformation to every atom in this segment
		Size i_start = std::max(segment_[i].start(), (Size)1);
		Size i_end   = std::min(segment_[i].end() , pose.total_residue());

		for ( Size j = i_start; j <= i_end; ++j ) {
			for ( Size k = 1; k <= pose.residue(j).natoms(); ++k ) {
				id::AtomID id( k, j );
				localX = global2local * ( pose.xyz(id) - origin );
				localRX = rotation * localX;
				globalRX = local2global * localRX + origin;

				//pose.set_xyz( id, globalRX );
				atm_ids.push_back( id );
				atm_xyzs.push_back( globalRX );
			}
		}
	}
	pose.batch_set_xyz( atm_ids, atm_xyzs );
}

/// @brief Apply an arbitrary translation
void  RBSegmentMover::applyTranslation( core::pose::Pose & pose, numeric::xyzVector< Real > translation ) {
	//	get transformation
	numeric::xyzVector< Real > origin;
	numeric::xyzMatrix< Real > local2global;

	getCoordinateTransformation( pose, origin, local2global);

	utility::vector1< id::AtomID > atm_ids;
	utility::vector1< numeric::xyzVector< core::Real> > atm_xyzs;

	numeric::xyzVector< Real > localX, localRX, globalRX;
	for (Size i=1; i<=segment_.nContinuousSegments(); ++i) {
		// apply to transformation to every atom in this segment
		Size i_start = std::max(segment_[i].start(), (Size)1);
		Size i_end   = std::min(segment_[i].end(), pose.total_residue());

		for ( Size j = i_start; j <= i_end; ++j ) {
			for ( Size k = 1; k <= pose.residue(j).natoms(); ++k ) {
				id::AtomID id( k, j );
				globalRX = pose.xyz(id) + local2global * translation;

				//pose.set_xyz( id, globalRX );
				atm_ids.push_back( id );
				atm_xyzs.push_back( globalRX );
			}
		}
	}
	pose.batch_set_xyz( atm_ids, atm_xyzs );
}

/// @brief Apply a rotation followed by a translation (does not recompute coordinate transformation between the two!)
void  RBSegmentMover::applyTransformation( core::pose::Pose & pose, numeric::xyzMatrix< Real > rotation , numeric::xyzVector< Real > translation ) {
	//	get transformation
	numeric::xyzVector< Real > origin;
	numeric::xyzMatrix< Real > local2global, global2local;

	getCoordinateTransformation( pose, origin, local2global);
	global2local = numeric::inverse( local2global );

	utility::vector1< id::AtomID > atm_ids;
	utility::vector1< numeric::xyzVector< core::Real> > atm_xyzs;

	numeric::xyzVector< Real > localX, localRX, globalRX;
	for (Size i=1; i<=segment_.nContinuousSegments(); ++i) {
		// apply to transformation to every atom in this segment
		Size i_start = std::max(segment_[i].start(), (Size)1);
		Size i_end   = std::min(segment_[i].end(), pose.total_residue());

		for ( Size j = i_start; j <= i_end; ++j ) {
			for ( Size k = 1; k <= pose.residue(j).natoms(); ++k ) {
				id::AtomID id( k, j );
				localX = global2local * ( pose.xyz(id) - origin );
				localRX = rotation * localX + translation;
				globalRX = local2global * localRX + origin;

				//pose.set_xyz( id, globalRX );
				atm_ids.push_back( id );
				atm_xyzs.push_back( globalRX );
			}
		}
	}
	pose.batch_set_xyz( atm_ids, atm_xyzs );
}

/// @brief Apply a spin of the specified angle (in degrees) about arbitrary axis
void  RBSegmentMover::applySpin( core::pose::Pose & pose, numeric::xyzVector< Real > rotationAxis, Real degrees ) {
	numeric::xyzMatrix< Real > rotation =  numeric::rotation_matrix(rotationAxis, numeric::conversions::radians(degrees));
	applyRotation( pose, rotation );
}

/// @brief (re)set the starting and ending residues of this transform
void RBSegmentMover::setResidueRange( RBSegment const & seg ) {
	segment_ = seg;
}

/// @brief Get the the starting and ending residues of transform
RBSegment const & RBSegmentMover::getResidueRange() {
	return segment_;
}


///////////////////////////////////////////
///  Apply a random rigid-body transformation to an arbitrary segment
///////////////////////////////////////////
void  GaussianRBSegmentMover::apply( core::pose::Pose & pose ) {
	// random rotation ...
	applyRotation( pose , sigma_rot*rbseg_RG.gaussian() , sigma_rot*rbseg_RG.gaussian() , sigma_rot*rbseg_RG.gaussian() );

	// ... and translation
	numeric::xyzVector< Real > trans(sigma_trans*rbseg_RG.gaussian() , sigma_trans*rbseg_RG.gaussian() , sigma_trans*rbseg_RG.gaussian());
	applyTranslation( pose , trans );
}

std::string
GaussianRBSegmentMover::get_name() const {
	return "GaussianRBSegmentMover";
}

/// @brief Get the fragment center / matrix that rotates global coordinates into local fragment coordinates.
/// Defined such that +z points to the C-terminal end of the helix axis,
/// +x from the helix axis to the N-terminal residue
void GaussianRBSegmentMover::getCoordinateTransformation(
								core::pose::Pose const & pose,
								Vector &rotationCenter,
								numeric::xyzMatrix< Real > &coordinateTransform
																) {
	// rotate about center-of-mass
	rotationCenter = getCoM( pose );

	// random rotation, so local coord frame irrelevant
	coordinateTransform.clear();
	coordinateTransform.xx(1);
	coordinateTransform.yy(1);
	coordinateTransform.zz(1);
}



///////////////////////////////////////////
///  Apply a "register shift" move to this fragment.
///////////////////////////////////////////
void SequenceShiftMover::apply( core::pose::Pose & pose ) {
	// pick a direction at random
	int dir = (rbseg_RG.random_range(0,1))? -1 : 1;
	int mag = rbseg_RG.random_range(1,magnitude_);
	core::Size nres = hybridization::get_num_residues_nonvirt( pose );

	if (offsets_.size() != nres) {
		offsets_.clear();
		offsets_.resize(nres,0);
	}

	TR.Debug << "SequenceShiftMover::apply() [" << dir*mag << "]" << std::endl;

	utility::vector1< id::AtomID > atm_ids;
	utility::vector1< numeric::xyzVector< core::Real> > atm_xyzs;

	for (Size ii=1; ii<=segment_.nContinuousSegments(); ++ii) {
		// apply to transformation to every atom in this segment
		int i_start = std::max(segment_[ii].start(), (Size)1);
		int i_end   = std::min(segment_[ii].end(), pose.total_residue());

		int nres_seg  = i_end - i_start + 1;

		numeric::xyzVector< Real > C1,N1,C2,N2, CA1, CA2;
		numeric::xyzMatrix< Real > R, R1;

		// avoid compiler warning
		R.xx(0.0);R.xy(0.0);R.xz(0.0);
		R.yx(0.0);R.yy(0.0);R.yz(0.0);
		R.zx(0.0);R.zy(0.0);R.zz(0.0);
		R1.xx(0.0);R1.xy(0.0);R1.xz(0.0);
		R1.yx(0.0);R1.yy(0.0);R1.yz(0.0);
		R1.zx(0.0);R1.zy(0.0);R1.zz(0.0);

		for ( int i = 0; i < nres_seg-mag; ++i ) {
			// "transform" r_i to r_j
			Size r_i = (dir==1)? i_start+i : i_end-i;
			Size r_j = r_i+dir*mag;

			CA1 = pose.residue(r_i).atom("CA").xyz();
			C1 = pose.residue(r_i).atom("C").xyz() - CA1;  // offset from CA
			N1 = pose.residue(r_i).atom("N").xyz() - CA1;  // offset from CA
			CA2 = pose.residue(r_j).atom("CA").xyz();
			C2 = pose.residue(r_j).atom("C").xyz() - CA2;  // offset from CA
			N2 = pose.residue(r_j).atom("N").xyz() - CA2;  // offset from CA

			// get rotation from (i+dir) to i
			R = numeric::alignVectorSets( C1,N1, C2,N2 );

			// apply transformation to every atom in this res
			for ( Size a_i = 1; a_i<=pose.residue(r_i).natoms(); ++a_i ) {
				id::AtomID id( a_i, r_i );
				numeric::xyzVector< Real > newX = R * (pose.xyz(id) - CA1) + CA2;

				//pose.set_xyz( id, newX );
				atm_ids.push_back( id );
				atm_xyzs.push_back( newX );
			}
		}

		// final 'mag+1' residues
		//    if there is a reference pose, use that
		//    otherwise, duplicate the last residues' xform
		for ( int i = 0; i < mag; ++i ) {
			int r_i = (dir==1)? i_end-i : i_start+i;
			int r_j = r_i + dir*mag + offsets_[r_i];

			bool crosses_cut = (ref_pose_.total_residue() == 0);
			for (int k = 0; k<=i && !crosses_cut; ++k) {
				Size r_k = (dir==1)? i_end-k : i_start+k;
				Size r_l = r_k + dir*mag + offsets_[r_k];
				if (r_l>nres || r_l < 1) continue;
				crosses_cut |= ref_pose_.fold_tree().is_cutpoint((dir==1)?r_l:r_l+1);
				//if (crosses_cut) std::cerr << "CROSS CUT: " << i << " " << r_k << " -> " << r_l << std::endl;
			}

			if ( r_j<=nres && r_j >= 1 && !crosses_cut) {
				CA1 = pose.residue(r_i).atom("CA").xyz();
				C1 = pose.residue(r_i).atom("C").xyz() - CA1;  // offset from CA
				N1 = pose.residue(r_i).atom("N").xyz() - CA1;  // offset from CA
				CA2 = ref_pose_.residue(r_j).atom("CA").xyz();
				C2 = ref_pose_.residue(r_j).atom("C").xyz() - CA2;  // offset from CA
				N2 = ref_pose_.residue(r_j).atom("N").xyz() - CA2;  // offset from CA

				// get rotation from (i+dir) to i
				R1 = numeric::alignVectorSets( C1,N1, C2,N2 );

				// apply transformation to every atom in this res
				for ( Size a_i = 1; a_i<=pose.residue(r_i).natoms(); ++a_i ) {
					id::AtomID id( a_i, r_i );
					numeric::xyzVector< Real > newX = R1 * (pose.xyz(id) - CA1) + CA2;

					//pose.set_xyz( id, newX );
					atm_ids.push_back( id );
					atm_xyzs.push_back( newX );
				}
			} else {
				for ( Size a_i = 1; a_i<=pose.residue(r_i).natoms(); ++a_i ) {
					id::AtomID id( a_i, r_i );
					numeric::xyzVector< Real > newX = R * (pose.xyz(id) - CA1) + CA2;
					atm_ids.push_back( id );
					atm_xyzs.push_back( newX );
				}
			}
		}
	}

	pose.batch_set_xyz( atm_ids, atm_xyzs );

	last_move_ = dir*mag;
}

void SequenceShiftMover::trigger_accept() {
	// apply shift to array and score
	for (int j=1; j<=segment_.nContinuousSegments(); ++j) {
		for (int k=(int)segment_[j].start(); k<=(int)segment_[j].end(); ++k) {
			offsets_[k] += last_move_;
		}
	}
}

int
SequenceShiftMover::score() {
	int score_aln=0;

	// # block shifts
	for (int j=2; j<=(int)offsets_.size(); ++j)
		score_aln += abs(offsets_[j] - offsets_[j-1]);

	// penalize "unaligned" residues
	utility::vector1<bool> residues_to_rebuild (offsets_.size(), false);
	for (int i=1; i<=offsets_.size(); ++i) {

		if (residues_to_rebuild[i]) continue;

		// case 1: alignment passes N- or C- term
		if (offsets_[i]+i<1 || offsets_[i]+i>offsets_.size()) {
			residues_to_rebuild[i]=true;
		}

		// case 2: residues overlap
		for (int j=i+1; j<=offsets_.size(); ++j) {
			if (offsets_[i]+i == offsets_[j]+j) {
				residues_to_rebuild[i]=residues_to_rebuild[j]=true;
			}
		}

		// case 3: alignment crosses cutpoint
		int dir = (offsets_[i]>0) ? 1:-1;
		for (int j=i+dir; j!=i+offsets_[i]+dir && !residues_to_rebuild[i]; j+=dir) {
			if (ref_pose_.fold_tree().is_cutpoint((dir==1)?j:j+1)) {
				residues_to_rebuild[i]=true;
			}
		}
	}
	for (int i=1; i<=offsets_.size(); ++i) {
		if (residues_to_rebuild[i]) score_aln++;
	}

	return score_aln;
}


loops::LoopsOP
SequenceShiftMover::get_residues_to_rebuild() {
	loops::LoopsOP retval = new loops::Loops;
	if (ref_pose_.total_residue() == 0)
		return retval;

	utility::vector1<bool> residues_to_rebuild (offsets_.size(), false);
	for (int i=1; i<=offsets_.size(); ++i) {

		if (residues_to_rebuild[i]) continue;

		// case 1: alignment passes N- or C- term
		if (offsets_[i]+i<1 || offsets_[i]+i>offsets_.size()) {
			residues_to_rebuild[i]=true;
		}

		// case 2: residues overlap
		for (int j=i+1; j<=offsets_.size(); ++j) {
			if (offsets_[i]+i == offsets_[j]+j) {
				residues_to_rebuild[i]=residues_to_rebuild[j]=true;
			}
		}

		// case 3: alignment crosses cutpoint
		int dir = (offsets_[i]>0) ? 1:-1;
		for (int j=i+dir; j!=i+offsets_[i]+dir && !residues_to_rebuild[i]; j+=dir) {
			if (ref_pose_.fold_tree().is_cutpoint((dir==1)?j:j-1)) {
				residues_to_rebuild[i]=true;
			}
		}
	}

	bool inloop=false;
	int loopstart,loopstop;
	for (int i=1; i<=residues_to_rebuild.size(); ++i) {
		if (!inloop && residues_to_rebuild[i]) {
			inloop=true;
			loopstart=i;
		}
		if (inloop && !residues_to_rebuild[i]) {
			inloop=false;
			loopstop = i-1;
			if (!ref_pose_.fold_tree().is_cutpoint(loopstart-1)) loopstart--;
			if (!ref_pose_.fold_tree().is_cutpoint(loopstop)) loopstop++;
			retval->add_loop( loopstart,loopstop, (int)std::floor(0.5*(loopstart+loopstop+1)) );
		}
	}
	if (inloop) {
		if (!ref_pose_.fold_tree().is_cutpoint(loopstart-1)) loopstart--;
		retval->add_loop( loopstart,residues_to_rebuild.size(),
		                  (int)std::floor(0.5*(loopstart+residues_to_rebuild.size()+1)) );
	}
	return retval;
}


std::string
SequenceShiftMover::get_name() const {
	return "SequenceShiftMover";
}

////////////////////////////////////////////
///   Helical-axis segment movers
////////////////////////////////////////////
void HelicalGaussianMover::apply( core::pose::Pose & pose )
{
	// if the segment is not simple (i.e. contains >1 continuous segment) output error msg
	Real displacement_Z( sigAxisT_*rbseg_RG.gaussian() );
	Real displacement_X( sigOffAxisT_*rbseg_RG.gaussian() );
	Real displacement_Y( sigOffAxisT_*rbseg_RG.gaussian() );

	Real displacement_alpha( sigAxisR_*rbseg_RG.gaussian() );
	Real displacement_beta ( sigOffAxisR_*rbseg_RG.gaussian() );
	Real displacement_gamma( sigOffAxisR_*rbseg_RG.gaussian() );

	TR.Debug << "HelicalGaussianMover::apply() ["
	         << displacement_X << "," << displacement_Y << "," << displacement_Z << ","
	         << displacement_alpha << "," << displacement_beta << "," << displacement_gamma << "]" << std::endl;
	TR.Debug << "                              ["
	         << sigAxisT_ << "," << sigOffAxisT_ << "," << sigAxisR_ << "," << sigOffAxisR_ << "]" << std::endl;

//	TR << "HelixAxisGaussianTransMover::apply()" << std::endl;
	if ( sigAxisT_ > 1e-6 || sigOffAxisT_ > 1e-6 ) {
//		TR << "  Translation w.r.t. Helical Axis [" << displacement_X << " , "
//																								<< displacement_Y << " , "
//																								<< displacement_Z << "]" << std::endl;
		numeric::xyzVector< Real > trans( displacement_X, displacement_Y, displacement_Z );
		applyTranslation( pose, trans );
	}

	if ( sigAxisR_ > 1e-6 || sigOffAxisR_ > 1e-6 ) {
//		TR << "  Rotation w.r.t. Helical Axis [" << displacement_alpha << " , "
//																						 << displacement_beta  << " , "
//																						 << displacement_gamma << "]" << std::endl;
		applyRotation( pose, displacement_alpha, displacement_beta, displacement_gamma );
	}
}

std::string
HelicalGaussianMover::get_name() const {
	return "HelicalGaussianMover";
}

void HelicalGaussianMover::getCoordinateTransformation(
								core::pose::Pose const & pose,
								Vector &rotationCenter,
								numeric::xyzMatrix< Real > &coordinateTransform
																)
{
	// if the segment is not simple (i.e. contains >1 continuous segment) output error msg
	if (!segment_.isSimple()) {
		TR << "[ ERROR ] HelicalGaussianMover::getCoordinateTransformationy() called on compound segment!" << std::endl;
		exit(1);
	}

	// the helix axis is defined based on the residues near the center of the helix
	// if the helix is longer than 7 residues compare midPt-3..midPt with midPt..midPt+3
	// otherwise compare start..start+3 with end-3..end
	// undefined if helix is less than 4 residues
	Vector helixAxisNterm(0,0,0), helixAxisCterm(0,0,0);
	int startRes = segment_[1].start();
	int endRes = segment_[1].end();
	int nres =  endRes - startRes + 1;

	if (nres < 4) {
		TR << "[WARNING]  Helical axis of helices less than four residues is not correctly computed" << std::endl;
		rotationCenter = getCoM( pose );

		coordinateTransform.clear();
		coordinateTransform.xx(1);
		coordinateTransform.yy(1);
		coordinateTransform.zz(1);
	} else {
		if (nres >= 7) {
			int midRes = startRes + (nres-1)/2;
			helixAxisNterm += ( 0.6/3.6 ) * pose.residue( midRes ).xyz( "CA" );
			helixAxisNterm += ( 1.0/3.6 ) * pose.residue( midRes - 1 ).xyz( "CA" );
			helixAxisNterm += ( 1.0/3.6 ) * pose.residue( midRes - 2 ).xyz( "CA" );
			helixAxisNterm += ( 1.0/3.6 ) * pose.residue( midRes - 3 ).xyz( "CA" );
			helixAxisCterm += ( 0.6/3.6 ) * pose.residue( midRes ).xyz( "CA" );
			helixAxisCterm += ( 1.0/3.6 ) * pose.residue( midRes + 1 ).xyz( "CA" );
			helixAxisCterm += ( 1.0/3.6 ) * pose.residue( midRes + 2 ).xyz( "CA" );
			helixAxisCterm += ( 1.0/3.6 ) * pose.residue( midRes + 3 ).xyz( "CA" );
		} else {
			helixAxisNterm += (1.0/3.6) * pose.residue(startRes).xyz("CA");
			helixAxisNterm += (1.0/3.6) * pose.residue(startRes+1).xyz("CA");
			helixAxisNterm += (1.0/3.6) * pose.residue(startRes+2).xyz("CA");
			helixAxisNterm += (0.6/3.6) * pose.residue(startRes+3).xyz("CA");
			helixAxisCterm += (1.0/3.6) * pose.residue(endRes).xyz("CA");
			helixAxisCterm += (1.0/3.6) * pose.residue(endRes-1).xyz("CA");
			helixAxisCterm += (1.0/3.6) * pose.residue(endRes-2).xyz("CA");
			helixAxisCterm += (0.6/3.6) * pose.residue(endRes-3).xyz("CA");
		}

		rotationCenter = 0.5 * helixAxisNterm + 0.5 * helixAxisCterm;
		coordinateTransform = coordTransformFromThreePoints(
																 rotationCenter,
																 helixAxisCterm,
																 pose.residue(startRes).xyz("CA") );
	}
}


////////////////////////
////////////////////////
void StrandTwistingMover::apply( core::pose::Pose & pose )
{
	Real displacement_Z( sigAxisT_*rbseg_RG.gaussian() );
	Real displacement_X( sigOffAxisT_*rbseg_RG.gaussian() );
	Real displacement_Y( sigOffAxisT_*rbseg_RG.gaussian() );

	Real displacement_alpha( sigAxisR_*rbseg_RG.gaussian() );
	Real displacement_beta ( sigOffAxisR_*rbseg_RG.gaussian() );
	Real displacement_gamma( sigOffAxisR_*rbseg_RG.gaussian() );

	TR.Debug << "StrandTwistingMover::apply() ["
	         << displacement_X << "," << displacement_Y << "," << displacement_Z << ","
	         << displacement_alpha << "," << displacement_beta << "," << displacement_gamma << "]" << std::endl;

	if ( sigAxisT_ > 1e-6 || sigOffAxisT_ > 1e-6 ) {
		numeric::xyzVector< Real > trans( displacement_X, displacement_Y, displacement_Z );
		applyTranslation( pose, trans );
	}

	if ( sigAxisR_ > 1e-6 || sigOffAxisR_ > 1e-6 ) {
		applyRotation( pose, displacement_alpha, displacement_beta, displacement_gamma );
	}
}

std::string
StrandTwistingMover::get_name() const {
	return "StrandTwistingMover";
}

void StrandTwistingMover::getCoordinateTransformation(
	core::pose::Pose const & pose,
	Vector &rotationCenter,
	numeric::xyzMatrix< Real > &coordinateTransform
)
{
	// if the segment is not simple (i.e. contains >1 continuous segment) output error msg
	if (!segment_.isSimple()) {
		TR << "[ ERROR ] HelicalGaussianMover::getCoordinateTransformationy() called on compound segment!" << std::endl;
		exit(1);
	}

	// undefined if strand is less than 3 residues
	Vector sAxisNterm(0,0,0), sAxisCterm(0,0,0);
	int startRes = segment_[1].start();
	int endRes = segment_[1].end();
	int nres =  endRes - startRes + 1;

	if (nres < 3) {
		TR << "[WARNING]  Strand axis of less than three residues is not correctly computed" << std::endl;
		rotationCenter = getCoM( pose );

		coordinateTransform.clear();
		coordinateTransform.xx(1);
		coordinateTransform.yy(1);
		coordinateTransform.zz(1);
	} else {
		sAxisNterm += (0.5) * pose.residue(startRes).xyz("CA");
		sAxisNterm += (0.5) * pose.residue(startRes+1).xyz("CA");
		sAxisCterm += (0.5) * pose.residue(endRes).xyz("CA");
		sAxisCterm += (0.5) * pose.residue(endRes-1).xyz("CA");

		rotationCenter = 0.5 * sAxisNterm + 0.5 * sAxisCterm;
		coordinateTransform = coordTransformFromThreePoints(
																 rotationCenter,
																 sAxisCterm,
																 pose.residue(startRes).xyz("CA") );
	}
}

}
}
