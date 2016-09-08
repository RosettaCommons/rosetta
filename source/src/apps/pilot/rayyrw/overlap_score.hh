#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh> // remove virtual residues

#include <core/import_pose/import_pose.hh>

#include <core/chemical/ChemicalManager.hh>

#include <utility/vector1.hh>
#include <utility/exit.hh>

#include <apps/pilot/rayyrw/rms_util.hh>
#include <apps/pilot/rayyrw/util.hh>

#include <basic/Tracer.hh>

#ifndef apps_pilot_rayyrw_overlap_score_hh
#define apps_pilot_rayyrw_overlap_score_hh

static basic::Tracer tr_overlap_score("overlap_score");

// read into 2 fragments pose and find the compatibility between them
// given a pair of residue-overlapped fragments,
// compat_score =
//  overlap_score (dist check for overlapped residues)
//  + clash_score (to prevent residues which are not overlapped intertwined with each other //  no clash is allowed! )

bool isOverlapping(
	core::pose::Pose const &pose_lower,
	core::Size const pos_lower,
	core::Size const pos_upper
	//core::pose::Pose const  &pose_upper,
){
	assert ( pos_lower < pos_upper );

	core::Size pose_lower_endrsn = pos_lower + pose_lower.size() - 1;
	if ( pose_lower_endrsn >= pos_upper ) { // "=" mean 1 overlap
		return true;
	} else {
		return false;
	}
}


// Problems found
// 1. some intermediate overlap just killed me - need to have an yes/no overlap.
//    we only want perfect structural overlap for overlapping residues
// 2. give penalty for residues where they are overlapping,
//    but do not have structural overlap
core::Real
soft_overlap_score(
	core::pose::Pose const &pose1,
	core::Size const pos1,
	core::pose::Pose const &pose2,
	core::Size const pos2,
	core::Real overlap_dist_cutoff=3.0,
	core::Real steepness_wt=8.0,
	// clash_score = when there is a clash after seq_sep 5, return 500 immediately
	int seq_sep=5,
	core::Real clash_dist=3.0,
	core::Real clash_return_score=10
){
	if ( pos1 == pos2 ) return 0;

	// because of the numbering matching issue - the pos2 must be further than pos1
	core::pose::Pose const &pose_lower = ( pos1 > pos2 ) ? pose2 : pose1;
	core::pose::Pose const &pose_upper = ( pos1 > pos2 ) ? pose1 : pose2;
	core::Size pos_lower = std::min( pos1, pos2 );
	core::Size pos_upper = std::max( pos1, pos2 );

	core::Size pose_lower_nrsds = pose_lower.size();
	core::Size pose_upper_nrsds = pose_upper.size();

	// check overlap or not
	//if ( ! isOverlapping( pose_lower, pos_lower, pose_upper, pos_upper ) )
	if ( ! isOverlapping( pose_lower, pos_lower, pos_upper ) ) {
		return 0;
	}

	core::Real overlap_score( 0.0 );
	core::Real sum_overlap_score( 0.0 );

	int offset = pos2 - pos1; // should therefore always be positive

	for ( core::Size irsd=1; irsd <= pose_lower_nrsds; ++irsd ) {
		numeric::xyzVector<core::Real> atm_i = pose_lower.residue( irsd ).atom(2).xyz();

		for ( core::Size jrsd=1; jrsd <= pose_upper_nrsds; ++jrsd ) {
			numeric::xyzVector<core::Real> atm_j = pose_upper.residue( jrsd ).atom(2).xyz();
			core::Real dist = ( atm_i - atm_j ).length();

			//tr_overlap_score << "irsd: " << pos1+irsd-1 << " jrsd: " << pos1+jrsd+offset-1 << " dist: " << dist ;
			if ( irsd == jrsd+offset ) {
				// at least having an insurance here saying this is crazy shouldn't happen
				if ( dist > 5.0 ) return clash_return_score;

				overlap_score = 2.0 / ( 1 + std::exp( -steepness_wt*(dist-overlap_dist_cutoff)) ) - 1;
				sum_overlap_score += overlap_score;

			} else { // residues that are overlapping
				// ij_seq_sep less than seq_sep don't do clash check
				int ij_seq_sep = std::abs( int( irsd - (jrsd+offset) ) );
				if ( ij_seq_sep <= seq_sep ) {
					//tr_overlap_score << " clash_score: seq_sep: " << ij_seq_sep << " less than seq_seq: " << seq_sep << " skipped " << std::endl;
					continue;
				}
				if ( dist <= clash_dist ) {
					//tr_overlap_score << "clash_score: clash found! " << std::endl;
					return clash_return_score;
				}
			} // non-overlap regions - should not have clashes
		} // pose2
	} // pose1

	//tr_overlap_score << " compat_score = " << sum_overlap_score << std::endl;
	return sum_overlap_score;
}



// 140221: this is being called right now
core::Real
overlap_score(
	core::pose::Pose const &pose1,
	core::Size const pos1,
	core::pose::Pose const &pose2,
	core::Size const pos2,
	// overlap_score = -1*overlap_wt*exp(-1*dist_ij*steepness_wt)
	core::Real overlap_wt=1.0,
	core::Real steepness_wt=1.0,
	int seq_sep=5,
	core::Real clash_dist=3.0,
	core::Real clash_return_score=500// clash_score = when there is a clash after seq_sep 5, return 500 immediately
){
	//tr_overlap_score << "overlap_score got called" << std::endl;
	if ( pos1 == pos2 ) {
		return 0;
	}

	// because of the numbering matching issue - the pos2 must be further than pos1
	core::pose::Pose const &pose_lower = ( pos1 > pos2 ) ? pose2 : pose1;
	core::pose::Pose const &pose_upper = ( pos1 > pos2 ) ? pose1 : pose2;
	core::Size pos_lower = std::min( pos1, pos2 );
	core::Size pos_upper = std::max( pos1, pos2 );

	core::Real overlap_score( 0.0 );
	core::Real sum_overlap_score( 0.0 );

	// this is for speed up
	core::Size pose_lower_nrsds = pose_lower.size();
	core::Size pose_upper_nrsds = pose_upper.size();

	int offset = pos_upper - pos_lower; // should therefore always be positive

	// return 0 if fragments are not overlapping
	// this should be check before calling this overlap function, just to be safe
	//if ( ! isOverlapping( pose_lower, pos_lower, pose_upper, pos_upper ) )
	if ( ! isOverlapping( pose_lower, pos_lower, pos_upper ) ) {
		return 0;
	}


	for ( core::Size irsd=1; irsd <= pose_lower_nrsds; ++irsd ) {
		numeric::xyzVector<core::Real> atm_i = pose_lower.residue( irsd ).atom(2).xyz();

		for ( core::Size jrsd=1; jrsd <= pose_upper_nrsds; ++jrsd ) {
			numeric::xyzVector<core::Real> atm_j = pose_upper.residue( jrsd ).atom(2).xyz();
			core::Real dist = ( atm_i - atm_j ).length();

			//tr_overlap_score << "irsd: " << pos_lower+irsd-1 << " jrsd: " << pos_lower+jrsd+offset-1 << " ";
			if ( irsd == jrsd+offset ) {
				overlap_score = -1*overlap_wt*std::exp( -1*dist*steepness_wt);
				sum_overlap_score += overlap_score;
				//overlap_natoms += 1;
				//tr_overlap_score << " overlap_score: " << overlap_score << " natoms: " << overlap_natoms << std::endl;
			} else { // residues overlap
				// ij_seq_sep less than seq_sep don't do clash check
				int ij_seq_sep = std::abs( int( irsd - (jrsd+offset) ) );
				if ( ij_seq_sep <= seq_sep ) {
					//tr_overlap_score << "clash_score: seq_sep: " << ij_seq_sep << " less than seq_seq: " << seq_sep << " skipped" << std::endl;
					continue;
				}
				if ( dist <= clash_dist ) {
					//tr_overlap_score << "clash_score: clash found! " << std::endl;
					return clash_return_score;
				}
			} // non-overlap regions
		} // pose2
	} // pose1

	//tr_overlap_score << " compat_score = " << sum_overlap_score << std::endl;
	return sum_overlap_score;
}



#endif
