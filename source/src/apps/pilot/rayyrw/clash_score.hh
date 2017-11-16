#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>

#include <apps/pilot/rayyrw/util.hh>
#include <apps/pilot/rayyrw/rms_util.hh>

#ifndef apps_pilot_rayyrw_clash_score_hh
#define apps_pilot_rayyrw_clash_score_hh

static basic::Tracer tr_clash_score("clash_score");



//  just count how many clashes for two fragments
core::Size
std_clash_check(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	core::Real const clash_dist_threshold,
	core::Size const n_clashes_allowed=10 // this is just silly right? fragment length is only 9
){
	//core::Real score(0.0);
	core::Size clash_natoms( 0 );
	core::Size pose2_nrsds = pose2.size();

	for ( core::Size irsd=1; irsd <= pose1.size(); ++irsd ) {
		numeric::xyzVector<core::Real> atm_i = pose1.residue( irsd ).atom( 2 ).xyz();

		for ( core::Size jrsd=1; jrsd <= pose2_nrsds; ++jrsd ) {
			numeric::xyzVector<core::Real> atm_j = pose2.residue( jrsd ).atom( 2 ).xyz();
			core::Real dist = ( atm_i - atm_j ).length();

			//tr_clash_score << "irsd: " << irsd << " jrsd: " << jrsd << " dist: " << dist << std::endl;

			if ( dist <= clash_dist_threshold ) {
				//score += cal_score( dist );
				clash_natoms += 1;
				//tr_clash_score << "Found clash at: irsd " << irsd << " jrsd " << jrsd << " dist: " << dist << std::endl;
			}

			if ( clash_natoms > n_clashes_allowed ) {
				return clash_natoms;
			}

		} // j
	} // i
	//return score;
	return clash_natoms;
}




// for gap_size == 1, do a soften clash check - skipping
core::Size
soften_clash_check(
	core::pose::Pose const & pose_lower,
	core::pose::Pose const & pose_upper,
	core::Real const clash_dist_threshold,
	core::Size const n_clashes_allowed=10
){
	core::Size clash_natoms( 0 );
	//core::Real score(0.0);

	// here are the differences
	// - skip checking the last residue of the further fragment and the first residue of the latter one
	for ( core::Size irsd=1; irsd <= pose_lower.size()-1; ++irsd ) {
		for ( core::Size jrsd=2; jrsd <= pose_upper.size(); ++jrsd ) {

			core::Real dist = cal_distance( pose_lower, irsd, pose_upper, jrsd );
			// tr_clash_score << "irsd: " << irsd << " jrsd: " << jrsd << " dist: " << dist << std::endl;

			if ( dist <= clash_dist_threshold ) {
				//score += cal_score( dist );
				clash_natoms += 1;
				// tr_clash_score << "Found clash at: irsd " << irsd << " jrsd " << jrsd << " dist: " << dist << std::endl;
			}

			if ( clash_natoms > n_clashes_allowed ) {
				return clash_natoms;
			}

		} // j
	} // i
	//return score;
	return clash_natoms;
}


// clash score that don't count clash scores from the terminus residues
// this is using soften_clash_check when the gap size is 1
core::Size
clash_score(
	core::pose::Pose const &pose1,
	core::Size const pos1,
	core::pose::Pose const &pose2,
	core::Size const pos2,
	core::Real const clash_dist_threshold
){
	//tr_clash_score << "clash_score got called" << std::endl;
	if ( pos1 == pos2 ) return 0; // this is probably not going to happen since the script setup has ingored it

	//  find pose_lower and pose_upper; important for soften_clash_check
	//  because we ignore the junction region for (pos_lower+pose_lower.n_residues()-1)-1 and pos_upper+1
	core::pose::Pose const &pose_lower = ( pos1 > pos2 ) ? pose2 : pose1;
	core::pose::Pose const &pose_upper = ( pos1 > pos2 ) ? pose1 : pose2;
	core::Size pos_lower = std::min( pos1, pos2 );
	core::Size pos_upper = std::max( pos1, pos2 );

	int gap_size = pos_upper - ( pos_lower + pose_lower.size() - 1 );

	if ( gap_size <= 0 ) { // from one overlap to...
		return 0;

	} else if ( gap_size == 1 ) {
		return soften_clash_check( pose_lower, pose_upper, clash_dist_threshold );

	} else {
		return std_clash_check( pose1, pose2, clash_dist_threshold ); // positions doesn't matter
	}
}


// given two fragments - do they have clashes?
core::Real
cal_score(
	core::Real dist
){
	core::Real score = 2*( std::pow( (( dist - 2.0 ) / 2.0 ), 2 ) );
	std::cout << "clash_score: " << score << std::endl;
	return std::min( 1.5, score );
}


core::Real
clash_score(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	core::Real const clash_dist_threshold=2,
	core::Size const n_clashes_allowed=10 // this will be useful when do domain assignments
){
	return std_clash_check( pose1, pose2, clash_dist_threshold, n_clashes_allowed );
}

#endif
