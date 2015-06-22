// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/make_rot_lib/MakeRotLibMover.hh
/// @brief Header file for MakeRotLibMover
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

#ifndef INCLUDED_protocols_make_rot_lib_MakeRotLibMover_hh
#define INCLUDED_protocols_make_rot_lib_MakeRotLibMover_hh

// unit headers
#include <protocols/make_rot_lib/MakeRotLibMover.fwd.hh>
#include <protocols/make_rot_lib/MakeRotLibOptionsData.hh>
#include <protocols/make_rot_lib/RotData.hh>

// protocols headers
#include <protocols/moves/Mover.hh>

// core headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

#include <core/scoring/ScoreFunction.hh>

// c++ headers
#include <string>
#include <ostream>

namespace protocols {
namespace make_rot_lib {

class MakeRotLibMover : public protocols::moves::Mover {
public:

	// ctor
	MakeRotLibMover();

	// dtor
	virtual ~MakeRotLibMover(){}

	// pose interface
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "MakeRotLibMover"; }

	// setup methods
	void init_centroids( CentroidRotNumVecVec const & centroid_data, core::Size num_chi, bool semirotameric );
	void init_rotamers( TorsionRangeVec const & chi_ranges, core::Size num_clusters, core::Real omg, utility::vector1< core::Real> bbs, utility::vector1<core::Size> bbids, core::Real eps, bool semirotameric );

	// minmization methods
	void minimize_rotamer( RotData & rd, core::pose::Pose & pose, utility::vector1< core::Real> bbs, utility::vector1<core::Size> bb_ids,  MakeRotLibPolymerType polymer_type );
	void minimize_all_rotamers( core::pose::Pose & pose, utility::vector1< core::Real> bbs, utility::vector1<core::Size> bb_ids, MakeRotLibPolymerType polymer_type );

	// cluster loop
	void seed_centroids( bool semirotameric ); // for k-medoids only
	void calc_all_dist();
	bool calc_rotamer_clusters();
	bool calc_medoids();
	bool calc_centroids();
	bool bbs_appropriate_for_definitions( utility::vector1< core::Real > bbs );

	// finalizing methods
	void calc_final_rotamers();
	void calc_final_rotamer_probs();
	void calc_standard_deviations( core::pose::Pose & pose,utility::vector1< core::Real> bbs, utility::vector1<core::Size> bb_ids, MakeRotLibPolymerType polymer_type );

	// symmetry methods
	void make_two_fold_symmetry_135_315( RotDataVec & rdv, core::Size chi_num );
	void make_two_fold_symmetry_0_180( RotDataVec & rdv, core::Size chi_num );
	void make_three_fold_symmetry_90_210_330( RotDataVec & rdv, core::Size chi_num );

	// logging
	void print_rot_data( RotData & rd, utility::vector1< core::Size > bb_ids, std::ostream & os );
	void print_rot_data_vec( RotDataVec & rdv, utility::vector1< core::Size > bb_ids, std::ostream & os );
	core::Real print_avg_cluster_centroid_dist( std::ostream & os );
	void print_dunbrack02_rotlib( core::Real omg, utility::vector1< core::Real> bbs, utility::vector1<core::Size> bbids, core::Real eps, MakeRotLibPolymerType polymer_type, std::ostream & os );
	void print_definitions( std::ostream & os );
	void print_dunbrack10_rotlib( core::Real omg, utility::vector1< core::Real> bbs, utility::vector1<core::Size> bbids, core::Real eps, MakeRotLibPolymerType polymer_type, std::ostream & os );

	// utility
	core::Real calc_dist( RotData & rd1, RotData & rd2 );
	core::Real angle_diff( core::Real a1, core::Real a2 );
	void calc_running_avg( core::Real angle_new, core::Real & angle_old, core::Size & count );

private:

	core::scoring::ScoreFunctionOP scrfxn_;
	core::Real KbT_;
	RotDataVec centroids_;
	RotDataVec rotamers_;
	RotDataVec final_rotamers_;

};

}//make_rot_lib
}//protocols

#endif //INCLUDED_protocols_make_rot_lib_MakeRotLibMover_HH
