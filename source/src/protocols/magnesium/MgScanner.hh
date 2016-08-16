// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/magnesium/MgScanner.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_magnesium_MgScanner_HH
#define INCLUDED_protocols_magnesium_MgScanner_HH

#include <protocols/moves/Mover.hh>
#include <protocols/magnesium/MgScanner.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace magnesium {

class MgScanner: public protocols::moves::Mover {

public:

	//constructor
	MgScanner();

	//destructor
	~MgScanner();

public:

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const{ return "MgScanner"; }

	void set_minimize( bool const & setting ){ minimize_ = setting; }
	bool minimize() const { return minimize_; }

	void set_hydrate( bool const & setting ){ hydrate_ = setting; }
	bool hydrate() const { return hydrate_; }

	void set_minimize_during_scoring( bool const & setting ){ minimize_during_scoring_ = setting; }
	bool minimize_during_scoring() const { return minimize_during_scoring_; }

	void set_tether_to_closest_res( bool const & setting ){ tether_to_closest_res_ = setting; }
	bool tether_to_closest_res() const { return tether_to_closest_res_; }

	void set_integration_test( bool const & setting ){ integration_test_ = setting; }
	bool integration_test() const { return integration_test_; }

	void set_score_cut( core::Real const & setting ){ score_cut_ = setting; }
	core::Real score_cut() const { return score_cut_; }

	void set_score_cut_PDB( core::Real const & setting ){ score_cut_PDB_ = setting; }
	core::Real score_cut_PDB() const { return score_cut_PDB_; }

	void set_xyz_step( core::Real const & setting ){ xyz_step_ = setting; }
	core::Real xyz_step() const { return xyz_step_; }

	void set_scorefxn( core::scoring::ScoreFunctionOP const & setting ){ scorefxn_ = setting; }
	core::scoring::ScoreFunctionOP scorefxn() const { return scorefxn_; }

	void set_silent_file( std::string const & setting ){ silent_file_ = setting; }
	std::string silent_file() const { return silent_file_; }

	void set_output_pdb( std::string const & setting ){ output_pdb_ = setting; }
	std::string output_pdb() const { return output_pdb_; }

	void set_input_scan_res( utility::vector1< core::Size > const & setting ){ input_scan_res_ = setting; }
	utility::vector1< core::Size > input_scan_res() const { return input_scan_res_; }

	void set_minimize_mg_coord_constraint_distance( core::Distance const & setting ){ minimize_mg_coord_constraint_distance_ = setting; }
	core::Distance minimize_mg_coord_constraint_distance() const { return minimize_mg_coord_constraint_distance_; }


private:

	void
	scan_magnesiums( core::pose::Pose & pose );

	core::Real
	get_score( core::pose::Pose & pose,
		core::Vector const & mg_position,
		core::scoring::ScoreFunctionCOP scorefxn,
		bool const hydrate_magnesium = false,
		bool const keep_waters = false,
		bool const minimize = false );

	void
	cluster_mg();

	void
	output_mg_to_silent_file( std::string const & silent_file );

	void
	output_mg_into_one_PDB( core::pose::Pose const & pose );

	Size
	get_unique_mg_res( core::pose::Pose const & mg_pose );

	core::pose::PoseOP
	get_single_mg_pose();

	core::Distance
	distance_to_closest_magnesium( core::Vector const & mg_position,
		core::pose::Pose const & reference_pose );


private:

	utility::vector1< core::pose::PoseOP > mg_poses;
	bool hydrate_;
	bool minimize_;
	bool minimize_during_scoring_;
	bool tether_to_closest_res_;
	bool integration_test_;
	core::Real score_cut_;
	core::Real score_cut_PDB_;
	core::Real xyz_step_;
	std::string silent_file_;
	std::string output_pdb_;
	utility::vector1< core::Size > input_scan_res_;
	core::Distance minimize_mg_coord_constraint_distance_;

	core::scoring::ScoreFunctionOP scorefxn_;
	utility::vector1< core::pose::PoseOP > mg_poses_;

};

} //magnesium
} //protocols

#endif
