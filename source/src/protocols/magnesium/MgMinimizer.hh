// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/magnesium/MgMinimizer.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_magnesium_MgMinimizer_HH
#define INCLUDED_protocols_magnesium_MgMinimizer_HH

#include <protocols/moves/MoverForPoseList.hh>
#include <protocols/magnesium/MgMinimizer.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

namespace protocols {
namespace magnesium {

class MgMinimizer: public protocols::moves::MoverForPoseList {

public:

	//constructor
	MgMinimizer();

	//destructor
	~MgMinimizer();

public:

	virtual void apply( core::pose::Pose & pose );

	using protocols::moves::MoverForPoseList::apply;

	virtual std::string get_name() const{ return "MgMinimizer"; }

	void set_mg_res( utility::vector1<Size> const & setting ){ mg_res_ = setting; }
	utility::vector1<Size> mg_res() const { return mg_res_; }

	void set_minimize_scorefxn( core::scoring::ScoreFunctionCOP const & setting ){ minimize_scorefxn_ = setting; }
	core::scoring::ScoreFunctionCOP minimize_scorefxn() const { return minimize_scorefxn_; }

	void set_mg_coord_cst_dist( core::Distance const & setting ){ mg_coord_cst_dist_ = setting; }
	core::Distance mg_coord_cst_dist() const { return mg_coord_cst_dist_; }

private:

	core::kinematics::MoveMap
	get_mg_hoh_minimize_move_map( core::pose::Pose const & pose,
		utility::vector1< Size > const & mg_res ) const;

private:

	utility::vector1< Size > mg_res_; // use post numbering
	core::scoring::ScoreFunctionCOP minimize_scorefxn_;
	core::Distance mg_coord_cst_dist_;

};

} //magnesium
} //protocols

#endif
