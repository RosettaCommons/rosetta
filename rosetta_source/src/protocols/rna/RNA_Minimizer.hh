// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file loopRNA_minimizer.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_RNA_Minimizer_HH
#define INCLUDED_protocols_rna_RNA_Minimizer_HH

#include <protocols/moves/Mover.hh>
#include <protocols/toolbox/AllowInsert.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

//// C++ headers
#include <string>

#include <core/kinematics/MoveMap.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace rna {

/// @brief The RNA de novo structure modeling protocol
class RNA_Minimizer: public protocols::moves::Mover {
public:
	/// @brief Construct the protocol object given
	/// the RNA fragment library to use.
	RNA_Minimizer();

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const {
		return new RNA_Minimizer(*this);
	}

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	void deriv_check( bool const setting ){ deriv_check_ = setting; }

	void use_coordinate_constraints( bool const setting ){ use_coordinate_constraints_ = setting; }

	void skip_o2star_trials( bool const setting ){ skip_o2star_trials_ = setting; }

	void set_perform_minimizer_run( bool const setting ){ perform_minimizer_run_ = setting; }

	void vary_bond_geometry( bool const setting ){ vary_bond_geometry_ = setting; }

	void set_include_default_linear_chainbreak( bool const setting){ include_default_linear_chainbreak_ = setting; }

	void set_verbose( bool const setting){ verbose_ = setting; }

	void set_do_dump_pdb( bool const setting){ do_dump_pdb_ = setting; }

	void set_min_type( std::string const setting){ min_type_ = setting; }

	void
	set_allow_insert(toolbox::AllowInsertOP allow_insert  );

	void
	set_score_function( core::scoring::ScoreFunctionOP const & scorefxn );

	core::scoring::ScoreFunctionOP const &
	score_function() const{ return scorefxn_; }

private:

	// Make this a Mover?
	void
	o2star_trials( core::pose::Pose & pose, core::scoring::ScoreFunctionOP const & scorefxn ) const;

	void
	setup_movemap( core::kinematics::MoveMap & mm, core::pose::Pose & pose );

	//refactored out to ConstrainToIdealMover
	/*void
	vary_bond_geometry(
		core::kinematics::MoveMap & mm,
		core::pose::Pose & pose,
		core::pose::Pose const & pose_reference );

	void
	create_pose_reference(
		core::pose::Pose const & pose,
		core::pose::Pose & pose_reference );*/


	bool deriv_check_;
	bool use_coordinate_constraints_;
	bool skip_o2star_trials_;
	bool perform_minimizer_run_;
	bool vary_bond_geometry_;
	bool include_default_linear_chainbreak_;
	bool verbose_;
	bool do_dump_pdb_;
	std::string min_type_;

	toolbox::AllowInsertOP allow_insert_;

	core::scoring::ScoreFunctionOP scorefxn_;


}; // class RNA_Minimizer

} //rna
} // protocols

#endif
