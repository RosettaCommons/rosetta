// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/RemodelWorkingSet.hh
/// @brief
/// @author Possu Huang (possu@u.washington.edu)
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_remodel_RemodelWorkingSet_hh
#define INCLUDED_protocols_forge_remodel_RemodelWorkingSet_hh

#include <ObjexxFCL/FArray2D.hh>

#include <utility/vector1.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/remodel/RemodelData.hh>


namespace protocols{
namespace forge{
namespace remodel{

// this class holds all the info for the model pose
// in the future a new variable might be added to delete a certain jump
// to create domain assembly type of fold-tree
class RemodelWorkingSet {

public:

	// default constructor
	RemodelWorkingSet();

	// copy constrctor
	RemodelWorkingSet( RemodelWorkingSet const & rval );

	// copy assignment
	RemodelWorkingSet & operator = ( RemodelWorkingSet const & rval );

	~RemodelWorkingSet(){};

	void workingSetGen( core::pose::Pose const & input_pose, protocols::forge::remodel::RemodelData const & data );
	void manualPackerTaskGen( core::pose::Pose const & built_pose, protocols::forge::remodel::RemodelData const & data );

	//void design_matrix_from_blueprint( std::vector<protocols::forge::remodel::LineObject>  blueprint ); //manual
	//void setup_auto_design_matrix( core::pose::Pose const & model_pose, std::vector<protocols::forge::remodel::LineObject> const & blueprint, bool const core, bool const boundary, bool surface );

	//void setup_repack_residues(core::pose::Pose & model_pose, std::vector<protocols::forge::remodel::LineObject> const & blueprint);
	//void createDisulfideBuildingData(core::pose::Pose const & model_pose, protocols::forge::remodel::RemodelData const & remodel_data);
	//void updatePoseWithARandomDisulfideJump(core::pose::Pose & model_pose);
	//void makeDisulfPairs(core::pose::Pose & model_pose);

	/// @brief If remodel loop setup is calling for n-terminus movement then return true, otherwise false.
	//bool moving_n_terminus() const;

	/// @brief build a fold tree for loop modeling using the defined loops
	/// @note  builds directly from internal loop data and does not randomize cutpoints,
	///        so inside a fully stochastic loop building routine you most likely
	///        *do not* want to use this function
	//core::kinematics::FoldTree standard_loop_fold_tree() const;


	// Revisit by Sachko 03/29/2013
	// Whey are these all public???
	///
	protocols::loops::Loops loops;
	int safe_root_;
	std::string sequence;
	std::string ss;
	std::string abego;
	std::map<int,int> translate_index;

	std::vector<int> begin;
	std::vector<int> end;
	std::vector<int> copy_begin;
	std::vector<int> copy_end;
	std::vector<int> src_begin;
	std::vector<int> src_end;

	std::string aa; // for assigning generic type used for building.
	bool hasInsertion;

	protocols::forge::build::BuildManager manager;
	core::pack::task::PackerTaskOP task;

	// disulfide building
	core::pose::Pose rvjump_pose;
	bool buildDisulfide;
	ObjexxFCL::FArray2D_int disulfide_jump_points;
	int disulfide_cutpoint;
	std::string disulfide_ss;

	int insertionStartIndex;
	int insertionEndIndex;
	ObjexxFCL::FArray2D_bool design_matrix;

};

class Segment {

public:
	std::vector<int> residues;
};


}
}
}

#endif
