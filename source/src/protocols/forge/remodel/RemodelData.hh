// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/remodel/RemodelData.hh
/// @brief
/// @author Possu Huang (possu@u.washington.edu)
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_remodel_RemodelData_hh
#define INCLUDED_protocols_forge_remodel_RemodelData_hh

#include <utility/vector1.hh>

#include <ObjexxFCL/FArray1D.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/chemical/AA.hh>

#include <protocols/loops/Loops.hh>

// C++ headers
#include <vector>

namespace protocols {
namespace forge {
namespace remodel {

struct LineObject {
	int index;
	int original_index;
	std::string resname;
	std::string sstype;
	std::string design_type;
	bool isDesignable;
	bool has_constraints;
	bool has_ncaa;
	std::vector< std::string > ncaaList;
	std::vector< std::string > constraint_definition;
	std::vector< core::chemical::AA > aminoAcidList;
};

// this class stores information in the blueprint file
class RemodelData {

public:

	//constructor
	RemodelData();

	void getLoopsToBuildFromBlueprint( std::string text_blueprint );

	void getLoopsToBuildFromFile( std::string filename );

	void splitString( std::string str, std::string delim, std::vector< std::string > & results );

	void updateWithDsspAssignment( ObjexxFCL::FArray1D_char & dsspSS );

	void collectInsertionPose();

	void translateDSSP_ABEGO(std::string & ss, std::string & abego);


	protocols::loops::Loops loops_to_build;
	std::string sequence;

	// need this for "." switch to find remodel regions
	std::string ss;

	std::string abego;

	// merge the dssp assignment with ss string, exclude ".",
	// gets the final dssp_updated_ss for fragment pick
	std::string dssp_updated_ss;

	// 1 is fully auto
	// 2 is semi-auto (only design rebuilt and no neighbors),
	// 3 is manual which require resfile like assignments
	bool has_design_info_; //essential

	int pdb_start;
	int pdb_stop;

	int design_mode; // maybe defunct

	bool auto_design; // maybe defunct
	bool design_neighbor; // maybe defunct

	core::kinematics::MoveMap natro_movemap_;

	std::vector< protocols::forge::remodel::LineObject > blueprint;

	std::vector<core::Size> disulfMobileRange;
	std::vector<core::Size> disulfLandingRange;

	std::string parsed_string_for_resfile;
	// insertion related variables below. these vars only used when domain insertion is being done with remodel.
	core::pose::Pose insertPose;
	int insertionSize;
	std::string insertionSS;
	/*
	core::pose::Pose insertPose2;
	int insertion2Size;
	std::string insertion2SS;
	*/
	float total_chain_break_score;

};


} // remodel
} // forge
} // protocols

#endif
