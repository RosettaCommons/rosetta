// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_forge_remodel_RemodelData_hh
#define INCLUDED_protocols_forge_remodel_RemodelData_hh


#include <protocols/loops/Loops.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/chemical/AA.hh>

#include <ObjexxFCL/FArray1D.fwd.hh>

#include <utility/vector1.hh>



namespace protocols{
namespace forge{
namespace remodel{

       struct LineObject
        {
          int index;
          int original_index;
          std::string resname;
          std::string sstype;
          std::string design_type;
          bool isDesignable;
				  bool has_constraints;
					std::vector<std::string> constraint_definition;
          std::vector<core::chemical::AA> aminoAcidList;
        };

		// this class stores information in the blueprint file
        class RemodelData
        {
        public:

					//constructor
						RemodelData();

            protocols::loops::Loops loops_to_build;
            std::string sequence;

            // need this for "." switch to find remodel regions
						std::string ss;

            // merge the dssp assignment with ss string, exclude ".",
						// gets the final dssp_updated_ss for fragment pick
						std::string dssp_updated_ss;

            // 1 is fully auto, 2 is semi-auto (only design rebuilt and no neighbors,
						// 3 is manual which require resfile like assignments
						bool has_design_info_; //essential

						int pdb_start;
						int pdb_stop;

						int design_mode; // maybe defunct

						bool auto_design; // maybe defunct
						bool design_neighbor; // maybe defunct

						core::kinematics::MoveMap natro_movemap_;

            std::vector<protocols::forge::remodel::LineObject> blueprint;

            std::vector<core::Size> disulfMobileRange;
            std::vector<core::Size> disulfLandingRange;

						std::string parsed_string_for_resfile;

						//insertion related below
            core::pose::Pose insertPose;
            int insertionSize;
            std::string insertionSS;

            float total_chain_break_score;

						void getLoopsToBuildFromFile();
            void splitString(std::string str, std::string delim, std::vector<std::string> & results);
            void updateWithDsspAssignment(ObjexxFCL::FArray1D_char & dsspSS);

            void collectInsertionPose();
        };
}
}
}

#endif
