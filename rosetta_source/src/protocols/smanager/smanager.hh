// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/smanager/smanager.hh
///
/// @brief  Score manager
/// @author Sergey Lyskov

#ifndef INCLUDED_protocols_smanager_smanager_hh
#define INCLUDED_protocols_smanager_smanager_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>


#include <vector>
#include <string>


namespace protocols {
namespace smanager {

typedef core::Real Real;

struct Score
{
	std::string file_name;
	Real score;
	Real RMS;
};

///
/// This is very first draf, just to see what we really need to do.
///
class ScoreManager
{
public:
	ScoreManager() {}

	/// Add given score line to score list
	void add(Score &s, core::pose::Pose const &pose);


	void write_score_file(std::string fname);

private:
	std::string get_file_name(int decoy_num);

	std::vector<Score> Scores;

};


} // smanager
} // protocols


#endif // INCLUDED_protocols_smanager_SManager_HH


