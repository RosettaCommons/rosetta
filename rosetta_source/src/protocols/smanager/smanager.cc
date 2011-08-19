// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/smanager/smanager.cc
///
/// @brief  Score manager
/// @author Sergey Lyskov

#include <protocols/smanager/smanager.hh>

#include <core/pose/Pose.hh>

#include <sstream>
#include <fstream>
#include <iomanip>

#ifdef WIN32
// apparently this is required for a Visual Studio build.
#include <core/conformation/Residue.hh>
#endif

namespace protocols {
namespace smanager {

void ScoreManager::add(Score &s, core::pose::Pose const &pose)
{
	s.file_name = get_file_name(Scores.size());

	Scores.push_back(s);

	if ( false ) {
		pose.dump_pdb( s.file_name );
	}
	write_score_file("scores.txt");
}

std::string ScoreManager::get_file_name(int decoy_num)
{
	std::ostringstream file_name;

	file_name << "D-" << std::setw(7) << std::setbase(16) << std::setfill('0') << decoy_num << ".pdb";
	//std::cout << "File name: " << file_name.str() << "\n";
	return file_name.str();
}


void ScoreManager::write_score_file(std::string fname)
{
	std::ostringstream S;

	S << "File_name     score    RMS\n";

	for(core::Size i=0; i<Scores.size(); i++) {
		S << Scores[i].file_name << "  " << Scores[i].score << "  " << Scores[i].RMS << "\n";
	}

	std::ofstream file(fname.c_str(), std::ios::out | std::ios::binary);
	if(!file) {
		std::cout << "Error: ScoreManager::write_score_file: Unable to open file:" << fname << " for writing!!!\n";
		return;
	}

	file.write( S.str().c_str(), S.str().size() );

	file.close();
}

} // smanager
} // protocols

