// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoringManager.hh
/// @brief  Scoring manager class header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_ScoringManager_fwd_hh
#define INCLUDED_core_scoring_ScoringManager_fwd_hh

// Unit headers
#include <string> // doh

namespace core {
namespace scoring {

class ScoringManager;

extern std::string const FA_STANDARD_DEFAULT;
extern std::string const FA_STANDARD_SOFT;

extern std::string const UNFOLDED_SCORE12;
extern std::string const UNFOLDED_MM_STD;
extern std::string const UNFOLDED_RNA;

extern std::string const UNFOLDED_SPLIT_TALARIS2013;
extern std::string const UNFOLDED_SPLIT_MM_STD;

extern std::string const UNFOLDED_SPLIT_USER_DEFINED;

extern std::string const SPLIT_UNFOLDED_ELE;
extern std::string const SPLIT_UNFOLDED_PDB;
extern std::string const SPLIT_UNFOLDED_ROSETTA;
extern std::string const SPLIT_UNFOLDED_MM;
extern std::string const SPLIT_UNFOLDED_UNIQUE;

extern std::string const SPLIT_UNFOLDED_MEAN;
extern std::string const SPLIT_UNFOLDED_MEDIAN;
extern std::string const SPLIT_UNFOLDED_MODE;
extern std::string const SPLIT_UNFOLDED_BOLTZ;

extern std::string const SPLIT_UNFOLDED_USER_DEFINED;

} // namespace core
} // namespace scoring


#endif
