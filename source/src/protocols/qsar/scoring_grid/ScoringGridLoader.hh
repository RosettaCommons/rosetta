// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/parser/ScoringGridLoader.hh
/// @brief  Declartion of the XML parser's ScoringGridLoader class for adding named ScoringGrids to the basic::datacache::DataMap
/// @author Sam DeLuca
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com) -- moved here from DockDesignParser.cc

#ifndef INCLUDED_protocols_qsar_scoring_grid_ScoringGridLoader_HH
#define INCLUDED_protocols_qsar_scoring_grid_ScoringGridLoader_HH

// Package Headers
#include <protocols/qsar/scoring_grid/ScoringGridLoader.fwd.hh>
#include <protocols/jd2/parser/DataLoader.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace qsar {
namespace scoring_grid {

/// @brief A class for loading ScoringGrids into the XML parser's basic::datacache::DataMap.
class ScoringGridLoader : public protocols::jd2::parser::DataLoader
{
public:
	ScoringGridLoader();
	virtual ~ScoringGridLoader();

	/// @brief The ScoringGridLoader will create named ScoringGrids and load them into the basic::datacache::DataMap
	virtual
	void load_data(
		core::pose::Pose const & pose,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) const;

};

} //namespace scoring_grid
} //namespace qsar
} //namespace protocols

#endif //INCLUDED_protocols_qsar_scoring_grid_ScoringGridLoader_HH
