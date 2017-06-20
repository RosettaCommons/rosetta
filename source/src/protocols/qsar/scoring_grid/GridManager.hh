// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/scoring_grid/GridManager.hh
/// @author Sam DeLuca
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_qsar_scoring_grid_GridManager_hh
#define INCLUDED_protocols_qsar_scoring_grid_GridManager_hh


#include <protocols/qsar/scoring_grid/GridManager.fwd.hh>
#include <protocols/qsar/scoring_grid/GridSet.fwd.hh>
#include <protocols/qsar/qsarMap.fwd.hh>
#include <protocols/qsar/scoring_grid/GridBase.fwd.hh>
#include <protocols/qsar/scoring_grid/ScoreNormalization.fwd.hh>

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/UltraLightResidue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>
#include <utility/json_spirit/json_spirit_writer.h>
#include <utility/json_spirit/json_spirit_reader.h>

// C++ Headers
#include <string>
#include <map>

namespace protocols {
namespace qsar {
namespace scoring_grid {

typedef std::map<std::string,core::Real> ScoreMap;

class GridManager : public utility::SingletonBase< GridManager >
{
public:
	friend class utility::SingletonBase< GridManager >;

	/// @brief Get a set of grids for the given pose, centered on the given center, including/excluding the given chains.
	GridSetCOP get_grids(
		GridSet const & prototype,
		core::pose::Pose const & pose,
		core::Vector const & center,
		std::string const & chain,
		bool exclude = true );

	/// @brief Get a set of grids for the given pose, centered on the given center, including/excluding the given chains.
	GridSetCOP get_grids(
		GridSet const & prototype,
		core::pose::Pose const & pose,
		core::Vector const & center,
		char chain,
		bool exclude = true );

	/// @brief Get a set of grids for the given pose, centered on the given center, including/excluding the given chains,
	/// using the given prototype as the prototype.
	GridSetCOP get_grids(
		GridSet const & prototype,
		core::pose::Pose const & pose,
		core::Vector const & center,
		utility::vector1< std::string > chains,
		bool exclude = true );

private:

	GridManager();
	GridManager(GridManager const &) = delete;
	GridManager const & operator = (GridManager const & ) = delete;

	/// @brief Insert the given GridSet into the grid_set_cache under the given index value.
	void insert_into_cache( std::string const & hash_val, GridSetOP const & grid_set );

private:

	std::map< std::string, GridSetCOP > grid_set_cache_;

};

}
}
}

#endif
