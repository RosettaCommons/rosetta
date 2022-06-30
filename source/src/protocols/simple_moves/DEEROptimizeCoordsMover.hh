// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/simple_moves/DEEROptimizeCoordsMover.hh
/// @brief
/// @author Diego del Alamo

#include <random>
#include <iostream>

#include <protocols/simple_moves/DEEROptimizeCoordsMover.fwd.hh>
#include <protocols/simple_moves/DEEROptimizeCoordsMoverCreator.hh>

#include <basic/datacache/DataMap.hh>

#include <core/scoring/epr_deer/metrics/DEERData.hh>
#include <core/scoring/epr_deer/DEERDataCache.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/filters/Filter.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>

#include <string>
#include <cstdlib>
#include <random>

#if defined(_OPENMP)
#include "omp.h"
#endif // OPENMPI

namespace protocols {
namespace simple_moves {

class DEEROptimizeCoordsMover : public moves::Mover {

	using PairSizeString = std::pair< core::Size, std::string >;
	using ResPair = std::pair< PairSizeString, PairSizeString >;
	using SizePair = std::pair< core::Size, core::Size >;
	using PseudoSL = std::pair< numeric::xyzVector< core::Real >, core::Real >;
	using CoordModel = std::map< PairSizeString, utility::vector1< core::Real > >;

public:

	DEEROptimizeCoordsMover() ;

	~DEEROptimizeCoordsMover() override;

	moves::MoverOP
	clone() const override;

	void
	apply(
		core::pose::Pose & pose
	) override;

	/// @brief  Initializes coords used throughout method
	/// @param  pose: Pose being labeled
	/// @return Coords that pass clash evaluation
	std::map< PairSizeString, utility::vector1< PseudoSL > >
	init_coords(
		core::pose::Pose & pose,
		core::scoring::epr_deer::DEERDataCache const & datacache
	);

	/// @brief Initializes the data saved in the Mover
	/// @param pose: Pose being modified
	/// @param  datacache: Datacache
	/// @detail ResPair is std::pair< PairSizeString, PairSizeString >
	/// @detail PairSizeString is std::pair< core::Size, std::string >
	/// @detail SizePair is std::pair< core::Size, core::Size >
	/// @detail Map used: data[ resPair ][ std::make_pair( res1_id, res2_id ) ]
	std::map< ResPair, std::map< SizePair, utility::vector1< core::Real > > >
	init_data(
		std::map< PairSizeString, utility::vector1< PseudoSL > > const & default_coords,
		core::scoring::epr_deer::DEERDataCache & datacache
	);

	/// @brief Initialize and run the optimization procedure
	/// @param  start_resolution: Starting resolution for modifications
	/// @param  end_resolution: End resolution for modifications
	/// @return Best coordinate model obtained (by sum of squared residuals)
	CoordModel
	init_model(
		std::map< PairSizeString, utility::vector1< PseudoSL > > const & coords,
		core::scoring::epr_deer::DEERDataCache & datacache
	);

	/// @brief Performs an initial search (for non-3D backgrounds)
	/// @param  model: Model to refine
	/// @param  resolution: Refinement resolution
	/// @param  trace_map: Map of DEER traces for calculation
	/// @param  datacache: Datacache with data
	void
	initial_search(
		CoordModel & model,
		core::Real & resolution,
		std::map< ResPair, std::map< SizePair, utility::vector1< core::Real > > > const & trace_map,
		core::scoring::epr_deer::DEERDataCache & datacache
	);

	/// @brief Minimizes the model to fit the data
	/// @param  model: Model to refine
	/// @param  resolution: Refinement resolution
	/// @param  trace_map: Map of DEER traces for calculation
	/// @param  datacache: Datacache with data
	/// @return Score of the model
	core::Real
	minimize(
		CoordModel & model,
		core::Real & resolution,
		std::map< ResPair, std::map< SizePair, utility::vector1< core::Real > > > const & trace_map,
		core::scoring::epr_deer::DEERDataCache & datacache
	);

	bool
	boltzmann(
		core::Real const & diff,
		core::Real const & temp = 1.0
	) const;

	/// @brief Recursive function to normalize the weights of a model
	/// @param model: The model being normalized
	/// @param res: Residue to normalize
	/// @detail If residue is zero, the function calls itself for each residue
	void
	normalize(
		CoordModel & model,
		PairSizeString const & res = std::make_pair( 0, "" )
	) const;

	core::Real
	score(
		CoordModel const & coords,
		std::map< ResPair, std::map< SizePair, utility::vector1< core::Real > > > const & trace_map,
		core::scoring::epr_deer::DEERDataCache & datacache,
		bool const print = false
	);

	void
	print_coords(
		CoordModel const & model,
		std::map< PairSizeString, utility::vector1< PseudoSL > > const & default_coords,
		core::Real const & current_score,
		core::scoring::epr_deer::DEERDataCache & datacache,
		std::map< ResPair, std::map< SizePair, utility::vector1< core::Real > > > const & trace_map,
		core::pose::Pose & pose
	);

	std::string
	get_name() const override;

	std::string
	mover_name() const;

	void
	provide_xml_schema(
		utility::tag::XMLSchemaDefinition & xsd
	);

	void
	parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &
	) override;

	std::string
	keyname() const;

	moves::MoverOP
	create_mover() const;

	void
	provide_xml_schema(
		utility::tag::XMLSchemaDefinition & xsd
	) const;


};

}
}
