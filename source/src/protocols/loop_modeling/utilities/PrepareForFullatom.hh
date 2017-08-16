// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_loop_modeling_utilities_PrepareForFullatom_HH
#define INCLUDED_protocols_loop_modeling_utilities_PrepareForFullatom_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/utilities/PrepareForFullatom.fwd.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// RosettaScripts headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

/// @brief Convert a pose to fullatom mode for high-resolution loop modeling.
///
/// @details Before we start high-resolution loop modeling, we need to make
/// sure that all the residues in the pose have sidechains.  There are three
/// ways this could happen.  This first is that the pose is already in fullatom
/// mode and we don't have to do anything.  This is important to check, because
/// we don't want to mess anything up if low-resolution loop modeling was
/// skipped.  The second is by copying sidechains from the original structure,
/// except for those in regions where the backbone was moved, which have to be
/// repacked.  This is the default action.  The third option, which has to be
/// explicitly requested, is to simply repack everything.

class PrepareForFullatom : public LoopMover {

public:

	/// @brief Default constructor.
	PrepareForFullatom();

	/// @copydoc LoopMover::parse_my_tag
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		Pose const & pose) override;

	/// @copydoc LoopMover::get_name

	/// @brief Set the original pose.  This must be called before apply().
	void set_original_pose(Pose const & pose);

	/// @brief Get the score function to be used on the next call to apply().
	core::scoring::ScoreFunctionOP get_score_function();

	/// @brief Set the score function to be used on the next call to apply().
	void set_score_function(core::scoring::ScoreFunctionOP score_function);

	/// @brief Guarantee the pose is repacked before the fullatom stage.
	void set_force_repack(bool value);

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


protected:

	/// @brief Convert the given pose to centroid mode.
	bool do_apply(Pose & pose) override;

private:

	Pose original_pose_;
	bool force_repack_;

};

}
}
}


#endif

