// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_loop_modeling_utilities_PrepareForCentroid_HH
#define INCLUDED_protocols_loop_modeling_utilities_PrepareForCentroid_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/utilities/PrepareForCentroid.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

/// @brief Convert a pose to centroid mode for low-resolution loop modeling.
class PrepareForCentroid : public LoopMover {

public:

	/// @brief Default constructor.
	PrepareForCentroid();

	/// @copydoc LoopMover::get_name

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

};

}
}
}


#endif

