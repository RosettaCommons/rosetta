// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/symmetric_docking/membrane/MPSymDockMover.hh
///
/// @brief      Membrane symmetric docking protocol
/// @details    Given an asymmetric starting pose, setup the pose for symmetry,
///             add a membrane representation for the full symmetric complex,
///             position the pose at the center of mass of all transmembrane spans,
///             and proceed with the standard symmetric docking protocol.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 2/9/15

#ifndef INCLUDED_protocols_symmetric_docking_MPSymDockMover_hh
#define INCLUDED_protocols_symmetric_docking_MPSymDockMover_hh

// Unit Headers
#include <protocols/symmetric_docking/membrane/MPSymDockMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package Headers
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <utility/tag/Tag.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace symmetric_docking {
namespace membrane {

class MPSymDockMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Defualt constructor for the membrane protein symmetric
	/// docking protocol
	MPSymDockMover();

	/// @brief Destructor
	~MPSymDockMover();

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP
	clone() const;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void
	parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	//////////////////////
	/// Mover Methods  ///
	//////////////////////

	/// @brief Return the name of this mover (MPSymDockMover)
	virtual
	std::string
	get_name() const;

	/// @brief Apply Method: Symmetric Docking in the membrane
	/// @details Setup pose for symmetry, add the membrane components to the total pose,
	/// and then perform symmetric docking
	virtual
	void
	apply( Pose & pose );

private:

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Position Pose at TopologyCOM
	/// @details Position the initial pose at the center of mass of
	/// predicted transmembrane spanning regions of the whole complex
	void
	position_membrane_at_topology_com( core::pose::Pose & pose );

};

} // membrane
} // symmetric_docking
} // protocols

#endif // INCLUDED_protocols_symmetric_docking_MPSymDockMover_hh


