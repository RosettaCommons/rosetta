// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/TiltMover.hh
/// @brief      Tilts a protein in the membrane (Rosetta Scripts Hook)
/// @details Tilts a span, protein or part of a pose in the membrane,
///    depending on the jump number. The tilt axis is the axis
///    perpendicular to the axis connecting the embedding centers of the
///    two partners;
///    BEWARE: CANNOT USE MEMBRANE JUMP AS JUMP NUMBER!!!
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_TiltMover_hh
#define INCLUDED_protocols_membrane_TiltMover_hh

// Unit Headers
#include <protocols/membrane/TiltMover.fwd.hh>
#include <protocols/membrane/TiltMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.fwd.hh>

namespace protocols {
namespace membrane {

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;
using namespace core::pose;
using namespace protocols::moves;

/// @brief Tilts the downstream partner along the axis between
///   the COMs of the partners.
class TiltMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details Defaults: jump = 1, angle = random, axis =
	/// axis perpendicular to axis connecting protein embedding centers
	TiltMover();

	/// @brief Custom Constructor
	/// @details User can specify jump number
	TiltMover( core::Size jump_num );

	/// @brief Custom constructor
	/// @details User can specify jump number and angle
	TiltMover( core::Size jump_num, core::Real angle );

	/// @brief Copy Constructor
	TiltMover( TiltMover const & src );

	/// @brief Assignment Operator
	TiltMover & operator = ( TiltMover const & src );

	/// @brief Destructor
	~TiltMover() override;

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (TiltMover)
	std::string get_name() const override;

	/// @brief Flip the downstream partner in the membrane
	void apply( core::pose::Pose & pose ) override;

	/// @brief Set Random tilt angle between -20 and 20 degrees to keep
	///   protein oriented in the membrane correctly
	void set_random_membrane_tilt_angle();

private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options from Command Line
	void register_options();

	/// @brief Set default values
	void set_defaults();

private: // data

	/// @brief Jump number
	core::Size jump_num_;

	/// @brief Rotation angle in degrees
	int angle_;

	/// @brief Random tilt angle between -20 and 20 degrees
	bool random_angle_;
};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_TiltMover_hh
