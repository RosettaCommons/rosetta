// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief      Does translation and rotation of a pose or the membrane
/// @details Translates along a vector and rotates a pose or membrane around
///    the axis perpendicular to both vectors. Works for both
///    fixed protein/flexible membrane or fixed membrane/flexible protein
///    Takes the jump number as input, default is the membrane jump
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_TranslationRotationMover_hh
#define INCLUDED_protocols_membrane_TranslationRotationMover_hh

// Unit Headers
#include <protocols/membrane/TranslationRotationMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace membrane {

/// @brief Translation vector can be defined in -mp:setup center
///   flag to translate the new pose to. The mover is a general mover
///   but used mainly on membrane proteins, that's why we use this flag here
///   The default jump is going to be the membrane jump, but you can also
///   specify your own
class TranslationMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor, assumes old center = (0, 0, 0)
	TranslationMover();

	/// @brief Custom Constructor
	/// @details User can specify a translation vector
	TranslationMover(
		core::Vector translation_vector
	);


	/// @brief Custom Constructor
	/// @details User can specify a translation vector and a jump number;
	///   operation happens on the downstream stub
	TranslationMover(
		core::Vector translation_vector,
		core::Size jumpnum
	);

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
		basic::datacache::DataMap &
	) override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (TranslationMover)
	std::string get_name() const override;

	/// @brief Translate the pose along the defined vector
	void apply( core::pose::Pose & pose ) override;

private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options from Command Line
	/// @details Register mover-relevant options with JD2 - includes
	/// mp:setup options: center
	void register_options();

	/// @brief Initialize Mover options from the commandline
	/// @details Initialize mover settings from the commandline
	/// using the mp:setup center flag
	void init_from_cmd();


private: // data

	// translation vector
	core::Vector translation_vector_;

	// jump number
	core::Size jumpnum_;

};

/// @brief Rotates the pose such that a vector in the old orientation will be
///  overlayed in the new orientation. Requires two vectors (old and new)
///  and a point on the new vector around which the rotation takes place.
///  The mover is a general mover but used mainly on membrane proteins.
///  For membrane proteins, the two vectors will be the old and new normal
///  and the point will be the new center. The default jump is going to be
///  the membrane jump, but you can also specify your own. The rotation
///  happens around the axis perpendicular to both vectors with an angle
///  enclosed by both vectors.
class RotationMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details Uses the -mp:setup center and normal flags
	RotationMover();

	/// @brief Custom Constructor
	/// @details User can specify an old normal, a new normal, and a new center
	///   around which the rotation takes place
	RotationMover(
		core::Vector old_normal,
		core::Vector new_normal,
		core::Vector rot_center
	);


	/// @brief Custom Constructor
	/// @details User can specify an old normal, a new normal, and a new center
	///   around which the rotation takes place on this particular jump;
	///   operation happens on the downstream stub
	RotationMover(
		core::Vector old_normal,
		core::Vector new_normal,
		core::Vector rot_center,
		core::Size jumpnum
	);

	/// @brief Custom Constructor
	/// @details User can specify the same normal twice as old_normal_ and new_normal_, and a new center
	///   around which the rotation takes place on this particular jump;
	///   operation happens on the downstream stub. azimuthal_angle is the angle to be rotated about the normal.
	RotationMover(
		core::Vector old_normal,
		core::Vector new_normal,
		core::Vector rot_center,
		core::Size jumpnum,
		core::Real azimuthal_angle
	);

	/// @brief Copy Constructor
	RotationMover( RotationMover const & src );

	/// @brief Destructor
	~RotationMover() override;

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
		basic::datacache::DataMap &
	) override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (RotationMover)
	std::string get_name() const override;

	/// @brief Rotate the pose - see above
	void apply( core::pose::Pose & pose ) override;

private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options from Command Line
	/// @details Register mover-relevant options with JD2 - includes
	/// mp:setup options: center and normal
	void register_options();

	/// @brief Initialize Mover options from the commandline
	/// @details Initialize mover settings from the commandline
	/// using the mp:setup center and normal flags
	void init_from_cmd();


private: // data

	// vectors and point defining the rotation
	core::Vector old_normal_;
	core::Vector new_normal_;
	core::Vector rot_center_;

	// jump number
	core::Size jumpnum_;
	core::Real azimuthal_angle_;

};

/// @brief Translation and Rotation of a pose. The new position can be defined by
///  the -mp:setup center and normal flags. The mover is a
///  general mover, but used mainly on membrane proteins, that's why we
///  use this flag here. The default jump is going to be the membrane jump,
///  but you can also specify your own. See above for the TranslationMover
///  and the RotationMover
class TranslationRotationMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details Uses the -mp:setup center and normal flags
	TranslationRotationMover();

	/// @brief Custom Constructor
	/// @details User can specify an old normal and center and a new normal
	///   and center
	TranslationRotationMover(
		core::Vector old_center,
		core::Vector old_normal,
		core::Vector new_center,
		core::Vector new_normal
	);


	/// @brief Custom Constructor
	/// @details User can specify an old normal and center and a new normal
	///   and center and a jump; the downstream stub will be rotated
	TranslationRotationMover(
		core::Vector old_center,
		core::Vector old_normal,
		core::Vector new_center,
		core::Vector new_normal,
		core::Size jumpnum
	);

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
		basic::datacache::DataMap &
	) override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (TranslationRotationMover)
	std::string get_name() const override;

	/// @brief Translate the pose along the defined vector
	void apply( core::pose::Pose & pose ) override;

private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options from Command Line
	/// @details Register mover-relevant options with JD2 - includes
	/// mp:setup center and normal flags
	void register_options();

	/// @brief Initialize Mover options from the commandline
	/// @details Initialize mover settings from the commandline
	/// using the mp:setup center and normal flags
	void init_from_cmd();


private: // data

	// TranslationRotation vector
	core::Vector old_center_;
	core::Vector old_normal_;
	core::Vector new_center_;
	core::Vector new_normal_;

	// jump number
	core::Size jumpnum_;

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_TranslationRotationMover_hh
