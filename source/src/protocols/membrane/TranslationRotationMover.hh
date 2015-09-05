// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <protocols/membrane/TranslationRotationMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/rigid/RigidBodyMover.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/membrane/MembraneInfo.fwd.hh>
#include <core/conformation/membrane/SpanningTopology.fwd.hh>
#include <core/conformation/membrane/Span.fwd.hh>
#include <core/conformation/membrane/LipidAccInfo.hh>
#include <protocols/membrane/geometry/EmbeddingDef.fwd.hh>
#include <protocols/membrane/util.hh>

// Package Headers
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.fwd.hh>

namespace protocols {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::moves;

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

	/// @brief Copy Constructor
	TranslationMover( TranslationMover const & src );

	/// @brief Destructor
	virtual ~TranslationMover();

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	virtual protocols::moves::MoverOP clone() const;

	/// @brief Create a Fresh Instance of this Mover
	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (TranslationMover)
	virtual std::string get_name() const;

	/// @brief Translate the pose along the defined vector
	virtual void apply( core::pose::Pose & pose );

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

	/// @brief Copy Constructor
	RotationMover( RotationMover const & src );

	/// @brief Destructor
	virtual ~RotationMover();

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	virtual protocols::moves::MoverOP clone() const;

	/// @brief Create a Fresh Instance of this Mover
	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (RotationMover)
	virtual std::string get_name() const;

	/// @brief Rotate the pose - see above
	virtual void apply( core::pose::Pose & pose );

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

	/// @brief Copy Constructor
	TranslationRotationMover( TranslationRotationMover const & src );

	/// @brief Destructor
	virtual ~TranslationRotationMover();

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	virtual protocols::moves::MoverOP clone() const;

	/// @brief Create a Fresh Instance of this Mover
	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (TranslationRotationMover)
	virtual std::string get_name() const;

	/// @brief Translate the pose along the defined vector
	virtual void apply( core::pose::Pose & pose );

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
