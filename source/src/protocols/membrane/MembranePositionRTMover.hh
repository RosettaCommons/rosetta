// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		protocols/membrane/MembranePositionRTMover.hh
///
/// @brief		Membrane Position Rotation/Translation Mover
///	@details	Apply a uniform rigid translation & rotation of the
///				membrane defined by the center/normal coordinates
///				stored in the membrane.
///				Last Modified: 6/28/14
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_MembranePositionRTMover_hh
#define INCLUDED_protocols_membrane_MembranePositionRTMover_hh

// Unit Headers
#include <protocols/membrane/MembranePositionRTMover.fwd.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/UniformPositionMover.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <numeric/xyzVector.hh>

namespace protocols {
namespace membrane {

using namespace core;
using namespace protocols::moves;

/// @brief Membrane Position Translation-Rotation Mover
class MembranePositionRTMover : public Mover {
	
public:
	
	////////////////////
	/// Constructors ///
	////////////////////
	
	/// @brief Construct a Default Membrane Position Mover
	MembranePositionRTMover();
	
	/// @brief Custom Constructor
	/// @details Specify a new membrane center and normal
	///	to move this position to
	MembranePositionRTMover( Vector center, Vector normal );
	
	/// @brief Copy Constructor
	/// @details Make a deep copy of this mover object
	MembranePositionRTMover( MembranePositionRTMover const & src );
	
	/// @brief Assignment Operator
	/// @details Make a deep copy of this mover object, overriding the assignment operator
	MembranePositionRTMover &
	operator=( MembranePositionRTMover const & src );
	
	/// @brief Destructor
	~MembranePositionRTMover();
	
	/////////////////////
	/// Mover Methods ///
	/////////////////////
	
	/// @brief Get the name of this mover
	virtual std::string get_name() const;
	
	/// @brief Apply Rotation/Translation to Membrane
	/// @brief Translate the membrane position in this pose
	/// to the new center position, and rotate to new normal
	virtual void apply( Pose & pose );
	
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
	
private:
	
	// Store new normal/center pair
	Vector center_;
	Vector normal_;
	
};

/// @brief Membrane Position Rotation Move
/// @details Rotate the orientation of the membrane position to a new
/// normal position. Applies rotation to jump
class MembranePositionRotationMover : public Mover {
	
public:
	
	////////////////////
	/// Constructors ///
	////////////////////
	
	/// @brief Construct a Default Membrane Position Mover
	MembranePositionRotationMover();
	
	/// @brief Custom Constructor
	/// @details Specify a new normal to rotate membranes to
	///	to move this position to
	MembranePositionRotationMover( Vector normal );
	
	/// @brief Copy Constructor
	/// @details Make a deep copy of this mover object
	MembranePositionRotationMover( MembranePositionRotationMover const & src );
	
	/// @brief Assignment Operator
	/// @details Make a deep copy of this mover object, overriding the assignment operator
	MembranePositionRotationMover &
	operator=( MembranePositionRotationMover const & src );
	
	/// @brief Destructor
	~MembranePositionRotationMover();
	
	/////////////////////
	/// Mover Methods ///
	/////////////////////
	
	/// @brief Get the name of this mover
	virtual std::string get_name() const;
	
	/// @brief Apply Rotation to Membrane
	/// @brief Rotate the membrane to the new normal position
	virtual void apply( Pose & pose );
	
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
	
private:
	
	// Store new normal
	Vector normal_;
	
};

/// @brief Membrane Position Translation Move
/// @details Translate the center of the membrane stub ot the specified position
class MembranePositionTranslationMover : public Mover {
	
public:
	
	////////////////////
	/// Constructors ///
	////////////////////
	
	/// @brief Construct a Default Membrane Position Mover
	MembranePositionTranslationMover();
	
	/// @brief Custom Constructor
	/// @details Specify a new center position to translate this stub to
	MembranePositionTranslationMover( Vector center );
	
	/// @brief Copy Constructor
	/// @details Make a deep copy of this mover object
	MembranePositionTranslationMover( MembranePositionTranslationMover const & src );
	
	/// @brief Assignment Operator
	/// @details Make a deep copy of this mover object, overriding the assignment operator
	MembranePositionTranslationMover &
	operator=( MembranePositionTranslationMover const & src );
	
	/// @brief Destructor
	~MembranePositionTranslationMover();
	
	/////////////////////
	/// Mover Methods ///
	/////////////////////
	
	/// @brief Get the name of this mover
	virtual std::string get_name() const;
	
	/// @brief Apply Translation to membrane position
	/// @brief Translate membrane position to new center
	virtual void apply( Pose & pose );
	
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
	
private:
	
	// Store new center
	Vector center_;
	
};

} // membrane
} // protocols

#endif // #ifndef INCLUDED_protocols_membrane_MembranePositionRTMover_hh
