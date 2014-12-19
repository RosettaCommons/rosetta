// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/FlipMoverCreator.hh
/// @brief      Flips a span or protein in the membrane (Rosetta Scripts Hook)
/// @details	Flips a span, protein or part of a pose in the membrane,
///				depending on the jump number.
///				ONLY FOR FIXED MEMBRANE AND FLEXIBLE PROTEIN
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_FlipMover_hh
#define INCLUDED_protocols_membrane_FlipMover_hh

// Unit Headers
#include <protocols/membrane/FlipMover.fwd.hh>
#include <protocols/membrane/FlipMoverCreator.hh>
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

using namespace core;
using namespace core::pose;
using namespace protocols::moves;
	  
/// @brief	Takes a pose and flips the downstream partner around the axis between
///			the COMs of the partners, projected into the membrane plane.
///			CAUTION: THIS MOVER ONLY WORKS FOR A FIXED MEMBRANE WHERE THE
///			MEMBRANE VIRTUAL RESIDUE IS AT THE ROOT OF THE FOLDTREE!!!
class FlipMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details Defaults: jump = membrane jump, angle = 180 deg, axis =
	///			 axis between COMs projected into the membrane plane
	FlipMover();
	
//	/// @brief Custom Constructor
//	/// @details User can specify jump number
//	FlipMover( Size jump_num );
//
//	/// @brief Custom constructor
//	/// @details User can specify jump number and rotation axis
//	FlipMover( Size jump_num, Vector axis );
//
//	/// @brief Custom constructor
//	/// @details User can specify jump number and angle
//	FlipMover( Size jump_num, Real angle );
//	
//	/// @brief Custom constructor
//	/// @details User can specify jump number and rotation axis
//	FlipMover( Size jump_num, Vector axis, Real angle );
//
//	/// @brief Copy Constructor
//	FlipMover( FlipMover const & src );

	/// @brief Assignment Operator
	FlipMover & operator = ( FlipMover const & src );
	
	/// @brief Destructor
	virtual ~FlipMover();
	
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
	
	/// @brief Get the name of this Mover (FlipMover)
	virtual std::string get_name() const;
		
	/// @brief Flip the downstream partner in the membrane
	virtual void apply( Pose & pose );
	
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
	Size jump_num_;
	
	/// @brief Rotation axis
	Vector axis_;

	/// @brief Rotation angle in degrees
	Real angle_;
};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_FlipMover_hh
