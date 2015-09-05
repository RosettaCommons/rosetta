// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/membrane/TransformIntoMembraneMover.hh
/// @brief  Transform a pose into a membrane coordinate frame
/// @author     JKLeman (julia.koehler1982@gmail.com)
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// Last Modified: 6/11/15
/// #RosettaMPMover

#ifndef INCLUDED_protocols_membrane_TransformIntoMembraneMover_hh
#define INCLUDED_protocols_membrane_TransformIntoMembraneMover_hh

// Unit Headers
#include <protocols/membrane/TransformIntoMembraneMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/membrane/geometry/EmbeddingDef.fwd.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <basic/Tracer.fwd.hh>

namespace protocols {
namespace membrane {

/// @brief Transform a pose into membrane coordinates based on the current
///   embedding of the protein
class TransformIntoMembraneMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Transform the protein into a membrane defined by MEM
	/// @details Transform the protein into membrane defined by MEM, protein
	/// embedding computed from structure and spanfile
	TransformIntoMembraneMover();

	// TODO: use and test this constructor
	// Use custom jump to transform protein into membrane
	// Using user-specified jump, transform the downstream partner
	// into default membrane, partner embedding computed from structure & spanfile
	TransformIntoMembraneMover( core::Size jump );

	/// @brief Transform the protein with a user-specified protein embedding into
	/// a default membrane (defined by MEM)
	/// @details Transform the protein with a user-defined embedding (might have
	/// been optimized before) into the default membrane
	TransformIntoMembraneMover( 
		protocols::membrane::geometry::EmbeddingDefOP current_embedding 
		);

	/// @brief Transform the protein into user-specified membrane coordinates
	/// @details Transform the protein into a user-defined membrane, protein
	/// embedding is computed from structure and spanfile
	TransformIntoMembraneMover( 
		core::Vector new_mem_cntr, 
		core::Vector new_mem_norm 
		);

	/// @brief Transform the protein into user-specified membrane coordinates
	/// @details Transform the protein into a user-defined membrane, protein
	/// embedding is computed from structure and spanfile
	TransformIntoMembraneMover( 
		protocols::membrane::geometry::EmbeddingDefOP current_embedding, 
		core::Vector new_mem_cntr, 
		core::Vector new_mem_norm 
		);

	/// @brief Copy Constructor
	TransformIntoMembraneMover( TransformIntoMembraneMover const & src );

	/// @brief Destructor
	virtual ~TransformIntoMembraneMover();

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

	/// @brief Use the default membrane (cntr 0,0,0 and normal 0,0,1) instead
	///   of the membrane from the MEM coordinates stored in MembraneInfo
	void use_default_membrane( bool truefalse );

	/// @brief Get the name of this Mover (TransformIntoMembraneMover)
	virtual std::string get_name() const;

	/// @brief Move the pose into membrane coordinate frame
	virtual void apply( core::pose::Pose & pose );

private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options from Command Line
	/// @details Register mover-relevant options with JD2 - includes
	/// mp, setup options: center, normal, spanfiles, jumpnum
	void register_options();

	/// @brief Initialize Mover options from the commandline
	/// @details Initialize mover settings from the commandline
	/// mainly in the mp, setup group: center, normal,
	/// spanfiles, jumpnum
	void init_from_cmd();


private: // data

	// Jump used for the move
	core::Size jump_;

	// new membrane coordinates to transform into
	core::Vector new_mem_cntr_;
	core::Vector new_mem_norm_;

	// Embedding of the protein prior to transformation
	protocols::membrane::geometry::EmbeddingDefOP current_embedding_;

	// use default membrane of (center 0,0,0 and normal 0,0,1)
	bool use_default_membrane_;

	// user-defined membrane
	bool user_defined_membrane_;

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_TransformIntoMembraneMover_hh
