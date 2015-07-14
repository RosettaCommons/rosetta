// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///	@file		protocols/membrane/TransformIntoMembraneMover.hh
/// @brief		Transform a pose into a membrane coordinate frame
/// @author		Rebecca Faye Alford (rfalford12@gmail.com)
/// @author     JKLeman (julia.koehler1982@gmail.com)
///				CAUTION: THIS MOVER ONLY WORKS FOR A FIXED MEMBRANE WHERE THE
///				MEMBRANE VIRTUAL RESIDUE IS AT THE ROOT OF THE FOLDTREE!!!
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

using namespace core;
using namespace core::pose;
using namespace protocols::membrane::geometry;
using namespace protocols::moves;
	  
/// @brief	Transform a pose into membrane coordinates based on the current
///			embedding of the protein
class TransformIntoMembraneMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Transform the protein into a defalt membrane
	/// @details Transform the protein into default membrane, current protein
    /// embedding computed from structure and spanfile
	TransformIntoMembraneMover();
    
    /// @brief Use custom jump to transform protein into membrane
    /// @details Using user-specified jump, transform the protein into default
    /// membrane, current protein computed from structure & spanfile
    TransformIntoMembraneMover( core::Size jump );
	
	/// @brief Transform the protein into a defalt membrane and user specified embedding
	/// @details Transform the protein with a user-defined embedding (might have
    /// been optimized before) into the default membrane
	TransformIntoMembraneMover( EmbeddingDefOP embedding );

	/// @brief Transform the protein into user-specified mmebrane coordinates
    /// @details Transform the protein with a user-defined embedding (might have
    /// been optimized before and is saved as MEM) into a user-defined membrane
    /// which is ARGV for the constructor
	TransformIntoMembraneMover(
		Vector center,
		Vector normal,
		bool from_current=false
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
	
	/// @brief Get the name of this Mover (TransformIntoMembraneMover)
	virtual std::string get_name() const;
		
	/// @brief Move the pose into membrane coordinate frame
	virtual void apply( Pose & pose );
	
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
    Size jump_;

	// center and normal coordinates with respect to fixed membrane
	Vector mem_center_;
	Vector mem_normal_;
	
	// Embedding of the protien prior to transformation
	EmbeddingDefOP embedding_;
	
	// Extract current membrane coordinates (prior to transformation)
	// from the membrane residue
	bool keep_current_protein_embedding_; 
	
};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_TransformIntoMembraneMover_hh
