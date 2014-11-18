// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	    protocols/membrane/MembranePositionFromTopologyMover.hh
///
/// @brief      Computes and sets the initial position of the membrane
/// @details	Computes and sets the initial position of the membrane from
///				sequence or structure (can be specified by the user at construction
///				or as a setup cmd flag).
///				CAUTION: ONLY FOR FLEXIBLE MEMBRANE AND FIXED PROTEIN!!!
///
///				NOTE: Requires a membrane pose!
///				NOTE: sequence not yet implemented
///				Last Modified: 6/21/14
///
/// @author		Rebecca Alford (rflaford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_MembranePositionFromTopologyMover_hh
#define INCLUDED_protocols_membrane_MembranePositionFromTopologyMover_hh

// Unit Headers
#include <protocols/membrane/MembranePositionFromTopologyMover.fwd.hh>

// Package headers
#include <core/pose/Pose.fwd.hh> 
#include <core/conformation/membrane/Span.hh>
#include <core/types.hh> 

// Project Headers
#include <protocols/moves/Mover.hh> 

namespace protocols {
namespace membrane {

using namespace core;
using namespace core::pose; 

/// @brief Compute the initial position of the membrane based upon sequence
/// or structure
class MembranePositionFromTopologyMover : public protocols::moves::Mover {
	
public:
	
	////////////////////
	/// Constructors ///
	////////////////////
	
	/// @brief Defualt Constructor
	/// @details Compute the embedding of the pose based on xyz coordinates
	/// and spanning topology provided in MembraneInfo
	MembranePositionFromTopologyMover();
	
	/// @brief Custom Constructor - for Pyrosetta
	/// @details Compute the embedding of the pose - if structure_based is
	/// true do this based on xyz coordinates. If structure_based is false,
	/// compute based on sequence.
	MembranePositionFromTopologyMover( bool structure_based );
	
	/// @brief Copy Constructor
	/// @details Make a deep copy of this mover
	MembranePositionFromTopologyMover( MembranePositionFromTopologyMover const & src );
	
	/// @brief Destructor
	~MembranePositionFromTopologyMover();
	
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
	
	/// @brief Update Membrane position in pose
	/// @details Compute membrane posiiton based on sequence or structure
	/// and then call pose.update_membrane_position() to update the membrane position
	virtual void apply( Pose & pose );
	
	/// @brief Get the name of this mover
	virtual std::string get_name() const;
		
private: // data

	// Structure or sequence based?
	bool structure_based_;
	
};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MembranePositionFromTopologyMover_hh
