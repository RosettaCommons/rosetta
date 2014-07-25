// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/MPDockingMover.fwd.hh
/// @brief      Dock two membrane proteins
/// @author     JKLeman (julia.koehler1982@gmail.com)
/// @note       Last Modified (6/24/14)

#ifndef INCLUDED_protocols_docking_membrane_MPDockingMover_hh
#define INCLUDED_protocols_docking_membrane_MPDockingMover_hh

// Unit Headers
#include <protocols/docking/membrane/MPDockingMover.fwd.hh>
#include <protocols/moves/Mover.hh> 

// Project Headers
#include <protocols/membrane/AddMembraneMover.fwd.hh>
#include <protocols/docking/DockMCMProtocol.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/PyMolMover.fwd.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh> 
#include <core/types.hh> 

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

namespace protocols {
namespace docking {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace protocols::moves;
	  
/// @brief Add Membrane components to the pose - includes
///	spanning topology, lips info, embeddings, and a membrane
/// virtual residue describing the membrane position
class MPDockingMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details Create a membrane pose setting the membrane center
	/// at center=(0, 0, 0), normal=(0, 0, 1) and loads in spans
	/// and lips from the command line interface.
	MPDockingMover();
		
	/// @brief Copy Constructor
	/// @details Create a deep copy of this mover
	MPDockingMover( MPDockingMover const & src );
	
	/// @brief Destructor
	virtual ~MPDockingMover();
	
public: // methods

	/// @brief Create a Clone of this mover
	virtual protocols::moves::MoverOP clone() const;
	
	/// @brief Create a Fresh Instance of this Mover
	virtual protocols::moves::MoverOP fresh_instance() const;
	
	/////////////////////
	/// Mover Methods ///
	/////////////////////
	
	/// @brief Get the name of this Mover (MPDockingMover)
	virtual std::string get_name() const;
		
	/// @brief Add Membrane Components to Pose
	/// @details Add membrane components to pose which includes
	///	spanning topology, lips info, embeddings, and a membrane
	/// virtual residue describing the membrane position
	virtual void apply( Pose & pose );
	
private: // methods
	
	// setup
	void setup();

private: // data

	// add membrane mover
	protocols::membrane::AddMembraneMoverOP add_membrane_mover_;
	
	// dock MCM protocol
	protocols::docking::DockMCMProtocolOP dock_mcm_protocol_;
	
	// sequence mover
	protocols::moves::RandomMoverOP random_mover_;
	
	// scorefunction
	core::scoring::ScoreFunctionOP scorefunction_;
	
	// kT
	Real kT_;

	// Membrane Center/Normal pair used for setup
	Vector center_;
	Vector normal_;

	// SpanningTopology
	std::string spanfile_;
	
};

} // membrane
} // docking
} // protocols

#endif // INCLUDED_protocols_docking_membrane_MPDockingMover_hh
