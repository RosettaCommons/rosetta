// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/docking/membrane/MPFindInterfaceMoverCreator.hh
/// @brief      Sample protein-protein interface in the membrane
/// @details	The foldtree after the mover is reset to the original one as it
///				was before - so it should work with both fixed and movable membrane
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_docking_membrane_MPFindInterfaceMover_hh
#define INCLUDED_protocols_docking_membrane_MPFindInterfaceMover_hh

// Unit Headers
#include <protocols/docking/membrane/MPFindInterfaceMover.fwd.hh>
#include <protocols/docking/membrane/MPFindInterfaceMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/conformation/membrane/SpanningTopology.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <protocols/scoring/Interface.fwd.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh> 
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers
#include <utility/file/file_sys_util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.fwd.hh>

namespace protocols {
namespace docking {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::scoring;
using namespace protocols::moves;
	  
class MPFindInterfaceMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details Takes topology from MembraneInfo, gets the jump from docking partners
	MPFindInterfaceMover();

	/// @brief Copy Constructor
	MPFindInterfaceMover( MPFindInterfaceMover const & src );
	
	/// @brief Destructor
	virtual ~MPFindInterfaceMover();
	
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
	
	/// @brief Get the name of this Mover (MPFindInterfaceMover)
	virtual std::string get_name() const;
		
	/// @brief Sampling protein-protein interface in the membrane
	virtual void apply( Pose & pose );
	
private: // methods
	
	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options with JD2
	void register_options();
	
	/// @brief Initialize Mover options from the commandline
	void init_from_cmd();

	/// @brief Finalize setup
	void finalize_setup( Pose & pose );
	
	/// @brief Superimpose upstream partner
	/// @details Superimpose upstream partner of the pose with the native
	void superimpose_upstream_partner( Pose & pose );
	
	/// @brief Calculate interface SASA
	/// @details Calculate SASA buried in the interface
	core::Real calculate_interface_SASA( Pose & pose, Interface & interface );
	
	/// @brief Calculate fraction of small residues in interface
	/// details See Andrew Bordner, 2009, BMC Bioinformatics
	std::pair< core::Real, core::Real > fractions_small_residues( Pose & pose, Interface & interface);
	
private: // data

	/// @brief Native pose
	Pose native_;

	/// @brief docking partners
	std::string partners_;
	
	/// @brief jumps
	int jump_;
	utility::vector1< int > jumps_;
	
	/// @brief SpanningTopology objects
	SpanningTopologyOP topo_; 		// full pose
	SpanningTopologyOP topo_up_;	// upstream partner
	SpanningTopologyOP topo_down_; 	// downstream partner
	
	/// @brief scorefunction
	core::scoring::ScoreFunctionOP sfxn_;
	
	/// @brief allow flips in the membrane?
	bool flips_;
	
	/// @brief Run quick relax before and after protocol
	/// @details Kind of like flexible backbone and sidechain docking
	bool flexible_bb_;
	
};

} // membrane
} // docking
} // protocols

#endif // INCLUDED_protocols_docking_membrane_MPFindInterfaceMover_hh
