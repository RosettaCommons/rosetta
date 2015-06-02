// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/AddMembraneMover.hh
///
/// @brief      Add Membrane Representation to the Pose
/// @details	Given a pose, setup membrane topology, lips info,
///				and a membrane virtual residue in the pose. All of this information
///				is coordinated via the MembraneInfo object maintained in
///				the Pose's conformation. After applying AddMembraneMover
///				to the pose, pose.conformation().is_membrane() should always
///				return true.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (7/10/14)

#ifndef INCLUDED_protocols_membrane_AddMembraneMover_hh
#define INCLUDED_protocols_membrane_AddMembraneMover_hh

// Unit Headers
#include <protocols/membrane/AddMembraneMover.fwd.hh>
#include <protocols/moves/Mover.hh> 

// Project Headers
#include <core/conformation/membrane/SpanningTopology.fwd.hh>
#include <core/conformation/membrane/Span.fwd.hh>
#include <core/conformation/membrane/LipidAccInfo.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh> 
#include <core/types.hh> 

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

// C++ Headers
#include <cstdlib>

namespace protocols {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::moves;
	  
/// @brief	Add Membrane components to the pose, sets up MembraneInfo object
///			includes:
///			spanning topology, lips info, embeddings, and a membrane
///			virtual residue describing the membrane position
class AddMembraneMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details Create a membrane pose setting the membrane center
	/// at center=(0, 0, 0), normal=(0, 0, 1) and loads in spans
	/// and lips from the command line interface.
	AddMembraneMover();
	
	/// @brief Custom Constructor - Create membrane pose from existing membrane resnum
	/// @details Loads in membrane information from the commandline. Just needs a
	/// membrane residue number (application is symmetry)
	AddMembraneMover( core::Size membrane_rsd );
	
	/// @brief Custom Constructor - mainly for PyRosetta
	/// @details Creates a membrane pose setting the membrane
	/// center at emb_center and normal at emb_normal and will load
	/// in spanning regions from list of spanfiles provided
	AddMembraneMover(
		std::string spanfile,
		core::Size membrane_rsd=0
	);
		
    /// @brief Custom Constructor - mainly for PyRosetta
    /// @details Creates a membrane pose setting the membrane
    /// center at emb_center and normal at emb_normal and will load
    /// in spanning regions from list of spanfiles provided
    AddMembraneMover(
        SpanningTopologyOP topology,
        core::Size anchor_rsd=1,
        core::Size membrane_rsd=0
        );
    
	/// @brief Custom Constructor using both anchor rsd and membrane rsd
	/// @details Creates a membrane pose setting the anchor residue and membrane residue; topology is read from cmd line
	AddMembraneMover(
		core::Size anchor_rsd,
		core::Size membrane_rsd
		);
	
    /// @brief Custorm Constructur with lips info - mainly for PyRosetta
	/// @details Creates a membrane pose setting the membrane
	/// center at emb_center and normal at emb_normal and will load
	/// in spanning regions from list of spanfiles provided. Will also
	/// load in lips info from lips_acc info provided
	AddMembraneMover(
		std::string spanfile,
		std::string lipsfile,
		core::Size membrane_rsd=0
	);
    
    /// @brief Custorm Constructur - Center/Normal init
    /// @details Creates a membrane protein, setting the initial
    /// membrane center and normal to the specified values, loads
    /// a spanfile from given path, and indicates a membrane
    /// residue if provided.
    AddMembraneMover(
         Vector init_center,
         Vector init_normal,
         std::string spanfile = "",
         core::Size membrane_rsd = 0
         );
	
	/// @brief Copy Constructor
	/// @details Create a deep copy of this mover
	AddMembraneMover( AddMembraneMover const & src );
	
	/// @brief Destructor
	virtual ~AddMembraneMover();
	
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

    ////////////////////////////////////////////////////
    /// Getters - For subclasses of AddMembraneMover ///
    ////////////////////////////////////////////////////

    /// @brief Return the current path to the spanfile held
    /// by this mover
    std::string get_spanfile() const;
    
    ////////////////////////////////////////////////////
    /// Setters - For subclasses of AddMembraneMover ///
    ////////////////////////////////////////////////////
    
    /// @brief Set Spanfile path
    /// @details Set the path to the spanfile
    void spanfile( std::string spanfile );

    /// @brief Set lipsfile path
    /// @details Set the path to the lipsfile
    void lipsfile( std::string lipsfile );

    /// @brief Set option for including lipophilicity data
    /// @details Incidate whether lipophilicity information should be read
    /// and used in MembraneInfo
    void include_lips( bool include_lips );
	
	/////////////////////
	/// Mover Methods ///
	/////////////////////
	
	/// @brief Get the name of this Mover (AddMembraneMover)
	virtual std::string get_name() const;
		
	/// @brief Add Membrane Components to Pose
	/// @details Add membrane components to pose which includes
	///	spanning topology, lips info, embeddings, and a membrane
	/// virtual residue describing the membrane position
	virtual void apply( Pose & pose );
	
	/// @brief Helper Method - Setup Membrane Virtual
	/// @details Create a new virtual residue of type MEM from
	/// the pose typeset (fullatom or centroid). Add this virtual
	/// residue by appending to the last residue of the pose. Then set
	/// this position as the root of the fold tree.
	virtual core::Size setup_membrane_virtual( Pose & pose );
    
    /// @brief Helper Method - Check for Membrane residue already in the PDB
    /// @details If there is an MEM residue in the PDB at the end of the pose
    /// with property MEMBRANE, return it's residue number. In the control flow of the
    /// apply method, if this returns non-zero, a new membrane residue will not be added.
    virtual utility::vector1< core::SSize > check_pdb_for_mem( Pose & pose );
    
private:
	
	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options from Command Line
	/// @details Register mover-relevant options with JD2 - includes
	/// mp, seutp options: center, normal, spanfiles and
	/// lipsfiles
	void register_options();
	
	/// @brief Initialize Mover options from the comandline
	/// @details Initialize mover settings from the commandline
	/// mainly in the mp, setup group: center, normal,
	/// spanfiles and lipsfiles paths
	void init_from_cmd();

private:

	// Pose residye typeset & include lips
	bool include_lips_;
	
	// SpanningTopology
	std::string spanfile_;
	SpanningTopologyOP topology_;

	// Lipid Accessibility Info - Lips Files
	std::string lipsfile_;
	
    // the membrane residue is anchored to this
    // residue position by jump
    core::Size anchor_rsd_;
    
	// Membrane residue number
	core::Size membrane_rsd_;
    
    // Initial Center/Normal Pair
    core::Vector center_;
    core::Vector normal_;
    
    // Set membrane thickness & steepness (FOR HIGHRES ONLY!!!)
    core::Real thickness_;
    core::Real steepness_;
    
};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_AddMembraneMover_hh
