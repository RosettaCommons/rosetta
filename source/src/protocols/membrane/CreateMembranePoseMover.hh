// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/CreateMembranePoseMover.hh
///
/// @brief      Create Membrane Pose - Mover Class
/// @details    This mover creates a multi-chain membrane pose from the membrane framework
///             using a series of JD2 resource manager initialized resources. The objective
///             of the mover is to call the membrane protein factory, as for a new membrane
///             pose regardless of what has been passed on the commandline, and override that
///             pose. This should also be the top level interface to the membrane protein
///             framework.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (2/22/14)

#ifndef INCLUDED_protocols_membrane_CreateMembranePoseMover_hh
#define INCLUDED_protocols_membrane_CreateMembranePoseMover_hh

// Unit Headers
#include <protocols/membrane/CreateMembranePoseMover.fwd.hh>
#include <protocols/membrane/CreateMembranePoseMoverCreator.hh>

// Project Headers
#include <core/membrane/MembraneProteinFactory.hh>
#include <core/conformation/membrane/Exceptions.hh>
#include <core/conformation/membrane/definitions.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Protocol Headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/filters/Filter.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/rosetta_scripts/util.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>

// Utility Headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <cstdlib>
#include <string>

namespace protocols {
namespace membrane {
        
    /// @brief  Create Membrane Pose Mover Class
    /// @details Mover class for initializing a pose as a membrane protein
    ///          in Rosetta.
    class CreateMembranePoseMover : public protocols::moves::Mover {
        
    public: // constructors
        
        /// @brief Standard Constructor
        CreateMembranePoseMover();
        
        /// @brief Copy Constructor
        CreateMembranePoseMover( CreateMembranePoseMover const & src );
        
        /// @brief Standard Destructor
        ~CreateMembranePoseMover();
        
    public: // options system
        
        /// @brief Register Relevant Options
        void register_options();
        
        /// @brief Initialize Options from the Command line
        void init_from_cmd();
        
    public: // Mover Methods
        
        /// @brief Returns the membrane pose created by the mover
        core::pose::PoseOP get_membrane_pose();
        
        // Mover Methods
        /// @brief Return the name of the Mover.
        virtual std::string get_name() const;
        
        /// @brief Clone Method
        virtual protocols::moves::MoverOP clone() const;
        
        /// @brief Fresh Instance
        virtual protocols::moves::MoverOP fresh_instance() const;
        
        /// @brief Apply the corresponding move to the pose
        virtual void apply( core::pose::Pose & pose );
			
				/// @brief RosettaScripts hook for XML scripting
				void parse_my_tag( utility::tag::TagCOP tag,
													 basic::datacache::DataMap &,
													 protocols::filters::Filters_map const &,
													 protocols::moves::Movers_map const &,
													 core::pose::Pose const & );
		
    private: // data
        
        core::pose::PoseOP membrane_pose_;
        
        // Options
        bool fullatom_;
        std::string chains_;
        
    }; // Create Membrane Pose Mover class
        
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_CreateMembranePoseMover_hh

