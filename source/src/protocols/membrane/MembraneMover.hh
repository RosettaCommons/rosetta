// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/MembraneMover.hh
///
/// @brief      Top-Level Unit Test for the Membrane Protein Factory (Mover Class)
/// @details    The purpose of this application is to test the membrane protein factory
///             initialization code including all external dependencies which cannot be tested
///             in JD2, Resource Manager, and the Pose cache. This can also serve as the integration
///             test for memrane protein initialization.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (1/3/14)

#ifndef INCLUDED_protocols_membrane_MembraneMover_hh
#define INCLUDED_protocols_membrane_MembraneMover_hh

// Unit Headers
#include <protocols/membrane/MembraneMover.fwd.hh>

// Project Headers
#include <core/membrane/MembraneProteinFactory.hh>
#include <core/membrane/util/Exceptions.hh>
#include <core/membrane/util/definitions.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Protocol Headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>

// Utility Headers
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>
#include <string>

namespace protocols {
namespace membrane {
    
    /// @brief Membrane Mover Class
    /// @details Helper Mover class for testing the membrane protein factory
    class MembraneMover : public protocols::moves::Mover {
    
    public: // constructors
        
        /// @brief Standard Constructor
        MembraneMover();
        
        /// @brief Copy Constructor
        MembraneMover( MembraneMover const & src );
        
        /// @brief Standard Destructor
        ~MembraneMover();
        
    public: // options
        
        /// @brief Register Relevant Options
        void register_options();
        
        /// @brief Initialize Options from the Command line
        void init_from_cmd();
        
    private: // non membrane pose construction methods
        
        /// @brief Initialize Chains
        /// @details Initialize Chains List from Resource Options
        std::map< core::Size, std::string > initialize_chains();
        
        /// @brief Make a Non membrane Pose
        /// @details Load in the non memrbane counterpart of a membrane pose
        core::pose::PoseOP initialize_non_membrane_pose( std::map< core::Size, std::string > chain_descriptions );
        
    public: // Methods
        
        /// @brief Return a Membrane pose
        core::pose::PoseOP get_membrane_pose();
        core::pose::PoseOP get_non_membrane_pose();
        
        // Mover Methods
        /// @brief Return the name of the Mover.
        virtual std::string get_name() const;
        
        /// @brief Clone Method
        virtual protocols::moves::MoverOP clone() const;
        
        /// @brief Fresh Instance
        virtual protocols::moves::MoverOP fresh_instance() const;
        
        /// @brief Apply the corresponding move to the pose
        virtual void apply( core::pose::Pose & pose);
        
    private: // data
        
        core::pose::PoseOP membrane_pose_;
        core::pose::PoseOP non_mp_pose_;
        
        // Options
        bool fullatom_;
        std::string chains_;
        
    }; // Membrane Mover class
    
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MembraneMover_hh

