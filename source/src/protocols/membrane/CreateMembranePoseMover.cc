// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/CreateMembranePoseMover.cc
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

#ifndef INCLUDED_protocols_membrane_CreateMembranePoseMover_cc
#define INCLUDED_protocols_membrane_CreateMembranePoseMover_cc

// Unit Headers
#include <protocols/membrane/CreateMembranePoseMover.hh>

// Project Headers
#include <core/membrane/MembraneProteinFactory.hh>
#include <core/membrane/util/Exceptions.hh>
#include <core/membrane/util/definitions.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>

// Protocol Headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <cstdlib>
#include <string>

using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.membrane.CreateMembranePoseMover");

namespace protocols {
namespace membrane {
        
    /// Constructors /////////////////////////////
    
    /// @brief Standard Constructor
    CreateMembranePoseMover::CreateMembranePoseMover() :
        protocols::moves::Mover(),
        membrane_pose_(NULL),
        fullatom_(true),
        chains_("")
    {}
    
    /// @brief Copy Constructor
    CreateMembranePoseMover::CreateMembranePoseMover( CreateMembranePoseMover const & src ) :
        protocols::moves::Mover(src),
        membrane_pose_(NULL),
        fullatom_(true),
        chains_("")
    {
        membrane_pose_ = src.membrane_pose_;
    }
    
    /// @brief Standard Destructor
    CreateMembranePoseMover::~CreateMembranePoseMover() {}
    
    //// Options System //////////
    
    /// @brief Register Relevant Options
    void CreateMembranePoseMover::register_options() {
        
        using namespace basic::options;
        
        option.add_relevant( OptionKeys::in::file::membrane_chains );
        option.add_relevant( OptionKeys::in::file::fullatom );
    }
    
    /// @brief Initialize Options from the Command line
    void CreateMembranePoseMover::init_from_cmd() {
        
        using namespace basic::options;
        
        // Check for chains specified
        if ( option[ OptionKeys::in::file::membrane_chains].user() ) {
            chains_ = option[ OptionKeys::in::file::membrane_chains ]();
        }
        
        // Check for fullatom option
        if ( option[ OptionKeys::in::file::fullatom ].user() ) {
            fullatom_ = option[ OptionKeys::in::file::fullatom ]();
        }
    }
    
    /// Mover Inherited Methods ////////////////////////
    
    /// @brief Return your poses
    core::pose::PoseOP CreateMembranePoseMover::get_membrane_pose() { return membrane_pose_; }
    
    // Mover Methods
    /// @brief Return the name of the Mover.
    std::string CreateMembranePoseMover::get_name() const { return "CreateMembranePoseMover"; }
    
    /// @brief Clone Method
    protocols::moves::MoverOP CreateMembranePoseMover::clone() const { return new CreateMembranePoseMover(*this); }
    
    /// @brief Fresh Instance
    protocols::moves::MoverOP CreateMembranePoseMover::fresh_instance() const { return new CreateMembranePoseMover(); }

    /// Mover Apply
    
    /// @brief Apply the corresponding move to the pose
    void CreateMembranePoseMover::apply( core::pose::Pose & )
    {
        using namespace core::membrane;
        using namespace basic::resource_manager;
        using namespace core::membrane::util;
        
        TR << "Loading in a membrane pose" << std::endl;
        MembraneProteinFactoryOP mpf = new MembraneProteinFactory( false, chains_, fullatom_ );
        membrane_pose_ = mpf->create_membrane_pose();

    }
        
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_CreateMembranePoseMover_cc

