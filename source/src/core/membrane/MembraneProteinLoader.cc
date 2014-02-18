// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/membrane/MembraneProteinLoader.cc
/// @brief  resource loader for membrane proteins
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneProteinLoader_cc
#define INCLUDED_core_membrane_MembraneProteinLoader_cc

// unit headers
#include <core/membrane/MembraneProteinLoader.hh>
#include <core/membrane/MembraneProteinLoaderCreator.hh>

// project headers
#include <core/membrane/MembraneProteinFactory.hh>
#include <core/membrane/MembraneProteinOptions.hh>

//utility headers
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

//C++ headers
#include <istream>

static basic::Tracer TR("core.membrane.MembraneProteinLoader");

namespace core {
namespace membrane {
    
    /// @brief Constructor
    MembraneProteinLoader::MembraneProteinLoader() {}
    
    /// @brief Destructor
    MembraneProteinLoader::~MembraneProteinLoader() {}
    
    /// @brief   Template Create Resource Method
    /// @details  Create membrane protein resource using RM resource loader method (doesn't use istream or locatorid)
    utility::pointer::ReferenceCountOP
    MembraneProteinLoader::create_resource(
                                     basic::resource_manager::ResourceOptions const & options,
                                     basic::resource_manager::LocatorID const &,
                                     std::istream &
                                     ) const
    {
        using namespace core::membrane;
        
        if ( ! dynamic_cast< MembraneProteinOptions const * > ( &options ) ) {
            throw utility::excn::EXCN_Msg_Exception( "MembraneProteinLoader expected to be given a MembraneProteinOptions object, " \
                                                    "but was given a non-MembraneProteinOptions object of type '" + options.type() + "', which has the name '" + options.name() + "'." );
        }
        //MembraneProteinOptions const & opts = static_cast< MembraneProteinOptions const & > ( options );
   
        // Create pose from membrane protein factory
//        MembraneProteinFactoryOP mpf = new MembraneProteinFactory( opts.include_lips(), opts.membrane_chains(), opts.fullatom() );
        MembraneProteinFactoryOP mpf = new MembraneProteinFactory( false, "protocols/membrane/chains.txt", true );
        core::pose::PoseOP pose = mpf->create_membrane_pose();
        TR << "I have created my membrane pose in loading" << std::endl;
        assert( pose );
        return pose;
    }
    
    /// @brief Return Default Options Object
    basic::resource_manager::ResourceOptionsOP
    MembraneProteinLoader::default_options() const
    {
        return new MembraneProteinOptions;
    }
    
    /// @brief Return an instance of teh loader to the loader creator
    basic::resource_manager::ResourceLoaderOP MembraneProteinLoaderCreator::create_resource_loader() const
    {
        return new MembraneProteinLoader();
    }
    
    /// @brief Return the Loader type (and tag for XML)
    std::string MembraneProteinLoaderCreator::loader_type() const
    {
        return "MembraneProtein";
    }
        
    
} // namespace membrane
} // namespace core

#endif // INCLUDED_core_membrane_MembraneProteinLoader_cc

