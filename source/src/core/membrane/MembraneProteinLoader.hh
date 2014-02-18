// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/membrane/MembraneProteinLoader.hh
/// @brief  resource loader for membrane proteins
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneProteinLoader_hh
#define INCLUDED_core_membrane_MembraneProteinLoader_hh

//unit headers
#include <core/membrane/MembraneProteinLoader.fwd.hh>

//project headers
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <istream>

namespace core {
namespace membrane {

    class MembraneProteinLoader : public basic::resource_manager::ResourceLoader
    {
    public:
        
        /// @brief Constructor
        MembraneProteinLoader();
        
        /// @brief Destructor
        virtual ~MembraneProteinLoader();
        
        /// @brief Tempalte create resource method for the resource manager
        virtual
        utility::pointer::ReferenceCountOP
        create_resource(
                        basic::resource_manager::ResourceOptions const & options,
                        basic::resource_manager::LocatorID const & locator_id,
                        std::istream & istream
                        ) const;
        
        /// @brief Return default resource options object
        virtual
        basic::resource_manager::ResourceOptionsOP
        default_options() const;
        
    };

        
} // namespace membrane
} // namespace core

#endif // INCLUDED_core_membrane_MembraneProteinLoader_hh

