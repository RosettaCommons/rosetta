// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/membrane/MembraneProteinOptions.cc
/// @brief  options header for membrane protein options
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneProteinOptions_cc
#define INCLUDED_core_membrane_MembraneProteinOptions_cc

//unit headers
#include <core/membrane/MembraneProteinOptions.hh>
#include <core/membrane/MembraneProteinOptionsCreator.hh>

//utility headers
#include <utility/tag/Tag.hh>

namespace core {
namespace membrane {
    
    /// @brief Return options type (string)
    std::string
    MembraneProteinOptionsCreator::options_type() const { return "MembraneProteinOptions"; }
    
    /// @brief Create a resource options class from parsed xml input
    basic::resource_manager::ResourceOptionsOP
    MembraneProteinOptionsCreator::create_options() const { return new MembraneProteinOptions; }
    
    /// @brief Constructor
    MembraneProteinOptions::MembraneProteinOptions() :
        basic::resource_manager::ResourceOptions(),
        fullatom_( true ),
        include_lips_( false ),
        membrane_chains_( "" )
    {}
    
    /// @brief Destructor
    MembraneProteinOptions::~MembraneProteinOptions() {}
    
    /// @brief Parse XML File for options
    void
    MembraneProteinOptions::parse_my_tag(
                                   utility::tag::TagCOP tag
                                   )
    {
        fullatom( tag->getOption< bool >( "fullatom", true ));
        include_lips( tag->getOption< bool >( "include_lips", true ));
        membrane_chains( tag->getOption< std::string >( "membrane_chains" ));
    }
    
    /// @brief Return membrane protein options class type
    std::string
    MembraneProteinOptions::type() const
    {
        return "MembraneProteinOptions";
    }
    
    /// @brief Specify a fullatom pose constructon to the membrane protein factory
    bool MembraneProteinOptions::fullatom() const { return fullatom_; }
    void MembraneProteinOptions::fullatom( bool setting ) { fullatom_ = setting; }
    
    /// @brief Specify including lipid accessibility data in membrane pose construction + scoring
    bool MembraneProteinOptions::include_lips() const { return include_lips_; }
    void MembraneProteinOptions::include_lips( bool setting ) { include_lips_ = setting; }
    
    /// @brief Specify a path to file containing membrane chain refs
    std::string MembraneProteinOptions::membrane_chains() const { return membrane_chains_; }
    void MembraneProteinOptions::membrane_chains( std::string setting ) { membrane_chains_ = setting; }

} // namespace membrane
} // namespace core


#endif // INCLUDED_core_membrane_MembraneProteinOptions_cc

