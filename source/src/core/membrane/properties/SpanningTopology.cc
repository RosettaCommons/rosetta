// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/properties/SpanningTopology.cc
///
/// @brief      Membrane Spanning Topology Data
/// @details    Stores information describing the membrane spanning
///             topology of a pose. This definition is a dependency for embedding definitions
///             and requires a spanningfile from OCTOPUS for initialization
///
/// @note       Last Modified: 1/1/14
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_properties_SpanningTopology_cc
#define INCLUDED_core_membrane_properties_SpanningTopology_cc

// Unit headers
#include <core/membrane/properties/SpanningTopology.hh>

// Project Headers
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

// Package Headers
#include <core/types.hh>

#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>
#include <utility/vector0.hh>

#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

// Platform headers
#include <platform/types.hh>

// C++ Headers
#include <cstddef>
#include <cstdlib>
#include <string>

static basic::Tracer TR("core.membrane.properties.SpanningTopology");

using namespace core;

/// @brief      Class: Membrane Spanning Topology
/// @details    Stores information describing the membrane spanning
///             topology of a pose. This definition is a dependency for embedding definitions
///             and requires a spanningfile from OCTOPUS for initialization

namespace core {
namespace membrane {
namespace properties {
    
    ///// Constructors /////////////////////////////
    
    /// @brief Constructor
    SpanningTopology::SpanningTopology() :
        utility::pointer::ReferenceCount(),
        total_residue_in_span_file_(0),
        total_tmhelix_(0),
        tmh_inserted_(0)
    {
        // Initialize Remaining Members
        helix_id_.clear();
        span_.clear();
        full_span_.clear();
        relative_tmh_ori_.clear();
        tmregion_.clear();
    }
    
    /// @brief Destructor
    SpanningTopology::~SpanningTopology() {}
    
    ////// Setters and Getters ////////////////////////////////
    
    /// @brief Getters
    // Spanning ranges for transmembrane helices
    ObjexxFCL::FArray1D< Size > SpanningTopology::helix_id() { return helix_id_; }
    ObjexxFCL::FArray2D< Size > SpanningTopology::span() { return span_; }
    ObjexxFCL::FArray2D< Size > SpanningTopology::full_span() { return full_span_; }
    ObjexxFCL::FArray2D< Size > SpanningTopology::relative_tmh_ori() { return relative_tmh_ori_; }
    
    // Needed for evaluating tm penalties in search methods
    utility::vector1< bool > SpanningTopology::allow_scoring() { return allow_scoring_; }
    utility::vector1< bool > SpanningTopology::allow_tmh_scoring() { return allow_tmh_scoring_; }
    
    // Info about sequence directly from spanfile
    core::Size SpanningTopology::total_residue_in_span_file() { return total_residue_in_span_file_; }
    core::Size SpanningTopology::total_tmhelix() { return total_tmhelix_; }
    core::Size SpanningTopology::tmh_inserted() { return tmh_inserted_; }
    
    // Vector of 0/1 specifying if tm region
    utility::vector1< bool > SpanningTopology::tmregion() { return tmregion_; }
    
    /// @brief Setters
    // Spanning ranges for transmembrane helices
    void SpanningTopology::set_helix_id( ObjexxFCL::FArray1D< Size > helix_id ) { helix_id_ = helix_id; }
    void SpanningTopology::set_span( ObjexxFCL::FArray2D< Size > span ) { span_ = span; }
    void SpanningTopology::set_full_span( ObjexxFCL::FArray2D< Size > full_span ) { full_span_ = full_span; }
    void SpanningTopology::set_relative_tmh_ori( ObjexxFCL::FArray2D< Size > relative_tmh_ori ) { relative_tmh_ori_ = relative_tmh_ori; }
    
    // Needed for evaluating tm penalties in search methods
    void SpanningTopology::set_allow_scoring( utility::vector1< bool > allow_scoring ) { allow_scoring_ = allow_scoring; }
    void SpanningTopology::set_allow_tmh_scoring( utility::vector1< bool > allow_tmh_scoring ) { allow_tmh_scoring_ = allow_tmh_scoring; }
    
    // Vector of 0/1 specifying if tm region
    void SpanningTopology::set_total_residue_in_spanfile( core::Size total ) { total_residue_in_span_file_ = total; }
    void SpanningTopology::set_total_tmhelix( core::Size total ) { total_tmhelix_ = total; }
    void SpanningTopology::set_tmh_inserted( core::Size inserted ) { tmh_inserted_ = inserted; }
    
    // Vector of 0/1 specifying if tm region
    void SpanningTopology::set_tmregion( utility::vector1< bool > tmregion ) { tmregion_ = tmregion; }
    
    //// Public Utility Methods /////////////////////////////
    
    /// @brief Reset # TMH Inserted
    core::Size
    SpanningTopology::reset_tmh_insert() {
        tmh_inserted_ = 0;
        return tmh_inserted_;
    }
    
    /// @brief Reset allowed scoring positions
    void
    SpanningTopology::reset_allowed_scoring() {
        
        // Set tmh_inserted to 0
        tmh_inserted_ = 0;
        
        // Loop through allow_tmh_scoring vector and set all values to false
        for ( Size i = 1; i <= allow_tmh_scoring_.size(); ++i ) {
            allow_tmh_scoring_[i] = false;
        }
        
        // Loop through allow_scoring vector and set all values to false
        for ( Size i = 1; i <= allow_scoring_.size(); ++i ) {
            allow_scoring_[i] = false;
        }
    }
    
    /// @brief Shift membrane spanning by shift factor
    void
    SpanningTopology::shift_span( Size shift_factor ) {
        
        // Loop throough span arrays and shift helices
        for (Size i = 1; i <= total_tmhelix_; ++i ) {
            
            // shift individual spans
            span_(i, 1) += shift_factor;
            span_(i, 2) += shift_factor;
            
            // Shift full span
            full_span_(i, 1) += shift_factor;
            full_span_(i, 2) += shift_factor;
            
            if ( full_span_(i, 1) <= 0 ) {
                utility_exit_with_message("Topology spanning out of bounds");
            }
        }
        
        // Done!
        return;
    }
    
    /// @brief Print tmh spanning info
    void
    SpanningTopology::show() {
        
        // Print Stuff
        TR << "Total Transmembrane Helices" << total_tmhelix_ << std::endl;
        for ( Size i = 1; i <= total_tmhelix_; ++i ) {
            TR << "SPAN" << i << " " << span_(i, 1) << " " << span_(i, 2) << std::endl;
        }
        
        // Done!
        return;
    }
    
    /// @brief Get a subset of helices from the topology object
    void
    SpanningTopology::get_subset( utility::vector1< Size > & TMH_list ) {
        
        TR << "Grabbing a subset of tm helices from the topology object" << std::endl;
        
        // Ensure the TMH list is in bounds of the number of helices
        if ( TMH_list.size() > total_tmhelix_ ) {
            utility_exit_with_message( "Transmembrane Helix list is too long to grab a subset" );
        }
        
        // Set Length
        Size const len( TMH_list.size() );
        
        // Set dimensions
        span_.dimension(len, 2);
        full_span_.dimension(len, 2);
        relative_tmh_ori_.dimension(len, len);
        
        // Read through the spans and copy to topology object
        for ( Size i = 1; i <= len; ++i ) {
            
            // Copy normal spanning
            span_(i, 1) = span_( TMH_list[i], 1);
            span_(i, 2) = span_( TMH_list[i], 2);
            
            /// Copy full spanning
            full_span_(i, 1) = full_span_(TMH_list[i], 1);
            full_span_(i, 2) = full_span_(TMH_list[i], 2);
            
            // Set Helix ID equal to list index
            helix_id_(i) = TMH_list[i];
            
            // Read through tmh relative list and copy
            for ( Size j = 1; j <= TMH_list.size(); ++j ) {
                relative_tmh_ori_(i, j) = relative_tmh_ori_( TMH_list[i], TMH_list[j] );
            }
        }
    }
    
    /// @brief Copy Data
    /// @details Copy Constructor Helper function
    void
    SpanningTopology::copy_data( SpanningTopology src, SpanningTopology copy ) {
        
        // copy data...
        src.helix_id_ = copy.helix_id_;
        src.span_ = copy.span_;
        src.full_span_ = copy.full_span_;
        src.relative_tmh_ori_ = copy.relative_tmh_ori_;
        
        src.allow_tmh_scoring_ = copy.allow_tmh_scoring_;
        src.allow_scoring_ = copy.allow_scoring_;
        
        src.total_residue_in_span_file_ = copy.total_residue_in_span_file_;
        src.total_tmhelix_ = copy.total_tmhelix_;
        src.tmh_inserted_ = copy.tmh_inserted_;
        
        src.tmregion_ = copy.tmregion_;
    }
    
} // properties
} // membrane
} // core

#endif // INCLUDED_core_membrane_properties_SpanningTopology_cc

