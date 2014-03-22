// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/conformation/membrane/SpanningTopology.hh
///
/// @brief      Membrane Spanning Topology Data
/// @details    Stores information describing the membrane spanning
///             topology of a pose. This definition is a dependency for embedding definitions
///             and requires a spanningfile from OCTOPUS for initialization
///
/// @note       Last Modified: 1/1/14
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_SpanningTopology_hh
#define INCLUDED_core_conformation_membrane_SpanningTopology_hh

// Unit headers
#include <core/conformation/membrane/SpanningTopology.fwd.hh>

// Package Headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

#include <utility/vector1.hh>
#include <utility/vector0.hh>

// Platform headers
#include <platform/types.hh>

// C++ Headers
#include <cstddef>
#include <cstdlib>
#include <string>

using namespace core;

namespace core {
namespace conformation {
namespace membrane {
    
    /// @brief      Class: Membrane Spanning Topology
    /// @details    Stores information describing the membrane spanning
    ///             topology of a pose. This definition is a dependency for embedding definitions
    ///             and requires a spanningfile from OCTOPUS for initialization
    class SpanningTopology : public utility::pointer::ReferenceCount {
    
    public: // constructors
        
        /// @brief Constructor
        SpanningTopology();
        
        /// @brief Destructor
        ~SpanningTopology();
        
    public: // getter/setter
        
        /// @brief Getters
        // Spanning ranges for transmembrane helices
        ObjexxFCL::FArray1D< Size > helix_id();
        ObjexxFCL::FArray2D< Size > span();
        ObjexxFCL::FArray2D< Size > full_span();
        ObjexxFCL::FArray2D< Size > relative_tmh_ori();
        
        // Needed for evaluating tm penalties in search methods
        utility::vector1< bool > allow_scoring();
        utility::vector1< bool > allow_tmh_scoring();
        
        // Info about sequence directly from spanfile
        core::Size total_residue_in_span_file();
        core::Size total_tmhelix();
        core::Size tmh_inserted();
        
        // Vector of 0/1 specifying if tm region
        utility::vector1< bool > tmregion();
        
        /// @brief Setters
        // Spanning ranges for transmembrane helices
        void set_helix_id( ObjexxFCL::FArray1D< Size > helix_id );
        void set_span( ObjexxFCL::FArray2D< Size > span );
        void set_full_span( ObjexxFCL::FArray2D< Size > full_span );
        void set_relative_tmh_ori( ObjexxFCL::FArray2D< Size > relative_tmh_ori );
        
        // Needed for evaluating tm penalties in search methods
        void set_allow_scoring( utility::vector1< bool > allow_scoring );
        void set_allow_tmh_scoring( utility::vector1< bool > allow_tmh_scoring );
        
        // Vector of 0/1 specifying if tm region
        void set_total_residue_in_spanfile( core::Size total );
        void set_total_tmhelix( core::Size total );
        void set_tmh_inserted( core::Size inserted );
        
        // Vector of 0/1 specifying if tm region
        void set_tmregion( utility::vector1< bool > );
        
    public: // methods
        
        /// @brief Reset # TMH Inserted
        core::Size reset_tmh_insert();
        
        /// @brief Reset allowed scoring positions
        void reset_allowed_scoring();
        
        /// @brief Shift membrane spanning by shift factor
        void shift_span( Size shift_factor );
        
        /// @brief Print tmh spanning info
        void show();
        
        /// @brief Get a subset of helices from the topology object
        void get_subset( utility::vector1< Size > & TMH_list );
        
    private: // helper methods
        
        /// @brief Copy Datag
        void copy_data( SpanningTopology src, SpanningTopology copy );
        
    private: // data
        
        // Spanning ranges for transmembrane helices
        ObjexxFCL::FArray1D< Size > helix_id_;
        ObjexxFCL::FArray2D< Size > span_;
        ObjexxFCL::FArray2D< Size > full_span_;
        ObjexxFCL::FArray2D< Size > relative_tmh_ori_;
        
        // Needed for evaluating tm penalties in search methods
        utility::vector1< bool > allow_scoring_;
        utility::vector1< bool > allow_tmh_scoring_;
        
        // Info about sequence directly from spanfile
        Size total_residue_in_span_file_;
        Size total_tmhelix_;
        Size tmh_inserted_;
        
        // Vector of 0/1 specifying if tm region
        utility::vector1< bool > tmregion_;

    }; // class SpanningTopology
    
} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_SpanningTopology_hh



