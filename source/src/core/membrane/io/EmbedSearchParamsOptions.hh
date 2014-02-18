// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedSearchParamsOptions.hh
///
/// @brief      Embedding Search Parameters Options class  - Contains options for membrane search and score
/// @details    Membrane proteins in rosetta use the membrane scoring function and an MCM embedidng search
///             which can be tuned and adjusted using the following options.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedSearchParamsOptions_hh
#define INCLUDED_core_membrane_io_EmbedSearchParamsOptions_hh

// Unit Headers
#include <core/membrane/io/EmbedSearchParamsOptions.fwd.hh>

// Project Headers
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/definitions_util.hh>
#include <core/membrane/util/Exceptions.hh>

#include <core/types.hh>

#include <basic/resource_manager/ResourceOptions.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace membrane {
namespace io {

/// @brief Embedding Search Param Options
/// @details Specify options for searching and scoring membrane protein embeddings
class EmbedSearchParamsOptions : public basic::resource_manager::ResourceOptions {

public:

    /// @brief Constructor
    EmbedSearchParamsOptions();
    
    /// @brief Destructor
    virtual ~EmbedSearchParamsOptions();

    /// @brief Parse .xml file for embedding opts
    virtual
    void
    parse_my_tag(
                 utility::tag::TagCOP tag
                 );

    /// @brief Return Options Class Type - Embedding Search Parameters
    virtual
    std::string
    type() const;

    /// @brief Setters and Getters for Options
    
    /// @brief Normal Search
    /// @details Search for membrane protein embedding normal
    bool normal_search() const;
    void set_normal_search( bool setting );

    /// @brief Normal Start Angle
    /// @details Search for membrane protein embedding start angle
    core::Real normal_start_angle() const;
    void set_normal_start_angle( core::Real setting );

    /// @brief Normal Max Angle
    /// @details Set max normal search angle (upper bound)
    core::Real normal_max_angle() const;
    void set_normal_max_angle( core::Real setting );

    /// @brief Normal Delta Angle
    /// @details Set allowed delta value for normal search (also bound on search)
    core::Real normal_delta_angle() const;
    void set_normal_delta_angle( core::Real setting );

    /// @brief Center Search
    /// @details Search for Center
    bool center_search() const;
    void set_center_search( bool setting );

    /// @brief Center Max Delta
    /// @details Center Max Delta for Search
    core::Real center_max_delta() const;
    void set_center_max_delta( core::Real setting );

    /// @brief Center Magnitude
    core::Size center_mag() const;
    void set_center_mag( core::Size setting );

    /// @brief Normal Mag
    core::Size normal_mag() const;
    void set_normal_mag( core::Size setting );

    /// @brief Normal Cycles
    /// @details Max Cycles allowed in normal search
    core::Size normal_cycles() const;
    void set_normal_cycles( core::Size setting );
    
    /// @brief Getter/Setter pair for penalties
    /// @details Param tyupe: bool
    bool penalties() const;
    void set_penalties( bool setting );

    /// @brief Getter/Setter pair for interpolating mpair term
    /// @details Param Type: bool
    bool no_interpolate_mpair() const;
    void set_no_interpolate_mpair( bool setting );
    
private: // data

    // Data
    bool normal_search_;
    core::Real normal_delta_angle_;
    core::Real normal_start_angle_;
    core::Real normal_max_angle_;

    bool center_search_;
    core::Real center_max_delta_;

    core::Size center_mag_;
    core::Size normal_mag_;
    core::Size normal_cycles_;
    
    bool penalties_;
    bool no_interpolate_mpair_;

}; // class EmbedSearchParamsOptions

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedSearchParamsOptions_hh


