// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedSearchParamsIO.cc
///
/// @brief      Embedding Search Parameters IO class  - Contains options for membrane search and score
/// @details    Membrane proteins in rosetta use the membrane scoring function and an MCM embedidng search
///             which can be tuned and adjusted using the following options.
///
/// @note       Resource Manager Component
/// @note       Container class - reads straight from options class
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedSearchParamsIO_cc
#define INCLUDED_core_membrane_io_EmbedSearchParamsIO_cc

// Unit Headers
#include <core/membrane/io/EmbedSearchParamsIO.hh>

// Project Headers
#include <core/conformation/membrane/definitions.hh>
#include <core/conformation/membrane/definitions_util.hh>
#include <core/conformation/membrane/Exceptions.hh>

#include <core/membrane/io/EmbedSearchParamsOptions.hh>

// Package Headers
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <algorithm>
#include <string>
#include <cstdlib>
#include <cmath>

using basic::Error;
using basic::Warning;

using std::stringstream;
using std::string;

static basic::Tracer TR( "core.membrane.io.EmbedSearchParamsIO" );

using namespace core::membrane;
using namespace core::conformation::membrane;

namespace core {
namespace membrane {
namespace io {

            /// @brief Constructor
            EmbedSearchParamsIO::EmbedSearchParamsIO() :
            utility::pointer::ReferenceCount()
            {}

            /// @brief Destructor
            EmbedSearchParamsIO::~EmbedSearchParamsIO()
            {}

            /// @brief Main IO Function for Reading in Search Parameters from an Options Obj
            EmbedSearchParamsOP
            EmbedSearchParamsIO::get_embed_params_from_file(
                                                         core::membrane::io::EmbedSearchParamsOptions const & opts )
            {

                using namespace core::membrane::io;

                TR << "Initializing embedding parameters info using Embed Param options" << std::endl;

                // Create a New object
                EmbedSearchParamsOP params = init_EmbedSearchParams();

                // Initializing Normal Params
                params->normal_search = opts.normal_search();
                params->normal_start_angle = opts.normal_start_angle();
                params->normal_max_angle = opts.normal_max_angle();
                params->normal_delta_angle = opts.normal_delta_angle();

                // Initializing Center params
                params->center_search = opts.center_search();
                params->center_max_delta = opts.center_max_delta();

                // Initializing Normal/Center mag params
                params->center_mag = opts.center_mag();
                params->normal_mag = opts.normal_mag();

                // Initializing normal cycles
                params->normal_cycles = opts.normal_cycles();
                
                // Set booleans
                params->penalties = opts.penalties();
                params->no_interpolate_Mpair = opts.no_interpolate_mpair();

                // Done!
                return params;
            }

} // geometry
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedSearchParamsIO_cc

