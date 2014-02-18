// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       SpanFileIO.hh
///
/// @brief      Generates per-chain membrane spanning data
/// @details    Generates membrane spanning topology data from Octopus topology
///             prediction information. Topology can be generated at http://octopus.cbr.su.se/
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_SpanFileIO_cc
#define INCLUDED_core_membrane_io_SpanFileIO_cc

// Unit Headers
#include <core/membrane/io/SpanFileIO.hh>

// Project Headers
#include <core/membrane/properties/SpanningTopology.hh>
#include <core/membrane/util/Exceptions.hh>

// Package Headers
#include <core/types.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>

using basic::Error;
using basic::Warning;

using namespace core;
using namespace core::membrane::properties;

static basic::Tracer TR( "core.membrane.io.SpanFileIO" );

namespace core {
namespace membrane {
namespace io {

    /// @brief Constructor
    /// @param [none]
    SpanFileIO::SpanFileIO() :
        utility::pointer::ReferenceCount()
    {}

    /// @brief Destructor
    /// @param [none]
    SpanFileIO::~SpanFileIO()
    {}


    /// @brief     Read Spanfile from Locator ID
    /// @details   Read in span file using izstream
    ///
    /// @param     spanfile
    ///             spanfile locoator ID
    void SpanFileIO::read_spanfile(
                    SpanningTopologyOP topology,
                    std::string spanfile
                    )
    {
        
        // Setup vars for reading spanfile using izstream
        std::string line;
        utility::io::izstream stream (spanfile);

        // Read file Header "TM region prediction for"
        getline(stream, line);
        TR << line << std::endl;

        // Read line which includes number of tm spans and total resnum
        getline(stream, line);
        TR << line << std::endl;
        std::istringstream l(line);
        core::Size total_tmhelix;
        core::Size total_residue_in_span_file;

        TR << line << std::endl;
        l >> total_tmhelix >> total_residue_in_span_file;
        
        topology->set_total_tmhelix( total_tmhelix );
        topology->set_total_residue_in_spanfile( total_residue_in_span_file );

        // If there is a negative number of total residues, throw exception
        if ( topology->total_residue_in_span_file() <= 0 ) {
          throw utility::excn::EXCN_Msg_Exception( "SpanFileIO: No residues in pose - check file format or data" );
        }

				// If there are no tm helices, tell the user to stop using a globular protein
				if ( topology->total_tmhelix() <= 0 ) {
					throw utility::excn::EXCN_Msg_Exception( "SpanFileIO: No TM helices in file - check file format or that you are actually using a membrane protein" );
			  }

        // Set info vars using total_tmhelix - make sizes of vectors
        Size max_tmhelix( topology->total_tmhelix() );
        ObjexxFCL::FArray1D< Size > helix_id;
        ObjexxFCL::FArray2D< Size > span;
        ObjexxFCL::FArray2D< Size > full_span;
        ObjexxFCL::FArray2D< Size > relative_tmh_ori;
        
        span.dimension( max_tmhelix, 2);
        full_span.dimension( max_tmhelix, 2);
        relative_tmh_ori.dimension( max_tmhelix, max_tmhelix);
        helix_id.dimension( max_tmhelix );

        // Read in antiparallel line
        getline(stream, line);
        TR << line << std::endl;

        // Read in n2c line
        getline(stream, line);
        TR << line << std::endl;

        // For each line of the file, get spanning region info
        for ( Size i = 1; i <= topology->total_tmhelix(); ++i ) {

            // Reading
            getline(stream, line);
            std::istringstream l(line);
            l >> span(i, 1) >> span(i, 2);

            // Check that spanning regions are in bounds
            if ( span(i, 1) > total_residue_in_span_file || span(i, 1) <= 0 ) {
                throw utility::excn::EXCN_Msg_Exception("SpanFileIO: Membrane spanning topology region out of bounds of total residue!");
            }

            if ( span(i, 2) > total_residue_in_span_file || span(i, 2) <= 0 ) {
                throw utility::excn::EXCN_Msg_Exception("SpanFileIO: Membrane spanning topology region out of bounds of total residue!");
            }

            TR << line << std::endl;

            // Set helix id
            helix_id(i) = i;
        }

        // Set arrays in the obj
        topology->set_span( span );
        topology->set_full_span( full_span );
        topology->set_helix_id( helix_id );
        topology->set_relative_tmh_ori( relative_tmh_ori );

        // Close the izstream
        stream.close();
        stream.clear();

        // Done
        return;
    }

    /// @brief Set up span and full span info
    /// @param Topology Info Object
    void SpanFileIO::setup_span_info( SpanningTopologyOP topology )
    {

        ObjexxFCL::FArray2D< Size > span = topology->span();
        ObjexxFCL::FArray2D< Size > full_span = topology->full_span();
        
        // Set full spanning info
        full_span(1,1)=1;
        full_span( topology->total_tmhelix(), 2 ) = topology->total_residue_in_span_file();

        // Loop through full_span_ and set equal to span
        for ( Size reg1 = 2; reg1 <= topology->total_tmhelix(); ++reg1 ) {
            full_span(reg1, 1) = span(reg1, 1);
            full_span(reg1-1, 2) = span(reg1, 1);
        }
        
        // Set spans
        topology->set_full_span( full_span );
        topology->set_span( span );

        // Done
        return;
    }

    /// @brief Set values in tmregion (true/false for is_tmregion)
    /// @param Topology Info Object
    void SpanFileIO::setup_tmregion( SpanningTopologyOP topology )
    {
 
        // Initialize vars
        utility::vector1< bool > allow_scoring = topology->allow_scoring();
        utility::vector1< bool > allow_tmh_scoring = topology->allow_tmh_scoring();
        utility::vector1< bool > tmregion = topology->tmregion();
        
        
        // Set member vars from span file info
        tmregion.resize( topology->total_residue_in_span_file(), false );
        allow_scoring.resize( topology->total_residue_in_span_file(), true );
        allow_tmh_scoring.resize( topology->total_tmhelix(), true );
        topology->set_tmh_inserted( topology->total_tmhelix() );

        for( Size i = 1; i <= topology->total_residue_in_span_file(); ++i ) {
            for ( Size reg1 = 1; reg1 <= topology->total_tmhelix(); ++reg1 ) {
                if( i >= topology->span()(reg1, 1) && i <= topology->span()(reg1, 2) ) {
                    tmregion[i] = true;
                    continue;
                }
            }
        }
        
        // Set final arrays
        topology->set_allow_scoring( allow_scoring );
        topology->set_allow_tmh_scoring( allow_tmh_scoring );
        topology->set_tmregion( tmregion );

        // Done
        return;
    }

    /// @brief Set relative tmh info
    /// @param Topology Info Object
    void SpanFileIO::setup_relative_tmh( SpanningTopologyOP topology )
    {

        ObjexxFCL::FArray2D< Size > relative_tmh_ori = topology->relative_tmh_ori();
        
        // Loop through relative_tmh_ori and set info
        for ( Size reg1 = 1; reg1 <= topology->total_tmhelix(); ++reg1 ) {
            for ( Size reg2 = 1; reg2 <= topology->total_tmhelix(); ++ reg2 ) {
                relative_tmh_ori( reg1, reg2 ) = 1;
            }
        }
        
        topology->set_relative_tmh_ori( relative_tmh_ori );

        // Done
        return;
    }

    /// @brief Get Topology from SpanFile
    SpanningTopologyOP SpanFileIO::get_topology_from_spanfile( std::string spanfile )
    {

        using namespace core::membrane::util;

        TR << "Building topology object from spanfile" << std::endl;

        // Initialize space for a new topology object
        SpanningTopologyOP topology = new SpanningTopology();

        // Read in spanfile
        read_spanfile(topology, spanfile);

        // Setup spanning info
        setup_span_info(topology);

        // Set up tmregion array
        setup_tmregion(topology);

        // Setup relative tmh
        setup_relative_tmh(topology);

        // All done!
        return topology;
    }


} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_SpanFileIO_cc


