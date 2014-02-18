// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       LipoFileIO.cc
///
/// @brief      Lipid Exposure Data from Membrane Spanning Topology
/// @details    Lipid exposure data calculated using run_lips.pl from octopus spanning
///             topology data. Requires blast.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_LipoFileIO_cc
#define INCLUDED_core_membrane_io_LipoFileIO_cc

// Unit Headers
#include <core/membrane/io/LipoFileIO.hh>

// Project Headers
#include <core/membrane/properties/LipidAccInfo.hh>
#include <core/membrane/util/Exceptions.hh>

// Package Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

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

static basic::Tracer TR( "core.membrane.io.LipoFileIO" );

using namespace core::membrane;
using namespace core::membrane::properties;

namespace core {
namespace membrane {
namespace io {

    /// @brief Constructor
    LipoFileIO::LipoFileIO() :
        utility::pointer::ReferenceCount()
    {}

    /// @brief Destructor
    LipoFileIO::~LipoFileIO()
    {}

    /// @brief Read and store data from lips file
    /// @param Lipid Info object, lipid info file, and pose
    void LipoFileIO::read_lips_data(
        LipidAccInfoOP lipid_exp,
        std::string lipsfile
        )
    {

        using namespace core;
        using namespace core::membrane::util;
        
        TR << "Initializing lips exposure info using " << lipsfile << std::endl;
        
        // Initialize local vars
        Size num_of_csts(0);
        Size resnum;
        Real exposure;

        // Initialize izstream
        std::string line;
        utility::io::izstream stream (lipsfile);

        // Open stream and start reading
        stream.open(lipsfile);
        if (stream) {
            
            // Grab the first line
            getline(stream, line);
            getline(stream, line);
            
            // Grab total residues
            int nres = std::atoi(line.c_str());
            TR << nres << std::endl;
            
            // Initialize vectors to match sizes
            utility::vector1< core::Real > lipid_exposure;
            utility::vector1< core::Real > lipid_burial;
            
            lipid_exposure.resize(nres, 0.0);
            lipid_burial.resize(nres, 0.0);

            while ( !stream.eof() ) {
                
                std::istringstream l(line);

                // Add exposure info to vectors
                l >> resnum;
                l >> exposure;
                
                if ( exposure > 0 ) {
                
                    lipid_exposure[ resnum ] = exposure;
                    num_of_csts++;
                
                } else {
                
                    lipid_burial[ resnum ] = abs( exposure );
                    num_of_csts++;
                }

                getline(stream, line);
            }
            
            // Set final
            lipid_exp->set_lipid_exposure( lipid_exposure );
            lipid_exp->set_lipid_burial( lipid_burial );
            
            stream.close();
            stream.clear();

            TR << num_of_csts << " exposure constraints read!" << std::endl;

        } else {
            throw new core::membrane::util::EXCN_Illegal_Arguments("Lips data file not found!");
        }

        // Done!
        return;
    }

    /// @brief Main IO Funciton - Reads Lipid Acc data from file
    /// @param [ lipofile ]
    LipidAccInfoOP LipoFileIO::get_lips_exp_from_lipofile( std::string lipofile )
    {

        using namespace core::membrane::properties;
        using namespace core::membrane::io;

        // CHecking for invalid input
        if ( lipofile.compare("") == 0 ) {
            throw new core::membrane::util::EXCN_Resource_Definition("Must specify a lips file to read!");
        }
        
        // Create new object
        LipidAccInfoOP lips_exp = new LipidAccInfo();
        
        // read in Lips data
        read_lips_data( lips_exp, lipofile );
        
        // Done!
        return lips_exp;
    }

} // geometry
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_LipoFileIO_cc



