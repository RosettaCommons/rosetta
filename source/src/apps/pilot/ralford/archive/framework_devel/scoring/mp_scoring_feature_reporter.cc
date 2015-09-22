// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       apps/pilot/ralford/mp_scoring_features_reporter.cc
///
/// @brief      Report Features of Residues 
/// @details    Report neighbor counts for 6, 10, and 12A radii, SASA scoring for 1.4A radii, 
///				residue layer assignemnts for 5 layer scheme, residue layer assignemnts for 2 layer
///				scheme
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (4/15/14)


// App headers
#include <devel/init.hh> 

// Rosetta Package Headers
#include <core/pose/Pose.hh> 
#include <core/import_pose/import_pose.hh> 
#include <core/types.hh> 

// Options System
#include <basic/options/option.hh> 
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Compute Features
#include <core/scoring/sasa.hh>
#include <core/scoring/sasa/SasaCalc.hh> 
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/conformation/PointGraphData.hh>
#include <core/graph/UpperEdgeGraph.hh>

#include <core/conformation/Residue.hh> 

// Computing Energy Terms
#include <core/scoring/MembranePotential.hh> 
#include <core/scoring/membrane/MPEnvEnergy.hh> 
#include <core/scoring/Energies.hh> 


// Basic headers
#include <basic/Tracer.hh> 
#include <numeric/xyzVector.hh> 

// C++ Headers
#include <cstdlib>
#include <algorithm> 

using namespace core::conformation;

static THREAD_LOCAL basic::Tracer TR( "mp_scoring_features_reporter" );

/// @brief Compute Sasa
core::Real compute_sasa( core::pose::PoseOP pose, core::Size resnum ) {

    using namespace core::scoring::sasa; 

    SasaCalcOP sasa = new SasaCalc();
    sasa->set_probe_radius(1.4);
    sasa->calculate(*pose);
    return sasa->get_residue_sasa()[ resnum ];

}

/// @brief Compute 5 Layer Designation
core::Size compute_5layer_designation( core::Real z_coord ) {

    // Inner Hydrophobic layer (+/- 12A)
    if ( z_coord > -12.0 && z_coord < 12.0 ) {
        return 1;

    // Outer Hydrophobic Layer (+/-)
    } else if ( z_coord > -18.0 && z_coord < 18.0 ) {
        return 2;

    // Interface (+/- 24A)   
    } else if ( z_coord > -24.0 && z_coord < 24.0 ) {
        return 3;

    // Polar (+/- #0A)
    } else if ( z_coord > -30.0 && z_coord < 30.0 ) { 
        return 4;
    
    // Water layer
    } else {
        return 5;
    }

}

/// @brief Compute Two Layer Designation
core::Size compute_2layer_designation( core::Real z_coord ) {
    
    // Hydrophobic layer (+/- 12.5A)
    if ( z_coord > -12.0 && z_coord < 12.0 ) {
        return 2;

    // water layer    
    } else { 
        return 1;
    } 
}

/// @brief Compute Per-Residue Neighbors within 6A Radii
utility::vector1< core::Size > compute_sixAneighbors( core::Size resnum, core::conformation::PointGraphOP pg ) {

    using namespace core::scoring;
    using namespace core::conformation;
    
    // Vector for storing sizA neighbors
    utility::vector1< core::Size > sixAneighbors;

    for ( PointGraph::UpperEdgeListConstIter
            i_iter     = pg->get_vertex( resnum ).const_upper_edge_list_begin(),
            i_end_iter = pg->get_vertex( resnum ).const_upper_edge_list_end();
            i_iter != i_end_iter; ++i_iter ) {  

        sixAneighbors.push_back( i_iter->upper_vertex() ); 
    }

    return sixAneighbors;  

}

/// @brief Compute Per-Residue Neighbors within 10A Radii
utility::vector1< core::Size > compute_tenAneighbors( core::Size resnum, core::conformation::PointGraphOP pg ) {

    using namespace core::scoring;
    using namespace core::conformation;
    

    // Vector for storing sizA neighbors
    utility::vector1< core::Size > tenAneighbors;

    for ( PointGraph::UpperEdgeListConstIter
            i_iter     = pg->get_vertex( resnum ).const_upper_edge_list_begin(),
            i_end_iter = pg->get_vertex( resnum ).const_upper_edge_list_end();
            i_iter != i_end_iter; ++i_iter ) {  

        tenAneighbors.push_back( i_iter->upper_vertex() ); 
    }
    
    return tenAneighbors;  

}

/// @brief Compute Per-Residue Neighbors within 10A Radii
utility::vector1< core::Size > compute_twelveAneighbors( core::Size resnum, core::conformation::PointGraphOP pg ) {

    using namespace core::scoring;
    using namespace core::conformation;

    // Vector for storing sizA neighbors
    utility::vector1< core::Size > twelveAneighbors;

    for ( PointGraph::UpperEdgeListConstIter
            i_iter     = pg->get_vertex( resnum ).const_upper_edge_list_begin(),
            i_end_iter = pg->get_vertex( resnum ).const_upper_edge_list_end();
            i_iter != i_end_iter; ++i_iter ) {  

        twelveAneighbors.push_back( i_iter->upper_vertex() ); 
    }
    
    return twelveAneighbors;  

}

/// @brief Main Function
int main( int argc, char* argv[] )
{
    try {
        
        using namespace basic::options;
        using namespace std;

        // Initialize Options System, RG, and All Factory_Registrators
        devel::init(argc, argv);
        
        TR << "Membrane Energy Funciton Feature Reporting" << std::endl;
        TR << "Compute membrane scoring features from Rosetta" << std::endl;
		
        // Read in user options
        std::string pdbfile = option[ OptionKeys::in::file::native ]();

        // Set up a pose from pdb
        core::pose::PoseOP pose = new core::pose::Pose(); 
        core::import_pose::pose_from_pdb( *pose, pdbfile ); 

        // Create six A neighbor point graph
        PointGraphOP sixPg( new PointGraph );
        residue_point_graph_from_conformation( pose->conformation(), *sixPg );
        find_neighbors<PointGraphVertexData, PointGraphEdgeData>( sixPg, 6.0 );

        // Create ten A neighbor point graph 
        PointGraphOP tenPg( new PointGraph );
        residue_point_graph_from_conformation( pose->conformation(), *tenPg );
        find_neighbors<PointGraphVertexData, PointGraphEdgeData>( tenPg, 10.0 );

        // Create twelve A neighbor point graph 
        PointGraphOP twelvePg( new PointGraph );
        residue_point_graph_from_conformation( pose->conformation(), *twelvePg );
        find_neighbors<PointGraphVertexData, PointGraphEdgeData>( twelvePg, 12.0 );


        std::cout << "Residue Z_coord Five_Layer Two_Layer SASA sixAnbr tenAnbr twelveAnbr" << std::endl;

        // For each residue, compute score features and print
        for ( core::Size i = 1; i <= pose->total_residue(); ++i ) {

            // Compute a score feature
            core::Real z_coord = pose->residue(i).xyz(2).z();

            // Update score feature object!
            core::Size fiveLayer = compute_5layer_designation( z_coord ); 
            core::Size twoLayer = compute_2layer_designation( z_coord ); 
            core::Real sasa = compute_sasa( pose, i ); 
            utility::vector1< core::Size > sixANeighbors = compute_sixAneighbors( i, sixPg ); 
            utility::vector1< core::Size > tenANeighbors = compute_tenAneighbors( i, tenPg ); 
            utility::vector1< core::Size > twelveANeighbors = compute_twelveAneighbors( i, twelvePg ); 

            std::cout << i << " " << z_coord << " " << fiveLayer << " " << twoLayer << " " << sasa;

            // Print 6 neighbor list
            std::cout << " [";
            for ( core::Size i = 1; i <= sixANeighbors.size(); i++ ) {
                
                if ( sixANeighbors[i] != sixANeighbors.back() ) {
                  cout << sixANeighbors[i] << ",";  
                } else {
                    cout << sixANeighbors[i];
                }
            }
            cout << "]";

            // Print 10 neighbor list
            std::cout << " [";
            for ( core::Size i = 1; i <= tenANeighbors.size(); i++ ) {
                
                if ( tenANeighbors[i] != tenANeighbors.back() ) {
                  cout << tenANeighbors[i] << ",";  
                } else {
                    cout << tenANeighbors[i];
                }
            }
            cout << "]";

            // Print 12 neighbor list
            std::cout << " [";
            for ( core::Size i = 1; i <= twelveANeighbors.size(); i++ ) {
                
                if ( twelveANeighbors[i] != twelveANeighbors.back() ) {
                  cout << twelveANeighbors[i] << ",";  
                } else {
                    cout << twelveANeighbors[i];
                }
            }
            cout << "]" << std::endl;


        }
        
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cout << "caught exception " << e.msg() << std::endl;
    }

}
