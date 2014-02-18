// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   	core/membrane/scoring/MembraneScoring.cc
///
/// @brief  	Membrane Scoring
/// @details	This class contains all of the independent scoring methods
///				to apply such methods to a membrane bound pose
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_scoring_MembraneScoring_hh
#define INCLUDED_core_membrane_scoring_MembraneScoring_hh

// Unit headers
#include <core/membrane/scoring/MembraneScoring.fwd.hh>
#include <core/scoring/EnvPairPotential.hh>

// Project Headers
#include <core/membrane/util/definitions.hh>

// Package Headers
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <basic/database/open.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility headers
#include <utility/io/izstream.hh>

using basic::Error;
using basic::Warning;

namespace core {
namespace membrane {
namespace scoring {

class MembraneScoring : public core::scoring::EnvPairPotential {
    
public: // methods
        
        /// @brief 		Standard Constructor (seriously world's longest constructor)
        /// @detail 	Initialize values for bins
        ///
        /// @param [none]
        /// @return MembraneScoring Object
        MembraneScoring();
    
        /// @brief    Evaluate Membrane Environemnt
        /// @details  <TBD>
        ///
        /// @param    pose
        ///             pose to score
        /// @param    residue
        ///             residue to score
        /// @param    membrane env score
        ///             ref to final score
        void evaluate_env(
              pose::Pose const & pose,
              conformation::Residue const & rsd,
              Real const MembraneDepth,
              Vector const & normal,
              Vector const & center,
              Real & membrane_env_score
              ) const;
    
            
        /// @brief    Evaluate Membrane Environment
        /// @details  Evaluate membrane environemnt given membrane depth
        ///
        /// @param    pose
        ///             pose to score
        /// @param    residue
        ///             residue env to score
        /// @param    membrane depth
        ///             given membrane depth
        /// @param    membrane environemnt score
        ///             pass nonconst ref to score
        void evaluate_env(
              pose::Pose const & pose,
              conformation::Residue const & rsd,
              Real const MembraneDepth,
              Real & membrane_env_score
              ) const;
    
        /// @brief  Evaluate cbeta term
        /// @details
        ///
        /// @param
        void evaluate_cbeta(
              pose::Pose const & pose,
              conformation::Residue const & rsd,
              Real & membrane_cb_score
              ) const;
    
        /// @brief      Evaluate Pair Term between two residues
        /// @details    Evaluate pair term
        ///
        /// @param      pose
        ///                 pose to score
        /// @param      residue 1
        ///                 first residue
        /// @param      residue 2
        ///                 second residue
        /// @param      pose cendlist
        /// @param      score
        void evaluate_pair(
                 pose::Pose const & pose,
                 conformation::Residue const & rsd1,
                 conformation::Residue const & rsd2,
                 core::Real const & depth1,
                 core::Real const & depth2,
                 Real const cendist,
                 Real & membrane_pair_score
                 ) const;
    
public: // change resources by new desc (possibly chain specific??)
    
        /// @brief Change default resources by description
        /// @details Give a resource description specific to a chain or job or whatever
        ///
        /// @param [none]
        /// @throws EXCN_Resource_Definition (Membrane)
        void change_resource_by_desc( std::string desc );
    
private: // methods
    
        /// @brief Helper Function for Constructor
        /// @detail Loads Database info
        ///
        /// @param [none]
        /// @return [none]
        void load_db(
                 Size const env_log_table_cen6_bins,
                 Size const env_log_table_cen10_bins,
                 Size const pair_log_table_size,
                 Size const cbeta_den_table_size,
                 Size const max_mem_layers,
                 Size const min_mem_layers
        );
    
        /// @brief Helper funciton to load required resources
        /// @details Loads required resources from the commandline
        ///
        /// @param [none]
        /// @throws EXCN_Resource_Definition (Membrane)
        void load_requried_resources();
    
private: // data
    
    // Membrane Database Info
    ObjexxFCL::FArray3D< Real > mem_env_log6_;
    ObjexxFCL::FArray3D< Real > mem_env_log10_;
    ObjexxFCL::FArray1D< Real > mem_cbeta_den6_;
    ObjexxFCL::FArray1D< Real > mem_cbeta_den12_;
    ObjexxFCL::FArray1D< Real > mem_cbeta_2TM_den6_;
    ObjexxFCL::FArray1D< Real > mem_cbeta_2TM_den12_;
    ObjexxFCL::FArray1D< Real > mem_cbeta_4TM_den6_;
    ObjexxFCL::FArray1D< Real > mem_cbeta_4TM_den12_;
    ObjexxFCL::FArray4D< Real > mem_pair_log_;
    
    // Cenlist calculated info
    bool calculated_;
    
    // Membrane Scoring Parameters
    core::membrane::util::EmbedSearchParamsOP params_;
    core::membrane::util::SpanningTopologyOP topology_;
    
    // Max allowed
    Size const max_aa;
    
    // Env Pair Potential Cendist Info
    Real const cen_dist5_pad;
    Real const cen_dist6_pad;
    Real const cen_dist7_pad;
    Real const cen_dist10_pad;
    Real const cen_dist12_pad;
    
    Real const cen_dist5_pad_plus ;
    Real const cen_dist6_pad_plus ;
    Real const cen_dist7_pad_plus ;
    Real const cen_dist10_pad_plus;
    Real const cen_dist12_pad_plus;
    
    Real const cen_dist5_pad_minus ;
    Real const cen_dist7_pad_minus ;
    Real const cen_dist10_pad_minus;
    Real const cen_dist12_pad_minus;
    
    Real const cen_dist5_pad_hinv ;
    Real const cen_dist6_pad_hinv ;
    Real const cen_dist7_pad_hinv ;
    Real const cen_dist10_pad_hinv;
    Real const cen_dist12_pad_hinv;
    
    
}; // class Membrane Scoring
    
} // scoring
} // membrane
} // core



#endif // INCLUDED_core_membrane_scoring_MembraneScoring_hh