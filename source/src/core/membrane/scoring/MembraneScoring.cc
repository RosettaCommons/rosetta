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

#ifndef INCLUDED_core_membrane_scoring_MembraneScoring_cc
#define INCLUDED_core_membrane_scoring_MembraneScoring_cc

// Unit headers
#include <core/membrane/scoring/MembraneScoring.hh>
#include <core/scoring/EnvPairPotential.hh>

// Project Headers
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

// Package Headers
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <basic/database/open.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

// Utility headers
#include <utility/io/izstream.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.membrane.scoring.MembraneScoring");
//static numeric::random::RandomGenerator RG(280628);  // <- Magic number, do not change it! <- :P bad magic number HA I WIN

namespace core {
namespace membrane {
namespace scoring {

        /// @brief 		Standard Constructor (seriously world's longest constructor)
        /// @details 	Initialize values for bins
        ///
        /// @param [none]
        /// @return MembraneScoring Object
        MembraneScoring::MembraneScoring() :
            core::scoring::EnvPairPotential(),
    
            calculated_(false),
            params_(NULL),
    
            max_aa(20),
            cen_dist5_pad( 0.5 ),
            cen_dist6_pad( 0.6 ),
            cen_dist7_pad( 0.65 ),
            cen_dist10_pad( 1.0 ),
            cen_dist12_pad( 1.2 ),

            cen_dist5_pad_plus ( cen_dist5_pad  + 25.0 ),
            cen_dist6_pad_plus( cen_dist6_pad + 36.0 ),
            cen_dist7_pad_plus ( cen_dist7_pad  + 56.25 ),
            cen_dist10_pad_plus( cen_dist10_pad + 100.0 ),
            cen_dist12_pad_plus( cen_dist12_pad + 144.0 ),
            
            cen_dist5_pad_minus ( cen_dist5_pad  - 25.0 ),
            cen_dist7_pad_minus ( cen_dist7_pad  - 56.25 ),
            cen_dist10_pad_minus( cen_dist10_pad - 100.0 ),
            cen_dist12_pad_minus( cen_dist12_pad - 144.0 ),
            
            cen_dist5_pad_hinv ( 0.5 / cen_dist5_pad ),
            cen_dist6_pad_hinv ( 0.5 / cen_dist6_pad ),
            cen_dist7_pad_hinv ( 0.5 / cen_dist7_pad ),
            cen_dist10_pad_hinv( 0.5 / cen_dist10_pad ),
            cen_dist12_pad_hinv( 0.5 / cen_dist12_pad )
        {
            
            using namespace core::membrane::util;
            
            // Load required resources (defualt resource descriptions)
            load_requried_resources();
            
            // Continue to Load DB Data
            Size const env_log_table_cen6_bins( 15 );
            Size const env_log_table_cen10_bins( 40 );
            Size const pair_log_table_size( 5 );
            Size const cbeta_den_table_size( 45 );
            Size const max_mem_layers( 3 );
            Size const min_mem_layers( 2 );
            
            // Load DB Data
            load_db( env_log_table_cen6_bins, env_log_table_cen10_bins, pair_log_table_size, cbeta_den_table_size, max_mem_layers, min_mem_layers );
        }
            
        ///////////////////////////////////////////////////////////////////////////////////////////////
    
        /// @brief    Evaluate Membrane Environemnt
        /// @details  <TBD>
        ///
        /// @param    pose
        ///             pose to score
        /// @param    residue
        ///             residue to score
        /// @param    membrane env score
        ///             ref to final score
        void
        MembraneScoring::evaluate_env(
                                  pose::Pose const & pose,
                                  conformation::Residue const & rsd,
                                  Real const MembraneDepth,
                                  Vector const & normal,
                                  Vector const & center,
                                  Real & membrane_env_score
                                  ) const
        {
            
            // Removing If Block for spanning data - should ALWAYS be present!
            Real termini_pen(0);
            
            evaluate_env( pose, rsd, MembraneDepth, membrane_env_score);
                
            // Evaluate membrane environment penalties
            if( params_->penalties && ( rsd.seqpos()==1 || rsd.seqpos()==pose.total_residue() ) ) {
                    
                Vector const & xyz( pose.residue(rsd.seqpos()).atom( 2 ).xyz() );
                Real depth = dot(xyz-center,normal)+30;
                if(depth>18 && depth<42) { termini_pen++; }
                membrane_env_score+=50*termini_pen;
            }
        }
            
            
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
        void
        MembraneScoring::evaluate_env(
                        pose::Pose const & pose,
                        conformation::Residue const & rsd,
                        Real const MembraneDepth,
                        Real & membrane_env_score
                        ) const
        {
            
            // Initialize Some Parameters
            Real t2 = 2.0;
            Real t3 = 2.0;
            int s2 = 14;
            int s3 = 14;
            int layer1, layer2, layer;
            Real f, z, zn, low;
            
            // Initialize weights (WHY ARE THESE FUCKING HARD CODED)
            Real const env6_weight=1.0;
            Real const env10_weight=1.0;
            
            // Grab Cenlist Data from pose
            Real fcen6  ( cenlist_from_pose( pose ).fcen6( rsd.seqpos() ) );
            Real fcen10 ( cenlist_from_pose( pose ).fcen10( rsd.seqpos() ) );
            
            // in rare cases, the density is over 15 within 6Å
            if (fcen6 > 15) fcen6 = 15 ;
            if (fcen10 > 40) fcen10 = 40;
            
            // Skip virtual residues
            if ( rsd.is_protein() ) {
                
                if ( ( MembraneDepth < 11.0 ) || ( MembraneDepth > 49.0 ) ) {
                    
                    //pure water layer
                    layer = 3;
                    Real score6 (env6_weight*mem_env_log6_( rsd.aa(), layer, static_cast< int >( fcen6 ) ));
                    Real score10 (env10_weight* mem_env_log10_( rsd.aa(), layer, static_cast< int >( fcen10 ) ) );
                    membrane_env_score = score6 + score10;
                    
                }
                
                else if ( ( MembraneDepth >= 11.0 && MembraneDepth <= 13.0 ) || ( MembraneDepth >= 47.0 && MembraneDepth <= 49.0 ) ) {
                    
                    //interpolate between water and interface phases
                    layer1 = 2; //interface layer
                    layer2 = 3; //water layer
                    
                    if ( MembraneDepth <= 13.0 ) { low = 13.0; }
                    else { low = 47.0; }
                    
                    z = 2*std::abs( (MembraneDepth - low) ) / t2;
                    zn = std::pow( z, s2 );
                    f = zn/(1 + zn);
                    
                    Real score6_layer2( env6_weight*mem_env_log6_( rsd.aa(), layer2, static_cast< int >( fcen6 ) ) );
                    Real score10_layer2( env10_weight* mem_env_log10_( rsd.aa(), layer2, static_cast< int >( fcen10 ) ) );
                    Real score6_layer1( env6_weight*mem_env_log6_( rsd.aa(), layer1, static_cast< int >( fcen6 ) ) );
                    Real score10_layer1( env10_weight*mem_env_log10_( rsd.aa(), layer1, static_cast< int >( fcen10 ) ) );
                    
                    membrane_env_score = f * ( score6_layer2 + score10_layer2 ) + ( 1 - f ) * ( score6_layer1 + score10_layer1 );
                    
                    if ( MembraneDepth <= 12.0 || MembraneDepth >= 48.0 ) { layer = 2; }
                    else { layer = 3; }
                    
                } else if ( ( MembraneDepth > 13.0 && MembraneDepth < 17.0 ) || ( MembraneDepth > 43.0 && MembraneDepth < 47.0 ) ) {
                    
                    //pure interface phase
                    layer = 2; //interface layer
                    
                    Real score6 ( env6_weight*mem_env_log6_( rsd.aa(), layer, static_cast< int >( fcen6 ) ) );
                    Real score10 ( env10_weight*mem_env_log10_( rsd.aa(), layer, static_cast< int >( fcen10 ) ) );
                    membrane_env_score = score6 + score10;
                
                } else if ( ( MembraneDepth >= 17.0 && MembraneDepth <= 19.0 ) || ( MembraneDepth >= 41.0 && MembraneDepth <= 43.0 ) ) {
                
                    //interpolate between interface and hydrophobic phases
                    layer1 = 1; //hydrophobic layer
                    layer2 = 2; //interface layer
                    
                    if ( MembraneDepth <= 19.0 ) { low = 19.0; }
                    else { low = 41.0; }
                    
                    z = 2*std::abs( (MembraneDepth - low) ) / t3;
                    zn = std::pow( z, s3 );
                    f = zn/(1 + zn);
                    
                    Real score6_layer2(env6_weight*mem_env_log6_( rsd.aa(), layer2, static_cast< int >( fcen6 )));
                    Real score10_layer2( env10_weight*mem_env_log10_( rsd.aa(), layer2, static_cast< int >( fcen10 ) ) );
                    Real score6_layer1( env6_weight*mem_env_log6_( rsd.aa(), layer1, static_cast< int >( fcen6 ) ) );
                    Real score10_layer1( env10_weight*mem_env_log10_( rsd.aa(), layer1, static_cast< int >( fcen10 ) ) );
                    
                    membrane_env_score = f * ( score6_layer2  + score10_layer2 ) + ( 1 - f ) * ( score6_layer1 + score10_layer1 );
                    
                    if ( MembraneDepth <= 18.0 || MembraneDepth >= 42.0 ) { layer = 2; }
                    else { layer = 1; }
                    
                } else {
                    
                    //pure hydrophobic phase
                    layer = 1;
                    
                    Real score6  (env6_weight *mem_env_log6_(  rsd.aa(), layer, static_cast< int >( fcen6  ) ));
                    Real score10 (env10_weight*mem_env_log10_( rsd.aa(), layer, static_cast< int >( fcen10 ) ));
                    
                    membrane_env_score = score6+score10;
                }
                
                membrane_env_score*=0.5; //bw membrane_embed_weight...
                
            } else { // amino acid check
                membrane_env_score = 0.0;
            }
        }
    
        ///////////////////////////////////////////////////////////////////////////////////////////////
        void
        MembraneScoring::evaluate_cbeta(
                                                 pose::Pose const & pose,
                                                 conformation::Residue const & rsd,
                                                 Real & membrane_cb_score
                                                 ) const
        {
            
            using namespace core::membrane::util;
            
            // Init cbeta score
            membrane_cb_score=0;
            
            // Grab cenlist info from pose
            Real const fcen6 ( cenlist_from_pose( pose ).fcen6( rsd.seqpos() ) );
            Real const fcen12 ( cenlist_from_pose( pose ).fcen12( rsd.seqpos() ) );
            
            Real membrane_cb_score6,membrane_cb_score12;
            
            // Gram Membrane Topology
            Size const TMHs = topology_->tmh_inserted;

            // interp1 rounds down to nearest (non-negative) integer.
            int const interp1 = static_cast< int >( fcen6 );
            int const interp3 = static_cast< int >( fcen12 );

            // note cen6 is always at least 1.0
            // fraction remainder after nearest lower integer is removed
            Real const interp2 = fcen6-interp1;
            Real const interp4 = fcen12-interp3;
            
            if ( TMHs <= 2 ) {
                membrane_cb_score6 =
                ( 1.0-interp2 ) * mem_cbeta_2TM_den6_( interp1 )+
                interp2         * mem_cbeta_2TM_den6_( interp1+1 );
                
            } else if ( TMHs <= 4 ) {
                membrane_cb_score6 =
                (1.0-interp2) * mem_cbeta_4TM_den6_( interp1 )+
                interp2       * mem_cbeta_4TM_den6_( interp1+1 );
                
            } else {
                membrane_cb_score6 =
                (1.0-interp2) * mem_cbeta_den6_( interp1 )+
                interp2       * mem_cbeta_den6_( interp1+1 );
            }
            membrane_cb_score12 =
            (1.0-interp4) * mem_cbeta_den12_( interp3 )+
            interp4       * mem_cbeta_den12_( interp3+1 );
            
            membrane_cb_score = (membrane_cb_score6+membrane_cb_score12);
        }
    
    
    
        ///////////////////////////////////////////////////////////////////////////////////////////////
        void
        MembraneScoring::evaluate_pair(
                                         pose::Pose const & pose,
                                         conformation::Residue const & rsd1,
                                         conformation::Residue const & rsd2,
                                         core::Real const & depth1,
                                         core::Real const & depth2,
                                         Real const cendist,
                                         Real & membrane_pair_score
                                         ) const
        {
            
            // Init Pair Score
            membrane_pair_score = 0.0;
            
            // Check residues are .aa residues and protein residues
            if ( !rsd1.is_protein() || !rsd2.is_protein() ) return;
            
            chemical::AA const aa1( rsd1.aa() );
            chemical::AA const aa2( rsd2.aa() );
            
            //CAR no pair score if a disulfide
            if (	aa1 == chemical::aa_cys && aa2 == chemical::aa_cys &&
                rsd1.is_bonded( rsd2 ) && rsd1.polymeric_sequence_distance( rsd2 ) > 1 &&
                rsd1.has_variant_type( chemical::DISULFIDE ) && rsd2.has_variant_type( chemical::DISULFIDE ) ) return;
            
            // no pair score for residues closer than 9 in sequence
            if ( rsd1.polymeric_sequence_distance( rsd2 ) /* j - i */ <= 8 ) return;
            
            //  we now try to find which bin the pair distance lies in
            //  I note this could in principle be calculated and updatded
            //  just like cen_dist is if there is a need for speed.
            //  this function interpolates between bins.
            //  An important(!) requirement on pair_log is that the
            //  value should approach zero as the radius increases.
            //  this fact permits us not to have to compute and score pairs are larger
            //  than cen_dist > cutoff.
            
            int icon = 5;
            Real interp2( 0.0 );

            Real const MembraneDepth1 (depth1);
            Real const MembraneDepth2 (depth2);
            
            int hydro_layer=1;  //1 not_hydrophobic_core 2 hydrophobic core
            Real AverageDepth=(MembraneDepth1+MembraneDepth2)/2;
            if(MembraneDepth1 > 18 &&
               MembraneDepth1 < 42 &&
               MembraneDepth2 >18 &&
               MembraneDepth2 <42) //bw currently both residues have to be in the hydrophobic core
            {
                hydro_layer=2;
            }
            
            if ( cendist > cen_dist10_pad_plus ) {
                icon = 4;
                interp2 = ( cendist + cen_dist12_pad_minus ) * cen_dist12_pad_hinv;
            } else {
                if ( cendist > cen_dist7_pad_plus ) {
                    icon = 3;
                    interp2 = ( cendist + cen_dist10_pad_minus ) * cen_dist10_pad_hinv;
                } else {
                    if ( cendist > cen_dist5_pad_plus ) {
                        icon = 2;
                        interp2 = ( cendist + cen_dist7_pad_minus ) * cen_dist7_pad_hinv;
                    } else {
                        icon = 1;
                        interp2 = ( cendist + cen_dist5_pad_minus ) * cen_dist5_pad_hinv;
                    }
                }
            }
            if ( interp2 < 0.0 ) interp2 = 0.0;
            
            // note in theory this will never happen but in practice round off
            // error can cause problem
            if ( interp2 > 1.0 ) interp2 = 1.0;
            
            // handle last bin specially since icon+1 would be past array end
            Real f(0);
            if ( icon != 5 ) {
                
                if ( !params_->no_interpolate_Mpair ) {
                    
                    if( std::abs( AverageDepth - 18 )< 4) {

                        f = 1 / ( 1 + std::exp( 1.5 * ( 18 - AverageDepth ) ) );
                        membrane_pair_score = ( ( 1.0f - interp2 ) * ((1-f)*mem_pair_log_( 1, icon  , aa1, aa2 ) + f*mem_pair_log_( 2, icon, aa1, aa2 )) +
                                               ( interp2 ) *  ((1-f)*mem_pair_log_( 1, icon+1, aa1, aa2 ) + f*mem_pair_log_( 2, icon+1, aa1, aa2 )));
                    
                    } else if( std::abs( AverageDepth - 42 ) < 4 ) {
                        
                        f = 1 / ( 1 + std::exp( 1.5 * ( AverageDepth - 42 ) ) );
                        membrane_pair_score = ( ( 1.0f - interp2 ) * ((1-f)*mem_pair_log_( 1, icon  , aa1, aa2 ) + f*mem_pair_log_( 2, icon, aa1, aa2 )) +
                                               ( interp2 ) *  ((1-f)*mem_pair_log_( 1, icon+1, aa1, aa2 ) + f*mem_pair_log_( 2, icon+1, aa1, aa2 )));
                        
                    } else {
                        
                        membrane_pair_score = ( ( 1.0f - interp2 ) * mem_pair_log_( hydro_layer, icon  , aa1, aa2 ) +
                                               ( interp2 ) *  mem_pair_log_( hydro_layer, icon+1, aa1, aa2 ));
                    }
                    
                } else {
                    
                    membrane_pair_score = ( ( 1.0f - interp2 ) * mem_pair_log_( hydro_layer, icon, aa1, aa2 ) +
                                           ( interp2 ) *  mem_pair_log_( hydro_layer, icon+1, aa1, aa2 ));
                
                }
            
            } else {
            
                membrane_pair_score = ( 1.0f - interp2 ) * mem_pair_log_( hydro_layer, icon, aa1, aa2 );
            
            }
            
            membrane_pair_score*=2.019; //why is this multiplied by 2.019???
            return;
        }


        /**
         * Private Helper Functions
         **/
    
    
        /// @brief Helper Function for Constructor
        /// @detail Loads Database info
        ///
        /// @param [none]
        /// @return [none]
        void
        MembraneScoring::load_db(
             Size const env_log_table_cen6_bins,
             Size const env_log_table_cen10_bins,
             Size const pair_log_table_size,
             Size const cbeta_den_table_size,
             Size const max_mem_layers,
             Size const min_mem_layers
        ) {

            // Initial Declarations
            std::string tag,line;
            chemical::AA aa;
            
            // Score - Load Mmembrane Environment Cen6
            {
                mem_env_log6_.dimension( max_aa, max_mem_layers, env_log_table_cen6_bins );
                utility::io::izstream stream;
                basic::database::open( stream, "scoring/score_functions/MembranePotential/CEN6_mem_env_log.txt" );
                for ( Size i=1; i<= max_aa; ++i ){
                    getline( stream, line );
                    std::istringstream l(line);
                    l >> tag >> aa;
                    if ( l.fail() || tag != "MEM_ENV_LOG_CEN6:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/CEN6_mem_env_log.txt (cen6)");
                    for( Size j=1;j<=max_mem_layers;++j) {
                        for( Size k=1;k<=env_log_table_cen6_bins;++k) {
                            l >> mem_env_log6_(aa,j,k);
                        }
                    }
                }
            }
            
            // Score - Load Membrane Environment Cen10
            {
                mem_env_log10_.dimension( max_aa, max_mem_layers,env_log_table_cen10_bins );
                
                utility::io::izstream stream;
                basic::database::open( stream, "scoring/score_functions/MembranePotential/CEN10_mem_env_log.txt" );
                for ( Size i=1; i<= max_aa; ++i ){
                    getline( stream, line );
                    std::istringstream l(line);
                    l >> tag >> aa;
                    if ( l.fail() || tag != "MEM_ENV_LOG_CEN10:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/CEN10_mem_env_log.txt (cen10)");
                    for( Size j=1;j<=max_mem_layers;++j) {
                        for( Size k=1;k<=env_log_table_cen10_bins;++k) {
                            l >> mem_env_log10_(aa,j,k);
                        }
                    }
                }
            }
            
            // cbeta_den_6/12
            {
                mem_cbeta_den6_.dimension( cbeta_den_table_size );
                mem_cbeta_den12_.dimension( cbeta_den_table_size );
                mem_cbeta_2TM_den6_.dimension( cbeta_den_table_size );
                mem_cbeta_2TM_den12_.dimension( cbeta_den_table_size );
                mem_cbeta_4TM_den6_.dimension( cbeta_den_table_size );
                mem_cbeta_4TM_den12_.dimension( cbeta_den_table_size );
                
                utility::io::izstream stream;
                basic::database::open( stream, "scoring/score_functions/MembranePotential/memcbeta_den.txt" );
                
                { // den6
                    getline( stream, line );
                    {
                        std::istringstream l(line);
                        l >>	tag;
                        for ( Size i=1; i<= cbeta_den_table_size; ++i ){
                            l >>mem_cbeta_den6_(i);
                        }
                        if ( l.fail() || tag != "MEMCBETA_DEN6:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (DEN6)");
                    }
                    getline( stream, line );
                    {
                        std::istringstream l(line);
                        l >> tag;
                        for ( Size i=1; i<= cbeta_den_table_size; ++i ){
                            l >> mem_cbeta_2TM_den6_(i);
                        }
                        if ( l.fail() || tag != "MEMCBETA_2TM_DEN6:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (2TM_DEN6)");
                    }
                    getline( stream, line );
                    {
                        std::istringstream l(line);
                        l >> tag;
                        for ( Size i=1; i<= cbeta_den_table_size; ++i ){
                            l >> mem_cbeta_4TM_den6_(i);
                        }
                        if ( l.fail() || tag != "MEMCBETA_4TM_DEN6:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (4TM_DEN6)");
                    }
                }
                
                {
                    { // den12
                        getline( stream, line );
                        std::istringstream l(line);
                        l >> tag;
                        for ( Size i=1; i<= cbeta_den_table_size; ++i ){
                            l >> mem_cbeta_den12_(i);
                        }
                        if ( l.fail() || tag != "MEMCBETA_DEN12:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (DEN12)");
                    }
                    getline( stream, line );
                    {
                        std::istringstream l(line);
                        l >> tag;
                        for ( Size i=1; i<= cbeta_den_table_size; ++i ){
                            l >> mem_cbeta_2TM_den12_(i);
                        }
                        if ( l.fail() || tag != "MEMCBETA_2TM_DEN12:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (2TM_DEN12)");
                    }
                    getline( stream, line );
                    {
                        std::istringstream l(line);
                        l >> tag;
                        for ( Size i=1; i<= cbeta_den_table_size; ++i ){
                            l >> mem_cbeta_4TM_den12_(i);
                        }
                        if ( l.fail() || tag != "MEMCBETA_4TM_DEN12:"  ) utility_exit_with_message("bad format for scoring/score_functions/MembranePotential/memcbeta_den.txt (4TM_DEN12)");
                    }
                    
                }
            }
            
            // pair_log
            {
                mem_pair_log_.dimension( min_mem_layers,pair_log_table_size,max_aa, max_aa );
                
                utility::io::izstream stream;
                basic::database::open( stream, "scoring/score_functions/MembranePotential/mem_pair_log.txt" );
                for ( Size i=1; i<= min_mem_layers;++i) {
                    for ( Size j=1; j<= pair_log_table_size; ++j ) {
                        for ( Size k=1; k<= max_aa; ++k ) {
                            getline( stream, line );
                            std::istringstream l(line);
                            Size ii,jj;
                            l >> tag >> ii >> jj >> aa;
                            assert( Size(aa) == k );
                            for ( Size n=1;n<=max_aa;++n)
                            {
                                l >> mem_pair_log_(i,j,aa,n);
                            }
                            if ( l.fail() || ii != i || jj != j || tag != "MEM_PAIR_LOG:"  ) utility_exit_with_message("scoring/score_functions/MembranePotential/mem_pair_log.txt");
                        }
                        
                    }
                }
            }
        }
    
    /// @brief Load Required resources (will load default by startstruct)
    /// @details Loads in the scoring parameters and the membrane topology object
    void MembraneScoring::load_requried_resources() {
        
        using namespace basic::resource_manager;
        using namespace core::membrane::util;
        
        // load in topology resource
        if ( !ResourceManager::get_instance()->has_resource_with_description("startstruct_topology") ) {
            throw new EXCN_Resource_Definition("Either a resource definition file or the command line option -in:file:spanfile must be specified for membrane proteins");
        }
        topology_ = get_resource< SpanningTopology >("startstruct_topology");
        
        // load in search params resource
        if ( !ResourceManager::get_instance()->has_resource_with_description("startstruct_params") ) {
            throw new EXCN_Resource_Definition("Either a resource definition file or the command line option -in:membrane:params must be specified for membrane protiens");
        }
        params_ = get_resource< EmbedSearchParams >("startstruct_params");
        
    }
    
    /// @brief Change default resources by description
    /// @details Give a resource description specific to a chain or job or whatever
    void MembraneScoring::change_resource_by_desc( std::string desc ) {
        
        using namespace basic::resource_manager;
        using namespace core::membrane::util;
        
        TR << "Changing default scoring resources - start struct params and topology to use resources with the base description " << desc << std::endl;
        TR << "Note: This assumes that for whatever base description you have specified, there are corresponding resources available that are paired correctly. If not, kill your process, fix your descriptions and try again" << std::endl;
        
        // Create resource specific descriptions from base description
        std::string st_desc = desc + "_topology";
        std::string param_desc = desc + "_params";
        
        // load in topology resource
        if ( !ResourceManager::get_instance()->has_resource_with_description(st_desc) ) {
            throw new EXCN_Resource_Definition("Either a resource definition file or the command line option -in:file:spanfile must be specified for membrane proteins");
        }
        topology_ = get_resource< SpanningTopology >(st_desc);
        
        // load in search params resource
        if ( !ResourceManager::get_instance()->has_resource_with_description(param_desc) ) {
            throw new EXCN_Resource_Definition("Either a resource definition file or the command line option -in:membrane:params must be specified for membrane protiens");
        }
        params_ = get_resource< EmbedSearchParams >(param_desc);
        
    }
    
} // scoring
} // membrane
} // core

#endif // INCLUDED_core_membrane_scoring_MembraneScoring_cc