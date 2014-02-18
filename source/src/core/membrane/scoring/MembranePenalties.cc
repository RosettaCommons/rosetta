// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file           MembranePenalties.cc
///
/// @brief          Membrane Penalties
/// @details		This class's methods evaluate scoring penalties for the membrane region
///                 This will later get instantiated in MembraneSearch and MembranePtoential
///
/// @note           Updated Documentation - Documented from Vladmir et al. 2006
///
/// @author         Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_scoring_MembranePenalties_cc
#define INCLUDED_core_membrane_scoring_MembranePenalties_cc

// Unit Headers
#include <core/membrane/scoring/MembranePenalties.hh>

// Project Headers
#include <core/membrane/util/definitions.hh>

// Package Headers
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>

#include <core/types.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

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

static basic::Tracer TR("core.membrane.scoring.MembranePenalties");

/// @brief 	Class: MembranePenalties
/// @detail This class contains methods for evaluating and applying penalties
///			for tm helix orientations relative to pose membrane center and normal.

namespace core {
namespace membrane {
namespace scoring {

	/// @brief      Constructor
	/// @details	Creates an instance of the Membrane Penalties object
    MembranePenalties::MembranePenalties() :
        utility::pointer::ReferenceCount()
    {}

	/// @brief      Evlauate TM Projection Penalty
	/// @details    Penalty for TM helixes that do not follow projection from
    ///             normal and center
	void
    MembranePenalties::tm_projection_penalty(
                                     std::string desc,
                                     std::string chain,
                                     Vector const & normal,
                                     Vector const & center,
                                     core::Real & tm_proj
                                     )
	{
        
        using namespace core::membrane::util;
        
		// Initialize TM Projection
		tm_proj = 0.0;
        
        // Get pose and topology (dummy lines)
        core::pose::PoseOP pose = new Pose();
        core::membrane::util::SpanningTopology topology;

		// Define vectors for inside and outside cap residue
		Vector inside(0);
		Vector outside(0);

		for( Size i = 1; i <= topology->total_tmhelix; ++i )
		{
			// Check that the given position allows scoring
			if( !topology->allow_tmh_scoring[i] ) continue;

			// Get start and end coordinates from the pose
			Vector const & start( pose->residue( topology->span(i, 1) ).atom( 2 ).xyz());
			Vector const & end( pose->residue( topology->span(i, 2) ).atom( 2 ).xyz());

			// Calculate helix length
			Real tm_length = std::abs( dot( start-center, normal ) - dot( end-center, normal ) );

			// Calculate ratio
			Real ratio = tm_length / ( topology->span(i, 2)-topology->span(i, 1)+1 );

			// Evaluate penalty
			if(tm_length<15) { tm_proj++; }
			if(ratio<1 || ratio > 1.5) { tm_proj++; }
		}

		tm_proj*=50;
	}

	/// @brief      Evaluate Non Membrane Helix Penalty
	/// @details    Penalty for Predicted TM Helixes that do not span the membrane
	void MembranePenalties::non_helix_in_membrane_penalty(
                                       std::string desc,
                                       std::string chain,
                                       Vector const & normal,
                                       Vector const & center,
                                       Real & non_helix_pen
                                       )
	{

        using namespace core::membrane::util;
        using namespace core::conformation::symmetry;
        
		// Initialize non_helix penalty
		non_helix_pen=0.0;

        // Get pose and topology
        core::pose::PoseOP pose = load_pose(desc, chain);
        SpanningTopologyOP topology = load_topology(desc, chain);

		// Loop through the posee and calculate penalty at each res
		for ( Size i = 0; i <= pose->total_residue(); ++i ) {

			Size rsdSeq(i);

			if ( core::pose::symmetry::is_symmetric( *pose ) ) {
				
				SymmetricConformation const & symm_conf (
													 dynamic_cast< SymmetricConformation const & > ( pose->conformation() ) );
				SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

				if ( !symm_info->bb_is_independent( pose->residue(i).seqpos()) ) {
					rsdSeq = symm_info->bb_follows( pose->residue(i).seqpos() );
				}

				if (symm_info->is_virtual(i)) { rsdSeq = 0; }
			}

            if (rsdSeq ==0 ) continue; // skip virtual residue
            if ( pose->residue(rsdSeq).aa() == core::chemical::aa_vrt ) continue;
            if( !topology->allow_scoring[rsdSeq] ) continue;

            if( topology->tmregion[rsdSeq] && pose->conformation().secstruct(i)!='H' ) {
                Vector const & xyz( pose->residue( i ).atom( 2 ).xyz());
                Real depth=dot(xyz-center,normal)+30; // const 30 - cannot readjust
                
                // Increment non helix penalty if depth out of bounds
                if( depth > 18 && depth < 42 ) { non_helix_pen++; } // const - also hard coded THANKS PATRICK :P
            }
		}

		non_helix_pen *= 10;
	}

	/// @brief      Evaluate Termini Penalty
	/// @details    Penalty for termini positions contradictory to membrane
    ///             spanning
	void
    MembranePenalties::termini_penalty(
                                       std::string desc,
                                       std::string chain,
                                       Vector const & normal,
                                       Vector const & center,
                                       Real & termini_pen
                                       )
	{
        
        using namespace core::membrane::util;
        using namespace core::conformation::symmetry;
        
		// Initialize termini penalty
		termini_pen=0.0;

        // Get pose and topology
        core::pose::PoseOP pose = load_pose(desc, chain);
        SpanningTopologyOP topology = load_topology(desc, chain);

		// Loop through the pose and calculate termini penalty per residue
		for(Size i=1;i<=pose->total_residue();++i) {

			if (!pose->residue(i).is_terminus()) continue;

			Size rsdSeq(i);
			if ( core::pose::symmetry::is_symmetric( *pose ) ) {
				
				SymmetricConformation const & symm_conf (
					dynamic_cast< SymmetricConformation const & > ( pose->conformation() ) );
				SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

                if (!symm_info->bb_is_independent(pose->residue(i).seqpos())) {
                    rsdSeq = symm_info->bb_follows(pose->residue(i).seqpos());
                }

                if (symm_info->is_virtual(i)) { rsdSeq = 0; }
            }
		
            if ( rsdSeq ==0 ) continue;
            if ( pose->residue(rsdSeq).aa() == core::chemical::aa_vrt ) continue;


            // if topology allowed scoring, go
            Vector const & xyz( pose->residue( i ).atom( 2 ).xyz());
            Real depth = dot( xyz-center, normal );
            if( depth > -12 && depth < 12) { termini_pen++; }
        }

        termini_pen *= 50;
    }

} // scoring
} // membrane
} // core

#endif // INCLUDED_core_membrane_scoring_MembranePenalties_cc
