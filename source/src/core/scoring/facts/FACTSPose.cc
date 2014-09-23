// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:f;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
 // (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// @file:   core/scoring/facts/FACTSPotential.cc
// @brief:  The definitions of 3 classes of the FACTS algorithm resides here (see devel/khorvash/FACTSPotential.hh
// @author: Hahnbeom Park

// Unit headers
#include <core/scoring/facts/FACTSPotential.fwd.hh>
#include <core/scoring/facts/FACTSResidue.hh>
#include <core/scoring/facts/FACTSPose.hh>

// Project headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/RotamerSetCacheableDataType.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>

#include <basic/prof.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <math.h>
#include <stdio.h>
//#include <time>

static thread_local basic::Tracer TR( "core.scoring.FACTSPoseInfo" );

using namespace std;

namespace core {
namespace scoring {

/**************************************************************************************************/
/*                                                                                                */
/*    @breif: The class FACTSPoseInfo                                                             */
/*                                                                                                */
/**************************************************************************************************/

// Constructor
FACTSPoseInfo::FACTSPoseInfo() :
  context_derivative_empty_(true)
{}

// Constructor
FACTSPoseInfo::FACTSPoseInfo( FACTSPoseInfo const & src ) : CacheableData()
{
  Size const src_size( src.size() );
  
  residue_info_.resize( src_size );
  placeholder_residue_.resize( src_size );
  placeholder_info_.resize( src_size );
  
  for ( Size i=1; i<= src_size; ++i ) {
    residue_info_[i] = src.residue_info_[i]->clone();
    if ( src.placeholder_residue_[i] ) {
      placeholder_residue_[i] = src.placeholder_residue_[i]->clone();
      placeholder_info_[i] = src.placeholder_info_[i]->clone();
    } else {
      placeholder_residue_[i] = 0;
      placeholder_info_[i] = 0;
    }
  }
  being_packed_ = src.being_packed_;
  context_derivative_empty_ = src.context_derivative_empty_;
}
  
void FACTSPoseInfo::initialize( pose::Pose const & pose, FACTSRsdTypeMap &rsdtypemap )
{
  Size const nres( pose.total_residue() ); // maybe faster for symm if it is made symm-aware
  
  residue_info_.resize( nres, 0 );
  placeholder_residue_.resize( nres, 0 );
  placeholder_info_.resize( nres, 0 );
  
  for ( Size i=1; i<= nres; ++i ) {
    if ( !residue_info_[i] ) {
      residue_info_[i] = FACTSResidueInfoOP( new FACTSResidueInfo() );
      residue_info_[i]->set_enumeration_shell( true );
    }
    
    // Initialize only if the residue is in enumeration_shell
    // otherwise keep information stored previously
    if( residue_info_[i]->enumeration_shell() ){
      // initialize residuetypeinfo if it has not been
      core::chemical::ResidueType const &rsdtype = pose.residue(i).type();
      FACTSRsdTypeMap::const_iterator it = rsdtypemap.find( &rsdtype );

      if ( it == rsdtypemap.end() ) {
				TR << "Adding new FACTS residue type info: " << rsdtype.name() << std::endl;
				FACTSRsdTypeInfoOP rsdtypeinfo( new FACTSRsdTypeInfo );
				rsdtypeinfo->create_info( rsdtype );
				rsdtypemap[ &rsdtype ] = rsdtypeinfo;
				it = rsdtypemap.find( &rsdtype );
      }
      residue_info_[i]->initialize( pose.residue(i), it->second );
    }
  }
}
  
void FACTSPoseInfo::set_placeholder( Size const i, ResidueOP rsd, FACTSResidueInfoOP info )
{
  placeholder_residue_[ i ] = rsd;
  placeholder_info_[ i ] = info;
}

void FACTSPoseInfo::set_repack_list( utility::vector1< bool > const & repacking_residues )
{
  being_packed_.resize( size(), false );
  for ( Size i=1; i<= size(); ++i ) {
    being_packed_[i] = repacking_residues[ i ];
  }
}

bool FACTSPoseInfo::is_changed( pose::Pose const &pose ){
  for( Size ires = 1; ires <= pose.total_residue(); ++ ires ){
    
    Size const natom( pose.residue(ires).natoms() );
    utility::vector1<Vector> const facts_xyz = residue_info( ires ).xyz();
    
    if( natom != facts_xyz.size() )
      return true;
    
    for( Size iatm = 1; iatm <= natom; ++ iatm ){
      Vector const dxyz = facts_xyz[iatm] - pose.residue(ires).xyz(iatm);
      Real const d2 = dxyz.dot(dxyz);
      if( d2 > 1.0e-6 ) return true;
    }
    
  }
  return false;
}

// Store change in residue level
void
FACTSPoseInfo::update_enumeration_shell( pose::Pose const &pose,
																				 bool const enumerate_second_shell ){

	Energies const & energies( pose.energies() );
	EnergyGraph const & energy_graph( energies.energy_graph() );
	
	// use minimization graph instead!
	

	// Old way - use coordinate
	/*
	// First check change in coordinate in residue level
	for( Size ires = 1; ires <= pose.total_residue(); ++ ires ){
		
		FACTSResidueInfo & facts1( residue_info( ires ) );
		facts1.set_changed( false );
		facts1.set_enumeration_shell( false );
		
		Size const natom( pose.residue(ires).natoms() );
		utility::vector1<Vector> const facts_xyz = residue_info( ires ).xyz();
		
		// Check residue conformation change by looking at coordinate
		// Is there better way than this?
		if( natom != facts_xyz.size() ){
			facts1.set_changed( true );
			facts1.set_enumeration_shell( true );
			
		} else {
			for( Size iatm = 1; iatm <= natom; ++ iatm ){
				Vector const dxyz = facts_xyz[iatm] - pose.residue(ires).xyz(iatm);
				Real const d2 = dxyz.dot(dxyz);
				if( d2 > 1.0e-6 ){
					facts1.set_changed( true );
					facts1.set_enumeration_shell( true );
					break;
				}
			}
		}
	}
	*/

	// Decide whether to expand enumeration to the second shell:
	// say A-B-C where A,B,C are residues and A-B < 10 A, B-C < 10 A, A-C > 10 A
	// if A's conformation changes, it will affect B's context, which will again change B-C interaction energy
	// In order to calculate energy induced by A's change, 
	// one should enumerate over second shell also otherwise B-C interaction would have error

	if( !enumerate_second_shell ) return;

	// Then iter over neighbors to expand shell where its structure change can
	// affect context difference
	for( Size res1 = 1; res1 <= pose.total_residue(); ++res1 ) {
		FACTSResidueInfo & facts1( residue_info( res1 ) );

		// Propagate change info into its first neighbor shell
		if( facts1.changed() ){
			facts1.set_enumeration_shell( true );
			for ( graph::Graph::EdgeListConstIter
							iru  = energy_graph.get_node( res1 )->const_edge_list_begin(),
							irue = energy_graph.get_node( res1 )->const_edge_list_end();
						iru != irue; ++iru ) {
				Size const res2( (*iru)->get_other_ind( res1 ) );
				FACTSResidueInfo & facts2( residue_info( res2 ) );
				facts2.set_enumeration_shell( true );
			}
		}
	}

}

} // namespace scoring
} // namespace core
