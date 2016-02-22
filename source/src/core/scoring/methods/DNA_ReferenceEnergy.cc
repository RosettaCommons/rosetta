// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/methods/DNA_ReferenceEnergy.cc
/// @brief  dna scoring
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/DNA_ReferenceEnergy.hh>
#include <core/scoring/methods/DNA_ReferenceEnergyCreator.hh>

// Package headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
//#include <core/scoring/EnergyGraph.hh>
//#include <core/scoring/Energies.hh>
//#include <core/scoring/dna/DNA_ReferencePotential.hh>
#include <core/scoring/dna/BasePartner.hh>
//#include <core/scoring/etable/count_pair/CountPair1BC4.hh>
#include <core/conformation/Residue.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>

// Project headers
//#include <core/pose/PoseCachedData.hh>
#include <core/pose/Pose.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
//#include <basic/options/util.hh> // HACK
#include <basic/Tracer.hh>
#include <utility/exit.hh> // HACK

// Utility headers


// C++


namespace core {
namespace scoring {
namespace methods {

using namespace dna; //////////////////// NOTE NOTE NOTE
static basic::Tracer TR("core.scoring.methods.DNA_ReferenceEnergy");



methods::EnergyMethodOP
DNA_ReferenceEnergyCreator::create_energy_method(
	EnergyMethodOptions const & options
) const {
	return EnergyMethodOP( new DNA_ReferenceEnergy( options ) );
}

ScoreTypes
DNA_ReferenceEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( dna_ref );
	return sts;
}





DNA_ReferenceEnergy::DNA_ReferenceEnergy( EnergyMethodOptions const & options ):
	parent( EnergyMethodCreatorOP( new DNA_ReferenceEnergyCreator ) )
{
	if ( !options.has_method_weights( dna_ref ) ) utility_exit_with_message( "dna_ref requires method weights!" );
	// order is aa,ac,ag,at,ca,cc,cg,ct,...
	utility::vector1< Real > wts( options.method_weights( dna_ref ) );
	base_step_reference_energies_.resize(4);
	runtime_assert( wts.size() == 16 );
	for ( Size i=1; i<= 4; ++i ) {
		base_step_reference_energies_[i].resize(4);
		for ( Size j=1; j<= 4; ++j ) {
			base_step_reference_energies_[i][j] = wts.front();
			wts.erase( wts.begin() );
		}
	}
	runtime_assert( wts.empty() );
	//
	// add_score_type( dna_ref );

	// {
	//  if ( options::option[ options::OptionKeys::dna::specificity::water_anchor_atom_weights ].user() ) {
	//   use_water_anchor_atom_weights_ = true;
	//   utility::vector1< std::string > wts
	//    ( options::option[ options::OptionKeys::dna::specificity::water_anchor_atom_weights ]() );
	//   for ( Size ii=1; ii<= wts.size(); ++ii ) {
	//    std::string const w( wts[ii] ); // looks like "t:O4:-0.1"
	//    if ( w[1] != ':' ) utility_exit_with_message("parse error: "+w);
	//    Size const pos( w.find_last_of( ':' ) );
	//    if ( !is_float( w.substr( pos+1 ) ) ) utility_exit_with_message("parse error: "+w);
	//    std::string const base_and_atom( w.substr( 0, pos ) );
	//    Real const weight( float_of( w.substr( pos+1 ) ) );
	//    water_anchor_atom_weights_[ base_and_atom ] = weight;
	//    TR.Trace << "water_anchor_atom_weights: " << base_and_atom << ' ' << weight << std::endl;
	//   }
	//  }
	// }
}

/// copy c-tor
DNA_ReferenceEnergy::DNA_ReferenceEnergy( DNA_ReferenceEnergy const & src ):
	parent( EnergyMethodCreatorOP( new DNA_ReferenceEnergyCreator ) ),
	base_step_reference_energies_( src.base_step_reference_energies_ )
	// use_water_anchor_atom_weights_( src.use_water_anchor_atom_weights_ ),
	// water_anchor_atom_weights_( src.water_anchor_atom_weights_ )
{
}


/// clone
EnergyMethodOP
DNA_ReferenceEnergy::clone() const
{
	return EnergyMethodOP( new DNA_ReferenceEnergy( *this ) );
}

/// are these really necessary??????????? move to scheme that doesnt depend on nbr calcn
///
// void
// DNA_ReferenceEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
// {
//  pose.update_residue_neighbors();
// }

// ///
// void
// DNA_ReferenceEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
// {
//  pose.update_residue_neighbors();
// }

// void
// DNA_ReferenceEnergy::setup_for_packing( pose::Pose & pose, pack::task::PackerTask const & ) const
// {
//  pose.update_residue_neighbors();
// }

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


/// same as dna::retrieve_base_partner_from_pose
inline
BasePartner const &
retrieve_base_partner_from_pose_inline( pose::Pose const & pose )
{
	assert( pose.data().has( core::pose::datacache::CacheableDataType::BASE_PARTNER ) );
	assert( dynamic_cast< BasePartner const *>( &( pose.data().get( core::pose::datacache::CacheableDataType::BASE_PARTNER ))));
	return ( static_cast< BasePartner const &>(    pose.data().get( core::pose::datacache::CacheableDataType::BASE_PARTNER )));
}


///
/// 04/29/11 === change to allow rsd1 and rsd2 in either sequence order
///
void
DNA_ReferenceEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	if ( rsd1.is_DNA() && rsd2.is_DNA() ) {
		Size const pos1( rsd1.seqpos() ), pos2( rsd2.seqpos() );
		if ( pos2 == pos1 + 1 && !rsd1.is_upper_terminus() &&
				count_pair_bs( pos1, pos2, retrieve_base_partner_from_pose_inline( pose ) ) ) {
			// this is not true during minimization, see MinimizationGraph/ScoreFunction
			//runtime_assert( std::fabs( emap[ dna_ref ] )<1e-3 );
			emap[ dna_ref ] += base_step_energy( rsd1.aa(), rsd2.aa() );
		} else if ( pos1 == pos2+1 && !rsd2.is_upper_terminus() &&
				count_pair_bs( pos2, pos1, retrieve_base_partner_from_pose_inline( pose ) ) ) {
			emap[ dna_ref ] += base_step_energy( rsd2.aa(), rsd1.aa() );
		}
	}
	// if ( use_water_anchor_atom_weights_ &&
	//    ( rsd1.is_DNA() || rsd2.is_DNA() ) &&
	//    ( rsd1.aa() == chemical::aa_h2o || rsd2.aa() == chemical::aa_h2o ) ) {
	//  using namespace pack::rotamer_set;
	//  Size waterpos,dnapos;
	//  if ( rsd1.aa() == chemical::aa_h2o ) { waterpos=rsd1.seqpos(); dnapos = rsd2.seqpos(); }
	//  else { waterpos=rsd2.seqpos(); dnapos = rsd1.seqpos(); }
	//  if ( has_water_packing_info( pose ) && get_water_packing_info( pose ).has( waterpos ) &&
	//     get_water_packing_info( pose )[ waterpos ].anchor_residue() == dnapos ) {
	//   char const name1( rsd1.is_DNA() ? rsd1.name1() : rsd2.name1() );
	//   std::string atom( get_water_packing_info(pose)[waterpos].anchor_atom_name());
	//   ObjexxFCL::strip_whitespace( atom );
	//   std::string const key( std::string(1,name1)+":"+atom );
	//   std::map< std::string, Real >::const_iterator it( water_anchor_atom_weights_.find( key ) );
	//   if ( it != water_anchor_atom_weights_.end() ) {
	//    emap[ dna_ref ] += it->second;
	//   }
	//  }
	// }
}



// void
// DNA_ReferenceEnergy::eval_atom_derivative(
//                      id::AtomID const &,
//                      pose::Pose const &,
//                      kinematics::DomainMap const &,
//                      ScoreFunction const &,
//                      EnergyMap const &,
//                      Vector &,
//                      Vector &
//                      ) const
// {
// }



/// @brief DNA_ReferenceEnergy distance cutoff
Distance
DNA_ReferenceEnergy::atomic_interaction_cutoff() const
{
	return 5.5; // -- temporary hack to allow us to use the standard neighbor array
}

/// @brief DNA_ReferenceEnergy
void
DNA_ReferenceEnergy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{
}



} // methods
} // scoring
} // core
