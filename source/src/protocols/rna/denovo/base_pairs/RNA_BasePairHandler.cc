// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/denovo/base_pairs/RNA_BasePairHandler.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/rna/denovo/base_pairs/RNA_BasePairHandler.hh>
#include <protocols/rna/denovo/setup/RNA_DeNovoParameters.hh>
#include <protocols/rna/denovo/libraries/BasePairStepLibrary.hh>
#include <protocols/rna/denovo/base_pairs/BasePairStep.hh>
#include <protocols/rna/util.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.rna.denovo.base_pairs.RNA_BasePairHandler" );

using namespace core;
using namespace core::pose::rna;
using namespace core::chemical::rna;
using namespace ObjexxFCL::format;
using namespace protocols::rna::denovo::setup;

namespace protocols {
namespace rna {
namespace denovo {
namespace base_pairs {

//Constructor
RNA_BasePairHandler::RNA_BasePairHandler( core::pose::Pose const & pose ) {
	core::kinematics::FoldTree const & f( pose.fold_tree() );
	for ( Size n = 1; n <= f.num_jump(); n++ ) {
		Size const jump_pos1( f.upstream_jump_residue( n ) );
		Size const jump_pos2( f.downstream_jump_residue( n ) );
		if ( pose.residue_type( jump_pos1 ).is_RNA() &&
				pose.residue_type( jump_pos2 ).is_RNA() ) {
			rna_pairing_list_.push_back( core::pose::rna::BasePair( jump_pos1, jump_pos2 ) /* note that edges are not specified*/ );
		}
	}

	for ( Size i = 1; i < pose.size(); i++ ) {
		if ( f.is_cutpoint( i ) ) {
			if ( pose.residue_type( i ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) {
				runtime_assert( pose.residue_type( i+1 ).has_variant_type( core::chemical::CUTPOINT_UPPER ) );
				// cutpoints_closed_.push_back( i );
			} else {
				cutpoints_open_.push_back( i );
			}
		}
	}
}

RNA_BasePairHandler::RNA_BasePairHandler( RNA_DeNovoParameters const & rna_params )
{
	rna_pairing_list_  = rna_params.rna_pairing_list();
	chain_connections_ = rna_params.chain_connections();
	cutpoints_open_    = rna_params.cutpoints_open();
}

//Destructor
RNA_BasePairHandler::~RNA_BasePairHandler()
{}

/////////////////////////////////////////////////////////
bool
RNA_BasePairHandler::check_base_pairs( pose::Pose & pose,
	protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map ) const
{
	using namespace core::pose::rna;
	static scoring::rna::RNA_LowResolutionPotential const rna_low_resolution_potential;

	for ( BasePair const & rna_pairing : rna_pairing_list_ ) {
		Size i( rna_pairing.res1() );
		Size j( rna_pairing.res2() );
		if ( i > j ) {
			i = rna_pairing.res2();
			j = rna_pairing.res1();
		}

		// Check for non-RNA residues
		if ( !pose.residue_type(i).is_RNA() ) continue;
		if ( !pose.residue_type(j).is_RNA() ) continue;

		// are these part of the pose that is actually being moved?
		if  ( !atom_level_domain_map->get( named_atom_id_to_atom_id( id::NamedAtomID( "C1'", i ), pose ) ) ) continue;
		if  ( !atom_level_domain_map->get( named_atom_id_to_atom_id( id::NamedAtomID( "C1'", j ), pose ) ) ) continue;

		if ( !rna_low_resolution_potential.check_forming_base_pair(pose,i,j) ) {
			TR << "MISSING BASE PAIR " << i << " " << j << std::endl;
			return false;
		}

		if ( rna_pairing.edge1() == WATSON_CRICK && rna_pairing.edge2() == WATSON_CRICK && rna_pairing.orientation() == ANTIPARALLEL ) {
			if ( is_cutpoint_open(pose, i) && is_cutpoint_open( pose, j-1) ) {
				if ( !rna_low_resolution_potential.check_clear_for_stacking( pose, i, +1 /* sign */) ) return false;
				if ( !rna_low_resolution_potential.check_clear_for_stacking( pose, j, -1 /* sign*/ ) ) return false;
			} else if ( is_cutpoint_open( pose, i-1 ) && is_cutpoint_open( pose, j) ) {
				if ( !rna_low_resolution_potential.check_clear_for_stacking( pose, i, -1 /* sign */) ) return false;
				if ( !rna_low_resolution_potential.check_clear_for_stacking( pose, j, +1 /* sign*/ ) ) return false;
			}
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////
void
RNA_BasePairHandler::setup_base_pair_constraints( core::pose::Pose & pose,
	protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
	Real const suppress_bp_constraint /* = 1.0 */ ) const
{
	using namespace core::id;

	utility::vector1< std::pair< Size, Size > > pairings;

	// Go through all pairings and define atom pair constraints that will bring together
	//  appropriate atoms on bases.
	for ( BasePair const & rna_pairing : rna_pairing_list_ ) {
		Size const & i( rna_pairing.res1() );
		Size const & j( rna_pairing.res2() );

		//Basic check that its canonical...
		if ( !( rna_pairing.edge1() == WATSON_CRICK && rna_pairing.edge2() == WATSON_CRICK && rna_pairing.orientation() == ANTIPARALLEL ) ) {
			TR.Debug<<  "skipping constraints for non-canonical base pair: " << I(3,i) << " " << I(3,j) << " " << rna_pairing.edge1() << " " << rna_pairing.edge2() << " " << rna_pairing.orientation() << std::endl;
			continue;
		}

		if ( !pose.residue_type(i).is_RNA() ) continue;
		if ( !pose.residue_type(j).is_RNA() ) continue;

		if ( !atom_level_domain_map->get( named_atom_id_to_atom_id( NamedAtomID( "C1'", i ), pose ) ) &&
				!atom_level_domain_map->get( named_atom_id_to_atom_id( NamedAtomID( "C1'", j ), pose ) ) ) continue; //assumed to be frozen, so no need to set up constraints?

		pairings.push_back( std::make_pair( i, j ) ) ;
	}

	// In RNA_ProtocolUtil.cc :
	protocols::rna::setup_base_pair_constraints( pose, pairings, suppress_bp_constraint );
}

/////////////////////////////////////////////////////////////////////
std::map< Size, Size >
RNA_BasePairHandler::connections() const
{
	std::map< Size, Size > connections_local;
	for ( BasePair const & rna_pairing : rna_pairing_list_ ) {
		connections_local[ rna_pairing.res1() ] = rna_pairing.res2();
		connections_local[ rna_pairing.res2() ] = rna_pairing.res1();
	}
	return connections_local;
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_BasePairHandler::figure_out_partner( std::map< Size, Size > & partner, bool const force_canonical ) const
{
	for ( BasePair const & rna_pairing : rna_pairing_list_ ) {
		Size i( rna_pairing.res1() );
		Size j( rna_pairing.res2() );

		bool const pair_is_canonical = rna_pairing.edge1() == WATSON_CRICK && rna_pairing.edge2() == WATSON_CRICK &&
			rna_pairing.orientation() == ANTIPARALLEL;
		if ( force_canonical && !pair_is_canonical ) continue;

		// if pairing is user-specified to be, say, H/S/A, we won't be able to handle it.
		bool const pair_is_ambiguous = rna_pairing.edge1() == ANY_BASE_EDGE && rna_pairing.edge2() == ANY_BASE_EDGE &&
			rna_pairing.orientation() == ANY_BASE_DOUBLET_ORIENTATION;
		//  TR << TR.Red << "PAIRING " << rna_pairing.res1() << " " << rna_pairing.res2() << " " << rna_pairing.edge1() << " " << rna_pairing.edge2() << " " << rna_pairing.orientation() << " --> "  << pair_is_ambiguous << std::endl;
		if ( !force_canonical && !pair_is_ambiguous ) continue;

		if ( partner.find( i ) != partner.end() )  {
			TR << "Warning: in base pair steps, already found a partner for " << i << " which is " << partner[ i ] << ", so not including additional pairing to " << j << std::endl;
		}
		if ( partner.find( j ) != partner.end() )  {
			TR << "Warning: in base pair steps, already found a partner for " << j << " which is " << partner[ j ] << ", so not including additional pairing to " << i << std::endl;
		}

		partner[ i ] = j;
		partner[ j ] = i;
	}
}

////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< BasePairStep >
RNA_BasePairHandler::get_base_pair_steps( bool const just_canonical ) const {

	// Parse base pair steps.
	std::map< Size, Size > partner;
	// Canonical pairs take precedence
	figure_out_partner( partner, true  /*force_canonical*/ );
	// Then any ambiguous pairs for residues that remain "partner free" (base edges=X, orientation=X)
	if ( !just_canonical ) figure_out_partner( partner, false /*force_canonical*/ );

	utility::vector1< BasePairStep > base_pair_steps;

	for ( auto const & elem : partner ) {
		Size const i = elem.first;
		if ( partner.find( i+1 ) == partner.end() || cutpoints_open_.has_value( i ) ) continue;

		Size const j = partner[ i+1 ];
		// In following, q = 0  means flush base pair steps. Strands are (i, i+1)  and (j, j+1 )
		// For q > 0, there are q unpaired residues ('bulges') in between j and j+q. Strands are (i, i+1) and ( j, j+1, ... j+q+1 ).
		for ( int q = 0; q <= libraries::MAX_BULGE_LENGTH; q++ ) {
			if ( j != partner[i] - q - 1 ) continue;

			bool has_cutpoint( false );
			for ( Size n = j; n < partner[ i ]; n++ ) {
				if ( cutpoints_open_.has_value( n ) ) has_cutpoint = true;
			}
			if ( has_cutpoint ) continue;
			if ( q == 0 && base_pair_steps.has_value( BasePairStep( j, j+1, i, i+1 )  ) ) { // don't put in the 'reverse' of the base pair step.
				continue;
			}
			base_pair_steps.push_back( BasePairStep( i, i+1, j, j+q+1 ) );
		}
	}

	return base_pair_steps;
}

////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< BasePairStep >
RNA_BasePairHandler::get_canonical_base_pair_steps() const {
	return get_base_pair_steps( true /* just_canonical*/ );
}

////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< BasePairStep >
RNA_BasePairHandler::get_noncanonical_base_pair_steps() const {
	utility::vector1< BasePairStep > canonical_base_pair_steps = get_base_pair_steps( true /* just_canonical*/ );
	utility::vector1< BasePairStep > general_base_pair_steps   = get_base_pair_steps( false /* just_canonical*/ );

	utility::vector1< BasePairStep > noncanonical_base_pair_steps;
	for ( Size n = 1; n <= general_base_pair_steps.size(); n++ ) {
		BasePairStep const & base_pair_step = general_base_pair_steps[ n ];
		if ( canonical_base_pair_steps.has_value( base_pair_step ) ) continue;
		noncanonical_base_pair_steps.push_back( base_pair_step );
	}

	return noncanonical_base_pair_steps;
}


//////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
RNA_BasePairHandler::get_stem_residues( core::pose::Pose const & pose ) const
{
	std::list< Size > in_stem;

	for ( BasePair const & rna_pairing : rna_pairing_list_ ) {
		if (  rna_pairing.edge1() == WATSON_CRICK &&
				rna_pairing.edge2() == WATSON_CRICK &&
				rna_pairing.orientation() == ANTIPARALLEL &&
				core::chemical::rna::possibly_canonical( pose.residue( rna_pairing.res1() ).aa(),
				pose.residue( rna_pairing.res2() ).aa() ) )  {
			in_stem.push_back( rna_pairing.res1() );
			in_stem.push_back( rna_pairing.res2() );
		}
	}

	in_stem.sort();
	in_stem.unique();

	utility::vector1< Size > stem_res;
	stem_res.assign( in_stem.begin(), in_stem.end() );

	return stem_res;
}



} //base_pairs
} //denovo
} //rna
} //protocols
