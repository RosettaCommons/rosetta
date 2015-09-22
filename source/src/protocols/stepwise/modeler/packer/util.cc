// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/packer/util.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/packer/util.hh>
#include <protocols/stepwise/modeler/packer/StepWisePacker.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/RNA_ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <utility/tools/make_vector1.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.packer.util" );

using namespace core::scoring;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace packer {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
StepWisePackerOP
get_packer(
	core::scoring::ScoreFunctionCOP pack_scorefxn,
	utility::vector1< core::Size > const & moving_res_list,
	protocols::stepwise::modeler::options::StepWiseModelerOptionsCOP options ) {

	using namespace core::scoring;
	// may want to put in proper handling of moving_partition_res -- i.e., pack any residues that make new contacts:
	StepWisePackerOP stepwise_packer( new StepWisePacker( moving_res_list ) );
	stepwise_packer->set_scorefxn( pack_scorefxn );
	stepwise_packer->set_use_packer_instead_of_rotamer_trials( options->use_packer_instead_of_rotamer_trials() );
	stepwise_packer->set_allow_virtual_side_chains( options->allow_virtual_side_chains() );
	if ( !options->o2prime_legacy_mode() ) {
		stepwise_packer->set_allow_virtual_o2prime_hydrogens( options->allow_virtual_o2prime_hydrogens() );
	}
	stepwise_packer->set_pack_all_side_chains( options->global_optimize() );
	return stepwise_packer;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Soon can be smart about terminal phosphates and 2'-OH for RNA and side-chain atom for protein.
//  -- should see big code speedup
///////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< core::Size >
figure_out_working_interface_res( core::pose::Pose const & pose,
	core::Size const working_moving_res ) {
	return figure_out_working_interface_res( pose, utility::tools::make_vector1( working_moving_res ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
figure_out_working_interface_res( core::pose::Pose const & pose,
	utility::vector1< Size > const & working_moving_res_list ){

	utility::vector1< bool > at_interface( pose.total_residue(), false );
	// could make a map to vecs to save memory:
	utility::vector1< utility::vector1< bool > > checked_pair( pose.total_residue(), at_interface );

	for ( Size n = 1; n <= working_moving_res_list.size(); n++ ) {
		figure_out_working_interface_res( pose, working_moving_res_list[n],
			at_interface, checked_pair );
	}

	utility::vector1< Size > interface_res;
	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( at_interface[i] ) interface_res.push_back( i );
	}

	//  TR << TR.Magenta << make_tag_with_dashes( interface_res ) << TR.Reset << std::endl;
	return interface_res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
bool
check_o2prime_contact( pose::Pose const & pose, Size const i, Size const j ) {

	Real const dist_cutoff2 = ( 4.0 * 4.0 );
	Vector const & o2prime_xyz = pose.residue( i ).xyz(  pose.residue_type( i ).RNA_type().o2prime_index() );

	core::conformation::Residue const & nbr_rsd = pose.residue( j );
	for ( Size n = 1; n <= nbr_rsd.nheavyatoms(); n++ ) {
		if ( nbr_rsd.is_virtual( n ) ) continue;
		Vector const & nbr_xyz = nbr_rsd.xyz( n );
		Real const rsd_dist2 = ( nbr_xyz - o2prime_xyz ).length_squared();
		if ( rsd_dist2 <= dist_cutoff2 ) return true;
	}

	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
void
figure_out_working_interface_res( core::pose::Pose const & pose,
	core::Size const working_moving_res,
	utility::vector1< bool > & interface_res /* save work here */,
	utility::vector1< utility::vector1< bool > > & checked_pair /* save work here */ ) {

	EnergyGraph const & energy_graph( pose.energies().energy_graph() ); // note -- pose must have been scored already.

	utility::vector1< Size > const moving_partition_res = figure_out_moving_partition_res( pose, working_moving_res );

	for ( Size n = 1; n <= moving_partition_res.size(); n++ ) {

		Size const i = moving_partition_res[ n ];
		if ( pose.residue_type(i).has_variant_type( core::chemical::VIRTUAL_RESIDUE_VARIANT ) ) continue;

		for ( graph::Graph::EdgeListConstIter
				iter = energy_graph.get_node( i )->const_edge_list_begin();
				iter != energy_graph.get_node( i )->const_edge_list_end();
				++iter ) {
			Size const j( (*iter)->get_other_ind( i ) );
			if ( pose.residue_type(j).has_variant_type( core::chemical::VIRTUAL_RESIDUE_VARIANT ) ) continue;

			if ( checked_pair[ i ][ j ] ) continue;
			if ( interface_res[ i ] && interface_res[ j ] ) continue;

			checked_pair[ i ][ j ] = true; // don't need to check again.
			checked_pair[ j ][ i ] = true;

			// special geometry checks for 2'-OH
			if ( pose.residue_type( i ).is_RNA() && pose.residue_type( j ) .is_RNA() ) {
				if ( check_o2prime_contact( pose, i, j ) ) interface_res[ i ] = true;
				if ( check_o2prime_contact( pose, j, i ) ) interface_res[ j ] = true;
			} else {
				interface_res[ i ] = true;
				interface_res[ j ] = true;
			}
		} // j
	} // i

}

} //packer
} //modeler
} //stepwise
} //protocols
