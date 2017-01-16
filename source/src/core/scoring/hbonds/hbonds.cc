// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author


// Unit headers
#include <core/scoring/hbonds/hbonds.hh>

// Package headers
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>

// // Project headers
#include <core/conformation/Residue.hh>
#include <utility/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AtomType.hh>
#include <ObjexxFCL/format.hh>
#include <core/chemical/rna/RNA_Info.hh>

//pba membrane specific hbond
#include <core/scoring/Membrane_FAPotential.hh>
#include <core/scoring/Membrane_FAPotential.fwd.hh>
#include <core/scoring/MembranePotential.hh>

// Hydrogen bonding for membrane proteins - using membrane framework
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>

// Boost Headers
#include <boost/unordered_map.hpp>

//Standard I/O
//#include <stdio.h>


//#include <core/scoring/Energies.hh>

// // Numeric headers
// #include <numeric/numeric.functions.hh>

using namespace ObjexxFCL;
//pba
using namespace basic::options;
using namespace OptionKeys;

namespace core {
namespace scoring {
namespace hbonds {

static THREAD_LOCAL basic::Tracer tr( "core.scoring.hbonds.hbonds" );

/**
This routine fills an hbond-set with hbonds. All hbonds are included,
even ones which might be excluded later based on the backbone-hbond
exclusion.

WARNING WARNING WARNING
The pose must have an update energies object, eg it must be scored.
WARNING WARNING WARNING
**/


void
fill_hbond_set(
	pose::Pose const & pose,
	bool const calculate_derivative,
	HBondSet & hbond_set,
	bool const exclude_bb  /* default false */,
	bool const exclude_bsc /* default false */,
	bool const exclude_scb /* default false */,
	bool const exclude_sc  /* default false */
) {
	SSWeightParameters sswt;
	fill_hbond_set( pose, calculate_derivative, hbond_set, sswt, exclude_bb, exclude_bsc, exclude_scb, exclude_sc);
}


void
fill_hbond_set(
	pose::Pose const & pose,
	bool const calculate_derivative,
	HBondSet & hbond_set,
	SSWeightParameters const & sswt,
	bool const exclude_bb  /* default false */,
	bool const exclude_bsc /* default false */,
	bool const exclude_scb /* default false */,
	bool const exclude_sc  /* default false */
)
{
	debug_assert( pose.energies().residue_neighbors_updated() );

	// clear old data
	hbond_set.clear();
	HBondDatabase const & database( * HBondDatabase::get_database(hbond_set.hbond_options().params_database_tag()));

	// need to know which residues are neighbors
	// and what the neighbor-numbers are for each residue since some of the
	// weights are environment-dependent.
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );

	// loop over all nbr-pairs
	for ( Size res1 = 1; res1 <= pose.size(); ++res1 ) {
		int const nb1 = tenA_neighbor_graph.get_node( res1 )->num_neighbors_counting_self_static();
		conformation::Residue const & rsd1( pose.residue( res1 ) );

		for ( utility::graph::Graph::EdgeListConstIter
				iru = energy_graph.get_node(res1)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(res1)->const_upper_edge_list_end();
				iru != irue; ++iru ) {

			int const res2( (*iru)->get_second_node_ind() );

			conformation::Residue const & rsd2( pose.residue( res2 ) );
			if ( hbond_set.hbond_options().exclude_DNA_DNA() && rsd1.is_DNA() && rsd2.is_DNA() ) continue;

			int const nb2 = tenA_neighbor_graph.get_node( res2 )->num_neighbors_counting_self_static();

			// membrane specific hbond
			if ( hbond_set.hbond_options().Mbhbond() ) {
				identify_hbonds_1way_membrane(
					database,
					rsd1, rsd2, nb1, nb2, calculate_derivative,
					exclude_bb, exclude_bsc, exclude_scb, exclude_sc, hbond_set, pose );

				identify_hbonds_1way_membrane(
					database,
					rsd2, rsd1, nb2, nb1, calculate_derivative,
					exclude_bb, exclude_bsc, exclude_scb, exclude_sc, hbond_set, pose);
			} else {
				core::Real ssdep_weight = get_ssdep_weight(rsd1, rsd2, pose, sswt);

				identify_hbonds_1way(
					database,
					rsd1, rsd2, nb1, nb2, calculate_derivative,
					exclude_bb, exclude_bsc, exclude_scb, exclude_sc, hbond_set, ssdep_weight);

				identify_hbonds_1way(
					database,
					rsd2, rsd1, nb2, nb1, calculate_derivative,
					exclude_bb, exclude_bsc, exclude_scb, exclude_sc, hbond_set, ssdep_weight);

			}
		} // nbrs of res1
		if ( !hbond_set.hbond_options().exclude_self_hbonds() && hbond_set.hbond_options().exclude_DNA_DNA() && rsd1.is_DNA() ) {
			//pba membrane specific hbond
			if ( pose.conformation().is_membrane() || hbond_set.hbond_options().Mbhbond() ) {
				identify_hbonds_1way_membrane(
					database,
					rsd1, rsd1, nb1, nb1, calculate_derivative,
					exclude_bb, exclude_bsc, exclude_scb, exclude_sc, hbond_set, pose);
			} else {
				identify_hbonds_1way(
					database,
					rsd1, rsd1, nb1, nb1, calculate_derivative,
					exclude_bb, exclude_bsc, exclude_scb, exclude_sc, hbond_set);
			}
		}
	} // res1

	fill_intra_res_hbond_set( pose, calculate_derivative, hbond_set,
		exclude_bb, exclude_bsc, exclude_scb, exclude_sc ); // don't worry, will skip proteins.
}

///////////////////////////////////////////////////
bool
calculate_intra_res_hbonds( conformation::Residue const & rsd,
	HBondOptions const & options ) {
	if ( rsd.is_protein() && options.exclude_intra_res_protein() /* default true */ ) return false;
	if ( rsd.is_DNA() && options.exclude_DNA_DNA() /* default true */ ) return false;
	if ( rsd.is_RNA() && options.exclude_intra_res_RNA() /* default false */ ) return false;
	return true;
}

///////////////////////////////////////////////////
void
fill_intra_res_hbond_set(
	pose::Pose const & pose,
	bool const calculate_derivative,
	HBondSet & hbond_set,
	bool const exclude_bb  /* default false */,
	bool const exclude_bsc /* default false */,
	bool const exclude_scb /* default false */,
	bool const exclude_sc  /* default false */
) {
	HBondDatabase const & database( * HBondDatabase::get_database(hbond_set.hbond_options().params_database_tag()));
	TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );

	for ( Size res1 = 1; res1 <= pose.size(); ++res1 ) {
		conformation::Residue const & rsd1( pose.residue( res1 ) );
		if ( !calculate_intra_res_hbonds( rsd1, hbond_set.hbond_options() ) ) continue;
		Size const rsd_nb = tenA_neighbor_graph.get_node( rsd1.seqpos() )->num_neighbors_counting_self_static();
		identify_intra_res_hbonds( database, rsd1, rsd_nb, calculate_derivative, hbond_set,
			exclude_bb, exclude_bsc, exclude_scb, exclude_sc );
	}
}


core::Real
get_ssdep_weight(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	SSWeightParameters const & ssdep
) {
	if ( ! ssdep.ssdep_ ) return 1.0;

	// check if they are in same helix
	bool connected=true;
	Size hstart = std::min( rsd1.seqpos(), rsd2.seqpos()), hend=std::max( rsd1.seqpos(), rsd2.seqpos());
	for ( Size i=hstart; i<=hend && connected; ++i ) {
		connected = (pose.secstruct(i) == 'H');
	}
	if ( !connected ) return 1.0;
	Real ss_len_scalefactor = 1.0;

	while ( hstart >= 1 && pose.secstruct(hstart) == 'H' ) hstart--;
	while ( hend <= pose.size() && pose.secstruct(hend) == 'H' ) hend++;

	Size hlen = hend-hstart-1;
	if ( pose.secstruct(hstart) == 'H' ) hlen++;
	if ( pose.secstruct(hend) == 'H' ) hlen++;

	if ( hlen<=ssdep.len_l_ ) ss_len_scalefactor = ssdep.l_;
	else if ( hlen>=ssdep.len_h_ ) ss_len_scalefactor = ssdep.h_;
	else {
		Real m = (ssdep.h_-ssdep.l_)/(ssdep.len_h_-ssdep.len_l_);
		Real b = ssdep.l_-ssdep.len_l_*m;
		ss_len_scalefactor = m*hlen+b;
	}
	//std::cerr << "hlen = " << hlen << "scale = " << ss_len_scalefactor << std::endl;

	return ss_len_scalefactor;
}

///////////////////////////////////////////////////

void
fill_hbond_set_by_AHdist_threshold(
	core::pose::Pose const & pose,
	Real const AHdist_threshold,
	HBondSet & hbond_set
) {

	hbond_set.clear();
	HBondDatabase const & database( * HBondDatabase::get_database(hbond_set.hbond_options().params_database_tag()));

	// need to know which residues are neighbors
	// and what the neighbor-numbers are for each residue since some of the
	// weights are environment-dependent.
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	// It's possible to miss some relevant neighbors at 10A
	// 2*pose_max_nbr_radius(pose) gives me 12.2A for talaris2013 with fa_standard
	//TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
	TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );

	// loop over all nbr-pairs
	for ( Size res1 = 1; res1 <= pose.size(); ++res1 ) {
		core::conformation::Residue const & rsd1( pose.residue( res1 ) );
		int const nb1 = tenA_neighbor_graph.get_node( rsd1.seqpos() )->num_neighbors_counting_self_static();

		for ( utility::graph::Graph::EdgeListConstIter
				iru = energy_graph.get_node(res1)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(res1)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			int const res2( (*iru)->get_second_node_ind() );
			core::conformation::Residue const & rsd2(pose.residue(res2));
			int const nb2 = tenA_neighbor_graph.get_node(res2)->num_neighbors_counting_self_static();


			if ( hbond_set.hbond_options().exclude_DNA_DNA() &&
					rsd1.is_DNA() && rsd2.is_DNA() ) continue;

			identify_hbonds_1way_AHdist(database, rsd1, rsd2, nb1, nb2, AHdist_threshold, hbond_set);
			identify_hbonds_1way_AHdist(database, rsd2, rsd1, nb1, nb2, AHdist_threshold, hbond_set);
		}
	}

}


/// @details identify_hbonds_1way is overloaded to either add HBond objects to
/// an HBondSet or to accumulate energy into a EnergyMap
/// object.  This is done for performance reasons.  The allocation of
/// the temporary HBondSet on the heap causes a substatial slow down.
void
identify_hbonds_1way(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const evaluate_derivative,
	bool const exclude_bb,  /* exclude if acc=bb and don=bb */
	bool const exclude_bsc, /* exclude if acc=bb and don=sc */
	bool const exclude_scb, /* exclude if acc=sc and don=bb */
	bool const exclude_sc,  /* exclude if acc=sc and don=sc */
	// output
	HBondSet & hbond_set,
	Real ssdep_weight_factor
)
{
	debug_assert( don_rsd.seqpos() != acc_rsd.seqpos() );

	// <f1,f2> -- derivative vectors
	//std::pair< Vector, Vector > deriv( Vector(0.0), Vector(0.0 ) );
	HBondDerivs derivs;

	for ( Size const hatm : don_rsd.Hpos_polar() ) {

		Size const datm(don_rsd.atom_base(hatm));
		bool datm_is_bb = don_rsd.atom_is_backbone(datm);
		if ( datm_is_bb ) {
			if ( exclude_bb && exclude_scb ) continue;
		} else {
			if ( exclude_sc && exclude_bsc ) continue;
		}
		Vector const & hatm_xyz(don_rsd.atom(hatm).xyz());
		Vector const & datm_xyz(don_rsd.atom(datm).xyz());

		for ( Size const aatm : acc_rsd.accpt_pos() ) {

			if ( acc_rsd.atom_is_backbone(aatm) ) {
				if ( datm_is_bb ) {
					if ( exclude_bb ) continue;
				} else {
					if ( exclude_bsc ) continue;
				}
			} else {
				if ( datm_is_bb ) {
					if ( exclude_scb ) continue;
				} else {
					if ( exclude_sc ) continue;
				}
			}

			// rough filter for existence of hydrogen bond
			if ( hatm_xyz.distance_squared( acc_rsd.xyz( aatm ) ) > MAX_R2 ) continue;

			Real unweighted_energy( 0.0 );

			HBEvalTuple hbe_type( datm, don_rsd, aatm, acc_rsd);

			int const base ( acc_rsd.atom_base( aatm ) );
			int const base2( acc_rsd.abase2( aatm ) );
			debug_assert( base2 > 0 && base != base2 );

			//std::cout << "about to evaluate " << acc_rsd.atom_name(aatm ) << " " << acc_rsd.atom_name(base ) << acc_rsd.atom_name(base2) << std::endl;
			//std::cout << "coords " << acc_rsd.atom(aatm ) << "\t" << acc_rsd.atom(base ) << "\t" << acc_rsd.atom(base2) << std::endl;
			//std::cout << "seqpos " << don_rsd.seqpos() << "\t" << acc_rsd.seqpos() << std::endl;

			hb_energy_deriv( database, hbond_set.hbond_options(),
				hbe_type, datm_xyz, hatm_xyz,
				acc_rsd.atom(aatm ).xyz(),
				acc_rsd.atom(base ).xyz(),
				acc_rsd.atom(base2).xyz(),
				unweighted_energy, evaluate_derivative, derivs);

			if ( unweighted_energy >= MAX_HB_ENERGY ) continue;

			Real environmental_weight
				(!hbond_set.hbond_options().use_hb_env_dep() ? 1 :
				get_environment_dependent_weight(hbe_type, don_nb, acc_nb, hbond_set.hbond_options()));

			Real ssdep_weight = (get_hbond_weight_type(hbe_type.eval_type())==hbw_SR_BB) ?  ssdep_weight_factor : 1.0;

			//////
			// now we have identified a hbond -> append it into the hbond_set
			hbond_set.append_hbond( hatm, don_rsd, aatm, acc_rsd,
				hbe_type, unweighted_energy, environmental_weight*ssdep_weight, derivs );

			//////

		} // loop over acceptors
	} // loop over donors
}

void
identify_hbonds_1way(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const evaluate_derivative,
	bool const exclude_bb,  /* exclude if acc=bb and don=bb */
	bool const exclude_bsc, /* exclude if acc=bb and don=sc */
	bool const exclude_scb, /* exclude if acc=sc and don=bb */
	bool const exclude_sc,  /* exclude if acc=sc and don=sc */
	HBondOptions const & options,
	// output
	EnergyMap & emap,
	Real ssdep_weight_factor
)
{
	debug_assert( don_rsd.seqpos() != acc_rsd.seqpos() );

	// <f1,f2> -- derivative vectors
	//std::pair< Vector, Vector > deriv( Vector(0.0), Vector(0.0 ) );
	HBondDerivs derivs;

	for ( Size const hatm : don_rsd.Hpos_polar() ) {

		Size const datm(don_rsd.atom_base(hatm));
		bool datm_is_bb = don_rsd.atom_is_backbone(datm);
		if ( datm_is_bb ) {
			if ( exclude_bb && exclude_scb ) continue;
		} else {
			if ( exclude_sc && exclude_bsc ) continue;
		}
		Vector const & hatm_xyz(don_rsd.atom(hatm).xyz());
		Vector const & datm_xyz(don_rsd.atom(datm).xyz());

		for ( Size const aatm : acc_rsd.accpt_pos() ) {

			if ( acc_rsd.atom_is_backbone(aatm) ) {
				if ( datm_is_bb ) {
					if ( exclude_bb ) continue;
				} else {
					if ( exclude_bsc ) continue;
				}
			} else {
				if ( datm_is_bb ) {
					if ( exclude_scb ) continue;
				} else {
					if ( exclude_sc ) continue;
				}
			}

			// rough filter for existence of hydrogen bond
			if ( hatm_xyz.distance_squared( acc_rsd.xyz( aatm ) ) > MAX_R2 ) continue;

			Real unweighted_energy( 0.0 );

			HBEvalTuple hbe_type( datm, don_rsd, aatm, acc_rsd);

			int const base ( acc_rsd.atom_base( aatm ) );
			int const base2( acc_rsd.abase2( aatm ) );
			debug_assert( base2 > 0 && base != base2 );

			//std::cout << "about to evaluate " << acc_rsd.atom_name(aatm ) << " " << acc_rsd.atom_name(base ) << acc_rsd.atom_name(base2) << std::endl;
			//std::cout << "coords " << acc_rsd.atom(aatm ) << "\t" << acc_rsd.atom(base ) << "\t" << acc_rsd.atom(base2) << std::endl;
			//std::cout << "seqpos " << don_rsd.seqpos() << "\t" << acc_rsd.seqpos() << std::endl;

			hb_energy_deriv( database, options, hbe_type, datm_xyz, hatm_xyz,
				acc_rsd.atom(aatm ).xyz(),
				acc_rsd.atom(base ).xyz(),
				acc_rsd.atom(base2).xyz(),
				unweighted_energy, evaluate_derivative, derivs);

			if ( unweighted_energy >= MAX_HB_ENERGY ) continue;

			Real environmental_weight
				(!options.use_hb_env_dep() ? 1 :
				get_environment_dependent_weight(hbe_type, don_nb, acc_nb, options));

			////////
			// now we have identified an hbond -> accumulate its energy

			Real hbE = unweighted_energy /*raw energy*/ * environmental_weight /*env-dep-wt*/;
			emap[hbond] += hbE;
			switch(get_hbond_weight_type(hbe_type.eval_type())){
			case hbw_NONE:
			case hbw_SR_BB :
				emap[hbond_sr_bb] += ssdep_weight_factor*hbE; break;
			case hbw_LR_BB :
				emap[hbond_lr_bb] += hbE; break;
			case hbw_SR_BB_SC :
				//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
				emap[hbond_bb_sc] += hbE;
				emap[hbond_sr_bb_sc] += hbE; break;
			case hbw_LR_BB_SC :
				//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
				emap[hbond_bb_sc] += hbE;
				emap[hbond_lr_bb_sc] += hbE; break;
			case hbw_SC :
				emap[hbond_sc] += hbE; break;
			default :
				tr << "Warning: energy from unexpected HB type ignored "
					<< hbe_type.eval_type() << std::endl;
				runtime_assert(false);
				break;
			}
			/////////

		} // loop over acceptors
	} // loop over donors
}

void
identify_hbonds_1way(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const evaluate_derivative,
	bool const exclude_bb,  /* exclude if acc=bb and don=bb */
	bool const exclude_bsc, /* exclude if acc=bb and don=sc */
	bool const exclude_scb, /* exclude if acc=sc and don=bb */
	bool const exclude_sc,  /* exclude if acc=sc and don=sc */
	HBondOptions const & options,
	// output
	EnergyMap & emap,
	boost::unordered_map<core::Size, core::Size> & num_hbonds,
	Real ssdep_weight_factor
)
{
	debug_assert( don_rsd.seqpos() != acc_rsd.seqpos() );

	// <f1,f2> -- derivative vectors
	//std::pair< Vector, Vector > deriv( Vector(0.0), Vector(0.0 ) );
	HBondDerivs derivs;

	for ( Size const hatm : don_rsd.Hpos_polar() ) {

		Size const datm(don_rsd.atom_base(hatm));
		bool datm_is_bb = don_rsd.atom_is_backbone(datm);

		if ( datm_is_bb ) {
			if ( exclude_bb && exclude_scb ) continue;
		} else {
			if ( exclude_sc && exclude_bsc ) continue;
		}
		Vector const & hatm_xyz(don_rsd.atom(hatm).xyz());
		Vector const & datm_xyz(don_rsd.atom(datm).xyz());

		for ( Size const aatm : acc_rsd.accpt_pos() ) {

			if ( acc_rsd.atom_is_backbone(aatm) ) {
				if ( datm_is_bb ) {
					if ( exclude_bb ) continue;
				} else {
					if ( exclude_bsc ) continue;
				}
			} else {
				if ( datm_is_bb ) {
					if ( exclude_scb ) continue;
				} else {
					if ( exclude_sc ) continue;
				}
			}

			// rough filter for existence of hydrogen bond
			if ( hatm_xyz.distance_squared( acc_rsd.xyz( aatm ) ) > MAX_R2 ) continue;

			Real unweighted_energy( 0.0 );

			HBEvalTuple hbe_type( datm, don_rsd, aatm, acc_rsd);

			int const base ( acc_rsd.atom_base( aatm ) );
			int const base2( acc_rsd.abase2( aatm ) );
			debug_assert( base2 > 0 && base != base2 );

			//std::cout << "about to evaluate " << acc_rsd.atom_name(aatm ) << " " << acc_rsd.atom_name(base ) << acc_rsd.atom_name(base2) << std::endl;
			//std::cout << "coords " << acc_rsd.atom(aatm ) << "\t" << acc_rsd.atom(base ) << "\t" << acc_rsd.atom(base2) << std::endl;
			//std::cout << "seqpos " << don_rsd.seqpos() << "\t" << acc_rsd.seqpos() << std::endl;

			hb_energy_deriv( database, options, hbe_type, datm_xyz, hatm_xyz,
				acc_rsd.atom(aatm ).xyz(),
				acc_rsd.atom(base ).xyz(),
				acc_rsd.atom(base2).xyz(),
				unweighted_energy, evaluate_derivative, derivs);

			if ( unweighted_energy >= MAX_HB_ENERGY ) continue;
			//std::cout << std::endl << "Found a hydrogen bond" << std::endl;
			num_hbonds[don_rsd.seqpos()]++;
			num_hbonds[acc_rsd.seqpos()]++;

			Real environmental_weight
				(!options.use_hb_env_dep() ? 1 :
				get_environment_dependent_weight(hbe_type, don_nb, acc_nb, options));

			////////
			// now we have identified an hbond -> accumulate its energy

			Real hbE = unweighted_energy /*raw energy*/ * environmental_weight /*env-dep-wt*/;
			emap[hbond] += hbE;
			switch(get_hbond_weight_type(hbe_type.eval_type())){
			case hbw_NONE:
			case hbw_SR_BB :
				emap[hbond_sr_bb] += ssdep_weight_factor*hbE; break;
			case hbw_LR_BB :
				emap[hbond_lr_bb] += hbE; break;
			case hbw_SR_BB_SC :
				//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
				emap[hbond_bb_sc] += hbE;
				emap[hbond_sr_bb_sc] += hbE; break;
			case hbw_LR_BB_SC :
				//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
				emap[hbond_bb_sc] += hbE;
				emap[hbond_lr_bb_sc] += hbE; break;
			case hbw_SC :
				emap[hbond_sc] += hbE; break;
			default :
				tr << "Warning: energy from unexpected HB type ignored "
					<< hbe_type.eval_type() << std::endl;
				runtime_assert(false);
				break;
			}
			/////////

		} // loop over acceptors
	} // loop over donors
}


void
identify_hbonds_1way_AHdist(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	Real const AHdist_threshold,
	HBondSet & hbond_set
)
{
	for ( Size const hatm : don_rsd.Hpos_polar() ) {
		Size const datm(don_rsd.atom_base(hatm));
		Vector const & hatm_xyz(don_rsd.atom(hatm).xyz());
		Vector const & datm_xyz(don_rsd.atom(datm).xyz());

		for ( Size const aatm : acc_rsd.accpt_pos() ) {
			if ( hatm_xyz.distance( acc_rsd.xyz( aatm )) > AHdist_threshold ) continue;

			HBEvalTuple hbe_type(datm, don_rsd, aatm, acc_rsd);
			int const base ( acc_rsd.atom_base( aatm ) );
			int const base2( acc_rsd.abase2( aatm ) );

			Real unweighted_energy( 0.0 );
			hb_energy_deriv(database, hbond_set.hbond_options(),
				hbe_type, datm_xyz, hatm_xyz,
				acc_rsd.atom(aatm ).xyz(),
				acc_rsd.atom(base ).xyz(),
				acc_rsd.atom(base2).xyz(),
				unweighted_energy, false, DUMMY_DERIVS);

			Real environmental_weight
				(!hbond_set.hbond_options().use_hb_env_dep() ? 1 :
				get_environment_dependent_weight(hbe_type, don_nb, acc_nb, hbond_set.hbond_options()));

			hbond_set.append_hbond(
				hatm, don_rsd, aatm, acc_rsd, hbe_type, unweighted_energy, environmental_weight, DUMMY_DERIVS );
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
identify_intra_res_hbonds(
	HBondDatabase const & database,
	conformation::Residue const & rsd,
	Size const rsd_nb,
	bool const evaluate_derivative,
	HBondSet & hbond_set,
	bool const exclude_bb,  /* exclude if acc=bb and don=bb */
	bool const exclude_bsc, /* exclude if acc=bb and don=sc */
	bool const exclude_scb, /* exclude if acc=sc and don=bb */
	bool const exclude_sc  /* exclude if acc=sc and don=sc */
)
{
	runtime_assert_string_msg( calculate_intra_res_hbonds( rsd, hbond_set.hbond_options() ), "Error in core::scoring::hbonds::identify_intra_res_hbonds(): This function was called, but the energy options are set to exclude intra-residue hydrogen bonds." );

	// <f1,f2> -- derivative vectors
	HBondDerivs derivs;

	for ( Size const hatm : rsd.Hpos_polar() ) {
		Size const datm(rsd.atom_base(hatm));
		bool datm_is_bb = rsd.atom_is_backbone(datm);
		if ( datm_is_bb ) {
			if ( exclude_bb && exclude_scb ) continue;
		} else {
			if ( exclude_sc && exclude_bsc ) continue;
		}

		for ( Size const aatm : rsd.accpt_pos() ) {
			int const base ( rsd.atom_base( aatm ) ); //"base" does not refer to RNA base here!
			int const base2( rsd.abase2( aatm ) );
			runtime_assert( base2 > 0 && base != base2 );

			if ( rsd.atom_is_backbone(aatm) ) {
				if ( datm_is_bb ) {
					if ( exclude_bb ) continue;
				} else {
					if ( exclude_bsc ) continue;
				}
			} else {
				if ( datm_is_bb ) {
					if ( exclude_scb ) continue;
				} else {
					if ( exclude_sc ) continue;
				}
			}

			if ( rsd.path_distance( aatm, datm ) < 4 ) continue; //utility_exit_with_message("rsd.path_distance(aatm, datm) < 4"); //consistency check

			Vector const & hatm_xyz( rsd.atom(hatm).xyz());
			Vector const & datm_xyz( rsd.atom(datm).xyz());

			// rough filter for existence of hydrogen bond
			if ( hatm_xyz.distance_squared( rsd.xyz( aatm ) ) > MAX_R2 ) continue;

			Real unweighted_energy( 0.0 );

			HBEvalTuple hbe_type( datm, rsd, aatm, rsd);

			hb_energy_deriv( database, hbond_set.hbond_options(),
				hbe_type, datm_xyz, hatm_xyz,
				rsd.atom(aatm ).xyz(),
				rsd.atom(base ).xyz(),
				rsd.atom(base2).xyz(),
				unweighted_energy, evaluate_derivative, derivs);

			if ( unweighted_energy >= MAX_HB_ENERGY ) continue;

			Real environmental_weight
				(!hbond_set.hbond_options().use_hb_env_dep() ? 1 :
				get_environment_dependent_weight(hbe_type, rsd_nb, rsd_nb, hbond_set.hbond_options()));

			// now we have identified a hbond -> put it into the hbond_set//////
			hbond_set.append_hbond( hatm, rsd, aatm, rsd, hbe_type, unweighted_energy, environmental_weight, derivs );

		} // loop over donors
	} // loop over acceptors
}


// Note: in other places calculating hbond_1way, code is duplicated for performance reasons.
//  This loop is not going to be time-critical, so decided to wrap it around the HBondSet setting code,
//  to avoid copying code yet again. --rhiju
void
identify_intra_res_hbonds(
	HBondDatabase const & database,
	conformation::Residue const & rsd,
	Size const rsd_nb,
	HBondOptions const & options,
	EnergyMap & emap){

	runtime_assert( calculate_intra_res_hbonds( rsd, options ) );
	HBondSet hbond_set( options );
	identify_intra_res_hbonds( database, rsd, rsd_nb, false /*evaluate_derivative*/, hbond_set );
	for ( Size n = 1; n <= hbond_set.nhbonds(); n++ ) {
		HBond const & hbond_ = hbond_set.hbond( n );
		Real const weighted_energy = hbond_.energy() * hbond_.weight();
		if ( options.put_intra_into_total() ) {
			emap[ hbond ]       += weighted_energy;
		} else {
			emap[ hbond_intra ] += weighted_energy;
		}
	}
}

//mjo this should be the only way to assign hbond energies.  If you
//feel the need to collect the energies from some of the bonds,
//consider filtering the hbond set instead.  Don't assign energies for
//"special cases" with other functions, add cases to the HBEvalType
//instead.

void
get_hbond_energies(
	HBondSet const & hbond_set,
	EnergyMap & emap
)
{
	for ( Size i = 1; i <= hbond_set.nhbonds(); ++i ) {
		if ( !hbond_set.allow_hbond(i) ) continue;

		HBond const & hbond_(hbond_set.hbond(i));
		HBEvalType const hbe_type = hbond_.eval_type();
		if ( hbe_type == hbe_UNKNOWN ) {
			tr.Error << "ERROR: Unknown HBEvalType for " << hbond_.eval_tuple() << std::endl;
			utility_exit_with_message("Can't get energy for hbond interactions.");
		}

		Real hbE = hbond_.energy() /*raw energy*/ * hbond_.weight() /*env-dep-wt*/;
		emap[hbond] += hbE;
		switch(get_hbond_weight_type(hbe_type)){
		case hbw_NONE:
		case hbw_SR_BB :
			emap[hbond_sr_bb] += hbE; break;
		case hbw_LR_BB :
			emap[hbond_lr_bb] += hbE; break;
		case hbw_SR_BB_SC :
			//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
			emap[hbond_bb_sc] += hbE;
			emap[hbond_sr_bb_sc] += hbE; break;
		case hbw_LR_BB_SC :
			//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
			emap[hbond_bb_sc] += hbE;
			emap[hbond_lr_bb_sc] += hbE; break;
		case hbw_SC :
			emap[hbond_sc] += hbE; break;
		default :
			tr << "Warning: energy from unexpected HB type ignored "
				<< hbond << std::endl;
			runtime_assert(false);
			break;
		}
	}
}

Real
hb_eval_type_weight(
	HBEvalType const & hbe_type,
	EnergyMap const & weights,
	bool const intra_res /*false*/,
	bool const put_intra_into_total /*true*/ )
{
	Real weight(0.0);

	if ( intra_res ) {
		if ( put_intra_into_total ) {
			weight += weights[hbond]; // typically zero if other weights are non-zero.
		} else {
			weight += weights[hbond_intra];
		}
	} else {

		weight += weights[hbond]; // typically zero if other weights are non-zero.

		switch(get_hbond_weight_type(hbe_type)){
		case hbw_NONE:
		case hbw_SR_BB :
			weight += weights[hbond_sr_bb]; break;
		case hbw_LR_BB :
			weight += weights[hbond_lr_bb]; break;
		case hbw_SR_BB_SC :
			//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
			weight += weights[hbond_bb_sc];
			weight += weights[hbond_sr_bb_sc]; break;
		case hbw_LR_BB_SC :
			//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
			weight += weights[hbond_bb_sc];
			weight += weights[hbond_lr_bb_sc]; break;
		case hbw_SC :
			weight += weights[hbond_sc]; break;
		default :
			tr << "Warning: Unexpected HBondWeightType " << get_hbond_weight_type(hbe_type) << std::endl;
			runtime_assert(false);
			break;
		}
	}

	return weight;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// functions for environment-dependent weighting
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

////////////////////////
// lin jiang's approach:

//JSS note: this is half of the weight from one atom;
// the burial weight is the sum from donor and acceptor.
inline core::Real
burial_weight(int const nb)
{
	if ( nb < 7 ) return 0.1;
	if ( nb > 24 ) return 0.5;
	return (nb-2.75)*(0.5/21.25);
}

core::Real
hb_env_dep_burial_lin(int const nb1, int const nb2)
{
	return (burial_weight(nb1) + burial_weight(nb2));
}

////////////////////////////
// tanja kortemme's approach

void
burial3class_weight_initializer( FArray2D_double & burial )
{
	burial( 1, 1 ) = 0.2 ; burial( 1, 2 ) = 0.2 ; burial( 1, 3 ) = 0.55;
	burial( 2, 1 ) = 0.2 ; burial( 2, 2 ) = 0.55; burial( 2, 3 ) = 1.0 ;
	burial( 3, 1 ) = 0.55; burial( 3, 2 ) = 1.0 ; burial( 3, 3 ) = 1.0 ;
}

inline int
get_burial_3(
	int const neighbors,
	int const threshold_1,
	int const threshold_3
)
{
	//tk get burial measure, three possible classes:
	//tk 1: exposed, 2: intermediate, 3:buried
	if ( neighbors > threshold_1 ) {
		if ( neighbors >= threshold_3 ) return 3;
		else return 2;
	} else return 1;
}

core::Real
hb_env_dep_burial_tk(int const nb1, int const nb2)
{
	//tk assign weight based on CB neighbors of two interacting residues

	// local
	int const exposed_threshold = { 10 };
	int const buried_threshold = { 20 };

	static FArray2D_double const burial3class_weight
		( 3, 3, burial3class_weight_initializer );

	return burial3class_weight(
		get_burial_3( nb1, exposed_threshold, buried_threshold ),
		get_burial_3( nb2, exposed_threshold, buried_threshold ) );
}

///////////////////////////////////////////////////////////////////////////////
Real
get_environment_dependent_weight(
	HBEvalTuple const & hbe_type,
	int const don_nb,
	int const acc_nb,
	HBondOptions const & options
)
{
	Real weight( 1.0 );

	// mjo why is this only applied to side chains!!!!
	if ( hbe_is_SC_type(hbe_type.eval_type()) ) {
		if ( options.smooth_hb_env_dep() ) {
			weight = hb_env_dep_burial_lin( acc_nb, don_nb );
		} else {
			weight = hb_env_dep_burial_tk( acc_nb, don_nb );
		}
	}

	// std::cout << "HB_ENV_WEIGHT: " << weight << std::endl;
	return weight;
}

///////////////////////////////////////////////////////////////////////////////

bool
nonzero_hbond_weight( ScoreFunction const & scorefxn )
{
	return ( scorefxn.has_nonzero_weight( hbond_lr_bb ) ||
		scorefxn.has_nonzero_weight( hbond_sr_bb ) ||
		scorefxn.has_nonzero_weight( hbond_bb_sc ) ||
		scorefxn.has_nonzero_weight( hbond_sr_bb_sc ) ||
		scorefxn.has_nonzero_weight( hbond_lr_bb_sc ) ||
		scorefxn.has_nonzero_weight( hbond_sc ) ||
		scorefxn.has_nonzero_weight( hbond_intra ) ||
		scorefxn.has_nonzero_weight( hbond ) );
}


///////////////////////////////////////////////////////////////////////////////
void
identify_hbonds_1way_membrane(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const evaluate_derivative,
	bool const exclude_bb,  /* exclude if acc=bb and don=bb */
	bool const exclude_bsc, /* exclude if acc=bb and don=sc */
	bool const exclude_scb, /* exclude if acc=sc and don=bb */
	bool const exclude_sc,  /* exclude if acc=sc and don=sc */
	// output
	HBondSet & hbond_set,
	pose::Pose const & pose
){
	assert( don_rsd.seqpos() != acc_rsd.seqpos() );

	HBondDerivs derivs;

	for ( chemical::AtomIndices::const_iterator
			hnum = don_rsd.Hpos_polar().begin(), hnume = don_rsd.Hpos_polar().end();
			hnum != hnume; ++hnum ) {
		Size const hatm( *hnum );
		Size const datm(don_rsd.atom_base(hatm));
		bool datm_is_bb = don_rsd.atom_is_backbone(datm);

		if ( datm_is_bb ) {
			if ( exclude_bb && exclude_scb ) continue;
		} else {
			if ( exclude_sc && exclude_bsc ) continue;
		}
		Vector const & hatm_xyz(don_rsd.atom(hatm).xyz());
		Vector const & datm_xyz(don_rsd.atom(datm).xyz());

		for ( chemical::AtomIndices::const_iterator
				anum = acc_rsd.accpt_pos().begin(), anume = acc_rsd.accpt_pos().end();
				anum != anume; ++anum ) {
			Size const aatm( *anum );
			if ( acc_rsd.atom_is_backbone(aatm) ) {
				if ( datm_is_bb ) {
					if ( exclude_bb ) continue;
				} else {
					if ( exclude_bsc ) continue;
				}
			} else {
				if ( datm_is_bb ) {
					if ( exclude_scb ) continue;
				} else {
					if ( exclude_sc ) continue;
				}
			}

			// rough filter for existance of hydrogen bond
			if ( hatm_xyz.distance_squared( acc_rsd.xyz( aatm ) ) > MAX_R2 ) continue;

			Real unweighted_energy( 0.0 );

			HBEvalTuple hbe_type( datm, don_rsd, aatm, acc_rsd);

			int const base ( acc_rsd.atom_base( aatm ) );
			int const base2( acc_rsd.abase2( aatm ) );
			assert( base2 > 0 && base != base2 );

			hb_energy_deriv( database, hbond_set.hbond_options(),
				hbe_type, datm_xyz, hatm_xyz,
				acc_rsd.atom(aatm ).xyz(),
				acc_rsd.atom(base ).xyz(),
				acc_rsd.atom(base2).xyz(),
				unweighted_energy, evaluate_derivative, derivs);

			if ( unweighted_energy >= MAX_HB_ENERGY ) continue;

			//pba membrane depth dependent correction to the environmental_weight
			Real environmental_weight(
				get_membrane_depth_dependent_weight(pose, don_nb, acc_nb, hatm_xyz,
				acc_rsd.atom(aatm ).xyz()));

			/// Add Membrane hbonds to the hydrogen bonds set
			hbond_set.append_hbond( hatm, don_rsd, aatm, acc_rsd,
				hbe_type, unweighted_energy, environmental_weight, derivs );

		} // loop over acceptors
	} // loop over donors
}

void
identify_hbonds_1way_membrane(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const evaluate_derivative,
	bool const exclude_bb,  /* exclude if acc=bb and don=bb */
	bool const exclude_bsc, /* exclude if acc=bb and don=sc */
	bool const exclude_scb, /* exclude if acc=sc and don=bb */
	bool const exclude_sc,  /* exclude if acc=sc and don=sc */
	HBondOptions const & options,
	EnergyMap & emap,
	pose::Pose const & pose
)
{
	assert( don_rsd.seqpos() != acc_rsd.seqpos() );

	HBondDerivs derivs;

	for ( chemical::AtomIndices::const_iterator
			hnum = don_rsd.Hpos_polar().begin(), hnume = don_rsd.Hpos_polar().end();
			hnum != hnume; ++hnum ) {
		Size const hatm( *hnum );
		Size const datm(don_rsd.atom_base(hatm));
		bool datm_is_bb = don_rsd.atom_is_backbone(datm);

		if ( datm_is_bb ) {
			if ( exclude_bb && exclude_scb ) continue;
		} else {
			if ( exclude_sc && exclude_bsc ) continue;
		}
		Vector const & hatm_xyz(don_rsd.atom(hatm).xyz());
		Vector const & datm_xyz(don_rsd.atom(datm).xyz());

		for ( chemical::AtomIndices::const_iterator
				anum = acc_rsd.accpt_pos().begin(), anume = acc_rsd.accpt_pos().end();
				anum != anume; ++anum ) {
			Size const aatm( *anum );
			if ( acc_rsd.atom_is_backbone(aatm) ) {
				if ( datm_is_bb ) {
					if ( exclude_bb ) continue;
				} else {
					if ( exclude_bsc ) continue;
				}
			} else {
				if ( datm_is_bb ) {
					if ( exclude_scb ) continue;
				} else {
					if ( exclude_sc ) continue;
				}
			}

			// rough filter for existance of hydrogen bond
			if ( hatm_xyz.distance_squared( acc_rsd.xyz( aatm ) ) > MAX_R2 ) continue;

			Real unweighted_energy( 0.0 );

			HBEvalTuple hbe_type( datm, don_rsd, aatm, acc_rsd);

			int const base ( acc_rsd.atom_base( aatm ) );
			int const base2( acc_rsd.abase2( aatm ) );
			assert( base2 > 0 && base != base2 );

			hb_energy_deriv( database, options,
				hbe_type, datm_xyz, hatm_xyz,
				acc_rsd.atom(aatm ).xyz(),
				acc_rsd.atom(base ).xyz(),
				acc_rsd.atom(base2).xyz(),
				unweighted_energy, evaluate_derivative, derivs);

			if ( unweighted_energy >= MAX_HB_ENERGY ) continue;

			//pba membrane depth dependent weight

			Real environmental_weight(
				get_membrane_depth_dependent_weight(pose, don_nb, acc_nb, hatm_xyz,
				acc_rsd.atom(aatm ).xyz()));

			////////
			// now we have identified an hbond -> accumulate its energy

			Real hbE = unweighted_energy /*raw energy*/ * environmental_weight /*env-dep-wt*/;
			emap[hbond] += hbE;
			switch(get_hbond_weight_type(hbe_type.eval_type())){
			case hbw_NONE:
			case hbw_SR_BB :
				emap[hbond_sr_bb] += hbE; break;
			case hbw_LR_BB :
				emap[hbond_lr_bb] += hbE; break;
			case hbw_SR_BB_SC :
				//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
				emap[hbond_bb_sc] += hbE;
				emap[hbond_sr_bb_sc] += hbE; break;
			case hbw_LR_BB_SC :
				//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
				emap[hbond_bb_sc] += hbE;
				emap[hbond_lr_bb_sc] += hbE; break;
			case hbw_SC :
				emap[hbond_sc] += hbE; break;
			default :
				tr << "Warning: energy from unexpected HB type ignored "
					<< hbe_type.eval_type() << std::endl;
				runtime_assert(false);
				break;
			}
			/////////

		} // loop over acceptors
	} // loop over donors
}

Real
get_membrane_depth_dependent_weight(
	pose::Pose const & pose,
	int const don_nb,
	int const acc_nb,
	Vector const & Hxyz, // proton
	Vector const & Axyz  // acceptor
)
{
	Real wat_weight(1.0), memb_weight(1.0), total_weight(1.0);

	// water phase smooth_hb_env_dep
	wat_weight = hb_env_dep_burial_lin( acc_nb, don_nb );


	// Initialize to defaults (taken from Barth, 2007)
	Vector normal( 0, 0, 1 );
	Vector center( 0, 0, 0 );
	Real thickness( 15.0 );
	Real steepness( 10.0 );


	if ( pose.conformation().is_membrane() ) {

		center = pose.conformation().membrane_info()->membrane_center(pose.conformation());
		normal = pose.conformation().membrane_info()->membrane_normal(pose.conformation());
		thickness = pose.conformation().membrane_info()->membrane_thickness();
		steepness = pose.conformation().membrane_info()->membrane_steepness();

	} else {

		normal = MembraneEmbed_from_pose( pose ).normal();
		center = MembraneEmbed_from_pose( pose ).center();
		thickness = Membrane_FAEmbed_from_pose( pose ).thickness();
		steepness = Membrane_FAEmbed_from_pose( pose ).steepness();
	}

	// Hdonor depth
	Real fa_depth_H = dot(Hxyz-center, normal); // non consistent z_position
	Real internal_product = std::abs(fa_depth_H);
	Real z = internal_product;
	z /= thickness;
	Real zn = std::pow( z, steepness );
	Real fa_proj_H = zn/(1 + zn);

	// Acc depth
	Real fa_depth_A = dot(Axyz-center, normal);
	internal_product = std::abs(fa_depth_A);
	z = internal_product;
	z /= thickness;
	zn = std::pow( z, steepness );
	Real fa_proj_A = zn/(1 + zn);

	Real fa_proj_AH = 0.5*(fa_proj_H+fa_proj_A);
	total_weight = fa_proj_AH * wat_weight + (1-fa_proj_AH) * memb_weight;

	return total_weight;
}

Real
get_membrane_depth_dependent_weight(
	Vector const & normal,
	Vector const & center,
	Real const & thickness,
	Real const & steepness,
	int const don_nb,
	int const acc_nb,
	Vector const & Hxyz, // proton
	Vector const & Axyz  // acceptor
)
{
	Real wat_weight(1.0), memb_weight(1.0), total_weight(1.0);

	// water phase smooth_hb_env_dep
	wat_weight = hb_env_dep_burial_lin( acc_nb, don_nb );

	// membrane phase dependent weight
	// Hdonor depth
	Real fa_depth_H = dot(Hxyz-center, normal);
	Real internal_product = std::abs(fa_depth_H);
	Real z = internal_product;
	z /= thickness;
	Real zn = std::pow( z, steepness );
	Real fa_proj_H = zn/(1 + zn);

	// Acc depth
	Real fa_depth_A = dot(Axyz-center, normal);
	internal_product = std::abs(fa_depth_A);
	z = internal_product;
	z /= thickness;
	zn = std::pow( z, steepness );
	Real fa_proj_A = zn/(1 + zn);

	Real fa_proj_AH = 0.5*(fa_proj_H+fa_proj_A);
	total_weight = fa_proj_AH * wat_weight + (1-fa_proj_AH) * memb_weight;

	return total_weight;
}

} // hbonds
} // scoring
} // core
