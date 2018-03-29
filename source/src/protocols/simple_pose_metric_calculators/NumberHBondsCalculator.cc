// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/metrics/NumberHBondsCalculator.cc
/// @brief  number of hbonds calculator class
/// @author Florian Richter
/// @author Scott Boyken ( sboyken@gmail.com); updated 2017 to be compatible with new unsat filters and generous way of h-bond counting

// Unit headers
#include <protocols/simple_pose_metric_calculators/NumberHBondsCalculator.hh>
#include <core/pose/Pose.hh>


// Utility headers
#include <utility>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>

// core headers
//#include <core/pose/util.tmpl.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

#include <utility/assert.hh>

#include <utility/vector1.hh>

using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;

#ifdef    SERIALIZATION
// Project serialization headers
#include <core/id/AtomID_Map.srlz.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/set.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace simple_pose_metric_calculators {


std::string
choose_hbond_parameter_set() {
	// use the "score12_params" set if the -restore_pre_talaris_2013_behavior flag is on the command line
	// and otherwise use the new and improved sp2_elec_params parameter set
	if ( basic::options::option[ basic::options::OptionKeys::mistakes::restore_pre_talaris_2013_behavior ] ) {
		return "score12_params";
	} else {
		return "sp2_elec_params";
	}
}

// default constructor
NumberHBondsCalculator::NumberHBondsCalculator() :
	use_generous_hbonds_( true ),
	special_region_( /* NULL */ ),
	all_Hbonds_( /* NULL */ ),
	special_region_Hbonds_( /* NULL */ ),
	atom_Hbonds_( /* NULL */ ),
	residue_Hbonds_( /* NULL */ ),
	ref_residue_total_energies_( /* NULL */ )
	//hbond_set_( /* NULL */ ),
	//hb_database_( core::scoring::hbonds::HBondDatabase::get_database( choose_hbond_parameter_set() ) )
{}

// constructor with options
NumberHBondsCalculator::NumberHBondsCalculator( bool const generous, std::set< core::Size > special_region /* NULL */ ) :
	use_generous_hbonds_( generous ),
	special_region_(special_region),
	all_Hbonds_( /* NULL */ ),
	special_region_Hbonds_( /* NULL */ ),
	atom_Hbonds_( /* NULL */ ),
	residue_Hbonds_( /* NULL */ ),
	ref_residue_total_energies_( /* NULL */ )
	//hbond_set_( /* NULL */ ),
	//hb_database_( core::scoring::hbonds::HBondDatabase::get_database( choose_hbond_parameter_set() ) )
{}

// desctructor
NumberHBondsCalculator::~NumberHBondsCalculator(){}

//core::pose::metrics::PoseMetricCalculatorOP NumberHBondsCalculator::clone() const { return core::pose::metrics::PoseMetricCalculatorOP( new NumberHBondsCalculator() ); }

void
NumberHBondsCalculator::lookup(
	std::string const & key,
	basic::MetricValueBase * valptr
) const
{

	if ( key == "all_Hbonds" ) {
		basic::check_cast( valptr, &all_Hbonds_, "all_Hbonds expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( all_Hbonds_ );

	} else if ( key == "special_region_Hbonds" ) {
		basic::check_cast( valptr, &special_region_Hbonds_, "special_region_Hbonds expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( special_region_Hbonds_ );

	} else if ( key == "atom_Hbonds" ) {
		basic::check_cast( valptr, &atom_Hbonds_, "atom_Hbonds expects to return a id::AtomID_Map< Size >" );
		(static_cast<basic::MetricValue<id::AtomID_Map< Size > > *>(valptr))->set( atom_Hbonds_ );

	} else if ( key == "residue_Hbonds" ) {
		basic::check_cast( valptr, &residue_Hbonds_, "residue_Hbonds expects to return a utility::vector1< Size >" );
		(static_cast<basic::MetricValue<utility::vector1< Size > > *>(valptr))->set( residue_Hbonds_ );

	} else {
		basic::Error() << "NumberHbondsCalculator cannot compute the requested metric " << key << std::endl;
		utility_exit();
	}

} //lookup


std::string
NumberHBondsCalculator::print( std::string const & key ) const
{

	if ( key == "all_Hbonds" ) {
		return utility::to_string( all_Hbonds_ );
	} else if ( key == "special_region_Hbonds" ) {
		return utility::to_string( special_region_Hbonds_ );
	} else if ( key == "atom_Hbonds" ) {
		basic::Error() << "id::AtomID_Map< Size > has no output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "residue_Hbonds" ) {
		return utility::to_string( residue_Hbonds_ );
	}

	basic::Error() << "NumberHbondsCalculator cannot compute metric " << key << std::endl;
	utility_exit();
	return "";

} //print


void
NumberHBondsCalculator::recompute( Pose const & this_pose )
{

	using namespace core::scoring;
	using namespace core::scoring::hbonds;

	//first we have to figure out which of the hbonds to (re)calculate
	utility::vector1< bool >res_to_recompute( this_pose.size(), false );
	//hbonds::HBondSet hb_set( this_pose.size() );

	//hbonds::HBondSet test_hb_set( this_pose.size() );
	//hbonds::fill_hbond_set( this_pose, false, test_hb_set, false);

	HBondSet temp_hbond_set; // we need to recompute HBondSet each time
	HBondOptions new_options( temp_hbond_set.hbond_options() ); // set default options
	if ( use_generous_hbonds_ ) {
		new_options.use_hb_env_dep(false);
		new_options.decompose_bb_hb_into_pair_energies(true);
		new_options.bb_donor_acceptor_check( false ); // don't use bb exclusion logic when penalizing unsatisfied -- ideally would only eclude N-H donors and not exclude C=O with only 1 h-bond
		// according to ALF will get rid of exclusion logic soon
		new_options.exclude_intra_res_protein( false ); // we want to count this for unsat calc, by default they are excluded
		new_options.exclude_DNA_DNA( false );
		new_options.exclude_intra_res_RNA( false );
		new_options.exclude_self_hbonds(false);
	}
	HBondSetOP hb_set( new HBondSet( new_options, this_pose.size() ) ); // we need to recompute HBondSet each time

	determine_res_to_recompute( this_pose, res_to_recompute );
	runtime_assert( residue_Hbonds_.size() == res_to_recompute.size() );

	//set the hbonds vector for all recompute residues to 0
	for ( Size i = 1; i <= res_to_recompute.size(); ++i ) {

		if ( res_to_recompute[ i ] ) residue_Hbonds_[ i ] = 0;
	}


	for ( Size i = 1; i <= res_to_recompute.size(); ++i ) {

		compute_Hbonds_for_residue( this_pose, i, res_to_recompute, *hb_set );

	} //loop over all residues


	//some double checking
	//std::cerr << "comp hbset has " << test_hb_set.nhbonds() << " hbonds, while other has " << hb_set.nhbonds() << std::endl;


	//now we have to setup the AtomID Map
	for ( Size i = 1; i <= res_to_recompute.size(); ++i ) {

		if ( !res_to_recompute[i] ) continue;

		conformation::Residue const & rsd = this_pose.residue( i );

		for ( Size at = 1; at <= rsd.nheavyatoms(); ++at ) {

			core::id::AtomID atid( at, i );

			if ( residue_Hbonds_[i] == 0 ) {
				atom_Hbonds_.set( atid, 0 );
				continue;
			}

			if ( rsd.atom_type( at ).is_acceptor() && rsd.atom_type( at ).is_donor() ) {
				Size hbonds_this_donor(0);
				// count donated hbonds
				for ( Size hcount = rsd.type().attached_H_begin( at ); hcount<= rsd.type().attached_H_end( at ); hcount++ ) {
					hbonds_this_donor = hbonds_this_donor + hb_set->atom_hbonds( core::id::AtomID (hcount, i ), !use_generous_hbonds_ /* include_only_allowed */ ).size();
					// sboyken added; want to store Hpol as +1 in atom map if makes an h-bond
					//    the !use_generous_hbonds_ results in giving us all h-bonds, which we want for generous case
					atom_Hbonds_.set( core::id::AtomID(hcount, i ), hb_set->atom_hbonds( core::id::AtomID (hcount, i ), !use_generous_hbonds_ /* include_only_allowed */ ).size() );
				}
				// for heavy atom, store total
				atom_Hbonds_.set( atid, hbonds_this_donor + hb_set->atom_hbonds( atid, !use_generous_hbonds_ /* include_only_allowed */ ).size() );

			} else if ( rsd.atom_type( at ).is_acceptor() ) {
				// for each heavy atom acceptor, store number of h-bonds it participates in
				atom_Hbonds_.set( atid, hb_set->atom_hbonds( atid, !use_generous_hbonds_ /* include_only_allowed */ ).size() );
			} else if ( rsd.atom_type( at ).is_donor() ) {
				Size hbonds_this_donor(0);

				for ( Size hcount = rsd.type().attached_H_begin( at ); hcount<= rsd.type().attached_H_end( at ); hcount++ ) {

					hbonds_this_donor = hbonds_this_donor + hb_set->atom_hbonds( core::id::AtomID (hcount, i ), !use_generous_hbonds_ /* include_only_allowed */ ).size();
					// sboyken added; want to store Hpol as +1 in atom map if makes an h-bond
					//    the !use_generous_hbonds_ results in giving us all h-bonds, which we want for generous case
					atom_Hbonds_.set( core::id::AtomID(hcount, i ), hb_set->atom_hbonds( core::id::AtomID (hcount, i ), !use_generous_hbonds_ /* include_only_allowed */ ).size() );
				}
				// for each heavy atom donor, store number of h-bonds it participates in
				atom_Hbonds_.set( atid, hbonds_this_donor );
			} else atom_Hbonds_.set( atid, 0 );
		}
	}


	all_Hbonds_ = hb_set->nhbonds();

	//hbond_set_ = hb_set; // store it for later access, but we need to recompute HBondSet each time

	//  this is not accurate, this assumes only inter h-bonds
	// //finally, we compute the total number of h-bonds
	// all_Hbonds_ = 0;
	// for ( utility::vector1< core::Size >::const_iterator resh_it = residue_Hbonds_.begin();
	//   resh_it != residue_Hbonds_.end(); ++resh_it ) {
	//  all_Hbonds_ = all_Hbonds_ + *resh_it;
	// }
	// all_Hbonds_ = all_Hbonds_ / 2; //remember not to overcount

} //recompute


/// @brief function to figure out which residues to recompute the hydrogen bonds for
/// @brief strategy: for each residue, we check whether the internally cached total energies
/// @brief           correspond to the energies found in the pose for that residue. If they do,
/// @brief           this means that the number of H-bonds hasn't changed.
void
NumberHBondsCalculator::determine_res_to_recompute(
	core::pose::Pose const & pose,
	utility::vector1< bool > & res_to_recompute)
{
	using namespace core::scoring;

	//check1: does the internal reference array have the same size as the pose has residues?
	//if not, means we have to recompute everything
	if ( ref_residue_total_energies_.size() != pose.size() ) {
		atom_Hbonds_.clear();
		atom_Hbonds_.resize( pose.size() );
		residue_Hbonds_.resize( pose.size() );
		ref_residue_total_energies_.resize( pose.size() );
		for ( Size i = 1; i <= pose.size(); ++i ) {
			res_to_recompute[ i ] = true;
			ref_residue_total_energies_[i] = pose.energies().residue_total_energies(i)[ total_score ];
		}

		return;
	}


	//check2: for each residue, check whether calculator cached residue energies have changed
	for ( Size i = 1; i <= pose.size(); ++i ) {

		if ( ref_residue_total_energies_[i] != pose.energies().residue_total_energies(i)[ total_score ] ) {

			ref_residue_total_energies_[i] = pose.energies().residue_total_energies(i)[ total_score ];
			res_to_recompute[ i ] = true;
			//std::cerr << "res " << i << " needs recomputing." << std::endl;
		}

		if ( special_region_.find( i ) != special_region_.end() ) res_to_recompute[ i ] = true;
	}

} //determine_res_to_recompute


void
NumberHBondsCalculator::compute_Hbonds_for_residue(
	core::pose::Pose const & pose,
	core::Size i,
	utility::vector1< bool > const & res_to_recompute,
	core::scoring::hbonds::HBondSet & hb_set)
{
	using namespace core::scoring;
	using namespace core::scoring::hbonds;


	conformation::Residue const & rsd1( pose.residue( i ) );
	int const nb1 = pose.energies().tenA_neighbor_graph().get_node( i )->num_neighbors_counting_self();

	// sboyken added
	// no need to store database, can get as needed
	core::scoring::hbonds::HBondDatabaseCOP hb_database( core::scoring::hbonds::HBondDatabase::get_database( choose_hbond_parameter_set() ) );

	if ( use_generous_hbonds_ ) {
		identify_intra_res_hbonds( *hb_database, rsd1, nb1, false /*evaluate derivative*/, hb_set);
	}
	// will not add these to residue_Hbonds_[ i ] for now

	//go over the neighbours of this residue
	//in the pose energy graph and see whether there
	//is an Hbond energy between them
	for ( utility::graph::EdgeListConstIterator egraph_it = pose.energies().energy_graph().get_node( i )->const_upper_edge_list_begin();
			egraph_it != pose.energies().energy_graph().get_node( i )->const_upper_edge_list_end(); ++egraph_it ) {

		Size other_node = (*egraph_it)->get_other_ind( i );

		if ( (!res_to_recompute[i]) && !(res_to_recompute[other_node] ) ) continue;

		//EnergyEdge const * Eedge = static_cast< EnergyEdge const * > (*egraph_it);

		//if( sum_Hbond_terms( Eedge->energy_map() ) <= 100 ){ //hardcoded cutoff for now

		//ok, we gotta compute the number of Hbonds between these two residues
		Size prev_no_hb = hb_set.nhbonds();

		conformation::Residue const & rsd2( pose.residue( other_node ) );

		int const nb2 = pose.energies().tenA_neighbor_graph().get_node( other_node )->num_neighbors_counting_self();

		// rsd1 as donor, rsd2 as acceptor
		hbonds::identify_hbonds_1way( *hb_database, rsd1, rsd2, nb1, nb2, false /*evaluate derivative*/,
			false, false, false, false, hb_set);
		//rsd2 as  donor, rsd1 as acceptor
		hbonds::identify_hbonds_1way( *hb_database, rsd2, rsd1, nb2, nb1, false /*evaluate derivative*/,
			false, false, false, false, hb_set);

		Size num_hb_these_two_res = hb_set.nhbonds() - prev_no_hb;

		if ( res_to_recompute[i] ) residue_Hbonds_[ i ] += num_hb_these_two_res;
		if ( res_to_recompute[ other_node ] ) residue_Hbonds_[ other_node ] += num_hb_these_two_res;

		if ( special_region_.find( i ) != special_region_.end() ) special_region_Hbonds_ += num_hb_these_two_res;
		else if ( special_region_.find( other_node ) != special_region_.end() ) special_region_Hbonds_ += num_hb_these_two_res;

		//std::cerr << "residues " << i << " and " << other_node << " make " << num_hb_this_two_res << " hbonds.\n";

	}//iterator over all neighbors for this residue


} //compute_Hbonds_for_residue


core::Real
NumberHBondsCalculator::sum_Hbond_terms(
	core::scoring::EnergyMap const & emap
)
{
	using namespace core::scoring;

	return emap[ hbond_sr_bb ] + emap[hbond_lr_bb] + emap[ hbond_bb_sc ] + emap[hbond_sc];

}

} //namespace simple_pose_metric_calculators
} //namespace protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::simple_pose_metric_calculators::NumberHBondsCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::EnergyDependentCalculator >( this ) );
	// arc( CEREAL_NVP( hb_database_ ) ); // core::scoring::hbonds::HBondDatabaseCOP

	// EXEMPT hb_database
	// Don't serialize the hbond database -- instead, plan to get the global HBondDatabase
	// during deserialization

	arc( CEREAL_NVP( use_generous_hbonds_ ) ); // bool
	arc( CEREAL_NVP( all_Hbonds_ ) ); // core::Size
	arc( CEREAL_NVP( special_region_Hbonds_ ) ); // core::Size
	arc( CEREAL_NVP( atom_Hbonds_ ) ); // core::id::AtomID_Map<core::Size>
	arc( CEREAL_NVP( residue_Hbonds_ ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( ref_residue_total_energies_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( special_region_ ) ); // std::set<core::Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::simple_pose_metric_calculators::NumberHBondsCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::EnergyDependentCalculator >( this ) );

	//hb_database_ = core::scoring::hbonds::HBondDatabase::get_database( choose_hbond_parameter_set() );

	arc( use_generous_hbonds_ ); // bool
	arc( all_Hbonds_ ); // core::Size
	arc( special_region_Hbonds_ ); // core::Size
	arc( atom_Hbonds_ ); // core::id::AtomID_Map<core::Size>
	arc( residue_Hbonds_ ); // utility::vector1<core::Size>
	arc( ref_residue_total_energies_ ); // utility::vector1<core::Real>
	arc( special_region_ ); // std::set<core::Size>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::simple_pose_metric_calculators::NumberHBondsCalculator );
CEREAL_REGISTER_TYPE( protocols::simple_pose_metric_calculators::NumberHBondsCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_simple_pose_metric_calculators_NumberHBondsCalculator )
#endif // SERIALIZATION
