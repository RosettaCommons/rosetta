// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/PoseMetricCalculator/NonlocalContactsCalculator.cc
/// @brief  calculator to compute nonlocal/tertiary contacts in a given pose
/// @author Florian Richter

// Unit headers
#include <protocols/toolbox/pose_metric_calculators/NonlocalContactsCalculator.hh>

//#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/chemical/ResidueType.hh>
#include <core/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/MetricValue.hh>


#include <utility/assert.hh>

#include <utility/vector1.hh>



using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;

static thread_local basic::Tracer TR( "protocols/toolbox/PoseMetricCalculators/NonlocalContactsCalculator" );

namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {


NonlocalContactsCalculator::NonlocalContactsCalculator(
	core::Size min_sequence_separation,
	core::Real contact_cutoffE
) : total_nlcontacts_(0),
		special_region1_nlcontacts_(0),
		special_region2_nlcontacts_(0),
		special_region1_intra_nlcontacts_(0),
		special_region1_to_other_nlcontacts_(0),
		region1_region2_nlcontacts_(0),
		nlcontacts_graph_( /* NULL */ ),
		min_seq_separation_(min_sequence_separation),
		cutoffE_(contact_cutoffE)
{
	residue_nlcontacts_.clear();
	residue_nlscore_.clear();
  special_region1_.clear();
	special_region2_.clear();
}


NonlocalContactsCalculator::NonlocalContactsCalculator(
	std::set< core::Size > const & special_region,
	core::Size min_sequence_separation,
	core::Real contact_cutoffE
) : total_nlcontacts_(0),
		special_region1_nlcontacts_(0),
		special_region2_nlcontacts_(0),
		special_region1_intra_nlcontacts_(0),
		special_region1_to_other_nlcontacts_(0),
		region1_region2_nlcontacts_(0),
		nlcontacts_graph_( /* NULL */ ),
		min_seq_separation_(min_sequence_separation),
		cutoffE_(contact_cutoffE),
		special_region1_(special_region)
{
	residue_nlcontacts_.clear();
	residue_nlscore_.clear();
	special_region2_.clear();
}


NonlocalContactsCalculator::NonlocalContactsCalculator(
	std::set< core::Size > const & special_region1,
	std::set< core::Size > const & special_region2,
	core::Size min_sequence_separation,
	core::Real contact_cutoffE
) : total_nlcontacts_(0),
		special_region1_nlcontacts_(0),
		special_region2_nlcontacts_(0),
		special_region1_intra_nlcontacts_(0),
		special_region1_to_other_nlcontacts_(0),
		region1_region2_nlcontacts_(0),
		nlcontacts_graph_( /* NULL */ ),
		min_seq_separation_(min_sequence_separation),
		cutoffE_(contact_cutoffE),
		special_region1_(special_region1),
		special_region2_(special_region2)
{
	residue_nlcontacts_.clear();
	residue_nlscore_.clear();
}

NonlocalContactsCalculator::~NonlocalContactsCalculator(){}


void
NonlocalContactsCalculator::lookup(
  std::string const & key,
  basic::MetricValueBase * valptr
) const
{

	if ( key == "total_nlcontacts" ) {
		basic::check_cast( valptr, &total_nlcontacts_, "total_nlcontacts expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( total_nlcontacts_ );

	} else if ( (key == "special_region_nlcontacts") || (key == "special_region1_nlcontacts") ) {
		basic::check_cast( valptr, &special_region1_nlcontacts_, "special_region_nlcontacts expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( special_region1_nlcontacts_ );

	} else if ( key == "special_region2_nlcontacts_" ) {
		basic::check_cast( valptr, &special_region2_nlcontacts_, "special_region2_nlcontacts expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( special_region2_nlcontacts_ );

	} else if ( key == "special_region1_intra_nlcontacts_" ) {
		basic::check_cast( valptr, &special_region1_intra_nlcontacts_, "special_region1_intra_nlcontacts expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( special_region1_intra_nlcontacts_ );

	} else if ( key == "special_region1_to_other_nlcontacts_" ) {
		basic::check_cast( valptr, &special_region1_to_other_nlcontacts_, "special_region1_to_other_nlcontacts expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( special_region1_to_other_nlcontacts_ );

	} else if ( key == "region1_region2_nlcontacts_" ) {
		basic::check_cast( valptr, &region1_region2_nlcontacts_, "regio1_region2_nlcontacts expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( region1_region2_nlcontacts_ );

	} else if ( key == "residue_nlcontacts" ) {
		basic::check_cast( valptr, &residue_nlcontacts_, "residue_nlcontacts expects to return a utility::vector1< Size >" );
		(static_cast<basic::MetricValue<utility::vector1< Size > > *>(valptr))->set( residue_nlcontacts_ );

	} else if ( key == "residue_nlscore" ) {
		basic::check_cast( valptr, &residue_nlscore_, "residue_nlscore expects to return a utility::vector1< Real >" );
		(static_cast<basic::MetricValue<utility::vector1< Real > > *>(valptr))->set( residue_nlscore_ );

	} else if ( key == "nlcontacts_graph" ) {
		basic::check_cast( valptr, &nlcontacts_graph_, "nlcontacts_graph expects to return a core::Graph::GraphOP" );
		(static_cast<basic::MetricValue< core::graph::GraphOP > *>(valptr))->set( nlcontacts_graph_ );

	}

	else {
		basic::Error() << "NonlocalContactsCalculator cannot compute the requested metric " << key << std::endl;
		utility_exit();
	}

} //lookup



std::string
NonlocalContactsCalculator::print( std::string const & key ) const
{


  basic::Error() << "NonlocalContactsCalculator cannot compute metric " << key << std::endl;
  utility_exit();
  return "";

} //print


void
NonlocalContactsCalculator::recompute( Pose const & this_pose )
{
	using namespace core::scoring;

	residue_nlcontacts_.clear();
	residue_nlscore_.clear();
	residue_nlcontacts_.resize( this_pose.total_residue(), 0 );
	residue_nlscore_.resize( this_pose.total_residue(), 0.0 );
	total_nlcontacts_ = 0;
	special_region1_nlcontacts_ = 0;
	special_region2_nlcontacts_ = 0;
	special_region1_intra_nlcontacts_ = 0;
	region1_region2_nlcontacts_ = 0;

	nlcontacts_graph_ = core::graph::GraphOP( new core::graph::Graph( this_pose.total_residue() ) );

	EnergyMap cur_weights = this_pose.energies().weights();

	for( core::Size i = 1; i <= this_pose.total_residue(); ++i ){

		if( ! this_pose.residue_type( i ).is_protein() ) continue;
		//get the node for this residue in the energy graph
		  for( graph::EdgeListConstIterator egraph_it = this_pose.energies().energy_graph().get_node( i )->const_upper_edge_list_begin();
       egraph_it != this_pose.energies().energy_graph().get_node( i )->const_upper_edge_list_end(); ++egraph_it){

				core::Size other_res = (*egraph_it)->get_other_ind( i );

				if( ( ( other_res - i ) <= min_seq_separation_ ) || !this_pose.residue_type( other_res ).is_protein() ) continue;

				//TR << other_res << " - " << i << " is bigger than " << min_seq_separation_ << std::endl;
				//downcast to energy edge
				EnergyEdge const * Eedge = static_cast< EnergyEdge const * > (*egraph_it);

				//to do: get the long range energies

				core::Real resresE( Eedge->dot( cur_weights ) );
				core::Real resresE_half( resresE / 2);

				residue_nlscore_[ i ] += resresE_half;
				residue_nlscore_[ other_res ] += resresE_half;

				if( resresE <= cutoffE_ ){

					TR.Debug << "residues " << i << " and " << other_res << " make nonlocal contact, interactionE is " << resresE << std::endl;

					total_nlcontacts_++;
					residue_nlcontacts_[i]++;
					residue_nlcontacts_[other_res]++;
					nlcontacts_graph_->add_edge( i, other_res );

					bool i_in_region1( special_region1_.find(i) != special_region1_.end() );
					bool i_in_region2( special_region2_.find(i) != special_region2_.end() );

					bool other_in_region1( special_region1_.find( other_res ) != special_region1_.end() );
					bool other_in_region2( special_region2_.find( other_res ) != special_region2_.end() );

					if( i_in_region1 || other_in_region1 ) special_region1_nlcontacts_++;

					if( i_in_region2 || other_in_region2 ) special_region2_nlcontacts_++;

					if( ( i_in_region1 && other_in_region2 ) || ( i_in_region2 && other_in_region1 ) ) region1_region2_nlcontacts_++;
					else if( i_in_region1 && other_in_region1 ) special_region1_intra_nlcontacts_++;


				} //if residues form nonlocal contact

			} //neighbors of this residue
	} //loop over residues

	special_region1_to_other_nlcontacts_ = special_region1_nlcontacts_ - special_region1_intra_nlcontacts_;

	//some optional debug output. done at the end to prevent tracer if evaluations in normal production runs
	if( basic::options::option[basic::options::OptionKeys::out::level] >= basic::t_debug ){
		for( core::Size i = 1; i <= this_pose.total_residue(); ++i) TR.Debug << "Residue " << i << " makes " << residue_nlcontacts_[i] << " nonlocal contacts and has a total nonlocal interaction energy of " << residue_nlscore_[i] << "." << std::endl;
	}
} //recompute


} //namespace pose_metric_calculators
} //namespace toolbox
} //namespace protocols
