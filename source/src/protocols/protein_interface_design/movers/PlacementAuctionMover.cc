// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/PlacementAuctionMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)
/// @author Lei Shi (shilei@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/PlacementAuctionMover.hh>
#include <protocols/protein_interface_design/movers/PlacementAuctionMoverCreator.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <core/scoring/constraints/BackboneStubConstraint.hh>
#include <core/scoring/constraints/BackboneStubLinearConstraint.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>

//#include <protocols/docking/DockingProtocol.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/hotspot_hashing/HotspotStub.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
// Unit Headers
#include <protocols/filters/Filter.hh>
#include <protocols/protein_interface_design/movers/PlaceUtils.hh>

#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AA.hh>
#include <numeric/xyzVector.hh>


// Unit Headers
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

// C++ headers
#include <map>
#include <boost/foreach.hpp>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_moves/DesignRepackMover.hh>


using namespace protocols::protein_interface_design;

static thread_local basic::Tracer TR( "protocols.protein_interface_design.movers.PlacementAuctionMover" );

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace protocols::moves;
using namespace core;

std::string
PlacementAuctionMoverCreator::keyname() const
{
	return PlacementAuctionMoverCreator::mover_name();
}

protocols::moves::MoverOP
PlacementAuctionMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new PlacementAuctionMover );
}

std::string
PlacementAuctionMoverCreator::mover_name()
{
	return "Auction";
}

PlacementAuctionMover::PlacementAuctionMover() :
	simple_moves::DesignRepackMover( PlacementAuctionMoverCreator::mover_name() )
{}

PlacementAuctionMover::~PlacementAuctionMover() {}

protocols::moves::MoverOP
PlacementAuctionMover::clone() const {
	return( protocols::moves::MoverOP( new PlacementAuctionMover( *this ) ) );
}

PlacementAuctionMover::ResidueAuction
PlacementAuctionMover::auction_results() const {
	return auction_results_;
}

void
PlacementAuctionMover::insert( PlacementAuctionMover::ResidueAuctionItem const & item )
{
	auction_results_.insert( item );
}

core::Size
PlacementAuctionMover::size() const
{
	return auction_results_.size();
}

void
PlacementAuctionMover::erase( PlacementAuctionMover::iterator it ){
	auction_results_.erase( it );
}

void
PlacementAuctionMover::clear(){
	auction_results_.clear();
}

PlacementAuctionMover::iterator
PlacementAuctionMover::begin(){
	return auction_results_.begin();
}

PlacementAuctionMover::iterator
PlacementAuctionMover::end(){
	return auction_results_.end();
}

PlacementAuctionMover::const_iterator
PlacementAuctionMover::begin() const{
	return auction_results_.begin();
}

PlacementAuctionMover::const_iterator
PlacementAuctionMover::end() const{
	return auction_results_.end();
}

void
PlacementAuctionMover::apply( core::pose::Pose & pose )
{
  core::pose::Pose const saved_pose( pose ); // the pose should not actually be changed within this function

	using namespace protocols::hotspot_hashing;
	using namespace core::scoring;

	core::Size const host_chain_begin( pose.conformation().chain_begin( host_chain_ ) );
	core::Size const host_chain_end  ( pose.conformation().chain_end(   host_chain_ ) );
	using namespace core::scoring;
	ScoreFunctionOP only_stub_scorefxn( new ScoreFunction );
	only_stub_scorefxn->reset();

  protocols::filters::FilterOP stf;
	if ( stub_energy_fxn_ == "backbone_stub_constraint" ) {
		only_stub_scorefxn->set_weight( backbone_stub_constraint, 1.0 );
		stf = protocols::filters::FilterOP( new protocols::simple_filters::ScoreTypeFilter(only_stub_scorefxn, backbone_stub_constraint, 1.0 ) );

	} else if ( stub_energy_fxn_ == "backbone_stub_linear_constraint" ) {
		only_stub_scorefxn->set_weight( backbone_stub_linear_constraint, 1.0 );
		stf = protocols::filters::FilterOP( new protocols::simple_filters::ScoreTypeFilter(only_stub_scorefxn, backbone_stub_linear_constraint, 1.0 ) );
	} else {
		utility_exit_with_message( "ERROR: unrecognized stub_energy_fxn_. Only support backbone_stub_constraint or backbone_stub_linear_constraint");
	}

	core::Size fixed_res(1);
	if( host_chain_ == 1 ) fixed_res = pose.total_residue();
	core::id::AtomID const fixed_atom_id = core::id::AtomID( pose.residue(fixed_res).atom_index("CA"), fixed_res );
/// ResidueAuction is keyed by energy => we select the residue,stub,stubset combination with the best energy for each stubset,stub combination
	typedef std::pair< HotspotStubSetOP, HotspotStubOP > StubsetStubPair;
	//typedef std::pair< core::Real, std::pair< core::Size, StubsetStubPair > > ResidueAuctionItem;
	typedef std::multimap< core::Real, std::pair< core::Size, StubsetStubPair > > ResidueAuction;
/// Preventing positions that have already been prevented throught task factory or through the prevent_repacking method
	using namespace core::pack::task;
	PackerTaskOP task;
	if( task_factory() )
		task = task_factory()->create_task_and_apply_taskoperations( pose );
	else
		task = TaskFactory::create_packer_task( pose );
	utility::vector1< core::Size > host_positions;
	for( core::Size host_position( host_chain_begin+1 ); host_position<=host_chain_end-1; ++host_position ){
		using namespace core::chemical;
		// exclude gly/pro and don't allow prevented residues
		if( std::find( prevent_repacking_.begin(), prevent_repacking_.end(), host_position ) == prevent_repacking_.end() && pose.residue( host_position ).aa() != aa_gly && pose.residue( host_position ).aa() != aa_pro ){
			if( task->nonconst_residue_task( host_position ).being_packed() )
				host_positions.push_back( host_position );
		}
	}// for host_position

	ResidueAuction saved_auction; /// auction_results_ will be depleted in the following. Then, if successful, I'll reinstate it.
	if ( stub_energy_fxn_ == "backbone_stub_linear_constraint" ) {
		saved_auction = auction_results() ; /// auction_results_ will be depleted in the following. Then, if successful, I'll reinstate it.
	}

	BOOST_FOREACH( StubSetStubPos const hs_set, stub_sets_ ){
		HotspotStubSetOP stub_set( hs_set.first );
    //TR << "Loop restart: " << std::endl; //loop through each library
		BOOST_FOREACH( HotspotStubSet::Hs_data const stub_pair, *stub_set ){
			HotspotStubOP stub( stub_pair.second.second );
			BOOST_FOREACH( core::Size const host_residue, host_positions )
			{
				core::Real const distance( pose.residue( host_residue ).xyz( "CB" ).distance( stub->residue()->xyz( "CB" ) ) );
				if( distance >= max_cb_cb_dist_ ) continue;
				//TR<< "residue: " << pose.residue( host_residue ).seqpos() << " " << pose.residue( host_residue ).name()<< " stub: " << stub->residue()->seqpos() << " " << stub->residue()->name() << " distance: " << distance  << " max_cb_cb_dist_: " << max_cb_cb_dist_ <<std::endl;
				core::Real const bonus( stub->bonus_value() );
				// I'm circumventing add_hotspot_constraints to pose and adding the constraint directly
				// since there's no ambiguity here and no need to switch to ala pose etc. And I don't
				// want all the quality control machinery to be applied to this stub; I know it's good.
				using namespace core::scoring::constraints;
				ConstraintCOPs stub_constraints;
				core::conformation::Residue const host_res( pose.conformation().residue( host_residue ) );
  				if ( stub_energy_fxn_ == "backbone_stub_constraint" ) {
						stub_constraints.push_back( core::scoring::constraints::ConstraintOP( new BackboneStubConstraint( pose, host_residue, fixed_atom_id, host_res, bonus, cb_force_ ) ) );
						stub_constraints = pose.add_constraints( stub_constraints );
						core::Real const bb_cst_score( stf->report_sm( pose ) );
						if( bb_cst_score <= -0.5 ) // take only residues that make some appreciable contribution
							insert( std::make_pair( bb_cst_score, std::make_pair( host_residue, std::make_pair( stub_set, stub ) ) ) );
				} else if ( stub_energy_fxn_ == "backbone_stub_linear_constraint" ) {
						stub_constraints.push_back( core::scoring::constraints::ConstraintOP( new BackboneStubLinearConstraint( pose, host_residue, fixed_atom_id, *(stub->residue()), bonus, cb_force_ ) ) );
						stub_constraints = pose.add_constraints( stub_constraints );
        		core::Real const bb_cst_score( stf->report_sm( pose ) );
						insert( std::make_pair( bonus+bb_cst_score, std::make_pair( host_residue, std::make_pair( stub_set, stub ) ) ) );
				}
			 pose = saved_pose;
			}// foreach host_residue
		}//foreach stub_pair

		//only for backbone_stub_linear_constraint
     if ( stub_energy_fxn_ == "backbone_stub_linear_constraint" ) {
				if( size() == 0 ){
					TR<<"No pairing found. Failing."<<std::endl;
					set_last_move_status( protocols::moves::FAIL_RETRY );
					return;
				}

      //Insert all auctioned postiions to be used by PlaceSimultaneousMover
       core::Real best_combined_energy=100;
       //TR<<"Total possiblitly found: " << size() <<std::endl;
       for( ResidueAuction::iterator lowest_energy = begin(); lowest_energy != end(); ++lowest_energy) {
              if (lowest_energy->first <= best_combined_energy ) {
                 saved_auction.insert(*lowest_energy);
                 core::Size const position( lowest_energy->second.first );
                 HotspotStubSetCOP stubset( lowest_energy->second.second.first );
                 HotspotStubOP stub( lowest_energy->second.second.second );

                 //assign matched positions to stub_sets_
                BOOST_FOREACH( StubSetStubPos & hs_set, stub_sets() ){
                 if ( position == hs_set.second.second ) // it has already being used!
                   break;
                 if( hs_set.first == stubset ){
                   hs_set.second.first = stub;
                   hs_set.second.second = position;
                   break;
                 }
               }
             }
           }
           clear(); //clear auction results for one hotspot
		} //end of backbone_stub_linear_constraint
	}//foreach hs_set

  if ( stub_energy_fxn_ == "backbone_stub_constraint" ) {
					if( size() == 0 ){
						TR<<"No pairing found. Failing."<<std::endl;
						set_last_move_status( protocols::moves::FAIL_RETRY );
						return;
					}

					saved_auction = auction_results() ; /// auction_results_ will be depleted in the following. Then, if successful, I'll reinstate it.
					while( size() > 0 ){
						// The auction ensures that each position is paired to the stubset that bids the lowest-energy stub on that
						/// position, but allows each stubset to bid multiple positions. This allows stubsets that lose on one position
						// to potentially succeed on another.
						PlacementAuctionMover::const_iterator lowest_energy( begin() );
						core::Size const position( lowest_energy->second.first );
						HotspotStubSetCOP stubset( lowest_energy->second.second.first );
						HotspotStubOP stub( lowest_energy->second.second.second );
						BOOST_FOREACH( StubSetStubPos & hs_set, stub_sets() ){
					// This is where the pairing takes place
							if( hs_set.first == stubset ){
								hs_set.second.first = stub;
								hs_set.second.second = position;
								break;
							}
						}

						for( ResidueAuction::iterator energy_set_pair = begin(); energy_set_pair != end(); /*incrementing done within the loop*/ ){
							ResidueAuction::iterator next_it = energy_set_pair;
							core::Size const erased_pos( energy_set_pair->second.first );
							HotspotStubSetCOP erased_stubset( energy_set_pair->second.second.first );

							if( position == erased_pos || stubset == erased_stubset ){
								++next_it;
								erase( energy_set_pair );
								energy_set_pair = next_it;
							}
							else ++energy_set_pair;
						}//for energy_set_pair
					}//while size()
	} //end of backbone_stub_constraint

	//check if all stub positions have been paired
	BOOST_FOREACH( StubSetStubPos const stubset_pos_pair, stub_sets() ){
		core::Size const pos( stubset_pos_pair.second.second );
		if( pos == 0 ){
			TR<<"Pairing failed"<<std::endl;
			set_last_move_status( protocols::moves::FAIL_RETRY );
			return;
		}
	}//foreach stubset_pos_pair

	// If all went well then retrieve the saved copy of the auction
	auction_results_ = saved_auction;
	TR<<"Pairing successful: "<< size() << std::endl;
	set_last_move_status( protocols::moves::MS_SUCCESS );
	pose=saved_pose;
	TR.flush();
}

std::string
PlacementAuctionMover::get_name() const {
	return PlacementAuctionMoverCreator::mover_name();
}

void
PlacementAuctionMover::host_chain( core::Size const hc ){
	host_chain_ = hc;
}

void
PlacementAuctionMover::max_cb_cb_dist( core::Real const c ){
	max_cb_cb_dist_ = c;
}

void
PlacementAuctionMover::stub_sets( utility::vector1< StubSetStubPos > const & stub_sets ){
	stub_sets_ = stub_sets;
}

void
PlacementAuctionMover::cb_force( core::Real const cb ){
	cb_force_ = cb;
}

std::string
PlacementAuctionMover::get_stub_scorefxn() const {
	return stub_energy_fxn_;
}

utility::vector1< protocols::protein_interface_design::movers::PlacementAuctionMover::StubSetStubPos > const &
PlacementAuctionMover::stub_sets() const{
	return stub_sets_;
}

utility::vector1< protocols::protein_interface_design::movers::PlacementAuctionMover::StubSetStubPos > &
PlacementAuctionMover::stub_sets() {
	return stub_sets_;
}

void
PlacementAuctionMover::parse_my_tag( TagCOP const tag,
		basic::datacache::DataMap &data,
		protocols::filters::Filters_map const &,
		Movers_map const &,
		core::pose::Pose const & pose )
{
	using namespace protocols::hotspot_hashing;
	using namespace protocols::filters;
	using namespace core::scoring;

	host_chain_ = tag->getOption<core::Size>( "chain_to_design", 2 );
	max_cb_cb_dist_ = tag->getOption< core::Real >( "max_cb_dist", 3.0 );
	cb_force_ = tag->getOption< core::Real >( "cb_force", 0.5 );
	stub_energy_fxn_ = tag->getOption<std::string>( "stubscorefxn", "backbone_stub_constraint" ) ;
	stub_sets_ = parse_stub_sets( tag, pose, host_chain_, data );
	runtime_assert( stub_sets_.size() );
	TR<<"max cb cb distance set to "<<max_cb_cb_dist_<<" and cb_force to "<<cb_force_<< " stub energy function" << stub_energy_fxn_ << '\n';
	TR<<"PlacementAuction mover on chain "<<host_chain_<<" with repack_non_ala set to "<<std::endl;
}

} //movers
} //protein_interface_design
} //protocols

