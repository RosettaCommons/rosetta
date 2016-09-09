// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ResidueSetChainEnergyFilter.cc
/// @brief
/// @author Chris King (chrisk1@uw.edu)


//Unit Headers
#include <protocols/simple_filters/ResidueSetChainEnergyFilter.hh>
#include <protocols/simple_filters/ResidueSetChainEnergyFilterCreator.hh>

//Project Headers
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/exit.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/format.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/excn/Exceptions.hh>
#include <core/pose/selection.hh>

namespace protocols {
namespace simple_filters {

using namespace core;
using namespace core::scoring;
using namespace ObjexxFCL::format;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.ResidueSetChainEnergyFilter" );

protocols::filters::FilterOP
ResidueSetChainEnergyFilterCreator::create_filter() const { return protocols::filters::FilterOP( new ResidueSetChainEnergyFilter ); }

std::string
ResidueSetChainEnergyFilterCreator::keyname() const { return "ResidueSetChainEnergy"; }

//default ctor
ResidueSetChainEnergyFilter::ResidueSetChainEnergyFilter() :
	protocols::filters::Filter( "ResidueSetChainEnergy" )
{}

//full ctor, default ctor defined in header file
ResidueSetChainEnergyFilter::ResidueSetChainEnergyFilter(
	core::scoring::ScoreFunctionCOP scorefxn,
	core::scoring::ScoreType const score_type,
	core::Real const score_type_threshold,
	std::string const string_resnums,
	core::Size const chain
) :
	protocols::filters::Filter( "ResidueSetChainEnergy" )
{
	score_type_ = score_type;
	score_type_threshold_ = score_type_threshold;
	scorefxn_ = scorefxn->clone();
	string_resnums_ = string_resnums;
	chain_ = chain;
}

ResidueSetChainEnergyFilter::~ResidueSetChainEnergyFilter() = default;

void
ResidueSetChainEnergyFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	using namespace core::scoring;

	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();

	score_type_ = core::scoring::score_type_from_name( tag->getOption<std::string>( "score_type", "total_score" ) );
	if ( ! tag->hasOption( "threshold" ) ) throw utility::excn::EXCN_RosettaScriptsOption("Must specify 'threshold' for ResidueSetChainEnergyFilter.");
	score_type_threshold_ = tag->getOption< core::Real >( "threshold" );

	TR<< "filter for score_type "<< score_type_ <<" with threshold "<< score_type_threshold_ << std::endl;

	string_resnums_ = tag->getOption< std::string >( "resnums" );// these are kept in memory until the pose is available (at apply time)
	chain_ = tag->getOption< core::Size >( "chain", 0 );
}

bool
ResidueSetChainEnergyFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const score( compute( pose ) );
	TR << "score " << core::scoring::ScoreTypeManager::name_from_score_type( score_type_ ) << " is " << score << ". ";
	if ( score <= score_type_threshold_ ) {
		TR<<"passing." << std::endl;
		return true;
	} else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

void
ResidueSetChainEnergyFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	out<<"Weighted score of "<<core::scoring::ScoreTypeManager::name_from_score_type( score_type_ )<<" "<<compute( pose )<<'\n';
}

core::Real
ResidueSetChainEnergyFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
ResidueSetChainEnergyFilter::compute(
	core::pose::Pose const & pose
) const {

	using namespace core::pose;
	using namespace core::scoring;

	//copy the pose so we can set a possibly new scorefxn in it (input to apply is const!)
	Pose in_pose( pose );

	// make sure that scoring weights are compatible with pose's residue type set
	// check centroid case
	//TODO: this is a hacky and not robust way of checking for centroid!
	if ( ( ( *scorefxn_ )[ fa_rep ] == 0.0 && ( *scorefxn_ )[ fa_atr ] == 0.0 ) // full atom terms are off
			&& ( ( *scorefxn_ )[ interchain_vdw ] > 0.0 || ( *scorefxn_ )[ vdw ] > 0.0 )  ) { // a centroid term is on
		if ( in_pose.is_fullatom() ) { // but pose is full atom
			core::util::switch_to_residue_type_set( in_pose, core::chemical::CENTROID );
		}
	} else { // full atom case
		if ( in_pose.is_centroid() ) { // but pose is centroid
			core::util::switch_to_residue_type_set( in_pose, core::chemical::FA_STANDARD );
		}
	}

	scorefxn_->score( in_pose );

	utility::vector1< core::Size > const res_set_vec( core::pose::get_resnum_list_ordered( string_resnums_, in_pose ) );
	// core::Size chain_seqpos_begin( in_pose.conformation().chain_begin( chain_ ) );
	// core::Size chain_seqpos_end( in_pose.conformation().chain_end( chain_ ) );
	// the neighbor/energy links
	EnergyGraph & energy_graph( in_pose.energies().energy_graph() );
	core::Real score( 0. );
	for ( core::Size i_res_vec = 1; i_res_vec <= res_set_vec.size(); ++i_res_vec ) {
		Size iseq1( res_set_vec[ i_res_vec ] );
		assert( iseq1 <= in_pose.size() );
		//TR<< "Summing energies: seqpos: " << iseq1 << std::endl;
		//search over energy edges
		for ( utility::graph::Graph::EdgeListIter el_iter  = energy_graph.get_node( iseq1 )->edge_list_begin();
				el_iter != energy_graph.get_node( iseq1 )->edge_list_end(); ++el_iter ) {
			EnergyEdge * edge( static_cast< EnergyEdge *> ( *el_iter ) );
			//the other seqpos connected to this edge
			core::Size iseq2( edge->get_first_node_ind() );
			if ( iseq2 == iseq1 ) iseq2 = edge->get_second_node_ind();
			//TR<< "set seqpos: " << iseq1 << " chain: " << in_pose.chain( iseq1 ) << " :: chain seqpos: " << iseq2 << " chain: " << in_pose.chain( iseq2 ) << std::endl;
			//skip if iseq2 is not in the chain we care about
			if ( static_cast< core::Size >( in_pose.chain( iseq2 ) ) != chain_ ) continue;

			// the pair energies cached in the link
			EnergyMap const & emap( edge->fill_energy_map() );
			if ( score_type_ == total_score ) score += emap.dot( scorefxn_->weights() );
			else score += emap[ ScoreType( score_type_ ) ];
		}
	}

	//get wt for unweighted score
	core::Real weight;
	if ( score_type_ == total_score ) weight = 1;
	else weight = scorefxn_->get_weight( ScoreType( score_type_ ) );
	core::Real const weighted_score( weight * score );
	return( weighted_score );
}

}
}
