// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/ScoreCutoffFilter.cc
/// @brief
/// @detailed
///	  Contains currently:
///
///
/// @author Florian Richter

// Unit Headers
#include <protocols/simple_filters/ScoreCutoffFilter.hh>
#include <protocols/simple_filters/ScoreCutoffFilterCreator.hh>

// Package Headers

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// Utility headers
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/scoring/EnergyGraph.hh>



//// C++ headers
static basic::Tracer tr("protocols.filters.ScoreCutoffFilter");

namespace protocols {
namespace simple_filters {


ScoreCutoffFilter::ScoreCutoffFilter() :
	parent("ScoreCutoffFilter"), cutoff_( 0.0 ), report_residue_pair_energies_( false ), total_score_( true ), unweighted_( false )
{
	score_types_.push_back( core::scoring::total_score );
	positions_.clear();
}

ScoreCutoffFilter::ScoreCutoffFilter( core::Real cutoff_in ) :
	parent("ScoreCutoffFilter"),cutoff_(cutoff_in), report_residue_pair_energies_( false ), total_score_( true ), unweighted_( false )
{
	score_types_.push_back( core::scoring::total_score );
	positions_.clear();
}


void
ScoreCutoffFilter::set_score_type( core::scoring::ScoreType scotype )
{
	score_types_.clear();
	add_score_type( scotype );
}

void
ScoreCutoffFilter::add_score_type( core::scoring::ScoreType scotype )
{
	using namespace core::scoring;

	if( scotype == total_score ){

		total_score_ = true;
		if ( score_types_.size() != 0 ) {
			if( (score_types_.size() > 1) || ( score_types_[1] != total_score ) ){

				score_types_.clear();
				tr << "WARNING: setting score_type to total score even though other score types have previously been set. these will be erased." << std::endl;
			}
		}
	}
	else{

		if( total_score_ == true ){

			score_types_.clear();
			total_score_ = false;
		}
	}

	score_types_.push_back( scotype );
}


bool
ScoreCutoffFilter::apply( core::pose::Pose const & pose ) const {

	core::Real cur_score = get_score( pose );

	if ( cur_score <= cutoff() ) return true;

	return false;

} // apply_filter


core::Real
ScoreCutoffFilter::get_score( core::pose::Pose const & pose ) const {

	using namespace core::scoring;
	EnergyMap emap;

	assert( score_types_.size() > 0 );

	if ( positions_.size() == 0 ) { //means we're interested in pose totals

		for ( utility::vector1< ScoreType >::const_iterator sco_it = score_types_.begin();
				 sco_it != score_types_.end(); ++sco_it ){

			emap[ *sco_it ] = pose.energies().total_energies()[ *sco_it ];
		}
	}

	else {
		for ( utility::vector1< core::Size >::const_iterator pos_it = positions_.begin(); pos_it != positions_.end(); ++pos_it ) {

			for ( utility::vector1< ScoreType >::const_iterator sco_it = score_types_.begin();
						 sco_it != score_types_.end(); ++sco_it ){
				emap[ *sco_it ] += pose.energies().residue_total_energies( *pos_it )[ *sco_it ];
			}
		}
	} // else

	//now accumulate into score: total, unweighted, or weighted
	if ( total_score_ ) return emap[ total_score ];
	if ( unweighted_ )  return emap.sum();
	return emap.dot( pose.energies().weights() );

}

void
ScoreCutoffFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & , protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const &  )
{
	if (tag->hasOption("report_residue_pair_energies")) {
		if( tag->getOption<core::Size>("report_residue_pair_energies",0) == 1 )	report_residue_pair_energies_ = true;
	}
	if (tag->hasOption("cutoff")) {
		cutoff_ = tag->getOption<core::Real>("cutoff", 10000.0 );
	}
}

void
ScoreCutoffFilter::report( std::ostream & ostr, core::pose::Pose const & pose ) const
{
	if( report_residue_pair_energies_ ) output_residue_pair_energies( ostr, pose );
	else {
		using namespace core::scoring;
		for ( utility::vector1< ScoreType >::const_iterator sco_it = score_types_.begin();
				sco_it != score_types_.end(); ++sco_it ) {
			EnergyMap emap;
			emap[ *sco_it ] = pose.energies().total_energies()[ *sco_it ];
			ostr << " Scoretype: " << *sco_it << " score: " << emap.sum() << ", cutoff: " << cutoff_ << std::endl;
		}
	}
}

void
ScoreCutoffFilter::output_residue_pair_energies( std::ostream & ostr, core::pose::Pose const & pose ) const
{
	using namespace core::scoring;
	core::Size field_width( 10 ), last_active_sr_st(0);
	utility::vector1< ScoreType > active_st;
	EnergyMap const & weights( pose.energies().weights() );
	for( core::Size i = 1; i <= n_score_types; ++i ){
		if( weights[ ScoreType(i)] != 0 ) active_st.push_back( ScoreType(i) );
		if( i == n_shortranged_2b_score_types ) last_active_sr_st = active_st.size();
	}
	utility::vector1< LREnergyContainerCOP > active_lr_e;
	for( Size lr = 1; lr <= methods::n_long_range_types; lr++){
		methods::LongRangeEnergyType lr_type = methods::LongRangeEnergyType( lr );
		LREnergyContainerCOP lrec = pose.energies().long_range_container( lr_type );
		if( !lrec ) continue;
		if(  !lrec->empty() ) active_lr_e.push_back( lrec );
	}

	ostr << "\nResResE " <<  ObjexxFCL::format::A( field_width, "Res1" ) << " " << ObjexxFCL::format::A(field_width, "Res2") << " ";
	for( core::Size i =1; i <= active_st.size(); ++i) ostr << ObjexxFCL::format::A( field_width , name_from_score_type( active_st[i] ) ) << " ";
	ostr << ObjexxFCL::format::A( field_width , "total" )<< " \n";

	ostr << "ResResE " << ObjexxFCL::format::A( field_width, "nonzero" ) << " " << ObjexxFCL::format::A( field_width, "weights" ) << " ";
	for( core::Size i =1; i <= active_st.size(); ++i) ostr << ObjexxFCL::format::F( field_width, 2, weights[ active_st[ i ] ] ) << " ";
	ostr << ObjexxFCL::format::A( field_width , "NA" ) << "\n";

	EnergyGraph const & egraph( pose.energies().energy_graph() );

	for( core::Size res1 = 1; res1 <= pose.total_residue(); ++res1 ){
		std::string res1name( pose.residue_type( res1 ).name1() + utility::to_string( res1 ) );
		std::map< core::Size, EnergyMap > upper_interactions;

		//1st, get short range enegires
		for( core::graph::EdgeListConstIterator egraph_it = egraph.get_energy_node( res1 )->const_upper_edge_list_begin();
				 egraph_it != egraph.get_energy_node( res1 )->const_upper_edge_list_end(); ++egraph_it){
			core::Size other_ind( (*egraph_it)->get_other_ind( res1 ) );
			upper_interactions.insert( std::pair< Size, EnergyMap >(other_ind, EnergyMap() ) );
			EnergyMap & this_emap( upper_interactions.find( other_ind )->second );
			EnergyEdge const * eedge( static_cast< EnergyEdge const * >(*egraph_it));

			for( core::Size i =1; i <= last_active_sr_st; ++i) this_emap[ active_st[i] ] = ((*eedge)[  active_st[i] ]);

		} //loop over all short range upper edges for res

		//2nd, get long range energies
		for( Size lr = 1; lr <= active_lr_e.size(); lr++){
			for ( ResidueNeighborConstIteratorOP rni( active_lr_e[lr]->const_upper_neighbor_iterator_begin( res1 ) );
						*rni != *( active_lr_e[lr]->const_upper_neighbor_iterator_end( res1 ) ); ++(*rni) ) {
						EnergyMap lr_emap;
						core::Size other_ind( rni->upper_neighbor_id() );
						rni->retrieve_energy( lr_emap );
						std::map< core::Size, EnergyMap >::iterator map_it = upper_interactions.find( other_ind );
						if( map_it == upper_interactions.end() ) upper_interactions.insert( std::pair<Size, EnergyMap >(other_ind, lr_emap ) );
						else map_it->second += lr_emap;
			} // loop over all lr interactions for this res
		} //loop over all lr energies

		//3. finally output shit
		for( std::map<  core::Size, EnergyMap >::const_iterator map_it = upper_interactions.begin(); map_it != upper_interactions.end(); ++map_it){
			core::Size res2 = map_it->first;
			EnergyMap const & this_emap( map_it->second );
			std::string res2name( pose.residue_type( res2 ).name1() + utility::to_string( res2 ) );
			ostr << "ResResE " << ObjexxFCL::format::A( field_width, res1name ) << " " << ObjexxFCL::format::A( field_width, res2name ) << " ";
			core::Real totscore(0.0);
			for( core::Size i =1; i <=  active_st.size(); ++i){
				core::Real score(  weights[ active_st[i] ] * this_emap[ active_st[i] ] );
				totscore += score;
				ostr << ObjexxFCL::format::F( field_width, 2, score) << " ";
			}
			ostr << ObjexxFCL::format::F( field_width, 2, totscore) << "\n";
		}
	} //loop over all residues

}

filters::FilterOP
ScoreCutoffFilterCreator::create_filter() const { return new ScoreCutoffFilter; }

std::string
ScoreCutoffFilterCreator::keyname() const { return "ScoreCutoffFilter"; }



} // filters
} // protocols
