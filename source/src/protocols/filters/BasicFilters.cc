// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/BasicFilters.cc
/// @brief
/// @details
///  Contains currently:
///
///
/// @author Florian Richter, Sarel Fleishman (sarelf@uw.edu), Rocco Moretti (rmoretti@u.washington.edu)

// Unit Headers
#include <protocols/filters/BasicFilters.hh>
#include <protocols/filters/BasicFilterCreators.hh>

#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

// Package Headers

// Project Headers
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

// ObjexxFCL Headers

// Utility headers
#include <numeric/random/random.hh>

// Boost Headers
#include <boost/foreach.hpp>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

//// C++ headers
static thread_local basic::Tracer TR( "protocols.filters.Filter" );

namespace protocols {
namespace filters {

using namespace core;
typedef std::pair< std::string const, FilterCOP > StringFilter_pair;
typedef utility::tag::TagCOP TagCOP;
typedef core::pose::Pose Pose;

////////////////////////////////////////////////////////////////////////////////////////////////////

StochasticFilter::StochasticFilter() : Filter( "Stochastic" ) {}
StochasticFilter::~StochasticFilter() {}

StochasticFilter::StochasticFilter( core::Real const confidence )
	: Filter( "Stochastic" ), confidence_( confidence )
{}

bool
StochasticFilter::apply( Pose const & ) const
{
	if( confidence_ >= 0.999 ) return true;

	core::Real const random_number( numeric::random::rg().uniform() );
	if( random_number <= confidence_ ) {
		TR<<"stochastic filter returning false"<<std::endl;
		return false;
	}
	TR<<"stochastic filter returning true"<<std::endl;
	return true;
}


FilterOP
StochasticFilter::clone() const
{
	return FilterOP( new StochasticFilter( *this ) );
}

FilterOP
StochasticFilter::fresh_instance() const
{
	return FilterOP( new StochasticFilter() );
}

void
StochasticFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const & )
{
	confidence_ = tag->getOption< core::Real >( "confidence", 1.0 );
	TR<<"stochastic filter with confidence "<<confidence_<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// @brief Used to define a compound logical statement involving other filters with
// AND, OR and XOR
CompoundFilter::CompoundFilter() :
		Filter( "CompoundStatement" ),
		invert_(false),
		reset_filters_(false)
 {}
CompoundFilter::~CompoundFilter() {}

CompoundFilter::CompoundFilter( CompoundStatement const & compound_statement ) :
	Filter( "CompoundStatement" ),
	compound_statement_( compound_statement ),
	invert_(false),
	reset_filters_(false)
{}

bool
CompoundFilter::apply( Pose const & pose ) const
{
	bool const value( compute( pose ) );

	TR<<"Compound logical statement is "<<value<<"."<<std::endl;
	return( value );
}

FilterOP
CompoundFilter::clone() const
{
	return FilterOP( new CompoundFilter( *this ) );
}

FilterOP
CompoundFilter::fresh_instance() const
{
	return FilterOP( new CompoundFilter() );
}

void
CompoundFilter::report( std::ostream & out, Pose const & pose ) const
{
	if( compound_statement_.size() == 2 ){
	//special case for filters that are defined with a confidence value. In that case, we want to report the value of the filter regardless of the stochastic filter
		bool confidence( false );
		CompoundStatement::const_iterator non_stochastic_filter;
		for( CompoundStatement::const_iterator it=compound_statement_.begin(); it!=compound_statement_.end(); ++it ){
			if( it->first->get_type() == "Stochastic" ) confidence = true;
			else non_stochastic_filter = it;
		}
		if( confidence ) non_stochastic_filter->first->report( out, pose );
	}
	bool const value( compute( pose ) );

	out<<"Compound filter returns: "<<value<<'\n';
}

core::Real
CompoundFilter::report_sm( Pose const & pose ) const
{
	if( compound_statement_.size() == 2 ){
	//special case for filters that are defined with a confidence value. In that case, we want to report the value of the filter regardless of the stochastic filter
		bool confidence( false );
		CompoundStatement::const_iterator non_stochastic_filter;
		for( CompoundStatement::const_iterator it=compound_statement_.begin(); it!=compound_statement_.end(); ++it ){
			if( it->first->get_type() == "Stochastic" ) confidence = true;
			else non_stochastic_filter = it;
		}
		if( confidence ) return( non_stochastic_filter->first->report_sm( pose ) );
	}
	bool const value( compute( pose ) );
	return( value );
}

bool
CompoundFilter::compute( Pose const & pose ) const
{

	bool value( true );

	for( CompoundStatement::const_iterator it=compound_statement_.begin(); it!=compound_statement_.end(); ++it ) {
		if( it - compound_statement_.begin() == 0 ){
			// first logical op may only be NOT
			// ANDNOT and ORNOT are also treated as NOT (with a warning)
			value = it->first->apply( pose );
			if (it->second == NOT) value = !value;
			if (it->second == ORNOT) {
				TR << "WARNING: CompoundFilter treating operator ORNOT as NOT" << std::endl;
				value = !value;
			}
			if (it->second == ANDNOT) {
				TR << "WARNING: CompoundFilter treating operator ANDNOT as NOT" << std::endl;
				value = !value;
			}
		} else {
			switch( it->second  ) {
				case ( AND ) : value = value && it->first->apply( pose ); break;
				case ( OR  ) : value = value || it->first->apply( pose ); break;
				case ( XOR ) : value = value ^ it->first->apply( pose ); break;
				case ( ORNOT ) : value = value || !it->first->apply( pose ); break;
				case ( ANDNOT ) : value = value && !it->first->apply( pose ); break;
				case ( NOR ) : value = !( value || it->first->apply( pose ) ); break;
				case (NAND ) : value = !( value && it->first->apply( pose ) ); break;
				case (NOT ) :
					TR << "WARNING: CompoundFilter treating operator NOT as ANDNOT" << std::endl;
					value = value && !it->first->apply( pose );
					break;
			}
		}
	}
	if( invert_ ) value = !value;
	return( value );
}

void
CompoundFilter::clear()
{
	compound_statement_.clear();
}

CompoundFilter::iterator
CompoundFilter::begin()
{
	return( compound_statement_.begin() );
}
CompoundFilter::const_iterator
CompoundFilter::begin() const
{
	return( compound_statement_.begin() );
}

CompoundFilter::iterator
CompoundFilter::end()
{
	return( compound_statement_.end() );
}

CompoundFilter::const_iterator
CompoundFilter::end() const
{
	return( compound_statement_.end() );
}

void
CompoundFilter::invert( bool const inv )
{
	invert_ = inv;
}

void
CompoundFilter::set_reset_filters( utility::vector1<FilterOP> const reset_filters ) 
{
	reset_filters_ = reset_filters;
}

void
CompoundFilter::reset_filters()
{
	CompoundStatement new_compound_statement;
	std::pair< FilterOP, boolean_operations > filter_pair;
	for( CompoundStatement::const_iterator it=compound_statement_.begin(); it!=compound_statement_.end(); ++it ) {
		filter_pair.first = it->first;
		filter_pair.second = it->second;
		BOOST_FOREACH( FilterOP filter, reset_filters_) { 
			if( filter->get_user_defined_name() == it->first->get_user_defined_name() ) {
				filter_pair.first = filter->clone();
				break;
			}
		}
		new_compound_statement.push_back( filter_pair );
	}
	compound_statement_.clear();
	compound_statement_ = new_compound_statement;
}

void
CompoundFilter::clear_reset_filters()
{
	reset_filters_.clear();
}

/// @details call the compound statement's constituent filters' set_resid
void
CompoundFilter::set_resid( core::Size const resid )
{
	for( iterator it( compound_statement_.begin() ); it!=compound_statement_.end(); ++it )
		protocols::moves::modify_ResId_based_object( it->first, resid );
}

void
CompoundFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const & filters,
	moves::Movers_map const &,
	Pose const & )
{
	TR<<"CompoundStatement"<<std::endl;
	invert_ = tag->getOption<bool>( "invert", false );

	BOOST_FOREACH(TagCOP cmp_tag_ptr, tag->getTags() ){
		std::string const operation( cmp_tag_ptr->getName() );
		std::pair< FilterOP, boolean_operations > filter_pair;
		if( operation == "AND" ) filter_pair.second = AND;
		else if( operation == "OR" ) filter_pair.second = OR;
		else if( operation == "XOR" ) filter_pair.second = XOR;
		else if( operation == "NOR" ) filter_pair.second = NOR;
		else if( operation == "NAND" ) filter_pair.second = NAND;
		else if( operation == "ORNOT" ) filter_pair.second = ORNOT;
		else if( operation == "ANDNOT" ) filter_pair.second = ANDNOT;
		else if( operation == "NOT" ) filter_pair.second = NOT;
		else {
			throw utility::excn::EXCN_RosettaScriptsOption( "Error: Boolean operation in tag is undefined." );
		}
		std::string const filter_name( cmp_tag_ptr->getOption<std::string>( "filter_name" ) );

		Filters_map::const_iterator find_filter( filters.find( filter_name ));
		bool const filter_found( find_filter!=filters.end() );
		if( filter_found )
			filter_pair.first = find_filter->second->clone();
		else {
			TR<<"***WARNING WARNING! Filter defined for CompoundStatement not found in filter_list!!!! Defaulting to truefilter***"<<std::endl;
			filter_pair.first = FilterOP( new filters::TrueFilter );
		}
		runtime_assert( filter_found );
		compound_statement_.push_back( filter_pair );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// @brief Used to combine multiple seperate filters into a single filter value
	CombinedFilter::CombinedFilter() : Filter( "CombinedValue" ), threshold_(0.0) {}
CombinedFilter::~CombinedFilter() {}

bool
CombinedFilter::apply( core::pose::Pose const & pose ) const
{
	core::Real const value( compute( pose ) );

	TR<<"CombinedFilter value is "<<value<<", threshold is "<< threshold_ << "."<<std::endl;
	return( value <= threshold_ );
}

FilterOP
CombinedFilter::clone() const
{
	return FilterOP( new CombinedFilter( *this ) );
}

FilterOP
CombinedFilter::fresh_instance() const
{
	return FilterOP( new CombinedFilter() );
}

void
CombinedFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	Real value( compute( pose ) );

	out<<"Combined filter returns: "<<value<<'\n';
}

core::Real
CombinedFilter::report_sm( core::pose::Pose const & pose ) const
{
	return compute(pose);
}

core::Real
CombinedFilter::compute( core::pose::Pose const & pose ) const
{
	core::Real value( 0.0 );

	BOOST_FOREACH(FilterWeightPair fw_pair, filterlist_){
		value += fw_pair.second * fw_pair.first->report_sm( pose );
	}
	return( value );
}

void
CombinedFilter::set_reset_filters( utility::vector1<FilterOP> const reset_filters ) 
{
	reset_filters_ = reset_filters;
}

void
CombinedFilter::reset_filters()
{
	FilterList new_filterlist;
	FilterWeightPair filter_pair;
	BOOST_FOREACH(FilterWeightPair fw_pair, filterlist_){
		filter_pair.first = fw_pair.first;
		filter_pair.second = fw_pair.second;
		BOOST_FOREACH( FilterOP filter, reset_filters_) { 
			if( filter->get_user_defined_name() == fw_pair.first->get_user_defined_name() ) {
				filter_pair.first = filter->clone();
				break;
			}
		}
		new_filterlist.push_back( filter_pair );
	}
	filterlist_.clear();
	filterlist_ = new_filterlist;
}

void
CombinedFilter::clear_reset_filters()
{
	reset_filters_.clear();
}

void
CombinedFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const & filters,
	moves::Movers_map const &,
	Pose const & )
{
	threshold_ = tag->getOption<core::Real>( "threshold", 0.0 );
	utility::vector1< TagCOP > const sub_tags( tag->getTags() );
	BOOST_FOREACH(TagCOP tag_ptr, sub_tags){
		core::Real weight(1.0);
		if (tag_ptr->hasOption("factor") ) {
			weight = tag_ptr->getOption<core::Real>( "factor" );
		} else if (tag_ptr->hasOption("temp") ) {
			weight = 1.0 / tag_ptr->getOption<core::Real>( "temp" );
		}

		FilterOP filter;
		std::string const filter_name( tag_ptr->getOption<std::string>( "filter_name" ) );
		Filters_map::const_iterator find_filter( filters.find( filter_name ));
		bool const filter_found( find_filter!=filters.end() );
		if( filter_found ) {
			filter = find_filter->second->clone();
		}
		else {
			TR.Warning<<"***WARNING WARNING! Filter " << filter_name << " defined for CombinedValue not found in filter_list!!!! ***"<<std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption("Filter "+filter_name+" not found in filter list.");
		}
		filterlist_.push_back( FilterWeightPair(filter, weight) );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Apply a sub-mover prior to calculating a filter value
MoveBeforeFilter::MoveBeforeFilter() : Filter( "MoveBeforeFilter" ) {}
MoveBeforeFilter::MoveBeforeFilter(moves::MoverOP mover, FilterCOP filter) : Filter( "MoveBeforeFilter" ), subfilter_(filter), submover_(mover) {}
MoveBeforeFilter::~MoveBeforeFilter() {}

bool
MoveBeforeFilter::apply( core::pose::Pose const & pose ) const
{
	core::pose::Pose p_mod(pose);
	submover_->apply(p_mod);
	return subfilter_->apply(p_mod);
}

FilterOP
MoveBeforeFilter::clone() const
{
	return FilterOP( new MoveBeforeFilter( *this ) );
}

FilterOP
MoveBeforeFilter::fresh_instance() const
{
	return FilterOP( new MoveBeforeFilter() );
}

void
MoveBeforeFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	core::pose::Pose p_mod(pose);
	submover_->apply(p_mod);
	subfilter_->report(out, p_mod);
}

core::Real
MoveBeforeFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::pose::Pose p_mod(pose);
	submover_->apply(p_mod);
	return subfilter_->report_sm(p_mod);
}

void
MoveBeforeFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const & filters,
	moves::Movers_map const & movers,
	Pose const & )
{
	std::string mover_name("");
	std::string filter_name("");
	if ( tag->hasOption("mover") ) mover_name = tag->getOption< std::string >( "mover" );
	if ( tag->hasOption("mover_name") ) mover_name = tag->getOption< std::string >( "mover_name" );
	if ( tag->hasOption("filter") ) filter_name = tag->getOption< std::string >( "filter" );
	if ( tag->hasOption("filter_name") ) filter_name = tag->getOption< std::string >( "filter_name" );

	moves::Movers_map::const_iterator  find_mover ( movers.find( mover_name ));
  Filters_map::const_iterator find_filter( filters.find( filter_name ));

  if( find_mover == movers.end() ) {
    TR.Error << "ERROR !! mover '"<<mover_name<<"' not found in map: \n" << tag << std::endl;
    runtime_assert( find_mover != movers.end() );
  }
  if( find_filter == filters.end() ) {
    TR.Error << "ERROR !! filter '"<<filter_name<<"' not found in map: \n" << tag << std::endl;
    runtime_assert( find_filter != filters.end() );
  }
	submover_ = find_mover->second;
	subfilter_ = find_filter->second;

	TR << "Setting MoveBeforeFilter for mover '"<<mover_name<<"' and filter '"<<filter_name<<"'"<<std::endl;

}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Evaluate to a value contingent on the evaluation of another filter.

IfThenFilter::IfThenFilter() :
		Filter( "IfThenFilter" ),
		elsevalue_(0),
		threshold_(0),
		floor_(false)
 {}
IfThenFilter::~IfThenFilter() {}

void
IfThenFilter::add_condition( FilterCOP testfilter, FilterCOP valuefilter, core::Real value, bool invert, core::Real weight) {
	runtime_assert( iffilters_.size() == thenfilters_.size() );
	runtime_assert( iffilters_.size() == values_.size() );
	runtime_assert( iffilters_.size() == weights_.size() );
	runtime_assert( iffilters_.size() == invert_.size() );

	iffilters_.push_back( testfilter );
	thenfilters_.push_back( valuefilter );
	values_.push_back( value );
	weights_.push_back( weight );
	invert_.push_back( invert );
}

void
IfThenFilter::set_else( FilterCOP elsefilter, core::Real value, core::Real elseweight ) {
	elsefilter_ = elsefilter;
	elsevalue_ = value;
	elseweight_ = elseweight;
}

bool
IfThenFilter::apply( core::pose::Pose const & pose ) const
{
	core::Real const value( compute( pose ) );

	TR<<"IfThenFilter value is "<<value<<", threshold is "<< threshold_;
	if( floor_ ) {
		TR << " (lower limit)" << std::endl;
		return value >= threshold_;
	} else {
		TR << " (upper limit)" << std::endl;
		return value <= threshold_;
	}
}

FilterOP
IfThenFilter::clone() const
{
	return FilterOP( new IfThenFilter( *this ) );
}

FilterOP
IfThenFilter::fresh_instance() const
{
	return FilterOP( new IfThenFilter() );
}

void
IfThenFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	Real value( compute( pose ) );

	out<<"IfThen filter returns: "<<value<<'\n';
}

core::Real
IfThenFilter::report_sm( core::pose::Pose const & pose ) const
{
	return compute(pose);
}

core::Real
IfThenFilter::compute( core::pose::Pose const & pose ) const
{
	assert( iffilters_.size() == thenfilters_.size() );
	assert( iffilters_.size() == values_.size() );
	assert( iffilters_.size() == weights_.size() );
	assert( iffilters_.size() == invert_.size() );

	for( core::Size ii(1); ii<= iffilters_.size(); ++ii ) {
		assert( iffilters_[ii] );
		if( invert_[ii] ^ iffilters_[ii]->apply( pose ) ) { // XOR: if invert_ is true, true becomes false and false becomes true
			if( thenfilters_[ii] ) {
				return weights_[ii] * thenfilters_[ii]->report_sm( pose );
			} else {
				return weights_[ii] * values_[ii];
			}
		}
	}
	if( iffilters_.size() == 0 ) {
		TR.Warning << "WARNING: No conditional filters specified for IfThenFilter. Using else values only." << std::endl;
	}
	if( elsefilter_ ) {
		return elseweight_ * elsefilter_->report_sm( pose );
	}

	return elseweight_ * elsevalue_;
}


void
IfThenFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const & filters,
	moves::Movers_map const &,
	Pose const & )
{
	threshold( tag->getOption<core::Real>( "threshold", 0.0 ) );
	set_lower_threshold( tag->getOption<bool>( "lower_threshold", false ) );
	if( tag->hasOption( "lower_threshold" ) && ! tag->hasOption( "threshold" ) ){
		TR.Warning << "WARNING: In IfThenFilter, lower_threshold set without setting threshold." << std::endl;
		TR.Warning << "WARNING: Note that lower_threshold is a true/false flag, not a real-valued setting." << std::endl;
	}
	utility::vector1< TagCOP > const sub_tags( tag->getTags() );
	BOOST_FOREACH(TagCOP tag_ptr, sub_tags){
		std::string const tagname = tag_ptr->getName();
		FilterOP valuefilter = 0; //default NULL
		if( tag_ptr->hasOption("valuefilter") ) {
			valuefilter = protocols::rosetta_scripts::parse_filter( tag_ptr->getOption<std::string>( "valuefilter" ), filters);
		}
		core::Real value( tag_ptr->getOption<core::Real>( "value", 0 ) );
		core::Real weight( tag_ptr->getOption<core::Real>( "weight", 1 ) );

		if( tagname == "IF" || tagname == "ELIF" ){
			if( ! tag_ptr->hasOption("testfilter") )  {
				TR.Error << "In IfThenFilter, If and ELIF require a tesfilter option." << std::endl;
				throw utility::excn::EXCN_RosettaScriptsOption("In IfThenFilter, If and ELIF require a tesfilter option.");
			}
			FilterOP testfilter = protocols::rosetta_scripts::parse_filter( tag_ptr->getOption<std::string>( "testfilter" ), filters);
			bool inverttest = tag_ptr->getOption< bool >( "inverttest", false );
			add_condition( testfilter, valuefilter, value, inverttest, weight );
		} else if ( tagname == "ELSE" ) {
			set_else( valuefilter, value, weight );
		} else {
			TR.Error << "Unknown subtag name in IfThenFilter: " << tagname << std::endl;
			TR.Error << "   Acceptable values are:   IF   ELIF   ELSE" << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption("Unknown subtag name in IfThenFilter: " + tagname );
		}
	}
	TR << "IfThenFilter defined with " << iffilters_.size() << " conditions and "
			<< ( elsefilter_ ? "an ": "no " ) << " else filter and a "
			<< ( floor_ ? "lower " : "upper " ) << "threshold of " << threshold_ << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// @brief FilterCreator methods

FilterOP
TrueFilterCreator::create_filter() const { return FilterOP( new TrueFilter ); }

std::string
TrueFilterCreator::keyname() const { return "TrueFilter"; }

FilterOP
FalseFilterCreator::create_filter() const { return FilterOP( new FalseFilter ); }

std::string
FalseFilterCreator::keyname() const { return "FalseFilter"; }

FilterOP
StochasticFilterCreator::create_filter() const { return FilterOP( new StochasticFilter ); }

std::string
StochasticFilterCreator::keyname() const { return "Stochastic"; }

FilterOP
CompoundFilterCreator::create_filter() const { return FilterOP( new CompoundFilter ); }

std::string
CompoundFilterCreator::keyname() const { return "CompoundStatement"; }

FilterOP
CombinedFilterCreator::create_filter() const { return FilterOP( new CombinedFilter ); }

std::string
CombinedFilterCreator::keyname() const { return "CombinedValue"; }

FilterOP
MoveBeforeFilterCreator::create_filter() const { return FilterOP( new MoveBeforeFilter ); }

std::string
MoveBeforeFilterCreator::keyname() const { return "MoveBeforeFilter"; }

FilterOP
IfThenFilterCreator::create_filter() const { return FilterOP( new IfThenFilter ); }

std::string
IfThenFilterCreator::keyname() const { return "IfThenFilter"; }


} // filters
} // protocols
