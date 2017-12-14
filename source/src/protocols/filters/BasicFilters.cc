// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

#include <utility>
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

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

//// C++ headers
static basic::Tracer TR( "protocols.filters.Filter" );

namespace protocols {
namespace filters {

using namespace core;
typedef std::pair< std::string const, FilterCOP > StringFilter_pair;
using TagCOP = utility::tag::TagCOP;
using Pose = core::pose::Pose;

////////////////////////////////////////////////////////////////////////////////////////////////////

StochasticFilter::StochasticFilter() : Filter( "Stochastic" ) {}
StochasticFilter::~StochasticFilter() = default;

StochasticFilter::StochasticFilter( core::Real const confidence, FilterOP subfilter, bool run_subfilter_on ):
	Filter( "Stochastic" ),
	confidence_( confidence ),
	subfilter_(std::move( subfilter )),
	run_subfilter_on_( run_subfilter_on )
{}

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

bool
StochasticFilter::calculate() const {
	if ( confidence_ < 0.999 ) {
		core::Real const random_number( numeric::random::rg().uniform() );
		return ( random_number <= confidence_ );
	} else {
		return true;
	}
}

bool
StochasticFilter::apply( Pose const & pose ) const
{
	bool result = calculate();

	if ( subfilter_ && result == run_subfilter_on_ ) {
		return subfilter_->apply( pose );
	} else {
		return result;
	}
}


core::Real
StochasticFilter::report_sm( core::pose::Pose const & pose ) const {
	if ( subfilter_ ) {
		return subfilter_->report_sm( pose ); // Always use the sub-filter's metric.
	} else {
		return Filter::report_sm( pose ); // Use the base class behavior.
	}
}

void
StochasticFilter::report( std::ostream & ostream, core::pose::Pose const & pose) const {
	if ( subfilter_ ) {
		subfilter_->report( ostream, pose ); // Always use the sub-filters behavior
	} else {
		// Do nothing.
	}
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

///////////////////////////////////////////////////////////////////////////////////////////////////
// @brief Used to define a compound logical statement involving other filters with
// AND, OR and XOR
CompoundFilter::CompoundFilter() :
	Filter( "CompoundStatement" ),
	invert_(false),
	reset_filters_(false)
{}
CompoundFilter::~CompoundFilter() = default;

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
	bool const value( compute( pose ) );

	out<<"Compound filter returns: "<<value<<'\n';
}

core::Real
CompoundFilter::report_sm( Pose const & pose ) const
{
	bool const value( compute( pose ) );
	return( value );
}

bool
CompoundFilter::compute( Pose const & pose ) const
{

	bool value( true );

	for ( auto it=compound_statement_.begin(); it!=compound_statement_.end(); ++it ) {
		if ( it - compound_statement_.begin() == 0 ) {
			// first logical op may only be NOT
			// ANDNOT and ORNOT are also treated as NOT (with a warning)
			value = compute_subfilter( it->first, pose );
			if ( it->second == NOT ) value = !value;
			if ( it->second == ORNOT ) {
				TR.Warning << "CompoundFilter treating operator ORNOT as NOT" << std::endl;
				value = !value;
			}
			if ( it->second == ANDNOT ) {
				TR.Warning << "CompoundFilter treating operator ANDNOT as NOT" << std::endl;
				value = !value;
			}
		} else {
			switch( it->second  ) {
			case ( AND ) : value = value && compute_subfilter( it->first, pose ); break;
			case ( OR  ) : value = value || compute_subfilter( it->first, pose ); break;
			case ( XOR ) : value = value ^ compute_subfilter( it->first, pose ); break;
			case ( ORNOT ) : value = value || !compute_subfilter( it->first, pose ); break;
			case ( ANDNOT ) : value = value && !compute_subfilter( it->first, pose ); break;
			case ( NOR ) : value = !( value || compute_subfilter( it->first, pose ) ); break;
			case (NAND ) : value = !( value && compute_subfilter( it->first, pose ) ); break;
			case (NOT ) :
				TR.Warning << "CompoundFilter treating operator NOT as ANDNOT" << std::endl;
				value = value && !compute_subfilter( it->first, pose );
				break;
			}
		}
	}
	if ( invert_ ) value = !value;
	TR << "CompoundFilter as a whole " << (value ? " succeeded." : " failed." ) << std::endl;
	return( value );
}

bool
CompoundFilter::compute_subfilter( FilterOP const & filter, Pose const & pose ) const {
	bool subfilter_value( filter->apply( pose ) );
	TR << "CompoundFilter subfilter "
		<< (filter->get_user_defined_name().empty() ? filter->get_type() : filter->get_user_defined_name() )
		<< (subfilter_value ? " succeeded." : " failed." ) << std::endl;
	return subfilter_value;
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
CompoundFilter::set_reset_filters( utility::vector1<FilterOP> const & reset_filters )
{
	reset_filters_ = reset_filters;
}

void
CompoundFilter::reset_filters()
{
	CompoundStatement new_compound_statement;
	std::pair< FilterOP, boolean_operations > filter_pair;
	for ( CompoundStatement::const_iterator it=compound_statement_.begin(); it!=compound_statement_.end(); ++it ) {
		filter_pair.first = it->first;
		filter_pair.second = it->second;
		for ( FilterOP filter : reset_filters_ ) {
			if ( filter->get_user_defined_name() == it->first->get_user_defined_name() ) {
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
	for ( auto & it : compound_statement_ ) {
		protocols::moves::modify_ResId_based_object( it.first, resid );
	}
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

	for ( TagCOP cmp_tag_ptr : tag->getTags() ) {
		std::string const operation( cmp_tag_ptr->getName() );
		std::pair< FilterOP, boolean_operations > filter_pair;
		if ( operation == "AND" ) filter_pair.second = AND;
		else if ( operation == "OR" ) filter_pair.second = OR;
		else if ( operation == "XOR" ) filter_pair.second = XOR;
		else if ( operation == "NOR" ) filter_pair.second = NOR;
		else if ( operation == "NAND" ) filter_pair.second = NAND;
		else if ( operation == "ORNOT" ) filter_pair.second = ORNOT;
		else if ( operation == "ANDNOT" ) filter_pair.second = ANDNOT;
		else if ( operation == "NOT" ) filter_pair.second = NOT;
		else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Error: Boolean operation in tag is undefined." );
		}
		std::string const filter_name( cmp_tag_ptr->getOption<std::string>( "filter_name" ) );

		auto find_filter( filters.find( filter_name ));
		bool const filter_found( find_filter!=filters.end() );
		if ( filter_found ) {
			filter_pair.first = find_filter->second->clone();
		} else {
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
CombinedFilter::~CombinedFilter() = default;

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

	for ( FilterWeightPair const & fw_pair : filterlist_ ) {
		value += fw_pair.second * fw_pair.first->report_sm( pose );
	}
	return( value );
}

void
CombinedFilter::set_reset_filters( utility::vector1<FilterOP> const & reset_filters )
{
	reset_filters_ = reset_filters;
}

void
CombinedFilter::reset_filters()
{
	FilterList new_filterlist;
	FilterWeightPair filter_pair;
	for ( FilterWeightPair const & fw_pair : filterlist_ ) {
		filter_pair.first = fw_pair.first;
		filter_pair.second = fw_pair.second;
		for ( FilterOP filter : reset_filters_ ) {
			if ( filter->get_user_defined_name() == fw_pair.first->get_user_defined_name() ) {
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
	set_threshold( tag->getOption<core::Real>( "threshold", 0.0 ) );
	utility::vector1< TagCOP > const sub_tags( tag->getTags() );
	for ( TagCOP tag_ptr : sub_tags ) {
		core::Real weight(1.0);
		if ( tag_ptr->hasOption("factor") ) {
			weight = tag_ptr->getOption<core::Real>( "factor" );
		} else if ( tag_ptr->hasOption("temp") ) {
			weight = 1.0 / tag_ptr->getOption<core::Real>( "temp" );
		}

		FilterOP filter;
		std::string const filter_name( tag_ptr->getOption<std::string>( "filter_name" ) );
		auto find_filter( filters.find( filter_name ));
		bool const filter_found( find_filter!=filters.end() );
		if ( filter_found ) {
			filter = find_filter->second->clone();
		} else {
			TR.Warning<<"***Filter " << filter_name << " defined for CombinedValue not found in filter_list!!!! ***"<<std::endl;
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Filter "+filter_name+" not found in filter list.");
		}
		add_filter( filter, weight, false /*filter was already cloned*/ );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Apply a sub-mover prior to calculating a filter value
MoveBeforeFilter::MoveBeforeFilter() : Filter( "MoveBeforeFilter" ) {}
MoveBeforeFilter::MoveBeforeFilter(moves::MoverOP mover, FilterCOP filter) : Filter( "MoveBeforeFilter" ), subfilter_(std::move(filter)), submover_(std::move(mover)) {}
MoveBeforeFilter::~MoveBeforeFilter() = default;

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

	auto  find_mover ( movers.find( mover_name ));
	auto find_filter( filters.find( filter_name ));

	if ( find_mover == movers.end() ) {
		TR.Error << "Mover '"<<mover_name<<"' not found in map: \n" << tag << std::endl;
		runtime_assert( find_mover != movers.end() );
	}
	if ( find_filter == filters.end() ) {
		TR.Error << "Filter '"<<filter_name<<"' not found in map: \n" << tag << std::endl;
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
IfThenFilter::~IfThenFilter() = default;

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
	if ( floor_ ) {
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
	debug_assert( iffilters_.size() == thenfilters_.size() );
	debug_assert( iffilters_.size() == values_.size() );
	debug_assert( iffilters_.size() == weights_.size() );
	debug_assert( iffilters_.size() == invert_.size() );

	for ( core::Size ii(1); ii<= iffilters_.size(); ++ii ) {
		debug_assert( iffilters_[ii] );
		if ( invert_[ii] ^ iffilters_[ii]->apply( pose ) ) { // XOR: if invert_ is true, true becomes false and false becomes true
			if ( thenfilters_[ii] ) {
				return weights_[ii] * thenfilters_[ii]->report_sm( pose );
			} else {
				return weights_[ii] * values_[ii];
			}
		}
	}
	if ( iffilters_.size() == 0 ) {
		TR.Warning << "No conditional filters specified for IfThenFilter. Using else values only." << std::endl;
	}
	if ( elsefilter_ ) {
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
	if ( tag->hasOption( "lower_threshold" ) && ! tag->hasOption( "threshold" ) ) {
		TR.Warning << "In IfThenFilter, lower_threshold set without setting threshold." << std::endl;
		TR.Warning << "Note that lower_threshold is a true/false flag, not a real-valued setting." << std::endl;
	}
	utility::vector1< TagCOP > const sub_tags( tag->getTags() );
	for ( TagCOP tag_ptr : sub_tags ) {
		std::string const tagname = tag_ptr->getName();
		FilterOP valuefilter = nullptr; //default NULL
		if ( tag_ptr->hasOption("valuefilter") ) {
			valuefilter = protocols::rosetta_scripts::parse_filter( tag_ptr->getOption<std::string>( "valuefilter" ), filters);
		}
		auto value( tag_ptr->getOption<core::Real>( "value", 0 ) );
		auto weight( tag_ptr->getOption<core::Real>( "weight", 1 ) );

		if ( tagname == "IF" || tagname == "ELIF" ) {
			if ( ! tag_ptr->hasOption("testfilter") )  {
				TR.Error << "In IfThenFilter, If and ELIF require a tesfilter option." << std::endl;
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "In IfThenFilter, If and ELIF require a tesfilter option.");
			}
			FilterOP testfilter = protocols::rosetta_scripts::parse_filter( tag_ptr->getOption<std::string>( "testfilter" ), filters);
			bool inverttest = tag_ptr->getOption< bool >( "inverttest", false );
			add_condition( testfilter, valuefilter, value, inverttest, weight );
		} else if ( tagname == "ELSE" ) {
			set_else( valuefilter, value, weight );
		} else {
			TR.Error << "Unknown subtag name in IfThenFilter: " << tagname << std::endl;
			TR.Error << "   Acceptable values are:   IF   ELIF   ELSE" << std::endl;
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Unknown subtag name in IfThenFilter: " + tagname );
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

void TrueFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TrueFilter::provide_xml_schema( xsd );
}

std::string TrueFilter::name() const {
	return class_name();
}

std::string TrueFilter::class_name() {
	return "TrueFilter";
}

void TrueFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist; // No attributes

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Filter that always returns true", attlist );
}

std::string FalseFilter::name() const {
	return class_name();
}

std::string FalseFilter::class_name() {
	return "FalseFilter";
}

void FalseFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist; // No attributes

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Filter that always returns false", attlist );
}

std::string FalseFilterCreator::keyname() const {
	return FalseFilter::class_name();
}

protocols::filters::FilterOP
FalseFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new FalseFilter );
}

void FalseFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FalseFilter::provide_xml_schema( xsd );
}


// XRW TEMP StochasticFilterCreator::create_filter() const { return FilterOP( new StochasticFilter ); }


std::string StochasticFilter::name() const {
	return class_name();
}

std::string StochasticFilter::class_name() {
	return "Stochastic";
}

void StochasticFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "confidence", xsct_real, "XRW TO DO", "1.0" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string StochasticFilterCreator::keyname() const {
	return StochasticFilter::class_name();
}

protocols::filters::FilterOP
StochasticFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new StochasticFilter );
}

void StochasticFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	StochasticFilter::provide_xml_schema( xsd );
}


// XRW TEMP FilterOP
// XRW TEMP CompoundFilterCreator::create_filter() const { return FilterOP( new CompoundFilter ); }

// XRW TEMP std::string
// XRW TEMP CompoundFilterCreator::keyname() const { return "CompoundStatement"; }

std::string CompoundFilter::name() const {
	return class_name();
}

std::string CompoundFilter::class_name() {
	return "CompoundStatement";
}

void CompoundFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "invert", xsct_rosetta_bool, "XRW TO DO", "false" );
	AttributeList subelement_attlist;
	subelement_attlist
		+ XMLSchemaAttribute::required_attribute( "filter_name", xs_string, "XRW TO DO" );
	XMLSchemaSimpleSubelementList subelements;
	subelements
		.add_simple_subelement( "AND", subelement_attlist, "XRW TO DO" )
		.add_simple_subelement( "OR", subelement_attlist, "XRW TO DO" )
		.add_simple_subelement( "XOR", subelement_attlist, "XRW TO DO" )
		.add_simple_subelement( "NOR", subelement_attlist, "XRW TO DO" )
		.add_simple_subelement( "NAND", subelement_attlist, "XRW TO DO" )
		.add_simple_subelement( "ORNOT", subelement_attlist, "XRW TO DO" )
		.add_simple_subelement( "ANDNOT", subelement_attlist, "XRW TO DO" )
		.add_simple_subelement( "NOT", subelement_attlist, "XRW TO DO" );

	protocols::filters::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, class_name(), "XRW TO DO", attlist, subelements );
}

std::string CompoundFilterCreator::keyname() const {
	return CompoundFilter::class_name();
}

protocols::filters::FilterOP
CompoundFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new CompoundFilter );
}

void CompoundFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CompoundFilter::provide_xml_schema( xsd );
}


// XRW TEMP FilterOP
// XRW TEMP CombinedFilterCreator::create_filter() const { return FilterOP( new CombinedFilter ); }

// XRW TEMP std::string
// XRW TEMP CombinedFilterCreator::keyname() const { return "CombinedValue"; }

std::string CombinedFilter::name() const {
	return class_name();
}

std::string CombinedFilter::class_name() {
	return "CombinedValue";
}

void CombinedFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "XRW TO DO", "0.0" );

	XMLSchemaSimpleSubelementList subelements;
	AttributeList subelement_attlist;
	subelement_attlist
		+ XMLSchemaAttribute( "factor", xsct_real, "XRW TO DO" )
		+ XMLSchemaAttribute( "temp", xsct_real, "XRW TO DO" )
		+ XMLSchemaAttribute::required_attribute( "filter_name", xs_string, "XRW TO DO" );

	subelements.add_simple_subelement( "Add", subelement_attlist, "XRW TO DO" );

	protocols::filters::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, class_name(), "XRW TO DO", attlist, subelements );
}

std::string CombinedFilterCreator::keyname() const {
	return CombinedFilter::class_name();
}

protocols::filters::FilterOP
CombinedFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new CombinedFilter );
}

void CombinedFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CombinedFilter::provide_xml_schema( xsd );
}


// XRW TEMP FilterOP
// XRW TEMP MoveBeforeFilterCreator::create_filter() const { return FilterOP( new MoveBeforeFilter ); }

// XRW TEMP std::string
// XRW TEMP MoveBeforeFilterCreator::keyname() const { return "MoveBeforeFilter"; }

std::string MoveBeforeFilter::name() const {
	return class_name();
}

std::string MoveBeforeFilter::class_name() {
	return "MoveBeforeFilter";
}

void MoveBeforeFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute(
		"mover", xs_string,
		"Name of the mover to be applied in advance of filter execution" )
		+ XMLSchemaAttribute(
		"mover_name", xs_string,
		"Name of the mover to be applied in advance of filter execution" )
		+ XMLSchemaAttribute(
		"filter", xs_string,
		"Filter succeeded by the mover" )
		+ XMLSchemaAttribute(
		"filter_name", xs_string,
		"Filter succeeded by the mover" );

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Apply a given mover to the pose before calculating the results from another "
		"filter. Note that, like all filters, MoveBeforeFilter cannot change the input "
		"pose - the results of the submover will only be used for the subfilter "
		"calculation and then discarded.",
		attlist );
}

std::string MoveBeforeFilterCreator::keyname() const {
	return MoveBeforeFilter::class_name();
}

protocols::filters::FilterOP
MoveBeforeFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new MoveBeforeFilter );
}

void MoveBeforeFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MoveBeforeFilter::provide_xml_schema( xsd );
}


// XRW TEMP FilterOP
// XRW TEMP IfThenFilterCreator::create_filter() const { return FilterOP( new IfThenFilter ); }

// XRW TEMP std::string
// XRW TEMP IfThenFilterCreator::keyname() const { return "IfThenFilter"; }

std::string IfThenFilter::name() const {
	return class_name();
}

std::string IfThenFilter::class_name() {
	return "IfThenFilter";
}

void IfThenFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "XRW TO DO", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "lower_threshold", xsct_rosetta_bool, "XRW TO DO", "false" );

	AttributeList if_elif_attlist;
	if_elif_attlist
		+ XMLSchemaAttribute( "valuefilter", xs_string, "XRW TODO" )
		+ XMLSchemaAttribute::attribute_w_default(
		"value", xsct_real,
		"Value of the filter used when returned true",
		"0.0" )
		+ XMLSchemaAttribute::attribute_w_default(
		"weight", xsct_real,
		"The value will be multiplied with weight in case of true",
		"1.0" )
		+ XMLSchemaAttribute::required_attribute(
		"testfilter", xs_string,
		"If testfilter returns true its value will be multiplied by the weight as final result" )
		+ XMLSchemaAttribute::attribute_w_default(
		"inverttest", xsct_rosetta_bool,
		"Inverses the filter result",
		"false" );

	AttributeList else_attlist;
	else_attlist
		+ XMLSchemaAttribute(
		"valuefilter", xs_string,
		"If all all the testfilters failed, the return value is the value of the "
		"valuefilter multiplied by its weight" )
		+ XMLSchemaAttribute::attribute_w_default(
		"value", xsct_real,
		"Value of the valuefilter",
		"0.0" )
		+ XMLSchemaAttribute::attribute_w_default(
		"weight", xsct_real,
		"Weight of the valuefilter",
		"1.0" );

	XMLSchemaSimpleSubelementList subelements;
	subelements
		.add_simple_subelement( "IF", if_elif_attlist,
		"If clause" )
		.add_simple_subelement( "ELIF", if_elif_attlist,
		"else if clause" )
		.add_simple_subelement( "ELSE", else_attlist,
		"else statement" );

	protocols::filters::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, class_name(),
		"Evaluate to a value contingent on the true/false value of other filters. "
		"Each of the IF clauses are evaluated in order. If the testfilter evaluates "
		"to true, the real-valued result of the IfThenFilter is the real-valued return "
		"value of the valuefilter, multiplied by the corresponding weight. (If "
		"inverttest is true, a false testfilter will cause valuefilter evaluation.) "
		"Alternatively, you can omit the valuefilter, and give a literal value with "
		"the value parameter (which will also be multiplied by the given weight). If "
		"none of the IF clauses return true for their testfilters, then the real-"
		"valued result of the ELSE clause valuefilter (or the corresponding literal "
		"value) multiplied by the weight is used as the value instead. For truth value "
		"testing, the default is to return true if the value is less than or equal to "
		"the given threshold. If lower_threshold is true, then IfThenFilter returns "
		"true if the value is greater than or equal to the threshold.",
		attlist, subelements );

}

std::string IfThenFilterCreator::keyname() const {
	return IfThenFilter::class_name();
}

protocols::filters::FilterOP
IfThenFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new IfThenFilter );
}

void IfThenFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	IfThenFilter::provide_xml_schema( xsd );
}



} // filters
} // protocols
