// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ScoreTypeFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)


//Unit Headers
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_filters/ScoreTypeFilterCreator.hh>

//Project Headers
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
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
//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <protocols/elscripts/util.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>
#include <core/pose/util.hh>

namespace protocols {
namespace simple_filters {

using namespace core;
using namespace core::scoring;
using namespace ObjexxFCL::format;

static THREAD_LOCAL basic::Tracer score_type_filter_tracer( "protocols.simple_filters.ScoreTypeFilter" );

protocols::filters::FilterOP
ScoreTypeFilterCreator::create_filter() const { return protocols::filters::FilterOP( new ScoreTypeFilter ); }

std::string
ScoreTypeFilterCreator::keyname() const { return "ScoreType"; }

/// @brief Constructor
///
ScoreTypeFilter::ScoreTypeFilter() :
	Filter( "ScoreType" ),
	score_type_threshold_(0.0),
	score_type_( core::scoring::score_type_from_name( "total_score" ) ),
	scorefxn_() //Null pointer by default.
{}

/// @brief Copy constructor
///
ScoreTypeFilter::ScoreTypeFilter( ScoreTypeFilter const &src ) :
	Filter( "ScoreType" ),
	score_type_threshold_(src.score_type_threshold_),
	score_type_( src.score_type_ ),
	scorefxn_( src.scorefxn_->clone() ) //CLONE the scorefunction, don't copy it.
{}

/// @brief Constructor with parameters
///
ScoreTypeFilter::ScoreTypeFilter( core::scoring::ScoreFunctionCOP scorefxn, core::scoring::ScoreType const score_type, core::Real const score_type_threshold ) : Filter( "ScoreType" ) {
	score_type_threshold_ = score_type_threshold;
	score_type_ = score_type;
	scorefxn_ = scorefxn->clone();
}

ScoreTypeFilter::~ScoreTypeFilter() {}

void
ScoreTypeFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	using namespace core::scoring;

	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();

	score_type_ = core::scoring::score_type_from_name( tag->getOption<std::string>( "score_type", "total_score" ) );
	if ( ! tag->hasOption( "threshold" ) ) throw utility::excn::EXCN_RosettaScriptsOption("Must specify 'threshold' for ScoreTypeFilter.");
	score_type_threshold_ = tag->getOption<core::Real>( "threshold" );

	score_type_filter_tracer<<"ScoreType filter for score_type "<<score_type_<<" with threshold "<<score_type_threshold_<<std::endl;
}
void ScoreTypeFilter::parse_def( utility::lua::LuaObject const & def,
	utility::lua::LuaObject const & score_fxns,
	utility::lua::LuaObject const & /*tasks*/ ) {
	using namespace core::scoring;

	if ( def["scorefxn"] ) {
		scorefxn_ = protocols::elscripts::parse_scoredef( def["scorefxn"], score_fxns );
	} else {
		scorefxn_ = score_fxns["score12"].to<ScoreFunctionSP>()->clone();
	}

	score_type_ = core::scoring::score_type_from_name( def["score_type"] ? def["score_type"].to<std::string>() : "total_score" );
	if ( ! def["threshold"] ) utility_exit_with_message("Must specify 'threshold' for ScoreTypeFilter.");
	score_type_threshold_ = def["threshold"].to<core::Real>();

	score_type_filter_tracer<<"ScoreType filter for score_type "<<score_type_<<" with threshold "<<score_type_threshold_<<std::endl;
}

bool
ScoreTypeFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const score( compute( pose ) );
	score_type_filter_tracer<<"score "<<core::scoring::ScoreTypeManager::name_from_score_type( score_type_ )<<" is "<<score<<". ";
	if ( score <= score_type_threshold_ ) {
		score_type_filter_tracer<<"passing." << std::endl;
		return true;
	} else {
		score_type_filter_tracer<<"failing."<<std::endl;
		return false;
	}
}

void
ScoreTypeFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	out<<"Weighted score of "<<core::scoring::ScoreTypeManager::name_from_score_type( score_type_ )<<" "<<compute( pose )<<'\n';
}

core::Real
ScoreTypeFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
ScoreTypeFilter::compute( core::pose::Pose const & pose ) const {
	using namespace core::pose;
	using namespace core::scoring;

	PoseOP in_pose( new Pose( pose ) );

	// make sure that scoring weights are compatible with pose's residue type set
	// check centroid case
	if ( ( (*scorefxn_)[fa_rep] == 0.0 && (*scorefxn_)[fa_atr] == 0.0 ) // full atom terms are off
			&& ( (*scorefxn_)[interchain_vdw] > 0.0 || (*scorefxn_)[vdw] > 0.0)  ) { // a centroid term is on
		if ( in_pose->is_fullatom() ) { // but pose is full atom
			core::util::switch_to_residue_type_set( *in_pose, core::chemical::CENTROID );
		}
	} else { // full atom case
		if ( in_pose->is_centroid() ) { // but pose is centroid
			core::util::switch_to_residue_type_set( *in_pose, core::chemical::FA_STANDARD );
		}
	}

	(*scorefxn_)( *in_pose );
	/// Now handled automatically.  scorefxn_->accumulate_residue_total_energies( *in_pose );
	core::Real const weight( (*scorefxn_)[ ScoreType( score_type_ ) ] );
	core::Real const score( in_pose->energies().total_energies()[ ScoreType( score_type_ ) ]);
	if ( score_type_ == total_score ) return( score );
	core::Real const weighted_score( weight * score );
	return( weighted_score );
}

}
}
