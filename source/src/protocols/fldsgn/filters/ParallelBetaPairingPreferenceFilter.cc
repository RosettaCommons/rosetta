// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/ParallelBetaPairingPreferenceFilter.cc
/// @brief filter structures by parallel beta pairing preference
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/ParallelBetaPairingPreferenceFilter.hh>
#include <protocols/fldsgn/filters/ParallelBetaPairingPreferenceFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/fldsgn/topology/util.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>

// Utility headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <boost/lexical_cast.hpp>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

//// C++ headers
static basic::Tracer TR( "protocols.fldsgn.filters.ParallelBetaPairingPreferenceFilter" );

namespace protocols {
namespace fldsgn {
namespace filters {

// @brief default constructor
ParallelBetaPairingPreferenceFilter::ParallelBetaPairingPreferenceFilter():
	Filter( "ParallelBetaPairingPreference" ),
	filter_value_( 0.0 ),
	verbose_( false )
{
	using namespace core::chemical;

	utility::io::izstream stream;
	basic::database::open( stream, "parallel_beta_pairing_preference" );

	score_pairmatrix_.resize( 20 );

	String line;
	while ( getline( stream, line ) ) {
		utility::vector1< String > tokens ( utility::string_split( line, '\t' ) );
		runtime_assert( tokens.size() == 21 );
		AA aa = core::chemical::aa_from_oneletter_code( tokens[ 1 ][ 0 ] );
		score_pairmatrix_[ Size ( aa ) ].resize( 20 );
		for ( Size ii=1; ii<=20; ii++ ) {
			Real value;
			if ( aa == core::chemical::aa_pro ) {
				value = 0.01;
			} else {
				value = boost::lexical_cast<Real>( tokens[ ii+1 ] );
			}
			score_pairmatrix_[ Size( aa ) ][ ii ] = -std::log( value );
			// TR << tokens[ 0 ][ 0 ] << " " << ii << " " << tokens[ ii ] << std::endl;
		}
	}
	stream.close();
}


// @brief copy constructor
ParallelBetaPairingPreferenceFilter::ParallelBetaPairingPreferenceFilter( ParallelBetaPairingPreferenceFilter const & rval ):
	Super( rval ),
	filter_value_( rval.filter_value_ ),
	score_pairmatrix_( rval.score_pairmatrix_ ),
	verbose_( rval.verbose_ )
{}


// @brief set filter value ( defalt 0 )
void
ParallelBetaPairingPreferenceFilter::filter_value( Real const value )
{
	filter_value_ = value;
}


/// @brief return filter value
ParallelBetaPairingPreferenceFilter::Real
ParallelBetaPairingPreferenceFilter::report_sm( Pose const & pose ) const
{
	return compute( pose );
}


/// @brief report filter results
void
ParallelBetaPairingPreferenceFilter::report( std::ostream & out, Pose const & pose ) const
{
	out << "ParallelBetaPairingPreference: " << compute( pose ) << "\n";
}


// @brief returns true if the given pose passes the filter, false otherwise.
bool
ParallelBetaPairingPreferenceFilter::apply( Pose const & pose ) const
{
	Real score = compute( pose );

	if ( filter_value_ < score ) {
		return true;
	} else {
		return false;
	}
} // apply


/// @brief refer score of residue pair
ParallelBetaPairingPreferenceFilter::Real
ParallelBetaPairingPreferenceFilter::score_pairmatrix( AA aa1, AA aa2 ) const
{
	return score_pairmatrix_[ Size( aa1 ) ][ Size( aa2 ) ];
}


// @brief compute filter value give a pose
ParallelBetaPairingPreferenceFilter::Real
ParallelBetaPairingPreferenceFilter::compute( Pose const & pose ) const
{
	using core::scoring::dssp::Dssp;
	using protocols::fldsgn::topology::calc_strand_pairing_set;
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;
	using protocols::fldsgn::topology::StrandPairingSet;
	using protocols::fldsgn::topology::StrandPairings;
	using protocols::fldsgn::topology::StrandPairing;

	float dssp_hbond_threshold = -0.5;

	Dssp dssp( pose );
	SS_Info2_OP ssinfo( new SS_Info2( pose, dssp.get_dssp_secstruct() ) );
	StrandPairingSet spairset = calc_strand_pairing_set( pose, ssinfo );

	Size num_pair( 0 );
	Real score( 0.0 );
	for ( StrandPairings::const_iterator it=spairset.begin(), ite=spairset.end(); it != ite; ++it ) {
		StrandPairing spair( **it );

		if ( spair.orient() == 'A' ) continue;

		for ( Size ires=spair.begin1(); ires<=spair.end1(); ires++ ) {

			Size jres( spair.residue_pair( ires ) );
			float score1 = dssp.bb_pair_score( ires, jres );
			float score2( 0.0 );
			if ( jres <= pose.size() ) {
				score2 = dssp.bb_pair_score( ires, jres+1 );
			}

			bool ires_HB( false );
			if ( score1 >= dssp_hbond_threshold && score2 < dssp_hbond_threshold ) {
				ires_HB = true;
			} else if ( score1 >= dssp_hbond_threshold && score2 >= dssp_hbond_threshold ) {
				ires_HB = false;
			} else {
				TR.Warning << ires << ", " << jres << " is making strange hbond. " << std::endl;
				continue;
			}

			Real sc;
			String ptn("");
			if ( ires_HB ) {
				ptn = "HB";
				sc = score_pairmatrix( pose.aa( ires ), pose.aa( jres ) );
			} else {
				ptn = "nHB";
				sc = score_pairmatrix( pose.aa( jres ), pose.aa( ires ) );
			}

			if ( verbose_ ) {
				TR << ires << " " << jres << " "
					<< core::chemical::oneletter_code_from_aa( pose.aa( ires ) ) << " "
					<< core::chemical::oneletter_code_from_aa( pose.aa( jres ) ) << " "
					<< ptn << " " << score2 << " " << sc << std::endl;
			}

			score += sc;
			num_pair ++;

		}
	}

	return score/Real( num_pair );

} // apply_filter

/// @brief parse xml
void
ParallelBetaPairingPreferenceFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	// set threshold
	filter_value( tag->getOption<Real>( "threshold", 0 ) );

	// set threshold
	verbose_ = tag->getOption<bool>( "verbose", 0 );
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP ParallelBetaPairingPreferenceFilterCreator::create_filter() const { return protocols::filters::FilterOP( new ParallelBetaPairingPreferenceFilter ); }

// XRW TEMP std::string
// XRW TEMP ParallelBetaPairingPreferenceFilterCreator::keyname() const { return "ParallelBetaPairingPreference"; }

std::string ParallelBetaPairingPreferenceFilter::name() const {
	return class_name();
}

std::string ParallelBetaPairingPreferenceFilter::class_name() {
	return "ParallelBetaPairingPreference";
}

void ParallelBetaPairingPreferenceFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "XRW TO DO", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "verbose", xsct_rosetta_bool, "XRW TO DO", "false" );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string ParallelBetaPairingPreferenceFilterCreator::keyname() const {
	return ParallelBetaPairingPreferenceFilter::class_name();
}

protocols::filters::FilterOP
ParallelBetaPairingPreferenceFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ParallelBetaPairingPreferenceFilter );
}

void ParallelBetaPairingPreferenceFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ParallelBetaPairingPreferenceFilter::provide_xml_schema( xsd );
}



} // filters
} // fldsgn
} // protocols
