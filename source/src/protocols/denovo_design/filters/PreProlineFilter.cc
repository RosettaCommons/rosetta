// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/filters/PreProlineFilter.cc
/// @brief Tom's Denovo design protocol
/// @details
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <protocols/denovo_design/filters/PreProlineFilter.hh>
#include <protocols/denovo_design/filters/PreProlineFilterCreator.hh>

// Project Headers
#include <protocols/denovo_design/util.hh>
#include <protocols/rosetta_scripts/util.hh>

// Core Headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/sequence/ABEGOManager.hh>

// Basic Headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

// ObjexxFCL Headers

//C++ Headers

static basic::Tracer TR( "protocols.denovo_design.PreProlineFilter" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace filters {

// XRW TEMP std::string
// XRW TEMP PreProlineFilterCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return PreProlineFilter::class_name();
// XRW TEMP }

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP PreProlineFilterCreator::create_filter() const {
// XRW TEMP  return protocols::filters::FilterOP( new PreProlineFilter() );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP PreProlineFilter::class_name()
// XRW TEMP {
// XRW TEMP  return "PreProline";
// XRW TEMP }

///  ---------------------------------------------------------------------------------
///  PreProlineFilter main code:
///  ---------------------------------------------------------------------------------
PreProlineFilter::PreProlineFilter() :
	Filter( "PreProlineFilter" ),
	threshold_( 0.0 ),
	use_statistical_potential_( false ),
	selector_(),
	spline_()
{
	setup_spline();
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
PreProlineFilter::~PreProlineFilter()
{}


/// Return a copy of ourselves
protocols::filters::FilterOP
PreProlineFilter::clone() const
{
	return protocols::filters::FilterOP( new PreProlineFilter(*this) );
}

protocols::filters::FilterOP
PreProlineFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new PreProlineFilter() );
}

void
PreProlineFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	threshold_ = tag->getOption< core::Real >( "threshold", threshold_ );

	if ( tag->hasOption( "use_statistical_potential" ) ) {
		use_statistical_potential_ = tag->getOption< bool >( "use_statistical_potential" );
	}

	selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, data );
}

std::string
PreProlineFilter::get_name() const
{
	return "PreProline";
}

void
PreProlineFilter::set_selector( core::select::residue_selector::ResidueSelectorCOP selector )
{
	selector_ = selector;
}

void
PreProlineFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out << "PreProlineFilter returning " << compute( pose ) << std::endl;
}

core::Real
PreProlineFilter::report_sm( core::pose::Pose const & pose ) const
{
	return compute( pose );
}

core::Real
PreProlineFilter::compute( core::pose::Pose const & pose ) const
{
	// find selected residues
	utility::vector1< bool > selection( pose.size(), true );
	if ( selector_ ) {
		selection = selector_->apply( pose );
	}

	if ( use_statistical_potential_ ) {
		return compute_spline( pose, selection );
	} else {
		return compute_simple( pose, selection );
	}
}

core::Real
PreProlineFilter::compute_spline(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & selection ) const
{
	core::Size pro_count = 0;
	core::Real potential_sum = 0.0;
	for ( core::Size i = 1; i < pose.size(); ++i ) {
		if ( ! selection[ i + 1 ] ) {
			continue;
		}

		if ( core::pose::is_upper_terminus( pose, i ) ) {
			continue;
		}

		if ( pose.residue( i + 1 ).name1() != 'P' ) {
			continue;
		}

		++pro_count;
		core::Real const splinescore = spline_.F( pose.phi( i ), pose.psi( i ) );
		TR << "Phi = " << pose.phi( i ) << " Psi = " << pose.psi( i )
			<< " spline score = " << splinescore << std::endl;
		potential_sum += splinescore;
	}
	return potential_sum;
}

core::Real
PreProlineFilter::compute_simple(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & selection ) const
{
	// as a simple first trial, only "B" and "E" torsion spaces should be allowed
	core::Size bad_count = core::Size( 0.0 );
	core::Size pro_count = core::Size( 0 );
	std::string const & sequence = pose.sequence();
	std::string const abegos = abego_str( core::sequence::get_abego( pose, 1 ) );
	core::Size resi = 1;
	for ( std::string::const_iterator a = abegos.begin(), s = sequence.begin() + 1;
			( a != abegos.end() ) && ( s != sequence.end() );
			++a, ++s, ++resi ) {
		// ignore if proline is not selected
		if ( ! selection[ resi + 1 ] ) {
			continue;
		}

		// ignore if this is the end of the chain
		if ( core::pose::is_upper_terminus( pose, resi ) ) {
			continue;
		}

		if ( *s == 'P' ) {
			++pro_count;
			TR.Debug << "Res " << resi << " " << pose.residue( resi ).name() << " " << pose.residue( resi + 1 ).name() << " " << *a << std::endl;
			if ( ( *a != 'B' ) && ( *a != 'E' ) ) {
				++bad_count;
			}
		}
	}
	TR << "Prolines in pose: " << pro_count << " Bad pre-proline torsions: " << bad_count << std::endl;
	return static_cast< core::Real >( bad_count );
}

/// @brief Does the PreProline Filtering
bool
PreProlineFilter::apply( core::pose::Pose const & pose ) const
{
	core::Real const value( compute( pose ) );
	if ( threshold_ < 0 || value <= threshold_ ) {
		return true;
	} else {
		return false;
	}
}

numeric::MathMatrix< core::Real >
parse_matrix(
	std::istream & instream,
	core::Size const npoints_x,
	core::Size const npoints_y )
{
	numeric::MathMatrix< core::Real > matrix( npoints_y, npoints_x, core::Real( 0.0 ) );
	for ( core::Size i = 1; i <= npoints_y; ++i ) {
		if ( !instream.good() ) {
			std::stringstream err;
			err << "Error parsing matrix stream before y =" << i << std::endl;
			throw CREATE_EXCEPTION(utility::excn::Exception,  err.str() );
		}
		for ( core::Size ii = 1; ii <= npoints_x; ++ii ) {
			// matrix is 0-indexed
			instream >> matrix[ matrix.get_number_rows() - i  ][ ii - 1 ];
		}
	}
	if ( !instream.good() ) {
		std::stringstream err;
		err << "Error parsing matrix stream in preproline filter!" << std::endl;
		throw CREATE_EXCEPTION(utility::excn::Exception,  err.str() );
	}
	return matrix;
}

void
PreProlineFilter::setup_spline()
{
	// read in training data if necessary
	static std::string const dbfile = "protocol_data/denovo_design/preproline_normalized.gz";
	utility::io::izstream infile( dbfile );
	if ( ! basic::database::open( infile, dbfile ) ) {
		std::stringstream err;
		err << "Pre-proline filter could not open database file " << dbfile << std::endl;
		throw CREATE_EXCEPTION(utility::excn::Exception,  err.str() );
	}

	core::Size npoints_phi, npoints_psi;
	infile >> npoints_phi >> npoints_psi;

	numeric::MathMatrix< core::Real > const data = parse_matrix( infile, npoints_phi, npoints_psi );

	// line after v matrix is start values for each dimension
	core::Real const start[2] = {
		core::Real( -180.0 ),
		core::Real( -180.0 )
		};

	// third line is delta values for each dimension
	core::Real const delta[2] = {
		core::Real( 360.0 ) / core::Real( npoints_phi ),
		core::Real( 360.0 ) / core::Real( npoints_psi )
		};

	TR << "X start: " << start[0] << " delta: " << delta[0] << std::endl;
	TR << "Y start: " << start[1] << " delta: " << delta[1] << std::endl;

	numeric::interpolation::spline::BorderFlag boundary[2] = {
		numeric::interpolation::spline::e_FirstDer,
		numeric::interpolation::spline::e_FirstDer
		};

	std::pair< core::Real, core::Real > firstbe[2] = {
		std::make_pair( 0.0, 0.0 ),
		std::make_pair( 0.0, 0.0 )
		};

	bool lincont[2] = { true, true };

	spline_.train( boundary, start, delta, data, lincont, firstbe );

	TR << "Spline for preproline residues has been trained from " << dbfile << "." << std::endl;
}

void
PreProlineFilter::set_use_statistical_potential( bool const use_stat )
{
	use_statistical_potential_ = use_stat;
}

std::string PreProlineFilter::name() const {
	return class_name();
}

std::string PreProlineFilter::class_name() {
	return "PreProline";
}

void PreProlineFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute(
		"threshold", xsct_real,
		"Returns true if the value of the spline function is less than or equal to the threshold" )
		+ XMLSchemaAttribute(
		"use_statistical_potential", xsct_rosetta_bool,
		"If true, the bicublic spline fit to the statistical potential of "
		"Ramachandran space will be used to evaluate the torsions. If false, "
		"residues in potentially bad torsion bins will be counted. "
		"(default = false)" );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist );

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Helical space (ABEGO type A) is strongly disfavored before a proline, "
		"and Rosetta energy may not capture this effect. In the default mode, "
		"the filter simply checks the ABEGO of all residues before prolines and "
		"counts the residues that have non-B, non-E abego types (lower = better)"
		". With the \"use_statistical_potential\" option, the filter uses a "
		"bicubic spline fit to the data attached to guess at the favorability of "
		"a preproline residue with the given phi/psi angles -- this favorability "
		"is then returned as the filter score (lower = better).",
		attlist );
}

std::string PreProlineFilterCreator::keyname() const {
	return PreProlineFilter::class_name();
}

protocols::filters::FilterOP
PreProlineFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new PreProlineFilter );
}

void PreProlineFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PreProlineFilter::provide_xml_schema( xsd );
}


} // namespace filters
} // namespace denovo_design
} // namespace protocols
