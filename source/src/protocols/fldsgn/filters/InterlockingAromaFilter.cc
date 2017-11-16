// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/InterlockingAromaFilter.cc
/// @brief filter structures by packstat score
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/InterlockingAromaFilter.hh>
#include <protocols/fldsgn/filters/InterlockingAromaFilterCreator.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pose/Pose.hh>
#include <protocols/parser/BluePrint.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

//// C++ headers
static basic::Tracer tr( "protocols.fldsgn.filters.InterlockingAromaFilter" );

namespace protocols {
namespace fldsgn {
namespace filters {

// @brief default constructor
InterlockingAromaFilter::InterlockingAromaFilter():
	Filter( "InterlockingAroma" ),
	filter_value_( 0 ),
	contact_dist2_( numeric::square( 5.5 ) ),
	input_ss_( "" ),
	verbose_( true )
{}


// @brief value constructor
InterlockingAromaFilter::InterlockingAromaFilter( String const & ss ):
	Filter( "InterlockingAroma" ),
	filter_value_( 0 ),
	contact_dist2_( numeric::square( 5.5 ) ),
	input_ss_( ss ),
	verbose_( true )
{}


// @brief copy constructor
InterlockingAromaFilter::InterlockingAromaFilter( InterlockingAromaFilter const & rval ):
	Super( rval ),
	filter_value_( rval.filter_value_ ),
	contact_dist2_( rval.contact_dist2_ ),
	input_ss_( rval.input_ss_ ),
	verbose_( rval.verbose_ )
{}


// @brief set filter value ( defalt 0 )
void
InterlockingAromaFilter::filter_value( Real const value )
{
	filter_value_ = value;
}


// @brief set distance to judge contact ( defalt 4.5 )
void
InterlockingAromaFilter::contact_distance( Real const value )
{
	contact_dist2_ = numeric::square( value );
}


// @brief set verbose on or off ( default on )
void
InterlockingAromaFilter::verbose( bool const b )
{
	verbose_ = b;
}


/// @brief
InterlockingAromaFilter::Real
InterlockingAromaFilter::report_sm( Pose const & pose ) const
{
	return compute( pose );
}


/// @brief
void
InterlockingAromaFilter::report( std::ostream & out, Pose const & pose ) const
{
	out << "InterlockingAroma: " <<  compute( pose ) << std::endl;
}

/// @brief
bool
InterlockingAromaFilter::compute( Size const & res, Pose const & pose, SS_Info2_COP const ssinfo ) const
{
	runtime_assert( pose.aa( res ) == core::chemical::aa_phe ||
		pose.aa( res ) == core::chemical::aa_tyr ||
		pose.aa( res ) == core::chemical::aa_trp );

	using core::Vector;
	using core::conformation::Residue;
	using core::scoring::TenANeighborGraph;

	Vector centroid;
	Residue residue( pose.residue( res ) );
	Size natom( 0 );
	centroid.zero();
	for ( Size iatm=1, iatm_end=residue.natoms(); iatm<=iatm_end; ++iatm ) {
		if ( residue.atom_type( int(iatm) ).atom_type_name() == "aroC" ||
				residue.atom_type( int(iatm) ).atom_type_name() == "Ntrp" ) {
			centroid += residue.atom( iatm ).xyz();
			natom++;
		}
	}
	runtime_assert( natom != 0 );
	centroid = centroid/Real( natom );

	Size contact( 0 );
	utility::vector1< Size > lockres;
	TenANeighborGraph const & energy_graph( pose.energies().tenA_neighbor_graph() );
	for ( utility::graph::Graph::EdgeListConstIter
			iru  = energy_graph.get_node( res )->const_edge_list_begin(),
			irue = energy_graph.get_node( res )->const_edge_list_end(); iru != irue; ++iru ) {

		Size stres( (*iru)->get_second_node_ind() );
		if ( res == stres ) stres = (*iru)->get_first_node_ind();

		if ( ssinfo->strand_id( stres ) == 0 ) continue;

		Vector cb( ssinfo->bb_pos().CB( stres ) );
		Real const dsq( cb.distance_squared( centroid ) );

		if ( dsq <= contact_dist2_ ) {
			contact++;
			lockres.push_back( stres );
		}


		//  if( stres+2 > pose.size() ) continue;
		//  if( ssinfo->strand_id( stres ) == 0 && ssinfo->strand_id( stres+2 ) == 0 ) continue;
		// Vector cb1( ssinfo->bb_pos().CB( stres ) );
		// Real const dsq1( distance_squared( cb1, centroid ));
		//
		// Vector cb2( ssinfo->bb_pos().CB( stres+2 ) );
		// Real const dsq2( distance_squared( cb2, centroid ));
		// tr << res << " " << stres << " " << stres+2 << " " << dsq1 << " " << dsq2 << std::endl;
		//
		// if ( dsq1 <= contact_dist2_ ) {
		//  Vector cb2( ssinfo->bb_pos().CB( stres+2 ) );
		//  Real const dsq2( distance_squared( cb2, centroid ));
		//  if( dsq2 <= contact_dist2_ ) {
		//   if( verbose_ ) {
		//    tr << "Residue " << res << " is interlocking with " << stres << " & " << stres+2 << std::endl;
		//   }
		//   return true;
		//  }
		// }
	}

	if ( contact >= 2 ) {

		tr << "Residue " << res << " is interlocking with ";
		for ( Size ii=1; ii<=lockres.size(); ii++ ) {
			tr << lockres[ ii ] << " ";
		}
		tr << std::endl;
		return true;

	} else {
		return false;
	}

}


/// @brief
InterlockingAromaFilter::Real
InterlockingAromaFilter::compute( Pose const & pose ) const
{
	using core::scoring::dssp::Dssp;
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;

	String ss("");
	if ( ! input_ss_.empty() ) {
		ss = input_ss_;
		runtime_assert( input_ss_.length() == pose.size() );
	} else {
		Dssp dssp( pose );
		ss = dssp.get_dssp_secstruct();
	}
	SS_Info2_OP ssinfo( new SS_Info2( pose, ss ) );

	// calc number of interlocking aromatic residues
	Size num_interlocked( 0 );
	for ( Size ii=1; ii<=pose.size(); ++ii ) {
		if ( ss.at( ii-1 ) == 'E' ) continue;
		if ( pose.aa( ii ) == core::chemical::aa_phe ||
				pose.aa( ii ) == core::chemical::aa_tyr ||
				pose.aa( ii ) == core::chemical::aa_trp ) {
			if ( compute( ii, pose, ssinfo ) ) num_interlocked++;
		}
	}

	return Real(num_interlocked);
}


// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool InterlockingAromaFilter::apply( Pose const & pose ) const
{
	Real value = compute( pose );
	if ( value > filter_value_ ) {
		tr << "Successfully filtered: " << value << std::endl;
		return true;
	} else {
		tr << "Filter failed current/threshold=" << value << "/" << filter_value_ << std::endl;
		return false;
	}
} // apply_filter

/// @brief parse xml
void
InterlockingAromaFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	String const blueprint = tag->getOption<String>( "blueprint", "" );
	if ( blueprint != "" ) {
		parser::BluePrint blue( blueprint );
		input_ss_ = blue.secstruct();
	}
	// set threshold
	filter_value( tag->getOption<Real>( "threshold", 0 ) );

	// set threshold
	contact_distance( tag->getOption<Real>( "dist", 5.5 ) );

	// set threshold
	verbose_ = tag->getOption<bool>( "verbose", 1 );
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP InterlockingAromaFilterCreator::create_filter() const { return protocols::filters::FilterOP( new InterlockingAromaFilter ); }

// XRW TEMP std::string
// XRW TEMP InterlockingAromaFilterCreator::keyname() const { return "InterlockingAroma"; }

std::string InterlockingAromaFilter::name() const {
	return class_name();
}

std::string InterlockingAromaFilter::class_name() {
	return "InterlockingAroma";
}

void InterlockingAromaFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "blueprint", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "XRW TO DO", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "dist", xsct_real, "XRW TO DO", "5.5" )
		+ XMLSchemaAttribute::attribute_w_default( "verbose", xsct_rosetta_bool, "XRW TO DO", "true" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string InterlockingAromaFilterCreator::keyname() const {
	return InterlockingAromaFilter::class_name();
}

protocols::filters::FilterOP
InterlockingAromaFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new InterlockingAromaFilter );
}

void InterlockingAromaFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InterlockingAromaFilter::provide_xml_schema( xsd );
}



} // filters
} // fldsgn
} // protocols
