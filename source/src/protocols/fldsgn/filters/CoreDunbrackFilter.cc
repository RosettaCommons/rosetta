// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/CoreDunbrackFilter.cc
/// @brief filter structures by packstat score
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/CoreDunbrackFilter.hh>
#include <protocols/fldsgn/filters/CoreDunbrackFilterCreator.hh>

// Project Headers

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/sasa.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>


// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
//// C++ headers
static THREAD_LOCAL basic::Tracer tr( "protocols.fldsgn.filters.CoreDunbrackFilter" );

namespace protocols {
namespace fldsgn {
namespace filters {

// @brief default constructor
CoreDunbrackFilter::CoreDunbrackFilter():
	Filter( "CoreDunbrack" ),
	filter_value_( 100.0 ),
	type_( "average" ),
	fa_dun_danger_( 1.0 )
{}

// @brief constructor with arguments
CoreDunbrackFilter::CoreDunbrackFilter( String const & type, Real const value ):
	Filter( "CoreDunbrack" ),
	filter_value_( value ),
	type_( type ),
	fa_dun_danger_( 1.0 )
{}

// @brief copy constructor
CoreDunbrackFilter::CoreDunbrackFilter( CoreDunbrackFilter const & rval ):
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	filter_value_( rval.filter_value_ ),
	type_( rval.type_ ),
	fa_dun_danger_( rval.fa_dun_danger_ )
{}

// @brief set filter value
void CoreDunbrackFilter::filter_value( Real const & value )
{
	filter_value_ = value;
}

// @brief set filter type
void CoreDunbrackFilter::filter_type( String const & value )
{
	type_ = value;
}

/// @brief
CoreDunbrackFilter::Real
CoreDunbrackFilter::report_sm( Pose const & pose ) const
{
	return  compute( pose );
}

/// @brief
void
CoreDunbrackFilter::report( std::ostream & out, Pose const & pose ) const
{
	out << "CoreDunbrack: " <<  compute( pose ) << std::endl;
}

/// @brief
CoreDunbrackFilter::Real
CoreDunbrackFilter::compute( Pose const & pose ) const
{

	// define atom_map for main-chain and CB
	core::id::AtomID_Map< bool > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose, false );
	for ( Size ir = 1; ir <= pose.total_residue(); ++ir ) {
		for ( Size j = 1; j<=5; ++j ) {
			core::id::AtomID atom( j, ir );
			atom_map.set( atom, true );
		}
	}

	// calc sasa
	Real pore_radius( 2.0 );
	core::id::AtomID_Map< Real > atom_sasa;
	utility::vector1< Real > rsd_sasa;
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, pore_radius, false, atom_map );

	// store working pose
	Pose copy_pose( pose );
	ScoreFunctionOP sfxn =  core::scoring::get_score_function();
	(*sfxn)( copy_pose );

	Real asa_core( 40.0 );
	Real frustration( 0.0 );
	Size num_frustrated_residue( 0 );
	Size rnum_core( 0 );
	// take dunbrack score
	for ( Size i=1; i<=pose.total_residue(); i++ ) {
		if ( rsd_sasa[ i ] < asa_core ) {
			if ( name_from_aa( pose.aa( i ) ) == "TRP" ||
					name_from_aa( pose.aa( i ) ) == "TYR" ||
					name_from_aa( pose.aa( i ) ) == "MET" ||
					name_from_aa( pose.aa( i ) ) == "PHE" ||
					name_from_aa( pose.aa( i ) ) == "ILE" ||
					name_from_aa( pose.aa( i ) ) == "LEU" ||
					name_from_aa( pose.aa( i ) ) == "VAL" ) {

				Real fa_dun = ( copy_pose.energies().residue_total_energies( i ) )[ core::scoring::fa_dun ];
				if ( fa_dun > fa_dun_danger_ ) {
					tr << "CAUTION high dubrack score " << i << " " << name_from_aa( pose.aa( i ) ) << " " << rsd_sasa[ i ] << " "
						<< fa_dun << std::endl;
					num_frustrated_residue ++;
				}

				frustration += fa_dun;
				rnum_core ++;

			}
		}
	}


	Real ave = Real( frustration )/Real( rnum_core );
	tr << "Frustration: " << frustration << " ,num_hydrophobic_in_core: " << rnum_core
		<< "average: " << ave << " ,num_frustrated_residue: " << num_frustrated_residue << std::endl;

	Real value( 0.0 );
	if ( type_ == "average" ) {
		value = ave;
	} else if ( type_ == "num_frustrated_residue" ) {
		value = num_frustrated_residue;
	} else if ( type_ == "total" ) {
		value = frustration;
	} else {
		tr << "improper type specification " << type_ << std::endl;
		runtime_assert( false );
	}

	return value;
}


// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool CoreDunbrackFilter::apply( Pose const & pose ) const
{
	Real value = compute( pose );
	if ( value < filter_value_ ) {
		tr << "Successfully filtered: " << type_ << " " << value << std::endl;
		return true;
	} else {
		tr << "Filter failed current/threshold=" << value << "/" << filter_value_ << std::endl;
		return false;
	}
} // apply_filter

/// @brief parse xml
void
CoreDunbrackFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & pose )
{
	// set filter type
	type_ = tag->getOption<String>( "type", "average" );

	// set threshold
	filter_value_ = tag->getOption<Real>( "threshold", Real( pose.total_residue() ) );
	if ( type_ == "average" ) {}
	else if ( type_ == "num_frustrated_residue" ) {}
	else if ( type_ == "total" ) {}
	else {
		tr << "invalid type specification, slect from average, num_frustrated_residue, or total" << std::endl;
		runtime_assert( false );
	}

}

protocols::filters::FilterOP
CoreDunbrackFilterCreator::create_filter() const { return protocols::filters::FilterOP( new CoreDunbrackFilter ); }

std::string
CoreDunbrackFilterCreator::keyname() const { return "CoreDunbrack"; }


} // filters
} // fldsgn
} // protocols
