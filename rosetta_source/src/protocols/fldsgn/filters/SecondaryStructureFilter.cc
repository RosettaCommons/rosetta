// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/SecondaryStructureFilter.cc
/// @brief filter structures by sheet topology
/// @detailed
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/SecondaryStructureFilter.hh>
#include <protocols/fldsgn/filters/SecondaryStructureFilterCreator.hh>

// Project Headers
#include <protocols/fldsgn/BluePrint.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/util/ABEGOManager.hh>
#include <protocols/fldsgn/topology/HSSTriplet.hh> // REQUIRED FOR WINDOWS

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

//Auto Headers


//// C++ headers
static basic::Tracer tr("protocols.fldsgn.filters.SecondaryStructureFilter");

namespace protocols {
namespace fldsgn {
namespace filters {

// @brief default constructor
SecondaryStructureFilter::SecondaryStructureFilter():
	Filter( "SecondaryStructure" ),
	filtered_ss_( "" ),
	use_abego_( false )
{}


// @brief constructor with arguments
SecondaryStructureFilter::SecondaryStructureFilter( String const & ss ):
	Filter( "SecondaryStructure" ),
	filtered_ss_( ss ),
	use_abego_( false )
{}

// @brief copy constructor
SecondaryStructureFilter::SecondaryStructureFilter( SecondaryStructureFilter const & rval ):
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	filtered_ss_( rval.filtered_ss_ ),
	filtered_abego_( rval.filtered_abego_ ),
	use_abego_( rval.use_abego_ )
{}

// @brief set filtered secondary structure
void SecondaryStructureFilter::filtered_ss( String const & s )
{
	filtered_ss_ = s;
}

// @brief set filtered secondary structure
void SecondaryStructureFilter::filtered_abego( String const & s )
{
	filtered_abego_.clear();
	for( Size ii=1; s.length(); ++ii ) {
		filtered_abego_.push_back( s.substr( ii-1, 1 ) );
	}
	use_abego_ = true;
}
	
// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool SecondaryStructureFilter::apply( Pose const & pose ) const
{

	if( pose.total_residue() != filtered_ss_.length() ) {
		utility_exit_with_message("Length of input ss is not same as total residue of pose.");
	}
	if( filtered_abego_.size() >= 1 &&  pose.total_residue() != filtered_abego_.size() ) {
		utility_exit_with_message("Length of input abego is not same as total residue of pose.");
	}

	core::util::ABEGOManager abego_manager;

	for( Size i=1; i<=pose.total_residue(); i++ ){
		String sec = filtered_ss_.substr( i-1, 1 );
		if( *sec.c_str() == 'D' ){
			continue;
		}else{
			if( *sec.c_str() != pose.secstruct( i ) ){
				tr << "SS filter fail: current/filterd = "
				   << pose.secstruct( i ) << '/' << sec << " at position " << i << std::endl;
				return false;
			}
		}

		if( filtered_abego_.size() >= 1 && use_abego_ ) {
			String abego( filtered_abego_[ i ] );
			bool flag( false );
			for( Size j=1; j<=abego.size(); j++ ) {
				if( abego_manager.check_rama( *abego.substr( j-1, 1 ).c_str(), pose.phi( i ), pose.psi( i ), pose.omega( i ) )) {
					flag = true;
				}
			}
			if( !flag ) {
				tr << "Abego filter fail current(phi,psi,omega)/filtered = "
				   << pose.phi( i ) << "," << pose.psi( i ) << "," << pose.omega( i ) << "/" << abego
				   << " at position " << i << std::endl;
				return false;
			}
		}
	}

	tr << "Successfully " << filtered_ss_ << " was filtered. " << std::endl;
	if( filtered_abego_.size() >= 1 && use_abego_ ) {
		tr << "Successfully " << abego_manager.get_abego_string( filtered_abego_ ) << " was filtred. " << std::endl;
	}

	return true;

} // apply_filter

/// @brief parse xml
void
SecondaryStructureFilter::parse_my_tag(
	TagPtr const tag,
	DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	using protocols::fldsgn::BluePrint;

 	filtered_ss_ = tag->getOption<String>( "ss", "" );

	String abego( tag->getOption<String>( "abego", "" ) );
	filtered_abego_.clear();
	for( Size ii=1; abego.length(); ++ii ) {
		filtered_abego_.push_back( abego.substr( ii-1, 1 ) );
	}
	if( filtered_abego_.size() >= 1 ) {
		use_abego_ = true;
	}

	// use blueprint as input
	String blueprint = tag->getOption<String>( "blueprint", "" );
	if( blueprint != "" ) {

		if( filtered_ss_ != "" || filtered_abego_.size() > 1 ) {
			tr << "ss and abego definition in xml will be ignored, and blueprint info is used. " << std::endl;
			filtered_ss_ = "";
			filtered_abego_.clear();
		}

		BluePrint blue( blueprint );
		filtered_ss_ = blue.secstruct();
		filtered_abego_ = blue.abego();
		use_abego_ = tag->getOption<bool>( "use_abego", 0 );
	}

	if( filtered_ss_ == "" ){
		tr.Error << "Error!,  option of topology is empty." << std::endl;
		runtime_assert( false );
	}
	tr << filtered_ss_ << " is filtred." << std::endl;

	core::util::ABEGOManager abego_manager;
	if( filtered_abego_.size() >= 1 && use_abego_ ) {
		tr << abego_manager.get_abego_string( filtered_abego_ ) << " is filtred " << std::endl;
	}
}

protocols::filters::FilterOP
SecondaryStructureFilterCreator::create_filter() const { return new SecondaryStructureFilter; }

std::string
SecondaryStructureFilterCreator::keyname() const { return "SecondaryStructure"; }


} // filters
} // fldsgn
} // protocols
