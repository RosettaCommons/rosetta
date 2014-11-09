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
/// @author Doo Nam Kim ( doonam.kim@gmail.com ) responsible for 'h' which is either 'E' or 'L'

// Unit Headers
#include <protocols/fldsgn/filters/SecondaryStructureFilter.hh>
#include <protocols/fldsgn/filters/SecondaryStructureFilterCreator.hh>

// Project Headers
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/util/ABEGOManager.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>
// AUTO-REMOVED #include <protocols/fldsgn/topology/HSSTriplet.hh> // REQUIRED FOR WINDOWS

// Utility headers
#include <basic/Tracer.hh>
#include <protocols/flxbb/utility.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


//// C++ headers
static thread_local basic::Tracer tr( "protocols.fldsgn.filters.SecondaryStructureFilter" );

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
	// count protein residues
	core::Size protein_residues( pose.total_residue() );
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		if ( ! pose.residue( i ).is_protein() ) {
			--protein_residues;
		}
	}

	if( protein_residues != filtered_ss_.length() ) {
		utility_exit_with_message("Length of input ss is not same as total residue of pose.");
	}
	if( filtered_abego_.size() >= 1 &&  protein_residues != filtered_abego_.size() ) {
		utility_exit_with_message("Length of input abego is not same as total residue of pose.");
	}


	return ( compute( pose ) >= protein_residues );
} // apply_filter

/// @brief parse xml
void
SecondaryStructureFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
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

		set_blueprint( blueprint );
		use_abego_ = tag->getOption<bool>( "use_abego", 0 );
	}

	if( filtered_ss_ == "" ){
		tr.Error << "Error!,  option of topology is empty." << std::endl;
		runtime_assert( false );
	}
	tr << filtered_ss_ << " is filtered." << std::endl;

	core::util::ABEGOManager abego_manager;
	if( filtered_abego_.size() >= 1 && use_abego_ ) {
		tr << abego_manager.get_abego_string( filtered_abego_ ) << " is filtered " << std::endl;
	}
}

/// @brief sets the blueprint file based on filename.  If a strand pairing is impossible (i.e. the structure has two strands, 5 and 6 residues, respectively, it sets the unpaired residues to 'h' so that they still match.
void
SecondaryStructureFilter::set_blueprint( std::string const blueprint_file )
{	
	using protocols::jd2::parser::BluePrint;	
	BluePrint blue( blueprint_file );
	std::string ss = blue.secstruct();
	filtered_abego_ = blue.abego();
	
	// now we check the strand pairing definitions in the BluePrint -- if some of them are impossible, DSSP will never return 'E',
	// so we will want these residues to be L/E (which is specified as 'h').
	protocols::fldsgn::topology::SS_Info2_COP ss_info;
 	ss_info = protocols::fldsgn::topology::SS_Info2_COP( protocols::fldsgn::topology::SS_Info2_OP( new protocols::fldsgn::topology::SS_Info2( blue.secstruct() ) ) );
	protocols::fldsgn::topology::StrandPairingSet spairs( blue.strand_pairings(), ss_info );
	protocols::fldsgn::topology::StrandPairings pairs( spairs.strand_pairings() );

	// initialize map that says which residues are paired
	std::map< core::Size, utility::vector1< bool > > paired_residues;
	for ( core::Size strand=1; strand<=pairs.size(); ++strand ) {
		paired_residues[pairs[strand]->s1()] = utility::vector1< bool >( pairs[strand]->size1(), false );
	}

	// determine paired residues
	for ( core::Size strand=1; strand<=pairs.size(); ++strand ) {
		core::Size const s1( pairs[strand]->s1() );
		core::Size const s2( pairs[strand]->s2() );
		int const shift( pairs[strand]->rgstr_shift() );
		core::Size const len1( pairs[strand]->size1() );
		core::Size const len2( pairs[strand]->size2() );
		for ( core::Size res=1; res<=len1; ++res ) {
			core::Size const res2( res-shift );
			if ( res2 >= 1 && res2 <= len2 ) {
				paired_residues[s1][res] = true;
				paired_residues[s2][res2] = true;
			} 
		}
	}	

	// determine unpaired residues and alter them in the string
	for ( core::Size strand=1; strand<=pairs.size(); ++strand ) {
		core::Size const s1( pairs[strand]->s1() );
		core::Size const s2( pairs[strand]->s2() );
		for ( core::Size res=1; res<=paired_residues[s1].size(); ++res ) {
			if ( !paired_residues[s1][res] ) {
				tr << "unpaired " << s1 << " : " << res << std::endl;
				// pairs has begin and end only on paired residues...
				core::Size const pose_res( ss_info->strand(s1)->begin() + res - 1 );
				runtime_assert( pose_res <= ss.size() );
				tr << "pose residue is " << pose_res << std::endl;
				ss[pose_res-1] = 'h'; 
			}
		}
		for ( core::Size res=1; res<=paired_residues[s2].size(); ++res ) {
			if ( !paired_residues[s2][res] ) {
				tr << "unpaired " << s2 << " : " << res << std::endl;
				core::Size const pose_res( ss_info->strand(s2)->begin() + res - 1 );
				runtime_assert( pose_res <= ss.size() );
				tr << "pose residue is " << pose_res << std::endl;
				ss[pose_res-1] = 'h';
			}
		}
	}
	filtered_ss_ = ss;
}

/// returns the number of residues in the protein that have matching secondary structure
/// WARNING: ignores abegos for now, since abegomanager only returns a bool
/// The filter score will report only secondary structures, but if you specify abego,
/// the 'apply' function will return false if abegos don't match
core::Size
SecondaryStructureFilter::compute( core::pose::Pose const & pose ) const {

	core::util::ABEGOManager abego_manager;
	core::Size match_count( 0 );
	for( Size i=1; i<=pose.total_residue(); i++ ){
		// if this residue is a ligand, ignore and move on
		if ( ! pose.residue( i ).is_protein() ) continue;

		String sec = filtered_ss_.substr( i-1, 1 );
		if( *sec.c_str() == 'D' ){
			++match_count;
			continue;
		}else if( *sec.c_str() == 'h' )	{
			// small h means secondary structure other than H, so either L/E, therefore as long as it is not capital H (helix), it passes
			if( pose.secstruct( i ) == 'H')	{
				// it is capital H (helix), so fails
				tr << "SS filter fail: current/filtered = "
				   << pose.secstruct( i ) << '/' << sec << " at position " << i << std::endl;
				continue;
			}
		}else{
			if( *sec.c_str() != pose.secstruct( i ) ){
				tr << "SS filter fail: current/filtered = "
				   << pose.secstruct( i ) << '/' << sec << " at position " << i << std::endl;
				continue;
				//return false;
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
				continue;
			}
		} // if abego
		// to reach this line, we have matched both sec struct and abego, so we can count this residue
		++match_count;
	} //for

	tr << filtered_ss_ << " was filtered with " << match_count << " residues matching " << pose.secstruct() << std::endl;
	if( filtered_abego_.size() >= 1 && use_abego_ ) {
		tr << abego_manager.get_abego_string( filtered_abego_ ) << " was filtered with " << match_count << " residues matching." << std::endl;
	}
	return match_count;
}

core::Real
SecondaryStructureFilter::report_sm( core::pose::Pose const & pose ) const {
	return compute( pose );
}

protocols::filters::FilterOP
SecondaryStructureFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SecondaryStructureFilter ); }

std::string
SecondaryStructureFilterCreator::keyname() const { return "SecondaryStructure"; }


} // filters
} // fldsgn
} // protocols
