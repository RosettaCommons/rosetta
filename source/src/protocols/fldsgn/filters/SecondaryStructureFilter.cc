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
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )
/// @author Doo Nam Kim ( doonam.kim@gmail.com ) responsible for 'h' which is either 'E' or 'L'

// Unit Headers
#include <protocols/fldsgn/filters/SecondaryStructureFilter.hh>
#include <protocols/fldsgn/filters/SecondaryStructureFilterCreator.hh>

// Project Headers
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/sequence/ABEGOManager.hh>
#include <core/types.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <protocols/flxbb/utility.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


//// C++ headers
static THREAD_LOCAL basic::Tracer tr( "protocols.fldsgn.filters.SecondaryStructureFilter" );

namespace protocols {
namespace fldsgn {
namespace filters {

// @brief default constructor
SecondaryStructureFilter::SecondaryStructureFilter():
	Filter( "SecondaryStructure" ),
	filtered_ss_( "" ),
	use_abego_( false ),
	use_pose_secstruct_( false ),
	threshold_( 1.0 )
{}


// @brief constructor with arguments
SecondaryStructureFilter::SecondaryStructureFilter( String const & ss ):
	Filter( "SecondaryStructure" ),
	filtered_ss_( ss ),
	use_abego_( false ),
	use_pose_secstruct_( false ),
	threshold_( 1.0 )
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
	for ( Size ii=1; s.length(); ++ii ) {
		filtered_abego_.push_back( s.substr( ii-1, 1 ) );
	}
	use_abego_ = true;
}

void SecondaryStructureFilter::set_use_pose_secstruct( bool const use_ss )
{
	use_pose_secstruct_ = use_ss;
}

// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool SecondaryStructureFilter::apply( Pose const & pose ) const
{
	core::Real const score = report_sm( pose );
	tr.Debug << "Score is " << score << " threshold " << threshold_ << std::endl;
	return ( score >= threshold_ );
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
	for ( Size ii=1; abego.length(); ++ii ) {
		filtered_abego_.push_back( abego.substr( ii-1, 1 ) );
	}
	if ( filtered_abego_.size() >= 1 ) {
		use_abego_ = true;
	}

	// use blueprint as input
	String blueprint = tag->getOption<String>( "blueprint", "" );
	if ( blueprint != "" ) {

		if ( filtered_ss_ != "" || filtered_abego_.size() > 1 ) {
			tr << "ss and abego definition in xml will be ignored, and blueprint info is used. " << std::endl;
			filtered_ss_ = "";
			filtered_abego_.clear();
		}

		set_blueprint( blueprint );
		use_abego_ = tag->getOption<bool>( "use_abego", 0 );
	}

	if ( filtered_ss_ == "" ) {
		tr.Warning << "Warning!,  option of topology is empty. SecondaryStructureFilter will attempt to determine it from StructureData information in the pose at runtime" << std::endl;
	} else {
		tr << filtered_ss_ << " is filtered." << std::endl;
	}

	core::sequence::ABEGOManager abego_manager;
	if ( filtered_abego_.size() >= 1 && use_abego_ ) {
		tr << abego_manager.get_abego_string( filtered_abego_ ) << " is filtered " << std::endl;
	}

	if ( tag->hasOption( "use_pose_secstruct" ) ) {
		set_use_pose_secstruct( tag->getOption< bool >( "use_pose_secstruct" ) );
	}

	threshold_ = tag->getOption< core::Real >( "threshold", threshold_ );
}

/// @brief sets the blueprint file based on filename.  If a strand pairing is impossible (i.e. the structure has two strands, 5 and 6 residues, respectively, it sets the unpaired residues to 'h' so that they still match.
void
SecondaryStructureFilter::set_blueprint( std::string const & blueprint_file )
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
	typedef std::map< core::Size, utility::vector1< bool > > PairedResiduesMap;
	PairedResiduesMap paired_residues;

	for ( topology::StrandPairings::const_iterator pair=pairs.begin(); pair!=pairs.end(); ++pair ) {
		debug_assert( *pair );
		PairedResiduesMap::iterator pr = paired_residues.find( (*pair)->s1() );
		if ( pr == paired_residues.end() ) {
			paired_residues.insert( std::make_pair( (*pair)->s1(), utility::vector1< bool >( (*pair)->size1(), false ) ) );
		} else {
			if ( pr->second.size() < (*pair)->size1() ) {
				pr->second.assign( (*pair)->size1(), false );
			}
		}

		pr = paired_residues.find( (*pair)->s2() );
		if ( pr == paired_residues.end() ) {
			paired_residues.insert( std::make_pair( (*pair)->s2(), utility::vector1< bool >( (*pair)->size2(), false ) ) );
		} else {
			if ( pr->second.size() < (*pair)->size2() ) {
				pr->second.assign( (*pair)->size2(), false );
			}
		}
	}

	for ( PairedResiduesMap::const_iterator pres=paired_residues.begin(); pres!=paired_residues.end(); ++pres ) {
		tr << " Strand is " << pres->first << " pairings " << pres->second << std::endl;
	}

	// determine paired residues
	for ( topology::StrandPairings::const_iterator pair=pairs.begin(); pair!=pairs.end(); ++pair ) {
		for ( core::Size s1_resid=1; s1_resid<=(*pair)->size1(); ++s1_resid ) {
			core::Size const s2_resid = s1_resid - (*pair)->rgstr_shift();
			if ( s2_resid < 1 ) continue;
			if ( s2_resid > (*pair)->size2() ) continue;
			paired_residues[ (*pair)->s1() ][ s1_resid ] = true;
			paired_residues[ (*pair)->s2() ][ s2_resid ] = true;
		}
	}

	// parse the created map looking for unpaired residues
	for ( PairedResiduesMap::const_iterator pres=paired_residues.begin(); pres!=paired_residues.end(); ++pres ) {
		core::Size res = 1;
		for ( utility::vector1< bool >::const_iterator r=pres->second.begin(); r!=pres->second.end(); ++r, ++res ) {
			if ( ! *r ) {
				tr << "unpaired " << pres->first << " : " << res << std::endl;
				// pairs has begin and end only on paired residues...
				core::Size const pose_res( ss_info->strand( pres->first )->begin() + res - 1 );
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

	core::sequence::ABEGOManager abego_manager;
	core::Size match_count( 0 );
	String pose_ss = pose.secstruct();
	if ( !use_pose_secstruct_ ) {
		core::scoring::dssp::Dssp dssp( pose );
		pose_ss = dssp.get_dssp_secstruct();
	}

	std::string filter_ss = filtered_ss_;
	if ( filtered_ss_.empty() ) {
		protocols::denovo_design::components::StructureDataOP sd =
			protocols::denovo_design::components::StructureData::create_from_pose( pose, "SSFilter" );
		filter_ss = sd->ss();
	}

	if ( pose.total_residue() != filter_ss.length() ) {
		utility_exit_with_message("Length of input ss is not same as total residue of pose.");
	}
	if ( filtered_abego_.size() >= 1 && pose.total_residue() != filtered_abego_.size() ) {
		utility_exit_with_message("Length of input abego is not same as total residue of pose.");
	}

	for ( Size i=1; i<=pose.total_residue(); ++i ) {
		// if this residue is a ligand, ignore and move on
		if ( ! pose.residue( i ).is_protein() ) continue;

		char const sec = filter_ss[ i - 1 ];
		if ( sec == 'D' ) {
			++match_count;
			continue;
		} else if ( sec == 'h' ) {
			// small h means secondary structure other than H, so either L/E, therefore as long as it is not capital H (helix), it passes
			if ( pose_ss[ i - 1 ] == 'H' ) {
				// it is capital H (helix), so fails
				tr << "SS filter fail: current/filtered = "
					<< pose_ss[ i - 1 ] << '/' << sec << " at position " << i << std::endl;
				continue;
			}
		} else {
			if ( sec != pose_ss[ i - 1 ] ) {
				tr << "SS filter fail: current/filtered = "
					<< pose_ss[ i - 1 ] << '/' << sec << " at position " << i << std::endl;
				continue;
			}
		}

		if ( ( filtered_abego_.size() >= 1 ) && use_abego_ ) {
			String abego( filtered_abego_[ i ] );
			bool flag( false );
			for ( Size j=1; j<=abego.size(); ++j ) {
				if ( abego_manager.check_rama( *abego.substr( j-1, 1 ).c_str(), pose.phi( i ), pose.psi( i ), pose.omega( i ) ) ) {
					flag = true;
				}
			}
			if ( !flag ) {
				tr << "Abego filter fail current(phi,psi,omega)/filtered = "
					<< pose.phi( i ) << "," << pose.psi( i ) << "," << pose.omega( i ) << "/" << abego
					<< " at position " << i << std::endl;
				continue;
			}
		} // if abego
		// to reach this line, we have matched both sec struct and abego, so we can count this residue
		++match_count;
	} //for

	tr << filter_ss << " was filtered with " << match_count << " residues matching " << pose_ss << std::endl;
	if ( ( filtered_abego_.size() >= 1 ) && use_abego_ ) {
		tr << abego_manager.get_abego_string( filtered_abego_ ) << " was filtered with " << match_count << " residues matching." << std::endl;
	}
	return match_count;
}

core::Real
SecondaryStructureFilter::report_sm( core::pose::Pose const & pose ) const
{
	// count protein residues
	core::Size protein_residues( pose.total_residue() );
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		if ( ! pose.residue( i ).is_protein() ) {
			--protein_residues;
		}
	}

	return core::Real( compute( pose ) ) / core::Real( protein_residues );
}

protocols::filters::FilterOP
SecondaryStructureFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SecondaryStructureFilter ); }

std::string
SecondaryStructureFilterCreator::keyname() const { return "SecondaryStructure"; }


} // filters
} // fldsgn
} // protocols
