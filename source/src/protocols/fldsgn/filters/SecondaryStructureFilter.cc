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
#include <protocols/denovo_design/components/SegmentPairing.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/residue_selectors/PairedSheetResidueSelector.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
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
	selector_( new core::select::residue_selector::TrueResidueSelector ),
	filtered_ss_( "" ),
	filtered_abego_(),
	strand_pairings_( "" ),
	use_abego_( false ),
	use_dssp_( false ),
	threshold_( 1.0 )
{}


// @brief constructor with arguments
SecondaryStructureFilter::SecondaryStructureFilter( String const & ss ):
	Filter( "SecondaryStructure" ),
	selector_( new core::select::residue_selector::TrueResidueSelector ),
	filtered_ss_( ss ),
	filtered_abego_(),
	strand_pairings_( "" ),
	use_abego_( false ),
	use_dssp_( false ),
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

/// @brief Sets sheet topology string
void SecondaryStructureFilter::set_strand_pairings( String const & s )
{
	strand_pairings_ = s;
}

void SecondaryStructureFilter::set_use_dssp( bool const use_ss )
{
	use_dssp_ = use_ss;
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
	basic::datacache::DataMap & data,
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

	set_use_dssp( tag->getOption< bool >( "use_dssp", use_dssp_ ) );

	threshold_ = tag->getOption< core::Real >( "threshold", threshold_ );

	core::select::residue_selector::ResidueSelectorCOP const selector = core::select::residue_selector::parse_residue_selector( tag, data );
	if ( selector ) set_residue_selector( *selector );
}

void
SecondaryStructureFilter::set_residue_selector( core::select::residue_selector::ResidueSelector const & selector )
{
	selector_ = selector.clone();
}

/// @brief Sets the filtered_ss, filtered_abego, and strand_pairings from the blueprint file given.
void
SecondaryStructureFilter::set_blueprint( std::string const & blueprint_file )
{
	using protocols::jd2::parser::BluePrint;
	BluePrint blue( blueprint_file );
	filtered_ss_ = blue.secstruct();
	filtered_abego_ = blue.abego();
	strand_pairings_ = blue.strand_pairings();
}

/// @brief If a strand pairing is impossible (i.e. the structure has two strands, 5 and 6 residues, respectively, it sets the unpaired residues to 'h' so that they still match.
/// @param[in]     pose      Input pose to be tested.
/// @param[in/out] secstruct Desired secondary structure. It will be modified in-place
/// @details If strand_pairings_ is user-specified, the specific residue pairings will be computed based on it.
///          If strand_pairings_ is empty, residue pairings will be taken from the cached StructureData in
///          the pose, if present.
///          If strand_pairings_ is empty and no cached StructureData was found, this will do nothing.
///          Any unpaired residues found in the sheet topology with 'E' secondary structure will be changed
///          to 'h' (to allow either loop or strand)
void
SecondaryStructureFilter::correct_for_incomplete_strand_pairings(
	core::pose::Pose const & pose,
	std::string & secstruct ) const
{
	using core::select::residue_selector::ResidueSubset;
	using protocols::denovo_design::residue_selectors::PairedSheetResidueSelector;

	PairedSheetResidueSelector selector;
	selector.set_secstruct( secstruct );
	if ( !strand_pairings_.empty() ) selector.set_sheet_topology( strand_pairings_ );

	ResidueSubset const paired = selector.apply( pose );

	if ( secstruct.size() != pose.total_residue() ) {
		std::stringstream msg;
		msg << "SecondaryStructureFilter: The secondary structure length ("
			<< secstruct.size() << ") does not match the size of the pose ("
			<< pose.total_residue() << ")" << " secstruct=" << secstruct << std::endl;
		utility_exit_with_message( msg.str() );
	}

	if ( secstruct.size() != paired.size() ) {
		std::stringstream msg;
		msg << "SecondaryStructureFilter: The secondary structure length ("
			<< secstruct.size() << ") does not match the size of the computed residue subset ("
			<< paired.size() << ")" << " secstruct=" << secstruct << " subset=" << paired << std::endl;
		utility_exit_with_message( msg.str() );
	}

	core::Size resid = 1;
	std::string::iterator ss = secstruct.begin();
	for ( ResidueSubset::const_iterator in_pair=paired.begin(); ( ss!=secstruct.end() ) && ( in_pair!=paired.end() ); ++ss, ++in_pair, ++resid ) {
		if ( *ss != 'E' ) continue;
		if ( *in_pair ) continue;
		tr.Debug << "Found unpaired E residue: " << resid << std::endl;
		*ss = 'h';
	}
}

/// @brief returns the number of residues in the protein that have matching secondary structure
/// @param[in] pose   Pose to be analyzed
/// @param[in] subset Subset of the pose to be analyzed. Only residues that are true in the subset
///                   will be included in the computation.
/// WARNING: ignores abegos for now, since abegomanager only returns a bool
/// The filter score will report only secondary structures, but if you specify abego,
/// the 'apply' function will return false if abegos don't match
/// @details If impossible strand pairings exist (i.e a two-strand pose where one strand has
///          6 residues and the other has 5), either 'E' or 'L' will be allowed at the unpaired
///          position
core::Size
SecondaryStructureFilter::compute(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & subset ) const
{
	using protocols::denovo_design::components::StructureData;
	using protocols::denovo_design::components::StructureDataFactory;

	core::sequence::ABEGOManager abego_manager;
	core::Size match_count( 0 );
	String pose_ss = pose.secstruct();
	if ( use_dssp_ ) {
		core::scoring::dssp::Dssp dssp( pose );
		pose_ss = dssp.get_dssp_secstruct();
	}

	std::string filter_ss = get_filtered_secstruct( pose );
	if ( !strand_pairings_.empty() ) correct_for_incomplete_strand_pairings( pose, filter_ss );

	if ( pose.total_residue() != filter_ss.length() ) {
		utility_exit_with_message("Length of input ss is not same as total residue of pose.");
	}
	if ( filtered_abego_.size() >= 1 && pose.total_residue() != filtered_abego_.size() ) {
		utility_exit_with_message("Length of input abego is not same as total residue of pose.");
	}

	for ( Size i=1; i<=pose.total_residue(); ++i ) {
		// if this residue is a ligand, ignore and move on
		if ( ! pose.residue( i ).is_protein() ) continue;
		if ( ! subset[ i ] ) continue;


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
	// count selected protein residues
	core::Size protein_residues( pose.total_residue() );
	core::select::residue_selector::ResidueSubset const subset = selector_->apply( pose );
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		if ( pose.residue( i ).is_protein() && subset[ i ] ) continue;
		--protein_residues;
	}

	return core::Real( compute( pose, subset ) ) / core::Real( protein_residues );
}

/// @brief Returns the desired secondary structure
/// @details  Rules for selecting this seconary structure:
///           1. If a user-specified filtered_ss_ is set, return that
///           2. If StructureData is cached in the pose, use the secondary structure of that
///           3. Use pose secondary structure, throwing and error if use_dssp is false
std::string
SecondaryStructureFilter::get_filtered_secstruct( core::pose::Pose const & pose ) const
{
	// 1. Use user-specified filtered_ss_ if non-empty
	if ( !filtered_ss_.empty() ) return filtered_ss_;

	// 2. Use StructureData if present
	// TODO: uncomment this once I've updated the 'Tomponent' classes
	//denovo_design::components::StructureDataFactory const & factory =
	//	*denovo_design::components::StructureDataFactory::get_instance();
	//if ( factory.has_cached_data( pose ) ) {
	//	tr.Debug << "Getting filtered_ss from StructureData" << std::endl;
	//	return factory.get_from_const_pose( pose ).ss();
	//}

	// 3. Use pose secstruct
	tr.Debug << "Getting filtered_ss from pose" << std::endl;
	if ( !use_dssp_ ) {
		std::stringstream msg;
		msg << "SecondaryStructureFilter::get_filtered_secstruct(): User has requested to use the "
			<< "secondary structure in the pose instead of DSSP as the test secondary structure. "
			<< "However, the filtered secondary structure is also being determined from the pose. "
			<< "You should either set use_dssp to true, and put the desired secondary structure into the pose, "
			<< "or specify a filtered secondary structure."
			<< std::endl;
		utility_exit_with_message( msg.str() );
	}
	return pose.secstruct();
}

protocols::filters::FilterOP
SecondaryStructureFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SecondaryStructureFilter ); }

std::string
SecondaryStructureFilterCreator::keyname() const { return "SecondaryStructure"; }


} // filters
} // fldsgn
} // protocols
