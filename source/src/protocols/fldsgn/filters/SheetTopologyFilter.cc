// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/SheetTopologyFilter.cc
/// @brief filter structures by sheet topology
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/SheetTopologyFilter.hh>
#include <protocols/fldsgn/filters/SheetTopologyFilterCreator.hh>

// Package Headers
#include <protocols/denovo_design/components/SegmentPairing.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/components/SegmentPairing.hh>
#include <protocols/fldsgn/topology/SheetFoldTypeManager.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/util.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/sequence/ABEGOManager.hh>
#include <protocols/parser/BluePrint.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

//// C++ headers


// TEMP
#include <core/pose/PDBInfo.hh>

static THREAD_LOCAL basic::Tracer tr( "protocols.fldsgn.filters.SheetTopologyFilter" );

namespace protocols {
namespace fldsgn {
namespace filters {

// @brief default constructor
SheetTopologyFilter::SheetTopologyFilter():
	Filter( "SheetTopology" ),
	secstruct_input_( "" ),
	ignore_register_shift_( false ),
	use_dssp_( true )
{}

// @brief constructor with arguments
SheetTopologyFilter::SheetTopologyFilter( StrandPairingSetOP const & sps ):
	Filter( "SheetTopology" ),
	secstruct_input_( "" ),
	ignore_register_shift_( false ),
	use_dssp_( true )
{
	filtered_sheet_topology_ = (*sps).name();
}

// @brief constructor with arguments
SheetTopologyFilter::SheetTopologyFilter( String const & sheet_topology ):
	Filter( "SheetTopology" ),
	filtered_sheet_topology_( sheet_topology ),
	secstruct_input_( "" ),
	ignore_register_shift_( false ),
	use_dssp_( true )
{}

// @brief copy constructor
SheetTopologyFilter::SheetTopologyFilter( SheetTopologyFilter const & rval ):
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	filtered_sheet_topology_( rval.filtered_sheet_topology_ ),
	secstruct_input_( rval.secstruct_input_ ),
	ignore_register_shift_( rval.ignore_register_shift_ ),
	use_dssp_( rval.use_dssp_ )
{}

// @brief set filtered sheet_topology by SrandPairingSetOP
void SheetTopologyFilter::filtered_sheet_topology( StrandPairingSetOP const & sps )
{
	filtered_sheet_topology_ = (*sps).name();
}


// @brief set filtered sheet_topology by SrandPairingSetOP
void SheetTopologyFilter::filtered_sheet_topology( String const & sheet_topology )
{
	filtered_sheet_topology_ = sheet_topology;
}

void SheetTopologyFilter::set_use_dssp( bool const use_dssp )
{
	use_dssp_ = use_dssp;
}

core::Size
compute_max_strand( std::string const & sheet_topology )
{
	core::Size max_strand = 0;
	utility::vector1< std::string > const pairs = utility::string_split( sheet_topology, ';' );
	for ( utility::vector1< std::string >::const_iterator p=pairs.begin(); p!=pairs.end(); ++p ) {
		std::string const strands = *( utility::string_split( *p, '.' ).begin() );
		utility::vector1< std::string > const strandlist = utility::string_split( strands, '-' );
		debug_assert( strandlist.size() == 2 );
		for ( utility::vector1< std::string >::const_iterator s=strandlist.begin(); s!=strandlist.end(); ++s ) {
			core::Size const strand = boost::lexical_cast< core::Size >( *s );
			if ( strand > max_strand ) max_strand = strand;
		}
		tr << "strands = " << strands << std::endl;
	}
	return max_strand;
}

/// @brief returns the fraction of pairings that pass the filter
/// @param[in] pose Pose to be checked
/// @details Pose secondary structure is determined by the user inputs, and
///          must match the pose length.
///
///          If the filtered sheet topology doesn't contain and strand pairings,
///          the value returned is 1.0 (i.e. all pairings OK)
///
///          If the pose doesn't contain strands, the value returned is 0.0
///          (i.e. all pairings bad)
///
///          If the pose is missing a strand, the value returned is 0.0 (i.e.
///          all pairings bad)
///
///          Otherwise, NP_actual/NP_filtered is returned, where NP_actual is
///          the number of good residue pairings in the structure, and NP_filtered
///          is total possible residue pairings in the sheet
core::Real
SheetTopologyFilter::compute( Pose const & pose ) const
{
	using core::scoring::dssp::Dssp;
	using protocols::fldsgn::topology::StrandPairings;
	using protocols::fldsgn::topology::NO_STRANDS;

	std::string const ss = get_secstruct( pose );
	if ( ss.size() != pose.size() ) {
		std::stringstream msg;
		msg << "SheetTopologyFilter::compute(): Length of desired secondary structure ("
			<< ss << "; " << ss.size() << ") does not match pose length ("
			<< pose.total_residue() << ")" << std::endl;
		utility_exit_with_message( msg.str() );
	}

	/*
	utility::vector1< std::string > const abego = get_abego( pose );
	if ( abego.size() != pose.total_residue() ) {
	std::stringstream msg;
	msg << "SheetTopologyFilter::compute(): Length of desired abego ("
	<< ss << "; " << abego.size() << ") does not match pose length ("
	<< pose.total_residue() << ")" << std::endl;
	utility_exit_with_message( msg.str() );
	}
	*/

	std::string sheet_topology = get_filtered_sheet_topology( pose );
	if ( sheet_topology.empty() ) {
		tr << "Sheet topology is empty -- all pairings are therefore satisfied" << std::endl;
		return core::Real( 1.0 );
	}

	if ( ignore_register_shift_ ) {
		sheet_topology = remove_register_shifts( sheet_topology );
	}
	tr << "Topology " << sheet_topology << " will be filtered" << std::endl;

	topology::SS_Info2_OP ss_info( new topology::SS_Info2( pose, ss ) );

	// check for missing strands in pose
	if ( !( ss_info->strands().size() > 0 ) ) {
		tr << "Structure does not include strands." << std::endl;
		return core::Real( 0.0 );
	}

	core::Size const max_strand = compute_max_strand( sheet_topology );
	if ( max_strand > ss_info->strands().size() ) {
		tr << "sheet topology contains a strand number (" << max_strand << ") than the structure contains ("
			<< ss_info->strands().size() << "). Not a matching topology." << std::endl;
		return core::Real( 0.0 );
	}

	StrandPairingSet spairset_filter( sheet_topology, ss_info, core::sequence::get_abego( pose ) );
	tr << "spairset_filter: "<< spairset_filter.name() << std::endl;
	tr.Debug << spairset_filter << std::endl;

	StrandPairingSet spairset = protocols::fldsgn::topology::calc_strand_pairing_set( pose, ss_info );
	tr << "spairset: "<< spairset.name() << std::endl;
	tr.Debug << spairset << std::endl;

	replace_register_shifts( spairset, spairset_filter );

	ResiduePairingSets const filtered_residue_pairs = compute_residue_pairings( spairset_filter, *ss_info );
	ResiduePairingSets const pose_residue_pairs = compute_residue_pairings( spairset, *ss_info );

	tr << "Filtered residue pairings: " << filtered_residue_pairs << std::endl;
	tr << "Residue pairings in pose: " << pose_residue_pairs << std::endl;

	// Iterate through associated filtered pairings and computed residue pairing sets
	core::Size good_pairings = 0;
	ResiduePairingSets::const_iterator res_pairset;
	StrandPairings::const_iterator spair;
	for ( spair=spairset_filter.begin(), res_pairset=filtered_residue_pairs.begin();
			( spair!=spairset_filter.end() ) && ( res_pairset!=filtered_residue_pairs.end() );
			++spair, ++res_pairset ) {

		core::Size const pose_pairing_idx = find_pairing_idx( spairset, (*spair)->s1(), (*spair)->s2() );
		if ( pose_pairing_idx == 0 ) {
			tr << "No pose strand pairings found between strands " << (*spair)->s1() << " and "
				<< (*spair)->s2() << std::endl;
			continue;
		}

		topology::StrandPairing const & pose_pairing = *spairset.strand_pairing( pose_pairing_idx );

		if ( pose_pairing.orient() != (*spair)->orient() ) {
			tr << "Orientation for filtered pairing " << **spair << " does not match the one found in the pose ("
				<< pose_pairing << ")" << std::endl;
			continue;
		}

		debug_assert( pose_pairing_idx <= pose_residue_pairs.size() );
		good_pairings += count_good_pairings( *res_pairset, pose_residue_pairs[ pose_pairing_idx ] );
	}

	core::Size const total_pairings = count_residue_pairings( filtered_residue_pairs );

	tr << "SheetTopology: Good / Total sheet residues = " << good_pairings << " / "
		<< total_pairings << std::endl;

	return core::Real( good_pairings ) / core::Real( total_pairings );
}

/// @brief Computes number of pairings in the given StrandPairing
/// @param[in] pairing  StrandPairing which contains residue pairing information
/// @param[in] ss_info  SS_Info2 object describing the secondary structure of the pose
/// @returns ResiduePairingSet containing pairs of residues
SheetTopologyFilter::ResiduePairingSet
SheetTopologyFilter::compute_paired_residues(
	topology::StrandPairing const & pairing,
	topology::SS_Info2 const & ss_info ) const
{
	ResiduePairingSet pairings;
	if ( pairing.rgstr_shift() == 99 ) {
		core::Size const s1_size = ss_info.strand( pairing.s1() )->length();
		core::Size const s2_size = ss_info.strand( pairing.s2() )->length();
		core::Size len = s1_size;
		if ( s2_size < len ) len = s2_size;
		for ( core::Size res=1; res<=len; ++res ) {
			pairings.insert( ResiduePairing( res, res ) );
		}
		tr.Debug << "paired residues in " << pairing << " = " << pairings << std::endl;
		return pairings;
	}

	runtime_assert( pairing.begin1() <= pairing.end1() );
	for ( core::Size res=pairing.begin1(); res<=pairing.end1(); ++res ) {
		core::Size const paired_res = pairing.has_paired_residue( res ) ? pairing.residue_pair( res ) : 0;
		tr.Debug << "Res " << res << " paired residue: " << paired_res;
		if ( paired_res == 0 ) continue;
		pairings.insert( ResiduePairing( res, paired_res ) );
	}
	return pairings;
}

/// @brief Counts total number of residue pairings present in the ResiduePairingSets
core::Size
SheetTopologyFilter::count_residue_pairings( ResiduePairingSets const & pair_sets ) const
{
	core::Size count = 0;
	for ( ResiduePairingSets::const_iterator pset=pair_sets.begin(); pset!=pair_sets.end(); ++pset ) {
		count += pset->size();
	}
	return count;
}

/// @brief Counts number of residue pairs in the filtered_pair_set are present in the pose_pair_set
core::Size
SheetTopologyFilter::count_good_pairings(
	ResiduePairingSet const & filtered_pair_set,
	ResiduePairingSet const & pose_pair_set ) const
{
	core::Size good_pairings = 0;
	for ( ResiduePairingSet::const_iterator filt_pair=filtered_pair_set.begin(); filt_pair!=filtered_pair_set.end(); ++filt_pair ) {
		if ( pose_pair_set.find( *filt_pair ) == pose_pair_set.end() ) {
			tr.Debug << "Filtered pairing " << *filt_pair << " not present in pose pairings." << std::endl;
			continue;
		}

		tr.Debug << "Filtered pairing " << *filt_pair << " found in pose pairings." << std::endl;
		++good_pairings;
	}
	return good_pairings;
}

core::Real
SheetTopologyFilter::report_sm( Pose const & pose ) const
{
	return -compute( pose );
}

// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool SheetTopologyFilter::apply( Pose const & pose ) const
{
	core::Real const good_pair_fraction = compute( pose );

	if ( good_pair_fraction >= 0.9999 ) {
		tr << "Sheet topology " << filtered_sheet_topology_ << " was successfully filtered. " << std::endl;
		return true;
	} else {
		tr << "Sheet topology " << filtered_sheet_topology_ << " is not present. Fraction of pairings present: "
			<< good_pair_fraction << std::endl;
		return false;
	}

} // apply_filter

/// @brief parse xml
void
SheetTopologyFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	if ( tag->hasOption( "ignore_register_shift" ) ) {
		ignore_register_shift_ = tag->getOption< bool >( "ignore_register_shift" );
	}

	filtered_sheet_topology_ = tag->getOption<String>( "topology", "" );
	set_use_dssp( tag->getOption< bool >( "use_dssp", use_dssp_ ) );

	// Blueprint is for giving secondary structure information, otherwise dssp will run for ss definition
	// SSPAIR line is read for the topology of strand pairings

	String const blueprint = tag->getOption<String>( "blueprint", "" );
	if ( blueprint != "" ) {
		protocols::parser::BluePrint blue( blueprint );
		set_secstruct( blue.secstruct() );

		if ( ! blue.strand_pairings().empty() ) {
			if ( filtered_sheet_topology_ == "" ) {
				StrandPairingSet spairset( blue.strand_pairings() );
				filtered_sheet_topology_ = spairset.name();
			} else {
				tr << " SSPAIR line in blueprint will be ignored " << std::endl;
			}
		}
	} //

	if ( filtered_sheet_topology_ == "" ) {
		tr.Warning << "option of topology is empty -- it will be computed from StructureData at runtime" << std::endl;
	} else {
		tr << filtered_sheet_topology_ << " is filtred " << std::endl;
		if ( ignore_register_shift_ ) {
			tr.Warning << "ignore_register_shift_ option will be superceded by information from blueprint or topology string." << std::endl;
		}
	}

}

/// @brief Returns the pose secondary structure to be used in computation
/// @details  Rules for selecting the secondary structure:
///           1. If a user-specified secstruct_input_ is set, return this
///           2. If use_dssp is true, determine secondary structure by DSSP
///           3. Return pose secondary stucture otherwise
std::string
SheetTopologyFilter::get_secstruct( core::pose::Pose const & pose ) const
{
	if ( ! secstruct_input_.empty() ) {
		tr.Debug << "Using user-specified secondary structure: " << secstruct_input_ << std::endl;
		return secstruct_input_;
	}

	if ( use_dssp_ ) {
		core::scoring::dssp::Dssp dssp( pose );
		tr.Debug << "Using DSSP-derived secondary structure: " << dssp.get_dssp_secstruct() << std::endl;
		return dssp.get_dssp_secstruct();
	}

	tr.Debug << "Using pose secondary structure: " << pose.secstruct() << std::endl;
	return pose.secstruct();
}

/// @brief helper function for replacing register shift of a pair with 99
std::string
remove_register_shift_single_pair( std::string const & pair_str )
{
	std::stringstream newstr;
	utility::vector1< std::string > const fields = utility::string_split( pair_str, '.' );
	debug_assert( fields.size() == 3 );
	newstr << fields[1] << '.' << fields[2] << '.' << 99;
	return newstr.str();
}

/// @brief helper function for replacing register shift of all pairs with 99
std::string
remove_register_shifts( std::string const & pair_str )
{
	std::stringstream newstr;
	utility::vector1< std::string > const pairs = utility::string_split( pair_str, ';' );
	for ( utility::vector1< std::string >::const_iterator p=pairs.begin(); p!=pairs.end(); ++p ) {
		if ( !newstr.str().empty() ) newstr << ';';
		newstr << remove_register_shift_single_pair( *p );
	}
	return newstr.str();
}

/// @brief Returns the desired strand pairing topology string
/// @details  Rules for selecting this topology string:
///           1. If a user-specified filtered_sheet_topology_ is set, return that
///           2. If StructureData is cached in the pose, determine pairings from that
///           3. throw error
std::string
SheetTopologyFilter::get_filtered_sheet_topology( core::pose::Pose const & pose ) const
{
	// 1. Use user-set topology string, if specified
	if ( !filtered_sheet_topology_.empty() ) return filtered_sheet_topology_;

	// 2. Get topology from StructureData if present
	denovo_design::components::StructureDataFactory const & factory =
		*denovo_design::components::StructureDataFactory::get_instance();
	if ( factory.has_cached_data( pose ) ) {
		denovo_design::components::StructureData const & sd = factory.get_from_const_pose( pose );
		std::string sheet_topology = protocols::denovo_design::components::SegmentPairing::get_strand_pairings( sd );
		return sheet_topology;
	}

	// Nowhere else to get sheet topology -- input error
	std::stringstream msg;
	msg << "SheetTopologyFilter::get_filtered_sheet_topology(): No sheet topology was specified by the user, and "
		<< "no StructureData object was found in the pose cache. You must specify a topology to filter, "
		<< "or attach StructureData containing desired pairings." << std::endl;
	utility_exit_with_message( msg.str() );
	return "";
}

/// @brief Given the filtered strand pairings, compute the number of residue pairings possible
/// @param[in] spairset The strand pairing set to be used to find residue pairings.  It is
///                     non-const because the pairings are stored as OPs, so begin() and end()
///                     are non-const
/// @param[in] ss_info  SS_Info2 object describing the secondary structure of the pose
/// @returns Vector of ResiduePairingSets, one for each strand pairing, in the same order as
///          in spairset
SheetTopologyFilter::ResiduePairingSets
SheetTopologyFilter::compute_residue_pairings(
	topology::StrandPairingSet & spairset,
	topology::SS_Info2 const & ss_info ) const
{
	using topology::StrandPairings;

	ResiduePairingSets pairing_sets;
	for ( StrandPairings::const_iterator pair=spairset.begin(); pair!=spairset.end(); ++pair ) {
		debug_assert( *pair );
		tr << "Computing residue pairings for " << **pair << std::endl;
		pairing_sets.push_back( compute_paired_residues( **pair, ss_info ) );
	}
	return pairing_sets;
}

/// @brief Replace register shift of pairings in pose_spairset with 99 if register shift in filtered_spairset
///        is 99
void
SheetTopologyFilter::replace_register_shifts(
	topology::StrandPairingSet & spairset,
	topology::StrandPairingSet & filtered_spairset ) const
{
	using topology::StrandPairings;
	topology::StrandPairingSet new_pairset;
	for ( StrandPairings::const_iterator pair=filtered_spairset.begin(); pair!=filtered_spairset.end(); ++pair ) {
		topology::StrandPairingOP pose_pair = find_pairing( spairset, (*pair)->s1(), (*pair)->s2() );
		if ( !pose_pair ) continue;
		if ( (*pair)->rgstr_shift() != 99 ) {
			new_pairset.push_back( pose_pair );
			continue;
		}
		topology::StrandPairingOP new_pair( new topology::StrandPairing(
			pose_pair->s1(), pose_pair->s2(), 99, pose_pair->orient() ) );
		new_pairset.push_back( new_pair );
	}
	spairset = new_pairset;
}

/// @brief Searches the StrandPairingSet for a pairing containing s1 and s2. Returns OP to it
topology::StrandPairingOP
find_pairing( topology::StrandPairingSet & spairset, core::Size const s1, core::Size const s2 )
{
	using topology::StrandPairings;

	for ( StrandPairings::const_iterator p=spairset.begin(); p!=spairset.end(); ++p ) {
		if ( (*p)->s1() != s1 ) continue;
		if ( (*p)->s2() != s2 ) continue;
		return *p;
	}
	return topology::StrandPairingOP();
}

/// @brief Searches the StrandPairingSet for a pairing containing s1 and s2. Returns its 1-based index
core::Size
find_pairing_idx( topology::StrandPairingSet & spairset, core::Size const s1, core::Size const s2 )
{
	using topology::StrandPairings;
	core::Size idx = 1;
	for ( StrandPairings::const_iterator p=spairset.begin(); p!=spairset.end(); ++p, ++idx ) {
		if ( (*p)->s1() != s1 ) continue;
		if ( (*p)->s2() != s2 ) continue;
		return idx;
	}
	return 0;
}

//protocols::filters::FilterOP
//SheetTopologyFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SheetTopologyFilter ); }

void
SheetTopologyFilter::set_secstruct( std::string const & ss )
{
	secstruct_input_ = ss;
}

std::string SheetTopologyFilter::name() const {
	return class_name();
}

std::string SheetTopologyFilter::class_name() {
	return "SheetTopology";
}

void SheetTopologyFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;



	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "ignore_register_shift", xsct_rosetta_bool, "XRW TO DO" )
		+ XMLSchemaAttribute( "topology", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "use_dssp", xsct_rosetta_bool, "XRW TO DO" )
		//If no blueprint is provided, dssp will be used
		+ XMLSchemaAttribute( "blueprint", xs_string, "XRW TO DO" );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string SheetTopologyFilterCreator::keyname() const {
	return SheetTopologyFilter::class_name();
}

protocols::filters::FilterOP
SheetTopologyFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new SheetTopologyFilter );
}

void SheetTopologyFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SheetTopologyFilter::provide_xml_schema( xsd );
}



} // filters
} // fldsgn
} // protocols
