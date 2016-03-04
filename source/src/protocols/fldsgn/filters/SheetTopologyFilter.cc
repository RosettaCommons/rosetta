// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/SheetTopologyFilter.cc
/// @brief filter structures by sheet topology
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/SheetTopologyFilter.hh>
#include <protocols/fldsgn/filters/SheetTopologyFilterCreator.hh>

// Package Headers
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/fldsgn/topology/SheetFoldTypeManager.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/util.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


//// C++ headers
static THREAD_LOCAL basic::Tracer tr( "protocols.fldsgn.filters.SheetTopologyFilter" );

namespace protocols {
namespace fldsgn {
namespace filters {

// @brief default constructor
SheetTopologyFilter::SheetTopologyFilter():
	Filter( "SheetTopology" ),
	secstruct_input_( "" ),
	ignore_register_shift_( false ),
	use_dssp_( true ),
	ssinfo_( SS_Info2_OP( new SS_Info2 ) )
{}

// @brief constructor with arguments
SheetTopologyFilter::SheetTopologyFilter( StrandPairingSetOP const & sps ):
	Filter( "SheetTopology" ),
	secstruct_input_( "" ),
	ignore_register_shift_( false ),
	use_dssp_( true ),
	ssinfo_( SS_Info2_OP( new SS_Info2 ) )
{
	filtered_sheet_topology_ = (*sps).name();
}

// @brief constructor with arguments
SheetTopologyFilter::SheetTopologyFilter( String const & sheet_topology ):
	Filter( "SheetTopology" ),
	filtered_sheet_topology_( sheet_topology ),
	secstruct_input_( "" ),
	ignore_register_shift_( false ),
	use_dssp_( true ),
	ssinfo_( SS_Info2_OP( new SS_Info2 ) )
{}

// @brief copy constructor
SheetTopologyFilter::SheetTopologyFilter( SheetTopologyFilter const & rval ):
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	filtered_sheet_topology_( rval.filtered_sheet_topology_ ),
	secstruct_input_( rval.secstruct_input_ ),
	ignore_register_shift_( rval.ignore_register_shift_ ),
	use_dssp_( rval.use_dssp_ ),
	ssinfo_( rval.ssinfo_ )
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
compute_paired_residues( topology::StrandPairingCOP filt_pair, topology::StrandPairingCOP pair )
{
	core::Size const s1_size = filt_pair->end1() - filt_pair->begin1() + 1;
	core::Size const s2_size = filt_pair->end2() - filt_pair->begin2() + 1;
	core::Size startres = filt_pair->begin1();
	core::Size endres = filt_pair->end1();
	if ( s2_size < s1_size ) {
		startres = filt_pair->begin2();
		endres = filt_pair->end2();
	}

	debug_assert( startres <= endres );
	core::Size paircount = 0;
	for ( core::Size res=startres; res<=endres; ++res ) {
		if ( filt_pair->rgstr_shift() == 99 ) {
			++paircount;
		} else {
			if ( pair->has_paired_residue( res ) ) {
				if ( pair->residue_pair( res ) == filt_pair->residue_pair( res ) ) {
					tr.Debug << "Good residue = " << res << " paired to " << filt_pair->residue_pair( res ) << std::endl;
					++paircount;
				} else {
					tr.Debug << "Bad residue = " << res << " paired to " << pair->residue_pair( res ) << " and not " << filt_pair->residue_pair( res ) << std::endl;
				}
			} else {
				tr.Debug << "Bad residue = " << res << " not paired to anything." << std::endl;
			}
		}
	}
	return paircount;
}

core::Size
compute_total_paired_residues( topology::StrandPairingCOP filt_pair )
{
	if ( filt_pair->rgstr_shift() == 99 ) {
		core::Size const s1_size = filt_pair->end1() - filt_pair->begin1() + 1;
		core::Size const s2_size = filt_pair->end2() - filt_pair->begin2() + 1;
		core::Size len = s1_size;
		if ( s2_size < len ) len = s2_size;
		tr.Debug << "paired residues in " << *filt_pair << " = " << len << std::endl;
		return len;
	} else {
		core::Size len = 0;
		for ( core::Size res=filt_pair->begin1(); res<=filt_pair->end1(); ++res ) {
			tr << "Res " << res << " has paired residue? " << filt_pair->has_paired_residue( res ) << std::endl;
			if ( filt_pair->has_paired_residue( res ) ) ++len;
		}
		return len;
	}
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
core::Real
SheetTopologyFilter::compute( Pose const & pose ) const
{
	using core::scoring::dssp::Dssp;
	using protocols::fldsgn::topology::StrandPairings;
	using protocols::fldsgn::topology::NO_STRANDS;

	std::string ss = "";
	if ( secstruct_input_.empty() ) {
		if ( use_dssp_ ) {
			Dssp dssp( pose );
			ss = dssp.get_dssp_secstruct();
		} else {
			ss = pose.secstruct();
		}
	} else {
		for ( core::Size res=1; res<=pose.total_residue(); ++res ) {
			if ( pose.residue( res ).is_protein() ) {
				ss += secstruct_input_[ res - 1 ];
			} else {
				ss += 'L';
			}
		}
	}
	ssinfo_->initialize( pose, ss );

	if ( !( ssinfo_->strands().size() > 0 ) ) {
		tr << "Structure does not include strands." << std::endl;
		return core::Real( 0.0 );
	}

	std::string sheet_topology = filtered_sheet_topology_;
	if ( sheet_topology.empty() ) {
		protocols::denovo_design::components::StructureDataOP sd =
			protocols::denovo_design::components::StructureData::create_from_pose( pose, "sheetfilter" );
		debug_assert( sd );
		sheet_topology = protocols::denovo_design::get_strandpairings( *sd, !ignore_register_shift_ );
		tr << "Topology " << sheet_topology << " will be filtered" << std::endl;
	}

	core::Size const max_strand = compute_max_strand( sheet_topology );
	if ( max_strand > ssinfo_->strands().size() ) {
		tr << "sheet topology contains a strand number (" << max_strand << ") than the structure contains ("
			<< ssinfo_->strands().size() << "). Not a matching topology." << std::endl;
		return core::Real( 0.0 );
	}

	StrandPairingSet spairset_filter( sheet_topology, ssinfo_ );
	tr << "spairset_filter: "<< spairset_filter.name() << std::endl;
	tr.Debug << spairset_filter << std::endl;

	StrandPairingSet spairset = protocols::fldsgn::topology::calc_strand_pairing_set( pose, ssinfo_ );
	tr << "spairset: "<< spairset.name() << std::endl;
	tr.Debug << spairset << std::endl;

	core::Size good_residues = 0;
	core::Size total_residues = 0;
	for ( StrandPairings::const_iterator filt_pair = spairset_filter.begin(); filt_pair != spairset_filter.end(); ++filt_pair ) {
		total_residues += compute_total_paired_residues( *filt_pair );
		for ( StrandPairings::const_iterator pair = spairset.begin(); pair != spairset.end(); ++pair ) {
			// Strand Pairing must match up
			if ( (*pair)->s1() != (*filt_pair)->s1() ) continue;
			if ( (*pair)->s2() != (*filt_pair)->s2() ) continue;

			tr.Debug << "Comparing pairings for strand pairing " << **pair << std::endl;

			// orientation must match or there are no matching residues
			if ( (*pair)->orient() != (*filt_pair)->orient() ) continue;

			good_residues += compute_paired_residues( *filt_pair, *pair );
		}
	}
	tr << "SheetTopology: Good / Total sheet residues = " << good_residues << " / " << total_residues << std::endl;
	return core::Real( good_residues ) / core::Real( total_residues );
	/*
	for ( Size ii = 1; ii <= spairset_filter.size(); ++ii ) {

	bool flag( false );
	for ( Size jj = 1; jj <= spairset.size(); ++jj ) {
	if ( spairset.strand_pairing( jj )->s1() == spairset_filter.strand_pairing( ii )->s1() &&
	spairset.strand_pairing( jj )->s2() == spairset_filter.strand_pairing( ii )->s2() ) {
	//tr << "ii: " << ii  << std::endl;
	//tr << "jj: " << jj  << std::endl;
	//tr << "spairset.strand_pairing( jj )->s1(): " << spairset.strand_pairing( jj )->s1() << std::endl;
	//tr << "spairset_filter.strand_pairing( ii )->s1(): " << spairset_filter.strand_pairing( ii )->s1()  << std::endl;
	//tr << " spairset.strand_pairing( jj )->s2(): " <<  spairset.strand_pairing( jj )->s2() << std::endl;
	//tr << "spairset_filter.strand_pairing( ii )->s2(): " << spairset_filter.strand_pairing( ii )->s2() << std::endl;
	//tr << "spairset.strand_pairing( jj )->orient(): " << spairset.strand_pairing( jj )->orient() << std::endl;
	//tr << "spairset_filter.strand_pairing( ii )->orient(): " << spairset_filter.strand_pairing( ii )->orient() << std::endl;
	//tr << "spairset.strand_pairing( jj )->rgstr_shift(): " << spairset.strand_pairing( jj )->rgstr_shift() << std::endl;
	//tr << "spairset_filter.strand_pairing( ii )->rgstr_shift(): " << spairset_filter.strand_pairing( ii )->rgstr_shift() << std::endl;
	if ( spairset.strand_pairing( jj )->orient() == spairset_filter.strand_pairing( ii )->orient() ) {
	if ( spairset_filter.strand_pairing( ii )->rgstr_shift() != 99 ) {
	if ( spairset.strand_pairing( jj )->rgstr_shift() == spairset_filter.strand_pairing( ii )->rgstr_shift() ) {
	flag = true;
	break;
	}
	} else {
	flag = true;
	break;
	} // register shift ?
	} // orient ?

	if ( !flag ) {
	tr << "Filtering failed - current pair/desired pair :  " << *spairset.strand_pairing( jj )
	<< " / " << *spairset_filter.strand_pairing( ii ) << std::endl;
	}
	} // spairset?
	} // jj
	if ( flag ) {
	good_strands += 1.0;
	}
	} // ii
	good_strands /= core::Real( spairset_filter.size() );
	return good_strands;
	*/
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
	core::Real const good_pairs = compute( pose );

	if ( good_pairs >= 0.9999 ) {
		tr << "Sheet topology " << filtered_sheet_topology_ << " was successfully filtered. " << std::endl;
		return true;
	} else {
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

	// Blueprint is for giving secondary structure information, otherwise dssp will run for ss definition
	// SSPAIR line is read for the topology of strand pairings:w

	String const blueprint = tag->getOption<String>( "blueprint", "" );
	if ( blueprint != "" ) {
		protocols::jd2::parser::BluePrint blue( blueprint );
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
		tr.Warning << "Warning!,  option of topology is empty -- it will be computed from StructureData at runtime" << std::endl;
	} else {
		tr << filtered_sheet_topology_ << " is filtred " << std::endl;
		if ( ignore_register_shift_ ) {
			tr.Warning << "WARNING: ignore_register_shift_ option will be superceded by information from blueprint or topology string." << std::endl;
		}
	}

}

protocols::filters::FilterOP
SheetTopologyFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SheetTopologyFilter ); }

std::string
SheetTopologyFilterCreator::keyname() const { return "SheetTopology"; }

void
SheetTopologyFilter::set_secstruct( std::string const & ss )
{
	secstruct_input_ = ss;
}


} // filters
} // fldsgn
} // protocols
