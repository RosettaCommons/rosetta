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
#include <protocols/fldsgn/topology/SheetFoldTypeManager.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/util.hh>

// Project Headers
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
static thread_local basic::Tracer tr( "protocols.fldsgn.filters.SheetTopologyFilter" );

namespace protocols {
namespace fldsgn {
namespace filters {

// @brief default constructor
SheetTopologyFilter::SheetTopologyFilter():
	Filter( "SheetTopology" ),
	secstruct_input_( false ),
	ssinfo_( SS_Info2_OP( new SS_Info2 ) )
{}

// @brief constructor with arguments
SheetTopologyFilter::SheetTopologyFilter( StrandPairingSetOP const & sps ):
	Filter( "SheetTopology" ),
	secstruct_input_( false ),
	ssinfo_( SS_Info2_OP( new SS_Info2 ) )
{
	filtered_sheet_topology_ = (*sps).name();
}

// @brief constructor with arguments
SheetTopologyFilter::SheetTopologyFilter( String const & sheet_topology ):
	Filter( "SheetTopology" ),
	filtered_sheet_topology_( sheet_topology ),
	secstruct_input_( false ),
	ssinfo_( SS_Info2_OP( new SS_Info2 ) )
{}

// @brief copy constructor
SheetTopologyFilter::SheetTopologyFilter( SheetTopologyFilter const & rval ):
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	filtered_sheet_topology_( rval.filtered_sheet_topology_ ),
	secstruct_input_( rval.secstruct_input_ ),
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

// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool SheetTopologyFilter::apply( Pose const & pose ) const
{
	using core::scoring::dssp::Dssp;
	using protocols::fldsgn::topology::NO_STRANDS;

	if( ! secstruct_input_ ) {
		Dssp dssp( pose );
		ssinfo_->initialize( pose, dssp.get_dssp_secstruct() );
	}

	if( !( ssinfo_->strands().size() > 0 ) ){
		tr << "Structure does not include strands." << std::endl;
		return false;
	}

	//tr << "ss of input pose=" << pose.secstruct() << std::endl;
	//tr << "ssinfo_ " << *ssinfo_ << std::endl;

	StrandPairingSet spairset = protocols::fldsgn::topology::calc_strand_pairing_set( pose, ssinfo_ );
	tr << "spairset: "<< spairset.name() << std::endl;
	StrandPairingSet spairset_filter( filtered_sheet_topology_ );
	tr << "spairset_filter: "<< spairset_filter.name() << std::endl;

	for( Size ii=1; ii<=spairset_filter.size(); ++ii ) {

		bool flag( false );
		for( Size jj=1; jj<=spairset.size(); jj++ ) {
			if( spairset.strand_pairing( jj )->s1() == spairset_filter.strand_pairing( ii )->s1() &&
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
				if( spairset.strand_pairing( jj )->orient() == spairset_filter.strand_pairing( ii )->orient() ) {
					if( spairset_filter.strand_pairing( ii )->rgstr_shift() != 99 ) {
						if( spairset.strand_pairing( jj )->rgstr_shift() == spairset_filter.strand_pairing( ii )->rgstr_shift() ) {
							flag = true;
							break;
						}
					} else {
						flag = true;
						break;
					} // register shift ?
				} else {
					break;					
				} // orient ?
			} // spairset 
		} // jj

		if( !flag ) {
			tr << "Filtering failed: current/filtered sheet topology, " << spairset.name() << '/' << filtered_sheet_topology_ << '.' << std::endl;
			return false;
		}
	} // ii


	tr << "Successfully " << spairset.name() << " sheet topology was filtered. " << std::endl;
	return true;

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
 	filtered_sheet_topology_ = tag->getOption<String>( "topology", "" );

	// Blueprint is for giving secondary structure information, otherwise dssp will run for ss definition
	// SSPAIR line is read for the topology of strand pairings
	String const blueprint = tag->getOption<String>( "blueprint", "" );
	if( blueprint != "" ) {
		protocols::jd2::parser::BluePrint blue( blueprint );
		ssinfo_->initialize( blue.secstruct() );
		secstruct_input_ = true;

		if( ! blue.strand_pairings().empty() ) {
			if( filtered_sheet_topology_ == "" ) {
				StrandPairingSet spairset( blue.strand_pairings() );
				filtered_sheet_topology_ = spairset.name();
			} else {
				tr << " SSPAIR line in blueprint will be igonared " << std::endl;
			}
		}
 	} //

	if( filtered_sheet_topology_ == "" ) {
		tr.Error << "Error!,  option of topology is empty." << std::endl;
		runtime_assert( false );
	}

	tr << filtered_sheet_topology_ << " is filtred " << std::endl;
}

protocols::filters::FilterOP
SheetTopologyFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SheetTopologyFilter ); }

std::string
SheetTopologyFilterCreator::keyname() const { return "SheetTopology"; }


} // filters
} // fldsgn
} // protocols
