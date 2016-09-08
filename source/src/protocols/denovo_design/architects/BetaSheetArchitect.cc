// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/BetaSheetArchitect.cc
/// @brief Architect that creates a beta sheet
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/architects/BetaSheetArchitect.hh>
#include <protocols/denovo_design/architects/BetaSheetArchitectCreator.hh>

// Protocol headers
#include <protocols/denovo_design/architects/StrandArchitect.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/SegmentPairing.hh>
#include <protocols/denovo_design/components/SheetDB.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/denovo_design/types.hh>

// Core headers
#include <core/pose/Pose.hh>

// Basic/Utililty headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.architects.BetaSheetArchitect" );

namespace protocols {
namespace denovo_design {
namespace architects {

BetaSheetArchitect::BetaSheetArchitect( std::string const & id_value ):
	DeNovoArchitect( id_value ),
	permutations_(),
	strands_(),
	extensions_(),
	sheetdb_( new components::SheetDB ),
	use_sheetdb_( false ),
	updated_( false )
{
}

BetaSheetArchitect::~BetaSheetArchitect()
{}

BetaSheetArchitect::DeNovoArchitectOP
BetaSheetArchitect::clone() const
{
	return DeNovoArchitectOP( new BetaSheetArchitect( *this ) );
}

std::string
BetaSheetArchitect::type() const
{
	return BetaSheetArchitect::class_name();
}

void
BetaSheetArchitect::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	std::string const sheet_db_str = tag->getOption< std::string >( "sheet_db", "" );
	if ( !sheet_db_str.empty() ) {
		debug_assert( sheetdb_ );
		sheetdb_->set_db_path( sheet_db_str );
		use_sheetdb_ = true;
	}
	std::string const extensions_str = tag->getOption< std::string >( "strand_extensions", "" );
	if ( tag->hasOption( "strand_extensions" ) ) set_strand_extensions( extensions_str );

	strands_.clear();
	// parse strands, ensure that only strands are present
	utility::tag::Tag::tags_t const & tags = tag->getTags();
	for ( utility::tag::Tag::tags_t::const_iterator subtag=tags.begin(); subtag!=tags.end(); ++subtag ) {
		if ( (*subtag)->getName() == "StrandArchitect" ) {
			StrandArchitect new_strand( "" );
			new_strand.parse_my_tag( *subtag, data );
			TR << "Registered " << new_strand.id() << " as type " << new_strand.type() << std::endl;
			add_strand( new_strand );
		} else {
			utility_exit_with_message( "In sheet " + id() + ", " + (*subtag)->getName() + " is not a strand." );
		}
	}

	if ( !updated_ ) enumerate_permutations();
}

BetaSheetArchitect::StructureDataOP
BetaSheetArchitect::design( core::pose::Pose const &, core::Real & random ) const
{
	if ( !updated_ ) {
		std::stringstream msg;
		msg << "BetaSheetArchitect:apply(): The class list of beta sheet permutations is not updated.  "
			<< "You probably need to call enumerate_permutations() before designing." << std::endl;
		utility_exit_with_message( msg.str() );
	}
	if ( permutations_.empty() ) {
		std::stringstream msg;
		msg << type() << "::apply() No valid Beta sheet permutations could be generated with the given"
			<< " inputs." << std::endl;
		utility_exit_with_message( msg.str() );
	}
	core::Size const idx = extract_int( random, 1, permutations_.size() );
	TR << "Selected permutation " << idx << " of " << permutations_.size() << std::endl;
	return StructureDataOP( new StructureData( *permutations_[ idx ] ) );
}

void
BetaSheetArchitect::set_strand_extensions( std::string const & extensions_str )
{
	typedef utility::vector1< std::string > Strings;
	Strings const extensions = csv_to_container< Strings >( extensions_str, ';' );
	for ( Strings::const_iterator s=extensions.begin(); s!=extensions.end(); ++s ) {
		Strings const name_length = csv_to_container< Strings >( *s, ',' );
		if ( name_length.size() != 2 ) {
			std::stringstream msg;
			msg << "BetaSheetArchitect::set_strand_extensions(): malformed extensions string ("
				<< extensions_str << ") -- the extension " << *s << " has " << name_length.size()
				<< " fields, but two are required (name and length)" << std::endl;
			utility_exit_with_message( msg.str() );
		}
		core::Size const len = boost::lexical_cast< core::Size >( *name_length.rbegin() );
		add_strand_extension( *name_length.begin(), len );
	}
}

void
BetaSheetArchitect::add_strand_extension( std::string const & strand_name, core::Size const length )
{
	extensions_[ strand_name ] = length;
}

void
BetaSheetArchitect::add_strand( StrandArchitect const & strand )
{
	StrandArchitectOP new_strand( new StrandArchitect( strand ) );
	strands_.push_back( new_strand );
	needs_update();
}

/// @brief generates and stores a vector of permutations based on strands
void
BetaSheetArchitect::enumerate_permutations()
{
	TR << "Enumerating permutations!" << std::endl;
	permutations_.clear();
	utility::vector1< components::StructureDataCOPs > perms;
	for ( StrandArchitectOPs::const_iterator s=strands_.begin(); s!=strands_.end(); ++s ) {
		(*s)->enumerate_permutations();
		TR << (*s)->id() << " has " << std::distance( (*s)->motifs_begin(), (*s)->motifs_end() ) << std::endl;
		perms.push_back( components::StructureDataCOPs( (*s)->motifs_begin(), (*s)->motifs_end() ) );
	}
	components::StructureDataCOPs chain;
	combine_permutations_rec( chain, perms );

	updated_ = true;
	TR << "Computed " << permutations_.size() << " permutations!" << std::endl;
}

/// @brief combines the given set of permutations with the current set
void
BetaSheetArchitect::combine_permutations_rec(
	components::StructureDataCOPs const & chain,
	utility::vector1< components::StructureDataCOPs > const & plist )
{
	core::Size const component_idx( chain.size()+1 );

	// stopping condition
	if ( component_idx > plist.size() ) {
		TR.Debug << "adding a chain of length " << chain.size() << std::endl;
		components::StructureDataCOP const p = combine_permutations( chain );
		debug_assert( p );
		// only save if it passes check
		try {
			check_permutation( *p );
			modify_and_add_permutation( *p );
		} catch ( EXCN_PreFilterFailed const & e ) {
			TR.Debug << "Skipping due to failed pre-filter check." << std::endl;
			e.show( TR.Debug );
		}
		return;
	}
	for ( core::Size perm=1; perm<=plist[component_idx].size(); ++perm ) {
		components::StructureDataCOPs new_chain( chain );
		new_chain.push_back( plist[component_idx][perm] );
		combine_permutations_rec( new_chain, plist );
	}
}

/// @brief merges a list of permutations
components::StructureDataOP
BetaSheetArchitect::combine_permutations( components::StructureDataCOPs const & chain ) const
{
	StructureDataOP new_perm( new StructureData( this->id() ) );
	debug_assert( chain.size() == strands_.size() );
	StrandArchitectOPs::const_iterator s=strands_.begin();
	for ( components::StructureDataCOPs::const_iterator m=chain.begin(); m!=chain.end(); ++m, ++s ) {
		new_perm->merge( **m );
	}
	return new_perm;
}

/// @brief modifies/stores data into a permutation and adds it
void
BetaSheetArchitect::modify_and_add_permutation( components::StructureData const & perm )
{
	if ( !use_sheetdb_ ) {
		StructureDataOP toadd( new StructureData( perm ) );
		store_sheet_idx( *toadd, 0 );
		store_strand_pairings( *toadd );
		permutations_.push_back( toadd );
	} else {
		components::SheetList const & list = sheetdb_->sheet_list(
			retrieve_lengths( perm ),
			retrieve_orientations( perm ),
			retrieve_register_shifts( perm ) );

		core::Size const total_size = perm.pose_length();
		core::Size sheet_idx = 1;
		for ( components::SheetList::const_iterator s=list.begin(); s!=list.end(); ++s, ++sheet_idx ) {
			if ( (*s)->size() != total_size ) {
				std::stringstream msg;
				msg << "Sheet from DB with length " << (*s)->size() << " and lengths=" << retrieve_lengths( perm ) << " does not match SD: " << perm << std::endl;
				throw utility::excn::EXCN_BadInput( msg.str() );
			}

			// add template pose for all strands
			StructureDataOP toadd( new StructureData( perm ) );
			core::Size cur_res = 1;
			for ( StrandArchitectOPs::const_iterator strand=strands_.begin(); strand!=strands_.end(); ++strand ) {
				core::Size const start = cur_res + 1;
				core::Size const stopres = start + toadd->segment( (*strand)->id() ).elem_length() - 1;
				TR << "Setting template pose residues " << start << " to " << stopres << std::endl;
				toadd->set_template_pose( (*strand)->id(), **s, start, stopres );
				cur_res += toadd->segment( (*strand)->id() ).length();
			}
			store_sheet_idx( *toadd, sheet_idx );
			store_strand_pairings( *toadd );
			permutations_.push_back( toadd );
		}
	}
}

core::Size
bulge_residue(
	components::StructureData const & sd,
	StrandArchitect const & strand )
{
	std::string const & segment = strand.id();

	StrandBulge bulge = strand.retrieve_bulge( sd );
	if ( bulge < 0 ) {
		bulge = sd.segment( segment ).elem_length() + bulge + sd.segment( segment ).start_local();
	}
	if ( ( bulge < 0 ) || ( bulge > (StrandBulge)sd.segment( segment ).length() ) ) {
		std::stringstream msg;
		msg << "bulge_residue: Error converting " << strand.retrieve_bulge( sd ) << " to a residue id in the segment "
			<< segment << " result=" << bulge << std::endl;
		throw utility::excn::EXCN_BadInput( msg.str() );
	}
	return static_cast< core::Size >( bulge );
}

void
BetaSheetArchitect::store_strand_pairings( StructureData & sd ) const
{
	if ( strands_.size() < 2 ) return;
	for ( StrandArchitectOPs::const_iterator prev=strands_.begin(), c=++strands_.begin(); c!=strands_.end(); ++c, ++prev ) {
		components::StrandPairing pairing(
			(*prev)->id(), (*c)->id(),
			(*prev)->retrieve_orientation( sd ),
			(*c)->retrieve_orientation( sd ),
			(*c)->retrieve_register_shift( sd ) );
		sd.add_pairing( pairing );
	}
}

void
BetaSheetArchitect::store_sheet_idx( StructureData & sd, core::Size const sheet_idx ) const
{
	sd.set_data_int( id(), "sheet_idx", sheet_idx );
}

/// @brief given a permutation, returns Lengths
Lengths
BetaSheetArchitect::retrieve_lengths( StructureData const & perm ) const
{
	Lengths lengths;
	for ( StrandArchitectOPs::const_iterator c=strands_.begin(); c!=strands_.end(); ++c ) {
		lengths.push_back( perm.segment( (*c)->id() ).elem_length() );
	}
	return lengths;
}

/// @brief given a permutation, returns register shifts
RegisterShifts
BetaSheetArchitect::retrieve_register_shifts( StructureData const & perm ) const
{
	RegisterShifts shifts;
	for ( StrandArchitectOPs::const_iterator c=strands_.begin(); c!=strands_.end(); ++c ) {
		shifts.push_back( (*c)->retrieve_register_shift( perm ) );
	}
	return shifts;
}

/// @brief given a permutation, returns Orientations
StrandOrientations
BetaSheetArchitect::retrieve_orientations( StructureData const & perm ) const
{
	StrandOrientations orients;
	for ( StrandArchitectOPs::const_iterator c=strands_.begin(); c!=strands_.end(); ++c ) {
		orients.push_back( (*c)->retrieve_orientation( perm ) );
	}
	return orients;
}

/// @brief looks up and returns extension length for a strand
core::Size
BetaSheetArchitect::extension_length( std::string const & strand ) const
{
	StrandExtensionsMap::const_iterator e = extensions_.find( strand );
	if ( e == extensions_.end() ) return 0;
	else return e->second;
}

/// @brief checks whether the given permutation forms a valid sheet
void
BetaSheetArchitect::check_permutation( components::StructureData const & perm ) const
{
	TR.Debug << "STRANDS: " << std::endl;
	// initialize the residue pairing lookup
	utility::vector1< utility::vector1< bool > > paired_res;
	for ( StrandArchitectOPs::const_iterator c=strands_.begin(); c!=strands_.end(); ++c ) {
		TR.Debug << perm.segment( (*c)->id() ).elem_length() << " "
			<< (*c)->retrieve_orientation( perm ) << " "
			<< (*c)->retrieve_register_shift( perm ) << " "
			<< (*c)->retrieve_bulge( perm ) << std::endl;
		core::Size const strand_length = perm.segment( (*c)->id() ).elem_length() + extension_length( (*c)->id() );
		paired_res.push_back( utility::vector1< bool >( strand_length, false ) );
	}

	utility::vector1< core::Size > bulges;
	for ( StrandArchitectOPs::const_iterator c=strands_.begin(); c!=strands_.end(); ++c ) {
		bulges.push_back( bulge_residue( perm, **c ) );
	}
	debug_assert( bulges.size() == paired_res.size() );

	// bulges will be automatically considered paired
	for ( core::Size strand=1; strand<=strands_.size(); ++strand ) {
		if ( bulges[ strand ] == 0 ) continue;
		paired_res[ strand ][ bulges[ strand ] ] = true;
	}

	// now determine pairing
	for ( core::Size strand=2; strand<=strands_.size(); ++strand ) {
		// for each residue in the first strand, determine what it is paired with
		std::string const & c1_name = strands_[ strand - 1 ]->id();
		std::string const & c2_name = strands_[ strand ]->id();
		int const shift = strands_[ strand ]->retrieve_register_shift( perm );
		core::Size const bulge1 = bulges[ strand  - 1 ];
		core::Size const bulge2 = bulges[ strand ];
		core::Size const length1 = perm.segment( c1_name ).elem_length() + extension_length( c1_name );
		core::Size const length2 = perm.segment( c2_name ).elem_length() + extension_length( c2_name );
		int res2_offset = 0;
		for ( core::Size res = 1; res<=length1; ++res ) {
			core::Size const res2 = res - shift + res2_offset;
			if ( ( res2 == 0 ) || ( res2 > length2 ) ) {
				continue;
			}

			if ( bulge1 == res ) {
				--res2_offset;
				continue;
			}
			if ( bulge2 == res2 ) {
				++res2_offset;
				--res;
				continue;
			}

			TR.Debug << "Strand " << strand-1 << ":" << res << " paired with " << strand << ":" << res2 << std::endl;
			paired_res[ strand - 1 ][ res ] = true;
			paired_res[ strand ][ res2 ] = true;
		}
	}

	// now check the table
	for ( core::Size i = 1; i <= paired_res.size(); ++i ) {
		for ( core::Size j = 1; j <= paired_res[ i ].size(); ++j ) {
			if ( !paired_res[ i ][ j ] ) {
				std::stringstream msg;
				msg << "Returning false: strand " << i << ":" << j << std::endl;
				throw EXCN_PreFilterFailed( msg.str() );
			}
		}
	}

	// now check the database
	Lengths const & lengths_t = retrieve_lengths( perm );
	StrandOrientations const & orientations_t = retrieve_orientations( perm );
	RegisterShifts const & shifts_t = retrieve_register_shifts( perm );
	if ( use_sheetdb_ && ( !sheetdb_->sheet_list( lengths_t, orientations_t, shifts_t ).size() ) ) {
		std::stringstream msg;
		msg << "No sheets were found in the database for lengths=" << lengths_t << " orientations="
			<< orientations_t << " shifts=" << shifts_t << std::endl;
		throw EXCN_PreFilterFailed( msg.str() );
	}
	TR.Debug << "Passing check_permutation" << std::endl;
}

void
BetaSheetArchitect::needs_update()
{
	updated_ = false;
}

///////////////////////////////////////////////////////////////////////////////
/// Creator
///////////////////////////////////////////////////////////////////////////////
std::string
BetaSheetArchitectCreator::keyname() const
{
	return BetaSheetArchitect::class_name();
}

DeNovoArchitectOP
BetaSheetArchitectCreator::create_architect( std::string const & architect_id ) const
{
	return DeNovoArchitectOP( new BetaSheetArchitect( architect_id ) );
}

} //protocols
} //denovo_design
} //architects
