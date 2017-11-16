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
#include <protocols/denovo_design/architects/DeNovoArchitectFactory.hh>
// Protocol headers
#include <protocols/denovo_design/architects/StrandArchitect.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/SegmentPairing.hh>
#include <protocols/denovo_design/components/SheetDB.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/denovo_design/types.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueVector.hh>

// Basic/Utililty headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
// Boost headers
#include <boost/assign.hpp>

static basic::Tracer TR( "protocols.denovo_design.architects.BetaSheetArchitect" );

namespace protocols {
namespace denovo_design {
namespace architects {

BetaSheetArchitect::BetaSheetArchitect( std::string const & id_value ):
	DeNovoArchitect( id_value ),
	permutations_(),
	strands_(),
	orientations_(),
	shifts_(),
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
BetaSheetArchitectCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	BetaSheetArchitect::provide_xml_schema( xsd );
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
	orientations_.clear();
	shifts_.clear();

	// parse strands, ensure that only strands are present
	for ( utility::tag::TagCOP const & subtag : tag->getTags() ) {
		if ( subtag->getName() == "StrandArchitect" ) {
			StrandArchitect new_strand( "" );
			new_strand.parse_my_tag( subtag, data );
			add_strand( new_strand );
			TR.Debug << "Registered " << new_strand.id() << " as type " << new_strand.type() << std::endl;

			// get orientation/register shifts
			add_orientations( subtag->getOption< std::string >( "orientation", "" ) );
			add_register_shifts( subtag->getOption< std::string >( "register_shift", "" ) );
		} else {
			utility_exit_with_message( "In sheet " + id() + ", " + subtag->getName() + " is not a strand." );
		}
	}

	if ( !updated_ ) enumerate_permutations();
}

void
BetaSheetArchitect::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "sheet_db", xs_string, "Path to detabase of beta sheets" )
		+ XMLSchemaAttribute( "strand_extensions", xs_string, "Semicolon-separated string specifying strand extensions. Each extension should be a comma-separated pair of name and length" );

	//The subelements are StrandArchitects
	StrandArchitect::provide_xml_schema( xsd );
	XMLSchemaSimpleSubelementList subelements;
	subelements
		.add_already_defined_subelement( "StrandArchitect", & DeNovoArchitectFactory::complex_type_name_for_architect );
	//StrandArchitect

	DeNovoArchitect::add_common_denovo_architect_attributes( attlist );

	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & DeNovoArchitectFactory::complex_type_name_for_architect )
		.element_name( class_name() )
		.description( "Architect used to construct beta sheets" )
		.add_attributes( attlist )
		.set_subelements_repeatable( subelements )
		.write_complex_type_to_schema( xsd );
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
	for ( std::string const & s : extensions ) {
		Strings const name_length = csv_to_container< Strings >( s, ',' );
		if ( name_length.size() != 2 ) {
			std::stringstream msg;
			msg << "BetaSheetArchitect::set_strand_extensions(): malformed extensions string ("
				<< extensions_str << ") -- the extension " << s << " has " << name_length.size()
				<< " fields, but two are required (name and length)" << std::endl;
			utility_exit_with_message( msg.str() );
		}
		core::Size const len = boost::lexical_cast< core::Size >( *name_length.rbegin() );
		add_strand_extension( *name_length.begin(), len );
	}
}

/// @brief set allowed register shifts from a string
void
BetaSheetArchitect::add_register_shifts( std::string const & val )
{
	if ( val.empty() ) {
		shifts_.push_back( RegisterShifts() );
		return;
	}

	RegisterShifts retval;
	utility::vector1< std::string > const str_shifts( utility::string_split( val, ',' ) );
	for ( utility::vector1< std::string >::const_iterator s=str_shifts.begin(); s!=str_shifts.end(); ++s ) {
		TR.Debug << *s << " " << val << std::endl;
		if ( s->empty() ) continue;
		utility::vector1< std::string > const ranges( utility::string_split( *s, ':' ) );
		if ( ranges.size() == 1 ) {
			retval.push_back( boost::lexical_cast< RegisterShift >( ranges[1] ) );
		} else if ( ranges.size() == 2 ) {
			RegisterShift const start( boost::lexical_cast< RegisterShift >( ranges[1] ) );
			RegisterShift const end( boost::lexical_cast< RegisterShift >( ranges[2] ) );
			for ( RegisterShift i=start; i<=end; ++i ) {
				retval.push_back( i );
			}
		} else {
			utility_exit_with_message( "Invalid register shift input: " + val );
		}
	}
	shifts_.push_back( retval );
}

void
BetaSheetArchitect::add_orientations( std::string const & orientations_str )
{
	if ( orientations_str.empty() ) {
		TR.Debug << "Adding blank orientations" << std::endl;
		orientations_.push_back( StrandOrientations() );
		return;
	}

	StrandOrientations retval;
	utility::vector1< std::string > const str_orients( utility::string_split( orientations_str, ',' ) );
	for ( utility::vector1< std::string >::const_iterator s=str_orients.begin(); s!=str_orients.end(); ++s ) {
		TR.Debug << *s << " " << orientations_str << std::endl;
		if ( s->empty() ) continue;
		utility::vector1< std::string > const ranges( utility::string_split( *s, ':' ) );
		// check input string to make sure only "A", "P" are specified
		for ( core::Size i=1; i<=ranges.size(); ++i ) {
			if ( ( ranges[i] != "U" ) && ( ranges[i] != "D" ) ) {
				utility_exit_with_message( "Invalid orientation character: " + ranges[i] );
			}
		}
		if ( ranges.size() == 1 ) {
			if ( ranges[1] == "U" ) {
				retval.push_back( components::UP );
			} else {
				retval.push_back( components::DOWN );
			}
		} else if ( ranges.size() == 2 ) {
			retval.push_back( components::UP );
			retval.push_back( components::DOWN );
		} else {
			utility_exit_with_message( "Invalid orientation input: " + orientations_str );
		}
	}
	TR.Debug << "Adding orientations: " << retval << std::endl;
	orientations_.push_back( retval );
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
	TR.Debug << "Enumerating permutations!" << std::endl;
	permutations_.clear();
	utility::vector1< components::StructureDataCOPs > perms;

	debug_assert( strands_.size() == orientations_.size() );
	debug_assert( strands_.size() == shifts_.size() );

	for ( StrandArchitectOPs::const_iterator s=strands_.begin(); s!=strands_.end(); ++s  ) {
		(*s)->enumerate_permutations();
		TR.Debug << (*s)->id() << " has " << std::distance( (*s)->motifs_begin(), (*s)->motifs_end() ) << std::endl;
		perms.push_back( components::StructureDataCOPs( (*s)->motifs_begin(), (*s)->motifs_end() ) );
	}

	components::StructureDataCOPs chain;
	combine_permutations_rec( chain, perms );
	TR.Debug << "After combining perms, got " << permutations_.size() << " total." << std::endl;
	TR.Debug << "Orientations: " << orientations_ << std::endl;
	TR.Debug << "Shifts      : " << shifts_ << std::endl;

	// Add pairings
	components::StructureDataCOPs perms2 = add_pairings( permutations_ );
	perms2 = filter_permutations( perms2 );

	permutations_.clear();
	for ( components::StructureDataCOP const & p : perms2 ) {
		modify_and_add_permutation( *p );
	}

	updated_ = true;
	TR << "Computed " << permutations_.size() << " permutations!" << std::endl;
}

BetaSheetArchitect::PairingsInfoVector
fill_orientation_info( BetaSheetArchitect::PairingsInfoVector const & by_strand )
{
	if ( by_strand.empty() ) return BetaSheetArchitect::PairingsInfoVector();

	// compute list for remaining strands
	BetaSheetArchitect::PairingsInfoVector::const_iterator remaining_it = by_strand.begin();
	++remaining_it;
	BetaSheetArchitect::PairingsInfoVector const remaining( remaining_it, by_strand.end() );
	BetaSheetArchitect::PairingsInfoVector const combos = fill_orientation_info( remaining );

	if ( combos.empty() ) {
		BetaSheetArchitect::PairingsInfoVector retval;
		for ( components::StrandPairingCOP const & pinfo : *by_strand.begin() ) {
			retval.push_back( boost::assign::list_of (pinfo) );
		}
		return retval;
	}

	// Tack on the first-strand possibilities
	BetaSheetArchitect::PairingsInfoVector retval;
	for ( components::StrandPairingCOP const & pinfo : *by_strand.begin() ) {
		for ( BetaSheetArchitect::PairingsInfo const & prev : combos ) {
			BetaSheetArchitect::PairingsInfo newval = boost::assign::list_of (pinfo);
			for ( components::StrandPairingCOP const & pinfo2 : prev ) {
				newval.push_back( pinfo2 );
			}
			retval.push_back( newval );
		}
	}
	return retval;
}

components::StructureDataCOPs
BetaSheetArchitect::add_pairings( components::StructureDataCOPs const & perms ) const
{
	components::StructureDataCOPs perms2;

	// compute strand-indexed combinations of orientations/shifts
	PairingsInfoVector pinfo_by_strand;
	for ( core::Size sidx=2; sidx<=strands_.size(); ++sidx ) {
		PairingsInfo pinfo;
		for ( StrandOrientation const & o1 : orientations_[ sidx - 1 ] ) {
			for ( StrandOrientation const & o2 : orientations_[ sidx ] ) {
				for ( RegisterShift const s : shifts_[ sidx ] ) {
					pinfo.push_back( components::StrandPairingCOP( new components::StrandPairing(
						strands_[ sidx - 1 ]->id(), strands_[ sidx ]->id(),
						o1, o2, s ) ) );
				}
			}
		}
		pinfo_by_strand.push_back( pinfo );
	}

	PairingsInfoVector const pinfo_perms = fill_orientation_info( pinfo_by_strand );

	for ( components::StructureDataCOP const & perm : perms ) {
		for ( PairingsInfo const & pinfo : pinfo_perms ) {
			components::StructureDataOP newsd( new components::StructureData(*perm) );
			for ( components::StrandPairingCOP const & pair : pinfo ) {
				newsd->add_pairing( *pair );
			}
			perms2.push_back( newsd );
		}
	}
	TR.Debug << "Found " << perms2.size() << " possible permutations (before filtering invalid ones)!" << std::endl;
	return perms2;
}


components::StructureDataCOPs
BetaSheetArchitect::filter_permutations( components::StructureDataCOPs const & perms ) const
{
	components::StructureDataCOPs good;
	for ( components::StructureDataCOP const & p : perms ) {
		try {
			check_permutation( *p );
			good.push_back( p );
		} catch ( EXCN_PreFilterFailed const & e ) {
			TR.Debug << "Permutation failed check -- " << e << std::endl;
		}
	}
	return good;
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
		// add all combinations to list -- they will be filtered later
		permutations_.push_back( p );
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
			permutations_.push_back( toadd );
		}
	}
}

BetaSheetArchitect::ResidueVector
bulge_residues(
	components::StructureData const & sd,
	StrandArchitect const & strand )
{
	BetaSheetArchitect::ResidueVector resids;
	StrandBulges bulges = strand.retrieve_bulges( sd );
	for ( SegmentResid const & b : bulges ) {
		if ( b != 0 ) {
			core::Size const pose_resid = sd.segment( strand.id() ).segment_to_pose( b );
			resids.push_back( static_cast< core::Size >( pose_resid ) );
		}
	}
	return resids;
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
BetaSheetArchitect::RegisterShifts
BetaSheetArchitect::retrieve_register_shifts( StructureData const & perm ) const
{
	RegisterShifts shifts;
	StrandArchitectOPs::const_iterator prev = strands_.end();
	for ( StrandArchitectOPs::const_iterator c=strands_.begin(); c!=strands_.end(); ++c ) {
		if ( c == strands_.begin() ) {
			shifts.push_back( 0 );
		} else {
			components::SegmentPairingCOP pair = perm.pairing( boost::assign::list_of ((*prev)->id())((*c)->id()) );
			if ( !pair ) {
				std::stringstream msg;
				msg << "Pairing not found!!!" << (*prev)->id() << " <--> " << (*c)->id() << std::endl;
				utility_exit_with_message( msg.str() );
			}
			debug_assert( utility::pointer::dynamic_pointer_cast< components::StrandPairing const >( pair ) );
			components::StrandPairingCOP spair = utility::pointer::static_pointer_cast< components::StrandPairing const >( pair );
			shifts.push_back( spair->shift() );
		}
		prev = c;
	}
	return shifts;
}

/// @brief given a permutation, returns Orientations
BetaSheetArchitect::StrandOrientations
BetaSheetArchitect::retrieve_orientations( StructureData const & perm ) const
{
	using components::SegmentPairing;
	using components::SegmentPairingCOPs;

	typedef std::map< std::string, StrandOrientation > OrientationMap;

	OrientationMap orients;
	for ( SegmentPairingCOPs::const_iterator p=perm.pairings_begin(); p!=perm.pairings_end(); ++p ) {
		if ( (*p)->type() != components::SegmentPairing::STRAND ) continue;
		debug_assert( utility::pointer::dynamic_pointer_cast< components::StrandPairing const >( *p ) );
		components::StrandPairingCOP spair = utility::pointer::static_pointer_cast< components::StrandPairing const >( *p );
		debug_assert( spair->segments().size() == 2 );
		std::string const & seg1 = *( spair->segments().begin() );
		std::string const & seg2 = *( spair->segments().rbegin() );
		OrientationMap::iterator o = orients.find( seg1 );
		if ( (o != orients.end()) && (o->second != spair->orient1()) ) {
			std::stringstream msg;
			msg << class_name() << ": Different orientations are specified in the pairings for "
				<< seg1 << std::endl;
			msg << "SD=" << perm << std::endl;
			utility_exit_with_message( msg.str() );
		}
		orients[ seg1 ] = spair->orient1();

		o = orients.find( seg2 );
		if ( (o != orients.end()) && (o->second != spair->orient2()) ) {
			std::stringstream msg;
			msg << class_name() << ": Different orientations are specified in the pairings for "
				<< seg2 << std::endl;
			msg << "SD=" << perm << std::endl;
			utility_exit_with_message( msg.str() );
		}
		orients[ seg2 ] = spair->orient2();
	}

	// convert map
	StrandOrientations orientations;
	for ( StrandArchitectOP const & s : strands_ ) {
		OrientationMap::const_iterator o = orients.find( s->id() );
		if ( o == orients.end() ) {
			std::stringstream msg;
			msg << class_name() << ": No pairing information was found for strand " << s->id() << std::endl;
			msg << "SD=" << perm << std::endl;
			utility_exit_with_message( msg.str() );
		}
		orientations.push_back( o->second );
	}

	return orientations;
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
BetaSheetArchitect::check_permutation( components::StructureData const & perm_noextend ) const
{
	using core::select::residue_selector::ResidueSubset;
	using protocols::denovo_design::components::ResiduePair;
	using protocols::denovo_design::components::ResiduePairs;

	// add extensions
	components::StructureData perm = perm_noextend;
	for ( StrandExtension const & ext : extensions_ ) {
		std::stringstream ss, abego;
		ss << 'L';
		abego << 'X';
		for ( core::Size i=1; i<=ext.second+perm.segment( ext.first ).elem_length(); ++i ) {
			ss << 'E';
			abego << 'B';
		}
		ss << 'L';
		abego << 'X';
		components::Segment newseg( ext.first, ss.str(), abego.str(), false, false );
		perm.replace_segment( ext.first, newseg );
		TR.Debug << "Replacing " << ext.first << " with version extended by " << ext.second << std::endl;
	}

	TR.Debug << "Extended version: " << perm << std::endl;
	// check paired resiudes
	ResidueSubset subset( perm.pose_length(), false );
	ResiduePairs const residue_pairs = components::SegmentPairing::get_strand_residue_pairs( perm );
	TR.Debug << "Paired: " << residue_pairs << std::endl;

	for ( ResiduePair const & pair : residue_pairs ) {
		subset[ pair.first ] = true;
		subset[ pair.second ] = true;
	}

	// consider bulges to be paired
	for ( StrandArchitectOP const & strand : strands_ ) {
		components::Segment const & segment = perm.segment( strand->id() );
		ResidueVector const bulges = bulge_residues( perm, strand->id() );
		std::set< core::Size > const bulgeset( bulges.begin(), bulges.end() );
		for ( core::Size resid=segment.start(); resid<=segment.stop(); ++resid ) {
			if ( subset[ resid ] ) continue;
			if ( bulgeset.find( resid ) != bulgeset.end() ) continue;
			// BAD pairing!!
			std::stringstream msg;
			msg << "Returning false: strand " << strand->id() << ": segment "
				<< segment.pose_to_segment( resid ) << ": pose "
				<< resid << std::endl;
			throw EXCN_PreFilterFailed( msg.str() );
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
