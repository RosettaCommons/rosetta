// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/components/StructureDataFactory.cc
/// @brief Singleton for creating StructureData objects
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/components/StructureDataFactory.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureDataObserver.hh>
#include <protocols/denovo_design/util.hh>

// Core headers
#include <core/conformation/Conformation.hh>
#include <core/io/Remarks.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/select/residue_selector/ResidueVector.hh>
#include <core/sequence/ABEGOManager.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/WriteableCacheableMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>

// Boost headers
#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>

static basic::Tracer TR( "protocols.denovo_design.components.StructureDataFactory" );

namespace protocols {
namespace denovo_design {
namespace components {

int const
	StructureDataFactory::REMARK_NUM = 994;

StructureDataFactory::StructureDataFactory() {
	if ( ! basic::options::option[ basic::options::OptionKeys::run::preserve_header ].user() ) {
		TR.Warning << "-run:preserve_header is required for using \"Tomponents\" -- setting it to true."
			<< " To avoid this message, include -run:preserve header true in your flags" << std::endl;
		// OUCH! Messing with the option system at this level is bad mojo. Don't do it!
		basic::options::option[ basic::options::OptionKeys::run::preserve_header ].value( true );
	}
}

/// @brief stores the data of this permutation into a pose for later retrieval
///        StructureData stored in a pose this way can be retrieved by calling
///        get_from_pose(), or get_from_const_pose()
void
StructureDataFactory::save_into_pose( core::pose::Pose & pose, StructureData const & sd ) const
{
	sd.check_pose_consistency( pose );

	set_cached_data( pose, sd );

	// wipe out pdb info except for remarks
	pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( pose, true ) ) );

	// read cached remarks and save them as actual remarks
	debug_assert( pose.pdb_info() );
	pose.pdb_info()->remarks( sd.retrieve_remarks(pose) );

	//save_into_pose_with_id( pdb_info, "" );
	std::stringstream ss;
	ss << sd;

	std::string line;
	while ( std::getline( ss, line ) ) {
		add_perm_remark( pose.pdb_info()->remarks(), line );
	}
}

void
StructureDataFactory::clear_from_pose( core::pose::Pose & pose ) const
{
	clear_cached_data( pose );
	pose.pdb_info()->remarks( core::io::Remarks() );
	detach_observer( pose );
}

/// @brief retrieves a StructureData object from the pose observable cache
///        utility_exit() if no StructureData object is found
StructureData const &
StructureDataFactory::get_from_const_pose( core::pose::Pose const & pose ) const
{
	return retrieve_cached_data( pose );
}

/// @brief retrieves a StructureData object from the pose observable cache
///        creates one if necessary
StructureData const &
StructureDataFactory::get_from_pose( core::pose::Pose & pose ) const
{
	if ( !has_cached_data( pose ) ) {
		save_into_pose( pose, create_from_pose( pose ) );
	}
	return get_from_const_pose( pose );
}

/// @brief retrieves a StructureData object from the pose observable cache
///        Creates a StructureData using the pose (but doesn't attach it) if
///        the cached StructureData could not be retrieved properly
StructureData
StructureDataFactory::create_from_pose( core::pose::Pose const & pose ) const
{
	return create_from_pose( pose, "" );
}

/// @brief retrieves a StructureData object from the pose data cache
///        Creates a StructureData using the pose (but doesn't store it) if
///        the cached StructureData is not present. Adds id_val as a prefix
///        for all new segments
/// @param[in]  pose    The input pose
/// @param[in]  prefix  A prefix to be added to the segments in the pose. For example, if
///                     prefix is empty, the first helix would be named 'H01', but if prefix
///                     is set to 'pose1', the first helix will be named 'pose.H01'
StructureData
StructureDataFactory::create_from_pose( core::pose::Pose const & pose, SegmentName const & prefix ) const
{
	// 1. If StructureData is cached, copy and return it
	if ( has_cached_data( pose ) ) {
		return retrieve_cached_data( pose );
	}

	// 2. Try to construct StructureData from pose and remarks
	StructureData newperm;
	core::io::Remarks remarks;
	core::pose::PDBInfoCOP pdb_info = pose.pdb_info();

	bool found_in_remarks = false;
	if ( pdb_info ) {
		TR << "No StructureData information was found in the pose -- Trying to infer info from remarks" << std::endl;
		try {
			newperm = create_from_remarks( pdb_info->remarks() );
			found_in_remarks = true;
		} catch ( EXCN_RemarksNotPresent const & e ) {
			e.show( TR );
			TR.flush();
		}
	}

	if ( !found_in_remarks ) {
		newperm = infer_from_pose( pose, prefix );
	}

	if ( pdb_info ) {
		TR.Debug << "Saving remarks into StructureData" << std::endl;
		newperm.save_remarks( pdb_info->remarks() );
	}

	newperm.check_pose_consistency( pose );
	return newperm;
}

SegmentCounts::SegmentCounts():
	counts_( NUM_SEGMENT_TYPES, 1 ),
	ligand_subset_()
{}

SegmentCounts::SegmentCounts( core::pose::Pose const & pose ):
	counts_( NUM_SEGMENT_TYPES, 1 ),
	ligand_subset_( pose.size(), false )
{
	for ( core::Size resid=1; resid<=pose.size(); ++resid ) {
		if ( pose.residue( resid ).is_ligand() ) {
			ligand_subset_[ resid ] = true;
		}
	}
}

std::string
SegmentCounts::new_segment_name( char const ss_type, core::Size const seqpos )
{
	SegmentType const t = type( ss_type, seqpos );
	std::stringstream name;
	name << ss_type;
	if ( t == LIGAND ) name << "IG";
	name << std::setw(2) << std::setfill('0') << counts_[ t ];
	++counts_[ t ];
	return name.str();
}

SegmentCounts::SegmentType
SegmentCounts::type( char const ss_type, core::Size const seqpos ) const
{
	if ( ss_type == 'H' ) return HELIX;
	if ( ss_type == 'E' ) return STRAND;
	if ( ss_type == 'L' ) {
		if ( !ligand_subset_.empty() && ligand_subset_[ seqpos ] ) return LIGAND;
		else return LOOP;
	}
	std::stringstream msg;
	msg << "SegmentCounts: unrecognized ss type: " << ss_type << std::endl;
	utility_exit_with_message( msg.str() );
	return NUM_SEGMENT_TYPES;
}

std::string
add_segment(
	StructureData & sd,
	Segment & seg,
	char const ss_type,
	SegmentCounts & counts )
{
	std::stringstream name;
	if ( !seg.id().empty() ) {
		name << seg.id() << PARENT_DELIMETER;
	}
	name << counts.new_segment_name( ss_type, sd.pose_length() + 1 );
	seg.set_id( name.str() );
	TR << "Creating new segment " << seg << std::endl;
	sd.add_segment( seg );
	return name.str();
}

/// @brief Adds segments to the given SD
/// @param[in]     id_val      Parent ID of the segments.  Can be empty
/// @param[in,out] sd          StructureData to be modified
/// @param[in]     chain_ss    Secondary structure of the full chain
/// @param[in]     chain_abego Abego for the full chain
/// @param[in,out] counts      Secondary structure element counters for naming segments
void
add_segments_for_chain(
	std::string const & id_val,
	StructureData & sd,
	std::string const & chain_ss,
	std::string const & chain_abego,
	SegmentCounts & counts )
{
	// if chain is < 3 residues, add it as a single segment
	if ( chain_ss.size() < 3 ) {
		// determine whether the termini need to be included
		// they are only included if the segment isn't long enough
		bool nterm_included = false;
		if ( chain_ss.size() < 3 ) {
			nterm_included = true;
		}
		bool cterm_included = false;
		if ( chain_ss.size() < 2 ) {
			cterm_included = true;
		}
		Segment newseg( id_val, chain_ss, chain_abego, nterm_included, cterm_included );
		add_segment( sd, newseg, 'L', counts );
		return;
	}
	// collect segment ss and abegos
	SegmentNames ss_segments, abego_segments;

	// collect non-terminal ss, abego
	std::string const elem_ss = chain_ss.substr( 1, chain_ss.size() - 2 );
	std::string const elem_abego = chain_abego.substr( 1, chain_abego.size() - 2 );

	std::string::const_iterator prev_ss = elem_ss.begin();
	std::string cur_ss = chain_ss.substr( 0, 2 );
	std::string cur_abego = chain_abego.substr( 0, 2 );
	for ( std::string::const_iterator s=++elem_ss.begin(), a=++elem_abego.begin(); s!=elem_ss.end(); ++s, ++a, ++prev_ss ) {
		if ( *s != *prev_ss ) {
			ss_segments.push_back( cur_ss );
			abego_segments.push_back( cur_abego );
			cur_ss = "";
			cur_abego = "";
		}
		cur_ss += *s;
		cur_abego += *a;
	}

	if ( !cur_ss.empty() ) {
		ss_segments.push_back( cur_ss + *chain_ss.rbegin() );
		abego_segments.push_back( cur_abego + *chain_abego.rbegin() );
	}
	TR << "Segments are " << ss_segments << " and " << abego_segments << std::endl;

	std::string prev_name = "";
	for ( SegmentNames::const_iterator ss=ss_segments.begin(), ab=abego_segments.begin();
			(ss!=ss_segments.end()) && (ab!=abego_segments.end()); ++ss, ++ab ) {
		bool const nterm_included = ( ss != ss_segments.begin() );
		auto next = ss;
		++next;
		bool const cterm_included = ( next != ss_segments.end() );

		Segment newseg( id_val, *ss, *ab, nterm_included, cterm_included );

		char const ss_type =  ( ss == ss_segments.begin() ) ? *(++(ss->begin())) : *ss->begin();
		std::string const name_str = add_segment( sd, newseg, ss_type, counts );
		if ( !prev_name.empty() ) {
			sd.mark_connected( prev_name, name_str );
		}
		prev_name = name_str;
	}
}

/// @brief creates a StructureData from a motif string
/// @param[in] motifs String of secstruct/abego motifs (e.g. "1LX-5EB-2LG-5EB-1LA-1LB-1LG-10HA-1LX" )
/// @returns StructureData containing one segment per motif with the specified
///          secondary structure and abego. Each motif is connected to
///          the previous and next motifs.  Segments are named by their
///          secondary structure (e.g. L01, E01, L02, E02, L03, H01, L04 )
StructureData
StructureDataFactory::create_from_motifs( std::string const & motif_str ) const
{
	return create_from_motifs( motif_str, "" );
}


/// @brief creates a StructureData from a motif string, optionally prefixing a string to the
///        name of each segment
/// @param[in] motifs String of secstruct/abego motifs (e.g. "1LX-5EB-2LG-5EB-1LA-1LB-1LG-10HA-1LX" )
/// @param[in] prefix Name to be prepended before the segment names. If empty,
///                   names will be as [ L01, E01, L02, ... ].  If set to "myprefix", names
///          will be as [ myprefix.L01, myprefix.E01, ... ]
/// @returns StructureData containing one segment per motif with the specified
///          secondary structure and abego. Each motif is connected to
///          the previous and next motifs.  Segments are named by their
///          secondary structure (e.g. L01, E01, L02, E02, L03, H01, L04 )
StructureData
StructureDataFactory::create_from_motifs( std::string const & motif_str, SegmentName const & prefix ) const
{
	StructureData sd( "StructureDataFactory::create_from_motifs()" );

	// determine SS for chain
	std::string chain_ss = "";
	std::string chain_abego = "";
	parse_motif_string( motif_str, chain_ss, chain_abego );

	// store SS counts
	SegmentCounts counts;

	// add segments to SD
	add_segments_for_chain( prefix, sd, chain_ss, chain_abego, counts );

	return sd;
}

/// @brief sets data from a given pose's information, not taking into account PDB remarks
/// @param[in]  pose    The input pose
/// @param[in]  prefix  A prefix to be added to the segments in the pose. For example, if
///                     prefix is empty, the first helix would be named 'H01', but if prefix
///                     is set to 'pose1', the first helix will be named 'pose.H01'
StructureData
StructureDataFactory::infer_from_pose( core::pose::Pose const & pose, SegmentName const & prefix ) const
{
	using core::select::residue_selector::ResidueVector;

	// the object we will add to
	StructureData sd( "StructureDataFactory::infer_from_pose()" );

	if ( !pose.size() ) {
		return sd;
	}

	// collect secondary structure
	core::scoring::dssp::Dssp dssp( pose );
	std::string const pose_ss = dssp.get_dssp_secstruct();
	debug_assert( pose_ss.size() == pose.size() );

	// collect ABEGOS
	std::string const pose_abego = abego_str( core::sequence::ABEGOManager().get_symbols( pose, 1 ) );
	debug_assert( pose_abego.size() == pose.size() );

	TR << "Pose ss = " << pose_ss << std::endl;
	TR << "Pose abego= " << pose_abego << std::endl;

	utility::vector1< core::Size > chain_endings = pose.conformation().chain_endings();
	chain_endings.push_back( pose.size() );
	core::Size chain_start = 1;
	core::Size cur_chain = 1;
	SegmentCounts counts( pose );
	for ( utility::vector1< core::Size >::const_iterator r=chain_endings.begin(); r!=chain_endings.end(); ++r, ++cur_chain ) {
		core::Size const chain_end = *r;

		// collect information about the residues from [ chain_start, chain_end ]
		core::Size const chain_length = chain_end - chain_start + 1;

		// look for non-polymers
		ResidueVector non_polymer;
		for ( core::Size resid=chain_start; resid<=chain_end; ++resid ) {
			if ( pose.residue( resid ).is_polymer() ) continue;
			non_polymer.push_back( resid );
		}
		TR.Debug << "Found non-polymer residues: " << non_polymer << std::endl;

		// get chain ss
		std::string chain_ss = pose_ss.substr( chain_start - 1, chain_length );
		debug_assert( chain_ss.size() == chain_length );

		// get chain abego
		std::string chain_abego = pose_abego.substr( chain_start - 1, chain_length );
		debug_assert( chain_abego.size() == chain_length );
		// often, last residue is "O" -- change this to X
		if ( chain_abego[ chain_length-1 ] == 'O' ) {
			chain_abego[ chain_length-1 ] = 'X';
		}

		// if non-polymer residues are present, segments must be separate
		core::Size prev_resid = chain_start;
		std::string::iterator ss = chain_ss.begin();
		std::string::iterator abego = chain_abego.begin();
		for ( ResidueVector::const_iterator non_poly=non_polymer.begin(); non_poly!=non_polymer.end(); ++non_poly ) {
			core::Size const local_resid = *non_poly - prev_resid + 1;
			std::string const local_ss( ss, ss + local_resid );
			std::string const local_abego( abego, abego + local_resid );
			add_segments_for_chain( prefix, sd, local_ss, local_abego, counts );
			TR << "Local ss = " << local_ss << " Chain ss = " << chain_ss << " Local_resid = " << local_resid << std::endl;
			prev_resid += local_ss.size();
			ss = ss + local_resid;
			abego = abego + local_resid;
		}

		std::string const local_ss( ss, chain_ss.end() );
		std::string const local_abego( abego, chain_abego.end() );
		if ( !local_ss.empty() && !local_abego.empty() ) {
			add_segments_for_chain( prefix, sd, local_ss, local_abego, counts );
		}
		chain_start = *r + 1;
	}

	// locate non-polymeric covalent bonds
	for ( core::Size res=1; res<=pose.size(); ++res ) {
		if ( pose.residue( res ).n_non_polymeric_residue_connections() ) {
			for ( core::Size conn=1; conn<=pose.residue( res ).n_possible_residue_connections(); ++conn ) {
				core::Size const other_res = pose.residue( res ).connected_residue_at_resconn( conn );
				if ( ( other_res == res + 1 ) || ( other_res + 1 == res ) ) continue;
				std::string const atom1 = pose.residue( res ).type().atom_name( pose.residue( res ).residue_connect_atom_index( conn ) );
				std::string atom2 = "";
				// find other connections
				for ( core::Size oconn=1; oconn<=pose.residue( other_res ).n_possible_residue_connections(); ++oconn ) {
					if ( pose.residue( other_res ).connected_residue_at_resconn( oconn ) == res ) {
						atom2 = pose.residue( other_res ).type().atom_name( pose.residue( other_res ).residue_connect_atom_index( oconn ) );
						break;
					}
				}
				if ( !other_res ) continue;
				TR << "Found non-polymeric connection between residues " << res << ":" << atom1 << " and " << other_res << ":" << atom2 << std::endl;
				debug_assert( other_res );
				debug_assert( !atom2.empty() );
				sd.add_covalent_bond( res, atom1, other_res, atom2 );
			}
		}
	}
	return sd;
}

StructureData
StructureDataFactory::create_from_cacheable_data( std::istream & raw_stream ) const
{
	std::istreambuf_iterator< char > eos;
	std::string xml_string( std::istreambuf_iterator< char >( raw_stream ), eos );

	TR.Debug << "Raw XML: " << xml_string << std::endl;
	clean_from_storage( xml_string );
	TR.Debug << "Clean XML: " << xml_string << std::endl;
	std::stringstream clean_xml( xml_string );
	return create_from_xml( clean_xml );
}

/// @brief creates a StructureData from an xml stringstream
StructureData
StructureDataFactory::create_from_xml( std::istream & xmltag ) const
{
	utility::tag::TagCOP tag = utility::tag::Tag::create( xmltag );
	StructureData newperm( "StructureDataFactory::create_from_xml()" );
	newperm.parse_tag( tag ); // ID is set via XML
	return newperm;
}

StructureData
StructureDataFactory::create_from_remarks( core::io::Remarks const & rem ) const
{
	TR.Debug << "Parsing " << rem.size() << " remarks!" << std::endl;
	// create list of strings
	utility::vector1< std::string > lines;
	for ( auto it_rem=rem.begin(); it_rem!=rem.end(); ++it_rem ) {
		TR.Debug << "checking " << *it_rem << std::endl;
		if ( it_rem->num != REMARK_NUM ) {
			continue;
		}
		lines.push_back( get_remark_line( it_rem, rem.end() ) );
	}

	if ( lines.empty() ) {
		throw CREATE_EXCEPTION(EXCN_RemarksNotPresent, "No StructureData remark lines found" );
	}

	// piece together full xml tag
	std::stringstream xmltag;
	xmltag << utility::join( lines, "\n" );
	TR.Debug << "XML tag: " << utility::join( lines, "\n" ) << std::endl;
	return create_from_xml( xmltag );
}

/// @brief checks for a StructureData in the pose observer cache
///        Returns true if one is found, false otherwise
bool
StructureDataFactory::observer_attached( core::pose::Pose const & pose ) const
{
	return pose.observer_cache().is_attached( core::pose::datacache::CacheableObserverType::STRUCTUREDATA_OBSERVER );
}

/// @brief returns observer pointer if the pose has one cached
///        pointer returned can be null if no StructureData is present in the cache
StructureDataObserverCOP
StructureDataFactory::retrieve_observer( core::pose::Pose const & pose ) const
{
	if ( !observer_attached( pose ) ) return StructureDataObserverCOP();

	core::pose::datacache::CacheableObserverCOP cache_obs =
		pose.observer_cache().get_const_ptr( core::pose::datacache::CacheableObserverType::STRUCTUREDATA_OBSERVER );
	StructureDataObserverCOP sd_obs = utility::pointer::static_pointer_cast< StructureDataObserver const >( cache_obs );
	return sd_obs;
}

/// @brief attaches cacheable observer.  Overwrites whatever was there before
void
StructureDataFactory::attach_observer( core::pose::Pose & pose ) const
{
	StructureDataOP sd_ptr( retrieve_cached_data_ptr( pose ) );

	StructureDataObserverOP new_sd_obs( new StructureDataObserver( sd_ptr ) );
	pose.observer_cache().set( core::pose::datacache::CacheableObserverType::STRUCTUREDATA_OBSERVER, new_sd_obs );
	TR.Debug << "Set StructureDataObserver" << std::endl;
}

void
StructureDataFactory::detach_observer( core::pose::Pose & pose ) const
{
	pose.observer_cache().clear( core::pose::datacache::CacheableObserverType::STRUCTUREDATA_OBSERVER );
}

basic::datacache::WriteableCacheableMap const &
get_writeable_cacheable_map( core::pose::Pose const & pose )
{
	using basic::datacache::WriteableCacheableMap;
	using core::pose::datacache::CacheableDataType;

	debug_assert( pose.data().has( CacheableDataType::WRITEABLE_DATA ) );
	basic::datacache::CacheableData const & cachable = pose.data().get( CacheableDataType::WRITEABLE_DATA );
	debug_assert( dynamic_cast< WriteableCacheableMap const * >(&cachable) == &cachable );
	auto const & writeable = static_cast< WriteableCacheableMap const & >( cachable );
	return writeable;
}

basic::datacache::WriteableCacheableMap &
get_writeable_cacheable_map( core::pose::Pose & pose )
{
	using basic::datacache::WriteableCacheableMap;
	using core::pose::datacache::CacheableDataType;

	debug_assert( pose.data().has( CacheableDataType::WRITEABLE_DATA ) );
	basic::datacache::CacheableData & cachable = pose.data().get( CacheableDataType::WRITEABLE_DATA );
	debug_assert( dynamic_cast< WriteableCacheableMap * >(&cachable) == &cachable );
	auto & writeable = static_cast< WriteableCacheableMap & >( cachable );
	return writeable;
}

/// @brief checks whether a pose contains a structuredata in its datacache
bool
StructureDataFactory::has_cached_data( core::pose::Pose const & pose ) const
{
	using basic::datacache::WriteableCacheableMap;
	using core::pose::datacache::CacheableDataType;

	if ( !pose.data().has( CacheableDataType::WRITEABLE_DATA ) ) return false;

	WriteableCacheableMap const & wmap = get_writeable_cacheable_map( pose );
	return ( wmap.map().find( StructureData::class_name() ) != wmap.map().end() );
}

StructureDataOP
StructureDataFactory::retrieve_cached_data_ptr( core::pose::Pose & pose ) const
{
	using basic::datacache::WriteableCacheableMap;
	using basic::datacache::WriteableCacheableDataOP;
	using core::pose::datacache::CacheableDataType;

	if ( !has_cached_data( pose ) ) {
		StructureData sd = create_from_pose( pose );
		save_into_pose( pose, sd );
	}

	if ( !pose.data().has( CacheableDataType::WRITEABLE_DATA ) ) {
		std::stringstream msg;
		msg << "StructureDataFactory::retrieve_cached_data_ptr(): No WriteableCacheableData was found in the pose!"
			<< " Be sure the StructureData is attached to the pose before trying to access it, or"
			<< " call get_from_pose() with a non-const pose instead." << std::endl;
		utility_exit_with_message( msg.str() );
	}

	WriteableCacheableMap const & wmap = get_writeable_cacheable_map( pose );
	auto dat = wmap.find( StructureData::class_name() );
	if ( dat == wmap.end() ) {
		std::stringstream msg;
		msg << "StructureDataFactory::retrieve_cached_data_ptr(): No StructureData object was found in the pose!"
			<< " Be sure the StructureData is attached to the pose before trying to access it, or"
			<< " call get_from_pose() with a non-const pose instead." << std::endl;
		utility_exit_with_message( msg.str() );
	}

	debug_assert( dat != wmap.end() );

	// a set of writeablecacheabledata pointers is stored in the map
	// we will enforce a set size of 1
	WriteableCacheableDataOP wdata = *dat->second.begin();
	debug_assert( utility::pointer::dynamic_pointer_cast< StructureData >( wdata ) == wdata );
	return utility::pointer::static_pointer_cast< StructureData >( wdata );
}

StructureData const &
StructureDataFactory::retrieve_cached_data( core::pose::Pose const & pose ) const
{
	using basic::datacache::WriteableCacheableMap;
	using basic::datacache::WriteableCacheableData;
	using core::pose::datacache::CacheableDataType;

	if ( !pose.data().has( CacheableDataType::WRITEABLE_DATA ) ) {
		std::stringstream msg;
		msg << "StructureDataFactory::retrieve_cached_data(): No WriteableCacheableData was found in the pose!"
			<< " Be sure the StructureData is attached to the pose before trying to access it, or"
			<< " call get_from_pose() with a non-const pose instead." << std::endl;
		utility_exit_with_message( msg.str() );
	}

	WriteableCacheableMap const & wmap = get_writeable_cacheable_map( pose );
	auto dat = wmap.find( StructureData::class_name() );
	if ( dat == wmap.end() ) {
		std::stringstream msg;
		msg << "StructureDataFactory::retrieve_cached_data(): No StructureData object was found in the pose!"
			<< " Be sure the StructureData is attached to the pose before trying to access it, or"
			<< " call get_from_pose() with a non-const pose instead." << std::endl;
		utility_exit_with_message( msg.str() );
	}
	debug_assert( dat != wmap.end() );

	// a set of writeablecacheabledata pointers is stored in the map
	// we will enforce a set size of 1
	WriteableCacheableData const & wdata = **dat->second.begin();
	debug_assert( dynamic_cast< StructureData const * >( &wdata ) == &wdata );
	return static_cast< StructureData const & >( wdata );
}

void
StructureDataFactory::clear_cached_data( core::pose::Pose & pose ) const
{
	using basic::datacache::WriteableCacheableMap;
	using core::pose::datacache::CacheableDataType;

	if ( !pose.data().has( CacheableDataType::WRITEABLE_DATA ) ) return;

	basic::datacache::CacheableData & cached = pose.data().get( CacheableDataType::WRITEABLE_DATA );
	debug_assert( dynamic_cast< WriteableCacheableMap * >( &cached ) == &cached );
	auto & wmap = static_cast< WriteableCacheableMap & >( cached );
	wmap.map().erase( StructureData::class_name() );
}

void
StructureDataFactory::set_cached_data( core::pose::Pose & pose, StructureData const & sd ) const
{
	using basic::datacache::WriteableCacheableDataOP;
	using basic::datacache::WriteableCacheableMap;
	using basic::datacache::WriteableCacheableMapOP;
	using core::pose::datacache::CacheableDataType;

	if ( !pose.data().has( CacheableDataType::WRITEABLE_DATA ) ) {
		pose.data().set( CacheableDataType::WRITEABLE_DATA, WriteableCacheableMapOP( new WriteableCacheableMap ) );
	}

	debug_assert( pose.data().has( CacheableDataType::WRITEABLE_DATA ) );
	basic::datacache::CacheableData & cached = pose.data().get( CacheableDataType::WRITEABLE_DATA );
	debug_assert( dynamic_cast< WriteableCacheableMap * >( &cached ) == &cached );
	auto & wmap = static_cast< WriteableCacheableMap & >( cached );

	// Create data set if necessary
	wmap[ StructureData::class_name() ].clear();
	wmap.insert( WriteableCacheableDataOP( new StructureData( sd ) ) );
}

/// @brief adds a remark to remarks object
void
StructureDataFactory::add_perm_remark( core::io::Remarks & remarks, std::string const & rem_value ) const
{
	add_remark( remarks, REMARK_NUM, rem_value );
}

void
clean_from_storage( std::string & st )
{
	// re-introduct newlines and tabs
	boost::replace_all( st, "%_S", " " );
	boost::replace_all( st, "%_T", "\t" );
	boost::replace_all( st, "%_N", "\n" );
}

void
clean_for_storage( std::string & ss )
{
	boost::replace_all( ss, "\n", "%_N" );
	boost::replace_all( ss, "\t", "%_T" );
	boost::replace_all( ss, " ", "%_S" );
}

} //protocols
} //denovo_design
} //components
