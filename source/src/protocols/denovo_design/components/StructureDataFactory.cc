// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

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
#include <core/sequence/ABEGOManager.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableStringMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>

// Boost headers
#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.components.StructureDataFactory" );

namespace protocols {
namespace denovo_design {
namespace components {

int const
StructureDataFactory::REMARK_NUM = 994;

std::string const
StructureDataFactory::DATA_NAME = "PERMUTATION";

StructureDataFactory *
StructureDataFactory::create_singleton_instance()
{
	if ( ! basic::options::option[ basic::options::OptionKeys::run::preserve_header ].user() ) {
		TR.Warning << "-run:preserve_header is required for using \"Tomponents\" -- setting it to true."
			<< " To avoid this message, include -run:preserve header true in your flags" << std::endl;
		basic::options::option[ basic::options::OptionKeys::run::preserve_header ].value( true );
	}
	return new StructureDataFactory;
}

/// @brief stores the data of this permutation into a pose for later retrieval
///        StructureData stored in a pose this way can be retrieved by calling
///        get_from_pose(), or get_from_const_pose()
void
StructureDataFactory::save_into_pose( core::pose::Pose & pose, StructureData const & sd ) const
{
	sd.check_pose_consistency( pose );

	std::stringstream ss;
	ss << sd;
	set_cached_string( pose, ss.str() );

	// wipe out pdb info except for remarks
	pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( pose, true ) ) );

	// read cached remarks and save them as actual remarks
	debug_assert( pose.pdb_info() );
	pose.pdb_info()->remarks( sd.retrieve_remarks(pose) );

	//save_into_pose_with_id( pdb_info, "" );
	std::string line;
	while ( std::getline( ss, line ) ) {
		add_perm_remark( pose.pdb_info()->remarks(), line );
	}

	// eventually, I would like to make the cached observer the only place data is stored
	attach_observer( pose, sd );
}

void
StructureDataFactory::clear_from_pose( core::pose::Pose & pose ) const
{
	clear_cached_string( pose );
	pose.pdb_info()->remarks( core::io::Remarks() );
	detach_observer( pose );
}

/// @brief retrieves a StructureData object from the pose observable cache
///        utility_exit() if no StructureData object is found
StructureData const &
StructureDataFactory::get_from_const_pose( core::pose::Pose const & pose ) const
{
	StructureDataObserverCOP sd_obs = retrieve_observer( pose );
	if ( !sd_obs ) {
		std::stringstream msg;
		msg << "StructureDataFactory::get_from_const_pose(): No StructureData object was found in the pose!"
			<< " Be sure the StructureData observer is attached to the pose before trying to access it, or"
			<< " call get_from_pose() with a non-const pose instead." << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return sd_obs->sd();
}

/// @brief retrieves a StructureData object from the pose observable cache
///        creates one if necessary
StructureData const &
StructureDataFactory::get_from_pose( core::pose::Pose & pose ) const
{
	return get_from_pose( pose, "" );
}

/// @brief retrieves a StructureData object from observable cache
///        creates one if necessary, and sets it in the pose
StructureData const &
StructureDataFactory::get_from_pose(
	core::pose::Pose & pose,
	std::string const & newid ) const
{
	if ( !observer_attached( pose ) ) {
		StructureDataOP new_sd = create_new_from_pose( pose, newid );
		save_into_pose( pose, *new_sd );
	}

	StructureDataObserverCOP sd_obs = retrieve_observer( pose );
	if ( !sd_obs ) {
		std::stringstream msg;
		msg << "StructureData::get_from_pose(): Error obtaining cached StructureData observer" << std::endl;
		utility_exit_with_message( msg.str() );
	}

	return sd_obs->sd();
}

/// @brief retrieves a StructureData object from the pose observable cache
///        Creates a StructureData using the pose (but doesn't attach it) if
///        the cached StructureData could not be retrieved properly
StructureDataOP
StructureDataFactory::create_from_pose( core::pose::Pose const & pose ) const
{
	return create_from_pose( pose, "" );
}

/// @brief retrieves a StructureData object from the pose observable cache
///        Creates a StructureData using the pose (but doesn't attach it) if
///        the cached StructureData could not be retrieved properly
StructureDataOP
StructureDataFactory::create_from_pose(
	core::pose::Pose const & pose,
	std::string const & newid ) const
{
	if ( ! observer_attached( pose ) ) {
		return create_new_from_pose( pose, newid );
	}

	StructureDataObserverCOP sd_obs = retrieve_observer( pose );
	if ( !sd_obs ) {
		std::stringstream msg;
		msg << "StructureData::get_from_pose(): Error obtaining cached StructureData observer" << std::endl;
		utility_exit_with_message( msg.str() );
	}

	return sd_obs->sd_ptr()->clone();
}

void
add_virtual_residues( StructureData & sd, core::pose::Pose const & pose )
{
	core::Size virt_count = 1;
	for ( core::Size resid=1; resid<=pose.total_residue(); ++resid ) {
		if ( pose.residue( resid ).aa() != core::chemical::aa_vrt ) continue;
		Segment virt_seg( "L", "-", true, true );

		std::stringstream segname;
		segname << "VIRT" << std::setw(5) << std::setfill('0') << virt_count;
		++virt_count;

		if ( resid <= sd.pose_length() ) sd.add_segment( segname.str(), virt_seg, sd.segment_name( resid ) );
		else sd.add_segment( segname.str(), virt_seg );
	}
}

/// @brief sets data from a given pose's pdb remark records. Only looks at remarks which are subcomponents of the given id
/// if the pdb remarks are not found, the permutation information is inferred from the information that can be gleaned from the pose
StructureDataOP
StructureDataFactory::create_new_from_pose( core::pose::Pose const & pose, std::string const & id ) const
{
	StructureDataOP newperm;
	core::io::Remarks remarks;
	if ( has_cached_string(pose) ) {
		std::stringstream ss;
		ss << retrieve_cached_string( pose );
		TR.Debug << "Found StructureData information in datacache. Creating from that." << std::endl;
		newperm = create_from_xml( ss );
		if ( !newperm ) {
			std::stringstream err;
			err << "StructureDataFactory::create_new_from_pose():" << id
				<< ": could not create StructureData from XML." << std::endl;
			err << "XML: " << ss.str() << " chains: " << pose.conformation().chain_endings()
				<< " pose length: " << pose.total_residue() << std::endl;
			throw utility::excn::EXCN_Msg_Exception( err.str() );
		}
	} else {
		core::pose::PDBInfoCOP pdb_info = pose.pdb_info();
		if ( !pdb_info )  {
			TR << "No StructureData information was found in the pose -- Trying to infer info" << std::endl;
			newperm = infer_from_pose( pose, id );
		} else {
			newperm = create_from_remarks( pdb_info->remarks() );
			if ( !newperm ) {
				newperm = infer_from_pose( pose, id );
			}
		}
	}
	if ( !newperm ) {
		return NULL;
	}

	if ( pose.pdb_info() ) {
		TR.Debug << "Saving remarks into StructureData" << std::endl;
		newperm->save_remarks( pose.pdb_info()->remarks() );
	}

	// there may be virtual residues in this case
	if ( newperm->pose_length() != pose.total_residue() )  {
		add_virtual_residues( *newperm, pose );
	}

	newperm->check_pose_consistency( pose );
	return newperm;
}

class SegmentCounts {
public:
	enum SegmentType {
		HELIX = 1,
		STRAND = 2,
		LOOP = 3,
		LIGAND = 4,
		NUM_SEGMENT_TYPES = 5
	};

	SegmentCounts( core::pose::Pose const & pose ):
		counts_( NUM_SEGMENT_TYPES, 1 ),
		ligand_subset_( pose.total_residue(), false )
	{
		for ( core::Size resid=1; resid<=pose.total_residue(); ++resid ) {
			if ( pose.residue( resid ).is_ligand() ) {
				ligand_subset_[ resid ] = true;
			}
		}
	}

	std::string
	new_segment_name( char const ss_type, core::Size const seqpos )
	{
		SegmentType const t = type( ss_type, seqpos );
		std::stringstream name;
		name << ss_type;
		if ( t == LIGAND ) name << "IG";
		name << std::setw(2) << std::setfill('0') << counts_[ t ];
		++counts_[ t ];
		return name.str();
	}

private:
	SegmentType
	type( char const ss_type, core::Size const seqpos ) const
	{
		if ( ss_type == 'H' ) return HELIX;
		if ( ss_type == 'E' ) return STRAND;
		if ( ss_type == 'L' ) {
			if ( ligand_subset_[ seqpos ] ) return LIGAND;
			else return LOOP;
		}
		std::stringstream msg;
		msg << "SegmentCounts: unrecognized ss type: " << ss_type << std::endl;
		utility_exit_with_message( msg.str() );
		return NUM_SEGMENT_TYPES;
	};

private:
	utility::vector1< core::Size > counts_;
	core::select::residue_selector::ResidueSubset ligand_subset_;
};

std::string
add_segment(
	std::string const & id_val,
	StructureData & sd,
	Segment const & seg,
	char const ss_type,
	SegmentCounts & counts )
{
	std::stringstream name;
	if ( !id_val.empty() ) {
		name << id_val << PARENT_DELIMETER;
	}
	name << counts.new_segment_name( ss_type, sd.pose_length() + 1 );
	TR << "Creating new segment " << NamedSegment( name.str(), seg ) << std::endl;
	sd.add_segment( name.str(), seg );
	return name.str();
}

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
		Segment newseg( chain_ss, chain_abego, nterm_included, cterm_included );
		add_segment( id_val, sd, newseg, 'L', counts );
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
		SegmentNames::const_iterator next = ss;
		++next;
		bool const cterm_included = ( next != ss_segments.end() );

		Segment newseg( *ss, *ab, nterm_included, cterm_included );

		char const ss_type =  ( ss == ss_segments.begin() ) ? *(++(ss->begin())) : *ss->begin();
		std::string const name_str = add_segment( id_val, sd, newseg, ss_type, counts );
		if ( !prev_name.empty() ) {
			sd.mark_connected( prev_name, name_str );
		}
		prev_name = name_str;
	}
}

/// @brief sets data from a given pose's information, not taking into account PDB remarks
StructureDataOP
StructureDataFactory::infer_from_pose( core::pose::Pose const & pose, std::string const & id_val ) const
{
	// the object we will add to
	StructureDataOP sd( new StructureData( id_val ) );
	debug_assert( sd );

	if ( !pose.total_residue() ) {
		return sd;
	}

	// collect secondary structure
	core::scoring::dssp::Dssp dssp( pose );
	std::string const pose_ss = dssp.get_dssp_secstruct();
	debug_assert( pose_ss.size() == pose.total_residue() );

	// collect ABEGOS
	std::string const pose_abego = abego_str( core::sequence::ABEGOManager().get_symbols( pose, 1 ) );
	debug_assert( pose_abego.size() == pose.total_residue() );

	TR << "Pose ss = " << pose_ss << std::endl;
	TR << "Pose abego= " << pose_abego << std::endl;

	utility::vector1< core::Size > chain_endings = pose.conformation().chain_endings();
	chain_endings.push_back( pose.total_residue() );
	core::Size chain_start = 1;
	SegmentCounts counts( pose );
	for ( utility::vector1< core::Size >::const_iterator r=chain_endings.begin(); r!=chain_endings.end(); ++r ) {
		core::Size const chain_end = *r;

		// collect information about the residues from [ chain_start, chain_end ]
		core::Size const chain_length = chain_end - chain_start + 1;

		// get chain ss
		std::string const chain_ss = pose_ss.substr( chain_start - 1, chain_length );
		debug_assert( chain_ss.size() == chain_length );

		// get chain abego
		std::string chain_abego = pose_abego.substr( chain_start - 1, chain_length );
		debug_assert( chain_abego.size() == chain_length );
		// often, last residue is "O" -- change this to X
		if ( chain_abego[ chain_length-1 ] == 'O' ) {
			chain_abego[ chain_length-1 ] = 'X';
		}

		// create/add segments for chain
		add_segments_for_chain( id_val, *sd, chain_ss, chain_abego, counts );
		chain_start = *r + 1;
	}

	// locate non-polymeric covalent bonds
	for ( core::Size res=1; res<=pose.total_residue(); ++res ) {
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
				sd->add_covalent_bond( res, atom1, other_res, atom2 );
			}
		}
	}
	return sd;
}

/// @brief creates a StructureData from an xml stringstream
StructureDataOP
StructureDataFactory::create_from_xml( std::istream & xmltag ) const
{
	StructureDataOP newperm = StructureDataOP( NULL );
	utility::tag::TagCOP tag = utility::tag::Tag::create( xmltag );
	newperm = StructureDataOP( new StructureData( "StructureDataFactory::create_from_xml()" ) );
	newperm->parse_tag( tag ); // ID is set via XML
	return newperm;
}

StructureDataOP
StructureDataFactory::create_from_remarks( core::io::Remarks const & rem ) const
{
	TR.Debug << "Parsing " << rem.size() << " remarks!" << std::endl;
	// create list of strings
	utility::vector1< std::string > lines;
	for ( core::io::Remarks::const_iterator it_rem=rem.begin(); it_rem!=rem.end(); ++it_rem ) {
		TR.Debug << "checking " << *it_rem << std::endl;
		if ( it_rem->num != REMARK_NUM ) {
			continue;
		}
		lines.push_back( get_remark_line( it_rem, rem.end() ) );
	}

	if ( lines.empty() ) {
		TR.Debug << "No StructureData remark lines found." << std::endl;
		return StructureDataOP();
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
StructureDataFactory::attach_observer( core::pose::Pose & pose, StructureData const & sd ) const
{
	StructureDataObserverOP new_sd_obs( new StructureDataObserver( sd ) );
	pose.observer_cache().set( core::pose::datacache::CacheableObserverType::STRUCTUREDATA_OBSERVER, new_sd_obs );
	TR.Debug << "Set StructureDataObserver" << std::endl;
}

void
StructureDataFactory::detach_observer( core::pose::Pose & pose ) const
{
	pose.observer_cache().clear( core::pose::datacache::CacheableObserverType::STRUCTUREDATA_OBSERVER );
}

/// @brief stores a string in the pose's datacache
bool
StructureDataFactory::has_cached_string( core::pose::Pose const & pose ) const
{
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::STRING_MAP ) ) {
		return false;
	}
	basic::datacache::CacheableData const & cachable = pose.data().get( core::pose::datacache::CacheableDataType::STRING_MAP );
	debug_assert( dynamic_cast< basic::datacache::CacheableStringMap const * >(&cachable) == &cachable );
	basic::datacache::CacheableStringMap const & stringcache =
		static_cast< basic::datacache::CacheableStringMap const & >( cachable );
	std::map< std::string, std::string > const & smap = stringcache.map();
	return ( smap.find( DATA_NAME ) != smap.end() );
}

/// @brief stores a string in the pose's datacache
std::string
StructureDataFactory::retrieve_cached_string( core::pose::Pose const & pose ) const
{
	return retrieve_cached_string( pose, DATA_NAME );
}

std::string
StructureDataFactory::retrieve_cached_string( core::pose::Pose const & pose, std::string const & data_name ) const
{
	debug_assert( pose.data().has( core::pose::datacache::CacheableDataType::STRING_MAP ) );
	basic::datacache::CacheableData const & cachable = pose.data().get( core::pose::datacache::CacheableDataType::STRING_MAP );
	debug_assert( dynamic_cast< basic::datacache::CacheableStringMap const * >(&cachable) == &cachable );
	basic::datacache::CacheableStringMap const & stringcache =
		static_cast< basic::datacache::CacheableStringMap const & >( cachable );
	std::map< std::string, std::string > const & smap = stringcache.map();
	std::map< std::string, std::string >::const_iterator dat = smap.find( data_name );
	debug_assert( dat != smap.end() );
	std::string st = dat->second;
	clean_from_storage(st);
	return st;
}

/// @brief stores a string in the pose's datacache
void
StructureDataFactory::set_cached_string( core::pose::Pose & pose, std::string const & ssorig ) const
{
	set_cached_string( pose, ssorig, DATA_NAME );
}

void
StructureDataFactory::clear_cached_string( core::pose::Pose & pose ) const
{
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::STRING_MAP ) ) return;

	basic::datacache::CacheableData & cached = pose.data().get( core::pose::datacache::CacheableDataType::STRING_MAP );
	debug_assert( dynamic_cast< basic::datacache::CacheableStringMap * >(&cached) == &cached );
	basic::datacache::CacheableStringMap & smap =
		static_cast< basic::datacache::CacheableStringMap & >( cached );
	smap.map().erase( DATA_NAME );
}

void
StructureDataFactory::set_cached_string(
	core::pose::Pose & pose,
	std::string const & ssorig,
	std::string const & data_name ) const
{
	// swap out newlines and tabs
	std::string ss = ssorig;
	clean_for_storage( ss );
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::STRING_MAP ) ) {
		pose.data().set( core::pose::datacache::CacheableDataType::STRING_MAP,
			basic::datacache::CacheableDataOP( new basic::datacache::CacheableStringMap() ) );
	}
	debug_assert( pose.data().has( core::pose::datacache::CacheableDataType::STRING_MAP ) );
	basic::datacache::CacheableData & cached = pose.data().get( core::pose::datacache::CacheableDataType::STRING_MAP );
	debug_assert( dynamic_cast< basic::datacache::CacheableStringMap * >(&cached) == &cached );
	basic::datacache::CacheableStringMap & smap =
		static_cast< basic::datacache::CacheableStringMap & >( cached );
	smap.map()[data_name] = ss;
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

// Singleton instance and mutex static data members
namespace utility {

using protocols::denovo_design::components::StructureDataFactory;

#if defined MULTI_THREADED && defined CXX11
template<> std::mutex utility::SingletonBase< StructureDataFactory >::singleton_mutex_{};
template<> std::atomic< StructureDataFactory * > utility::SingletonBase< StructureDataFactory >::instance_( NULL );
#else
template<> StructureDataFactory * utility::SingletonBase< StructureDataFactory >::instance_( NULL );
#endif

} // namespace utility

