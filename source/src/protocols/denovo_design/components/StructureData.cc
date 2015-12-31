// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/components/StructureData.cc
/// @brief StructureData of a segment -- basically interfaces between segment and pose
/// @detailed
/// @author Tom Linsky

//Unit Headers
#include <protocols/denovo_design/components/StructureData.hh>

//Project Headers
#include <protocols/denovo_design/components/FoldGraph.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/util.hh>

//Protocol Headers
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/loops/Loop.hh>

//Core Headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/io/pdb/file_data.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/OptCysHG.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Remarks.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/sequence/ABEGOManager.hh>
#include <core/util/SwitchResidueTypeSet.hh>

//Basic/Utility/Numeric Headers
#include <basic/Tracer.hh>
#include <basic/datacache/CacheableStringMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

// Boost/ObjexxFCL Headers
#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/lexical_cast.hpp>

//C++ Headers

static basic::Tracer TR( "protocols.denovo_design.components.StructureData" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace components {

int const StructureData::REMARK_NUM = 994;
std::string const StructureData::DATA_NAME = "PERMUTATION";
char const StructureData::DATA_DELIMETER = '#';

StructureData::StructureData( std::string const & id_val ) :
	pose_(),
	id_( id_val ),
	ss_( "" ),
	pose_length_( 0 ),
	length_( 0 )
{
	segments_.clear();
	aliases_.clear();
	covalent_bonds_.clear();
	segment_order_.clear();
	abego_.clear();
}

/// @brief copy constructor -- poseOP is cloned
StructureData::StructureData( StructureData const & perm ) :
	utility::pointer::ReferenceCount(),
	pose_(),
	id_( perm.id_ ),
	ss_( perm.ss_ ),
	abego_( perm.abego_ ),
	pose_length_( perm.pose_length_ ),
	length_( perm.length_ ),
	data_int_( perm.data_int_ ),
	data_real_( perm.data_real_ ),
	data_str_( perm.data_str_ ),
	segments_( perm.segments_ ),
	aliases_( perm.aliases_ ),
	covalent_bonds_( perm.covalent_bonds_ ),
	segment_order_( perm.segment_order_ )
{
	if ( perm.pose() ) {
		pose_ = perm.pose()->clone();
	}
}

/// @brief destructors
StructureData::~StructureData()
{}

SingleChainStructureData::SingleChainStructureData( std::string const & id_val ) :
	StructureData( id_val )
{
}

SingleChainStructureData::SingleChainStructureData(
	std::string const & id_val,
	core::Size const ASSERT_ONLY(length_val),
	core::Size const pose_len_val,
	bool const is_loop,
	std::string const & ss_val,
	utility::vector1< std::string > const & abego_val ) :
	StructureData( id_val )
{
	bool nterm_included = false;
	core::Size n_anchor_res = 2;
	if ( pose_len_val <= 2 ) {
		n_anchor_res = 1;
		nterm_included = true;
	}
	debug_assert( n_anchor_res >= 1 );
	debug_assert( n_anchor_res <= pose_len_val );

	bool cterm_included = false;
	core::Size c_anchor_res = pose_len_val-1;
	if ( pose_len_val == 1 ) {
		c_anchor_res = 1;
		cterm_included = true;
	}
	debug_assert( c_anchor_res >= 1 );
	debug_assert( c_anchor_res <= pose_len_val );

	Segment resis(
		pose_len_val,                         // length
		( n_anchor_res + c_anchor_res ) / 2,  // safe residue
		0,                                    // cutpoint
		1,                                    // movable group
		is_loop,                              // is_loop
		nterm_included,                       // n terminus included in segment
		cterm_included,                       // c terminus included in segment
		"",                                   // lower connected segment
		"",                                   // upper connected segment
		ss_val,                               // secondary structure
		abego_val );                          // abego
	add_segment( id(), resis, "" );

	debug_assert( segment(id()).nterm_resi() == 1 );
	debug_assert( segment(id()).cterm_resi() == pose_length() );
	debug_assert( ss_val.size() == pose_length() );
	debug_assert( segment(id()).length() == pose_length() );
	debug_assert( length() == length_val );
}

SingleChainStructureData::~SingleChainStructureData()
{}

MultiChainStructureData::~MultiChainStructureData()
{}


/// @brief sets data from a given pose's information, not taking into account PDB remarks
StructureDataOP
StructureData::infer_from_pose( core::pose::Pose const & pose, std::string const & id_val )
{
	// the object we will add to
	StructureDataOP sd( new MultiChainStructureData( id_val ) );
	debug_assert( sd );

	if ( !pose.total_residue() ) {
		sd->set_pose( pose );
		return sd;
	}

	// collect secondary structure
	core::scoring::dssp::Dssp dssp( pose );
	std::string const pose_ss = dssp.get_dssp_secstruct();
	debug_assert( pose_ss.size() == pose.total_residue() );

	// collect ABEGOS
	StringVec const pose_abego = core::sequence::ABEGOManager().get_symbols( pose, 1 );
	debug_assert( pose_abego.size() == pose.total_residue() );

	TR << "Pose ss = " << pose_ss << std::endl;
	TR << "Pose abego= " << abego_str( pose_abego ) << std::endl;

	// find chains
	utility::vector1< core::Size > chain_endings = pose.conformation().chain_endings();
	// we will count the end of the pose as a chain ending too
	chain_endings.push_back( pose.total_residue() );

	Size chain_start = 1;
	Size const movable_group = 1;
	Size chain_num = 1;
	for ( utility::vector1< core::Size >::const_iterator r = chain_endings.begin();
			r != chain_endings.end();
			++r, ++chain_num ) {
		Size const chain_end = *r;

		// collect information about the residues from [ chain_start, chain_end ]
		Size const chain_length = chain_end - chain_start + 1;
		Size const local_saferes = ( (chain_end - chain_start) / 2 ) + 1;

		// determine whether the termini need to be included
		// they are only included if the segment isn't long enough
		bool nterm_included = false;
		if ( chain_length < 3 ) {
			nterm_included = true;
		}
		bool cterm_included = false;
		if ( chain_length < 2 ) {
			cterm_included = true;
		}

		std::string const chain_ss = pose_ss.substr( chain_start - 1, chain_length );
		debug_assert( chain_ss.size() == chain_length );

		StringVec chain_abego;
		for ( Size res = chain_start; res <= chain_end; ++res ) {
			chain_abego.push_back( pose_abego[ res ] );
		}
		debug_assert( chain_abego.size() == chain_length );
		// often, last residue is "O" -- change this to X
		if ( chain_abego[ chain_length ] == "O" ) {
			chain_abego[ chain_length ] = "X";
		}

		// name of new segment
		std::stringstream segname;
		if ( id_val.empty() ) {
			segname << chain_num;
		} else {
			// one chain, don't include chain number in name
			// multiple chains: name = ID.chain
			if ( pose.conformation().num_chains() == 1 ) {
				segname << id_val;
			} else {
				segname << id_val << PARENT_DELIMETER << chain_num;
			}
		}

		TR << "Creating new segment " << segname.str()
			<< ", length=" << chain_length << " safe=" << local_saferes
			<< " cutpoint=" << 0 << " mg=" << movable_group << " nterm_inc=" << nterm_included
			<< " cterm_inc=" << cterm_included << " ss=" << chain_ss << " abego=" << chain_abego
			<< std::endl;

		Segment new_segment(
			chain_length,                         // length
			local_saferes,                        // safe residue
			0,                                    // cutpoint
			movable_group,                        // movable group
			false,                                // is_loop
			nterm_included,                       // n terminus included in segment
			cterm_included,                       // c terminus included in segment
			"",                                   // lower connected segment
			"",                                   // upper connected segment
			chain_ss,                             // secondary structure
			chain_abego );                        // abego
		sd->add_segment( segname.str(), new_segment );
		chain_start = *r + 1;
	}

	// locate non-polymeric covalent bonds
	for ( core::Size res=1; res<=pose.total_residue(); ++res ) {
		if ( pose.residue( res ).n_non_polymeric_residue_connections() ) {
			for ( core::Size conn=1; conn<=pose.residue( res ).n_residue_connections(); ++conn ) {
				core::Size const other_res = pose.residue( res ).connected_residue_at_resconn( conn );
				if ( ( other_res == res + 1 ) || ( other_res + 1 == res ) ) continue;
				std::string const atom1 = pose.residue( res ).type().atom_name( pose.residue( res ).residue_connect_atom_index( conn ) );
				std::string atom2 = "";
				// find other connections
				for ( core::Size oconn=1; oconn<=pose.residue( other_res ).n_residue_connections(); ++oconn ) {
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
	sd->set_pose( pose );
	debug_assert( sd->pose() );
	debug_assert( sd->ss().size() == pose.total_residue() );
	return sd;
}

/// @brief sets data from a given pose's pdb remark records. Only looks at remarks which are subcomponents of the given id
/// if the pdb remarks are not found, the permutation information is inferred from the information that can be gleaned from the pose
StructureDataOP
StructureData::create_from_pose( core::pose::Pose const & pose, std::string const & id )
{
	StructureDataOP newperm;
	core::pose::Remarks remarks;
	if ( has_cached_string(pose) ) {
		std::stringstream ss;
		ss << cached_string( pose );
		TR.Debug << "Found StructureData information in datacache. Creating from that." << std::endl;
		newperm = create_from_xml( ss, id );
		if ( pose.total_residue() != newperm->pose_length() ) {
			std::stringstream err;
			err << newperm->id() << ": Size of StructureData does not match size of pose. ";
			err << "XML: " << ss.str() << " chains: " << pose.conformation().chain_endings() << " pose length: " << pose.total_residue() << std::endl;
			utility::vector1< core::Size > test;
			test[2] = 2;
			throw utility::excn::EXCN_Msg_Exception( err.str() );
		}
		newperm->set_pose( pose );
	} else {
		core::pose::PDBInfoCOP pdb_info = pose.pdb_info();
		if ( !pdb_info )  {
			TR << "No StructureData information was found in the pose -- Trying to infer info" << std::endl;
			newperm = infer_from_pose( pose, id );
		} else {
			newperm = create_from_remarks( pdb_info->remarks(), id );
			if ( newperm ) {
				newperm->set_pose( pose );
			} else {
				newperm = infer_from_pose( pose, id );
			}
			TR.Debug << "Saving remarks to datacache" << std::endl;
			newperm->save_remarks_to_datacache( pdb_info->remarks() );
		}
	}
	if ( !newperm ) {
		return NULL;
	}
	debug_assert( newperm->pose() );
	debug_assert( !newperm->pose()->total_residue() || newperm->ss().size() );
	debug_assert( !newperm->pose()->total_residue() || newperm->abego().size() );
	newperm->save_into_pose();
	return newperm;
}

StructureDataOP
StructureData::create_from_remarks( core::pose::Remarks const & rem, std::string const & newid )
{
	StructureDataOP newperm = NULL;
	core::Size read_lines = 0;
	for ( core::pose::Remarks::const_iterator it_rem=rem.begin(); it_rem != rem.end(); ++it_rem ) {
		if ( ( it_rem->num != 991 ) && ( it_rem->num != REMARK_NUM ) ) {
			continue;
		}
		++read_lines;
		// "new" xml method
		if ( it_rem->num == REMARK_NUM ) {
			newperm = parse_remarks( rem, newid );
			break;
		}
		// fallback to "old" serialize method
		std::string line = get_remark_line( it_rem, rem.end() );
		utility::vector1< std::string > fields = utility::string_split( line, ':' );
		debug_assert( fields.size() >= 2 );
		std::string & name = fields[1];
		std::string & type = fields[2];
		utility::vector1< std::string > names = utility::string_split( name, '.' );
		if ( names.size() == 1 ) {
			if ( type == "multi" ) {
				if ( fields[3] == "1" ) {
					debug_assert( !newperm );
					newperm = StructureDataOP( new MultiChainStructureData( name ) );
					newperm->load_pdb_info_old( rem, "" );
				} else if ( fields[3] == "0" ) {
					debug_assert( !newperm );
					newperm = StructureDataOP( new SingleChainStructureData( name ) );
					newperm->load_pdb_info_old( rem, "" );
				} else {
					TR.Error << "input error reading permutation data from pose. Line=" << line << std::endl;
					debug_assert( false );
					return NULL;
				}
			}
		}
	}
	return newperm;
}

/// @brief creates a permutation from pdb remarks
StructureDataOP
StructureData::parse_remarks( core::pose::Remarks const & rem, std::string const & newid )
{
	TR << "Parsing remarks!" << std::endl;
	// create list of strings
	utility::vector1< std::string > lines;
	for ( core::pose::Remarks::const_iterator it_rem=rem.begin(), it_end=rem.end(); it_rem != it_end; ++it_rem ) {
		if ( it_rem->num != REMARK_NUM ) {
			continue;
		}
		lines.push_back( get_remark_line( it_rem, it_end ) );
	}
	// create full xml tag
	std::stringstream xmltag;
	xmltag << utility::join( lines, "\n" );
	TR.Debug << "XML tag: " << utility::join( lines, "\n" ) << std::endl;
	return create_from_xml( xmltag, newid );
}

/// @brief creates a StructureData from an xml stringstream
StructureDataOP
StructureData::create_from_xml( std::istream & xmltag, std::string const & newid )
{
	StructureDataOP newperm = StructureDataOP( NULL );
	utility::tag::TagOP tag = utility::tag::Tag::create( xmltag );

	debug_assert( tag->getName() == "StructureData" );
	if ( tag->hasOption( "multi" ) ) {
		if ( tag->getOption< bool >( "multi" ) ) {
			newperm = StructureDataOP( new MultiChainStructureData( newid ) );
		} else {
			newperm = StructureDataOP( new SingleChainStructureData( newid ) );
		}
		for ( utility::vector0< utility::tag::TagCOP >::const_iterator t = tag->getTags().begin(), end = tag->getTags().end(); t != end; ++t ) {
			utility::tag::TagCOP subtag = *t;
			if ( subtag->getName() == "ResidueRange" ) {
				subtag->write( TR );
				debug_assert( subtag->hasOption( "name" ) );
				Segment newresis;
				newresis.parse_tag( *t );
				newperm->add_segment( subtag->getOption< std::string >( "name" ), newresis );
			} else if ( subtag->getName() == "Int" ) {
				debug_assert( subtag->hasOption( "name" ) );
				debug_assert( subtag->hasOption( "value" ) );
				newperm->set_data_int( subtag->getOption< std::string >( "name" ), subtag->getOption< int >( "value" ) );
			} else if ( subtag->getName() == "Real" ) {
				debug_assert( subtag->hasOption( "name" ) );
				debug_assert( subtag->hasOption( "value" ) );
				newperm->set_data_real( subtag->getOption< std::string >( "name" ), subtag->getOption< core::Real >( "value" ) );
			} else if ( subtag->getName() == "Str" ) {
				debug_assert( subtag->hasOption( "name" ) );
				if ( subtag->hasOption( "value" ) ) {
					newperm->set_data_str( subtag->getOption< std::string >( "name" ), subtag->getOption< std::string >( "value" ) );
				} else {
					newperm->set_data_str( subtag->getOption< std::string >( "name" ), "" );
				}
			} else if ( subtag->getName() == "Alias" ) {
				debug_assert( subtag->hasOption( "name" ) );
				debug_assert( subtag->hasOption( "segment" ) );
				debug_assert( subtag->hasOption( "res" ) );
				newperm->set_resnum_alias( subtag->getOption< std::string >( "name" ),
					subtag->getOption< std::string >( "segment" ),
					subtag->getOption< core::Size >( "res" ) );
			} else if ( subtag->getName() == "CovalentBond" ) {
				debug_assert( subtag->hasOption( "segment1" ) );
				debug_assert( subtag->hasOption( "segment2" ) );
				debug_assert( subtag->hasOption( "residue1" ) );
				debug_assert( subtag->hasOption( "residue2" ) );
				debug_assert( subtag->hasOption( "atom1" ) );
				debug_assert( subtag->hasOption( "atom2" ) );
				newperm->add_covalent_bond(
					subtag->getOption< std::string >( "segment1" ),
					subtag->getOption< core::Size >( "residue1" ),
					subtag->getOption< std::string >( "atom1" ),
					subtag->getOption< std::string >( "segment2" ),
					subtag->getOption< core::Size >( "residue2" ),
					subtag->getOption< std::string >( "atom2" )
				);
			} else {
				subtag->write( TR.Error );
				throw utility::excn::EXCN_RosettaScriptsOption( "Unknown tag in permutation: " + subtag->getName() );
			}
		}
	}
	return newperm;
}

/// @brief Saves remarks of the given pose into the pose's datacache -- changes enzdes residues to segment name/number
void
StructureData::save_remarks_to_datacache( core::pose::Remarks const & remarks )
{
	debug_assert( pose_ );
	core::Size remcount = 1;
	for ( core::pose::Remarks::const_iterator r=remarks.begin(), endr=remarks.end(); r != endr; ++r ) {
		if ( r->num == REMARK_NUM ) {
			continue;
		}

		if ( r->num == 666 ) {
			TR << "Processing enzdes header " << r->value << std::endl;
			//enzdes header
			utility::vector1< std::string > fields = utility::string_split_simple( r->value, ' ' );
			core::Size const resid1 = boost::lexical_cast< core::Size >( fields[5] );
			std::stringstream ss;
			ss << fields[1] << " " << fields[2] << " " << fields[3] << " " << fields[4] << " ";
			if ( resid1 ) {
				std::string const seg = segment_name( resid1 );
				debug_assert( resid1 >= segment(seg).start() );
				debug_assert( resid1 <= segment(seg).stop() );
				core::Size const localres = resid1 - segment(seg).start() + 1;
				ss << "%%" << seg << DATA_DELIMETER << localres << "%% ";
			} else {
				ss << fields[5] << " ";
			}
			ss << fields[6] << " " << fields[7] << " " << fields[8] << " " << fields[9] << " ";

			core::Size const resid2 = boost::lexical_cast< core::Size >( fields[10] );
			if ( resid2 ) {
				std::string const seg2 = segment_name( resid2 );
				debug_assert( resid2 >= segment(seg2).start() );
				debug_assert( resid2 <= segment(seg2).stop() );
				core::Size const localres2 = resid2 - segment(seg2).start() + 1;
				ss << "%%" << seg2 << DATA_DELIMETER << localres2 << "%% ";
			} else {
				ss << fields[10] << " ";
			}
			ss << fields[11] << " " << fields[12];
			set_cached_string( *pose_,
				boost::lexical_cast< std::string >(r->num) + ' ' + ss.str(),
				DATA_NAME + '.' + boost::lexical_cast< std::string >(remcount) );
		} else {
			set_cached_string( *pose_,
				boost::lexical_cast< std::string >(r->num) + ' ' + r->value,
				DATA_NAME + '.' + boost::lexical_cast< std::string >(remcount) );
		}
		++remcount;
	}
}

/// @brief loads data from pdb remarks into this permutation
void
StructureData::load_pdb_info_old(
	core::pose::Remarks const & rem,
	std::string const & prefix )
{
	for ( core::pose::Remarks::const_iterator it_rem=rem.begin(); it_rem!=rem.end(); ++it_rem ) {
		if ( it_rem->num != 991 ) {
			continue;
		}
		std::string line = get_remark_line( it_rem, rem.end() );
		utility::vector1< std::string > fields = utility::string_split( line, ':' );
		debug_assert( fields.size() >= 3);
		std::string const & name = fields[1];
		std::string const & type = fields[2];
		std::string const & field_name = fields[3];

		// check for sub-permutation of this one
		core::Size const idx = name.find( prefix + id() + "." );
		if ( idx != std::string::npos ) {
			continue;
		}

		if ( name != prefix + id() ) {
			TR.Debug << "ID not at the end. prefix=" << prefix << " id=" << id() << " name=" << name << std::endl;
			continue;
		}

		debug_assert( name == prefix + id() );
		if ( type == "multi" ) {
			continue;
		} else if ( type == "int" ) {
			debug_assert( fields.size() == 4 );
			TR << "Read id=" << name << " type=" << type << " param=" << field_name << " int=" << fields[4] << std::endl;
			int val = boost::lexical_cast< int >( fields[4] );
			if ( field_name == "pose_length" ) {
				// pose length will be automatically determined
				continue;
			} else if ( field_name == "length" ) {
				continue;
			} else {
				set_data_int( field_name, val );
			}
		} else if ( type == "real" ) {
			debug_assert( fields.size() == 4 );
			TR << "Read id=" << name << " type=" << type << " param=" << field_name << " real=" << fields[4] << std::endl;
			core::Real val = boost::lexical_cast< core::Real >( fields[4] );
			set_data_real( field_name, val );
		} else if ( type == "str" ) {
			debug_assert( fields.size() == 4 );
			TR << "Read id=" << name << " type=" << type << " param=" << field_name << " str=" << fields[4] << std::endl;
			set_data_str( field_name, fields[4] );
		} else if ( type == "ss" ) {
			// secondary structure is now automatically determined
			continue;
		} else if ( type == "residues" ) {
			debug_assert( fields.size() == 14 );
			core::Size startval = boost::lexical_cast< core::Size >( fields[4] );
			core::Size endval = boost::lexical_cast< core::Size >( fields[5] );
			core::Size const safe_resi = boost::lexical_cast< core::Size >( fields[6] );
			core::Size const movable_group = boost::lexical_cast< core::Size >( fields[7] );
			bool const is_loop = boost::lexical_cast< bool >( fields[8] );
			bool const nterm_included = boost::lexical_cast< bool >( fields[9] );
			bool const cterm_included = boost::lexical_cast< bool >( fields[10] );
			boost::erase_all( fields[11], " " );
			boost::erase_all( fields[12], " " );
			std::string const lower_segment = fields[11];
			std::string const upper_segment = fields[12];
			std::string const ss = fields[13];
			std::string const abego = fields[14];
			utility::vector1< std::string > abegovec;
			for ( core::Size idx=0; idx<abego.size(); ++idx ) {
				std::string news = "";
				news += abego[idx];
				abegovec.push_back( news );
			}
			if ( !abegovec.size() ) {
				for ( core::Size i=1; i<=ss.size(); ++i ) {
					std::string tmp = "";
					if ( ss[i-1] == 'E' ) {
						tmp += 'B';
					} else if ( ss[i-1] == 'H' ) {
						tmp += 'A';
					} else {
						tmp += 'X';
					}
					abegovec.push_back(tmp);
				}
			}
			debug_assert( abegovec.size() == ss.size() );
			TR << "Read id=" << name << " type=" << type << " param=" << name + '.' + field_name << " start=" << fields[4] << " end=" << fields[5] << " is_loop=" << is_loop << " nterm_inc=" << nterm_included << " cterm_inc=" << cterm_included << " mg=" << movable_group << " saferes=" << safe_resi <<  " lower_segment=" << lower_segment << " upper_segment=" << upper_segment << " ss=" << ss << " abego=" << abego << std::endl;
			if ( !nterm_included ) {
				--startval;
			}
			if ( !cterm_included ) {
				++endval;
			}
			core::Size len_val = endval - startval + 1;
			Segment resis( len_val, safe_resi - startval + 1, 0, movable_group, is_loop, nterm_included, cterm_included, lower_segment, upper_segment, ss, abegovec );
			// identify segment before this one
			StringList::iterator c, end;
			for ( c = segment_order_.begin(); c != segment_order_.end(); ++c ) {
				if ( segment( *c ).nterm_resi() >= startval ) {
					break;
				}
			}
			add_segment( field_name, resis, c );
		} else if ( type == "aliases" ) {
			debug_assert( fields.size() == 5 );
			core::Size const resi = boost::lexical_cast< core::Size >( fields[5] );
			set_resnum_alias( field_name, fields[4], resi );
			TR << "Set alias " << field_name << " --> " << fields[4] << ":" << resi << std::endl;
		} else {
			TR << "Unknown type: " << type << std::endl;
			debug_assert( false );
		}
	}
}

/// @brief stores the data in the permutation
void
StructureData::save_into_pose()
{
	if ( pose_ ) save_into_pose( *pose_ );
}

/// @brief stores the data of this permutation into a pose pdb remarks
void
StructureData::save_into_pose( core::pose::Pose & pose ) const
{
	check_consistency( pose );

	// eventually, I would like to make the cached string the only place permutation is stored
	std::stringstream ss;
	ss << *this;
	set_cached_string( pose, ss.str() );

	// wipe out pdb info except for remarks
	pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( pose, true ) ) );

	// read cached remarks and save them as actual remarks
	debug_assert( pose.pdb_info() );
	pose.pdb_info()->remarks( cached_remarks(pose) );

	//save_into_pose_with_id( pdb_info, "" );
	std::string line;
	while ( std::getline( ss, line ) ) {
		add_perm_remark( pose.pdb_info()->remarks(), line );
	}
}

utility::vector1< BondInfo >::const_iterator
StructureData::non_peptidic_bond( std::string const & seg1, std::string const & seg2 ) const
{
	for ( utility::vector1< BondInfo >::const_iterator bi=covalent_bonds_begin(); bi!=covalent_bonds_end(); ++bi ) {
		if ( ( ( seg1 == bi->seg1 ) && ( seg2 == bi->seg2 ) ) ||
				( ( seg2 == bi->seg1 ) && ( seg1 == bi->seg2 ) ) ) {
			return bi;
		}
	}
	return covalent_bonds_end();
}

utility::vector1< BondInfo >::const_iterator
StructureData::covalent_bonds_begin() const
{
	return covalent_bonds_.begin();
}

utility::vector1< BondInfo >::const_iterator
StructureData::covalent_bonds_end() const
{
	return covalent_bonds_.end();
}

/// @brief stores a string in the pose's datacache
bool
StructureData::has_cached_string() const
{
	debug_assert( pose_ );
	return has_cached_string( *pose_ );
}

/// @brief stores a string in the pose's datacache
bool
StructureData::has_cached_string( core::pose::Pose const & pose )
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
StructureData::cached_string() const
{
	debug_assert( pose_ );
	return cached_string( *pose_ );
}

/// @brief retrieves cached remarks from pose datacache
core::pose::Remarks
StructureData::cached_remarks() const
{
	debug_assert( pose_ );
	return cached_remarks( *pose_ );
}

void
clean_from_storage( std::string & st )
{
	// re-introduct newlines and tabs
	boost::replace_all( st, "%_S", " " );
	boost::replace_all( st, "%_T", "\t" );
	boost::replace_all( st, "%_N", "\n" );
}

/// @brief retrieves cached remarks from pose datacache
core::pose::Remarks
StructureData::cached_remarks( core::pose::Pose const & pose ) const
{
	debug_assert( pose.data().has( core::pose::datacache::CacheableDataType::STRING_MAP ) );
	basic::datacache::CacheableData const & cachable = pose.data().get( core::pose::datacache::CacheableDataType::STRING_MAP );
	debug_assert( dynamic_cast< basic::datacache::CacheableStringMap const * >(&cachable) == &cachable );
	basic::datacache::CacheableStringMap const & stringcache =
		static_cast< basic::datacache::CacheableStringMap const & >( cachable );
	std::map< std::string, std::string > const & smap = stringcache.map();
	core::pose::Remarks retval;
	for ( std::map< std::string, std::string >::const_iterator rem=smap.begin(); rem != smap.end(); ++rem ) {
		if ( boost::starts_with( rem->first, DATA_NAME + '.' ) ) {
			std::string val = rem->second;
			clean_from_storage( val );
			boost::trim_right(val);
			TR.Debug << "Cached line = " << val << std::endl;
			utility::vector1< std::string > fields = utility::string_split_simple( val );
			core::pose::RemarkInfo me;
			me.num = boost::lexical_cast< int >( fields[1] );
			std::stringstream ss;
			if ( me.num == 666 ) { // enzdes header
				for ( core::Size i=2; i<=fields.size(); ++i ) {
					if ( i>2 ) {
						ss << " ";
					}
					if ( i == 4 ) {
						std::stringstream res_stream( fields[6] );
						core::Size const resid_chain = boost::lexical_cast< core::Size >( substitute_variables( res_stream ) );
						if ( resid_chain ) {
							char chain = 'A' + pose.chain( resid_chain ) - 1;
							if ( pose.pdb_info() ) {
								char const pdbchain = pose.pdb_info()->chain( resid_chain );
								if ( pdbchain != '^' ) {
									chain = pdbchain;
								}
							}
							ss << chain;
						} else {
							ss << fields[6];
						}
					} else if ( i == 9 ) {
						std::stringstream res_stream( fields[11] );
						core::Size const resid_chain = boost::lexical_cast< core::Size >( substitute_variables( res_stream ) );
						if ( resid_chain ) {
							char chain = 'A' + pose.chain( resid_chain ) - 1;
							if ( pose.pdb_info() ) {
								char const pdbchain = pose.pdb_info()->chain( resid_chain );
								if ( pdbchain != '^' ) {
									chain = pdbchain;
								}
							}
							ss << chain;
						} else {
							ss << fields[11];
						}
					} else {
						ss << fields[i];
					}
				}
			} else { // non-enzdes header
				for ( core::Size i=2; i<=fields.size(); ++i ) {
					if ( i>2 ) {
						ss << " ";
					}
					ss << fields[i];
				}
			}
			me.value = substitute_variables( ss );
			retval.push_back( me );
		}
	}
	return retval;
}

/// @brief stores a string in the pose's datacache
std::string
StructureData::cached_string( core::pose::Pose const & pose )
{
	return cached_string( pose, DATA_NAME );
}

std::string
StructureData::cached_string( core::pose::Pose const & pose, std::string const & data_name )
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
StructureData::set_cached_string( std::string const & ss )
{
	debug_assert( pose_ );
	set_cached_string( *pose_, ss );
}

/// @brief stores a string in the pose's datacache
void
StructureData::set_cached_string( core::pose::Pose & pose, std::string const & ssorig )
{
	set_cached_string( pose, ssorig, DATA_NAME );
}

void
clean_for_storage( std::string & ss )
{
	boost::replace_all( ss, "\n", "%_N" );
	boost::replace_all( ss, "\t", "%_T" );
	boost::replace_all( ss, " ", "%_S" );
}

void
StructureData::set_cached_string( core::pose::Pose & pose, std::string const & ssorig, std::string const & data_name )
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
StructureData::add_perm_remark( core::pose::Remarks & remarks, std::string const & rem_value ) const
{
	add_remark( remarks, REMARK_NUM, rem_value );
}

/// @brief given an input stream, substitute all variables
/// @details variables are of the form: %%SEGMENTNAME#residue%%
/// SEGMENTNAME = name of the segment
/// residue = local residue number within the segment
/// The substituted value will be an core::Size corresponding to the pose residue
std::string
StructureData::substitute_variables( std::istream & input ) const
{
	std::stringstream sub_str;
	core::Size linecount = 0;
	// File format: %%SEGMENT_NAME#resid%% will give an integer resid
	while ( input.good() ) {
		std::string line = "";
		std::getline( input, line );
		if ( line[0] == '#' ) continue;
		++linecount;
		TR.Debug << "Line=" << line << std::endl;
		core::Size next_sub = line.find("%%");
		while ( next_sub != std::string::npos ) {
			core::Size second_sub = line.find("%%", next_sub+1);
			if ( second_sub == std::string::npos ) {
				throw utility::excn::EXCN_BadInput( "Malformed line in constraint file : " + line );
			}
			debug_assert( second_sub - next_sub >= 5 );
			std::string const variable = line.substr( next_sub+2, second_sub-next_sub-2 );
			utility::vector1< std::string > fields = utility::string_split( variable, DATA_DELIMETER );
			core::Size local_resid = 0;
			bool from_start = true;
			if ( fields.size() > 1 ) {
				std::string res_str = fields[ 2 ];
				if ( *(res_str.begin()) == '-' ) {
					from_start = false;
					res_str = res_str.substr( 1, std::string::npos );
				}
				local_resid = boost::lexical_cast< core::Size >( res_str );
			}
			core::Size actual_resid = 0;
			if ( from_start ) {
				actual_resid = pose_residue( fields[ 1 ], local_resid );
			} else {
				actual_resid = segment( fields[ 1 ] ).stop() - local_resid + 1;
				debug_assert( actual_resid >= segment( fields[ 1 ] ).start() );
			}

			line = line.substr(0,next_sub) + boost::lexical_cast< std::string >(actual_resid) + line.substr(second_sub+2,std::string::npos);
			next_sub = line.find("%%");
		}
		TR.Debug << "New line=" << line << std::endl;
		if ( linecount > 1 ) {
			sub_str << std::endl;
		}
		if ( !line.empty() ) {
			sub_str << line;
		}
	}
	return sub_str.str();
}

/// @brief returns the chain number of the given residue
core::Size
StructureData::chain( core::Size const resid ) const
{
	debug_assert( pose_ );
	return pose_->chain( resid );
}

/// @brief finds a segment in the segment_order list and returns an iterator to it
StringList::iterator
StructureData::find_segment( std::string const & segname )
{
	StringList::iterator c = segment_order_.begin();
	while ( c != segment_order_.end() ) {
		if ( *c == segname ) {
			break;
		}
		++c;
	}
	return c;
}

/// @brief finds a segment in the segment_order list and returns an iterator to it
StringList::const_iterator
StructureData::find_segment( std::string const & segname ) const
{
	StringList::const_iterator c = segment_order_.begin();
	while ( c != segment_order_.end() ) {
		if ( *c == segname ) {
			break;
		}
		++c;
	}
	return c;
}

/// @brief returns an ordered list of segments which are all connected containing seg
StringList
StructureData::connected_segments( std::string const & seg ) const
{
	StringList segmentlist;
	std::set< std::string > visited;
	std::stack< std::string > nodes;

	// go forward from segment
	nodes.push( seg );
	while ( nodes.size() ) {
		std::string const cur = nodes.top();
		nodes.pop();
		if ( visited.find( cur ) != visited.end() ) {
			continue;
		}
		visited.insert( cur );
		segmentlist.push_back( cur );
		Segment const & res = segment( cur );
		if ( !res.has_free_upper_terminus() ) {
			nodes.push( res.upper_segment() );
		}
	}

	// now go backward from segment
	if ( !segment( seg ).has_free_lower_terminus() ) {
		nodes.push( segment( seg ).lower_segment() );
	}
	while ( nodes.size() ) {
		std::string const cur = nodes.top();
		nodes.pop();
		if ( visited.find( cur ) != visited.end() ) {
			continue;
		}
		visited.insert( cur );
		segmentlist.push_front( cur );
		Segment const & res = segment( cur );
		if ( !res.has_free_lower_terminus() ) {
			nodes.push( res.lower_segment() );
		}
	}
	return segmentlist;
}

/// @brief marks the given segments as covanlently connected
void
StructureData::mark_connected(
	std::string const & lower_seg,
	std::string const & upper_seg )
{
	SegmentMap::iterator low = segments_.find( lower_seg );
	debug_assert( low != segments_.end() );
	SegmentMap::iterator up = segments_.find( upper_seg );
	debug_assert( up != segments_.end() );
	debug_assert( low->second.upper_segment() == "" );
	low->second.set_upper_segment( upper_seg );
	debug_assert( up->second.lower_segment() == "" );
	up->second.set_lower_segment( lower_seg );
	save_into_pose();
}

/// @brief determine pose chains based on termini
void
StructureData::chains_from_termini()
{
	if ( pose_ ) {
		pose_->conformation().chains_from_termini();
	}
}

/// @brief adds upper terminal variant to a residue
void
StructureData::add_upper_terminus_variant_type( core::Size const resi )
{
	debug_assert( pose_ );
	if ( !pose_->residue( resi ).is_upper_terminus() ) {
		core::pose::add_variant_type_to_pose_residue( *pose_, core::chemical::UPPER_TERMINUS_VARIANT, resi );
	}
}

/// @brief adds lower terminal variant to a residue
void
StructureData::add_lower_terminus_variant_type( core::Size const resi )
{
	debug_assert( pose_ );
	if ( !pose_->residue( resi ).is_lower_terminus() ) {
		core::pose::add_variant_type_to_pose_residue( *pose_, core::chemical::LOWER_TERMINUS_VARIANT, resi );
	}
}

/// @brief removes upper terminal variant from a residue
void
StructureData::remove_upper_terminus_variant_type( core::Size const resi )
{
	debug_assert( pose_ );
	if ( pose_->residue( resi ).is_upper_terminus() ) {
		core::pose::remove_variant_type_from_pose_residue( *pose_, core::chemical::UPPER_TERMINUS_VARIANT, resi );
	}
}

/// @brief removes lower terminal variant from a residue
void
StructureData::remove_lower_terminus_variant_type( core::Size const resi )
{
	debug_assert( pose_ );
	if ( pose_->residue( resi ).is_lower_terminus() ) {
		core::pose::remove_variant_type_from_pose_residue( *pose_, core::chemical::LOWER_TERMINUS_VARIANT, resi );
	}
}

/// @brief returns true of the last residue of segment1 contains a covalent bond to the first residue of segment2
bool
StructureData::polymer_bond_exists( std::string const & segment1, std::string const & segment2 ) const
{
	debug_assert( pose_ );
	return pose_->residue( segment(segment1).cterm_resi() ).is_polymer_bonded( segment(segment2).nterm_resi() );
}

/// @brief removes jump and cutpoint between the two segments to create a single polymer chain
void
StructureData::delete_jump_and_intervening_cutpoint( std::string const & segment1, std::string const & segment2 )
{
	if ( segment(segment1).cutpoint() ) {
		segments_[segment1].set_cutpoint( 0 );
	} else {
		if ( segment(segment2).cutpoint() ) {
			segments_[segment2].set_cutpoint( 0 );
		}
	}
	save_into_pose();
	if ( pose_ ) {
		core::Size const res1 = segment(segment1).cterm_resi();
		core::Size const res2 = segment(segment2).nterm_resi();
		int const jnum = find_jump(segment2);
		// check for cyclic case -- same jump for both segments.  Exit if cyclic
		if ( ( jnum == find_jump(segment1) ) || ( segment( segment1 ).cterm_resi() + 1 != segment( segment2 ).nterm_resi() ) ) {
			TR << "Can't delete jump/cutpoint between " << segment1 << ", res " << segment( segment1 ).safe() << " and " << segment2 << ", res " << segment( segment2 ).safe() << " because the polymer is apparently cyclic. Jump1=" << find_jump( segment1 ) << " Jump2=" << jnum << " FT=" << pose_->fold_tree() << " perm=" << *this << std::endl;
			return;
		}
		delete_jump_and_intervening_cutpoint( jnum, res1, res2 );
	}
}

/// @brief given a jump and a cutpoint residue, move the jump around the cut residue and delete to form a single edge
void
StructureData::delete_jump_and_intervening_cutpoint(
	int const jnum,
	core::Size const cut_resi1,
	core::Size const cut_resi2 )
{
	debug_assert( pose_ );

	// remove cutpoint variants from pose
	TR.Debug << "Removing cutpoints from " << cut_resi1 << " and " << cut_resi2 << std::endl;
	if ( pose_->residue(cut_resi1).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) {
		TR.Debug << "Removing lower cutpoint from " << cut_resi1 << std::endl;
		core::pose::remove_variant_type_from_pose_residue( *pose_, core::chemical::CUTPOINT_LOWER, cut_resi1 );
	}
	if ( pose_->residue(cut_resi2).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
		TR.Debug << "Removing upper cutpoint from " << cut_resi2 << std::endl;
		core::pose::remove_variant_type_from_pose_residue( *pose_, core::chemical::CUTPOINT_UPPER, cut_resi2 );
	}

	core::kinematics::FoldTree ft = pose()->fold_tree();
	// move jump so it surrounds cutpoint
	debug_assert( jnum <= static_cast< int >(ft.num_jump()) );
	debug_assert( jnum > 0 );
	TR.Debug << "Sliding jump " << jnum << " to " << cut_resi1 << "__" << cut_resi2 << std::endl;
	TR.Debug << pose_->fold_tree() << std::endl;
	ft.slide_jump( jnum, cut_resi1, cut_resi2 );
	// delete jump+cutpoint
	ft.delete_jump_and_intervening_cutpoint( jnum );
	debug_assert( ft.check_fold_tree() );
	pose_->fold_tree( ft );
}

void
StructureData::set_fold_tree( core::kinematics::FoldTree const & ft )
{
	debug_assert( ft.check_fold_tree() );
	if ( pose_ ) pose_->fold_tree( ft );
}

/// @brief switches residue type set of the contained pose
void
StructureData::switch_residue_type_set( std::string const & typeset )
{
	debug_assert( pose_ );
	core::util::switch_to_residue_type_set( *pose_, typeset );
}

/// @brief aligns the upper-terminal residue of segment1 to the start "anchor" residue of segment2
void
StructureData::align_segments(
	std::string const & segment1,
	std::string const & segment2 )
{
	core::Size const align_res_target( segment(segment1).cterm_resi() );
	core::Size const align_res_movable( segment(segment2).start() );
	core::Size const align_torsions = align_res_movable;
	int const jump_idx = find_jump( segment2 );
	TR.Debug << "aligning " << align_res_movable << " to " << align_res_target << " with torsions from " << align_torsions << std::endl;
	TR.Debug << "perm=" << *this << std::endl;
	align_residues( jump_idx, align_res_target, align_res_movable, align_torsions );
}

/// @brief aligns the lower-terminal residue of segment1 to the end "anchor" residue of segment2
void
StructureData::align_segments_rev(
	std::string const & segment1,
	std::string const & segment2 )
{
	core::Size const align_res_movable( segment(segment1).cterm_resi() );
	core::Size const align_res_target( segment(segment2).start() );
	core::Size const align_torsions = align_res_target;
	int const jump_idx = find_jump( segment2 );
	TR.Debug << "aligning " << align_res_movable << " to " << align_res_target << " with torsions from " << align_torsions << std::endl;
	TR.Debug << "perm=" << *this << std::endl;
	align_residues( jump_idx, align_res_target, align_res_movable, align_torsions );
}

/// @brief aligns the template_res'th residue of template_seg to the movable_res'th residue of movable_seg
/// re-arranges fold tree if necessary
void
StructureData::align_residues(
	std::string const & template_seg,
	core::Size const template_res,
	std::string const & movable_seg,
	core::Size const movable_res )
{
	debug_assert( template_res );
	debug_assert( template_res <= segment(template_seg).length() );
	debug_assert( movable_res );
	debug_assert( movable_res <= segment(movable_seg).length() );
	core::Size const align_res_target = segment(template_seg).resid(template_res);
	core::Size const align_res_movable = segment(movable_seg).resid(movable_res);
	core::Size const align_torsions = 0;

	TR.Debug << "aligning " << align_res_movable << " to " << align_res_target << " with torsions from " << align_torsions << std::endl;
	utility::vector1< std::string > roots;
	roots.push_back( template_seg );
	roots.push_back( movable_seg );
	consolidate_movable_groups( roots );
	int jump_idx = find_jump( movable_seg );
	debug_assert( jump_idx );

	align_residues( jump_idx, align_res_target, align_res_movable, align_torsions );
}

/// @brief replaces one residue with another
void
StructureData::replace_residue(
	std::string const & target_segment,
	core::Size const target_res,
	core::conformation::Residue const & res_in )
{
	core::Size const resnum = segment(target_segment).resid(target_res);
	replace_residue( resnum, res_in, true );
}

/// @brief replaces one residue with another
void
StructureData::replace_residue(
	core::Size const resnum,
	core::conformation::Residue const & res_in,
	bool const orient_bb = true )
{
	debug_assert( pose_ );
	bool is_upper_term = pose_->residue(resnum).is_upper_terminus();
	bool is_lower_term = pose_->residue(resnum).is_lower_terminus();
	pose_->replace_residue( resnum, res_in, orient_bb );
	pose_->energies().clear();
	pose_->conformation().update_polymeric_connection( resnum-1 );
	pose_->conformation().update_polymeric_connection( resnum );
	if ( !is_lower_term && pose_->residue(resnum).is_lower_terminus() ) {
		core::pose::remove_lower_terminus_type_from_pose_residue( *pose_, resnum );
	}
	if ( !is_upper_term && pose_->residue(resnum).is_upper_terminus() ) {
		core::pose::remove_upper_terminus_type_from_pose_residue( *pose_, resnum );
	}
}

/// @brief moves jump pointing to the first segment so that it will move with the second segment
void
StructureData::slide_jump(
	std::string const & child_segment,
	std::string const & parent_segment )
{
	if ( !pose_ ) {
		return;
	}
	int const childjump = find_jump(child_segment);
	debug_assert( childjump );
	int const parentjump = find_jump(parent_segment);
	if ( parentjump == childjump ) {
		TR.Warning << "Parent jump and child jump are the same while sliding jump for " << child_segment << " to have parent " << parent_segment << " --  not doing anything." << std::endl;
		return;
	}
	core::kinematics::FoldTree ft = pose_->fold_tree();
	core::kinematics::Edge jedge = ft.jump_edge(childjump);
	ft.slide_jump( childjump, segment(parent_segment).safe(), jedge.stop() );
	debug_assert( ft.check_fold_tree() );
	pose_->fold_tree( ft );
}

/// @brief aligns two residues so their backbones are completely superimposed
void
StructureData::align_residues(
	core::Size const jump_idx,
	core::Size const align_res_target,
	core::Size const align_res_movable,
	core::Size const res_with_torsions )
{
	debug_assert( pose_ );
	// move the jump temporarily
	core::kinematics::FoldTree ft = pose_->fold_tree();
	core::kinematics::FoldTree const ft_backup = ft;
	TR << "Sliding jump " << jump_idx << " to " << align_res_target << " --> " << align_res_movable << " in " << ft << std::endl;
	core::kinematics::Edge jedge = ft.jump_edge( jump_idx );
	ft.slide_jump( jump_idx, align_res_target, align_res_movable );
	debug_assert( ft.check_fold_tree() );
	pose_->fold_tree( ft );

	//setup torsions
	if ( res_with_torsions ) {
		pose_->set_phi( align_res_target, pose_->phi(res_with_torsions) );
		pose_->set_psi( align_res_target, pose_->psi(res_with_torsions) );
		pose_->set_phi( align_res_movable, pose_->phi(res_with_torsions) );
		pose_->set_psi( align_res_movable, pose_->psi(res_with_torsions) );
	}

	core::id::StubID stub1(
		core::id::AtomID(pose_->residue(align_res_target).type().atom_index("CA"), align_res_target),
		core::id::AtomID(pose_->residue(align_res_target).type().atom_index("N"), align_res_target),
		core::id::AtomID(pose_->residue(align_res_target).type().atom_index("C"), align_res_target) );

	core::id::StubID stub2(
		core::id::AtomID(pose_->residue(align_res_movable).type().atom_index("CA"), align_res_movable),
		core::id::AtomID(pose_->residue(align_res_movable).type().atom_index("N"), align_res_movable),
		core::id::AtomID(pose_->residue(align_res_movable).type().atom_index("C"), align_res_movable) );

	pose_->conformation().set_stub_transform( stub1, stub2, core::kinematics::RT() );
	debug_assert( ft_backup.check_fold_tree() );
	pose_->fold_tree( ft_backup );
}

/// @brief re-arranges residues such that segment 2 follows segment 1 in sequence
/// @details copies segment 2 to immediately after segment 1 and adds a jump
/// original segment 2 is then deleted.
void
StructureData::move_segment(
	std::string const & segment1_n,
	std::string const & segment1_c,
	std::string const & segment2_n,
	std::string const & segment2_c )
{
	if ( segment(segment1_c).cterm_resi()+1 != segment(segment2_n).nterm_resi() ) {
		core::Size s1 = segment(segment1_n).nterm_resi();
		core::Size e1 = segment(segment1_c).cterm_resi();
		core::Size s2 = segment(segment2_n).nterm_resi();
		core::Size e2 = segment(segment2_c).cterm_resi();
		TR << "moving segment " << segment2_n << "(" << s2 << ") <--> " << segment2_c << "(" << e2 << ")" << " to after " << segment1_c << "(" << e1 << ")" << std::endl;
		TR << "Before moving SD: " << *this << std::endl;
		try {
			check_consistency();
		} catch ( EXCN_PoseInconsistent const & e ) {
			TR.Error << "Input to move_segment() is not consistent with pose -- dumping to predump.pdb" << std::endl;
			if ( pose_ ) pose_->dump_pdb( "predump.pdb" );
			throw e;
		}

// find all involved segments
		StringList::iterator segment2_segment_begin = find_segment( segment2_n );
		StringList::iterator segment2_segment_end = find_segment( segment2_c );
		debug_assert( segment2_segment_begin != segment_order_.end() );
		debug_assert( segment2_segment_end != segment_order_.end() );
		++segment2_segment_end;

		if ( pose_ ) {
			move_jumps_to_safety();
			move_segment_in_pose( s1, e1, s2, e2 );
		}

		StringList segment2_segments( segment2_segment_begin, segment2_segment_end );
		debug_assert( segment_order_.size() == segments_.size() );
		segment_order_.erase( segment2_segment_begin, segment2_segment_end );
		StringList::iterator insert_pos = find_segment( segment1_c );
		debug_assert( insert_pos != segment_order_.end() );
		++insert_pos;
		segment_order_.insert( insert_pos, segment2_segments.begin(), segment2_segments.end() );
		debug_assert( segment_order_.size() == segments_.size() );
		update_numbering();
		try {
			check_consistency();
		} catch ( EXCN_PoseInconsistent const & e ) {
			TR.Error << "SD inconsistent with pose. Dumping to dump.pdb. SD=" << *this << std::endl;
			if ( pose_ ) pose_->dump_pdb( "dump.pdb" );
			throw e;
		}
		move_jumps_to_safety();
	}
}

/// @brief moves a segment of the pose such that the segment from start2<=res<=end2 is moved so that it starts at end1+1
/// returns the jump number that is to be ignored in fold tree searches
void
StructureData::move_segment_in_pose(
	core::Size start1,
	core::Size end1,
	core::Size start2,
	core::Size end2 )
{
	TR << "moving segment " << start2 << "->" << end2 << " to after residues " << start1 << "->" << end1 << std::endl;

	debug_assert( pose_ );
	int const len = end2 - start2 + 1;

	insert_after_residue_in_pose( start1, end1, start2, end2 );
	debug_assert( len >= 0 );
	if ( end1 < start2 ) {
		start2 += len;
		end2 += len;
	}
	delete_residues_in_pose( start2, end2 );
}

/// @brief copies and inserts a segment into the permutation after the given segment
void
StructureData::insert_after_residue_in_pose(
	core::Size segment1_start,
	core::Size segment1_end,
	core::Size segment2_start,
	core::Size segment2_end )
{
	debug_assert( pose_ );
	TR.Debug << "Inserting segment -- original segment2: " << segment2_start << "-->" << segment2_end << ", segment1_end=" << segment1_end << std::endl;

	// add segment2 residues
	for ( core::Size cseg_res=segment2_start; cseg_res<=segment2_end; ++cseg_res ) {
		TR.Debug << "copying " << cseg_res << " to " << segment1_end+1 << std::endl;
		// new residue will be attached by jump
		if ( cseg_res != segment2_start ) {
			if ( pose()->residue(cseg_res).is_lower_terminus() ) {
				pose_->insert_residue_by_jump( pose_->residue( cseg_res ), segment1_end+1, pose_->fold_tree().root() );
			} else {
				//remove_lower_terminus_variant_type( cseg_res );
				pose_->append_polymer_residue_after_seqpos( pose_->residue(cseg_res), segment1_end, false );
			}
		} else {
			int upstream_res = pose_->fold_tree().root();
			TR.Debug << "Adding jump from " << upstream_res << " after " << segment1_end << " to residue copied from " << cseg_res << std::endl;
			core::kinematics::FoldTree const & ft = pose_->fold_tree();
			for ( core::Size i=1; i<=pose_->total_residue(); ++i ) {
				TR.Debug << "Res " << i << " : " << ft.is_cutpoint(i) << std::endl;
			}
			pose_->insert_residue_by_jump( pose_->residue(cseg_res), segment1_end+1, upstream_res );
			core::pose::add_lower_terminus_type_to_pose_residue( *pose_, segment1_end+1 );
		}
		if ( segment1_end < segment2_start ) {
			cseg_res += 1;
			segment2_start += 1;
			segment2_end += 1;
		}
		segment1_end += 1;
	}
	TR.Debug << "after adding, segment1_end=" << segment1_end << " segment2_start=" << segment2_start << " segment2_end=" << segment2_end << std::endl;

	// fix the fold tree -- jumps without label "jump" which are in c2 need to be reset
	core::kinematics::FoldTree ft( pose()->fold_tree() );
	TR.Debug << "Initial foldtree: " << ft << std::endl;
	utility::vector1< core::kinematics::Edge > new_edges;
	for ( int j=1; j<=static_cast<int>(ft.num_jump()); ++j ) {
		int const u = ft.upstream_jump_residue(j);
		int const d = ft.downstream_jump_residue(j);
		// these are important for the casting below
		debug_assert( u >= 0 );
		debug_assert( d >= 0 );
		TR.Debug << "Looking at jump " << j << " upstream " << u << " downstream " << d << " Comp2=" << segment2_start << "->" << segment2_end << std::endl;
		bool const u_in_segment2 = ( segment2_start <= static_cast< core::Size >(u) ) && ( static_cast< core::Size >(u) <= segment2_end );
		bool const d_in_segment1 = ( segment1_start <= static_cast< core::Size >(d) ) && ( static_cast< core::Size >(d) <= segment1_end );
		core::kinematics::Edge & new_edge = ft.jump_edge(j);
		debug_assert( new_edge.start() == u );
		debug_assert( new_edge.stop() == d );
		// 4 cases
		if ( u_in_segment2 && d_in_segment1 ) {
			// 1. upstream in cut segment, downstream in new segment
			TR.Debug << "reversing jump " << j << " from newup=" << d << " to " << " newdown=" << u << std::endl;
			new_edge.start() = d;
			new_edge.stop() = u;
		} else if ( u_in_segment2 && (!d_in_segment1) ) {
			// 2. upstream in cut segment, downstream not in new segment
			core::Size const target_res = segment1_end - ( segment2_end - ft.upstream_jump_residue(j) );
			TR.Debug << "sliding up jump " << j << " from " << u << " to " << target_res << " down=" << d << std::endl;
			new_edge.start() = target_res;
			insert_peptide_edges( ft, new_edge );
		} else {
			// 4. Neither residue is involved in this addition -- don't touch
			TR.Debug << "ignoring jump " << j << " up=" << u << " down=" << d << std::endl;
		}
	}

	TR.Debug << "fold tree after insert : " << ft << std::endl;
	ft.delete_extra_vertices();
	TR.Debug << "cleaned fold tree      : " << ft << std::endl;
	debug_assert( ft.check_fold_tree() );
	pose_->fold_tree( ft );
}

/// @brief expands the segment so that the trailing pad residue(s) become part of the segment
void
StructureData::engulf_leading_residues( std::string const & seg )
{
	SegmentMap::iterator r = segments_.find( seg );
	debug_assert( r != segments_.end() );
	r->second.engulf_leading_residues();
	update_numbering();
}

/// @brief expands the segment so that the trailing pad residue(s) become part of the segment
void
StructureData::engulf_trailing_residues( std::string const & seg )
{
	SegmentMap::iterator r = segments_.find( seg );
	debug_assert( r != segments_.end() );
	r->second.engulf_trailing_residues();
	update_numbering();
}

/// @brief deletes the residues between the segment N terminus and the N anchor point
void
StructureData::delete_leading_residues( std::string const & seg )
{
	SegmentMap::iterator r = segments_.find( seg );
	debug_assert( r != segments_.end() );
	core::Size const start = r->second.nterm_resi();
	core::Size const stop = segment(seg).start() - 1;
	TR.Debug << "Lead delete start=" << start << " stop=" << stop << std::endl;
	if ( start == r->second.start() ) {
		return;
	}
	if ( pose_ ) {
		delete_residues_in_pose( start, stop );
	}
	r->second.delete_leading_residues();
	update_numbering();
}

/// @brief deletes the residues between the segment C terminus and the C anchor point
void
StructureData::delete_trailing_residues( std::string const & seg )
{
	SegmentMap::iterator r = segments_.find( seg );
	debug_assert( r != segments_.end() );
	core::Size const start = r->second.stop() + 1;
	core::Size const stop = r->second.cterm_resi();
	TR.Debug << "Trail delete start=" << start << " stop=" << stop << std::endl;
	if ( start > stop ) {
		return;
	}
	if ( pose_ ) {
		delete_residues_in_pose( start, stop );
	}
	r->second.delete_trailing_residues();
	update_numbering();
}

void
StructureData::delete_segment_nopose( std::string const & seg_val, SegmentMap::iterator r )
{
	std::string const segment1_c = r->second.lower_segment();
	std::string const segment2_n = r->second.upper_segment();
	if ( !r->second.has_free_lower_terminus() ) {
		SegmentMap::iterator r2 = segments_.find( segment1_c );
		debug_assert( r2 != segments_.end() );
		r2->second.set_upper_segment( "" );
	}
	if ( !r->second.has_free_upper_terminus() ) {
		SegmentMap::iterator r2 = segments_.find( segment2_n );
		debug_assert( r2 != segments_.end() );
		r2->second.set_lower_segment( "" );
	}
	StringList::iterator remove_me = find_segment( seg_val );
	debug_assert( remove_me != segment_order_.end() );
	segment_order_.erase( remove_me );
	core::Size const removed_mg = r->second.movable_group;
	segments_.erase( r );
	debug_assert( segment_order_.size() == segments_.size() );

	update_movable_groups_after_deletion( removed_mg );
}

/// @brief removes all traces of the given segment from the object
void
StructureData::delete_segment( std::string const & seg_val )
{
	TR.Debug << "delete_segment(" << seg_val << ")" << std::endl;
	std::string seg = seg_val;
	SegmentMap::iterator r = segments_.find( seg );
	if ( r == segments_.end() ) {
		seg = id() + PARENT_DELIMETER + seg_val;
		r = segments_.find( seg );
	}
	debug_assert( r != segments_.end() );

	// disconnect
	core::Size start_del = r->second.nterm_resi();
	core::Size stop_del = r->second.cterm_resi();
	std::string const segment1_c = r->second.lower_segment();
	std::string const segment2_n = r->second.upper_segment();
	if ( !r->second.has_free_lower_terminus() ) {
		SegmentMap::iterator r2 = segments_.find( segment1_c );
		debug_assert( r2 != segments_.end() );
		if ( start_del <= stop_del ) {
			start_del += 1;
			r2->second.set_cterm_included( false );
		}
	}
	if ( !r->second.has_free_upper_terminus() ) {
		SegmentMap::iterator r2 = segments_.find( segment2_n );
		debug_assert( r2 != segments_.end() );
		if ( start_del <= stop_del ) {
			debug_assert( stop_del > 1 );
			stop_del -= 1;
			r2->second.set_nterm_included( false );
		}
	}
	if ( pose_  && ( start_del <= stop_del ) ) {
		core::Size const resi_count = pose_->total_residue();
		delete_residues_in_pose( start_del, stop_del );
		debug_assert( resi_count - (stop_del-start_del+1) == pose_->total_residue() );
	}
	delete_segment_nopose( seg, r );
	update_numbering();

	if ( pose_ ) {
		if ( !segment1_c.empty() ) {
			add_upper_terminus_variant_type( segment( segment1_c ).cterm_resi() );
		}
		if ( !segment2_n.empty() ) {
			add_lower_terminus_variant_type( segment( segment2_n ).nterm_resi() );
		}
		chains_from_termini();
		utility::vector1< std::string > roots;
		if ( !segment1_c.empty() ) roots.push_back( segment1_c );
		if ( !segment2_n.empty() ) roots.push_back( segment2_n );
		if ( roots.empty() ) {
			move_jumps_to_safety();
		} else {
			consolidate_movable_groups( roots );
		}
	}
	TR << "Deleted segment " << seg << std::endl;
}

/// @brief deletes the given residues from the pose
void
StructureData::delete_residues_in_pose(
	core::Size const start,
	core::Size const end )
{
	debug_assert( pose_ );
	pose_->conformation().buffer_signals();
	TR.Debug << "Deleting pose residues " << start << " to " << end << std::endl;
	if ( start > end ) {
		std::stringstream ss;
		ss << id() << ": deleting residue from " << start << " to " << end
			<< ", pose length = " << pose_->total_residue() << std::endl;
		ss << "Bad start and end given to delete_residues_in_pose()" << std::endl;
		throw utility::excn::EXCN_BadInput( ss.str() );
	}
	pose_->conformation().delete_residue_range_slow( start, end );
	pose_->conformation().unblock_signals();
}

/// @brief sets jump with index jumpidx to datastructure j
void
StructureData::set_jump( int const jumpidx, core::kinematics::Jump const & j )
{
	if ( pose_ ) {
		pose_->set_jump( jumpidx, j );
	}
}

/// @brief sets phi for a specific residue
void
StructureData::set_phi( core::Size const seqpos, core::Real const phi_val )
{
	pose_->set_phi( seqpos, phi_val );
}

/// @brief sets psi for a specific residue
void
StructureData::set_psi( core::Size const seqpos, core::Real const psi_val )
{
	pose_->set_psi( seqpos, psi_val );
}

/// @brief sets omega for a specific residue
void
StructureData::set_omega( core::Size const seqpos, core::Real const omega_val )
{
	pose_->set_omega( seqpos, omega_val );
}

/// @brief sets bond length in pose
void
StructureData::set_bond_length(
	std::string const & segmentname,
	core::Size const resid,
	std::string const & atom1,
	std::string const & atom2,
	core::Real const newlength )
{
	debug_assert( pose_ );
	debug_assert( newlength >= 0 );
	debug_assert( has_segment( segmentname ) );
	Segment const & res = segment(segmentname);
	debug_assert( resid );
	debug_assert( resid <= res.length() );
	core::Size const poseres = res.resid(resid);
	debug_assert( poseres );
	debug_assert( poseres <= pose_->total_residue() );
	debug_assert( pose_->residue(poseres).type().has( atom1 ) );
	debug_assert( pose_->residue(poseres).type().has( atom2 ) );
	core::id::AtomID const atom1id( pose_->residue(poseres).type().atom_index(atom1), poseres );
	core::id::AtomID const atom2id( pose_->residue(poseres).type().atom_index(atom2), poseres );
	pose_->conformation().set_bond_length( atom1id, atom2id, newlength );

}

/// @brief (re-)detect disulfides, will convert CYD to CYS if disulfide bond is lost
void
StructureData::detect_disulfides( core::scoring::ScoreFunctionOP sfx )
{
	debug_assert( pose_ );
	debug_assert( sfx );
	if ( basic::options::option[ basic::options::OptionKeys::in::detect_disulf ].user() ?
			basic::options::option[ basic::options::OptionKeys::in::detect_disulf ]() : // detect_disulf true
			pose_->is_fullatom() // detect_disulf default but fa pose
			) {
		pose_->conformation().detect_disulfides();
	}

	// fix HG of CYS to relieve clashes of any newly converted CYS
	core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
	tf->push_back( core::pack::task::operation::TaskOperationOP( new core::pack::task::operation::OptCysHG() ) );
	core::pack::pack_rotamers( *pose_, *sfx, tf->create_task_and_apply_taskoperations( *pose_ ) );

	// safety
	pose_->energies().clear();
}

/// @brief the total number of free termini in the permutation
core::Size
StructureData::num_chains() const
{
	core::Size lower_count = 0;
	core::Size upper_count = 0;
	for ( SegmentMap::const_iterator r=segments_.begin(); r != segments_.end(); ++r ) {
		if ( r->second.has_free_lower_terminus() ) {
			++lower_count;
		}
		if ( r->second.has_free_upper_terminus() ) {
			++upper_count;
		}
	}

	if ( lower_count != upper_count ) {
		TR.Error << "Error in counting number of chains..." << lower_count << " vs. " << upper_count << std::endl;
		TR.Error << *this << std::endl;
	}
	debug_assert( lower_count == upper_count );
	return lower_count;
}

/// @brief returns the actual residue number of the given name and res #
core::Size
StructureData::pose_residue( std::string const & segment_name, core::Size const local_res ) const
{
	// first look in residues
	SegmentMap::const_iterator r = segments_.find( segment_name );
	if ( r == segments_.end() ) {
		if ( has_alias( segment_name ) ) {
			return alias_resnum( segment_name );
		}
	}
	if ( r == segments_.end() ) {
		std::stringstream err;
		err << "Can't resolve the segment name " << segment_name
			<< "into a valid segment in the permutation. Valid segments are: " << segments_ << std::endl;
		throw utility::excn::EXCN_BadInput( err.str() );
	}
	return r->second.resid(local_res);
}

/// @brief returns segments which have free lower termini in sequence order
utility::vector1< std::string >
StructureData::available_lower_termini() const
{
	utility::vector1< std::string > retval;
	for ( StringList::const_iterator r = segments_begin(); r != segments_end(); ++r ) {
		if ( segment( *r ).has_free_lower_terminus() ) {
			TR.Debug << *r << " has free lower terminus: " << segment( *r ).lower_segment() << std::endl;
			retval.push_back( *r );
		}
	}
	return retval;
}

/// @brief returns segments which have free upper termini in sequence order
utility::vector1< std::string >
StructureData::available_upper_termini() const
{
	utility::vector1< std::string > retval;
	for ( StringList::const_iterator r = segments_begin(); r != segments_end(); ++r ) {
		if ( segment( *r ).has_free_upper_terminus() ) {
			TR.Debug << *r << " has free upper terminus: " << segment( *r ).upper_segment() << std::endl;
			retval.push_back( *r );
		}
	}
	return retval;
}

/// @brief gives a quick yes-no answer as to whether it might be possible to connect these termini with an nres-residue loop
/// pose_ MUST BE SET if use_distance=true!!
bool
StructureData::are_connectable(
	std::string const & id1,
	std::string const & id2,
	core::Size const nres,
	bool const use_distance,
	bool const connection_performs_orientation,
	bool const allow_cyclic = false,
	core::Real const bond_dist = 1.5,
	core::Real const max_dist_per_res = 3.8 ) const
{
	// if one id or the other is already connected to something, this isn't connectable
	if ( ! segment(id1).has_free_upper_terminus() ) {
		TR.Debug << id1 << " and " << id2 << " len " << nres << " not connectable due to no free c anchor." << std::endl;
		TR.Debug << " available upper termini= " << available_upper_termini() << std::endl;
		return false;
	}
	if ( ! segment(id2).has_free_lower_terminus() ) {
		TR.Debug << id1 << " and " << id2 << " len " << nres << " not connectable due to no free n anchor." << std::endl;
		TR.Debug << " available lower termini= " << available_lower_termini() << std::endl;
		return false;
	}

	// ensure the two segments are on different chains
	// connecting a segment to itself is illegal unless allow_cyclic is true
	std::pair< std::string, std::string > termini1 = termini( id1 );
	std::pair< std::string, std::string > termini2 = termini( id2 );
	if ( !allow_cyclic && ( termini1 == termini2 ) ) {
		TR.Debug << id1 << " and " << id2 << " len " << nres << " not connectable because it would create a cyclic peptide." << std::endl;
		return false;
	}

	// ensure the two segments aren't in the same movable jump if performing orientation
	if ( connection_performs_orientation &&
			( movable_group( id1 ) == movable_group( id2 ) ) ) {
		TR.Debug << id1 << " and " << id2 << " len " << nres << " not connectable because they are in the same movable group : " << movable_group( id1 ) << "." << std::endl;
		return false;
	}

	if ( !use_distance ) {
		return true;
	}

	if ( !pose() ) {
		std::stringstream err;
		err << id() << ": no pose is present in this StructureData object, but a pose is required for StructureData::are_connectable() when use_distance is true." << std::endl;
		throw utility::excn::EXCN_Msg_Exception( err.str() );
	}
	// if the distance is > the fully extended Ca-Ca distance of an nres residue insert, this connection is physically impossible
	core::Size const res1( segment(id1).stop() ), res2( segment(id2).start() );
	// if there isn't a connect atom, find its name and use that
	core::chemical::ResidueType const & rtype = pose()->residue(res1).residue_type_set()->name_map( pose()->residue(res1).name3() );
	std::string const & aname = rtype.atom_name( rtype.upper_connect_atom() );
	core::Size const atom1 = pose()->residue(res1).type().atom_index( aname );
	debug_assert( atom1 );
	core::chemical::ResidueType const & rtype2 = pose()->residue(res2).residue_type_set()->name_map( pose()->residue(res2).name3() );
	std::string const & aname2 = rtype2.atom_name( rtype2.lower_connect_atom() );
	core::Size const atom2 = pose()->residue(res2).type().atom_index( aname2 );
	debug_assert( atom2 );
	core::Real const dist( pose()->residue(res1).xyz( atom1 ).distance(
		pose()->residue(res2).xyz( atom2 ) ) );
	// max 3.8 angstroms per residue, plus ~1.5 angstroms for N-C bond
	debug_assert( bond_dist >= -0.0000001 );
	if ( dist > (max_dist_per_res*nres + bond_dist) ) {
		TR << id1 << " and " << id2 << " not connectable due to nres=" << nres << " distance=" << dist << std::endl;
		return false;
	}
	return true;
}

/// @brief names a set of residues -- start_resid MUST be an nterm_resi() for a segment
void
StructureData::add_segment(
	std::string const & id_val,
	StringList::iterator insert_before_pos,
	core::Size const segment_length,
	core::Size const local_safe_residue,
	core::Size const local_cutpoint,
	core::Size const movable_group,
	bool const is_loop,
	bool const nterm_included,
	bool const cterm_included,
	std::string const & lower_conn,
	std::string const & upper_conn,
	std::string const & ss,
	utility::vector1< std::string > const & abego )
{
	Segment res( segment_length, local_safe_residue, local_cutpoint, movable_group, is_loop, nterm_included, cterm_included, lower_conn, upper_conn, ss, abego );

	add_segment( id_val, res, insert_before_pos );
}

/// @brief adds a residues segment to the end of the list
void
StructureData::add_segment( std::string const & id_val, Segment const & resis )
{
	TR.Debug << "Adding " << NamedSegment( id_val, resis ) << " to end of list" << std::endl;
	add_segment( id_val, resis, segment_order_.end() );
}

/// @brief adds a residues segment -- final will be inserted before segment named "insert_segment"
/// if "after_seg" is empty, or not found in the list of segments, the new segment will be inserted at the end
void
StructureData::add_segment( std::string const & id_val, Segment const & resis, std::string const & insert_before_segment )
{
	add_segment( id_val, resis, find_segment( insert_before_segment ) );
}

/// @brief adds a residues segment -- final will be inserted before the iterator given as insert_pose
void
StructureData::add_segment(
	std::string const & id_val,
	Segment const & resis,
	StringList::iterator insert_pos,
	core::pose::PoseCOP residues )
{
	NamedSegment toadd( id_val, resis );
	if ( insert_pos != segment_order_.end() ) {
		TR.Debug << "Adding " << toadd << " to after " << *insert_pos << std::endl;
	} else {
		TR.Debug << "Adding " << toadd << " to end of list" << std::endl;
	}
	debug_assert( segment_order_.size() == segments_.size() );
	SegmentMap::iterator s = segments_.find( id_val );
	if ( s == segments_.end() ) {
		if ( !segments_.insert( toadd ).second ) {
			throw utility::excn::EXCN_Msg_Exception( "failed to insert segment " + id_val );
		}
		segment_order_.insert( insert_pos, id_val );
		debug_assert( segments_.size() == segment_order_.size() );
	} else {
		s->second = resis;
		debug_assert( std::count( segment_order_.begin(), segment_order_.end(), id_val ) == 1 );
	}
	TR.Debug << "Segment order is now " << segment_order_ << std::endl;
	if ( pose_ && residues ) {
		core::Size insert_resid = pose_->total_residue() + 1;
		if ( insert_pos != segment_order_.end() ) {
			insert_resid = segment( *insert_pos ).nterm_resi();
		}
		TR << "Inserting new residues before position " << insert_resid << std::endl;
		pose_->conformation().insert_conformation_by_jump(
			residues->conformation(),
			insert_resid,
			pose_->num_jump()+2,
			1,
			pose_->num_jump()+1 );
	}
	update_numbering();
}

void
StructureData::add_segment(
	std::string const & id_val,
	Segment const & resid,
	StringList::iterator insert_pos )
{
	add_segment( id_val, resid, insert_pos, core::pose::PoseOP() );
}

/// @brief finds all segments that are loops
utility::vector1< std::string >
StructureData::loops() const
{
	utility::vector1< std::string > looplist;
	for ( SegmentMap::const_iterator s = segments_.begin(); s != segments_.end(); ++s ) {
		if ( s->second.is_loop ) {
			looplist.push_back( s->first );
		}
	}
	return looplist;
}

/// @brief returns n and c terminal segments of the chain which includes res
std::pair< std::string, std::string >
StructureData::termini( std::string const & seg ) const
{
	// go forward
	std::set< std::string > visited;
	std::string lower_segment = find_lower_terminus( visited, seg );
	TR.Debug << "Lower segment is " << lower_segment << " visited=" << visited << std::endl;
	visited.clear();
	std::string upper_segment = find_upper_terminus( visited, seg );
	TR.Debug << "Upper segment is " << upper_segment << " visited=" << visited << std::endl;
	return std::make_pair( lower_segment, upper_segment );
}

/// @brief performs dfs in lower direction looking for termini
std::string
StructureData::find_lower_terminus( std::set< std::string > & visited, std::string const & seg ) const
{
	// if we've already seen this segment, we are cyclic
	if ( visited.find(seg) != visited.end() ) {
		return "";
	}
	visited.insert(seg);
	if ( segment(seg).has_free_lower_terminus() ) {
		return seg;
	}
	debug_assert( segment(seg).lower_segment() != "" );
	return find_lower_terminus( visited, segment(seg).lower_segment() );
}

/// @brief performs dfs in upper direction looking for termini
std::string
StructureData::find_upper_terminus( std::set< std::string > & visited, std::string const & seg ) const
{
	// if we've already seen this segment, we are cyclic
	if ( visited.find(seg) != visited.end() ) {
		return "";
	}
	visited.insert(seg);
	if ( segment(seg).has_free_upper_terminus() ) {
		return seg;
	}
	debug_assert( segment(seg).upper_segment() != "" );
	return find_upper_terminus( visited, segment(seg).upper_segment() );
}

/// @brief sets an "alias" for a particular residue inside a segment which allows for it to be easily accessed
void
StructureData::set_resnum_alias(
	std::string const & alias_name,
	std::string const & segment_name,
	core::Size const resi )
{
	if ( has_segment( segment_name ) ) {
		debug_assert( resi );
		debug_assert( resi <= segment(segment_name).length() );
		aliases_[ alias_name ] = Alias( segment_name, resi );
	} else if ( has_alias( segment_name ) ) {
		aliases_[ alias_name ] = aliases_[ segment_name ];
	} else {
		std::stringstream ss;
		ss << "Segment named " << segment_name
			<< " does not exist in the permutation as an alias or residue segment. Perm="
			<< *this << std::endl;
		throw utility::excn::EXCN_BadInput( ss.str() );
	}
	save_into_pose();
}

/// @brief sets an "alias" for a particular residue which allows for it to be easily accessed
void
StructureData::set_resnum_alias(
	std::string const & alias_name,
	core::Size const resi )
{
	debug_assert( resi );
	debug_assert( resi <= pose_length() );
	std::string const segmentname = segment_name( resi );
	core::Size const localresi = resi - segment( segmentname ).start() + 1;
	TR << "Going to set alias with segment = " << segmentname << " Residue = " << localresi << std::endl;
	set_resnum_alias( alias_name, segmentname, localresi );
}

/// @brief given a residue alias, returns a pose residue number
core::Size
StructureData::alias_resnum( std::string const & alias ) const
{
	std::map< std::string, Alias >::const_iterator alias_it = aliases_.find( alias );
	if ( alias_it == aliases_.end() ) {
		TR.Warning << " alias " << alias << " not found!!  returning 0" << std::endl;
	}

	debug_assert( has_segment( alias_it->second.first ) );
	Segment const & resis = segment( alias_it->second.first );
	return resis.resid( alias_it->second.second );
}

/// @brief renames a residue segment and updates all connections. If the prefix is already there for a given segment name, it does nothing
void
StructureData::add_prefix_to_segments( std::string const & prefix )
{
	add_prefix_to_segments( prefix, PARENT_DELIMETER );
}

/// @brief renames a residue segment and updates all connections. If the prefix is already there for a given segment name, it does nothing
void
StructureData::add_prefix_to_segments( std::string const & prefix, char const delimeter )
{
	std::string const prepend_str = prefix + delimeter;
	SegmentMap newmap;
	for ( SegmentMap::iterator r = segments_.begin(); r != segments_.end(); ++r ) {
		if ( ! r->second.has_free_lower_terminus() ) {
			if ( !boost::starts_with( r->second.lower_segment(), prepend_str ) ) {
				r->second.set_lower_segment( prepend_str + r->second.lower_segment() );
			}
		}
		if ( ! r->second.has_free_upper_terminus() ) {
			if ( !boost::starts_with( r->second.upper_segment(), prepend_str ) ) {
				r->second.set_upper_segment( prepend_str + r->second.upper_segment() );
			}
		}
		std::string newname;
		if ( boost::starts_with( r->first, prepend_str ) ) {
			newname = r->first;
		} else {
			newname = prepend_str + r->first;
		}
		newmap[newname] = r->second;
		if ( r->first != newname ) {
			StringList::iterator c = find_segment( r->first );
			debug_assert( c != segment_order_.end() );
			*c = newname;
		}
	}
	// fix aliases too
	for ( std::map< std::string, Alias >::iterator a=aliases_.begin(); a!=aliases_.end(); ++a ) {
		if ( !boost::starts_with( a->second.first, prepend_str ) ) {
			a->second.first = prepend_str + a->second.first;
		}
	}
	segments_ = newmap;
	save_into_pose();
}

/// @brief renames a residue segment and updates all connections
void
StructureData::rename_segment( std::string const & old_name, std::string const & new_name )
{
	// find and rename the original segment
	SegmentMap::iterator r = segments_.find( old_name );
	debug_assert( r != segments_.end() );
	Segment resis = r->second;
	segments_.erase( r );
	segments_[new_name] = resis;

	StringList::iterator c = find_segment( old_name );
	debug_assert( c != segment_order_.end() );
	*c = new_name;

	for ( SegmentMap::iterator r2 = segments_.begin(); r2 != segments_.end(); ++r2 ) {
		if ( r2->second.lower_segment() == old_name ) {
			r2->second.set_lower_segment( new_name );
		}
		if ( r2->second.upper_segment() == old_name ) {
			r2->second.set_upper_segment( new_name );
		}
	}

	// fix aliases too
	for ( std::map< std::string, Alias >::iterator a=aliases_.begin(); a!=aliases_.end(); ++a ) {
		if ( a->second.first == old_name ) a->second.first = new_name;
	}
	save_into_pose();
}

/// @brief true if this permutation contains a residue segment named seg
bool
StructureData::has_segment( std::string const & seg ) const
{
	return ( ( segments_.find(seg) != segments_.end() ) ||
		( segments_.find( id() + PARENT_DELIMETER + seg ) != segments_.end() ) );
}

/// @brief start of segments list
StringList::const_iterator
StructureData::segments_begin() const
{
	return segment_order_.begin();
}

/// @brief end of segment list
StringList::const_iterator
StructureData::segments_end() const
{
	return segment_order_.end();
}


core::Size
StructureData::choose_new_movable_group() const
{
	std::set< core::Size > const mgs = movable_groups();
	core::Size mg = 1;
	while ( mgs.find( mg ) != mgs.end() ) {
		++mg;
	}
	return mg;
}

/// @brief merge all data and segments from "other" into this StructureData
void
StructureData::merge( StructureData const & other )
{
	merge( other, StringList( other.segments_begin(), other.segments_end() ) );
}


/// @brief merge given data and segments from "other" into this StructureData
void
StructureData::merge( StructureData const & other, StringList const & segments )
{
	// this should only be called on multi-chain permutations
	runtime_assert( is_multi() );

	core::Size const movable_group_new = choose_new_movable_group();
	for ( StringList::const_iterator s=segments.begin(); s!=segments.end(); ++s ) {
		Segment newseg = other.segment( *s );
		newseg.movable_group = movable_group_new;
		core::pose::PoseOP newpose;
		if ( other.pose() ) {
			newpose = core::pose::PoseOP(
				new core::pose::Pose(
				*(other.pose()),
				other.segment( *s ).nterm_resi(),
				other.segment( *s ).cterm_resi() ) );
		}
		add_segment( *s, newseg, segment_order_.end(), newpose );
		TR.Debug << "Added " << NamedSegment( *s, newseg ) << std::endl;
	}

	copy_data( other );
}

void
StructureData::merge_before( StructureData const & other, std::string const & position, StringList const & segments )
{
	// this should only be called on multi-chain permutations
	runtime_assert( is_multi() );

	core::Size const movable_group_new = choose_new_movable_group();
	for ( StringList::const_iterator c=segments.begin(); c!=segments.end(); ++c ) {
		Segment newseg = other.segment( *c );
		newseg.movable_group = movable_group_new;
		core::pose::PoseOP newpose;
		if ( other.pose() ) {
			newpose = core::pose::PoseOP(
				new core::pose::Pose(
				*(other.pose()),
				other.segment( *c ).nterm_resi(),
				other.segment( *c ).cterm_resi() ) );
		}
		add_segment( *c, newseg, find_segment( position ), newpose );
		TR.Debug << "Added " << NamedSegment( *c, newseg ) << std::endl;
	}

	copy_data( other );
}

void
StructureData::merge_before( StructureData const & other, std::string const & position )
{
	merge_before( other, position, StringList( other.segments_begin(), other.segments_end() ) );
}

/// @brief computes and returns a set of segments which are in the given movable group
utility::vector1< std::string >
StructureData::segments_in_movable_group( core::Size const group ) const
{
	utility::vector1< std::string > segments;
	for ( SegmentMap::const_iterator r = segments_.begin(); r != segments_.end(); ++r ) {
		if ( r->second.movable_group == group ) {
			segments.push_back( r->first );
		}
	}
	return segments;
}

/// @brief computes and returns a set of movable groups
std::set< core::Size >
StructureData::movable_groups() const
{
	std::set< core::Size > groups;
	for ( SegmentMap::const_iterator res_it = segments_.begin(); res_it != segments_.end(); ++res_it ) {
		groups.insert( res_it->second.movable_group );
	}
	return groups;
}

/// @brief moves around jumps so they are located in the center of rigid residue groups
void
StructureData::move_jumps_to_safety()
{
	if ( !pose() ) return;

	core::kinematics::FoldTree ft = pose()->fold_tree();

	for ( core::Size i = 1; i <= ft.num_jump(); ++i ) {
		core::kinematics::Edge e = ft.jump_edge(i);
		TR << "Jump " << e << " connects " << segment_name( e.start() ) << " to " << segment_name( e.stop() ) << std::endl;
		std::string const upname = segment_name( e.start() );
		e.start() = segment( upname ).safe();
		std::string const downname = segment_name( e.stop() );
		e.stop() = segment( downname ).safe();
		TR.Debug << "FT: " << ft << std::endl;
		TR.Debug << "Sliding jump " << e << std::endl;
		ft.slide_jump( e.label(), e.start(), e.stop() );
	}
	if ( !ft.check_fold_tree() ) {
		TR << "FAILING fold tree check:" << *this << std::endl;
		debug_assert( ft.check_fold_tree() );
	}
	pose_->fold_tree( ft );
}

/// @brief moves around jumps so that movable groups will all move together during folding
void
StructureData::consolidate_movable_groups( utility::vector1< std::string > const & root_segments )
{
	if ( pose_ ) {
		consolidate_movable_groups( pose_, root_segments );
	}
}

/// @brief moves around jumps so that movable groups will all move together during folding
void
StructureData::consolidate_movable_groups(
	core::pose::PoseOP pose,
	utility::vector1< std::string > const & root_segments )
{
	debug_assert( pose );

	// make fold graph
	FoldGraph fg( *this, pose );
	core::kinematics::FoldTree ft = fg.fold_tree( *this, root_segments );
	if ( ft.nres() != pose->total_residue() ) {
		TR.Error << "FT nres " << ft.nres() << " != pose nres " << pose->total_residue() << " FT: " << ft << std::endl;
		throw utility::excn::EXCN_Msg_Exception( "Bad nres" );
	}
	pose->fold_tree( ft );
	TR.Debug << "Created fold tree: " << pose->fold_tree() << std::endl;
	TR.Debug << "FT Perm = " << *this << std::endl;
}

/// @brief declares a covalent bond between the specified atoms
void
StructureData::declare_covalent_bond(
	std::string const & seg1, core::Size const res1, std::string const & atom1,
	std::string const & seg2, core::Size const res2, std::string const & atom2 )
{
	declare_covalent_bond( pose_residue( seg1, res1 ), atom1, pose_residue( seg2, res2 ), atom2 );
}

void
StructureData::add_covalent_bond(
	core::Size const res1, std::string const & atom1,
	core::Size const res2, std::string const & atom2 )
{
	std::string const & seg1 = segment_name( res1 );
	core::Size const localres1 = res1 - segment( seg1 ).start() + 1;
	std::string const & seg2 = segment_name( res2 );
	core::Size const localres2 = res2 - segment( seg2 ).start() + 1;
	if ( localres1 && localres2 ) {
		add_covalent_bond( seg1, localres1, atom1, seg2, localres2, atom2 );
	} else {
		TR << "Warning: connection between residues " << res1 << " and " << res2 << " may be lost since one/both of the residues are present as \"padding\" residues." << std::endl;
	}
}

void
StructureData::add_covalent_bond(
	std::string const & seg1, core::Size const res1, std::string const & atom1,
	std::string const & seg2, core::Size const res2, std::string const & atom2 )
{
	add_covalent_bond( BondInfo( seg1, seg2, res1, res2, atom1, atom2 ) );
}

void
StructureData::add_covalent_bond( BondInfo const & bi )
{
	// look for this bond, stop if it's found
	for ( utility::vector1< BondInfo >::const_iterator bi_it=covalent_bonds_.begin(); bi_it != covalent_bonds_.end(); ++bi_it ) {
		if ( bi == *bi_it ) {
			TR.Debug << "Skipping adding existing bond info " << *bi_it << std::endl;
			return;
		}
	}
	covalent_bonds_.push_back( bi );
	save_into_pose();
}

void
StructureData::update_covalent_bonds_in_pose()
{
	if ( !pose_ ) return;
	for ( utility::vector1< BondInfo >::const_iterator p = covalent_bonds_.begin(); p != covalent_bonds_.end(); ++p ) {
		TR.Debug << "Updating covalent bond " << p->seg1 << "-->" << p->seg2 << std::endl;
		declare_covalent_bond_in_pose( pose_residue( p->seg1, p->res1 ), p->atom1, pose_residue( p->seg2, p->res2 ), p->atom2 );
	}
}

/// @brief declares a covalent bond using pose residues
void
StructureData::declare_covalent_bond(
	core::Size const res1, std::string const & atom1,
	core::Size const res2, std::string const & atom2 )
{
	TR.Debug << "Creating covalent bond between " << res1 << " and " << res2 << std::endl;
	if ( ( res1 != res2 + 1 ) && ( res1 + 1 != res2 ) ) {
		add_covalent_bond( res1, atom1, res2, atom2 );
	}
	declare_covalent_bond_in_pose( res1, atom1, res2, atom2 );
}

void
StructureData::declare_covalent_bond_in_pose(
	core::Size const res1, std::string const & atom1,
	core::Size const res2, std::string const & atom2 )
{
	if ( !pose_ ) return;
	if ( pose_->residue( res1 ).has( atom1 ) && pose_->residue( res2 ).has( atom2 ) ) {
		pose_->conformation().declare_chemical_bond( res1, atom1, res2, atom2 );
		//Rebuild the connection atoms:
		for ( core::Size i = 1; i <= pose_->conformation().residue(res1).n_residue_connections(); ++i ) {
			pose_->conformation().rebuild_residue_connection_dependent_atoms( res1, i );
		}
		for ( core::Size i = 1; i <= pose_->conformation().residue(res2).n_residue_connections(); ++i ) {
			pose_->conformation().rebuild_residue_connection_dependent_atoms( res2, i );
		}
	} else {
		TR.Debug << "Skipping covalent bond declaration between " << res1 << ":" << atom1
			<< " and " << res2 << ":" << atom2 << " due to missing atoms. The residue identity may have changed." << std::endl;
	}
}

/// @brief declares a covalent bond between the c-terminal atom of segment 1 and the n-terminal atom of segment 2
void
StructureData::declare_polymer_bond( std::string const & segment1, std::string const & segment2 )
{
	core::Size const pos1 = segment(segment1).cterm_resi();
	core::Size const pos2 = segment(segment2).nterm_resi();
	debug_assert( pose_ );
	debug_assert( pose_->residue(pos1).is_polymer() );
	debug_assert( pose_->residue(pos2).is_polymer() );
	debug_assert( ! pose_->residue(pos1).is_polymer_bonded( pos2 ) );

	declare_covalent_bond( pos1, pose_->residue(pos1).atom_name( pose_->residue(pos1).upper_connect_atom() ),
		pos2, pose_->residue(pos2).atom_name( pose_->residue(pos2).lower_connect_atom() ) );

	debug_assert( pose_->residue(pos1).is_polymer_bonded( pos2 ) );
}

/// @brief marks the resi-th residue of seg as a cutpoint
void
StructureData::set_cutpoint( std::string const & seg, core::Size const resi )
{
	SegmentMap::iterator r = segments_.find( seg );
	debug_assert( r != segments_.end() );
	if ( r != segments_.end() ) {
		r->second.set_cutpoint( resi );
	}
	save_into_pose();
}

/// @brief connects the given chains together, doesn't update anything -- don't call this on its own unless you know what you're doing.
void
StructureData::connect_segments(
	std::string const & segment1_c,
	std::string const & segment2_n )
{
	// connected segments
	if ( has_free_upper_terminus( segment1_c ) && has_free_lower_terminus( segment2_n ) ) {
		mark_connected( segment1_c, segment2_n );
	}
	debug_assert( segment(segment1_c).upper_segment() == segment2_n );
	debug_assert( segment(segment2_n).lower_segment() == segment1_c );

	// pose manipulation
	if ( pose_ ) {
		remove_upper_terminus_variant_type( segment(segment1_c).cterm_resi() );
		remove_lower_terminus_variant_type( segment(segment2_n).nterm_resi() );
		if ( ! polymer_bond_exists( segment1_c, segment2_n ) ) {
			declare_polymer_bond( segment1_c, segment2_n );
		}
		chains_from_termini();
	}
	TR.Debug << "Marked " << segment1_c << " and " << segment2_n << " as connected." << std::endl;
	save_into_pose();
}

void
StructureData::mark_disconnected(
	std::string const & seg1,
	std::string const & seg2 )
{
	SegmentMap::iterator s1 = segments_.find( seg1 );
	SegmentMap::iterator s2 = segments_.find( seg2 );
	debug_assert( s1 != segments_.end() );
	debug_assert( s2 != segments_.end() );

	s1->second.set_upper_segment( "" );
	s2->second.set_lower_segment( "" );

	TR.Debug << "Marked " << s1->first << " and " << s2->first << " as disconnected." << std::endl;
	save_into_pose();
}

/// @brief disconnects the given chains doesn't update anything -- don't call this on its own unless you know what you're doing.
void
StructureData::disconnect_segments(
	std::string const & segment1_c,
	std::string const & segment2_n )
{
	mark_disconnected( segment1_c, segment2_n );

	// pose manipulation
	if ( pose_ ) {
		/*
		if ( polymer_bond_exists( segment1_c, segment2_n ) )
		pas
		*/

		add_upper_terminus_variant_type( segment( segment1_c ).cterm_resi() );
		add_lower_terminus_variant_type( segment( segment2_n ).nterm_resi() );
		chains_from_termini();
	}
}

/// @brief merges two segments into one that has the name new_name. They must be next to each other in sequence.
void
StructureData::merge_segments(
	std::string const & segment1,
	std::string const & segment2,
	std::string const & new_name )
{
	Segment c1 = segment(segment1);
	Segment c2 = segment(segment2);
	utility::vector1< std::string > c1_grp = segments_in_movable_group( c1.movable_group );
	utility::vector1< std::string > c2_grp = segments_in_movable_group( c2.movable_group );
	core::Size new_group = c1.movable_group;
	//core::Size del_group = c2.movable_group;
	debug_assert( ( c1_grp.size() <= 1 ) || ( c2_grp.size() <= 1 ) );
	if ( c2_grp.size() > c1_grp.size() ) {
		new_group = c2.movable_group;
		//del_group = c1.movable_group;
	}

	// remove terminal variants if necessary
	if ( pose_ ) {
		remove_upper_terminus_variant_type( segment(segment1).cterm_resi() );
		remove_lower_terminus_variant_type( segment(segment2).nterm_resi() );
		chains_from_termini();
	}
	debug_assert( c1.cterm_resi()+1 == c2.nterm_resi() );

	std::string const newss = segment(segment1).ss() + segment(segment2).ss();
	utility::vector1< std::string > const & c1_abego = segment(segment1).abego();
	utility::vector1< std::string > const & c2_abego = segment(segment2).abego();
	utility::vector1< std::string > new_abego( c1_abego.size() + c2_abego.size() );
	std::copy( c1_abego.begin(), c1_abego.end(), new_abego.begin() );
	std::copy( c2_abego.begin(), c2_abego.end(), new_abego.begin() + c1_abego.size() );
	debug_assert( newss.size() == c1.length() + c2.length() );
	debug_assert( new_abego.size() == c1.length() + c2.length() );

	// handle merging of residues segments
	std::string const low_seg = segment( segment1 ).lower_segment();
	std::string const upp_seg = segment( segment2 ).upper_segment();
	Segment resis(
		c1.length() + c2.length(), // length
		( ( c1.start() + c2.stop() ) / 2 ) - c1.nterm_resi() + 1, // safe residue
		0, // cut residue
		new_group,
		( c1.is_loop && c2.is_loop ),
		c1.nterm_included(),
		c2.cterm_included(),
		"",
		"",
		newss,
		new_abego );

	TR << "Removing old segments " << NamedSegment( *segments_.find( segment1 ) ) << " and " << NamedSegment( *segments_.find( segment2 ) )
		<< " and adding new merged segment " << NamedSegment( new_name, resis ) << std::endl;
	StringList::const_iterator after_me = find_segment( segment2 );
	if ( after_me != segment_order_.end() ) ++after_me;
	std::string const after_str = *after_me;
	delete_segment_nopose( segment2, segments_.find( segment2 ) );
	core::Size const num_groups = movable_groups().size();
	delete_segment_nopose( segment1, segments_.find( segment1 ) );
	if ( movable_groups().size() < num_groups ) {
		resis.movable_group = movable_groups().size() + 1;
	}
	add_segment( new_name, resis, find_segment( after_str ) );
	if ( !low_seg.empty() ) {
		connect_segments( low_seg, new_name );
	}
	if ( !upp_seg.empty() ) {
		connect_segments( new_name, upp_seg );
	}
	/*
	// erase old residues objects
	segments_.erase( segment1 );
	segments_.erase( segment2 );

	// erase old strings in order
	StringList::iterator pos = find_segment( segment1 );
	debug_assert( pos != segment_order_.end() );
	pos = segment_order_.erase( pos );
	debug_assert( pos != segment_order_.end() );
	debug_assert( *pos == segment2 );
	pos = segment_order_.erase( pos );

	//rename within residues
	for ( SegmentMap::iterator r=segments_.begin(); r!=segments_.end(); ++r ) {
	if ( ( r->second.upper_segment() == segment2 ) || ( r->second.upper_segment() == segment1 ) ) {
	r->second.set_upper_segment( new_name );
	}
	if ( ( r->second.lower_segment() == segment2 ) || ( r->second.lower_segment() == segment1 ) ) {
	TR << "Renaming lower segment of " << NamedSegment(*r) << " from " << r->second.lower_segment() << " to " << new_name << std::endl;
	r->second.set_lower_segment( new_name );
	}
	}
	std::set< core::Size > mgs = movable_groups();
	if ( mgs.find( mgs.size() + 1 ) != mgs.end() ) {
	renumber_movable_group( mgs.size()+1, del_group );
	}
	add_segment( new_name, resis, pos );
	*/
}

/// @brief sets the pose and does checks to ensure data is consistent
void
StructureData::set_pose( core::pose::Pose const & new_pose )
{
	set_pose( new_pose.clone() );
}

/// @brief checks pose vs. StructureData info
/// @throw EXCN_BadInput if things don't match
void
StructureData::check_pose( core::pose::Pose const & pose ) const
{
	if ( ! pose.total_residue() ) return;
	if ( pose.total_residue() != pose_length() ) {
		std::stringstream msg;
		msg << id() << ": pose length does not match StructureData.  Pose length = "
			<< pose.total_residue() << " SD length = " << pose_length() << std::endl;
		msg << *this << std::endl;
		//debug_assert( pose.total_residue() == pose_length() );
		throw EXCN_PoseInconsistent( msg.str() );
	}
	if ( pose.total_residue() != ss().size()  ) {
		std::stringstream msg;
		msg << id() << ": pose length does not match StructureData secstruct length. Pose length = "
			<< pose.total_residue() << " SS length = " << ss().size() << " SD: " << *this << std::endl;
		//debug_assert( pose.total_residue() == pose_length() );
		throw EXCN_PoseInconsistent( msg.str() );
	}
}

/// @brief sets the pose and does checks to ensure data is consistent
/// @throws EXCN_PoseInconsistent if the pose doesn't match up with SD
void
StructureData::set_pose( core::pose::PoseOP new_pose )
{
	if ( new_pose && new_pose->total_residue() ) {
		TR.Debug << id() << ": new_pose len=" << new_pose->total_residue() << " pose residue =" << pose_length() << " perm len=" << length() << std::endl;

		check_consistency( *new_pose );

		// put terminal variants onto chain endings
		utility::vector1< core::Size > const conformation_chainend = new_pose->conformation().chain_endings();
		utility::vector1< core::Size > chain_endings = conformation_chainend;
		chain_endings.push_back( new_pose->total_residue() );
		utility::vector1< core::Size > chain_starts;
		if ( new_pose->total_residue() ) {
			chain_starts.push_back( 1 );
		}
		for ( utility::vector1< core::Size >::const_iterator r = conformation_chainend.begin(); r != conformation_chainend.end(); ++r ) {
			if ( *r == new_pose->total_residue() ) continue;
			chain_starts.push_back( *r + 1 );
		}
		TR.Debug << "Chain endings: " << chain_starts << chain_endings << std::endl;

		for ( utility::vector1< core::Size >::const_iterator r = chain_starts.begin(); r != chain_starts.end(); ++r ) {
			debug_assert( *r > 0 );
			debug_assert( *r <= new_pose->total_residue() );
			if ( new_pose->residue( *r ).is_polymer() &&
					( ! new_pose->residue( *r ).is_lower_terminus() ) &&
					( ! new_pose->residue( *r ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) &&
					( ! new_pose->residue( *r ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) ) {
				core::pose::add_lower_terminus_type_to_pose_residue( *new_pose, *r );
				TR.Debug << "Added lower terminus type to residue " << *r << std::endl;
			} else {
				TR.Debug << "Residue " << new_pose->residue( *r ).name() << *r << " is already a lower terminus or non-polymer." << std::endl;
			}
		}
		for ( utility::vector1< core::Size >::const_iterator r = chain_endings.begin(); r != chain_endings.end(); ++r ) {
			debug_assert( *r > 0 );
			debug_assert( *r <= new_pose->total_residue() );
			if ( new_pose->residue( *r ).is_polymer() &&
					( ! new_pose->residue( *r ).is_upper_terminus() ) &&
					( ! new_pose->residue( *r ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) &&
					( ! new_pose->residue( *r ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) ) {
				core::pose::add_upper_terminus_type_to_pose_residue( *new_pose, *r );
				TR.Debug << "Added upper terminus type to residue " << *r << std::endl;
			} else {
				TR.Debug << "Residue " << new_pose->residue( *r ).name() << *r << " is already a upper terminus or non-polymer." << std::endl;
			}
		}
		new_pose->conformation().chains_from_termini();
	}
	pose_ = new_pose;
	save_into_pose();
}

/// @brief based on the pose, determines the jump number for the jump pointing to the segment specified
int
StructureData::find_jump( std::string const & seg ) const
{
	debug_assert( pose_ );
	int const j = find_jump_rec( pose_->fold_tree(), segment( seg ).safe() );
	debug_assert( j >= 0 );
	return j;
}

/// @brief returns n-terminal residue of the chain represented by given string
core::Size
StructureData::lower_anchor( std::string const & id_val ) const
{
	return segment(id_val).start();
}

/// @brief returns c-terminal residue of the chain represented by given string
core::Size
StructureData::upper_anchor( std::string const & id_val ) const
{
	return segment(id_val).stop();
}

/// @brief returns non-const residue range of the segment represented by given string
Segment &
StructureData::segment_nonconst( std::string const & id_val )
{
	SegmentMap::iterator it = segments_.find( id_val );
	if ( it == segments_.end() ) {
		it = segments_.find( id() + PARENT_DELIMETER + id_val );
	}
	if ( it == segments_.end() ) {
		std::stringstream err;
		err << id() << ": Segment not found in residue lists! ";
		err << "Search term is: " << id_val << "; Segment map is: " << segments_ << std::endl;
		debug_assert( false );
		throw utility::excn::EXCN_Msg_Exception( err.str() );
	}
	return it->second;
}

/// @brief returns residue range of the segment represented by given string
Segment const &
StructureData::segment( std::string const & id_val ) const
{
	SegmentMap::const_iterator it = segments_.find( id_val );
	if ( it == segments_.end() ) {
		it = segments_.find( id() + PARENT_DELIMETER + id_val );
	}
	if ( it == segments_.end() ) {
		std::stringstream err;
		err << id() << ": Segment not found in residue lists! ";
		err << "Search term is: " << id_val << "; Segment map is: " << segments_ << std::endl;
		debug_assert( false );
		throw utility::excn::EXCN_Msg_Exception( err.str() );
	}
	return it->second;
}

bool
StructureData::has_data_int( std::string const & segment_id, std::string const & data_name ) const
{
	return has_data_int( segment_id + DATA_DELIMETER + data_name );
}

/// @brief check for real number data
bool
StructureData::has_data_real( std::string const & segment_id, std::string const & data_name ) const
{
	return has_data_real( segment_id + DATA_DELIMETER + data_name );
}

/// @brief gets real number data
bool
StructureData::has_data_str( std::string const & segment_id, std::string const & data_name ) const
{
	return has_data_str( segment_id + DATA_DELIMETER + data_name );
}

/// @brief sets integer number data
void
StructureData::set_data_int( std::string const & segment_id, std::string const & data_name, int const val )
{
	set_data_int( segment_id + DATA_DELIMETER + data_name, val );
}

/// @brief sets integer number data
void
StructureData::set_data_int( std::string const & data_name, int const val )
{
	std::map< std::string, int >::iterator it = data_int_.find(data_name);
	if ( it == data_int_.end() ) {
		data_int_[data_name] = val;
	} else {
		it->second = val;
	}
	save_into_pose();
}

/// @brief sets real number data
void
StructureData::set_data_real( std::string const & segment_id, std::string const & data_name, core::Real const val )
{
	set_data_real( segment_id + DATA_DELIMETER + data_name, val );
}

/// @brief sets real number data
void
StructureData::set_data_real( std::string const & data_name, core::Real const val )
{
	std::map< std::string, core::Real >::iterator it = data_real_.find(data_name);
	if ( it == data_real_.end() ) {
		data_real_[data_name] = val;
	} else {
		it->second = val;
	}
	save_into_pose();
}

/// @brief sets string data
void
StructureData::set_data_str( std::string const & segment_id, std::string const & data_name, std::string const & val )
{
	return set_data_str( segment_id + DATA_DELIMETER + data_name, val );
}

/// @brief sets string data
void
StructureData::set_data_str( std::string const & data_name, std::string const & val )
{
	std::map< std::string, std::string >::iterator it = data_str_.find(data_name);
	if ( it == data_str_.end() ) {
		data_str_[data_name] = val;
	} else {
		it->second = val;
	}
	save_into_pose();
}

/// @brief gets int number data
int
StructureData::get_data_int( std::string const & segment_id, std::string const & data_name ) const
{
	return get_data_int( segment_id + DATA_DELIMETER + data_name );
}

/// @brief gets real number data
int
StructureData::get_data_int( std::string const & data_name ) const
{
	std::map< std::string, int >::const_iterator it = data_int_.find(data_name);
	if ( it == data_int_.end() ) {
		std::stringstream err( "In StructureData: " );
		err << id() << ": " << data_name << " not found in data map." << std::endl;
		for ( it = data_int_.begin(); it != data_int_.end(); ++it ) {
			err << it->first << " : " << it->second << std::endl;
		}
		throw utility::excn::EXCN_Msg_Exception( err.str() );
	}
	return it->second;
}

/// @brief gets real number data
core::Real
StructureData::get_data_real( std::string const & segment_id, std::string const & data_name ) const
{
	return get_data_real( segment_id + DATA_DELIMETER + data_name );
}

/// @brief gets real number data
core::Real
StructureData::get_data_real( std::string const & data_name ) const
{
	std::map< std::string, core::Real >::const_iterator it = data_real_.find(data_name);
	if ( it == data_real_.end() ) {
		std::stringstream err( "In StructureData: " );
		err << id() << ": " << data_name << " not found in data map." << std::endl;
		for ( it = data_real_.begin(); it != data_real_.end(); ++it ) {
			err << it->first << " : " << it->second << std::endl;
		}
		throw utility::excn::EXCN_Msg_Exception( err.str() );
	}
	return it->second;
}

/// @brief gets string data
std::string const &
StructureData::get_data_str( std::string const & segment_id, std::string const & data_name ) const
{
	return get_data_str( segment_id + DATA_DELIMETER + data_name );
}

/// @brief gets string data
std::string const &
StructureData::get_data_str( std::string const & data_name ) const
{
	std::map< std::string, std::string >::const_iterator it = data_str_.find(data_name);
	if ( it == data_str_.end() ) {
		std::stringstream err( "In StructureData: " );
		err << id() << ": " << data_name << " not found in data map." << std::endl;
		for ( it = data_str_.begin(); it != data_str_.end(); ++it ) {
			err << it->first << " : " << it->second << std::endl;
		}
		throw utility::excn::EXCN_Msg_Exception( err.str() );
	}
	return it->second;
}

/// @brief copies user data fields from one permutation to this one -- existing data is overridden
void
StructureData::copy_data( StructureData const & perm )
{
	copy_data( perm, true );
}

template< class T >
void
set_map_data(
	std::map< std::string, T > & map,
	std::string const & key,
	T const & value,
	bool const overwrite )
{
	if ( overwrite ) {
		map[ key ] = value;
	} else {
		if ( map.insert( std::make_pair( key, value ) ).second == false ) {
			TR.Debug << "skipping new map data with different value for key " << key
				<< " : " << map.find( key )->second << " exists vs. new value " << value << std::endl;
		}
	}
}

/// @brief copies user data fields from one permutation to this one -- gives the option to not overwrite data
void
StructureData::copy_data( StructureData const & perm, bool const overwrite )
{
	for ( std::map< std::string, core::Real >::const_iterator d=perm.data_real().begin(); d!=perm.data_real().end(); ++d ) {
		set_map_data( data_real_, d->first, d->second, overwrite );
	}
	for ( std::map< std::string, int >::const_iterator d=perm.data_int().begin(); d!=perm.data_int().end(); ++d ) {
		set_map_data( data_int_, d->first, d->second, overwrite );
	}
	for ( std::map< std::string, std::string >::const_iterator d=perm.data_str().begin(); d!=perm.data_str().end(); ++d ) {
		set_map_data( data_str_, d->first, d->second, overwrite );
	}
	for ( std::map< std::string, Alias >::const_iterator a=perm.aliases().begin(); a!=perm.aliases().end(); ++a ) {
		set_map_data( aliases_, a->first, a->second, overwrite );
	}
	for ( utility::vector1< BondInfo >::const_iterator bi=perm.covalent_bonds_begin(); bi!=perm.covalent_bonds_end(); ++bi ) {
		add_covalent_bond( *bi );
	}
	save_into_pose();
}
/*
void
StructureData::copy_data( StructureData const & perm, bool const do_rename )
{
for ( std::map< std::string, core::Real >::const_iterator d=perm.data_real().begin(), end=perm.data_real().end(); d != end; ++d ) {
std::pair< std::string, core::Real > p = *d;
if ( do_rename && ( !boost::starts_with( p.first, perm.id() ) ) ) {
p.first = perm.id() + PARENT_DELIMETER + p.first;
}
data_real_[ p.first ] = p.second;
if ( data_real_.insert( p ).second == false ) {
if ( p.second != d->second ) {
TR.Debug << "overwriting data with different value " << p.first << " : " << p.second << " vs " << d->first << " : " << d->second << std::endl;
data_real_[ p.first ] = p.second;
}
}
}
for ( std::map< std::string, int >::const_iterator d=perm.data_int().begin(), end=perm.data_int().end(); d != end; ++d ) {
std::pair< std::string, int > p = *d;
if ( do_rename && ( !boost::starts_with( p.first, perm.id() ) ) ) {
TR << "Renaming " << p.first << " with id " << perm.id() << " to " << perm.id() + PARENT_DELIMETER + p.first << std::endl;
p.first = perm.id() + PARENT_DELIMETER + p.first;
}
data_int_[ p.first ] = p.second;
if ( data_int_.insert( p ).second == false ) {
if ( p.second != d->second ) {
TR.Debug << "overwriting data with different value " << p.first << " : " << p.second << " vs " << d->first << " : " << d->second << std::endl;
data_int_[ p.first ] = p.second;
}
}
}
for ( std::map< std::string, std::string >::const_iterator d=perm.data_str().begin(), end=perm.data_str().end(); d != end; ++d ) {
std::pair< std::string, std::string > p = *d;
if ( do_rename && ( !boost::starts_with( p.first, perm.id() ) ) ) {
p.first = perm.id() + PARENT_DELIMETER + p.first;
}
data_str_[ p.first ] = p.second;
if ( data_str_.insert( p ).second == false ) {
if ( p.second != d->second ) {
TR.Debug << "overwriting data with different value " << p.first << " : " << p.second << " vs " << d->first << " : " << d->second << std::endl;
data_str_[ p.first ] = p.second;
}
}
}
for ( std::map< std::string, std::pair< std::string, core::Size > >::const_iterator a=perm.aliases().begin(), enda=perm.aliases().end(); a != enda; ++ a ) {
std::pair< std::string, std::pair< std::string, core::Size > > p = *a;
if ( do_rename && ( !boost::starts_with( p.second.first, perm.id() ) ) ) {
p.second.first = perm.id() + PARENT_DELIMETER + p.second.first;
}
aliases_[ p.first ] = p.second;
if ( aliases_.insert( p ).second == false ) {
if ( p.second != a->second ) {
TR.Debug << "overwriting data with different value " << p.first << " : " << p.second.first << "__" << p.second.second << " vs " << a->first << " : " << a->second.first << "__" << a->second.second << std::endl;
aliases_[ p.first ] = p.second;
}
}
}
for ( utility::vector1< BondInfo >::const_iterator bi=perm.covalent_bonds_begin(); bi!=perm.covalent_bonds_end(); ++bi ) {
add_covalent_bond( *bi );
}
save_into_pose();
} */

/// @brief just return the pose (doesn't build pose if it doesn't exist)
core::pose::PoseCOP
StructureData::pose() const
{
	return pose_;
}

/// @brief returns segment which includes residue number res
std::string const &
StructureData::segment_name( core::Size const res ) const
{
	for ( SegmentMap::const_iterator it=segments_.begin(); it != segments_.end(); ++it ) {
		if ( ( res >= it->second.nterm_resi() ) && ( res <= it->second.cterm_resi() ) ) {
			return it->first;
		}
	}
	std::stringstream err;
	err << " Residue " << res << " was not found in the residues map!" << std::endl;
	throw utility::excn::EXCN_Msg_Exception( err.str() );
}

/// @brief finds and returns the loop residues
utility::vector1< bool >
StructureData::loop_residues() const
{
	utility::vector1< bool > is_loop( pose_length(), false );
	for ( SegmentMap::const_iterator it=segments_.begin(); it != segments_.end(); ++it ) {
		if ( it->second.is_loop ) {
			TR << "Loop Range: " << it->second.start() << " -> " << it->second.stop() << std::endl;
			for ( core::Size i = it->second.start(); i <= it->second.stop(); ++i ) {
				debug_assert( i > 0 );
				debug_assert( i <= is_loop.size() );
				is_loop[i] = true;
			}
		}
	}

	return is_loop;
}

/// @brief updates numbering based on the saved order of Segment objects
void
StructureData::update_numbering()
{
	debug_assert( segment_order_.size() == segments_.size() );
	core::Size cur_num = 1;
	core::Size non_dummy_count = 0;
	std::string new_ss = "";
	utility::vector1< std::string > new_abego;

	for ( StringList::const_iterator c = segments_begin(); c != segments_end(); ++c ) {
		SegmentMap::iterator r = segments_.find( *c );
		debug_assert( r != segments_.end() );
		r->second.set_pose_start( cur_num );
		cur_num += r->second.length();
		non_dummy_count += ( r->second.stop() - r->second.start() + 1 );
		new_ss += r->second.ss();
		for ( utility::vector1< std::string >::const_iterator a = r->second.abego().begin();
				a != r->second.abego().end(); ++a ) {
			new_abego.push_back( *a );
		}
	}
	pose_length_= cur_num-1;
	length_ = non_dummy_count;
	if ( new_ss.size() != pose_length_ ) {
		std::stringstream err;
		err << id() << ": StructureData pose size doesn't match secondary structure string size.  this is probably an internal bug that needs to be fixed." << std::endl;
		err << *this << std::endl;
		err << "new ss= " << new_ss << std::endl;
		throw utility::excn::EXCN_BadInput( err.str() );
	}
	debug_assert( new_ss.size() == pose_length_ );
	ss_ = new_ss;
	debug_assert( new_abego.size() == pose_length_ );
	abego_ = new_abego;
	if ( pose_ ) {
		update_covalent_bonds_in_pose();
		save_into_pose();
	}
	TR.Debug << "Numbering updated - new pose length = " << pose_length_ << " new ss = " << new_ss << std::endl;
}

/// @brief updates movable group numbering after a deletion -- deleted mg is passed
void
StructureData::update_movable_groups_after_deletion( core::Size const /*mg_old*/ )
{
	/* TL: commented out because I don't think this is necessary anymore
	std::set< core::Size > mgs = movable_groups();
	core::Size replaced_mg = 0;

	for ( SegmentMap::iterator r = segments_.begin(); r != segments_.end(); ++r ) {
		if ( r->second.movable_group == mg_old ) {
			replaced_mg = r->second.movable_group;
			r->second.movable_group = mg_old;
		}
	}

	TR.Debug << "Renumbered movable group " << replaced_mg << " to " << mg_old << std::endl;
	*/
}

/// @brief renumbers movable group "oldg" to have new number "newg"
void
StructureData::renumber_movable_group( core::Size const oldg, core::Size const newg )
{
	for ( SegmentMap::iterator r = segments_.begin(); r != segments_.end(); ++r ) {
		if ( r->second.movable_group == oldg ) {
			r->second.movable_group = newg;
		}
	}
	update_movable_groups_after_deletion( oldg );
	save_into_pose();
}

/// @brief returns true if this object has a group of segments with the given name
bool
StructureData::has_segment_group( std::string const & sname ) const
{
	if ( has_segment( sname ) ) return true;
	std::string const match_str = sname + PARENT_DELIMETER;
	for ( StringList::const_iterator c=segments_begin(); c!=segments_end(); ++c ) {
		if ( boost::starts_with( *c, match_str ) ) {
			return true;
		}
	}
	return false;
}

/// @brief returns true if this object has a group of segments with the given name
StringList
StructureData::segment_group( std::string const & sname ) const
{
	std::string const match_str = sname + PARENT_DELIMETER;
	StringList retval;

	// check for perfect match first
	if ( has_segment( sname ) ) {
		retval.push_back( sname );
		return retval;
	}

	for ( StringList::const_iterator c = segments_begin(); c != segments_end(); ++c ) {
		if ( boost::starts_with( *c, match_str ) ) {
			retval.push_back( *c );
		}
	}
	return retval;
}

/// @brief checks consistency of the data
/// @throws EXCN_PoseInconsistent if there is a problem
void
StructureData::check_consistency() const
{
	if ( pose_ ) {
		check_consistency( *pose_ );
	} else {
		check_residues();
		check_movable_groups();
	}
}

void
StructureData::check_consistency( core::pose::Pose const & pose ) const
{
	check_residues();
	check_movable_groups();
	check_pose( pose );
	check_chain_endings( pose );
	check_chain_beginnings( pose );
	check_improper_termini( pose );
}

void
StructureData::check_residues() const
{
	std::stringstream msg;
	utility::vector1< bool > accounted_for( pose_length(), false );
	for ( SegmentMap::const_iterator r = segments_.begin(); r != segments_.end(); ++r ) {
		for ( core::Size i=r->second.nterm_resi(); i<=r->second.cterm_resi(); ++i ) {
			if ( !accounted_for[i] ) {
				accounted_for[i] = true;
			} else {
				msg << NamedSegment( r->first, r->second ) << " overlaps with something else at position " << i << std::endl << *this << std::endl;
				throw EXCN_PoseInconsistent( msg.str() );
			}
		}
		if ( r->second.lower_segment() != "" ) {
			SegmentMap::const_iterator r2 = segments_.find( r->second.lower_segment() );
			if ( r2 == segments_.end() ) {
				msg << "Lower segment of " << NamedSegment( r->first, r->second ) << " does not exist. SD=" << *this << std::endl;
				throw EXCN_PoseInconsistent( msg.str() );
			}
		}
		if ( r->second.upper_segment() != "" ) {
			SegmentMap::const_iterator r2 = segments_.find( r->second.upper_segment() );
			if ( r2 == segments_.end() ) {
				msg << "Upper segment of " << NamedSegment( r->first, r->second ) << " does not exist. SD=" << *this << std::endl;
				throw EXCN_PoseInconsistent( msg.str() );
			}
		}

	}
	for ( core::Size i=1; i<=pose_length(); ++i ) {
		if ( !accounted_for[i] ) {
			msg << " Residue " << i << " is not accounted for by the permutation. " << *this << std::endl;
			throw EXCN_PoseInconsistent( msg.str() );
		}
	}
}

void
StructureData::check_improper_termini( core::pose::Pose const & pose ) const
{
	for ( SegmentMap::const_iterator r=segments_.begin(); r!=segments_.end(); ++r ) {
		for ( core::Size i=r->second.nterm_resi()+1; i<r->second.cterm_resi(); ++i ) {
			if ( pose.residue(i).is_terminus() ) {
				std::stringstream msg;
				msg << " Residue " << i << " has a terminal variant but is inside segment " << NamedSegment( r->first, r->second ) << "." << std::endl;
				throw EXCN_PoseInconsistent( msg.str() );
			}
		}
	}
}

void
StructureData::check_chain_endings( core::pose::Pose const & pose ) const
{
	// check chain endings -- they should be on segment boundaries
	utility::vector1< core::Size > chain_endings = pose.conformation().chain_endings();
	if ( pose.total_residue() ) {
		chain_endings.push_back( pose.total_residue() );
	}
	for ( utility::vector1< core::Size >::const_iterator c_end=chain_endings.begin(); c_end!=chain_endings.end(); ++c_end ) {
		bool found = false;
		for ( SegmentMap::const_iterator s=segments_.begin(); s!=segments_.end(); ++s ) {
			if ( *c_end == s->second.cterm_resi() ) {
				found = true;
				break;
			}
		}
		if ( !found ) {
			std::stringstream msg;
			msg << "Chain ending " << *c_end << " does not match up with the upper terminus of any segment. " << *this << std::endl;
			throw EXCN_PoseInconsistent( msg.str() );
		}
	}
}

void
StructureData::check_chain_beginnings( core::pose::Pose const & pose ) const
{
	utility::vector1< core::Size > chain_beginnings( 1, 1 );
	for ( utility::vector1< core::Size >::const_iterator c_end=pose.conformation().chain_endings().begin(); c_end!=pose.conformation().chain_endings().end(); ++c_end ) {
		if ( chain_beginnings.size() == pose.conformation().chain_endings().size() ) break;
		chain_beginnings.push_back( *c_end + 1 );
	}
	for ( utility::vector1< core::Size >::const_iterator c_begin=chain_beginnings.begin(); c_begin!=chain_beginnings.end(); ++c_begin ) {
		bool found = false;
		for ( SegmentMap::const_iterator s=segments_.begin(); s!=segments_.end(); ++s ) {
			if ( *c_begin == s->second.nterm_resi() ) {
				found = true;
				break;
			}
		}
		if ( !found ) {
			std::stringstream msg;
			msg << "Chain beginning " << *c_begin << " does not match up with the lower terminus of any segment. " << *this << std::endl;
			throw EXCN_PoseInconsistent( msg.str() );
		}
	}
}

void
StructureData::check_movable_groups() const
{
	// check movable groups
	std::set< core::Size > const mg = movable_groups();
	TR.Debug << "Movable groups set is " << mg << std::endl;
	for ( std::set< core::Size >::const_iterator g = mg.begin(); g != mg.end(); ++g ) {
		if ( *g <= 0 ) {
			std::stringstream msg;
			msg << " StructureData has a movable group <= 0 : " << *this << std::endl;
			throw EXCN_PoseInconsistent( msg.str() );
		}
		/*
		if ( *g > mg.size() ) {
		std::stringstream msg;
		msg << " StructureData has a movable group (" << *g << ") >= size (" << mg.size() << ") : " << *this << std::endl;
		throw EXCN_PoseInconsistent( msg.str() );
		}
		*/
	}
}

/// @brief counts and returns the number of movable residue groups
core::Size
StructureData::movable_group( std::string const & id ) const
{
	return segment( id ).movable_group;
}

/// @brief sets movable group of a segment
void
StructureData::set_movable_group( std::string const & segid, core::Size const mg )
{
	SegmentMap::iterator s = segments_.find( segid );
	if ( s == segments_.end() ) {
		std::stringstream err;
		err << "Error in StructureData::set_movable_group( " << segid << ", " << mg << "):"
			<< "segment not found! perm=" << *this << std::endl;
		throw utility::excn::EXCN_Msg_Exception( err.str() );
	}
	s->second.movable_group = mg;
	save_into_pose();
}

/// @brief tells if the segment given has an available lower terminus
bool
StructureData::has_free_lower_terminus( std::string const & id_val ) const
{
	if ( segments_.find( id_val ) == segments_.end() ) {
		return false;
	}
	return segment( id_val ).has_free_lower_terminus();
}

/// @brief tells if the segment given has an available lower terminus
bool
StructureData::has_free_upper_terminus( std::string const & id_val ) const
{
	if ( segments_.find( id_val ) == segments_.end() ) {
		return false;
	}
	return segment( id_val ).has_free_upper_terminus();
}

/// @brief add lower cutpoint to residue cut and upper cutpoint to residue cut+1
void
StructureData::add_cutpoint_variants( core::Size const cut_res )
{
	debug_assert( pose_ );
	protocols::forge::methods::add_cutpoint_variants( *pose_, cut_res );
}

/// @brief removes cutpoint variants from residues cut and cut+1
void
StructureData::remove_cutpoint_variants( core::Size const cut_res )
{
	debug_assert( pose_ );
	debug_assert( cut_res + 1 <= pose_->total_residue() );

	if ( pose_->residue( cut_res ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) {
		core::pose::remove_variant_type_from_pose_residue( *pose_, core::chemical::CUTPOINT_LOWER, cut_res );
	}

	if ( pose_->residue( cut_res + 1 ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
		core::pose::remove_variant_type_from_pose_residue( *pose_, core::chemical::CUTPOINT_UPPER, cut_res + 1 );
	}
}

/// @brief applies an arbitrary mover to the contained pose
void
StructureData::apply_mover( protocols::moves::Mover & mover )
{
	core::Size const orig_len = pose_->total_residue();
	mover.apply( *pose_ );
	debug_assert( orig_len == pose_->total_residue() );
}

/// @brief applies an arbitrary mover to the contained pose
void
StructureData::apply_mover( protocols::moves::MoverOP mover )
{
	debug_assert( mover );
	apply_mover( *mover );
}

/// @brief removes constraints added by the given RCG
void
StructureData::remove_constraints_from_pose( protocols::forge::remodel::RemodelConstraintGeneratorOP rcg )
{
	debug_assert( pose_ );
	debug_assert( rcg );
	rcg->remove_remodel_constraints_from_pose( *pose_ );
}

/// @brief creates a new jump and cutpoint to build the given loop object
int
StructureData::new_jump_and_cutpoint( protocols::loops::Loop const & loop, core::Size const loop_overlap )
{
	core::Size const startres = loop_start_without_overlap( *pose_, loop.start(), loop_overlap );
	core::Size const stopres = loop_stop_without_overlap( *pose_, loop.stop(), loop_overlap );
	core::Size const cutres = loop.cut();
	core::Size saferes1_dist = 0;
	core::Size saferes1 = 0;
	// find the segment containing the cutpoint and set it
	for ( SegmentMap::iterator res = segments_.begin(); res != segments_.end(); ++res ) {
		if ( res->second.contains( cutres ) ) {
			res->second.set_cutpoint( cutres - res->second.nterm_resi() + 1 );
			break;
		}
	}
	// find the closest two safe residues to the loop and add a jump between them.
	for ( SegmentMap::const_iterator res = segments_.begin(); res != segments_.end(); ++res ) {
		Segment const & resis = res->second;
		if ( resis.safe() >= startres ) {
			continue;
		}
		if ( !saferes1 || ( startres-resis.safe() < saferes1_dist ) ) {
			saferes1 = resis.safe();
			saferes1_dist = startres-saferes1;
		}
	}
	core::Size saferes2_dist = 0;
	core::Size saferes2 = 0;
	for ( SegmentMap::const_iterator res = segments_.begin(); res != segments_.end(); ++res ) {
		Segment const & resis = res->second;
		if ( resis.safe() <= stopres ) {
			continue;
		}
		if ( !saferes2 || ( resis.safe()-stopres < saferes2_dist ) ) {
			saferes2 = resis.safe();
			saferes2_dist = saferes2-stopres;
		}
	}

	debug_assert( saferes1 < cutres );
	debug_assert( saferes2 > cutres );
	TR << "Inserting jump " << saferes1 << "__" << saferes2 << " with cut at " << cutres << std::endl;

	// modify fold tree
	core::kinematics::FoldTree ft = pose()->fold_tree();
	int const newjump = ft.new_jump( saferes1, saferes2, cutres );
	if ( !ft.check_fold_tree() ) {
		TR.Error << "FOLDTREE=" << ft << std::endl;
		TR.Error << "StructureData=" << *this << std::endl;
	}
	debug_assert( ft.check_fold_tree() );
	pose_->fold_tree( ft );

	// remove terminal variants (if applicable)
	if ( pose_->residue(cutres).is_upper_terminus() ) {
		remove_upper_terminus_variant_type( cutres );
	}
	if ( pose_->residue(cutres+1).is_lower_terminus() ) {
		remove_lower_terminus_variant_type( cutres+1 );
	}
	// add cutpoint variant types
	add_cutpoint_variants( cutres );
	rebuild_missing_atoms( *pose_, cutres );
	rebuild_missing_atoms( *pose_, cutres+1 );

	chains_from_termini();

	// they should ALWAYS be on the same chain at this point
	TR.Debug << pose_->fold_tree() << std::endl;
	TR.Debug << "chain of " << saferes1 << " = " << pose_->chain( saferes1 ) << " chain of " << saferes2 << " = " << pose_->chain( saferes2 ) << std::endl;
	for ( core::Size i=1, endi=pose_->total_residue(); i<=endi; ++i ) {
		TR.Debug << i << " " << pose_->residue(i).name() << std::endl;
	}
	debug_assert( pose_->chain( saferes1 ) == pose_->chain( saferes2 ) );

	return newjump;
}

/// output
std::ostream &
operator<<( std::ostream & os, StructureData const & perm )
{
	os << "<StructureData name=\"" << perm.id()
		<< "\" multi=\"" << perm.is_multi()
		<< "\" length=\"" << perm.length()
		<< "\" pose_length=\"" << perm.pose_length()
		<< "\" >" << std::endl;

	// residues
	for ( StringList::const_iterator c = perm.segment_order_.begin(); c != perm.segment_order_.end(); ++c ) {
		SegmentMap::const_iterator res = perm.segments_.find( *c );
		if ( res == perm.segments_.end() ) {
			std::stringstream err;
			err << perm.id() << ": segment not found = " << *c << std::endl;
			err << "Segment List = " << perm.segment_order_ << std::endl;
			throw utility::excn::EXCN_Msg_Exception( err.str() );
		}
		os << "\t" << NamedSegment(*res) << std::endl;
	}
	// int data
	std::map< std::string, int >::const_iterator dat;
	for ( dat = perm.data_int_.begin(); dat != perm.data_int_.end(); ++dat ) {
		os << "\t<Int name=\"" << dat->first << "\" value=\"" << dat->second << "\" />" << std::endl;
	}
	// real data
	std::map< std::string, core::Real >::const_iterator datr;
	for ( datr = perm.data_real_.begin(); datr != perm.data_real_.end(); ++datr ) {
		os << "\t<Real name=\"" << datr->first << "\" value=\"" << datr->second << "\" />" << std::endl;
	}
	// string data
	std::map< std::string, std::string >::const_iterator dats;
	for ( dats = perm.data_str_.begin(); dats != perm.data_str_.end(); ++dats ) {
		os << "\t<Str name=\"" << dats->first << "\" value=\"" << dats->second << "\" />" << std::endl;
	}
	// aliases
	for ( std::map< std::string, Alias >::const_iterator a_it=perm.aliases_.begin(); a_it!=perm.aliases_.end(); ++a_it ) {
		os << "\t<Alias name=\"" << a_it->first << "\" " << a_it->second << " />" << std::endl;
	}
	// covalent bonds
	for ( utility::vector1< BondInfo >::const_iterator b = perm.covalent_bonds_.begin(); b != perm.covalent_bonds_.end(); ++b ) {
		os << *b << std::endl;
	}
	os << "</StructureData>";
	return os;
}

/// dump contents of residues map
std::ostream &
operator<<( std::ostream & os, SegmentMap const & resmap )
{
	for ( SegmentMap::const_iterator r=resmap.begin(), endr=resmap.end(); r!=endr; ++r ) {
		os << NamedSegment( *r );
	}
	return os;
}

std::ostream &
operator<<( std::ostream & os, Alias const & alias )
{
	os << "segment=\"" << alias.first << "\" res=\"" << alias.second << "\"";
	return os;
}

} // namespace components
} // namespace denovo_design
} // namespace protocols
