// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/StartFrom.cc
/// @brief  a mover to place a ligand at defined position(s)
/// @author Gordon Lemmon (glemmon@gmail.com)
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/StartFrom.hh>
#include <protocols/ligand_docking/StartFromCreator.hh>
#include <protocols/ligand_docking/util.hh>

// Package headers
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/io/ResidueInformation.hh>
#include <core/io/AtomInformation.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/io/pdb/pdb_reader.hh>
#include <core/io/pose_from_sfr/PoseFromSFRBuilder.hh>

#include <protocols/jd2/util.hh>

// Numeric header
#include <numeric/random/random.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>

#include <utility/json_spirit/json_spirit_reader.h>
#include <utility/excn/Exceptions.hh>

// Basic header
#include <basic/Tracer.hh>

// C++ header
#include <fstream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static basic::Tracer TR( "protocols.ligand_docking.StartFrom" );

// XRW TEMP std::string
// XRW TEMP StartFromCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return StartFrom::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP StartFromCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new StartFrom );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP StartFrom::mover_name()
// XRW TEMP {
// XRW TEMP  return "StartFrom";
// XRW TEMP }

/// @brief
StartFrom::StartFrom():
	//utility::pointer::ReferenceCount(),
	Mover("StartFrom"),
	chains_(),
	use_nbr_(false),
	starting_positions_(){}


StartFrom::StartFrom(StartFrom const & )= default;

StartFrom::~StartFrom() = default;

protocols::moves::MoverOP StartFrom::clone() const {
	return protocols::moves::MoverOP( new StartFrom( *this ) );
}

protocols::moves::MoverOP StartFrom::fresh_instance() const {
	return protocols::moves::MoverOP( new StartFrom );
}

// XRW TEMP std::string StartFrom::get_name() const{
// XRW TEMP  return "StartFrom";
// XRW TEMP }

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
StartFrom::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*datamap*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "StartFrom" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This should be impossible");
	}
	if ( ! tag->hasOption("chain") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'StartFrom' mover requires chain tag");


	std::string const all_chains_str = tag->getOption<std::string>("chain");
	chains_ = utility::string_split(all_chains_str, ',');

	use_nbr( tag->getOption<bool>("use_nbr", false) );

	for ( utility::tag::TagCOP child_tag : tag->getTags() ) {
		std::string name= child_tag->getName();
		if ( name == "features" || name == "Features" ) {
			TR << "found features tag with type '" << child_tag->getOption<std::string>("type") << "'" << std::endl;
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "The 'Features' sub-tag in the StartFrom mover is not currently supported.");
		} else if ( name == "Coordinates" ) {
			if ( ! child_tag->hasOption("x") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'StartFrom' mover Coordinates tag requires 'x' coordinates option");
			if ( ! child_tag->hasOption("y") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'StartFrom' mover Coordinates tag requires 'y' coordinates option");
			if ( ! child_tag->hasOption("z") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'StartFrom' mover Coordinates tag requires 'z' coordinates option");

			std::string pdb_tag( child_tag->getOption<std::string>("pdb_tag","default") );

			core::Vector v(
				child_tag->getOption<core::Real>("x"),
				child_tag->getOption<core::Real>("y"),
				child_tag->getOption<core::Real>("z")
			);
			add_coords(v,pdb_tag);

		} else if ( name == "File" ) {
			std::string location_id_mode = child_tag->getOption<std::string>("struct_identifier","");
			if ( location_id_mode.size() != 0 ) {
				TR.Debug << "No need to specify 'struct_identifier' in tag - pdb_id/hash will now be taken care of automatically." << std::endl;
			}

			if ( !child_tag->hasOption("filename") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'StartFrom' mover File tag requires 'filename' coordinates option");

			parse_startfrom_file(child_tag->getOption<std::string>("filename"));

		} else if ( name == "PDB" ) {
			if ( !child_tag->hasOption("filename") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'StartFrom' mover File tag requires 'filename' coordinates option");
			std::string atom_name = child_tag->getOption<std::string>("atom_name","");
			std::string pdb_tag( child_tag->getOption<std::string>("pdb_tag","default") );

			parse_pdb_file(child_tag->getOption<std::string>("filename"), atom_name, pdb_tag);
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'StartFrom' mover doesn't understand the sub-tag '"+name+"'" );
		}
	}
}

/// @brief Add the given coordinates as a valid starting position for the given pdb_tag
void StartFrom::add_coords(core::Vector const & coords,std::string const & pdb_tag /*="default"*/)
{
	// std::map::operator[]  will create an empty item if the key does not exist
	starting_positions_[ pdb_tag ].push_back(coords);
}

/// @brief Add the given coordinates as a valid starting position for the given hash value
void StartFrom::add_coords_hash(core::Vector const & coords,std::string const & hash_value)
{
	// std::map::operator[]  will create an empty item if the key does not exist
	hash_starting_positions_[ hash_value ].push_back(coords);
}

void StartFrom::apply(core::pose::Pose & pose){

	for ( std::string chain : chains_ ) {
		if ( !core::pose::has_chain(chain,pose) ) {
			utility_exit_with_message("StartFrom mover cannot find the chain " +chain+ " in the current pose.");
		}

		utility::vector1<core::Vector> possible_centroids;

		// If we've already stored a startfrom, add that to the values to use.
		std::map< std::string, core::Real > string_real_data(protocols::jd2::get_string_real_pairs_from_current_job());
		if ( string_real_data.find("start_x") != string_real_data.end() ) {
			core::Vector start_coords;
			start_coords.x(string_real_data["start_x"]);
			start_coords.y(string_real_data["start_y"]);
			start_coords.z(string_real_data["start_z"]);
			possible_centroids.push_back( start_coords );
			TR.Debug << "Using starting coords stored in job" <<std::endl;
		}

		std::map< std::string, utility::vector1< core::Vector > >::iterator position_id;
		bool specific_tag_found( false );

		// Add hash-tagged data, if any.
		if ( ! hash_starting_positions_.empty() ) {
			std::string hash = core::pose::get_sha1_hash_excluding_chain(chain[0],pose);
			position_id = hash_starting_positions_.find(hash);
			if ( position_id != hash_starting_positions_.end() ) {
				specific_tag_found = true;
				possible_centroids.insert( possible_centroids.end(), position_id->second.begin(), position_id->second.end() );
				TR.Debug << "Using starting positions for hash " << hash << std::endl;
			}
		}

		// Now add all the job-tagged data.
		if ( ! starting_positions_.empty() && ! specific_tag_found ) {
			std::string const job_tag( protocols::jd2::current_input_tag());
			utility::vector1<std::string> const input_filenames(utility::split(job_tag));
			for ( std::string filename : input_filenames ) {
				position_id = starting_positions_.find(filename);
				if ( position_id != starting_positions_.end() ) {
					specific_tag_found = true;
					possible_centroids.insert( possible_centroids.end(), position_id->second.begin(), position_id->second.end() );
					TR.Debug << "Using starting positions for path " << filename << std::endl;
					break; // Should we stack multiple values?
				} else {
					// Try finding just the basename, rather than the whole path.
					utility::file::FileName file_data(filename);
					std::string base_name(file_data.base());
					position_id = starting_positions_.find(base_name);
					if ( position_id != starting_positions_.end() ) {
						specific_tag_found = true;
						possible_centroids.insert( possible_centroids.end(), position_id->second.begin(), position_id->second.end() );
						TR.Debug << "Using starting positions for base filename " << base_name << std::endl;
						break;
					}
				}
			}
		}

		// Try the defaults, if nothing else better came along.
		if ( ! specific_tag_found ) {
			position_id = starting_positions_.find("default");
			if ( position_id != starting_positions_.end() ) {
				possible_centroids.insert( possible_centroids.end(), position_id->second.begin(), position_id->second.end() );
				TR.Debug << "Using default starting position set." << std::endl;
			}
		}

		if ( possible_centroids.size() == 0 ) {
			utility_exit_with_message("The current pose does not have behavior specified, and there are no default starting coordinates specified in the StartFrom mover");
		}

		core::Size const picked_centroid( numeric::random::rg().random_range(1, possible_centroids.size()) );
		core::Vector const desired_centroid = possible_centroids[picked_centroid];
		if ( use_nbr_ ) {
			move_ligand_neighbor_to_desired_position(chains_, desired_centroid, pose);
		} else {
			move_ligand_to_desired_centroid(chains_, desired_centroid, pose);
		}
	}
}


void StartFrom::parse_startfrom_file(std::string const & filename)
{
	if ( !utility::file::file_exists(filename) ) {
		utility_exit_with_message("cannot parse "+filename+" because it does not exist");
	}

	utility::io::izstream infile;
	infile.open(filename.c_str(),std::ifstream::in);
	utility::json_spirit::mValue startfrom_data;
	if ( !utility::json_spirit::read(infile,startfrom_data) ) {
		infile.close();
		utility_exit_with_message("cannot parse JSON in file "+ filename);
	}

	infile.close();

	//The format is something like this:
	/*
	[

	{
	"file_name" : "infile.pdb",
	"x" : 0.0020,
	"y" : -0.004,
	"z" : 0.0020,
	"hash" : "aa2aff055d19bc32e483df7ff4ae08361a768931"
	}

	]
	*/
	// Only one of 'file_name' or 'hash' needs to be present (hash will take precident if both are present)
	// "input_tag" is considered a synonym for "file_name"
	// If you just have the x/y/z coordinates, they will be added to the default set.

	utility::json_spirit::mArray start_positions = startfrom_data.get_array();
	for ( auto & start_position : start_positions ) {
		utility::json_spirit::mObject position_data(start_position.get_obj());

		core::Real x = position_data["x"].get_real();
		core::Real y = position_data["y"].get_real();
		core::Real z = position_data["z"].get_real();
		core::Vector coords(x,y,z);

		std::string identifier;
		if ( position_data.find("hash") != position_data.end() ) {
			identifier = position_data["hash"].get_str();
			add_coords_hash( coords, identifier );
		} else if ( position_data.find("file_name") != position_data.end() ) {
			identifier = position_data["file_name"].get_str();
			add_coords( coords, identifier );
		} else if ( position_data.find("input_tag") != position_data.end() ) {
			identifier = position_data["input_tag"].get_str();
			add_coords( coords, identifier );
		} else {
			add_coords( coords, "default" );
		}

	}

}

/// @brief Parse a PDB file, grabbing the positions from the heavy atom coordinates.
/// If atom_name is not empty, only grab coordinates from the specified atom name.
void StartFrom::parse_pdb_file(std::string const & filename, std::string const & atom_name /*= "" */, std::string const & tag /*="default"*/) {
	// We don't read the PDB file into a Pose because we don't know/care if all the residues
	// are something Rosetta knows about.
	utility::io::izstream file( filename );
	if ( ! file ) {
		utility_exit_with_message("Unable to read file '"+filename+"' when reading PDB file for StartFrom mover.");
	}
	std::string file_content;
	utility::slurp( file, file_content );
	core::io::StructFileRepOptions options;
	core::io::StructFileRep sfr( core::io::pdb::create_sfr_from_pdb_file_contents( file_content ) );
	utility::vector1< core::io::ResidueInformation > rinfos;
	core::io::pose_from_sfr::create_working_data( options, sfr, rinfos );

	for ( core::Size rr(1); rr <= rinfos.size(); ++rr ) {
		for ( core::Size aa(1); aa <= rinfos[rr].atoms().size(); ++aa ) {
			core::io::AtomInformation const & atom( rinfos[rr].atoms()[aa] );
			std::string element( utility::stripped_whitespace( atom.element ) );
			std::string name( utility::stripped_whitespace( atom.name ) );
			if ( element.size() == 0 ) {
				element = name[0]; // If the element column is empty, heuristic from the atom name field
			}
			if ( element == "H" ) { continue; }
			if ( atom_name.size() != 0 && name != utility::stripped_whitespace(atom_name) ) { continue; }

			core::Vector v( atom.x, atom.y, atom.z );
			add_coords(v,tag);
			TR.Debug << "Added position from PDB file to StartFrom mover: " << v.x() << "   " << v.y() << "   " << v.z() << std::endl;
		}
	}
}

std::string StartFrom::get_name() const {
	return mover_name();
}

std::string StartFrom::mover_name() {
	return "StartFrom";
}

void StartFrom::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("chain", xs_string, "Chain ID")
		+ XMLSchemaAttribute::attribute_w_default("use_nbr", xsct_rosetta_bool,
		"Place neighbor atom (the atom which is superimposed during "
		"conformer repacking) at the given location", "false");

	// subelements
	XMLSchemaSimpleSubelementList subelement_list;
	AttributeList Coordinates_attributes;
	Coordinates_attributes
		+ XMLSchemaAttribute::required_attribute("x", xsct_real, "x coordinate of the desired start coordinate.")
		+ XMLSchemaAttribute::required_attribute("y", xsct_real, "y coordinate of the desired start coordinate.")
		+ XMLSchemaAttribute::required_attribute("z", xsct_real, "z coordinate of the desired start coordinate.")
		+ XMLSchemaAttribute::attribute_w_default("pdb_tag", xs_string, "Job-tag of a pdb.", "default");
	subelement_list.add_simple_subelement("Coordinates", Coordinates_attributes, "Set x,y and z coordinates, where to put the chain.");

	AttributeList File_attributes;
	File_attributes
		+ XMLSchemaAttribute("struct_identifier", xs_string, "pdb_tag, file_name, or hash")
		+ XMLSchemaAttribute::required_attribute("filename", xs_string, "Name of SON formatted file containing starting positions");
	subelement_list.add_simple_subelement("File", File_attributes, "Provide a JSON formatted file containing starting positions");

	AttributeList PDB_attributes;
	PDB_attributes
		+ XMLSchemaAttribute::required_attribute("filename", xs_string, "Name of pdb file")
		+ XMLSchemaAttribute("atom_name", xs_string, "Use only the positions of atoms matching that atom name.")
		+ XMLSchemaAttribute::attribute_w_default("pdb_tag", xs_string, "XRW TO DO", "default");
	subelement_list.add_simple_subelement("PDB", PDB_attributes,
		"The starting coordinates are provided as the heavy atom positions in a PDB file.");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(),
		"Move a chain (normally a ligand) to a specified coordinate."
		"By default the centroid of the specified chain (the average position "
		"of all atoms/residues) is centered on the given coordinate",
		attlist, subelement_list );
}

std::string StartFromCreator::keyname() const {
	return StartFrom::mover_name();
}

protocols::moves::MoverOP
StartFromCreator::create_mover() const {
	return protocols::moves::MoverOP( new StartFrom );
}

void StartFromCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	StartFrom::provide_xml_schema( xsd );
}


} //namespace ligand_docking
} //namespace protocols
