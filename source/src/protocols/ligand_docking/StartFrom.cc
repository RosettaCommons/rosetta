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
#include <core/io/ResidueInformation.hh>
#include <core/io/AtomInformation.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/io/pdb/pdb_reader.hh>
#include <core/io/pose_from_sfr/PoseFromSFRBuilder.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

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

// Boost header
#include <boost/foreach.hpp>

// C++ header
#include <fstream>

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static THREAD_LOCAL basic::Tracer TR( "protocols.ligand_docking.StartFrom" );

std::string
StartFromCreator::keyname() const
{
	return StartFromCreator::mover_name();
}

protocols::moves::MoverOP
StartFromCreator::create_mover() const {
	return protocols::moves::MoverOP( new StartFrom );
}

std::string
StartFromCreator::mover_name()
{
	return "StartFrom";
}

/// @brief
StartFrom::StartFrom():
	Mover("StartFrom"),
	chain_(""),
	use_nbr_(false)
{}

StartFrom::StartFrom(StartFrom const & that):
	protocols::moves::Mover( that ),
	chain_(that.chain_),
	use_nbr_(that.use_nbr_),
	starting_positions_(that.starting_positions_),
	hash_starting_positions_(that.hash_starting_positions_)
{}

StartFrom::~StartFrom() {}

protocols::moves::MoverOP StartFrom::clone() const {
	return protocols::moves::MoverOP( new StartFrom( *this ) );
}

protocols::moves::MoverOP StartFrom::fresh_instance() const {
	return protocols::moves::MoverOP( new StartFrom );
}

std::string StartFrom::get_name() const{
	return "StartFrom";
}

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
		throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible");
	}
	if ( ! tag->hasOption("chain") ) throw utility::excn::EXCN_RosettaScriptsOption("'StartFrom' mover requires chain tag");

	chain( tag->getOption<std::string>("chain") );
	use_nbr( tag->getOption<bool>("use_nbr", false) );

	BOOST_FOREACH ( utility::tag::TagCOP child_tag, tag->getTags() ) {
		std::string name= child_tag->getName();
		if ( name == "features" || name == "Features" ) {
			TR << "found features tag with type '" << child_tag->getOption<std::string>("type") << "'" << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption("The 'Features' sub-tag in the StartFrom mover is not currently supported.");
		} else if ( name == "Coordinates" ) {
			if ( ! child_tag->hasOption("x") ) throw utility::excn::EXCN_RosettaScriptsOption("'StartFrom' mover Coordinates tag requires 'x' coordinates option");
			if ( ! child_tag->hasOption("y") ) throw utility::excn::EXCN_RosettaScriptsOption("'StartFrom' mover Coordinates tag requires 'y' coordinates option");
			if ( ! child_tag->hasOption("z") ) throw utility::excn::EXCN_RosettaScriptsOption("'StartFrom' mover Coordinates tag requires 'z' coordinates option");

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

			if ( !child_tag->hasOption("filename") ) throw utility::excn::EXCN_RosettaScriptsOption("'StartFrom' mover File tag requires 'filename' coordinates option");
			parse_startfrom_file(child_tag->getOption<std::string>("filename"));

		} else if ( name == "PDB" ) {
			if ( !child_tag->hasOption("filename") ) throw utility::excn::EXCN_RosettaScriptsOption("'StartFrom' mover File tag requires 'filename' coordinates option");
			std::string atom_name = child_tag->getOption<std::string>("atom_name","");
			std::string pdb_tag( child_tag->getOption<std::string>("pdb_tag","default") );

			parse_pdb_file(child_tag->getOption<std::string>("filename"), atom_name, pdb_tag);
		} else {
			throw utility::excn::EXCN_RosettaScriptsOption("'StartFrom' mover doesn't understand the sub-tag '"+name+"'" );
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

void StartFrom::apply(core::pose::Pose & pose) {

	if ( !core::pose::has_chain(chain_,pose) ) {
		utility_exit_with_message("StartFrom mover cannot find the chain " +chain_+ " in the current pose.");
	}

	utility::vector1<core::Vector> possible_centroids;

	// If we've already stored a startfrom, add that to the values to use.
	jd2::Job::StringRealPairs string_real_data(jd2::JobDistributor::get_instance()->current_job()->get_string_real_pairs());
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
		std::string hash = core::pose::get_sha1_hash_excluding_chain(chain_[0],pose);
		position_id = hash_starting_positions_.find(hash);
		if ( position_id != hash_starting_positions_.end() ) {
			specific_tag_found = true;
			possible_centroids.insert( possible_centroids.end(), position_id->second.begin(), position_id->second.end() );
			TR.Debug << "Using starting positions for hash " << hash << std::endl;
		}
	}

	// Now add all the job-tagged data.
	if ( ! starting_positions_.empty() && ! specific_tag_found ) {
		std::string const job_tag(jd2::JobDistributor::get_instance()->current_job()->input_tag());
		utility::vector1<std::string> const input_filenames(utility::split(job_tag));
		BOOST_FOREACH ( std::string filename, input_filenames ) {
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
		move_ligand_neighbor_to_desired_position(chain_, desired_centroid, pose);
	} else {
		move_ligand_to_desired_centroid(chain_, desired_centroid, pose);
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
	for ( utility::json_spirit::mArray::iterator start_it = start_positions.begin(); start_it != start_positions.end(); ++start_it ) {
		utility::json_spirit::mObject position_data(start_it->get_obj());

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

} //namespace ligand_docking
} //namespace protocols
