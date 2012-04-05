// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/import_pose/import_pose.cc
///
/// @brief  various functions to construct Pose object(s) from PDB(s)
/// @author Sergey Lyskov

// Unit headers
#include <core/import_pose/import_pose.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/chemical/residue_io.hh>
#include <core/kinematics/FoldTree.hh>

// AUTO-REMOVED #include <core/scoring/Energies.hh>

// option key includes
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>


// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>
// AUTO-REMOVED #include <core/id/AtomID_Mask.hh>

// AUTO-REMOVED #include <ObjexxFCL/format.hh>
// AUTO-REMOVED #include <ObjexxFCL/ObjexxFCL.hh>

#include <core/pack/pack_missing_sidechains.hh>
#include <core/pack/optimizeH.hh>
// AUTO-REMOVED #include <core/io/raw_data/DisulfideFile.hh>
#include <core/pose/PDBInfo.hh>

#include <core/io/pdb/pdb_dynamic_reader.hh>
#include <core/io/pdb/file_data.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/string.functions.hh>

// option key includes

#include <basic/options/keys/inout.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>

#include <utility/vector1.hh>
#include <boost/foreach.hpp>

//Auto Headers
#define foreach BOOST_FOREACH
/// A temporary copy of the pose_from_pdb code from the demo directory.
/// Will be phased out in favor of file_data routines soon.
///

namespace core {
namespace import_pose {

using core::Size;
using core::SSize;

using basic::T;
using basic::Error;
using basic::Warning;

using namespace core::io::pdb;
using namespace ObjexxFCL;

basic::Tracer TR("core.import_pose.import_pose");

using utility::vector1;

//
void
read_additional_pdb_data(
	std::string const & s,
	pose::Pose & pose,
	io::pdb::FileData const &,
	bool read_fold_tree
)
{

	if ( (!read_fold_tree) && (!basic::options::option[ basic::options::OptionKeys::inout::fold_tree_io ].user()) ) return;

	// split on newlines
	utility::vector1< std::string > lines;
	Size start=0, i=0;
	while(start < s.size()) {
		if( s[i] == '\n' || s[i] == '\r' /* || i==s.size()-1 */) {
			lines.push_back( std::string(s.begin()+start, s.begin()+i) );
			start = i+1;
		}
		i++;
		if( i == s.size() ) {
			lines.push_back( std::string(s.begin()+start, s.begin()+i) );
			break;
		}
	}

	//
	for ( Size i=1; i<= lines.size(); ++i ) {
		std::string const & line( lines[i] );

		// look for fold_tree info
		if ( line.size() >= 16 && line.substr(0,16) == "REMARK FOLD_TREE" ) {
			std::istringstream l( line );
			std::string tag;
			kinematics::FoldTree f;
			l >> tag >> f;
			if ( !l.fail() && Size(f.nres()) == pose.total_residue() ) {
				TR << "setting foldtree from pdb file: " << f << std::endl;
				pose.fold_tree( f );
			} else {
				TR.Fatal << "pose_io:: foldtree io failure: " << line << ' ' << pose.total_residue()
					<< ' ' << f << std::endl;
				utility_exit();
			}
		}
	}

}



pose::PoseOP pose_from_pdb(std::string const & filename, bool read_fold_tree)
{
	pose::PoseOP pose = new pose::Pose();
	pose_from_pdb( *pose, filename, read_fold_tree);
	return pose;
}


pose::PoseOP pose_from_pdb(chemical::ResidueTypeSet const & residue_set, std::string const & filename, 	bool read_fold_tree)
{
	pose::PoseOP pose = new pose::Pose();
	pose_from_pdb( *pose, residue_set, filename, read_fold_tree);
	return pose;
}


void
pose_from_pdb(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filenames_string,
	bool read_fold_tree
)
{
	utility::vector1<std::string> filenames= utility::split(filenames_string);


	std::string res;

	foreach(std::string filename, filenames){
		utility::io::izstream file( filename );
		if (!file) {
			TR.Error << "PDB File:" << filename << " not found!" << std::endl;
			utility_exit_with_message( "Cannot open PDB file \"" + filename + "\"" );
		} else {
			TR.Debug << "read file: " << filename << std::endl;
		}
		utility::slurp( file, res );
	}

	//fpd If the conformation is not of type core::Conformation, reset it
	if ( !pose.conformation().same_type_as_me( conformation::Conformation(), true ) ) {
		pose.set_new_conformation( new conformation::Conformation() );
	}

	io::pdb::FileData fd = core::import_pose::PDB_DReader::createFileData(res);
	if ( fd.filename == "" ) {
		fd.filename = utility::join(filenames, "_");
	}
	build_pose(fd, pose, residue_set);

	// set secondary structure for centroid PDBs
	if ( residue_set.name() == core::chemical::CENTROID ) {
		core::pose::set_ss_from_phipsi( pose );
	}

	// check for foldtree info
	core::import_pose::read_additional_pdb_data( res, pose, fd, read_fold_tree);
}

void
pose_from_pdb(
	core::pose::Pose & pose,
	std::string const & filename,
	bool read_fold_tree
) {
	using basic::options::option;
	using namespace chemical;
	using namespace basic::options::OptionKeys;

	bool want_centroid( option[ in::file::centroid_input ]()
		|| option[ in::file::centroid ]()
		|| ( option[ in::file::fullatom ].user() && !option[ in::file::fullatom ]() )
		|| ( option[ in::file::residue_type_set ].user() && option[ in::file::residue_type_set ]() == "centroid" ) );

	if ( want_centroid &&
		( option[ in::file::fullatom ]()
			|| ( option[ in::file::residue_type_set ].user() && option[ in::file::residue_type_set ]() == "fa_standard" ) ) ) {
		TR.Warning << "conflicting command line flags for centroid/full-atom input. Choosing fullatom!" << std::endl;
		want_centroid = false;
	}
	ResidueTypeSetCAP residue_set( want_centroid ?
		ChemicalManager::get_instance()->residue_type_set( CENTROID ) :
		ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )
	);

	if ( option[ in::file::residue_type_set ].user() && option[ in::file::residue_type_set]()  == "rna" ) residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	core::import_pose::pose_from_pdb( pose, *residue_set, filename, read_fold_tree );
}

utility::vector1< core::pose::PoseOP > poseOPs_from_pdbs(
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree
) {
	using basic::options::option;
	using namespace chemical;
	using namespace basic::options::OptionKeys;

	bool want_centroid( option[ in::file::centroid_input ]()
		|| option[ in::file::centroid ]()
		|| ( option[ in::file::fullatom ].user() && !option[ in::file::fullatom ]() )
		|| ( option[ in::file::residue_type_set ].user() && option[ in::file::residue_type_set ]() == "centroid" ) );
	if ( want_centroid &&
		( option[ in::file::fullatom ]()
			|| ( option[ in::file::residue_type_set ].user() && option[ in::file::residue_type_set ]() == "fa_standard" ) ) ) {
		TR.Warning << "conflicting command line flags for centroid/full-atom input. Choosing fullatom!" << std::endl;
		want_centroid = false;
	}
	ResidueTypeSetCAP residue_set( want_centroid ?
		ChemicalManager::get_instance()->residue_type_set( CENTROID ) :
		ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )
	);

	return core::import_pose::poseOPs_from_pdbs( *residue_set, filenames, read_fold_tree );
}

utility::vector1< core::pose::Pose > poses_from_pdbs(
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree
) {
	using namespace chemical;
	ResidueTypeSetCAP residue_set
		( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

	return core::import_pose::poses_from_pdbs( *residue_set, filenames, read_fold_tree );
}

utility::vector1< core::pose::Pose > poses_from_pdbs(
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree
) {
	using namespace chemical;

	using std::string;
	using utility::vector1;
	using core::pose::Pose;

	vector1< Pose > poses;
	typedef vector1< string >::const_iterator vec_it;
	for ( vec_it it = filenames.begin(), end = filenames.end(); it != end; ++it ) {
		Pose pose;
		core::import_pose::pose_from_pdb( pose, residue_set, *it, read_fold_tree );
		poses.push_back( pose );
	}

	return poses;
}

void
read_all_poses(const utility::vector1<std::string>& filenames,
							 utility::vector1<core::pose::Pose>* poses) {
	using core::pose::Pose;
	using std::string;
	using utility::vector1;
	assert(poses);

	for (vector1<string>::const_iterator i = filenames.begin(); i != filenames.end(); ++i) {
		Pose pose;
		pose_from_pdb(pose, *i);
		poses->push_back(pose);
	}
}

utility::vector1< core::pose::PoseOP > poseOPs_from_pdbs(
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree
) {
	using namespace chemical;

	using std::string;
	using utility::vector1;
	using core::pose::Pose;

	vector1< pose::PoseOP > poses;
	typedef vector1< string >::const_iterator vec_it;
	for ( vec_it it = filenames.begin(), end = filenames.end(); it != end; ++it ) {
		pose::PoseOP pose = new pose::Pose;
		core::import_pose::pose_from_pdb( *pose, residue_set, *it, read_fold_tree );
		poses.push_back( pose );
	}

	return poses;
}


void
pose_from_pdb(
	utility::vector1< pose::Pose > & poses,
	std::string const & filename,
  bool read_fold_tree
)
{
	using namespace chemical;
	ResidueTypeSetCAP residue_set(
		ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )
	);
	core::import_pose::pose_from_pdb( poses, *residue_set, filename, read_fold_tree );
}


void
pose_from_pdb(
	utility::vector1< pose::Pose > & poses,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	bool read_fold_tree
)
{
	// Size fsize;
	pose::Pose pose;
	std::string all_lines, sub_lines;

	utility::io::izstream file( filename );
	if (!file) {
		TR.Error << "File:" << filename << " not found!" << std::endl;
		utility_exit_with_message( "Cannot open file " + filename );
	} else {
		TR.Debug << "read file: " << filename << std::endl;
	}

	utility::slurp( file, all_lines );
	// Hacky way to make sure that we find all models. I blame society.
	all_lines = "\n" + all_lines;

	Size pos1 = 0;
	Size pos2 = 0;
	Size n_models = 0;

	// count the number of poses will be reading
	while( pos1 != std::string::npos )
	{
		pos1 = all_lines.find( "\nMODEL ", pos1 );
		if ( pos1 != std::string::npos )
		{
			++n_models;
			++pos1;
		}
	}

	TR.Debug << "Reading " << n_models << " poses." << std::endl;
	// make space for all of our poses
	if ( n_models == 0 ) n_models = 1;
	poses.reserve( poses.size() + n_models );

	pos1 = 0;

	pos1 = all_lines.find( "\nMODEL ", pos1 );
	if ( pos1 != std::string::npos ) {
		pos2 = 0;
		while( pos2 != std::string::npos ) {
			// set pos1 == "M" in MODEL
			++pos1;

			// pos2 = position of newline character, start somewhere after pos1
			pos2 = all_lines.find( "\nMODEL ", pos1);
			sub_lines = all_lines.substr( pos1, pos2-pos1 ) ;
			pos1 = pos2;

			io::pdb::FileData fd = core::import_pose::PDB_DReader::createFileData( sub_lines );
			fd.filename = filename;
			build_pose( fd, pose, residue_set);

			// check for foldtree info
			core::import_pose::read_additional_pdb_data( sub_lines, pose, fd, read_fold_tree);
			poses.push_back( pose );
		}
	} else {
		FileData fd = core::import_pose::PDB_DReader::createFileData( all_lines );
		if ( fd.filename == "" ) {
			fd.filename = filename;
		}
		core::import_pose::build_pose( fd, pose, residue_set);
		// check for foldtree info
		core::import_pose::read_additional_pdb_data( all_lines, pose, fd, read_fold_tree);
		poses.push_back( pose );
	}
}

void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	std::string const & filename
)
{
	io::pdb::FileData fd = import_pose::PDB_DReader::createFileData( pdbcontents );
	fd.filename = filename;
	chemical::ResidueTypeSetCAP residue_set
		( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );
	core::import_pose::build_pose( fd, pose, *residue_set);

}

void
centroid_pose_from_pdb(
	pose::Pose & pose,
	std::string const & filename,
	bool read_fold_tree
)
{
	using namespace chemical;
	ResidueTypeSetCAP residue_set
		( ChemicalManager::get_instance()->residue_type_set( CENTROID ) );

	core::import_pose::pose_from_pdb( pose, *residue_set, filename, read_fold_tree);
}

///
/// @details Build mini Rosetta pose object from FileData.
///
void build_pose(io::pdb::FileData & fd, pose::Pose & pose, chemical::ResidueTypeSet const & residue_set)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	TR.Debug << "build_pose..." << std::endl;
	build_pose_as_is( fd, pose, residue_set);
	TR.Debug << "build_pose... Ok." << std::endl;
}

//void build_pose_as_is1( io::pdb::FileData & fd, pose::Pose & pose, chemical::ResidueTypeSet const & residue_set, id::AtomID_Mask & missing );
void build_pose_as_is2( io::pdb::FileData & fd, pose::Pose & pose, chemical::ResidueTypeSet const & residue_set, id::AtomID_Mask & missing );

/// @details
/// trying to Build pose object from pdb 'as-is'. PDB file must be _really_ clean.
///
///////////////////////////////////////////////////////////////////////////////
// "super-simple" (C) by Phil
//
void build_pose_as_is( io::pdb::FileData & fd, pose::Pose & pose, chemical::ResidueTypeSet const & residue_set)
{
	id::AtomID_Mask missing( false );

	io::pdb::build_pose_as_is1( fd, pose, residue_set, missing );
	build_pose_as_is2( fd, pose, residue_set, missing );
}

void build_pose_as_is2( io::pdb::FileData & /*fd*/, pose::Pose & pose, chemical::ResidueTypeSet const & residue_set, id::AtomID_Mask & missing )
{
	using namespace chemical;
	using namespace conformation;
	typedef std::map< std::string, double > ResidueTemps;
	typedef std::map< std::string, ResidueTemps > Temps;
	typedef std::map< std::string, Vector > ResidueCoords;
	typedef std::map< std::string, ResidueCoords > Coords;
	typedef utility::vector1< std::string > Strings;

	core::pose::PDBInfoOP pdb_info( new core::pose::PDBInfo(*pose.pdb_info()) );

	if ( !(basic::options::option[ basic::options::OptionKeys::run::skip_set_reasonable_fold_tree ].value()) ) {
		set_reasonable_fold_tree( pose );
	}

	/// optimize H if using a fullatom residue type set, and no_optH is not specified
	if ( residue_set.name() == FA_STANDARD ) {
		//if pack_missing_density specified, repack residues w/ missing density
		if( basic::options::option[ basic::options::OptionKeys::packing::pack_missing_sidechains ].value() ) {
			pack::pack_missing_sidechains( pose, missing );
		}
		/// optimize H if using a fullatom residue type set, and no_optH is not specified
		if( !basic::options::option[ basic::options::OptionKeys::packing::no_optH ]() ) {
			pack::optimize_H_and_notify( pose, missing );
		}
	}

	pose.pdb_info( pdb_info );
}

/// @details
/// all ligand residues have been appended by a jump -- create a fold tree without jumps through
/// ligands.
void
set_reasonable_fold_tree(
	pose::Pose & pose
)
{
	// An empty pose doesn't have jumps through ligands
	// (Will encounter a SegFault otherwise)
	if( pose.total_residue() == 0 ) {
		return;
	}

	using namespace core::kinematics;
	TR.Debug << "original fold tree: " << pose.fold_tree() << std::endl;

	FoldTree const & origft = pose.fold_tree();
	FoldTree newft;

	// Copy the original fold tree edges to the new fold tree, possibly replacing
	// some jumps to ligand residues with ligand-ligand chemical bonds.
	// As a result, jumps must be renumbered.
	Size last_jump_id = 0;
	for( FoldTree::const_iterator i = origft.begin(), i_end = origft.end(); i != i_end; ++i ) {
		Edge e = *i;
		// Jump to a ligand residue?
		if( e.is_jump() && !pose.residue_type(e.stop()).is_polymer() ) {
			Size const ii = e.stop(); // the ligand residue at the end of the jump
			conformation::Residue const & ii_res = pose.residue(ii);
			// Now we'll prepare a chemical edge by first finding the connecting atoms
			// between the two residues
			bool found_connection_residue_in_fold_tree( false );
			for ( Size jj = 1; jj <= pose.residue_type( ii ).n_residue_connections(); ++jj ) {
				if(ii_res.connection_incomplete(jj)) continue;  // allow incomplete connections for design
				Size jj_res_ID= ii_res.connect_map( jj ).resid();
				if ( jj_res_ID < ii){
					core::conformation::Residue const &  jj_res=pose.residue(jj_res_ID);
					if( jj_res.chain() == ii_res.chain() ) {
						int ii_connect_ID= ii_res.connect_atom(jj_res);
						int jj_connect_ID= jj_res.connect_atom(ii_res);

						std::string ii_connector= ii_res.atom_name(ii_connect_ID);
						std::string jj_connector= jj_res.atom_name(jj_connect_ID);

						newft.add_edge( jj_res_ID, ii, jj_connector, ii_connector);

						found_connection_residue_in_fold_tree = true;
						break;
					}
				}
			}
			// If we couldn't find a chemical bond to make, then we should probably just keep the jump as-is.
			if( ! found_connection_residue_in_fold_tree ) {
				if( pose.residue_type( ii ).n_residue_connections() > 0 ) {
					TR.Warning << "Can't find a chemical connection for residue " << ii << " " << pose.residue_type(ii).name() << std::endl;
				}
				if( e.is_jump() ) e.label() = ++last_jump_id;
				newft.add_edge(e);
			}
		} else { // no, just a normal peptide edge, inter-chain jump, etc.
			if( e.is_jump() ) e.label() = ++last_jump_id;
			newft.add_edge(e);
		}
	}

	runtime_assert( newft.size() > 0 || pose.total_residue() == 0 ); // valid fold tree must have 1+ edges

	TR.Debug << "new fold tree " << newft << std::endl;
	pose.fold_tree( newft );
}

} // namespace import_pose
} // namespace core
