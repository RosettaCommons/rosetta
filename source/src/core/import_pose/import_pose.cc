// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/import_pose.cc
/// @brief  various functions to construct Pose object(s) from PDB(s)
/// @details A temporary copy of the pose_from_pdb code from the demo directory.
/// Will be phased out in favor of file_data routines soon.
/// @author Sergey Lyskov

// Unit headers
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose_options.hh>

// Project headers
#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/carbohydrates/util.hh>

#include <core/util/metalloproteins_util.hh>

#include <core/pack/pack_missing_sidechains.hh>
#include <core/pack/optimizeH.hh>

#include <core/io/pdb/pdb_reader.hh>
#include <core/io/StructFileReaderOptions.hh>  // TODO: Rename after refactor.
#include <core/io/pose_from_sfr/PoseFromSFRBuilder.hh>
#include <core/io/mmcif/cif_reader.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>

#include <core/io/pdb/pdb_reader.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <cifparse/CifFile.h>
#include <cifparse/CifParserBase.h>
typedef utility::pointer::shared_ptr< CifFile > CifFileOP;
typedef utility::pointer::shared_ptr< CifParser > CifParserOP;


namespace core {
namespace import_pose {

using core::Size;
using core::SSize;

using basic::T;
using basic::Error;
using basic::Warning;

using namespace core::io;
using namespace ObjexxFCL;

static THREAD_LOCAL basic::Tracer TR( "core.import_pose.import_pose" );

using utility::vector1;

void
read_all_poses(
	utility::vector1< std::string > const & filenames,
	utility::vector1< core::pose::PoseOP > & poses
)
{
	using core::pose::Pose;
	using std::string;
	using utility::vector1;

	for ( auto const & filename : filenames ) {
		core::pose::PoseOP pose( new Pose );
		pose_from_file(*pose, filename, core::import_pose::PDB_file);
		poses.push_back(pose);
	}
}


void
read_additional_pdb_data(
	std::string const & s,
	pose::Pose & pose,
	io::StructFileRepCOP ,
	bool read_fold_tree
)
{
	ImportPoseOptions options;
	read_additional_pdb_data( s, pose, options, read_fold_tree );
}

void
read_additional_pdb_data(
	std::string const & s,
	pose::Pose & pose,
	ImportPoseOptions const & options,
	bool read_fold_tree
)
{
	// split additional_pdb_data s into newlines
	utility::vector1< std::string > lines;
	Size start=0, i=0;

	while ( start < s.size() ) {
		if ( s[i] == '\n' || s[i] == '\r' /* || i==s.size()-1 */ ) {
			lines.push_back( std::string(s.begin()+start, s.begin()+i) );
			start = i+1;
		}
		i++;
		if ( i == s.size() ) {
			lines.push_back( std::string(s.begin()+start, s.begin()+i) );
			break;
		}
	}

	//Added by Daniel-Adriano Silva, used to read PDBinfo-LABEL
	//TR.Debug << "Setting PDBinfo-labels from PDB file." << std::endl;
	for ( i=1; i<= lines.size(); ++i ) {
		std::string const & line( lines[i] );
		if ( line.size() > 21 && line.substr(0,21) == "REMARK PDBinfo-LABEL:" ) {
			//Parse and split string
			utility::vector1 < std::string > remark_values;
			utility::vector1 < std::string > tmp_remark_values = utility::string_split(line, ' ');
			//Copy non-empty (i.e. !' ') elements to remark_values
			if ( tmp_remark_values .size() > 3 ) {
				for ( Size j=3; j<= tmp_remark_values.size(); ++j ) {
					if ( tmp_remark_values[j] != "" ) {
						remark_values.push_back(tmp_remark_values[j]);
					}
				}
			}
			//Check that we have at least two elements left ([1]=index, [2-n]=PDBinfo-labels)
			if ( remark_values.size() > 1 ) {
				core::Size tmp_ndx=atoi(remark_values[1].c_str());
				if ( tmp_ndx <= pose.size() ) {
					for ( Size j=2; j<= remark_values.size(); ++j ) {
						pose.pdb_info()->add_reslabel(tmp_ndx,remark_values[j]);
					}
				} else {
					TR.Fatal << "pose_io:: PDBinfo-LABEL io failure: " << line << ' ' << pose.size()  << std::endl;
				}
			} else {
				TR.Fatal << "pose_io:: PDBinfo-LABEL io failure: " << line << ' ' << pose.size() << std::endl;
			}
		}
	}

	if ( (!read_fold_tree) && (!options.fold_tree_io()) ) return;

	for ( i=1; i<= lines.size(); ++i ) {
		std::string const & line( lines[i] );
		// Look for fold_tree info
		if ( line.size() >= 16 && line.substr(0,16) == "REMARK FOLD_TREE" ) {
			std::istringstream l( line );
			std::string tag;
			kinematics::FoldTree f;
			l >> tag >> f;
			if ( !l.fail() && Size(f.nres()) == pose.size() ) {
				TR << "setting foldtree from pdb file: " << f << std::endl;
				pose.fold_tree( f );
			} else {
				TR.Fatal << "pose_io:: foldtree io failure: " << line << ' ' << pose.size()
					<< ' ' << f << std::endl;
				utility_exit();
			}
		}
	}
}


pose::PoseOP pose_from_file( std::string const & filename, bool read_fold_tree, FileType type )
{
	pose::PoseOP pose( new pose::Pose() );
	pose_from_file( *pose, filename, read_fold_tree, type);
	return pose;
}

pose::PoseOP
pose_from_file(
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree,
	FileType type
)
{
	pose::PoseOP pose( new pose::Pose() );
	pose_from_file( *pose, filename, options, read_fold_tree, type);
	return pose;
}


pose::PoseOP pose_from_file(chemical::ResidueTypeSet const & residue_set, std::string const & filename,  bool read_fold_tree, FileType type)
{
	pose::PoseOP pose( new pose::Pose() );
	pose_from_file( *pose, residue_set, filename, read_fold_tree, type);
	return pose;
}

FileType
determine_file_type( std::string const &contents_of_file) {

	//See if this is a pdb file
	//code to determine the type of file
	utility::vector1< io::pdb::Record> records( io::pdb::create_records_from_pdb_file_contents( contents_of_file ) );
	for ( core::Size ii=1; ii<= records.size(); ++ii ) {
		if ( records[ ii ][ "type" ].value != "UNKNOW" ) {
			return PDB_file;
		}
	}
	//See if this is a CIF file
	std::string diagnostics;
	CifFileOP cifFile( new CifFile );

	CifParserOP cifParser( new CifParser( cifFile.get() ) );

	cifParser->ParseString( contents_of_file, diagnostics );
	if ( diagnostics.empty() ) {
		return CIF_file;
	}

	return Unknown_file;
}


void
pose_from_file(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filenames_string,
	ImportPoseOptions const & options,
	bool read_fold_tree,
	FileType file_type
) {
	utility::vector1< std::string > filenames = utility::split(filenames_string);

	std::string contents_of_file;

	for ( std::string const & filename : filenames ) {
		utility::io::izstream file( filename );
		if ( !file ) {
			TR.Error << "File: " << filename << " not found!" << std::endl;
			utility_exit_with_message( "Cannot open file \"" + filename + "\"" );
		} else {
			TR.Debug << "read file: " << filename << std::endl;
		}
		utility::slurp( file, contents_of_file );
	}

	if ( file_type == Unknown_file ) {
		file_type = determine_file_type( contents_of_file);
	}

	if ( file_type == Unknown_file ) {
		utility_exit_with_message( "Cannot determine file type. Current supported types are: PDB, CIF");
	} else if ( file_type == PDB_file ) {
		//fpd If the conformation is not of type core::Conformation, reset it
		conformation::ConformationOP conformation_op( new conformation::Conformation() );
		if ( !pose.conformation().same_type_as_me( *conformation_op, true ) ) {
			pose.set_new_conformation( conformation_op );
		}

		io::StructFileRepOP sfr( io::pdb::create_sfr_from_pdb_file_contents(contents_of_file, options).clone() );
		if ( sfr->filename() == "" ) {
			sfr->filename() = utility::join(filenames, "_");
		}

		build_pose( sfr, pose, residue_set, options );

		// set secondary structure for centroid PDBs
		if ( residue_set.mode() == core::chemical::CENTROID_t ) {
			core::pose::set_ss_from_phipsi( pose );
		}

		// check for foldtree info
		read_additional_pdb_data( contents_of_file, pose, options, read_fold_tree );

	} else if ( file_type == CIF_file ) {
		std::string diagnostics;
		CifFileOP cifFile( new CifFile );
		CifParserOP cifParser( new CifParser( cifFile.get() ) );
		cifParser->ParseString( contents_of_file, diagnostics );
		io::StructFileRepOP sfr ( io::mmcif::create_sfr_from_cif_file_op( cifFile, options ) );
		build_pose( sfr, pose, residue_set, options );

	}


}

void
pose_from_file(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filenames_string,
	bool read_fold_tree,
	FileType type
)
{
	ImportPoseOptions options;
	pose_from_file(pose, residue_set, filenames_string, options, read_fold_tree, type);
}



void
pose_from_file(
	core::pose::Pose & pose,
	std::string const & filename,
	bool read_fold_tree,
	FileType type
) {
	ImportPoseOptions options;
	pose_from_file( pose, filename, options, read_fold_tree , type);
}

void
pose_from_file(
	pose::Pose & pose,
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree,
	FileType type
) {
	using namespace chemical;

	ResidueTypeSetCOP residue_set( options.centroid() ?
		pose.residue_type_set_for_pose( CENTROID_t ) :
		pose.residue_type_set_for_pose( FULL_ATOM_t )
	);

	core::import_pose::pose_from_file( pose, *residue_set, filename, options, read_fold_tree, type);
}

utility::vector1< core::pose::PoseOP >
poseOPs_from_files(
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree,
	FileType type
) {
	ImportPoseOptions options;
	return poseOPs_from_files( filenames, options, read_fold_tree, type );
}

utility::vector1< core::pose::PoseOP >
poseOPs_from_files(
	utility::vector1< std::string > const & filenames,
	ImportPoseOptions const & options,
	bool read_fold_tree,
	FileType type
) {
	using namespace chemical;

	ResidueTypeSetCOP residue_set( options.centroid() ?
		ChemicalManager::get_instance()->residue_type_set( CENTROID ) :
		ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )
	);

	return core::import_pose::poseOPs_from_files( *residue_set, filenames, options, read_fold_tree, type);
}

/// @details Only returns full-atom poses
utility::vector1< core::pose::Pose >
poses_from_files(
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree,
	FileType type
) {
	using namespace chemical;
	ResidueTypeSetCOP residue_set
		( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

	return core::import_pose::poses_from_files( *residue_set, filenames, read_fold_tree, type);
}

utility::vector1< core::pose::Pose >
poses_from_files(
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree,
	FileType type
) {
	using namespace chemical;

	using std::string;
	using utility::vector1;
	using core::pose::Pose;

	ImportPoseOptions options;

	vector1< Pose > poses;
	for ( auto const & filename : filenames ) {
		Pose pose;
		core::import_pose::pose_from_file( pose, residue_set, filename, options, read_fold_tree, type);
		poses.push_back( pose );
	}

	return poses;
}


utility::vector1< core::pose::PoseOP >
poseOPs_from_files(
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< std::string > const & filenames,
	ImportPoseOptions const & options,
	bool read_fold_tree,
	FileType type
) {
	using namespace chemical;

	using std::string;
	using utility::vector1;
	using core::pose::Pose;

	vector1< pose::PoseOP > poses;
	for ( auto const & filename : filenames ) {
		pose::PoseOP pose( new pose::Pose );
		core::import_pose::pose_from_file( *pose, residue_set, filename, options, read_fold_tree, type);
		poses.push_back( pose );
	}

	return poses;
}


void
pose_from_file(
	utility::vector1< pose::Pose > & poses,
	std::string const & filename,
	bool read_fold_tree,
	FileType type
)
{
	using namespace chemical;
	ResidueTypeSetCOP residue_set(
		ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )
	);
	core::import_pose::pose_from_file( poses, *residue_set, filename, read_fold_tree, type);
}


void
pose_from_file(
	utility::vector1< pose::Pose > & poses,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	bool read_fold_tree,
	FileType type
)
{
	ImportPoseOptions options;
	pose_from_file( poses, residue_set, filename, options, read_fold_tree, type);
}

void
pose_from_file(
	utility::vector1< pose::Pose > & poses,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree,
	FileType type
)
{
	// Size fsize;
	pose::Pose pose;
	std::string all_lines, sub_lines;

	utility::io::izstream file( filename );
	if ( !file ) {
		TR.Error << "File:" << filename << " not found!" << std::endl;
		utility_exit_with_message( "Cannot open file " + filename );
	} else {
		TR.Debug << "read file: " << filename << std::endl;
	}

	utility::slurp( file, all_lines );


	if ( type == PDB_file || type == Unknown_file ) {
		// Hacky way to make sure that we find all models. I blame society.
		all_lines = "\n" + all_lines;

		Size pos1 = 0;
		Size pos2 = 0;
		Size n_models = 0;

		// count the number of poses will be reading
		while ( pos1 != std::string::npos )
				{
			pos1 = all_lines.find( "\nMODEL ", pos1 );
			if ( pos1 != std::string::npos ) {
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
			while ( pos2 != std::string::npos ) {
				// set pos1 == "M" in MODEL
				++pos1;

				// pos2 = position of newline character, start somewhere after pos1
				pos2 = all_lines.find( "\nMODEL ", pos1);
				sub_lines = all_lines.substr( pos1, pos2-pos1 ) ;
				pos1 = pos2;

				io::StructFileRepOP sfr( io::pdb::create_sfr_from_pdb_file_contents( sub_lines, options ).clone() );
				sfr->filename() = filename;
				build_pose( sfr, pose, residue_set, options );

				// check for foldtree info
				core::import_pose::read_additional_pdb_data( sub_lines, pose, options, read_fold_tree);
				poses.push_back( pose );
			}
		} else {
			StructFileRepOP sfr( io::pdb::create_sfr_from_pdb_file_contents( all_lines, options ).clone() );
			if ( sfr->filename() == "" ) {
				sfr->filename() = filename;
			}
			core::import_pose::build_pose( sfr, pose, residue_set, options );
			// check for foldtree info
			core::import_pose::read_additional_pdb_data( all_lines, pose, options, read_fold_tree);
			poses.push_back( pose );
		}
	}
}

void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	ImportPoseOptions const & options,
	std::string const & filename
)
{
	io::StructFileRepOP sfr( io::pdb::create_sfr_from_pdb_file_contents( pdbcontents, options ).clone() );
	sfr->filename() = filename;
	chemical::ResidueTypeSetCOP residue_set
		( pose.residue_type_set_for_pose( chemical::FULL_ATOM_t ) );
	core::import_pose::build_pose( sfr, pose, *residue_set, options );

}

void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	std::string const & filename
)
{
	ImportPoseOptions options;
	pose_from_pdbstring( pose, pdbcontents, options, filename );
}


void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename
){
	ImportPoseOptions options;
	pose_from_pdbstring( pose, pdbcontents, residue_set, options, filename );
}

void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options,
	std::string const & filename
){
	io::StructFileRepOP sfr( io::pdb::create_sfr_from_pdb_file_contents( pdbcontents, options ).clone() );
	sfr->filename() = filename;
	core::import_pose::build_pose( sfr, pose, residue_set, options);
}

void pose_from_pdb_stream(
	pose::Pose & pose,
	std::istream & pdb_stream,
	std::string const & filename,
	ImportPoseOptions const & options
){
	std::string pdb_file_contents;
	utility::slurp( pdb_stream, pdb_file_contents );

	chemical::ResidueTypeSetCOP residue_set( pose.residue_type_set_for_pose( options.residue_type_set_mode() ) );
	pose_from_pdbstring( pose, pdb_file_contents, *residue_set, options, filename);
	read_additional_pdb_data( pdb_file_contents, pose, options, options.read_fold_tree() );
}

void
centroid_pose_from_pdb(
	pose::Pose & pose,
	std::string const & filename,
	bool read_fold_tree
)
{
	using namespace chemical;
	ResidueTypeSetCOP residue_set
		( pose.residue_type_set_for_pose( CENTROID_t ) );

	core::import_pose::pose_from_file( pose, *residue_set, filename, read_fold_tree, core::import_pose::PDB_file);
}

void build_pose(
	io::StructFileRepOP fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set
)
{
	ImportPoseOptions options; // read from the command line
	build_pose( fd, pose, residue_set, options);
}

/// @brief Build Rosetta 3 Pose object from StructFileRep.
void build_pose(
	io::StructFileRepOP fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options
)
{
	TR.Debug << "build_pose..." << std::endl;
	build_pose_as_is( fd, pose, residue_set, options);
	TR.Debug << "build_pose... Ok." << std::endl;
}

// "super-simple" (C) by Phil
/// @brief Try to Build pose object from pdb 'as-is'. PDB file must be _really_ clean.
void build_pose_as_is(
	io::StructFileRepOP sfr,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options
)
{
	id::AtomID_Mask missing( false );

	//io::pdb::build_pose_as_is1( fd, pose, residue_set, missing, options );
	io::pose_from_sfr::PoseFromSFRBuilder builder( residue_set.get_self_ptr(), options );
	builder.build_pose( *sfr, pose );
	missing = builder.missing_atoms();
	build_pose_as_is2( sfr, pose, residue_set, missing, options );
}

void build_pose_as_is2(
	io::StructFileRepCOP /*fd*/,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	id::AtomID_Mask & missing,
	ImportPoseOptions const & options
)
{
	using namespace chemical;
	using namespace conformation;

	core::pose::PDBInfoOP pdb_info( new core::pose::PDBInfo(*pose.pdb_info()) );

	if ( !options.skip_set_reasonable_fold_tree() ) {
		set_reasonable_fold_tree( pose );
	}

	// optimize H if using a full-atom residue type set, and no_optH is not specified
	if ( residue_set.mode() == FULL_ATOM_t ) {
		//if pack_missing_density specified, repack residues w/ missing density
		if ( options.pack_missing_sidechains()  && ! options.membrane() ) {
			pack::pack_missing_sidechains( pose, missing );
		}
		// optimize H if using a fullatom residue type set, and no_optH is not specified
		if ( !options.no_optH() ) {
			pack::optimize_H_and_notify( pose, missing );
		}
	}


	// If pose contains carbohydrate residues, assure that their virtual atoms have the correct coordinates.
	if ( pose.conformation().contains_carbohydrate_residues() ) {
		for ( uint i = 1, n_residues = pose.size(); i <= n_residues; ++i ) {
			ResidueType const & res_type = pose.residue_type(i);
			if ( res_type.is_carbohydrate() ) {
				pose::carbohydrates::align_virtual_atoms_in_carbohydrate_residue(pose, i);
			}
		}
	}

	// If the user has set appropriate flags, check whether the pose contains metal ions,
	// and automatically set up covalent bonds and constraints to them.
	if ( options.set_up_metal_bonds() ) {
		core::util::auto_setup_all_metal_bonds(pose, options.metal_bond_LJ_multiplier(), true);
		if ( options.set_up_metal_constraints() ) {
			core::util::auto_setup_all_metal_constraints( pose, options.metal_bond_dist_constraint_multiplier(), options.metal_bond_angle_constraint_multiplier() );
		}
	}

	pose.pdb_info( pdb_info );
}

} // namespace import_pose
} // namespace core
