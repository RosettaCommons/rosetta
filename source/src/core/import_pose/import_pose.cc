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
#include <core/import_pose/FullModelPoseBuilder.hh>
#include <core/import_pose/import_pose_options.hh>

// Project headers
#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>

#include <core/util/metalloproteins_util.hh>

#include <core/pack/pack_missing_sidechains.hh>
#include <core/pack/optimizeH.hh>

#include <core/io/pdb/pdb_reader.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/pose_from_sfr/PoseFromSFRBuilder.hh>
#include <core/io/mmcif/cif_reader.hh>
#include <core/io/mmtf/mmtf_reader.hh>
#include <core/io/StructFileRep.hh>


#include <core/fragment/rna/RNA_MatchType.hh>

#include <core/scoring/rna/RNA_CentroidInfo.hh>

#include <core/id/TorsionID.hh>

// Basic headers
#include <basic/Tracer.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/magnesium.OptionKeys.gen.hh>
#include <basic/options/option.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

#include <ObjexxFCL/FArray2D.hh>
// External headers
#include <ObjexxFCL/string.functions.hh>
#include <cifparse/CifFile.h>
#include <cifparse/CifParserBase.h>

#include <tuple>
#include <boost/algorithm/string/predicate.hpp>

#include <core/pose/full_model_info/FullModelParameters.hh> // MANUAL IWYU

#include <utility/stream_util.hh> // AUTO IWYU For operator<<

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/details/helpers.hpp>
#endif

using CifFileOP = utility::pointer::shared_ptr<CifFile>;
using CifParserOP = utility::pointer::shared_ptr<CifParser>;

static basic::Tracer TR( "core.import_pose.util" );

namespace core {
namespace import_pose {

using namespace kinematics;

using core::Size;
using core::SSize;

using namespace core::io;
using namespace core::fragment::rna;
using namespace ObjexxFCL;
using namespace core::pose;
using namespace core::pose::full_model_info;

static basic::Tracer TR( "core.import_pose.import_pose" );

using utility::vector1;

std::ostream & operator<<( std::ostream & stream, FileType type ) {
	switch( type ) {
	case PDB_file :
		stream << "PDB";
		break;
	case CIF_file :
		stream << "mmCIF";
		break;
	case MMTF_file :
		stream << "MMTF";
		break;
	case SRLZ_file :
		stream << "SRLZ";
		break;
	default :
		stream << "UNKNOWN";
		break;
	}
	return stream;
}

/// @brief Given a filetype, return the string for the extension.
/// @details Extensions are in lowercase (e.g. "pdb", "cif", "mmtf"), etc.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
std::string
extension_from_filetype(
	FileType const filetype
) {
	switch( filetype ) {
	case PDB_file :
		return "pdb";
	case CIF_file :
		return "cif";
	case MMTF_file :
		return "mmtf";
	case SRLZ_file :
		return "srlz";
	default :
		utility_exit_with_message( "Error in core::import_pose::extension_from_filetype(): Invalid filetype provided!" );
		break;
	}
	return "ERROR"; //To keep compiler happy, though we should never reach here.
}

FileType
filetype_from_extension(
	std::string const & filename
) {
	utility::vector1< std::string > split = utility::string_split( filename, '.');

	std::string ext = split[split.size()-1];
	if ( ext == "gz" ) {
		ext = split[split.size()-2];
	}
	ext = utility::lower(ext);

	if ( ext == "pdb" ) {
		return PDB_file;
	} else if ( ext == "cif" ) {
		return CIF_file;
	} else if ( ext == "mmcif" ) {
		return CIF_file;
	} else if ( ext == "mmtf" ) {
		return MMTF_file;
	} else if ( ext == "srlz" ) {
		return SRLZ_file;
	} else {
		return Unknown_file;
	}
}

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
	pose.pdb_info()->parse_pdbinfo_labels( lines, pose );

	if ( (!read_fold_tree) && (!options.fold_tree_io()) ) return;

	for ( std::string const & line : lines ) {
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
	utility::vector1< std::string > lines( utility::split_by_newlines( contents_of_file ) );

	// The mmCIF format has a large number of initial underscores
	// (We put this test first as the "has ATOM record" test will pass on standard mmCIF files,
	// as they have a table that begins with "ATOM")
	core::Size n_initial_under(0);
	for ( std::string const& line: lines ) {
		if ( line.size() > 0 && line[0] == '_' ) {
			++n_initial_under;
		}
	}

	if ( n_initial_under > 1 ) {
		//See if this is a CIF file
		std::string diagnostics;
		CifFileOP cifFile( new CifFile );

		CifParserOP cifParser( new CifParser( cifFile.get() ) );

		cifParser->ParseString( contents_of_file, diagnostics );
		if ( diagnostics.empty() ) {
			return CIF_file;
		} else if ( TR.Debug.visible() ) {
			TR.Debug << "Attempted to read file as an mmCIF file. The mmCIF parser didn't like it, saying:" << std::endl;
			TR.Debug <<  diagnostics << std::endl;
		}
	}

	// See if this is a pdb file - Do we have proper ATOM/HETATM records?
	utility::vector1< io::pdb::Record> records( io::pdb::create_records_from_pdb_lines( lines ) );
	core::Size n_atom_records( 0 );
	for ( io::pdb::Record const & record: records ) {
		if ( record.count( "type" ) &&
				(record.at( "type" ).value == "ATOM  " ||
				record.at( "type" ).value == "HETATM" ) ) {
			++n_atom_records;
		}
	}
	if ( n_atom_records > 0 ) {
		return PDB_file;
	}

	return Unknown_file;
}


void
pose_from_file(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree,
	FileType file_type
) {

	//mmTF must go first as its contents are binary and opened by the mmTF machinery, IE not IZstream.
	if ( file_type == MMTF_file ||
			boost::algorithm::ends_with( utility::lower(filename), ".mmtf") ) {
		core::io::StructFileRepOP sfr( core::io::mmtf::create_sfr_from_mmtf_filename( filename, options ) );
		build_pose( sfr, pose, residue_set, options );
		return;
	}

	if ( boost::algorithm::ends_with( filename, ".mmtf.gz") ) {
		utility_exit_with_message(filename + " mmtf.gz reading is currently not supported.  You must first decompress you file(s)");
	}




	//JAB 7/23/19 - There can be paired ligands and receptors that can be concatonated on the fly using this code.
	//  This is done for large-scale screening jobs. For example: Some users do -s "some_ligand.pdb some_receptor.pdb"
	//  In the interest of not break existing use cases, I am leaving this functionality in after MUCH confusion.
	utility::vector1< std::string > filenames = utility::split(filename);
	std::string contents_of_file;
	for ( std::string const & fname : filenames ) {
		utility::io::izstream file( fname );
		if ( !file ) {
			TR.Error << "File: " << filename << " not found!" << std::endl;
			utility_exit_with_message( "Cannot open file \"" + filename + "\"" );
		} else {
			TR.Debug << "read file: " << filename << std::endl;
			utility::slurp( file, contents_of_file );
		}
	}


	if ( file_type == Unknown_file && boost::algorithm::ends_with( filenames[1], ".srlz" ) ) {
		file_type = SRLZ_file;
	}

	if ( file_type == Unknown_file ) {
		file_type = determine_file_type( contents_of_file );
		TR << "File '" << filename << "' automatically determined to be of type " << file_type << " from contents." << std::endl;
	}

	if ( file_type == Unknown_file ) {
		// Attempt to use extension as a proxy. Note that for merged filenames, this will only trigger on the last filename extension
		file_type =  filetype_from_extension( filename );
		TR << "Determining filetype from extension. File '" << filename << "' is typed as " << file_type << std::endl;
	}

	if ( file_type == Unknown_file ) {
		utility_exit_with_message( "Cannot determine file type. Current supported types are: PDB, CIF, SRLZ, MMTF");
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
		if ( !diagnostics.empty() ) {
			TR.Warning << "mmCIF parser reports issues with file '" << filename << "' : " << std::endl;
			TR.Warning << diagnostics << std::endl;
		}
		io::StructFileRepOP sfr ( io::mmcif::create_sfr_from_cif_file_op( cifFile, options ) );
		// We need to get rid of the cif parser IMMEDIATELY because of the arcanity in
		// the library that there can only be one at once... even though you can only have one per
		// file.
		cifParser.reset();
		build_pose( sfr, pose, residue_set, options );
	} else if ( file_type == SRLZ_file ) {
#ifdef SERIALIZATION
		std::istringstream iss( contents_of_file );
		cereal::BinaryInputArchive arc( iss );
		arc( pose );
#else
		utility_exit_with_message( "Cannot deserialize pose with non-serialization builds. Please compile with extras=serialization.");
#endif
	}
}

void
pose_from_file(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	bool read_fold_tree,
	FileType type
)
{
	ImportPoseOptions options;
	pose_from_file(pose, residue_set, filename, options, read_fold_tree, type);
}



void
pose_from_file(
	core::pose::Pose & pose,
	std::string const & filename,
	bool read_fold_tree,
	FileType type
) {
	ImportPoseOptions options;
	pose_from_file( pose, filename, options, read_fold_tree, type);
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

	// If in:file:residue_type_set is nondefault, use that.
	if ( options.residue_type_set_mode() != FULL_ATOM_t ) {
		TR.Info << "Specified explicit residue type set mode; disregarding -in:file:centroid and -in:file:fullatom." << std::endl;
		core::import_pose::pose_from_file( pose, *pose.residue_type_set_for_pose( options.residue_type_set_mode() ), filename, options, read_fold_tree, type );
	} else {
		ResidueTypeSetCOP residue_set( options.centroid() ?
			pose.residue_type_set_for_pose( CENTROID_t ) :
			pose.residue_type_set_for_pose( FULL_ATOM_t )
		);

		core::import_pose::pose_from_file( pose, *residue_set, filename, options, read_fold_tree, type);
	}
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

	// If in:file:residue_type_set is nondefault, use that.
	if ( options.residue_type_set_mode() != FULL_ATOM_t ) {
		TR.Info << "Specified explicit residue type set mode; disregarding -in:file:centroid and -in:file:fullatom." << std::endl;
		return core::import_pose::poseOPs_from_files( *ChemicalManager::get_instance()->residue_type_set( string_from_type_set_mode( options.residue_type_set_mode() ) ), filenames, options, read_fold_tree, type );
	} else {
		ResidueTypeSetCOP residue_set( options.centroid() ?
			ChemicalManager::get_instance()->residue_type_set( CENTROID ) :
			ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )
		);

		return core::import_pose::poseOPs_from_files( *residue_set, filenames, options, read_fold_tree, type);
	}
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
	chemical::ResidueTypeSetCOP residue_set
		( pose.residue_type_set_for_pose( chemical::FULL_ATOM_t ) );

	pose_from_pdbstring(pose, pdbcontents, *residue_set, options, filename );
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
	read_additional_pdb_data( pdbcontents, pose, options, options.read_fold_tree() );
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

	// Given this fold tree, look for residues that are not connected in the FT
	// but yet are chemically bonded upper to lower. (Chief example: N-to-C cyclization.)
	// To these residues, add CUTPOINT_LOWER and CUTPOINT_UPPER variants.
	//
	// AMW TODO: the above; for now, just see if the upper and lower of each edge (kinda
	// assumes each one is threading N to C) need cutpoint variants (to be scored by chainbreak).
	//
	// AMW: just did something a HAIR fancier; now we check ALL starts and stops that might
	// be chemically bonded to each other. Note that there is no danger of accidentally
	// melding termini together; we are after 'fixup termini' so only things that
	// have already been recognized as Weird Termini (due to, e.g., link_map() presence) are
	// gonna be eligible. The reason we check any starts and stops is because it's possible that
	// you have a slightly messed up chain with a chain-break in the middle, but it's cyclized
	// overall -- despite having been chopped into two Edges.
	for ( auto const &edge1 : pose.fold_tree() ) {
		for ( auto const &edge2 : pose.fold_tree() ) {
			if ( !edge1.is_polymer() || !edge2.is_polymer() ) continue;

			if ( !pose.residue_type( edge1.start() ).is_polymer() ) continue;
			if ( !pose.residue_type( edge2.stop()  ).is_polymer() ) continue;

			if ( pose.residue_type( edge1.start() ).lower_connect_id() == 0 ) continue;
			if ( pose.residue_type( edge2.stop()  ).upper_connect_id() == 0 ) continue;

			// AMW: 2KTQ fails here for reasons unclear; it seems to find a cutpoint
			// to be connected between two very distant residues.
			id::AtomID lower( pose.residue_type( edge1.start() ).lower_connect_atom(), edge1.start() );
			id::AtomID upper( pose.residue_type( edge2.stop()  ).upper_connect_atom(),  edge2.stop() );
			if ( pose.conformation().atoms_are_bonded( upper, lower ) ) {
				pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, edge1.start() );
				pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, edge2.stop()  );
			}
		}
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

			ResidueType const & res_type = pose.residue_type( i );
			if ( res_type.is_carbohydrate() ) {
				pose::carbohydrates::align_virtual_atoms_in_carbohydrate_residue( pose, i );
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


// Initialization from command line -- everything that doesn't require protocols
///////////////////////////////////////////////////////////////////////////////////////



pose::PoseOP
initialize_pose_and_other_poses_from_command_line( core::chemical::ResidueTypeSetCOP rsd_set ) {

	return initialize_pose_and_other_poses_from_options( rsd_set, basic::options::option );
}

PoseOP
initialize_pose_and_other_poses_from_options( core::chemical::ResidueTypeSetCOP rsd_set, utility::options::OptionCollection const & options ) {
	FullModelPoseBuilder builder;
	builder.set_options( options );
	builder.initialize_input_poses_from_options( rsd_set );
	builder.initialize_further_from_options();
	return builder.build();
}

PoseOP
initialize_pose_and_other_poses_from_options_and_input_poses(
	core::chemical::ResidueTypeSetCOP rsd_set,
	utility::options::OptionCollection const & options,
	utility::vector1< pose::PoseOP > & input_poses
) {
	// AMW: This function predates the FullModelPoseBuilder (developed above) but I think
	// it is still reasonable to use it, as we may simply set its input_poses manually.

	using namespace basic::options::OptionKeys;
	using namespace core::io::silent;
	using namespace core::pose;
	using namespace utility;

	std::tuple< vector1< Size >, vector1< char >, vector1< std::string > > const & input_resnum_and_chain_and_segid = options[ in::file::input_res ].resnum_and_chain();
	vector1< Size > const & input_res_list = std::get<0>( input_resnum_and_chain_and_segid );
	if ( input_res_list.size() ) {
		vector1< char > input_chain_list = std::get<1>( input_resnum_and_chain_and_segid );
		vector1< std::string > input_segid_list = std::get<2>( input_resnum_and_chain_and_segid );
		Size input_res_count = 0;
		for ( Size n = 1; n <= input_poses.size(); n++ ) {
			Pose & pose = *input_poses[ n ];
			PDBInfoOP pdb_info( new PDBInfo( pose ) );
			vector1< Size > input_res_for_pose;
			vector1< char > input_chain_for_pose;
			vector1< std::string > input_segid_for_pose;
			for ( Size k = 1; k <= pose.size(); k++ ) {
				input_res_count++;
				runtime_assert( input_res_count <= input_res_list.size() );
				Size const & number_in_full_model = input_res_list[ input_res_count ];
				input_res_for_pose.push_back( number_in_full_model );
				input_chain_for_pose.push_back( input_chain_list[ input_res_count ] );
				input_segid_for_pose.push_back( input_segid_list[ input_res_count ] );
			}
			pdb_info->set_numbering( input_res_for_pose );
			pdb_info->set_chains(    input_chain_for_pose );
			pdb_info->set_segment_ids(    input_segid_for_pose );
			pose.pdb_info( pdb_info );
		}
		runtime_assert( input_res_count == input_res_list.size() );
	}

	if ( input_poses.size() == 0 ) input_poses.push_back( utility::pointer::make_shared< Pose >() ); // just a blank pose for now.

	if ( options[ full_model::other_poses ].user() ) {
		get_other_poses( input_poses, options[ full_model::other_poses ](), rsd_set );
	}

	/*
	fill_full_model_info_from_options( input_poses, options );  //FullModelInfo (minimal object needed for add/delete)
	return input_poses[1];
	*/
	FullModelPoseBuilder builder;
	builder.set_options( options );
	builder.set_input_poses( input_poses );
	builder.initialize_further_from_options();
	return builder.build();
}

///////////////////////////////////////////////////////////////////////////////////////
void
get_other_poses( utility::vector1< pose::PoseOP > & other_poses,
	utility::vector1< std::string > const & other_files,
	core::chemical::ResidueTypeSetCOP rsd_set ) {

	for ( Size n = 1; n <= other_files.size(); n++ ) {
		other_poses.push_back( get_pdb_and_cleanup( other_files[ n ], rsd_set ) );
	}
}

//////////////////////////////////////////////////////////////////////////////////////
core::pose::PoseOP
get_pdb_with_full_model_info( std::string const & input_file,
	core::chemical::ResidueTypeSetCOP rsd_set ) {

	core::pose::PoseOP pose = get_pdb_and_cleanup( input_file, rsd_set );

	// a bit wasteful, since this checks command-line options over and over again, but hey this works.
	//  fill_full_model_info_from_command_line( *pose );

	return pose;
}

//////////////////////////////////////////////////////////////////////////////////////
// might be better to move these into core (e.g., core::pose::full_model_info ),
// or into a new protocols/full_model_setup/ directory.
core::pose::PoseOP
get_pdb_and_cleanup( std::string const & input_file )
{
	using namespace core::chemical;
	return get_pdb_and_cleanup( input_file, ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
}

core::pose::PoseOP
get_pdb_and_cleanup( std::string const & input_file,
	core::chemical::ResidueTypeSetCOP rsd_set )
{
	using namespace core::pose;
	PoseOP input_pose( new Pose );
	core::chemical::ResidueTypeSetCOP rsd_set_op( rsd_set );
	import_pose::pose_from_file( *input_pose, *rsd_set_op, input_file , core::import_pose::PDB_file);
	tag_into_pose( *input_pose, input_file );
	cleanup( *input_pose );
	core::pose::full_model_info::make_sure_full_model_info_is_setup( *input_pose );
	return input_pose;
}


//////////////////////////////////////////////////////////////////////////////////////
// currently have stuff we need for RNA... put any protein cleanup here too.
void
cleanup( pose::Pose & pose, bool const force_cut_at_rna_chainbreak /* = false */  ) {
	pose::rna::figure_out_reasonable_rna_fold_tree( pose, force_cut_at_rna_chainbreak);
	pose::rna::virtualize_5prime_phosphates( pose );
	pose.conformation().detect_disulfides();
}


///////////////////////////////////////////////////////////////////////////////////////
FullModelParametersOP
get_sequence_information(
	std::string const & fasta_file,
	vector1< Size > & cutpoint_open_in_full_model,
	bool const append_virtual /*=false*/ )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	vector1< core::sequence::SequenceOP > fasta_sequences = core::sequence::read_fasta_file( fasta_file );

	// calebgeniesse: setup for edensity scoring, if map is provided via cmd-line
	if ( option[ edensity::mapfile ].user() || append_virtual ) {
		// update fasta_sequences accordingly
		fasta_sequences.push_back( core::sequence::SequenceOP( new core::sequence::Sequence( "X" /*seq*/, " z:1" /*id*/) ) );
		//Size idx = fasta_sequences.size();
		//fasta_sequences[idx]->append_char('X');
		//fasta_sequences[idx]->id( fasta_sequences[idx]->id() + " z:1" );
	}

	look_for_dna( fasta_sequences );
	std::map< Size, std::string > non_standard_residue_map  = core::sequence::parse_out_non_standard_residues( fasta_sequences /*will reduce to one-letter*/ );
	if ( option[ magnesium::hydrate ]() ) setup_water_bank_for_magnesiums( non_standard_residue_map, fasta_sequences );
	std::string const desired_sequence           = core::sequence::get_concatenated_sequence( fasta_sequences );

	FullModelParametersOP full_model_parameters( new FullModelParameters( desired_sequence ) );
	vector1< char > conventional_chains;
	vector1< std::string > conventional_segids;
	vector1< int  > conventional_numbering;
	core::sequence::get_conventional_chains_and_numbering( fasta_sequences, conventional_chains, conventional_numbering, conventional_segids );

	full_model_parameters->set_conventional_numbering( conventional_numbering );
	full_model_parameters->set_conventional_segids( conventional_segids );
	full_model_parameters->set_conventional_chains( conventional_chains );

	full_model_parameters->set_non_standard_residue_map( non_standard_residue_map );
	cutpoint_open_in_full_model  = get_cutpoints( fasta_sequences, non_standard_residue_map,
		conventional_chains, conventional_numbering, conventional_segids );
	return full_model_parameters;
}

/////////////////////////////////////////////////////////////////////////////
//calebgeniesse: setup for density scoring
void
setup_for_density_scoring( core::pose::Pose & pose ) {
	// add virtual root for density scoring
	TR << "Adding virtual residue as root" << std::endl;
	pose::addVirtualResAsRoot( pose );
	// Note from rhiju -- should we update full_model_info & full_model_parameters?
	// I put a new function in core/pose/full_model_info/util.hh called append_virtual_residue_to_full_model_info().
}


///////////////////////////////////////////////////////////////////////////////////////
vector1< Size >
get_cutpoints_from_numbering( vector1< core::sequence::SequenceCOP > const & fasta_sequences,
	vector1< char > const & conventional_chains,
	vector1< std::string > const & conventional_segids,
	vector1< int  > const & conventional_numbering ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	vector1< Size > cutpoints;
	Size ntot( 0 );
	// explicit chain boundaries
	for ( Size n = 1; n <= fasta_sequences.size(); n++ ) {
		std::string const sequence = core::pose::rna::remove_bracketed( fasta_sequences[ n ]->sequence() );
		vector1< Size > const & spacer_pos = fasta_sequences[ n ]->spacer_positions();
		for ( Size q = 1; q <= spacer_pos.size(); q++ ) cutpoints.push_back( ntot + spacer_pos[ q ] );
		ntot += sequence.size();
		if ( n != fasta_sequences.size() /*very end is not 'cutpoint'*/ ) cutpoints.push_back( ntot );
	}
	runtime_assert( ntot == conventional_chains.size() );
	runtime_assert( ntot == conventional_numbering.size() );
	for ( Size k = 1; k < ntot; k++ ) {
		if ( cutpoints.has_value( k ) ) continue;
		if ( conventional_chains[ k ] != conventional_chains[ k+1 ] ) {
			cutpoints.push_back( k );
		} else if ( conventional_segids[ k ] != conventional_segids[ k+1 ] ) {
			cutpoints.push_back( k );
		} else { // break within chain.
			if ( conventional_numbering[ k+1 ] > conventional_numbering[ k ] + 1  ) {
				if ( !option[ full_model::allow_jump_in_numbering ]() ) {
					TR << TR.Red << "There appears to be a break in numbering at " << conventional_numbering[ k ] << " so adding to cutpoint list. If that is wrong, use flag -allow_jump_in_numbering." << TR.Reset << std::endl;
					cutpoints.push_back( k );
				}
			}
		}
	}

	std::sort( cutpoints.begin(), cutpoints.end() );
	return cutpoints;
}


/////////////////////////////////////////////////////////////////////////////////////////
// Kind of a hack. Would be better to somehow figure out if residue is not polymeric --
//  perhaps try to instantiate with ResidueTypeSet?
void
get_extra_cutpoints_from_names( Size const nres,
	vector1< Size > & cutpoints,
	std::map< Size, std::string > const & non_standard_residue_map )
{
	for ( Size n = 1; n <= nres; n++ ) {
		if ( cutpoints.has_value( n ) ) continue;
		if ( !stepwise_addable_residue( n, non_standard_residue_map ) ) {
			cutpoints.push_back( n );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////
vector1< Size >
get_cutpoints( vector1< core::sequence::SequenceCOP > const & fasta_sequences,
	std::map< Size, std::string > const & non_standard_residue_map,
	vector1< char > const & conventional_chains,
	vector1< int  > const & conventional_numbering,
	vector1< std::string > const & conventional_segids ) {
	vector1< Size > cutpoints = get_cutpoints_from_numbering( fasta_sequences, conventional_chains, conventional_segids, conventional_numbering );
	get_extra_cutpoints_from_names( conventional_numbering.size(), cutpoints, non_standard_residue_map );
	return cutpoints;
}


///////////////////////////////////////////////////////////////////////////////////////
// Need a 'water bank' to hold up to 6 waters for each Mg(2+) during stepwise modeling
//  with water hydration. Some of these waters will get virtualized, but its important
//  that number of residues stays constant.
//
// Following is based on an assumption that the *only* waters to be modeled in stepwise
//  would be ones associated with Mg(2+)... obviously that is not very general.
//
// It might be smarter to look inside the PDBs (esp. any native PDBs) to count
//   waters associated with each Mg(2+), one by one, and supplement water bank to
//   ensure 6 waters per Mg(2+).
void
setup_water_bank_for_magnesiums( std::map< Size, std::string > & non_standard_residue_map,
	vector1< core::sequence::SequenceOP > & fasta_sequences ) {
	using namespace core::sequence;

	// how many magnesiums are there? how many waters are there?
	Size offset( 0 ), num_magnesiums( 0 ), num_waters( 0 );
	for ( Size n = 1; n <= fasta_sequences.size(); n++ ) {
		std::string const sequence = fasta_sequences[n]->sequence();
		for ( Size i = 1; i <= sequence.size(); i++ ) {
			if ( sequence[ i - 1 ] == 'Z' ) {
				std::map< Size, std::string >::const_iterator it = non_standard_residue_map.find( offset+i );
				if ( it != non_standard_residue_map.end() && it->second == "MG" )  num_magnesiums++;
			} else if ( sequence[ i - 1 ] == 'w' ) {
				std::map< Size, std::string >::const_iterator it = non_standard_residue_map.find( offset+i );
				if ( it != non_standard_residue_map.end() && it->second == "HOH" )  num_waters++;
			}
		}
		offset += sequence.size();
	}

	runtime_assert( num_waters <= 6 * num_magnesiums ); // currently do not model any non-Mg(2+) waters.

	// are there at least six waters for each magnesium? if not, let's put them in.
	// Use chain 'w' semi-arbitrarily, and numbering from 1001 onward. -- this could be a disaster.
	std::string sequence;
	std::stringstream id;
	Size extra_water_number( 1000 );
	for ( Size n = num_waters; n <= 6 * num_magnesiums; n++ ) {
		sequence += 'w';
		offset++; // number of residues.
		extra_water_number++; // 1001, 1002, ...
		non_standard_residue_map[ offset ] = "HOH";
		id << " w:" << extra_water_number;
	}

	fasta_sequences.push_back( utility::pointer::make_shared< Sequence >( sequence, id.str() ) );
}

void
look_for_dna( vector1< core::sequence::SequenceOP > & fasta_sequences )
{
	using namespace core::sequence;
	std::map<char, std::string> rna_name = { {'a',"RAD"}, {'c',"RCY"}, {'g',"RGU"}, {'t',"5MU"} };
	for ( Size n = 1; n <= fasta_sequences.size(); n++ ) {
		auto fasta_sequence = fasta_sequences[ n ];
		auto id = fasta_sequences[ n ]->id();
		std::stringstream ss( fasta_sequences[n]->id() );
		bool is_DNA( false );
		while ( ss.good() ) {
			std::string tag;
			ss >> tag;
			if ( tag == "DNA" ) {
				is_DNA = true; break;
			}
		}
		auto sequence = fasta_sequences[ n ]->sequence();
		utility::vector1< std::string > fullname_list; // a vector of non-standard full names
		std::vector< Size > oneletter_to_fullname_index; // for each one-letter sequence, zero means no fullname given
		std::string one_letter_sequence;
		parse_sequence( sequence, fullname_list, oneletter_to_fullname_index, one_letter_sequence );
		if ( !is_DNA ) {
			for ( Size k = 1; k <= one_letter_sequence.size(); k++ ) {
				if ( one_letter_sequence[ k-1 ] == 't' ) {
					Size const pos = oneletter_to_fullname_index[ k - 1 ];
					if ( pos == 0 || fullname_list[ pos ].substr( 3, 14 ) != ":deoxy_O2prime" ) {
						utility_exit_with_message( "Seeing a 't' in sequence -- if you want this to be DNA, include the tag 'DNA' in the FASTA id for this sequence. if you want this to be RNA, use X[5MU] as the character." );
					}
				}
			}
		}
		if ( is_DNA ) {
			auto fasta_sequence_new = fasta_sequence;
			vector1< std::string > fullname_list_new;
			for ( Size k = 1; k <= oneletter_to_fullname_index.size(); k++ ) {
				Size const pos = oneletter_to_fullname_index[ k - 1 ];
				std::string fullname_new;
				if ( pos > 0 ) {
					fullname_new = fullname_list[ pos ];
				} else {
					char q = one_letter_sequence[ k - 1 ];
					if ( !rna_name.count( q ) ) utility_exit_with_message( "The character "+std::string(&q)+" is not a recognized DNA nucleotide. Use a,c,g,t." );
					fullname_new = rna_name[ q ];
				}
				if ( fullname_new.substr(3,14) != ":deoxy_O2prime" ) {
					fullname_new = fullname_new.substr(0,3) + ":deoxy_O2prime" + fullname_new.substr(3,std::string::npos);
				}
				fullname_list_new.push_back( fullname_new );
			}
			std::string sequence_new;
			for ( Size k = 1; k <= one_letter_sequence.size(); k++ ) {
				sequence_new.push_back( one_letter_sequence[ k-1 ] );
				sequence_new.push_back( '[' );
				sequence_new.append( fullname_list_new[ k ] );
				sequence_new.push_back( ']' );
			}
			fasta_sequences[ n ] = utility::pointer::make_shared< Sequence >( sequence_new, id );
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////
// In principle, all the parameters being sent in to this function could be
// encapsulated into a FullModelParameters object. However, its important to be conscious
// of what is const and what can be updated at the level of the individual variables.
void
setup_fold_trees( vector1< Pose * > & pose_pointers,
	vector1< Size > & cutpoint_open_in_full_model /* can be updated here*/,
	vector1< Size > & fixed_domain_map /* domain colors can be updated here */,
	vector1< Size > const & cutpoint_closed,
	vector1< Size > const & extra_minimize_res,
	vector1< Size > const & extra_minimize_jump_res,
	vector1< Size > const & sample_res,
	vector1< Size > const & working_res,
	vector1< Size > const & jump_res,
	vector1< Size > const & preferred_root_res,
	vector1< Size > const & virtual_sugar_res,
	FullModelParameters const & full_model_parameters,
	vector1< vector1< Size > > const & pose_res_lists ) {

	for ( Size n = 1; n <= pose_pointers.size(); n++ ) {
		Pose & pose = *pose_pointers[n];
		vector1< Size > const & res_list = pose_res_lists[ n ];
		for ( Size i = 1; i < pose.size(); i++ ) {
			if ( (res_list[ i+1 ] > res_list[ i ] + 1) && !pose.fold_tree().is_cutpoint(i) ) {
				put_in_cutpoint( pose, i );
			}
			if ( cutpoint_open_in_full_model.has_value( res_list[ i ]) ) continue;
			if ( (res_list[ i+1 ] == res_list[ i ] + 1) &&
					pose.fold_tree().is_cutpoint( i ) &&
					!pose.residue_type( i   ).has_variant_type( core::chemical::CUTPOINT_LOWER ) &&
					!pose.residue_type( i+1 ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
				TR << TR.Red << "There appears to be a strand boundary at " << res_list[ i ] << " so adding to cutpoint_in_full_model." << TR.Reset << std::endl;
				cutpoint_open_in_full_model.push_back( res_list[ i ] ); continue;
			}
			if ( pose.residue_type( i ).is_RNA() != pose.residue_type( i+1 ).is_RNA()
					&& !pose.residue_type( i+1 ).is_virtual_residue() ) {
				cutpoint_open_in_full_model.push_back( res_list[ i ] ); continue;
			}
			if ( pose.residue_type( i ).is_protein() != pose.residue_type( i+1 ).is_protein()
					&& !pose.residue_type( i+1 ).is_virtual_residue() ) {
				cutpoint_open_in_full_model.push_back( res_list[ i ] ); continue;
			}
		}

		add_cutpoint_closed( pose, res_list, cutpoint_closed );

		update_pose_fold_tree( pose,
			res_list,
			extra_minimize_res, sample_res, jump_res,
			full_model_parameters /* set preferences for chain connections*/ );
		add_virtual_sugar_res( pose, res_list,
			virtual_sugar_res ); // this was for checks -- no longer in use?

		update_fixed_domain_from_extra_minimize_jump_res( fixed_domain_map, pose, res_list, extra_minimize_jump_res );

		vector1< Size> root_partition_res;
		for ( Size nn = 1; nn <= pose.size(); ++nn ) root_partition_res.push_back( nn );
		core::pose::reroot( pose, root_partition_res, res_list, preferred_root_res, fixed_domain_map,
			cutpoint_open_in_full_model, working_res );

	}
}

////////////////////////////////////////////////////////////////////////////////////
void
update_pose_fold_tree( pose::Pose & pose,
	vector1< Size > const & res_list,
	vector1< Size > const & extra_min_res,
	vector1< Size > const & sample_res,
	vector1< Size > const & jump_res,
	core::pose::full_model_info::FullModelParameters const & full_model_parameters ) {

	if ( pose.size() == 0 ) return;

	vector1< vector1< Size > > all_res_in_chain, all_fixed_res_in_chain;
	vector1< Size > moveable_res = sample_res;
	for ( Size n = 1; n <= extra_min_res.size(); n++ ) moveable_res.push_back( extra_min_res[n] );
	define_chains( pose, all_res_in_chain, all_fixed_res_in_chain, res_list, moveable_res );
	Size nchains = all_res_in_chain.size();

	vector1< Size > jump_partners1, jump_partners2, cuts, blank_vector;
	vector1< pair< Size, Size > > chain_connections;
	setup_user_defined_jumps( jump_res, jump_partners1, jump_partners2,
		chain_connections, res_list, all_res_in_chain );
	runtime_assert( jump_partners1.size() < nchains );

	// needed to determine preferences for chain connections
	std::tuple< vector1< int >, vector1< char >, vector1< std::string > > const resnum_and_chain_in_pose = full_model_parameters.full_to_conventional_resnum_and_chain_and_segid( res_list );

	// choose fixed_res as jump partners first.
	setup_jumps( pose, jump_partners1, jump_partners2, chain_connections, all_fixed_res_in_chain, resnum_and_chain_in_pose );
	// use moving res if necessary
	setup_jumps( pose, jump_partners1, jump_partners2, chain_connections, all_res_in_chain, resnum_and_chain_in_pose );
	runtime_assert( jump_partners1.size() == (nchains - 1) );

	for ( Size n = 1; n < nchains; n++ ) cuts.push_back( all_res_in_chain[n][ all_res_in_chain[n].size() ] );
	FoldTree f = get_tree( pose, cuts, jump_partners1, jump_partners2 );
	pose.fold_tree( f );
}

///////////////////////////////////////////////////////////////////////////////////////
void
define_chains( pose::Pose const & pose,
	vector1< vector1< Size > > & all_res_in_chain,
	vector1< vector1< Size > > & all_fixed_res_in_chain,
	vector1< Size > const & res_list,
	vector1< Size > const & moveable_res ) {

	Size chain_start( 1 ), chain_end( 0 );
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( !pose.fold_tree().is_cutpoint( n ) &&
				n < pose.size() ) continue;
		chain_end = n;
		vector1< Size > res_in_chain, fixed_res_in_chain;
		for ( Size i = chain_start; i <= chain_end; i++ ) {
			res_in_chain.push_back( i );
			if ( !moveable_res.has_value( res_list[ i ] ) ) fixed_res_in_chain.push_back( i );
		}
		all_res_in_chain.push_back( res_in_chain );
		all_fixed_res_in_chain.push_back( fixed_res_in_chain );
		chain_start = chain_end + 1;
	}
}


////////////////////////////////////////////////////////////////////////////////
void
setup_user_defined_jumps( vector1< Size > const & jump_res,
	vector1< Size > & jump_partners1,
	vector1< Size > & jump_partners2,
	vector1< pair< Size, Size > > & chain_connections,
	vector1< Size > const & res_list,
	vector1< vector1< Size > > const & all_res_in_chain ) {

	vector1< Size > connection_domains = get_connection_domains( chain_connections, all_res_in_chain.size() );
	// Figure out jump res.
	for ( Size n = 1; n <= jump_res.size()/2; n++ ) {
		if ( res_list.has_value( jump_res[ 2*n - 1 ] )  &&
				res_list.has_value( jump_res[ 2*n ] ) ) {
			Size const i = res_list.index( jump_res[ 2*n - 1 ] );
			Size const j = res_list.index( jump_res[ 2*n ] );
			jump_partners1.push_back( i );
			jump_partners2.push_back( j );
			Size const chain_i( get_chain( i, all_res_in_chain ) );
			Size const chain_j( get_chain( j, all_res_in_chain ) );
			runtime_assert( connection_domains[ chain_i ] != connection_domains[ chain_j ] );
			chain_connections.push_back( make_pair( chain_i, chain_j ) );
			connection_domains = get_connection_domains( chain_connections, all_res_in_chain.size() );
		}
	}

}

////////////////////////////////////////////////////////////////////////////////
Size
get_chain( Size const i, vector1< vector1< Size > > const & all_res_in_chain ) {
	for ( Size n = 1; n <= all_res_in_chain.size(); n++ ) {
		if ( all_res_in_chain[ n ].has_value( i ) ) return n;
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////////
void
setup_jumps( core::pose::Pose const & pose,
	vector1< Size > & jump_partners1,
	vector1< Size > & jump_partners2,
	vector1< pair< Size, Size > > & chain_connections,
	vector1< vector1< Size > > const & all_res_in_chain,
	std::tuple< vector1< int >, vector1< char >, vector1< std::string > > const & resnum_and_chain_and_segid_in_pose ) {

	using namespace std;

	Size const num_chains = all_res_in_chain.size();
	if ( jump_partners1.size() == num_chains - 1 ) return;

	// order of preference for chain connections:
	//
	// 1. start with chains that have maximum number of contacts, then
	// 2. choose chains that are *closest in sequence*, as defined
	//      by end of chain i and start of chain j.
	// 3. break remaining ties based on chain order.
	//
	// For jump connection residues, prefer residues that have the most contacts.
	//
	vector1< int  > const & conventional_numbering_in_pose = std::get< 0 >( resnum_and_chain_and_segid_in_pose );
	vector1< char > const & conventional_chains_in_pose = std::get< 1 >( resnum_and_chain_and_segid_in_pose );
	vector1< std::string > const & conventional_segids_in_pose = std::get< 2 >( resnum_and_chain_and_segid_in_pose );

	// Data structure set up for sorting...
	//
	// chain_pair = ( -num_contacts, sequence_separation, ( ( chain_idx i, chain_idx j ), ( res in chain i, res in chain j ) ) )
	//
	utility::vector1< pair< int, pair< Size, pair< pair< Size, Size >, pair< Size, Size > > > > > chain_pairs;
	Size const max_seq_separation = max( conventional_numbering_in_pose ) - min( conventional_numbering_in_pose ) + 1;
	for ( Size i = 1; i < num_chains; i++ ) {
		vector1< Size > const & fixed_res_in_chain_i = all_res_in_chain[ i ];
		if ( fixed_res_in_chain_i.size() == 0 ) continue;

		for ( Size j = i+1; j <= num_chains; j++ ) {
			pair< Size, Size > jump_res_pair( 0, 0 );
			vector1< Size > const & fixed_res_in_chain_j = all_res_in_chain[ j ];
			if ( fixed_res_in_chain_j.size() == 0 ) continue;

			// num_contacts
			static Distance const CONTACT_DIST_CUTOFF( 4.0 );
			int num_contacts( 0 ), boundary_pair_contacts( 0 );
			vector1< pair< int, pair< Size, Size > > > num_contacts_pairwise;
			for ( Size m = 1; m <= fixed_res_in_chain_i.size(); m++ ) {
				core::conformation::Residue rsd_i = pose.residue( fixed_res_in_chain_i[ m ] );
				for ( Size n = 1; n <= fixed_res_in_chain_j.size(); n++ ) {
					core::conformation::Residue rsd_j = pose.residue( fixed_res_in_chain_j[ n ] );
					Size count( 0 );
					for ( Size mm = 1; mm <= rsd_i.nheavyatoms(); mm++ ) {
						for ( Size nn = 1; nn <= rsd_j.nheavyatoms(); nn++ ) {
							if ( ( rsd_i.xyz( mm ) - rsd_j.xyz( nn ) ).length() < CONTACT_DIST_CUTOFF ) {
								count++;
							}
						}
					}
					// this used to be a bug -- should save -count to *maximize* number of contacts
					num_contacts_pairwise.push_back( make_pair( -count, make_pair( fixed_res_in_chain_i[m],
						fixed_res_in_chain_j[n] ) ) );
					num_contacts += count;
					if ( m == fixed_res_in_chain_i.size() && n == 1 ) boundary_pair_contacts = count;
				}
			}

			// find residue pair that best represents this chain-to-chain connection
			if ( boundary_pair_contacts > 0 && ( j == i+1 ) ) {
				// use the 'boundary pair', which is the default in figure_out_reasonable_rna_fold_tree().
				// End of first chain, beginning of second chain
				jump_res_pair = make_pair( fixed_res_in_chain_i[ fixed_res_in_chain_i.size() ], fixed_res_in_chain_j[ 1 ] );
			} else {
				// sort based on number of contacts
				std::sort( num_contacts_pairwise.begin(), num_contacts_pairwise.end() );
				if ( num_contacts_pairwise.size() > 0 && num_contacts > 0 ) {
					jump_res_pair = num_contacts_pairwise[ 1 ].second;
				}
			}

			// figure out sequence separation
			Size sequence_separation( max_seq_separation );
			Size const chain_i_end   = fixed_res_in_chain_i[ fixed_res_in_chain_i.size() ];
			Size const chain_j_begin = fixed_res_in_chain_j[ 1 ];
			if ( conventional_chains_in_pose[ chain_i_end ] == conventional_chains_in_pose[ chain_j_begin ]
					&& conventional_segids_in_pose[ chain_i_end ] == conventional_segids_in_pose[ chain_j_begin ] ) {
				//runtime_assert( conventional_numbering_in_pose[ chain_j_begin ] > conventional_numbering_in_pose[ chain_i_end ] );
				sequence_separation = std::abs( conventional_numbering_in_pose[ chain_j_begin ] - conventional_numbering_in_pose[ chain_i_end ] );
			}
			if ( jump_res_pair == make_pair( Size(0), Size(0) ) ) jump_res_pair = make_pair( chain_i_end, chain_j_begin);

			chain_pairs.push_back( make_pair( -num_contacts /*want to maximize contact*/,
				make_pair( sequence_separation,
				make_pair( make_pair( i, j ), jump_res_pair ) ) ) );
		}
	}

	std::sort( chain_pairs.begin(), chain_pairs.end() );

	utility::vector1< Size > connection_domains = get_connection_domains( chain_connections, all_res_in_chain.size() );
	for ( Size q = 1; q <= chain_pairs.size(); q++ ) {
		Size const & i = chain_pairs[ q ].second.second.first.first;
		Size const & j = chain_pairs[ q ].second.second.first.second;
		if ( connection_domains[ i ] == connection_domains[ j ] ) continue;
		Size const & jump_res_i = chain_pairs[ q ].second.second.second.first;
		Size const & jump_res_j = chain_pairs[ q ].second.second.second.second;
		jump_partners1.push_back( jump_res_i );
		jump_partners2.push_back( jump_res_j );
		chain_connections.push_back( make_pair( i, j ) );
		connection_domains = get_connection_domains( chain_connections, all_res_in_chain.size() );
		if ( jump_partners1.size() == num_chains - 1 ) break;
	}

}

/////////////////////////////////////////////////////////////////////////////////
FoldTree
get_tree( pose::Pose const & pose,
	vector1< Size > const & cuts,
	vector1< Size > const & jump_partners1,
	vector1< Size > const & jump_partners2 ) {

	vector1< std::string > jump_atoms1, jump_atoms2;
	for ( Size n = 1; n <= jump_partners1.size(); n++ ) {
		jump_atoms1.push_back( chemical::rna::default_jump_atom( pose.residue_type( jump_partners1[n] ) ) );
		jump_atoms2.push_back( chemical::rna::default_jump_atom( pose.residue_type( jump_partners2[n] ) ) );
	}
	return get_tree( pose.size(), cuts, jump_partners1, jump_partners2, jump_atoms1, jump_atoms2 );
}


/////////////////////////////////////////////////////////////////////////////////
FoldTree
get_tree( Size const nres,
	vector1< Size > const & cuts,
	vector1< Size > const & jump_partners1,
	vector1< Size > const & jump_partners2,
	vector1< std::string > const & jump_atoms1,
	vector1< std::string > const & jump_atoms2 ) {

	Size const num_cuts = cuts.size();

	FoldTree f;
	ObjexxFCL::FArray2D< Size > jump_point_( 2, num_cuts );
	ObjexxFCL::FArray1D< Size > cuts_( num_cuts );
	for ( Size i = 1; i <= num_cuts; i++ ) {
		jump_point_( 1, i ) = std::min( jump_partners1[ i ], jump_partners2[ i ] );
		jump_point_( 2, i ) = std::max( jump_partners1[ i ], jump_partners2[ i ] );
		cuts_( i ) = cuts[ i ];
	}
	f.tree_from_jumps_and_cuts( nres, num_cuts, jump_point_, cuts_ );

	bool const KeepStubInResidue( true );
	for ( Size i = 1; i <= num_cuts; i++ ) {
		Size const n = f.jump_nr( jump_partners1[ i ], jump_partners2[ i ] );
		f.set_jump_atoms( n,
			jump_partners1[ i ], jump_atoms1[ i ],
			jump_partners2[ i ], jump_atoms2[ i ], KeepStubInResidue );
	}
	f.reassign_atoms_for_intra_residue_stubs(); // it seems silly that we need to do this separately.

	return f;
}

///////////////////////////////////////////////////////////////////////////////////////
void
update_fixed_domain_from_extra_minimize_jump_pairs( utility::vector1< Size > & fixed_domain,
	pose::Pose const & pose,
	utility::vector1< Size > const & res_list,
	utility::vector1< pair< Size, Size > > const & extra_minimize_jump_pairs ) {


	for ( Size i = 1; i <= extra_minimize_jump_pairs.size(); i++ ) {
		Size const res1_full = extra_minimize_jump_pairs[ i ].first;
		Size const res2_full = extra_minimize_jump_pairs[ i ].second;

		if ( !res_list.has_value( res1_full ) ) continue;
		if ( !res_list.has_value( res2_full ) ) continue;

		runtime_assert( fixed_domain[ res1_full ] > 0 );
		runtime_assert( fixed_domain[ res2_full ] > 0 );
		if ( fixed_domain[ res1_full ] != fixed_domain[ res2_full ] ) continue; // don't need to do anything.

		Size const res1 = res_list.index( res1_full );
		Size const res2 = res_list.index( res2_full );
		runtime_assert( pose.fold_tree().jump_exists( res1, res2 ) );
		Size jump_nr = pose.fold_tree().jump_nr( res1, res2 ); // need to check both orders!
		if ( jump_nr == 0 ) jump_nr = pose.fold_tree().jump_nr( res2, res1 );
		runtime_assert( jump_nr > 0 );

		utility::vector1< bool > partition_definition = pose.fold_tree().partition_by_jump( jump_nr );
		utility::vector1< Size > residues_to_update;
		Size const downstream_res = pose.fold_tree().downstream_jump_residue( jump_nr );
		for ( Size n = 1; n <= pose.size(); n++ ) {
			if ( partition_definition[ n ] == partition_definition[ downstream_res ] &&
					fixed_domain[ res_list[ n ] ] ==         fixed_domain[ res_list[ downstream_res ] ] ) {
				residues_to_update.push_back( n );
			}
		}
		runtime_assert( residues_to_update.size() > 0 );
		Size const new_domain = max( fixed_domain ) + 1;
		for ( Size n = 1; n <= residues_to_update.size(); n++ ) {
			fixed_domain[ res_list[ residues_to_update[ n ] ] ] = new_domain;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////
void
update_fixed_domain_from_extra_minimize_jump_res( vector1< Size > & fixed_domain,
	pose::Pose const & pose,
	vector1< Size > const & res_list,
	vector1< Size > const & extra_minimize_jump_res ) {
	utility::vector1< pair< Size, Size > > extra_minimize_jump_pairs;
	for ( Size n = 1; n <= extra_minimize_jump_res.size()/2; n++ ) {
		extra_minimize_jump_pairs.push_back( make_pair( extra_minimize_jump_res[ 2*n - 1 ],
			extra_minimize_jump_res[ 2*n     ] ) );
	}
	update_fixed_domain_from_extra_minimize_jump_pairs( fixed_domain, pose, res_list, extra_minimize_jump_pairs );
}

/////////////////////////////////////////////////////////////////////////////////////
void
add_cutpoint_closed( pose::Pose & pose,
	vector1< Size > const & res_list,
	vector1< Size > const & cutpoint_closed ) {
	using namespace core::pose::rna;
	for ( Size n = 1; n <= cutpoint_closed.size(); n++ ) {
		if ( !res_list.has_value( cutpoint_closed[ n ] ) ) continue;
		Size const i = res_list.index( cutpoint_closed[ n ] );
		// could be useful in general -- share this with TransientCutpointHandler?
		vector1< pair< core::id::TorsionID, Real > > const suite_torsion_info = get_suite_torsion_info( pose, i );
		put_in_cutpoint( pose, i );
		correctly_add_cutpoint_variants( pose, i );
		apply_suite_torsion_info( pose, suite_torsion_info );
	}
}


/////////////////////////////////////////////////////////////////////////////////////
void
put_in_cutpoint( pose::Pose & pose, Size const i ) {
	core::kinematics::FoldTree f = pose.fold_tree();
	Size const new_jump = f.new_jump( i, i+1, i );
	if ( pose.residue_type( i ).is_RNA() && pose.residue_type( i+1 ).is_RNA() ) {
		f.set_jump_atoms( new_jump,
			chemical::rna::default_jump_atom( pose.residue_type( f.upstream_jump_residue( new_jump ) ) ),
			chemical::rna::default_jump_atom( pose.residue_type( f.downstream_jump_residue( new_jump ) ) ) );
	}
	pose.fold_tree( f );
}

/////////////////////////////////////////////////////////////////////////////////////
void
add_virtual_sugar_res( pose::Pose & pose,
	vector1< Size > const & res_list,
	vector1< Size > const & virtual_sugar_res )
{
	for ( Size const res : virtual_sugar_res ) {
		if ( !res_list.has_value( res ) ) continue;
		Size const i = res_list.index( res );
		runtime_assert( i == 1 || pose.fold_tree().is_cutpoint( i - 1 ) );
		runtime_assert( i == pose.size() || pose.fold_tree().is_cutpoint( i ) );
		add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_RIBOSE, i );
	}
}

/////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
figure_out_working_res( utility::vector1< Size > const & input_domain_map,
	utility::vector1< Size > const & sample_res ) {
	vector1< Size > working_res;
	for ( Size n = 1; n <= input_domain_map.size(); n++ ) {
		if ( input_domain_map[ n ] == 0  && !sample_res.has_value( n ) ) continue;
		working_res.push_back( n );
	}
	return working_res;
}

/////////////////////////////////////////////////////////////////////////////////////
// 'default': sample any residues that are not inputted as fixed PDBs.
utility::vector1< Size >
figure_out_sample_res( utility::vector1< Size > const & input_domain_map,
	utility::vector1< Size > const & working_res ) {
	vector1< Size > sample_res;
	for ( Size n = 1; n <= input_domain_map.size(); n++ ) {
		if ( working_res.size() > 0 && !working_res.has_value( n ) ) continue;
		if ( input_domain_map[ n ] == 0  ) sample_res.push_back( n );
	}
	return sample_res;
}

/////////////////////////////////////////////////////////////////////////////////////
void
check_working_res( utility::vector1< Size > const & working_res,
	utility::vector1< Size > const & input_domain_map,
	utility::vector1< Size > const & sample_res ) {
	for ( Size i = 1; i <= sample_res.size(); i++ ) {
		runtime_assert( working_res.has_value( sample_res[ i ] ) );
	}
	for ( Size n = 1; n <= input_domain_map.size(); n++ ) {
		if ( input_domain_map[ n ] > 0 ) {
			if ( !working_res.has_value( n ) /* && !reference_pose */ ) {
				utility_exit_with_message( "Working res does not have input_domain_map residue "+ObjexxFCL::string_of(n) );
			}
			// then how can you sample something from input?
			//if ( sample_res.has_value( n ) ) utility_exit_with_message( "Sample res should not have "+ObjexxFCL::string_of(n) );
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// motif_mode:
//
//  extra_min_res -- minimize residues in input PDBs immediately neighboring loops to be built.
//  terminal_res  -- don't stack/pair loops on any cutpoint residues in input PDBs (except right next to loops)
//
// originally worked out in python (setup_stepwise_benchmark.py) by rhiju, 2014.
//
void
figure_out_motif_mode( utility::vector1< Size > & extra_min_res,
	utility::vector1< Size > & terminal_res,
	utility::vector1< Size > const & working_res,
	utility::vector1< Size > const & input_domain_map,
	utility::vector1< Size > const & cutpoint_open_in_full_model ) {
	// for double-checking work against any user-inputted extra_min_res & terminal_res
	utility::vector1< Size > const extra_min_res_save = extra_min_res;
	utility::vector1< Size > const terminal_res_save = terminal_res;
	extra_min_res.clear();
	terminal_res.clear();

	for ( Size m = 1; m <= input_domain_map.size(); m++ ) {
		if ( input_domain_map[ m ] == 0 ) continue;
		bool const right_before_chainbreak = ( m == input_domain_map.size() || cutpoint_open_in_full_model.has_value( m ) ||
			( working_res.size() > 0 && !working_res.has_value( m + 1 ) ) );
		bool const right_after_chainbreak  = ( m == 1                 || cutpoint_open_in_full_model.has_value( m - 1 ) ||
			( working_res.size() > 0 && !working_res.has_value( m - 1 )  ) );
		bool const prev_moving = !right_after_chainbreak && m > 1 &&
			( ( input_domain_map[ m - 1 ] == 0  && ( working_res.size() == 0 || working_res.has_value( m - 1 ) ) ) ||
			( input_domain_map[ m - 1 ] != input_domain_map[ m ] ) ) ;
		bool const next_moving = m < input_domain_map.size() && !right_before_chainbreak &&
			( ( input_domain_map[ m + 1 ] == 0 &&  ( working_res.size() == 0 || working_res.has_value( m + 1 ) ) ) ||
			( input_domain_map[ m + 1 ] != input_domain_map[ m ] ) );
		if ( ( right_after_chainbreak  && !right_before_chainbreak && !next_moving ) ||
				( right_before_chainbreak && !right_after_chainbreak  && !prev_moving ) ) {
			terminal_res.push_back( m );
		}
		if ( ( prev_moving && !next_moving && !right_before_chainbreak ) ||
				( next_moving && !prev_moving && !right_after_chainbreak ) ) {
			extra_min_res.push_back( m );
		}
	}

	//  if ( looks_like_reference_pose( input_domain_map ) ) return;
	if ( extra_min_res_save.size() > 0 ) {
		if ( extra_min_res != extra_min_res_save ) {
			TR << TR.Red << "Auto-computed EXTRA_MIN_RES " << extra_min_res <<      TR.Reset << std::endl;
			TR << TR.Red << "User-input    EXTRA_MIN_RES " << extra_min_res_save << TR.Reset << std::endl;
		}
		runtime_assert( extra_min_res_save.size() == 0 || extra_min_res == extra_min_res_save );
	}
	if ( terminal_res_save.size() > 0 ) {
		if ( terminal_res != terminal_res_save ) {
			TR << TR.Red << "Auto-computed TERMINAL_RES " << terminal_res <<      TR.Reset << std::endl;
			TR << TR.Red << "User-input    TERMINAL_RES " << terminal_res_save << TR.Reset << std::endl;
		}
	}
	runtime_assert( terminal_res_save.size() == 0 || terminal_res == terminal_res_save );

}

void
update_jump_res( utility::vector1< Size > & jump_res,
	utility::vector1< Size > const & extra_minimize_jump_res ) {

	runtime_assert( jump_res.size() % 2 == 0 );
	runtime_assert( extra_minimize_jump_res.size() % 2 == 0 );

	for ( Size n = 1; n <= extra_minimize_jump_res.size()/2; n++ ) {
		Size const & res1 = extra_minimize_jump_res[ 2*n - 1 ];
		Size const & res2 = extra_minimize_jump_res[ 2*n     ];
		bool matches_existing( false );
		for ( Size q = 1; q <= jump_res.size()/2; q++ ) {
			Size const & existing_res1 = jump_res[ 2*n - 1 ];
			Size const & existing_res2 = jump_res[ 2*n     ];
			if ( ( res1 == existing_res1 && res2 == existing_res2 ) ||
					( res1 == existing_res2 && res2 == existing_res1 ) ) {
				matches_existing = true;
				break;
			}
		}
		if ( matches_existing ) continue;
		jump_res.push_back( res1 );
		jump_res.push_back( res2 );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// An 'advanced' version of terminal_res -- places 'repulsion atoms' above or below base faces to
// prevent stacking during both sampling and minimization.
void
add_block_stack_variants( vector1< pose::Pose * > const & pose_pointers,
	vector1< vector1< Size > > const & pose_res_lists,
	vector1< Size > const & block_stack_above_res,
	vector1< Size > const & block_stack_below_res ) {
	using namespace core::chemical;
	for ( Size n = 1; n <= pose_pointers.size(); n++ ) {
		pose::Pose & pose = *( pose_pointers[ n ] ) ;
		vector1 < Size > const & pose_res_list = pose_res_lists[ n ] ;
		for ( Size i = 1; i <= pose_res_list.size(); i++ ) {
			if ( block_stack_above_res.has_value( pose_res_list[ i ] ) ) {
				runtime_assert( pose.residue_type( i ).is_RNA() );
				add_variant_type_to_pose_residue( pose, BLOCK_STACK_ABOVE, i );
			}
			if ( block_stack_below_res.has_value( pose_res_list[ i ] ) ) {
				runtime_assert( pose.residue_type( i ).is_RNA() );
				add_variant_type_to_pose_residue( pose, BLOCK_STACK_BELOW, i );
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
void
check_extra_minimize_res_are_input( utility::vector1< core::Size > const & extra_minimize_res,
	utility::vector1< core::Size > const & input_domain_map ) {
	for ( Size n = 1; n <= extra_minimize_res.size(); n++ ) {
		runtime_assert( input_domain_map[ extra_minimize_res[ n ] ] > 0 );
	}
}

/////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
figure_out_fixed_domain_map( utility::vector1< Size > const & input_domain_map,
	utility::vector1< Size > const & extra_minimize_res ) {
	vector1< Size > fixed_domain_map = input_domain_map;
	for ( Size n = 1; n <= extra_minimize_res.size(); n++ ) {
		fixed_domain_map[ extra_minimize_res[ n ] ] = 0;
	}
	return fixed_domain_map;
}

/////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
figure_out_dock_domain_map( utility::vector1< Size > & cutpoint_open_in_full_model,
	utility::vector1< utility::vector1< Size > > const & pose_res_lists,
	utility::vector1< Size > const & working_res,
	utility::vector1< Size > const & sample_res,
	Size const nres ) {

	// supplement cutpoint_open with breaks in working_res.
	vector1< Size > chains( nres, 0 ); // 0 outside working_res
	Size chain_number( 0 );
	for ( Size n = 1; n <= nres; n++ ) {
		if ( !working_res.has_value( n ) ) continue;
		bool new_chain( false );
		if ( ( n == 1 || !working_res.has_value( n - 1 )) && working_res.has_value( n ) ) new_chain = true;
		if ( n > 1 && cutpoint_open_in_full_model.has_value( n - 1 ) ) new_chain = true;
		if ( new_chain ) chain_number++;
		chains[ n ] = chain_number; // every working_res gets a chain number.
	}

	// which segments are connected by input poses?
	vector1< pair< Size, Size > > chain_connections;
	for ( vector1< Size > const & res_list : pose_res_lists ) {
		std::set< Size > chains_in_pose;
		for ( Size const seqpos : res_list ) {
			if ( sample_res.has_value( seqpos ) ) continue;
			if ( !working_res.has_value( seqpos ) ) continue;
			chains_in_pose.insert( chains[ seqpos ] );
		}
		for ( auto it1 = chains_in_pose.begin(), end = chains_in_pose.end(); it1 != end; ++it1 ) {
			for ( auto it2 = it1; it2 != end; ++it2 ) {
				if ( it1 != it2 ) chain_connections.push_back( make_pair( *it1, *it2 ) );
			}
		}
	}
	// find connected clusters
	utility::vector1< Size > connection_domains = get_connection_domains( chain_connections, max( chains ) );

	vector1< Size > dock_domain_map( nres, 0 );
	for ( Size n = 1; n <= nres; n++ ) {
		// to get from 1,2,... (connection_domains) to 0,1,... (convention in FullModelParameters)
		if ( !working_res.has_value( n ) ) continue;
		dock_domain_map[ n ] = connection_domains[ chains[ n ] ];
	}

	for ( Size k = 1; k < working_res.size(); k++ ) {
		Size const n = working_res[ k ];
		Size const n_next = working_res[ k+1 ];
		if (  dock_domain_map[ n ] != dock_domain_map[ n_next ] &&
				!cutpoint_open_in_full_model.has_value( n ) ) {
			TR << TR.Red << "There appears to be a dock boundary at " << n << " so adding to cutpoint_in_full_model." << TR.Reset << std::endl;
			cutpoint_open_in_full_model.push_back( n );
		}
	}

	return dock_domain_map;
}

/////////////////////////////////////////////////////////////////////////////
void
reorder_pose( pose::Pose & pose, utility::vector1< Size > & res_list ) {
	utility::vector1< Size > res_list_ordered = res_list;
	std::sort( res_list_ordered.begin(), res_list_ordered.end() );
	if ( res_list == res_list_ordered ) return;

	vector1< Size > slice_res;
	for ( Size n = 1; n <= res_list_ordered.size(); n++ ) slice_res.push_back( res_list.index( res_list_ordered[ n ] ) );
	pose::pdbslice( pose, slice_res );  // will do the rearrangement.
	res_list = res_list_ordered;
}

/////////////////////////////////////////////////////////////////////////////
// helper function to see if this is just an RNA modeling problem, as defined
// in fasta file.
bool
just_modeling_RNA( utility::vector1< std::string > const & fasta_files ) {
	runtime_assert( fasta_files.size() >= 1 );
	std::string sequence = core::sequence::read_fasta_file_return_str( fasta_files[1] );
	core::sequence::parse_out_non_standard_residues( sequence );
	return ( pose::just_modeling_RNA( sequence ) );
}


// AMW TODO: version that takes an OptionsCollection
// AMW: for now, create version that leaves
// alone if it has already been initialized.
///////////////////////////////////////////////////////////////
void
initialize_native_and_align_pose( PoseOP & native_pose,
	PoseOP & align_pose,
	core::chemical::ResidueTypeSetCOP rsd_set,
	PoseCOP start_pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( ( !native_pose || native_pose->size() == 0 ) && option[ in::file::native ].user() )  {
		// If native is on the command line, stuff all that in align_pose!
		align_pose = native_pose = core::import_pose::get_pdb_with_full_model_info( option[ in::file::native ](), rsd_set );
	} else if ( ( !align_pose || align_pose->size() == 0 ) && native_pose ) {
		// If native wasn't specified on the command line but it is already set up
		// nonetheless...
		align_pose = native_pose;
	}

	if ( option[ OptionKeys::stepwise::new_align_pdb ].user() ) {
		align_pose = core::import_pose::get_pdb_with_full_model_info(  option[ OptionKeys::stepwise::new_align_pdb ](), rsd_set );
	} else if ( ( !align_pose || align_pose->size() == 0 )  && option[ OptionKeys::stepwise::align_pdb ].user() ) {
		align_pose = core::import_pose::get_pdb_with_full_model_info(  option[ OptionKeys::stepwise::align_pdb ](), rsd_set );
	}

	if ( align_pose == nullptr && option[ in::file::s ].user() ) {
		align_pose = start_pose->clone();
	}
	if ( option[ OptionKeys::stepwise::virtualize_free_moieties_in_native ]() ) { // could generalize to proteins
		if ( native_pose != nullptr )  pose::rna::virtualize_free_rna_moieties( *native_pose );
		if ( align_pose  != nullptr ) pose::rna::virtualize_free_rna_moieties( *align_pose );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
process_input_file( std::string const & input_file,
	utility::vector1< pose::PoseOP > & pose_list,
	bool is_pdb /*= false*/,
	bool coarse_rna /* = false */)
{
	using namespace core::io::silent;

	core::chemical::ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	if ( is_pdb ) {

		pose::PoseOP pose_op( new pose::Pose );
		core::import_pose::pose_from_file( *pose_op, *rsd_set, input_file , core::import_pose::PDB_file);
		//   ensure_phosphate_nomenclature_matches_mini( *pose_op );
		core::pose::rna::figure_out_reasonable_rna_fold_tree( *pose_op );
		pose_list.push_back( pose_op );

	} else { //its a silent file.
		SilentFileOptions opts;
		SilentFileData silent_file_data( opts );
		silent_file_data.read_file( input_file );
		for ( core::io::silent::SilentFileData::iterator iter = silent_file_data.begin(),
				end = silent_file_data.end(); iter != end; ++iter ) {
			pose::PoseOP pose_op( new pose::Pose );
			iter->fill_pose( *pose_op );
			pose_list.push_back( pose_op );
		}

	}

	// further cleanup.
	for ( Size n = 1; n <= pose_list.size(); n++ ) {

		pose::PoseOP pose_op = pose_list[ n ];

		remove_cutpoints_closed( *pose_op );

		if ( coarse_rna && !pose_op->residue(1).is_coarse() ) {
			pose::Pose coarse_pose;
			make_coarse_pose( *pose_op, coarse_pose );
			*pose_op = coarse_pose;
		}

		core::pose::rna::virtualize_5prime_phosphates( *pose_op );
	}

	if ( pose_list.size() < 1 )  {
		utility_exit_with_message(  "No structure found in input file  " + input_file );
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////
bool
compare_RNA_char( char const char1, char const char2 ) {
	//Man this is silly, there must be a more elegant way to do this.
	if ( char1 == char2 ) return true;
	if ( char1 == 'n' || char2 == 'n' ) return true;
	if ( char1 == 'r' && (char2 == 'a' || char2 == 'g') ) return true;
	if ( char1 == 'y' && (char2 == 'c' || char2 == 'u') ) return true;
	if ( char2 == 'r' && (char1 == 'a' || char1 == 'g') ) return true;
	if ( char2 == 'y' && (char1 == 'c' || char1 == 'u') ) return true;
	return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
bool
compare_RNA_secstruct( char const char1, char const char2 ) {
	if ( char1 == char2 ) return true;
	if ( char1 == 'X' || char2 == 'X' ) return true;
	if ( char1 == 'L' && ( char2 == 'N' || char2 == 'P') ) return true;
	if ( char2 == 'L' && ( char1 == 'N' || char1 == 'P') ) return true;
	return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
std::string const
convert_based_on_match_type( std::string const & RNA_string, Size const type ){

	// AMW: this can probably be made more concise using a functional approach.

	std::string RNA_string_local = RNA_string;

	Size const size = RNA_string.length();

	static bool print_warning( false );

	//Obey orders to match exactly, match pyrimidine/purine, or match all.
	if ( type == MATCH_ALL ) {
		for ( Size i = 0; i < size; i++ )  RNA_string_local[ i ] = 'n';
	} else if ( type == MATCH_YR ) {
		for ( Size i = 0; i < size; i++ ) {
			if ( RNA_string[ i ] == 'g' || RNA_string[ i ] == 'a' ) {
				RNA_string_local[ i ] = 'r';
			} else {
				runtime_assert( RNA_string[ i ] == 'u' || RNA_string[ i ] == 'c' || RNA_string[ i ] == 't' );
				RNA_string_local[ i ] = 'y';
			}
		}
	} else {
		for ( Size i = 0; i < size; i++ )  {
			if ( RNA_string[ i ] == 't' ) {
				if ( !print_warning ) {
					TR.Warning << TR.Red << "Requesting an RNA fragment for t. Instead choosing fragment based on u!" << std::endl;
					print_warning = true;
				}
				RNA_string_local[ i ] = 'u';
			}
		}
	}

	return RNA_string_local;
}


/////////////////////////////////////////////////////////////////////////////////////
Size
get_anchor_rsd( pose::Pose const & pose )
{
	utility::vector1< Size > const rigid_body_jumps = pose::rna::get_rigid_body_jumps( pose );
	if ( rigid_body_jumps.size() == 0 ) return 0;

	Size const nres = pose.size(); // This better be a virtual residue -- checked in get_rigid_body_jumps() above.

	Size anchor_rsd = pose.fold_tree().downstream_jump_residue( rigid_body_jumps[1] );
	if ( anchor_rsd == nres ) anchor_rsd = pose.fold_tree().upstream_jump_residue( rigid_body_jumps[1] );
	return anchor_rsd;
}

////////////////////////////////////////////////////////////////////////////////////////
bool
involved_in_phosphate_torsion( std::string atomname )
{
	utility::vector1< std::string > const & atoms_involved = core::chemical::rna::atoms_involved_in_phosphate_torsion;

	for ( Size n = 1; n <= atoms_involved.size(); n++ ) {
		if (  atomname == atoms_involved[ n ] ) return true;
	}
	return false;
}



////////////////////////////////////////////////////////
void
remove_cutpoint_closed( pose::Pose & pose, Size const i ) {

	using namespace core::chemical;
	using namespace core::kinematics;

	remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, i );
	remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, i+1 );

	// using namespace protocols::forge::methods;
	//  FoldTree f( pose.fold_tree() );
	//  remove_cutpoint( i, f );
	//  pose.fold_tree( f );

	utility::vector1< Size > const & cutpoints = pose.fold_tree().cutpoints();

	Size const num_jump =  pose.fold_tree().num_jump();
	utility::vector1< Size > upstream_pos, downstream_pos;
	for ( Size n = 1; n <= num_jump; n++ ) {
		upstream_pos.push_back( pose.fold_tree().upstream_jump_residue( n ) );
		downstream_pos.push_back( pose.fold_tree().downstream_jump_residue( n ) );
	}

	ObjexxFCL::FArray1D< Size > cuts( num_jump-1 );
	Size count( 0 );
	for ( Size n = 1; n <= num_jump; n++ ) {
		if ( cutpoints[n] == i ) continue;
		count++;
		cuts( count ) = cutpoints[ n ];
	}

	Size const nres( pose.size() );

	// Just brute-force iterate through to find a jump we can remove.
	Size jump_to_remove( 0 );
	for ( Size k = 1; k <= num_jump; k++ ) {
		FArray1D< bool > partition_definition( nres, false );
		pose.fold_tree().partition_by_jump( k, partition_definition );
		if ( partition_definition( i ) != partition_definition( i+1 ) ) {
			jump_to_remove = k; break;
		}
	}

	bool success( false );

	ObjexxFCL::FArray2D< Size > jump_point( 2, num_jump-1 );

	count = 0;
	for ( Size n = 1; n <= num_jump; n++ ) {
		if ( n == jump_to_remove ) continue;
		count++;
		if ( upstream_pos[ n ] < downstream_pos[ n ] ) {
			jump_point( 1, count ) = upstream_pos[ n ];
			jump_point( 2, count ) = downstream_pos[ n ];
		} else {
			jump_point( 1, count ) = downstream_pos[ n ];
			jump_point( 2, count ) = upstream_pos[ n ];
		}
	}

	FoldTree f( nres );

	success = f.tree_from_jumps_and_cuts( nres, num_jump-1, jump_point, cuts, 1, false /*verbose*/ );

	if ( !success ) utility_exit_with_message( "FAIL to remove cutpoint "+string_of( i ) );

	pose.fold_tree( f );
}

////////////////////////////////////////////////////////
void
remove_cutpoints_closed( pose::Pose & pose ){
	// Make a list of each cutpoint_closed.
	for ( Size i = 1; i < pose.size(); i++ ) {
		if ( pose.fold_tree().is_cutpoint( i ) &&
				pose.residue_type( i   ).has_variant_type( chemical::CUTPOINT_LOWER ) &&
				pose.residue_type( i+1 ).has_variant_type( chemical::CUTPOINT_UPPER ) ) {
			remove_cutpoint_closed( pose, i ); // this will cycle through to find a jump that is removable.
		}
	}
}



/////////////////////////////////////////////////
void
make_extended_coarse_pose( pose::Pose & coarse_pose, std::string const & full_sequence ){

	using namespace core::chemical;
	using namespace core::id;

	ResidueTypeSetCOP rsd_set_coarse = ChemicalManager::get_instance()->residue_type_set( COARSE_RNA );

	make_pose_from_sequence( coarse_pose, full_sequence, *rsd_set_coarse );

	for ( Size n = 1; n <= coarse_pose.size(); n++ ) {
		coarse_pose.set_torsion( TorsionID( n, BB, 1 ), -150.0 );
		coarse_pose.set_torsion( TorsionID( n, BB, 2 ),  180.0 );
	}
}


///////////////////////////////////////////////////////////////////////////
void
make_coarse_pose( pose::Pose const & pose, pose::Pose & coarse_pose ){

	using namespace core::chemical;
	using namespace core::id;
	using namespace core::scoring::rna;

	ResidueTypeSetCOP rsd_set, rsd_set_coarse;
	RNA_CentroidInfo rna_centroid_info;

	rsd_set_coarse = ChemicalManager::get_instance()->residue_type_set( COARSE_RNA );
	make_extended_coarse_pose( coarse_pose, pose.sequence() );

	for ( Size n = 1; n <= pose.size(); n++ ) {
		coarse_pose.set_xyz(  NamedAtomID( " P  ", n ),  pose.xyz( NamedAtomID( " P  ", n )) );

		//coarse_pose.set_xyz(  NamedAtomID( " S  ", n ),  pose.xyz( NamedAtomID( " C4'", n )) );
		Vector sugar_centroid = pose::rna::get_sugar_centroid( pose.residue( n ) );
		coarse_pose.set_xyz(  NamedAtomID( " S  ", n ),  sugar_centroid );

		Vector base_centroid = rna_centroid_info.get_base_centroid( pose.residue( n ) );
		kinematics::Stub stub = rna_centroid_info.get_base_coordinate_system( pose.residue( n ), base_centroid );

		coarse_pose.set_xyz(  NamedAtomID( " CEN", n ),  base_centroid );
		coarse_pose.set_xyz(  NamedAtomID( " X  ", n ),  base_centroid + stub.M.col_x() );
		coarse_pose.set_xyz(  NamedAtomID( " Y  ", n ),  base_centroid + stub.M.col_y() );
	}

	core::kinematics::FoldTree f( pose.fold_tree() );
	for ( Size n = 1; n <= pose.num_jump(); n++ ) {
		f.set_jump_atoms( n, " Y  ", " Y  " );
	}
	coarse_pose.fold_tree( f );
}


} // namespace import_pose
} // namespace core
