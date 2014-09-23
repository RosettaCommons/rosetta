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
#include <core/pose/full_model_info/FullModelInfo.hh>

#include <core/util/metalloproteins_util.hh>

#include <core/pack/pack_missing_sidechains.hh>
#include <core/pack/optimizeH.hh>

#include <core/io/pdb/pdb_dynamic_reader.hh>
#include <core/io/pdb/pdb_dynamic_reader_options.hh>
#include <core/io/pdb/file_data.hh>
#include <core/io/pdb/file_data_options.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <boost/foreach.hpp>


namespace core {
namespace import_pose {

using core::Size;
using core::SSize;

using basic::T;
using basic::Error;
using basic::Warning;

using namespace core::io::pdb;
using namespace ObjexxFCL;

static thread_local basic::Tracer TR( "core.import_pose.import_pose" );

using utility::vector1;


void
read_additional_pdb_data(
	std::string const & s,
	pose::Pose & pose,
	io::pdb::FileData const &,
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

	//Added by Daniel-Adriano Silva, used to read PDBinfo-LABEL
	//TR.Debug << "Setting PDBinfo-labels from PDB file." << std::endl;
	for ( i=1; i<= lines.size(); ++i ) {
		std::string const & line( lines[i] );
		if( line.size() > 21 && line.substr(0,21) == "REMARK PDBinfo-LABEL:" ){
			//Parse and split string
			utility::vector1 < std::string > remark_values;
			utility::vector1 < std::string > tmp_remark_values = utility::string_split(line, ' ');
			//Copy non-empty (i.e. !' ') elements to remark_values
			if ( tmp_remark_values .size() > 3){
				for ( Size j=3; j<= tmp_remark_values.size(); ++j ) {
					if(tmp_remark_values[j] != ""){
						remark_values.push_back(tmp_remark_values[j]);
					}
				}
			}
			//Check that we have at least two elements left ([1]=index, [2-n]=PDBinfo-labels)
			if (remark_values.size() > 1){
				core::Size tmp_ndx=atoi(remark_values[1].c_str());
				if ( tmp_ndx <= pose.total_residue() ) {
					for ( Size j=2; j<= remark_values.size(); ++j ) {
						pose.pdb_info()->add_reslabel(tmp_ndx,remark_values[j]);
					}
				}else{
					TR.Fatal << "pose_io:: PDBinfo-LABEL io failure: " << line << ' ' << pose.total_residue()  << std::endl;
				}
			}else{
				TR.Fatal << "pose_io:: PDBinfo-LABEL io failure: " << line << ' ' << pose.total_residue() << std::endl;
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
	pose::PoseOP pose( new pose::Pose() );
	pose_from_pdb( *pose, filename, read_fold_tree);
	return pose;
}


pose::PoseOP pose_from_pdb(chemical::ResidueTypeSet const & residue_set, std::string const & filename, 	bool read_fold_tree)
{
	pose::PoseOP pose( new pose::Pose() );
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
	ImportPoseOptions options;
	if ( residue_set.name() == core::chemical::FA_RNA ) options.set_rna( true );
	pose_from_pdb(pose, residue_set, filenames_string, options, read_fold_tree);
}

void
pose_from_pdb(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filenames_string,
	ImportPoseOptions const & options,
	bool read_fold_tree
)
{
	utility::vector1<std::string> filenames = utility::split(filenames_string);

	std::string res;

	BOOST_FOREACH(std::string filename, filenames){
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
	conformation::ConformationOP conformation_op( new conformation::Conformation() );
	if ( !pose.conformation().same_type_as_me( *conformation_op, true ) ) {
		pose.set_new_conformation( conformation_op );
	}

	io::pdb::FileData fd = core::import_pose::PDB_DReader::createFileData(res, options);
	if ( fd.filename == "" ) {
		fd.filename = utility::join(filenames, "_");
	}
	build_pose(fd, pose, residue_set, options);

	// set secondary structure for centroid PDBs
	if ( residue_set.name() == core::chemical::CENTROID ) {
		core::pose::set_ss_from_phipsi( pose );
	}

	// check for foldtree info
	read_additional_pdb_data( res, pose, options, read_fold_tree );
}


void
pose_from_pdb(
	core::pose::Pose & pose,
	std::string const & filename,
	bool read_fold_tree
) {
	ImportPoseOptions options;
	pose_from_pdb( pose, filename, options, read_fold_tree );
}

void
pose_from_pdb(
	pose::Pose & pose,
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree
) {
	using namespace chemical;

	ResidueTypeSetCOP residue_set( options.centroid() ?
		ChemicalManager::get_instance()->residue_type_set( CENTROID ) :
		ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )
	);

	if ( options.rna() ) residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	core::import_pose::pose_from_pdb( pose, *residue_set, filename, options, read_fold_tree );
}

utility::vector1< core::pose::PoseOP >
poseOPs_from_pdbs(
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree
) {
	ImportPoseOptions options;
	return poseOPs_from_pdbs( filenames, options, read_fold_tree );
}

utility::vector1< core::pose::PoseOP >
poseOPs_from_pdbs(
	utility::vector1< std::string > const & filenames,
	ImportPoseOptions const & options,
	bool read_fold_tree
) {
	using namespace chemical;

	ResidueTypeSetCOP residue_set( options.centroid() ?
		ChemicalManager::get_instance()->residue_type_set( CENTROID ) :
		ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )
	);

	return core::import_pose::poseOPs_from_pdbs( *residue_set, filenames, options, read_fold_tree );
}

/// @details Only returns full-atom poses
utility::vector1< core::pose::Pose >
poses_from_pdbs(
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree
) {
	using namespace chemical;
	ResidueTypeSetCOP residue_set
		( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

	return core::import_pose::poses_from_pdbs( *residue_set, filenames, read_fold_tree );
}

utility::vector1< core::pose::Pose >
poses_from_pdbs(
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree
) {
	using namespace chemical;

	using std::string;
	using utility::vector1;
	using core::pose::Pose;

	ImportPoseOptions options;
	if ( residue_set.name() == core::chemical::FA_RNA ) options.set_rna( true );

	vector1< Pose > poses;
	typedef vector1< string >::const_iterator vec_it;
	for ( vec_it it = filenames.begin(), end = filenames.end(); it != end; ++it ) {
		Pose pose;
		core::import_pose::pose_from_pdb( pose, residue_set, *it, options, read_fold_tree );
		poses.push_back( pose );
	}

	return poses;
}

void
read_all_poses(
	const utility::vector1<std::string>& filenames,
	utility::vector1<core::pose::Pose>* poses
)
{
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

utility::vector1< core::pose::PoseOP >
poseOPs_from_pdbs(
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< std::string > const & filenames,
	ImportPoseOptions const & options,
	bool read_fold_tree
) {
	using namespace chemical;

	using std::string;
	using utility::vector1;
	using core::pose::Pose;

	vector1< pose::PoseOP > poses;
	typedef vector1< string >::const_iterator vec_it;
	for ( vec_it it = filenames.begin(), end = filenames.end(); it != end; ++it ) {
		pose::PoseOP pose( new pose::Pose );
		core::import_pose::pose_from_pdb( *pose, residue_set, *it, options, read_fold_tree );
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
	ResidueTypeSetCOP residue_set(
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
	ImportPoseOptions options;
	if ( residue_set.name() == core::chemical::FA_RNA ) options.set_rna( true );
	pose_from_pdb( poses, residue_set, filename, options, read_fold_tree );
}

void
pose_from_pdb(
	utility::vector1< pose::Pose > & poses,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	ImportPoseOptions const & options,
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

			io::pdb::FileData fd = core::import_pose::PDB_DReader::createFileData( sub_lines, options );
			fd.filename = filename;
			build_pose( fd, pose, residue_set, options );

			// check for foldtree info
			core::import_pose::read_additional_pdb_data( sub_lines, pose, options, read_fold_tree);
			poses.push_back( pose );
		}
	} else {
		FileData fd = core::import_pose::PDB_DReader::createFileData( all_lines, options );
		if ( fd.filename == "" ) {
			fd.filename = filename;
		}
		core::import_pose::build_pose( fd, pose, residue_set, options );
		// check for foldtree info
		core::import_pose::read_additional_pdb_data( all_lines, pose, options, read_fold_tree);
		poses.push_back( pose );
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
	io::pdb::FileData fd = import_pose::PDB_DReader::createFileData( pdbcontents, options );
	fd.filename = filename;
	chemical::ResidueTypeSetCOP residue_set
		( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );
	core::import_pose::build_pose( fd, pose, *residue_set, options );

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
	if ( residue_set.name() == core::chemical::FA_RNA ) options.set_rna( true );
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
	io::pdb::FileData fd = import_pose::PDB_DReader::createFileData( pdbcontents, options );
	fd.filename = filename;
	core::import_pose::build_pose( fd, pose, residue_set, options);
}

void pose_from_pdb_stream(
	pose::Pose & pose,
	std::istream & pdb_stream,
	std::string const & filename,
	ImportPoseOptions const & options
){
	std::string pdb_file_contents;
	utility::slurp( pdb_stream, pdb_file_contents );

	chemical::ResidueTypeSetCOP residue_set( chemical::ChemicalManager::get_instance()->residue_type_set( options.residue_type_set()) );
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
		( ChemicalManager::get_instance()->residue_type_set( CENTROID ) );

	core::import_pose::pose_from_pdb( pose, *residue_set, filename, read_fold_tree);
}

void build_pose(
	io::pdb::FileData & fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set
)
{
	ImportPoseOptions options; // read from the command line
	if ( residue_set.name() == core::chemical::FA_RNA ) options.set_rna( true );
	build_pose( fd, pose, residue_set, options);
}

/// @brief Build Rosetta 3 Pose object from FileData.
void build_pose(
	io::pdb::FileData & fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options
)
{
	TR.Debug << "build_pose..." << std::endl;
	build_pose_as_is( fd, pose, residue_set, options);
	TR.Debug << "build_pose... Ok." << std::endl;
}

void build_pose_as_is2(
		io::pdb::FileData & fd,
		pose::Pose & pose,
		chemical::ResidueTypeSet const & residue_set,
		id::AtomID_Mask & missing,
		ImportPoseOptions const & options );

// "super-simple" (C) by Phil
//
/// @brief Try to Build pose object from pdb 'as-is'. PDB file must be _really_ clean.
void build_pose_as_is(
	io::pdb::FileData & fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options
)
{
	id::AtomID_Mask missing( false );

	io::pdb::build_pose_as_is1( fd, pose, residue_set, missing, options );
	build_pose_as_is2( fd, pose, residue_set, missing, options );
}

void build_pose_as_is2(
	io::pdb::FileData & /*fd*/,
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
	if ( residue_set.name() == FA_STANDARD ) {
		//if pack_missing_density specified, repack residues w/ missing density
		if( options.pack_missing_sidechains() ) {
			pack::pack_missing_sidechains( pose, missing );
		}
		// optimize H if using a fullatom residue type set, and no_optH is not specified
		if( !options.no_optH() ) {
			pack::optimize_H_and_notify( pose, missing );
		}
	}


	// If pose contains carbohydrate residues, assure that their virtual atoms have the correct coordinates.
	if (pose.conformation().contains_carbohydrate_residues()) {
		for (uint i = 1, n_residues = pose.total_residue(); i <= n_residues; ++i) {
			ResidueType const & res_type = pose.residue_type(i);
			if (res_type.is_carbohydrate()) {
				pose::carbohydrates::align_virtual_atoms_in_carbohydrate_residue(pose, i);
			}
		}
	}

	// If the user has set appropriate flags, check whether the pose contains metal ions,
	// and automatically set up covalent bonds and constraints to them.
	if( options.set_up_metal_bonds() ) {
		core::util::auto_setup_all_metal_bonds(pose, options.metal_bond_LJ_multiplier(), true);
		if( options.set_up_metal_constraints() ) {
			core::util::auto_setup_all_metal_constraints( pose, options.metal_bond_dist_constraint_multiplier(), options.metal_bond_angle_constraint_multiplier() );
		}
	}

	pose.pdb_info( pdb_info );
	pose::full_model_info::make_sure_full_model_info_is_setup( pose );
}

/// @details
/// All ligand residues and polymer branches have been appended by a jump.  This method creates a new fold tree without
/// jumps through ligands, using CHEMICAL edges instead.
void
set_reasonable_fold_tree(pose::Pose & pose)
{
	// An empty pose doesn't have jumps through ligands
	// (Will encounter a SegFault otherwise)
	if( pose.total_residue() == 0 ) {
		return;
	}

	using namespace core::chemical;
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
		// Jump to a ligand residue or polymer branch?
		if( e.is_jump() &&
				(pose.residue_type(e.stop()).has_variant_type(BRANCH_LOWER_TERMINUS_VARIANT) ||
				!pose.residue_type(e.stop()).is_polymer() ) ) {
			Size const ii = e.stop(); // the residue at the end of the jump
			conformation::Residue const & ii_res = pose.residue(ii);
			// Now we'll prepare a chemical edge by first finding the connecting atoms
			// between the two residues
			bool found_connection_residue_in_fold_tree( false );
			for ( Size jj = 1; jj <= pose.residue_type( ii ).n_residue_connections(); ++jj ) {
				if(ii_res.connection_incomplete(jj)) continue;  // allow incomplete connections for design
				Size jj_res_ID= ii_res.connect_map( jj ).resid();
				if ( jj_res_ID < ii){
					core::conformation::Residue const &  jj_res=pose.residue(jj_res_ID);
					// Ensure that the connection is either a polymer branching or a ligand of the same chain.
					if ((jj_res.has_variant_type(BRANCH_POINT_VARIANT) && ii_res.has_variant_type(BRANCH_LOWER_TERMINUS_VARIANT)) ||
							(jj_res.chain() == ii_res.chain())) {
						int ii_connect_ID = ii_res.connect_atom(jj_res);
						int jj_connect_ID = jj_res.connect_atom(ii_res);

						std::string ii_connector = ii_res.atom_name(ii_connect_ID);
						std::string jj_connector = jj_res.atom_name(jj_connect_ID);

						newft.add_edge(jj_res_ID, ii, jj_connector, ii_connector);

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
