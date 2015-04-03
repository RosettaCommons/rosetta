// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file motif_utils.cc
/// @brief Motif helper/conversion/io functions
/// @author havranek, sthyme (sthyme@gmail.com)

// Unit Headers
#include <protocols/motifs/motif_utils.hh>

// Package Headers
#include <protocols/motifs/BuildPosition.hh>
#include <protocols/motifs/Motif.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <protocols/motifs/SingleMotif.hh>

// Project Headers (protocols)
#include <protocols/dna/DnaDesignDef.hh>
#include <protocols/dna/util.hh>
#include <protocols/loops/Loops.hh>

// Project Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.functions.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/graph/Graph.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/dna/setup.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/tools/make_map.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>

// C++ Headers

// Option Key Includes
#include <basic/options/option.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>

namespace protocols {
namespace motifs {

static thread_local basic::Tracer mu_tr( "protocols.motifs.motif_utils", basic::t_info );

//this function uses the constructor that accepts two residues, rather than input atoms, because if the constructors that allow you to setup the atoms explicitly are necessary you are likely to have the information setup in a motif_file (or you can simply make another matching function that allows input atoms as well as a pdb, rather than a istream)
SingleMotifOP
single_motif_from_filename(
	utility::file::FileName const & motif_filename
)
{
	core::pose::PoseOP pose( new core::pose::Pose );
	core::import_pose::pose_from_pdb( *pose, motif_filename );
	core::conformation::Residue res1( pose->residue(1) ); //the motifs in my motif pdbs are setup so that residue #1 is the amino acid and #2 is the base
	core::conformation::Residue res2( pose->residue(2) ); //this function at some point include some checks to see what residue is base or amino acid or both bases, both amino acids (or something else, like a ligand, etc)

	SingleMotifOP retval( new SingleMotif( res1, res2 ) );
	std::string this_remark( motif_filename.base() ); //I am using the remark to keep track of the pdb name for output purposes (I don't know what Jim is using the remark for, since he made it)
	std::string this_path( motif_filename.path() ); //I am using the remark to keep track of the pdb name for output purposes (I don't know what Jim is using the remark for, since he made it)
	retval->store_remark( this_remark );
	retval->store_path( this_path );

	return retval;
}

SingleMotifOP
single_motif_from_stream(
	utility::io::izstream & motif_info
)
{
	std::string res1;
	std::string r1atom1, r1atom2, r1atom3;
	std::string res2;
	std::string r2atom1, r2atom2, r2atom3;
	core::kinematics::Jump motif_jump;

	core::Size const max_line_len(1024);
	char first_line[max_line_len];
	motif_info.getline( first_line, sizeof(first_line) );
	std::istringstream line_in( first_line );

	line_in >> res1;
	line_in >> r1atom1;
	line_in >> r1atom2;
	line_in >> r1atom3;
	line_in >> res2;
	line_in >> r2atom1;
	line_in >> r2atom2;
	line_in >> r2atom3;

	char rest_of_line[max_line_len];
	line_in.getline( rest_of_line, sizeof(rest_of_line) );
	std::string rest( rest_of_line );
	std::string::size_type index( rest.find("REMARK ") );
	std::string this_remark;
	bool has_remark( false );
	if( index != std::string::npos ) {
		this_remark = rest.substr( index + 7 );
		has_remark = true;
	} else {
	}

	// Scan to end of line to skip comments
	// motif_info.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	// motif_info.get( remark_in, max_remark );

	motif_info >> motif_jump;

	SingleMotifOP retval( new SingleMotif( res1, r1atom1, r1atom2, r1atom3, res2, r2atom1, r2atom2, r2atom3, motif_jump ) );
	if( has_remark ) {
		retval->store_remark( this_remark );
	}

	return retval;
}

SingleMotifOP
single_motif_from_stream(
	std::istream & motif_info
)
{
	std::string res1;
	std::string r1atom1, r1atom2, r1atom3;
	std::string res2;
	std::string r2atom1, r2atom2, r2atom3;
	core::kinematics::Jump motif_jump;

	core::Size const max_line_len(1024);
	char first_line[max_line_len];
	motif_info.getline( first_line, sizeof(first_line) );
	std::istringstream line_in( first_line );

	line_in >> res1;
	line_in >> r1atom1;
	line_in >> r1atom2;
	line_in >> r1atom3;
	line_in >> res2;
	line_in >> r2atom1;
	line_in >> r2atom2;
	line_in >> r2atom3;

	char rest_of_line[max_line_len];
	line_in.getline( rest_of_line, sizeof(rest_of_line) );
	//std::cout << "Rest of line:" << rest_of_line << std::endl;
	std::string rest( rest_of_line );
	std::string::size_type index( rest.find("REMARK ") );
	std::string this_remark;
	bool has_remark( false );
	if( index != std::string::npos ) {
		this_remark = rest.substr( index + 7 );
		has_remark = true;
		//	std::cout << "Remark is " << this_remark << std::endl;
		} else {
		//	std::cout << "No remark in motif" << std::endl;
	}

	// Scan to end of line to skip comments
	// motif_info.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	// motif_info.get( remark_in, max_remark );

	motif_info >> motif_jump;

	SingleMotifOP retval( new SingleMotif( res1, r1atom1, r1atom2, r1atom3, res2, r2atom1, r2atom2, r2atom3, motif_jump ) );
	if( has_remark ) {
		retval->store_remark( this_remark );
	}

	return retval;
}

SingleMotifOP
single_ligand_motif_from_stream(
	std::istream & motif_info
)
{
	std::string res1;
	std::string r1atom1, r1atom2, r1atom3;
	std::string res2;
	std::string r2atom1, r2atom2, r2atom3;
	core::kinematics::Jump motif_jump;

	core::Size const max_line_len(1024);
	char first_line[max_line_len];
	motif_info.getline( first_line, sizeof(first_line) );
	std::istringstream line_in( first_line );

	line_in >> res1;
	line_in >> r1atom1;
	line_in >> r1atom2;
	line_in >> r1atom3;
	line_in >> res2;
	line_in >> r2atom1;
	line_in >> r2atom2;
	line_in >> r2atom3;

	char rest_of_line[max_line_len];
	line_in.getline( rest_of_line, sizeof(rest_of_line) );
	std::string rest( rest_of_line );
	std::string::size_type index( rest.find("REMARK ") );
	std::string this_remark;
	bool has_remark( false );
	if( index != std::string::npos ) {
		this_remark = rest.substr( index + 7 );
		has_remark = true;
		} else {
	}
	motif_info >> motif_jump;
	SingleMotifOP retval( new SingleMotif( res1, r1atom1, r1atom2, r1atom3, r2atom1, r2atom2, r2atom3, motif_jump ) );
	if( has_remark ) {
		retval->store_remark( this_remark );
	}

	return retval;
}

core::Real
parallel_base_test(
	core::conformation::Residue const & pose_dna,
	core::conformation::Residue const & motif_dna
)
{
	typedef numeric::xyzVector< core::Real > xyzVec;
	xyzVec p( core::scoring::dna::get_z_axis( motif_dna, core::scoring::dna::get_y_axis( motif_dna, 1 ) ) );
	// The second argument is the chain and if the value is changed from 1 to 2 the resulting value becomes negative
	xyzVec q( core::scoring::dna::get_z_axis( pose_dna, core::scoring::dna::get_y_axis( pose_dna, 1 ) ) );
	core::Real dotproduct( p.dot(q));
	// The dot product should equal 1.0 if the motif base is parallel to the dna that it will be connected with - a reasonable cutoff is 0.95
	return dotproduct;
}

core::Real
backbone_stub_match(
	core::conformation::Residue const & r1,
	core::conformation::Residue const & r2
)
{
	core::Real retval( 0.0 );

	retval += r1.xyz( "N" ).distance_squared( r2.xyz( "N" ) );
	retval += r1.xyz( "CA" ).distance_squared( r2.xyz( "CA" ) );
	retval += r1.xyz( "C" ).distance_squared( r2.xyz( "C" ) );

	// For the fourth atom, if either residue is glycine, you need to use HA, else use CB
	if( r1.type().aa() == core::chemical::aa_gly || r2.type().aa() == core::chemical::aa_gly ) {
		core::Size index1( r1.type().aa() == core::chemical::aa_gly ? r1.atom_index( "1HA" ) : r1.atom_index( "HA" ) );
		core::Size index2( r2.type().aa() == core::chemical::aa_gly ? r2.atom_index( "1HA" ) : r2.atom_index( "HA" ) );
		retval += r1.xyz( index1 ).distance_squared( r2.xyz( index2 ) );
	} else {
		retval += r1.xyz( "CB" ).distance_squared( r2.xyz( "CB" ) );
	}

	retval = std::sqrt( 0.25 * retval );

	return retval;
}

void
add_motif_bb_constraints(
	core::scoring::constraints::ConstraintSetOP cst_set,
	core::pose::Pose & pose,
  core::Size this_pos,
  core::conformation::Residue const & inv_rotamer
)
{
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	FuncOP fx1( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	cst_set->add_constraint( ConstraintCOP( ConstraintOP( new CoordinateConstraint( core::id::AtomID( pose.residue( this_pos ).atom_index( "CA" ), this_pos ),
																											core::id::AtomID( pose.residue( 1 ).atom_index( "CA" ), 1 ),
																											inv_rotamer.xyz( "CA" ),
																											fx1 ) ) ) );

	FuncOP fx2( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	cst_set->add_constraint( ConstraintCOP( ConstraintOP( new CoordinateConstraint( core::id::AtomID( pose.residue( this_pos ).atom_index( "C" ), this_pos ),
																											core::id::AtomID( pose.residue( 1 ).atom_index( "CA" ), 1 ),
																											inv_rotamer.xyz( "C" ),
																											fx2 ) ) ) );

	FuncOP fx3( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	cst_set->add_constraint( ConstraintCOP( ConstraintOP( new CoordinateConstraint( core::id::AtomID( pose.residue( this_pos ).atom_index( "N" ), this_pos ),
																											core::id::AtomID( pose.residue( 1 ).atom_index( "CA" ), 1 ),
																											inv_rotamer.xyz( "N" ),
																											fx3 ) ) ) );

	// For the fourth atom, if either residue is glycine, you need to use HA, else use CB
	if( pose.residue( this_pos ).type().aa() == core::chemical::aa_gly || inv_rotamer.type().aa() == core::chemical::aa_gly ) {
		core::Size index1( pose.residue( this_pos ).type().aa() == core::chemical::aa_gly ?
					pose.residue( this_pos ).atom_index( "1HA" ) : pose.residue( this_pos ).atom_index( "HA" ) );
		core::Size index2( inv_rotamer.type().aa() == core::chemical::aa_gly ?
					inv_rotamer.atom_index( "1HA" ) : inv_rotamer.atom_index( "HA" ) );

		FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
		cst_set->add_constraint( ConstraintCOP( ConstraintOP( new CoordinateConstraint( core::id::AtomID( index1, this_pos ),
																											core::id::AtomID( pose.residue( 1 ).atom_index( "CA" ), 1 ),
																											inv_rotamer.xyz( index2 ),
																											fx ) ) ) );
	} else {
		FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
		cst_set->add_constraint( ConstraintCOP( ConstraintOP( new CoordinateConstraint( core::id::AtomID( pose.residue( this_pos ).atom_index( "CB" ), this_pos ),
																											core::id::AtomID( pose.residue( 1 ).atom_index( "CA" ), 1 ),
																											inv_rotamer.xyz( "CB" ),
																											fx ) ) ) );
	}

	return;
}

void
add_motif_sc_constraints(
	core::scoring::constraints::ConstraintSetOP cst_set,
	core::pose::Pose & pose,
  core::Size this_pos,
  core::conformation::Residue const & inv_rotamer,
	MotifCOP this_motif,
	bool const is_it_forward
)
{
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	core::Size const index1( pose.residue( this_pos ).atom_index(
		( is_it_forward ? this_motif->res2_atom1_name() : this_motif->res1_atom1_name() )
	) );

	core::Size const index2( pose.residue( this_pos ).atom_index(
		( is_it_forward ? this_motif->res2_atom2_name() : this_motif->res1_atom2_name() )
	) );

	core::Size const index3( pose.residue( this_pos ).atom_index(
		( is_it_forward ? this_motif->res2_atom3_name() : this_motif->res1_atom3_name() )
	) );

	// This section is necessary in case pose.residue(1) is not a protein residue, because a DNA residue would not have a CA atom
	core::Size first_protein_resi(1);
	core::Size nres( pose.total_residue() );
	for ( core::Size i(1); i <= nres; ++i ) {
		if ( pose.residue_type(i).is_protein() ) {
			first_protein_resi = i;
			break;
		}
	}

	FuncOP fx1( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	cst_set->add_constraint( ConstraintCOP( ConstraintOP( new CoordinateConstraint( core::id::AtomID( index1, this_pos ),
																											core::id::AtomID( pose.residue( first_protein_resi ).atom_index( "CA" ), 1 ),
																											inv_rotamer.xyz( index1 ),
																											fx1 ) ) ) );

	FuncOP fx2( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	cst_set->add_constraint( ConstraintCOP( ConstraintOP( new CoordinateConstraint( core::id::AtomID( index2, this_pos ),
																											core::id::AtomID( pose.residue( first_protein_resi ).atom_index( "CA" ), 1 ),
																											inv_rotamer.xyz( index2 ),
																											fx2 ) ) ) );

	FuncOP fx3( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	cst_set->add_constraint( ConstraintCOP( ConstraintOP( new CoordinateConstraint( core::id::AtomID( index3, this_pos ),
																											core::id::AtomID( pose.residue( first_protein_resi ).atom_index( "CA" ), 1 ),
																											inv_rotamer.xyz( index3 ),
																											fx3 ) ) ) );

	return;
}

void mutate_position_vector_for_search(
	core::pose::Pose & pose,
	utility::vector1< core::Size > & trim_positions
)
{
	// Get a score function
	core::scoring::ScoreFunctionOP score_fxn( core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ) );

	// Set up the packer task
	core::pack::task::PackerTaskOP alanize_task( core::pack::task::TaskFactory::create_packer_task( pose ) );

	utility::vector1< bool > allow_vector( core::chemical::num_canonical_aas, false );
	allow_vector[ core::chemical::aa_ala ] = true;

	for( core::Size ires = 1, end_i = pose.total_residue() ; ires <= end_i ; ++ires ) {

		if( find( trim_positions.begin(), trim_positions.end(), ires ) != trim_positions.end() ) {
			if( pose.residue( ires ).aa() != core::chemical::aa_pro &&
					pose.residue( ires ).aa() != core::chemical::aa_gly ) {
				// Add it to the packer task
				alanize_task->nonconst_residue_task( ires ).restrict_absent_canonical_aas( allow_vector );
			} else {
				// Let prolines and glycines stay
				alanize_task->temporarily_set_pack_residue( ires, false );
			}
		} else {
			// Leave non-loop residues alone
			alanize_task->temporarily_set_pack_residue( ires, false );
		}

	}

	// Do the actual packing
	core::pack::pack_rotamers( pose, *score_fxn, alanize_task );
}

void mutate_loops_for_search(
	core::pose::Pose & pose,
	protocols::loops::Loops & flex_regions
)
{
	// Get a score function
	core::scoring::ScoreFunctionOP score_fxn( core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ) );

	// Set up the packer task
	core::pack::task::PackerTaskOP alanize_task( core::pack::task::TaskFactory::create_packer_task( pose ) );

	utility::vector1< bool > allow_vector( core::chemical::num_canonical_aas, false );
	allow_vector[ core::chemical::aa_ala ] = true;

	for( core::Size ires = 1, end_res = pose.total_residue() ; ires <= end_res ; ++ires ) {
		if( flex_regions.is_loop_residue( ires ) ) {
			if( pose.residue( ires ).aa() != core::chemical::aa_pro &&
					pose.residue( ires ).aa() != core::chemical::aa_gly ) {
				// Add it to the packer task
				alanize_task->nonconst_residue_task( ires ).restrict_absent_canonical_aas( allow_vector );
			} else {
				// Let prolines and glycines stay
				alanize_task->temporarily_set_pack_residue( ires, false );
			}
		} else {
			// Leave non-loop residues alone
			alanize_task->temporarily_set_pack_residue( ires, false );
		}
	}

	// Do the actual packing
	core::pack::pack_rotamers( pose, *score_fxn, alanize_task );
}

// The following functions are all related to user input and setup for the MotifSearch
MotifLibrary const
get_MotifLibrary_user()
{
	using namespace basic::options;
	// Make MotifLibrary from list of motif pdbs
	if ( option[ OptionKeys::motifs::list_motifs ].user() ) {
		utility::vector1< utility::file::FileName > listnames( option[ OptionKeys::motifs::list_motifs ]().vector() );
		utility::vector1< utility::file::FileName > motifnames( get_filenames( listnames ) );
		MotifLibrary motifs( motifnames );
		return motifs;
		// Make MotifLibrary from a file containing motifs
	} else if ( option[ OptionKeys::motifs::motif_filename ].user() ) {
		std::string motif_filename( option[ OptionKeys::motifs::motif_filename ]() );
		MotifLibrary motifs;
		motifs.add_from_file( motif_filename );
		return motifs;
	} else {
		// (not anymore) Program terminates if no MotifLibrary is given by the command line options
		// utility_exit_with_message( "ERROR! User must provide input motifs as either a list of PDBs or as a motif file" );
		mu_tr << "User did not provide input motifs via cmd line, but could be coming in as BuildPosition specific data" << std::endl;
		MotifLibrary motifs;
		return motifs;
	}
}

MotifLibrary const
get_LigandMotifLibrary_user()
{
	using namespace basic::options;
	// Make MotifLibrary from list of motif pdbs
	if ( option[ OptionKeys::motifs::list_motifs ].user() ) {
		utility::vector1< utility::file::FileName > listnames( option[ OptionKeys::motifs::list_motifs ]().vector() );
		mu_tr << "In get_LigandMotifLibrary, it's working" << std::endl;
		utility::vector1< utility::file::FileName > motifnames( get_filenames( listnames ) );
		mu_tr << "" << std::endl;
		MotifLibrary motifs( motifnames );
		mu_tr << "" << std::endl;
		return motifs;
		// Make MotifLibrary from a file containing motifs
	} else if ( option[ OptionKeys::motifs::motif_filename ].user() ) {
		std::string motif_filename( option[ OptionKeys::motifs::motif_filename ]() );
		mu_tr << "Got filename" << std::endl;
		MotifLibrary motifs;
		mu_tr << "Made motiflibrary" << std::endl;
		motifs.add_ligand_from_file( motif_filename ); //This is main change to continue on ligand path
		mu_tr << "added motifs from file" << std::endl;
		return motifs;
		//mu_tr << "returned motifs" << std::endl;  // This code will never be reached. ~Labonte
	} else {
		// (not anymore) Program terminates if no MotifLibrary is given by the command line options
		// utility_exit_with_message( "ERROR! User must provide input motifs as either a list of PDBs or as a motif file" );
		mu_tr << "User did not provide input motifs via cmd line, but could be coming in as BuildPosition specific data" << std::endl;
		MotifLibrary motifs;
		return motifs;
	}
}

utility::vector1< core::conformation::ResidueOP > const
get_targetconformers_user()
{
	using namespace basic::options;
	// Load a list of conformer PDBs, DNA in my case
	// At some point there will be input options other than PDB format
	// Currently this fxn will add all the conformers to one vector, even
	// if they are of multiple ResidueTypes, ie all 4 DNA bases
	utility::vector1< core::conformation::ResidueOP > conformerOPs;
	if ( option[ OptionKeys::motifs::list_dnaconformers ].user() ) {
		utility::vector1< utility::file::FileName > listnames( option[ OptionKeys::motifs::list_dnaconformers ]().vector() );
		utility::vector1< utility::file::FileName > conformernames( get_filenames( listnames ) );
		for ( utility::vector1< utility::file::FileName >::const_iterator filename( conformernames.begin() );
				filename != conformernames.end(); ++filename ) {
			if ( !utility::file::file_exists( *filename ) ) {
				continue;
			}
			core::pose::PoseOP pose( new core::pose::Pose );
			core::import_pose::pose_from_pdb( *pose, *filename );
			if ( pose->total_residue() > 1 ) {
				std::cerr << "WARNING!!! Conformer PDB contains more than one residue, loading all residues in PDB as conformers." << std::endl;
			}
			for ( core::Size i(1); i <= pose->total_residue(); ++i ) {
				conformerOPs.push_back( pose->residue(i).clone() );
			}
		}
	} else {
	}
	return conformerOPs;
}

std::map< std::string, utility::vector1< core::conformation::ResidueOP > > const
setup_conformer_map(
	utility::vector1< core::conformation::ResidueOP > const & conformerOPs
)
{
	using namespace core::conformation;
	std::map< std::string, ResidueOPs > conformer_map;
	for	( ResidueOPs::const_iterator itr = conformerOPs.begin(), end_itr = conformerOPs.end();
			itr != end_itr; ++itr ) {
		std::string name( (*itr)->name3() );
		conformer_map[name].push_back( *itr );
	}
	return conformer_map;
}

// This function (and related functions) does not belong in the motifs namespace
// It also need to be more general, only uses DnaDefs currently
// Make resfile compatible, see what options already exist (JA TARGET strategy)
// Defs assume PDB numbering!!
utility::vector1< core::Size >
get_target_positions_make_dna_mutations(
	core::pose::Pose & pose
)
{
	using namespace basic::options;
	using namespace protocols::dna;
	utility::vector1< core::Size > target_positions(0);
	if ( option[ OptionKeys::motifs::target_dna_defs ].user() &&
				option[ OptionKeys::dna::design::dna_defs ].user() ) {
		// This situation would only arise if you want to make many
		// DNA mutations, but only focus motif targeting on a subset
		// The option target_dna_defs should only be used in this situation
		// At least as I can see now, I could imagine input mutations in another manner
		// And use that def input to direct motifs
		// A second situation that should involve target_dna_defs is the one where
		// multiple bases are allowed a single position, so the Def is in the form of
		// 409.X.N or 409.X.V, using the standard code for degeneracy
		DnaDesignDefOPs mutated_dna;
		DnaDesignDefOPs targeted_dna;
		load_dna_design_defs_from_options( mutated_dna );
		make_dna_mutations( pose, mutated_dna );
		utility::vector1< std::string > str_def( option[ OptionKeys::motifs::target_dna_defs ]().vector() );
		load_dna_design_defs_from_strings( targeted_dna, str_def );
		target_positions = defs2vector( pose, targeted_dna );
	} else if ( option[ OptionKeys::dna::design::dna_defs ].user() ) {
		DnaDesignDefOPs mutated_dna;
		load_dna_design_defs_from_options( mutated_dna );
		make_dna_mutations( pose, mutated_dna );
		target_positions = defs2vector( pose, mutated_dna );
	} else if ( option[ OptionKeys::motifs::target_dna_defs ].user() ) {
		DnaDesignDefOPs targeted_dna;
		mu_tr << "DNA is not being mutated, but alternative bases are allowed, so mutation may occur at a later point, use dna::dna_defs to input positions to mutate in the straightforward manner." << std::endl;
		utility::vector1< std::string > str_def( option[ OptionKeys::motifs::target_dna_defs ]().vector() );
		load_dna_design_defs_from_strings( targeted_dna, str_def );
		target_positions = defs2vector( pose, targeted_dna );
	} else {
		mu_tr << "No input given for DNA target positions, will identify target positions based on input motif building positions." << std::endl;
	}
	return target_positions;
}

// 12-29-09, this is the updated version that can keep track of allowed types, as well as a build position number
// parts of this are sort of ridiculous, like if you are only doing the mutating part, then you don't really need to make a pair vector, since it's really just a vector, same goes for the situation where everything should remain wild-type . . .
// ie, the first two if statements are really not worthwhile, since you wouldn't be making the single mutations if you were allowing all possiblities, doesn't make sense
// this will be modified later, at least right now it is functional and always produces the vector of pairs
// the real problem is that the vector of strings is not useful (and thus the purpose of the pair), and the output should really just be a vector of sizes if you are going to allow for degeneracy at some positions . . . although maybe some positions allow degeneracy, some must be mutated to single position, who knows, there are all sorts of input possibilities, at least I know right now that this should function, even if it isn't perfect . . .
//utility::vector1< std::pair< core::Size, utility::vector1< std::string > > >
std::map< core::Size, std::set< std::string > >
get_target_position_map_make_dna_mutations(
	core::pose::Pose & pose
)
{
	using namespace basic::options;
	using namespace protocols::dna;
	//utility::vector1< std::pair< core::Size, utility::vector1< std::string > > > target_positions(0);
	std::map< core::Size, std::set< std::string > > target_positions;
	if ( option[ OptionKeys::motifs::target_dna_defs ].user() &&
				option[ OptionKeys::dna::design::dna_defs ].user() ) {
		DnaDesignDefOPs mutated_dna;
		DnaDesignDefOPs targeted_dna;
		load_dna_design_defs_from_options( mutated_dna );
		make_dna_mutations( pose, mutated_dna );
		utility::vector1< std::string > str_def( option[ OptionKeys::motifs::target_dna_defs ]().vector() );
		load_dna_design_defs_from_strings( targeted_dna, str_def );
		//target_positions = defs2allowedtypes( pose, targeted_dna );
		target_positions = defs2map( pose, targeted_dna );
	} else if ( option[ OptionKeys::dna::design::dna_defs ].user() ) {
		DnaDesignDefOPs mutated_dna;
		load_dna_design_defs_from_options( mutated_dna );
		make_dna_mutations( pose, mutated_dna );
		//target_positions = defs2allowedtypes( pose, mutated_dna );
		target_positions = defs2map( pose, mutated_dna );
	} else if ( option[ OptionKeys::motifs::target_dna_defs ].user() ) {
		DnaDesignDefOPs targeted_dna;
		mu_tr << "DNA is not being mutated, but alternative bases are allowed, so mutation may occur at a later point, use dna::dna_defs to input positions to mutate in the straightforward manner." << std::endl;
		utility::vector1< std::string > str_def( option[ OptionKeys::motifs::target_dna_defs ]().vector() );
		load_dna_design_defs_from_strings( targeted_dna, str_def );
		//target_positions = defs2allowedtypes( pose, targeted_dna );
		target_positions = defs2map( pose, targeted_dna );
	} else {
		mu_tr << "No input given for DNA target positions, will identify target positions based on input motif building positions." << std::endl;
	}
	return target_positions;
}

void
make_dna_mutations(
	core::pose::Pose & pose
)
{
	using namespace basic::options;
	using namespace protocols::dna;
	if ( option[ OptionKeys::motifs::target_dna_defs ].user() &&
				option[ OptionKeys::dna::design::dna_defs ].user() ) {
		DnaDesignDefOPs mutated_dna;
		DnaDesignDefOPs targeted_dna;
		load_dna_design_defs_from_options( mutated_dna );
		make_dna_mutations( pose, mutated_dna );
		utility::vector1< std::string > str_def( option[ OptionKeys::motifs::target_dna_defs ]().vector() );
		load_dna_design_defs_from_strings( targeted_dna, str_def );
	} else if ( option[ OptionKeys::dna::design::dna_defs ].user() ) {
		DnaDesignDefOPs mutated_dna;
		load_dna_design_defs_from_options( mutated_dna );
		make_dna_mutations( pose, mutated_dna );
	} else if ( option[ OptionKeys::motifs::target_dna_defs ].user() ) {
		DnaDesignDefOPs targeted_dna;
		mu_tr << "DNA is not being mutated, but alternative bases are allowed, so mutation may occur at a later point, use dna::dna_defs to input positions to mutate in the straightforward manner." << std::endl;
		utility::vector1< std::string > str_def( option[ OptionKeys::motifs::target_dna_defs ]().vector() );
		load_dna_design_defs_from_strings( targeted_dna, str_def );
	} else {
			mu_tr << "No input given for DNA target positions, will identify target positions based on input motif building positions." << std::endl;
	}
}

void
make_dna_mutations(
	core::pose::Pose & pose,
	protocols::dna::DnaDesignDefOPs const & target
)
{
	using namespace protocols::dna;
	core::pose::PDBPoseMap const & pdb_pose_map( pose.pdb_info()->pdb2pose() );
	core::scoring::dna::set_base_partner( pose );
	for ( DnaDesignDefOPs::const_iterator def( target.begin() );
			def != target.end(); ++def ) {
		// SHOULD INCLUDE JA checks to ensure that the input is DNA and is the correct strand
		core::Size index( pdb_pose_map.find( (*def)->chain, (*def)->pdbpos ) );
		if ( ! (*def)->name3.empty() ) {
			std::string basepairID( (*def)->name3 );
			//if ( basepairID.length() == 3 ) {
				if( protocols::dna::dna_full_name3( pose.residue(index).name3() ) != basepairID ) {
					make_base_pair_mutation( pose, index, core::chemical::aa_from_name( basepairID ) ); //what if it has an X instead of a name3
				}
		//	} else {
			// in the future this will end up checking to make sure that it = N or whatever other one-letter identifier you want
		//		mu_tr << "DNA is set to be mutated to multiple DNA bases during the motif search" << std::endl;
		//	}
		// Actually, cannot add this functionality to this function, because the input from the dna_defs will affect the later steps with design, so I will need to find a new place to add this functionality, probably have to use target_defs if I want something other than simple one-to-one dna mutations to be made
		} else {
			mu_tr << "DNA was not mutated because input Def did not include a type!" << std::endl;
			// Future: DEF Does not include a type, therefore I should mutate later in protocol to all types of bases for the RMSD comparison
		}
	}
}

utility::vector1< core::Size >
defs2vector(
	core::pose::Pose const & pose,
	protocols::dna::DnaDesignDefOPs const & targets
)
{
	using namespace protocols::dna;
	core::pose::PDBPoseMap const & pdb_pose_map( pose.pdb_info()->pdb2pose() );
	utility::vector1< core::Size > positions;
	for ( DnaDesignDefOPs::const_iterator def( targets.begin() );
			def != targets.end(); ++def ) {
		core::Size index( pdb_pose_map.find( (*def)->chain, (*def)->pdbpos ) );
		positions.push_back( index );
	}
	return positions;
}

utility::vector1< std::pair< core::Size, utility::vector1< std::string > > >
defs2allowedtypes(
	core::pose::Pose const & pose,
	protocols::dna::DnaDesignDefOPs const & targets
)
{
	using namespace protocols::dna;
	core::pose::PDBPoseMap const & pdb_pose_map( pose.pdb_info()->pdb2pose() );
	utility::vector1< std::pair< core::Size, utility::vector1< std::string > > > positions;
	for ( DnaDesignDefOPs::const_iterator def( targets.begin() );
			def != targets.end(); ++def ) {
		core::Size index( pdb_pose_map.find( (*def)->chain, (*def)->pdbpos ) );
		utility::vector1< std::string > names;
		std::string name( (*def)->name3 );
		//if name3 something, then names.push_back();
		if ( name.length() > 1 ) {
			positions.push_back( std::make_pair( index, names ) );
		} else if ( name.length() == 1 ) {
			// this map should maybe go elsewhere (almost certainly)
			std::map < std::string, utility::vector1< std::string > > degeneracycodes(
				utility::tools::make_map(
				std::string("N"), utility::tools::make_vector1( std::string("ADE"), std::string("CYT"), std::string("CYT"), std::string("THY") )
			) );
			names = degeneracycodes[name];
			positions.push_back( std::make_pair( index, names ) );
		} else {
			mu_tr << "All target positions will remain wild-type." << std::endl; // need to make sure code checks that in motif search and really does use wild-type if no other types are considered
		}
	}
	return positions;
}

// THIS FUNCTION CURRENTLY ISN'T GREAT, WON'T WORK FOR THE BUILDPOS AND MULTIPLE AAs YET
std::map< core::Size, std::set< std::string > >
defs2map(
	core::pose::Pose const & pose,
	protocols::dna::DnaDesignDefOPs const & targets
)
{
	using namespace protocols::dna;
	core::pose::PDBPoseMap const & pdb_pose_map( pose.pdb_info()->pdb2pose() );
	std::map< core::Size, std::set< std::string > > positions;
	for ( DnaDesignDefOPs::const_iterator def( targets.begin() );
			def != targets.end(); ++def ) {
		core::Size index( pdb_pose_map.find( (*def)->chain, (*def)->pdbpos ) );
		std::set< std::string > names;
		std::string name( (*def)->name3 );
		if ( name.length() > 1 ) {
			//positions.push_back( std::make_pair( index, names ) );
			names.insert( name );
		} else if ( name.length() == 1 ) {
			// THIS NEEDS TO BE FIXED FOR SET VERSUS VECTOR1 IF I WANT TO USE IT
			// this map should maybe go elsewhere (almost certainly)
			std::map < std::string, utility::vector1< std::string > > degeneracycodes(
				utility::tools::make_map(
				std::string("N"), utility::tools::make_vector1( std::string("ADE"), std::string("CYT"), std::string("CYT"), std::string("THY") )
			) );
			//names = degeneracycodes[name];
		} else {
			mu_tr << "All target positions will remain wild-type." << std::endl; // need to make sure code checks that in motif search and really does use wild-type if no other types are considered
		}
		positions[index] = names;
	}
	return positions;
}

std::map< core::Size, std::set< std::string > >
bpdefs2map(
	core::pose::Pose const & pose,
	protocols::dna::DnaDesignDefOPs const & targets
)
{
	using namespace protocols::dna;
	core::pose::PDBPoseMap const & pdb_pose_map( pose.pdb_info()->pdb2pose() );
	std::map< core::Size, std::set< std::string > > positions;
	for ( DnaDesignDefOPs::const_iterator def( targets.begin() );
			def != targets.end(); ++def ) {
		core::Size index( pdb_pose_map.find( (*def)->chain, (*def)->pdbpos ) );
		std::set< std::string > names;
		std::string name( (*def)->name3 );
		for ( core::Size c(0); c < name.size(); ++c ) {
			mu_tr << "Allowing AAtype " << name[c] << " for motif search." << std::endl;
			// put all canonical AAs in it
			std::stringstream name1;
			name1 << name[c];
			if ( name1.str() == "X" ) {
				utility::vector1< std::string > allAA(utility::tools::make_vector1(std::string("A"), std::string("C"), std::string("D"), std::string("E"), std::string("F"), std::string("H"), std::string("I"), std::string("K"), std::string("L"), std::string("M"), std::string("N"), std::string("P"), std::string("Q"), std::string("R"), std::string("S"), std::string("T"), std::string("V"), std::string("W"), std::string("Y") ) );
				for( core::Size x(1); x <= allAA.size(); ++x ) {
					names.insert( name3_from_oneletter( allAA[x] ) );
				};
			};
			std::string name3( name3_from_oneletter( name1.str() ) );
			names.insert( name3 );
			if (name3 == "GLY") {
				mu_tr << "There are no such thing as glycine motifs, check your build_position_defs and remove G" << std::endl;
			}

		// This function could be also potentially be converted to take a 3 letter code for some residue type that is not a canonical aa with a 1-letter
		}
		positions[index] = names;
	}
	return positions;
}

std::string
name3_from_oneletter(
	std::string const & oneletter
)
{
	std::map< std::string, std::string > name3_name1( utility::tools::make_map(
		std::string("A"), std::string("ALA"),
		std::string("C"), std::string("CYS"),
		std::string("D"), std::string("ASP"),
		std::string("E"), std::string("GLU"),
		std::string("F"), std::string("PHE"),
		std::string("G"), std::string("GLY"),
		std::string("H"), std::string("HIS"),
		std::string("I"), std::string("ILE"),
		std::string("K"), std::string("LYS"),
		std::string("L"), std::string("LEU"),
		std::string("M"), std::string("MET"),
		std::string("N"), std::string("ASN"),
		std::string("P"), std::string("PRO"),
		std::string("Q"), std::string("GLN"),
		std::string("R"), std::string("ARG"),
		std::string("S"), std::string("SER"),
		std::string("T"), std::string("THR"),
		std::string("V"), std::string("VAL"),
		std::string("W"), std::string("TRP"),
		std::string("Y"), std::string("TYR")
	) );
	return name3_name1[ oneletter ];
}

utility::vector1< core::Size >
get_motif_build_positions_user(
	core::pose::Pose const & pose
)
{
	using namespace basic::options;
	using namespace protocols::dna;
	utility::vector1< core::Size > motif_build_positions;
	if ( option[ OptionKeys::motifs::motif_build_defs ].user() ) {
		DnaDesignDefOPs build_positions;  // it's not a DnaDef anymore, since I'm using it for protein . .
		utility::vector1< std::string > str_def( option[ OptionKeys::motifs::motif_build_defs ]().vector() );
		load_dna_design_defs_from_strings( build_positions, str_def );
		motif_build_positions = defs2vector( pose, build_positions );
	} else {
		mu_tr << "No build positions specified by user input, will identify build positions based on proximity to target positions" << std::endl;
	}
	return motif_build_positions;
}

utility::vector1< protocols::dna::DnaDesignDefOP >
get_motif_build_position_defs_user()
{
	using namespace basic::options;
	utility::vector1< protocols::dna::DnaDesignDefOP > motif_build_positions(0);
	using namespace protocols::dna;
	if ( option[ OptionKeys::motifs::motif_build_defs ].user() ) {
		DnaDesignDefOPs build_positions;  // it's not a DnaDef anymore, since I'm using it for protein . .
		utility::vector1< std::string > str_def( option[ OptionKeys::motifs::motif_build_defs ]().vector() );
		load_dna_design_defs_from_strings( build_positions, str_def );
		motif_build_positions = build_positions;
	} else {
		mu_tr << "No build positions specified by user input, will identify build positions based on proximity to target positions" << std::endl;
	}
	return motif_build_positions;
}

//Ligand motif load_build_position_data
//This is a pretty stupid, hacky way to fix the problem.  I should just fix the other load_build_position_data with an _optional_ ligand_marker which allows me to divert ligand motifs into conditional statements.
void
load_build_position_data(
	BuildPosition & bp,
	std::string const & filename,
	core::pose::Pose & pose,
	core::Size const
)
{
	bool keep_one_motif( true );
	std::map< std::string, SingleMotifOP > single_motifs;
	utility::io::izstream data_file( filename.c_str() );
	std::string key_in;
	data_file >> key_in;
	if( key_in == "POSITION" ) {
		while( data_file >> key_in ) {
			std::stringstream bpseqpos;
			bpseqpos << bp.seqpos();
			if( key_in == bpseqpos.str() ) {
				std::string key2_in;
				while( data_file >> key2_in ) {
					if( key2_in == "POSITION" ) {
						std::string key3_in;
						data_file >> key3_in;
						break;
					} else if( key2_in == "SINGLE" ) {
						SingleMotifOP new_motif = single_ligand_motif_from_stream( data_file );
						if( ! keep_one_motif ) {
							bp.keep_motif( *new_motif );//best_motifs()->add_to_library( new_motif );
						} else {
							single_motifs[ new_motif->remark() ] = new_motif;
						}
					} else if( key2_in == "RESIDUE" ) {
					//	std::cout << "Pos2 id: " << key_in << std::endl;
						core::conformation::ResidueOP rsd;
						rsd = ( single_residue_from_stream( data_file )->clone() );
						rsd->seqpos( bp.seqpos() );
						utility::vector1<core::Real> mainchains( pose.residue(bp.seqpos()).mainchain_torsions() );
						rsd->mainchain_torsions( mainchains );
						set_chi_according_to_coordinates( *rsd );
						rsd->copy_residue_connections( pose.residue(bp.seqpos()) );
						bp.keep_rotamer( *rsd );
					}
			}
			break;
			} else {
					continue;
			}
		}
			if( keep_one_motif ) {
				for( std::map< std::string, SingleMotifOP >::iterator mot( single_motifs.begin() ),
						end( single_motifs.end()); mot != end; ++mot ) {
					bp.keep_motif( *(mot->second) );
				}
			}
		} else {
			mu_tr << "This file doesn't have any positions in it" << std::endl;
		}
}

void
load_build_position_data(
	BuildPosition & bp,
	std::string const & filename,
	core::pose::Pose & pose
)
{
	bool keep_one_motif( true );
	std::map< std::string, SingleMotifOP > single_motifs;
	utility::io::izstream data_file( filename.c_str() );
	std::string key_in;
	data_file >> key_in;
	//std::cout << "Pos id: " << key_in << std::endl;
	if( key_in == "POSITION" ) {
		while( data_file >> key_in ) {
		//		std::string seqpos = key_in;
			//data_file >> seqpos;
		//		std::cout << "ASeqpos: " << seqpos << std::endl;
			std::stringstream bpseqpos;
			bpseqpos << bp.seqpos();
			if( key_in == bpseqpos.str() ) {
				//std::cout << "SEQPOS: " << bpseqpos.str() << std::endl;
				//std::cout << "SEQPOS: " << key_in << std::endl;
				std::string key2_in;
				while( data_file >> key2_in ) {
					//std::cout << "Pos2 id: " << key_in << std::endl;
					if( key2_in == "POSITION" ) {
						std::string key3_in;
						data_file >> key3_in;
				//		std::cout << "READY TO BREAK ON: " << key2_in << " " << key3_in << std::endl;
						break;
					} else if( key2_in == "SINGLE" ) {
						SingleMotifOP new_motif = single_motif_from_stream( data_file );
						if( ! keep_one_motif ) {
							bp.keep_motif( *new_motif );//best_motifs()->add_to_library( new_motif );
						} else {
							single_motifs[ new_motif->remark() ] = new_motif;
						}
					} else if( key2_in == "RESIDUE" ) {
					//	std::cout << "Pos2 id: " << key_in << std::endl;
						core::conformation::ResidueOP rsd;
						rsd = ( single_residue_from_stream( data_file )->clone() );
						rsd->seqpos( bp.seqpos() );
						//core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( rsd->type(), *rsd, pose.conformation() ) );
						utility::vector1<core::Real> mainchains( pose.residue(bp.seqpos()).mainchain_torsions() );
						rsd->mainchain_torsions( mainchains );
						set_chi_according_to_coordinates( *rsd );
						rsd->copy_residue_connections( pose.residue(bp.seqpos()) );
						bp.keep_rotamer( *rsd );
					}
			}
			break;
			} else {
					continue;
			}
		}
			if( keep_one_motif ) {
				for( std::map< std::string, SingleMotifOP >::iterator mot( single_motifs.begin() ),
						end( single_motifs.end()); mot != end; ++mot ) {
					bp.keep_motif( *(mot->second) );
				}
			}
		} else {
			mu_tr << "This file doesn't have any positions in it" << std::endl;
		}
}

utility::vector1< utility::file::FileName >
get_filenames(
	utility::vector1< utility::file::FileName >  const & listnames
)
{
	utility::vector1< utility::file::FileName >  names;
	for ( utility::vector1< utility::file::FileName >::const_iterator filename( listnames.begin() );
			filename != listnames.end(); ++filename ) {
		utility::io::izstream list( (*filename).name().c_str() );
		while ( list ) {
			std::string name;
			list >> name;
			names.push_back( name );
		}
	}
	return names;
}


core::conformation::ResidueOP
single_residue_from_stream(
	utility::io::izstream & residue_info
)
{
	std::string resname;
	residue_info >> resname;
	//std::cout << "READING IN: " << resname << std::endl;
	std::string firstline;
	getline( residue_info, firstline );
	core::conformation::ResidueOP rsd = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map( resname ) );
	core::Size nlines( rsd->natoms() );
	for( core::Size i(1); i <= nlines; ++i ) {
		//std::string atomline;
		std::string atomname;
		residue_info >> atomname;
		//getline( residue_info, atomline );
		utility::vector1< std::string > atomwords( utility::string_split( atomname, ':' ) );
		std::string atomname2( atomwords.front() );
	//	std::cout << "atomname: " << atomname << std::endl;
	//	std::cout << "atomname2: " << atomname2 << std::endl;
		std::string skip, x, y, z;
		if( atomname == atomname2 ) {
			residue_info >> skip;
			residue_info >> x;
			residue_info >> y;
			residue_info >> z;
		//	std::cout << skip << x << y << z << std::endl;
		} else {
				residue_info >> x;
				residue_info >> y;
				residue_info >> z;
		}
		core::Vector atomxyz( utility::string2float(x), utility::string2float(y), utility::string2float(z) );
	//numeric::xyzVector atomxyz;
//	residueinfo >> atomxyz;
		rsd->set_xyz( atomname2, atomxyz );
	}
//	std::cout << "TESTING\n" << *rsd << std::endl;
	return rsd;
}

utility::vector1< bool >
bools_from_sizes(
	core::Size const nres,
	utility::vector1< core::Size > const & v
)
{
	utility::vector1< bool > b( nres, false );
	for ( utility::vector1< core::Size >::const_iterator pos= v.begin(), epos= v.end(); pos != epos; ++pos ) b[ *pos ] = true;
	return b;
}

void
make_base_pair_mutation(
	core::pose::Pose & pose,
	core::Size const seqpos,
	core::chemical::AA const & na
)
{
using namespace core::chemical;
using namespace core::conformation;
using namespace core::scoring::dna;

ResidueTypeSet const & residue_set( pose.residue(1).residue_type_set() );
BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );

for ( int r=1; r<= 2; ++r ) {
	core::Size const pos( r == 1 ? seqpos : partner[seqpos] );
	if ( pos == 0 ) continue; // unpaired
	AA const aa( r == 1 ? na : protocols::dna::dna_base_partner( na ) );

	Residue const & existing_residue( pose.residue( pos ) );
	assert( existing_residue.is_DNA() );

	// search for the matching residue type
	ResidueTypeCOPs rsd_types
		( ResidueSelector().set_aa( aa ).match_variants( existing_residue.type() ).select( residue_set ) );
	if ( rsd_types.size() != 1 ) {
		utility_exit_with_message("couldnt find residuetype for basepair mutation!");
	}

	ResidueOP rsd = ResidueFactory::create_residue( *(rsd_types[1]), existing_residue, pose.conformation() );
	rsd->set_chi( 1, existing_residue.chi(1) );

	pose.replace_residue( pos, *rsd, false );
}
}

core::Real
atom_specific_rms(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	utility::vector1< core::Size > const & atoms
)
{
	using namespace core;
	using namespace core::chemical;
	using namespace core::conformation;
	// Doesn't do automorphisms like original core/scoring/rms_util.cc functions!!
	//if( rsd1.type().name3() != rsd2.type().name3() ) utility_exit_with_message("Residue type name3 mismatch");
	//if( rsd1.nheavyatoms()  != rsd2.nheavyatoms()  ) utility_exit_with_message("Residue number-of-heavy-atoms mismatch");
	core::Real best_rms = 1e99;
	// Make atom-number translation table
	core::Real sum2( 0.0 );
	core::Size natoms( 0 );
	for( core::Size j = 1; j <= atoms.size(); ++j ) {
		core::Vector diff = rsd1.xyz( atoms[j] ) - rsd2.xyz( atoms[j] );
		sum2 += diff.length_squared();
		natoms +=1;
	}
	core::Real const curr_rms = std::sqrt(sum2 / natoms);

	// Check vs. minimum rmsd
	if( curr_rms < best_rms ) {
		best_rms = curr_rms;
	}
	return best_rms;
}

core::Real
atom_specific_rms(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	utility::vector1< std::string > const & atoms
)
{
	using namespace core;
	using namespace core::chemical;
	using namespace core::conformation;
	// Doesn't do automorphisms like original core/scoring/rms_util.cc functions!!
	//if( rsd1.type().name3() != rsd2.type().name3() ) utility_exit_with_message("Residue type name3 mismatch");
	//if( rsd1.nheavyatoms()  != rsd2.nheavyatoms()  ) utility_exit_with_message("Residue number-of-heavy-atoms mismatch");
	core::Real best_rms = 1e99;
	// Make atom-number translation table
	core::Real sum2( 0.0 );
	core::Size natoms( 0 );
	for( core::Size j = 1; j <= atoms.size(); ++j ) {
		core::Vector diff = rsd1.xyz( atoms[j] ) - rsd2.xyz( atoms[j] );
		sum2 += diff.length_squared();
		natoms +=1;
	}
	core::Real const curr_rms = std::sqrt(sum2 / natoms);

	// Check vs. minimum rmsd
	if( curr_rms < best_rms ) {
		best_rms = curr_rms;
	}
	return best_rms;
}
core::pack::rotamer_set::RotamerSetOP
build_rotamers_lite(
	core::pose::Pose & pose,
	core::Size const rotamer_build_position,
	utility::vector1< bool > aa_info,
	core::Size const ex_,
	bool bump_check
)
{
	using namespace core::pack::rotamer_set;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::pack;
	using namespace core::pack::task;

	RotamerSetFactory rsf;
	RotamerSetOP rotset = rsf.create_rotamer_set( pose.residue( rotamer_build_position ) );
	rotset->set_resid( rotamer_build_position );

	// Gather up the many things needed to build rotamers
	ScoreFunction scorefxn;
	scorefxn.set_weight( fa_atr, 1.00 );
	scorefxn.set_weight( fa_rep, 1.00 );

	PackerTaskOP task = TaskFactory::create_packer_task( pose );
	task->set_bump_check( bump_check );
	task->temporarily_fix_everything();
	task->temporarily_set_pack_residue( rotamer_build_position, true );

	task->nonconst_residue_task( rotamer_build_position ).restrict_absent_canonical_aas( aa_info );

	//pose.update_residue_neighbors(); //This was necessary before, but let's see if it is now
	//task->nonconst_residue_task( rotamer_build_position ).or_include_current( true );
	// You can't "include_current" because there is no current or native residue
	//task->nonconst_residue_task( rotamer_build_position).or_exrandom_sample_level(core::pack::task::NO_EXTRA_CHI_SAMPLES);
	if( ex_ > 0 ) task->nonconst_residue_task( rotamer_build_position ).or_ex1( true );
	if( ex_ > 1 ) task->nonconst_residue_task( rotamer_build_position ).or_ex2( true );
	if( ex_ > 2 ) task->nonconst_residue_task( rotamer_build_position ).or_ex3( true );
	if( ex_ > 3 ) task->nonconst_residue_task( rotamer_build_position ).or_ex4( true );
	if( ex_ > 4 ) task->nonconst_residue_task( rotamer_build_position ).or_ex1_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
	if( ex_ > 5 ) task->nonconst_residue_task( rotamer_build_position ).or_ex2_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
	if( ex_ > 6 ) task->nonconst_residue_task( rotamer_build_position ).or_ex3_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
	if( ex_ > 7 ) task->nonconst_residue_task( rotamer_build_position ).or_ex4_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);

	scorefxn( pose );
	core::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, scorefxn, task );

	rotset->build_rotamers( pose, scorefxn, *task, packer_neighbor_graph, false );

	return rotset;
}


} // namespace motifs
} // namespace protocols
