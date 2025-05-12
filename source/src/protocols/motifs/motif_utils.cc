// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.functions.hh>
#include <core/conformation/ResidueFactory.hh>
#include <utility/graph/Graph.hh>
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
#include <core/pose/util.hh>
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
#include <core/pack/task/ResidueLevelTask.hh>

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
#include <utility/io/ozstream.hh>

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

static basic::Tracer mu_tr( "protocols.motifs.motif_utils", basic::t_info );

//this function uses the constructor that accepts two residues, rather than input atoms, because if the constructors that allow you to setup the atoms explicitly are necessary you are likely to have the information setup in a motif_file (or you can simply make another matching function that allows input atoms as well as a pdb, rather than a istream)
SingleMotifOP
single_motif_from_filename(
	utility::file::FileName const & motif_filename
)
{
	core::pose::PoseOP pose( new core::pose::Pose );
	core::import_pose::pose_from_file( *pose, motif_filename , core::import_pose::PDB_file);
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

	/*
	//char array, old (keeping as backup)
	core::Size const max_line_len(1024);
	char first_line[max_line_len];
	motif_info.getline( first_line, sizeof(first_line) );
	*/

	std::string first_line;
	motif_info.getline(first_line);
	std::istringstream line_in( first_line );

	line_in >> res1;
	line_in >> r1atom1;
	line_in >> r1atom2;
	line_in >> r1atom3;
	line_in >> res2;
	line_in >> r2atom1;
	line_in >> r2atom2;
	line_in >> r2atom3;

	//char array, old (keeping as backup)
	/*
	char rest_of_line[max_line_len];
	line_in.getline( rest_of_line, sizeof(rest_of_line) );
	*/

	std::string rest_of_line;
	std::getline(line_in, rest_of_line);

	std::string rest( rest_of_line );
	std::string::size_type index( rest.find("REMARK ") );
	std::string this_remark;
	bool has_remark( false );
	if ( index != std::string::npos ) {
		this_remark = rest.substr( index + 7 );
		has_remark = true;
	}

	// Scan to end of line to skip comments
	// motif_info.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	// motif_info.get( remark_in, max_remark );

	motif_info >> motif_jump;

	SingleMotifOP retval( new SingleMotif( res1, r1atom1, r1atom2, r1atom3, res2, r2atom1, r2atom2, r2atom3, motif_jump ) );
	if ( has_remark ) {
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

	/*
	//char array, old (keeping as backup)
	core::Size const max_line_len(1024);
	char first_line[max_line_len];
	motif_info.getline( first_line, sizeof(first_line) );
	*/

	std::string first_line;
	std::getline(motif_info, first_line);
	std::istringstream line_in( first_line );

	line_in >> res1;
	line_in >> r1atom1;
	line_in >> r1atom2;
	line_in >> r1atom3;
	line_in >> res2;
	line_in >> r2atom1;
	line_in >> r2atom2;
	line_in >> r2atom3;

	//char array, old (keeping as backup)
	/*
	char rest_of_line[max_line_len];
	line_in.getline( rest_of_line, sizeof(rest_of_line) );
	*/

	std::string rest_of_line;
	std::getline(line_in, rest_of_line);
	//mu_tr << "Rest of line:" << rest_of_line << std::endl;
	std::string rest( rest_of_line );
	std::string::size_type index( rest.find("REMARK ") );
	std::string this_remark;
	bool has_remark( false );
	if ( index != std::string::npos ) {
		this_remark = rest.substr( index + 7 );
		has_remark = true;
		// mu_tr << "Remark is " << this_remark << std::endl;
	} else {
		// mu_tr << "No remark in motif" << std::endl;
	}

	// Scan to end of line to skip comments
	// motif_info.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	// motif_info.get( remark_in, max_remark );

	motif_info >> motif_jump;

	SingleMotifOP retval( new SingleMotif( res1, r1atom1, r1atom2, r1atom3, res2, r2atom1, r2atom2, r2atom3, motif_jump ) );
	if ( has_remark ) {
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

	/*
	//char array, old (keeping as backup)
	core::Size const max_line_len(1024);
	char first_line[max_line_len];
	motif_info.getline( first_line, sizeof(first_line) );
	*/

	std::string first_line;
	std::getline(motif_info, first_line);
	std::istringstream line_in( first_line );

	line_in >> res1;
	line_in >> r1atom1;
	line_in >> r1atom2;
	line_in >> r1atom3;
	line_in >> res2;
	line_in >> r2atom1;
	line_in >> r2atom2;
	line_in >> r2atom3;

	//char array, old (keeping as backup)
	/*
	char rest_of_line[max_line_len];
	line_in.getline( rest_of_line, sizeof(rest_of_line) );
	*/

	std::string rest_of_line;
	std::getline(line_in, rest_of_line);
	std::string rest( rest_of_line );
	std::string::size_type index( rest.find("REMARK ") );
	std::string this_remark;
	bool has_remark( false );
	if ( index != std::string::npos ) {
		this_remark = rest.substr( index + 7 );
		has_remark = true;
	}
	motif_info >> motif_jump;
	SingleMotifOP retval( new SingleMotif( res1, r1atom1, r1atom2, r1atom3, r2atom1, r2atom2, r2atom3, motif_jump ) );
	if ( has_remark ) {
		retval->store_remark( this_remark );
	}

	return retval;
}

SingleMotifOP
single_ligand_motif_from_stream(
	std::istream & motif_info, bool check_for_bad_motifs
)
{

	/*
	//test for >> overload
	mu_tr << "Calling >> overload" << std::endl;
	SingleMotifOP retval;
	motif_info >> retval;
	return retval;
	*/
	std::string res1;
	std::string r1atom1, r1atom2, r1atom3;
	std::string res2;
	std::string r2atom1, r2atom2, r2atom3;
	core::kinematics::Jump motif_jump;

	/*
	//char array, old (keeping as backup)
	core::Size const max_line_len(1024);
	char first_line[max_line_len];
	motif_info.getline( first_line, sizeof(first_line) );
	*/

	std::string first_line;
	std::getline(motif_info, first_line);
	std::istringstream line_in( first_line );

	line_in >> res1;
	line_in >> r1atom1;
	line_in >> r1atom2;
	line_in >> r1atom3;
	line_in >> res2;
	line_in >> r2atom1;
	line_in >> r2atom2;
	line_in >> r2atom3;

	//char array, old (keeping as backup)
	/*
	char rest_of_line[max_line_len];
	line_in.getline( rest_of_line, sizeof(rest_of_line) );
	*/

	std::string rest_of_line;
	std::getline(line_in, rest_of_line);

	std::string rest( rest_of_line );
	std::string::size_type index( rest.find("REMARK ") );
	std::string this_remark;
	bool has_remark( false );
	if ( index != std::string::npos ) {
		this_remark = rest.substr( index + 7 );
		has_remark = true;
	}


	motif_info >> motif_jump;
	SingleMotifOP retval( new SingleMotif( res1, r1atom1, r1atom2, r1atom3, r2atom1, r2atom2, r2atom3, motif_jump ) );

	//bool indicate whether to recursively call this function if a motif is encountered that fails the orthogonality check
	bool call_next_motif = false;

	//mu_tr << "Stream state: " << motif_info.rdstate() << std::endl;
	//indicate to user that the motif failed the orthogonality test and what the motif is
	if ( motif_info.rdstate() != 0 ) {
		if ( basic::options::option[ basic::options::OptionKeys::motifs::verbosity ]() ) {
			mu_tr << "RT failed orthogonality check!" << std::endl;
			mu_tr << "Failed motif is:" << std::endl;
			mu_tr << res1 << "   " << r1atom1 << "   " << r1atom2 << "   " << r1atom3 << "   " << res2 << "   " << r2atom1 << "   " << r2atom2 << "   " << r2atom3 << "   " << rest << std::endl;
			mu_tr << motif_jump << std::endl;
		}
		call_next_motif = true;
	}

	motif_info.clear();
	//set the state of stream to good
	motif_info.setstate( std::ios_base::goodbit );

	if ( has_remark ) {
		retval->store_remark( this_remark );
	}

	//override call_next_motif behavior, using check_for_bad_motifs variable; if false then return anyway
	if ( check_for_bad_motifs == false ) {
		return retval;
	}

	//recursively call this function if a failed motif was encountered. This will skip the failed motif and return the next one in the file
	if ( call_next_motif ) {
		//read in the next "SINGLE"
		std::string key_in;
		motif_info >> key_in;

		retval = single_ligand_motif_from_stream(motif_info, check_for_bad_motifs);
	}

	return retval;
}

core::Real
parallel_base_test(
	core::conformation::Residue const & pose_dna,
	core::conformation::Residue const & motif_dna
)
{
	using xyzVec = numeric::xyzVector<core::Real>;
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
	if ( r1.type().aa() == core::chemical::aa_gly || r2.type().aa() == core::chemical::aa_gly ) {
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
	cst_set->add_constraint( ConstraintCOP( utility::pointer::make_shared< CoordinateConstraint >( core::id::AtomID( pose.residue( this_pos ).atom_index( "CA" ), this_pos ),
		core::id::AtomID( pose.residue( 1 ).atom_index( "CA" ), 1 ),
		inv_rotamer.xyz( "CA" ),
		fx1 ) ) );

	FuncOP fx2( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	cst_set->add_constraint( ConstraintCOP( utility::pointer::make_shared< CoordinateConstraint >( core::id::AtomID( pose.residue( this_pos ).atom_index( "C" ), this_pos ),
		core::id::AtomID( pose.residue( 1 ).atom_index( "CA" ), 1 ),
		inv_rotamer.xyz( "C" ),
		fx2 ) ) );

	FuncOP fx3( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	cst_set->add_constraint( ConstraintCOP( utility::pointer::make_shared< CoordinateConstraint >( core::id::AtomID( pose.residue( this_pos ).atom_index( "N" ), this_pos ),
		core::id::AtomID( pose.residue( 1 ).atom_index( "CA" ), 1 ),
		inv_rotamer.xyz( "N" ),
		fx3 ) ) );

	// For the fourth atom, if either residue is glycine, you need to use HA, else use CB
	if ( pose.residue( this_pos ).type().aa() == core::chemical::aa_gly || inv_rotamer.type().aa() == core::chemical::aa_gly ) {
		core::Size index1( pose.residue( this_pos ).type().aa() == core::chemical::aa_gly ?
			pose.residue( this_pos ).atom_index( "1HA" ) : pose.residue( this_pos ).atom_index( "HA" ) );
		core::Size index2( inv_rotamer.type().aa() == core::chemical::aa_gly ?
			inv_rotamer.atom_index( "1HA" ) : inv_rotamer.atom_index( "HA" ) );

		FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
		cst_set->add_constraint( ConstraintCOP( utility::pointer::make_shared< CoordinateConstraint >( core::id::AtomID( index1, this_pos ),
			core::id::AtomID( pose.residue( 1 ).atom_index( "CA" ), 1 ),
			inv_rotamer.xyz( index2 ),
			fx ) ) );
	} else {
		FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
		cst_set->add_constraint( ConstraintCOP( utility::pointer::make_shared< CoordinateConstraint >( core::id::AtomID( pose.residue( this_pos ).atom_index( "CB" ), this_pos ),
			core::id::AtomID( pose.residue( 1 ).atom_index( "CA" ), 1 ),
			inv_rotamer.xyz( "CB" ),
			fx ) ) );
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
	core::Size nres( pose.size() );
	for ( core::Size i(1); i <= nres; ++i ) {
		if ( pose.residue_type(i).is_protein() ) {
			first_protein_resi = i;
			break;
		}
	}

	FuncOP fx1( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	cst_set->add_constraint( ConstraintCOP( utility::pointer::make_shared< CoordinateConstraint >( core::id::AtomID( index1, this_pos ),
		core::id::AtomID( pose.residue( first_protein_resi ).atom_index( "CA" ), 1 ),
		inv_rotamer.xyz( index1 ),
		fx1 ) ) );

	FuncOP fx2( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	cst_set->add_constraint( ConstraintCOP( utility::pointer::make_shared< CoordinateConstraint >( core::id::AtomID( index2, this_pos ),
		core::id::AtomID( pose.residue( first_protein_resi ).atom_index( "CA" ), 1 ),
		inv_rotamer.xyz( index2 ),
		fx2 ) ) );

	FuncOP fx3( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
	cst_set->add_constraint( ConstraintCOP( utility::pointer::make_shared< CoordinateConstraint >( core::id::AtomID( index3, this_pos ),
		core::id::AtomID( pose.residue( first_protein_resi ).atom_index( "CA" ), 1 ),
		inv_rotamer.xyz( index3 ),
		fx3 ) ) );

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

	for ( core::Size ires = 1, end_i = pose.size() ; ires <= end_i ; ++ires ) {

		if ( find( trim_positions.begin(), trim_positions.end(), ires ) != trim_positions.end() ) {
			if ( pose.residue( ires ).aa() != core::chemical::aa_pro &&
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

	for ( core::Size ires = 1, end_res = pose.size() ; ires <= end_res ; ++ires ) {
		if ( flex_regions.is_loop_residue( ires ) ) {
			if ( pose.residue( ires ).aa() != core::chemical::aa_pro &&
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



MotifLibrary const
get_LigandMotifLibrary_user(bool check_for_bad_motifs,  utility::vector1< std::string > const & ligand_atom_names)
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
		motifs.add_ligand_from_file( motif_filename, check_for_bad_motifs, ligand_atom_names ); //This is main change to continue on ligand path
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
			core::import_pose::pose_from_file( *pose, *filename , core::import_pose::PDB_file);
			if ( pose->size() > 1 ) {
				std::cerr << "WARNING!!! Conformer PDB contains more than one residue, loading all residues in PDB as conformers." << std::endl;
			}
			for ( core::Size i(1); i <= pose->size(); ++i ) {
				conformerOPs.push_back( pose->residue(i).clone() );
			}
		}
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
	for ( auto const & conformerOP : conformerOPs ) {
		std::string name( conformerOP->name3() );
		conformer_map[name].push_back( conformerOP );
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
	for ( auto const & def : target ) {
		// SHOULD INCLUDE JA checks to ensure that the input is DNA and is the correct strand
		core::Size index( pdb_pose_map.find( def->chain, def->pdbpos ) );
		if ( ! def->name3.empty() ) {
			std::string basepairID( def->name3 );
			//if ( basepairID.length() == 3 ) {
			if ( protocols::dna::dna_full_name3( pose.residue(index).name3() ) != basepairID ) {
				make_base_pair_mutation( pose, index, core::chemical::aa_from_name( basepairID ) ); //what if it has an X instead of a name3
			}
			// } else {
			// in the future this will end up checking to make sure that it = N or whatever other one-letter identifier you want
			//  mu_tr << "DNA is set to be mutated to multiple DNA bases during the motif search" << std::endl;
			// }
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
	for ( auto const & target : targets ) {
		core::Size index( pdb_pose_map.find( target->chain, target->pdbpos ) );
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
	for ( auto const & target : targets ) {
		core::Size index( pdb_pose_map.find( target->chain, target->pdbpos ) );
		utility::vector1< std::string > names;
		std::string name( target->name3 );
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
	for ( auto const & target : targets ) {
		core::Size index( pdb_pose_map.find( target->chain, target->pdbpos ) );
		std::set< std::string > names;
		std::string name( target->name3 );
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
	for ( auto const & target : targets ) {
		core::Size index( pdb_pose_map.find( target->chain, target->pdbpos ) );
		std::set< std::string > names;
		std::string name( target->name3 );
		for ( char c : name ) {
			mu_tr << "Allowing AAtype " << c << " for motif search." << std::endl;
			// put all canonical AAs in it
			std::stringstream name1;
			name1 << c;
			if ( name1.str() == "X" ) {
				utility::vector1< std::string > allAA(utility::tools::make_vector1(std::string("A"), std::string("C"), std::string("D"), std::string("E"), std::string("F"), std::string("H"), std::string("I"), std::string("K"), std::string("L"), std::string("M"), std::string("N"), std::string("P"), std::string("Q"), std::string("R"), std::string("S"), std::string("T"), std::string("V"), std::string("W"), std::string("Y") ) );
				for ( core::Size x(1); x <= allAA.size(); ++x ) {
					names.insert( name3_from_oneletter( allAA[x] ) );
				};
			};
			std::string name3( name3_from_oneletter( name1.str() ) );
			names.insert( name3 );
			if ( name3 == "GLY" ) {
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

void
load_build_position_data(
	BuildPosition & bp,
	std::string const & filename,
	core::pose::Pose & pose,
	core::Size const ligand_marker
)
{
	//mu_tr << "In 4 parameter load_build_position_data" << std::endl;
	//mu_tr << "Ligand_marker is" << ligand_marker << std::endl;

	std::map< std::string, SingleMotifOP > single_motifs;
	utility::io::izstream data_file( filename.c_str() );
	std::string key_in;
	data_file >> key_in;
	if ( key_in == "POSITION" ) {
		bool keep_one_motif( true );
		while ( data_file >> key_in ) {
			std::stringstream bpseqpos;
			bpseqpos << bp.seqpos();
			if ( key_in == bpseqpos.str() ) {
				std::string key2_in;
				while ( data_file >> key2_in ) {
					if ( key2_in == "POSITION" ) {
						std::string key3_in;
						data_file >> key3_in;
						break;
					} else if ( key2_in == "SINGLE" ) {

						SingleMotifOP new_motif = ( ligand_marker == LIGAND ) ?
							single_ligand_motif_from_stream( data_file ) :
							single_motif_from_stream( data_file );

						if ( ! keep_one_motif ) {
							bp.keep_motif( *new_motif );
						} else {
							single_motifs[ new_motif->remark() ] = new_motif;
						}
					} else if ( key2_in == "RESIDUE" ) {
						core::conformation::ResidueOP rsd;
						rsd = ( single_residue_from_stream( data_file )->clone() );
						rsd->seqpos( bp.seqpos() );
						utility::vector1<core::Real> mainchains( pose.residue(bp.seqpos()).mainchain_torsions() );
						rsd->mainchain_torsions( mainchains );
						set_chi_according_to_coordinates( *rsd );
						rsd->copy_residue_connections( pose.residue(bp.seqpos()) );
						//mu_tr << "Keeping rotamer" << std::endl;
						bp.keep_rotamer( *rsd );

						/*
						core::pack::rotamer_set::Rotamers bp_best_rotamers1( bp.best_rotamers() );
						core::Size bp_rots1( bp_best_rotamers1.size() );
						for ( core::Size r(1); r <= bp_rots1; ++r ) {
						mu_tr << "Rotamer " << r << ": " << bp_best_rotamers1[r]->name3() << std::endl;
						}
						*/

					}
				}
				break;
			} else {
				continue;
			}
		}
		if ( keep_one_motif ) {
			for ( auto & single_motif : single_motifs ) {
				bp.keep_motif( *(single_motif.second) );
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
	for ( auto const & listname : listnames ) {
		utility::io::izstream list( listname.name().c_str() );
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

	//Keep reading in from stream until a residue 3 letter code is read in
	bool is_residue = false;

	while ( is_residue == false )
			{
		residue_info >> resname;
		//mu_tr << "READING IN: " << resname << std::endl;

		//substring handling for D amino acids
		//entire process breaks if you don't feed in a canonical residue
		if ( resname.substr(0,3) == "ALA" || resname.substr(0,3) == "ARG" || resname.substr(0,3) == "ASN" || resname.substr(0,3) == "ASP" || resname.substr(0,3) == "CYS" ||
				resname.substr(0,3) == "GLN" || resname.substr(0,3) == "GLU" || resname.substr(0,3) == "GLY" || resname.substr(0,3) == "HIS" || resname.substr(0,3) == "ILE" ||
				resname.substr(0,3) == "LEU" || resname.substr(0,3) == "LYS" || resname.substr(0,3) == "MET" || resname.substr(0,3) == "PHE" || resname.substr(0,3) == "PRO" ||
				resname.substr(0,3) == "SER" || resname.substr(0,3) == "THR" || resname.substr(0,3) == "TRP" || resname.substr(0,3) == "TYR" || resname.substr(0,3) == "VAL" ) {
			is_residue  = true;
		}

	}


	std::string firstline;
	getline( residue_info, firstline );
	//mu_tr << "firstline: " << firstline << std::endl;
	core::conformation::ResidueOP rsd = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map( resname ) );
	core::Size nlines( rsd->natoms() );

	//iterate through non-rotamer text until the string "Coordinates:" is encountered
	bool before_coordinates = false;
	std::string next_word;
	while ( before_coordinates == false )
			{
		residue_info >> next_word;
		//mu_tr << "next_word: " << next_word << std::endl;
		if ( next_word == "Coordinates:" ) {
			before_coordinates = true;
		}
	}

	for ( core::Size i(1); i <= nlines; ++i ) {
		//std::string atomline;
		//mu_tr << "At new line" << std::endl;
		std::string atomname;
		residue_info >> atomname;

		//get next word if atomname is "(virtual)"
		if ( atomname == "(virtual)" ) {
			residue_info >> atomname;
		}

		//getline( residue_info, atomline );
		utility::vector1< std::string > atomwords( utility::string_split( atomname, ':' ) );
		std::string atomname2( atomwords.front() );
		// mu_tr << "atomname: " << atomname << std::endl;
		// mu_tr << "atomname2: " << atomname2 << std::endl;
		std::string skip, x, y, z;
		if ( atomname == atomname2 ) {
			residue_info >> skip;
			residue_info >> x;
			residue_info >> y;
			residue_info >> z;
			// mu_tr << skip << x << y << z << std::endl;
		} else {
			residue_info >> x;
			residue_info >> y;
			residue_info >> z;
		}
		//mu_tr << "At new line" << std::endl;
		//mu_tr << "atomname: " << atomname << std::endl;
		//mu_tr << "atomname2: " << atomname2 << std::endl;
		//mu_tr << "x: " << x << std::endl;
		//mu_tr << "y: " << y << std::endl;
		//mu_tr << "z: " << z << std::endl;

		core::Vector atomxyz( utility::string2float(x), utility::string2float(y), utility::string2float(z) );
		//numeric::xyzVector atomxyz;
		// residueinfo >> atomxyz;
		//mu_tr << "Before set_xyz call" << std::endl;
		rsd->set_xyz( atomname2, atomxyz );
		//mu_tr << "After set_xyz call" << std::endl;
	}
	// mu_tr << "TESTING\n" << *rsd << std::endl;
	return rsd;
}

utility::vector1< bool >
bools_from_sizes(
	core::Size const nres,
	utility::vector1< core::Size > const & v
)
{
	utility::vector1< bool > b( nres, false );
	for ( core::Size pos : v ) b[ pos ] = true;
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

	ResidueTypeSetCOP residue_set( pose.residue_type_set_for_pose() );
	BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );

	for ( int r=1; r<= 2; ++r ) {
		core::Size const pos( r == 1 ? seqpos : partner[seqpos] );
		if ( pos == 0 ) continue; // unpaired
		AA const aa( r == 1 ? na : protocols::dna::dna_base_partner( na ) );

		Residue const & existing_residue( pose.residue( pos ) );
		debug_assert( existing_residue.is_DNA() );

		// search for the matching residue type
		ResidueTypeCOP rsd_type( residue_set->get_representative_type_aa( aa, existing_residue.type().variant_types() ) );
		if ( rsd_type == nullptr ) {
			utility_exit_with_message("couldnt find residuetype for basepair mutation!");
		}

		ResidueOP rsd = ResidueFactory::create_residue( *rsd_type, existing_residue, pose.conformation() );
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
	for ( core::Size j = 1; j <= atoms.size(); ++j ) {
		core::Vector diff = rsd1.xyz( atoms[j] ) - rsd2.xyz( atoms[j] );
		sum2 += diff.length_squared();
		natoms +=1;
	}
	core::Real const curr_rms = std::sqrt(sum2 / natoms);

	// Check vs. minimum rmsd
	if ( curr_rms < best_rms ) {
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
	for ( core::Size j = 1; j <= atoms.size(); ++j ) {
		core::Vector diff = rsd1.xyz( atoms[j] ) - rsd2.xyz( atoms[j] );
		sum2 += diff.length_squared();
		natoms +=1;
	}
	core::Real const curr_rms = std::sqrt(sum2 / natoms);

	// Check vs. minimum rmsd
	if ( curr_rms < best_rms ) {
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

	RotamerSetOP rotset = RotamerSetFactory::create_rotamer_set( pose );
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
	if ( ex_ > 0 ) task->nonconst_residue_task( rotamer_build_position ).or_ex1( true );
	if ( ex_ > 1 ) task->nonconst_residue_task( rotamer_build_position ).or_ex2( true );
	if ( ex_ > 2 ) task->nonconst_residue_task( rotamer_build_position ).or_ex3( true );
	if ( ex_ > 3 ) task->nonconst_residue_task( rotamer_build_position ).or_ex4( true );
	if ( ex_ > 4 ) task->nonconst_residue_task( rotamer_build_position ).or_ex1_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
	if ( ex_ > 5 ) task->nonconst_residue_task( rotamer_build_position ).or_ex2_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
	if ( ex_ > 6 ) task->nonconst_residue_task( rotamer_build_position ).or_ex3_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
	if ( ex_ > 7 ) task->nonconst_residue_task( rotamer_build_position ).or_ex4_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);

	scorefxn( pose );
	utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, scorefxn, task );

	rotset->build_rotamers( pose, scorefxn, *task, packer_neighbor_graph, false );

	return rotset;
}

std::istream & operator >> (
	std::istream & motif_info, SingleMotifOP & retval
)
{

	std::string res1;
	std::string r1atom1, r1atom2, r1atom3;
	std::string res2;
	std::string r2atom1, r2atom2, r2atom3;
	core::kinematics::Jump motif_jump;

	/*
	//char array, old (keeping as backup)
	core::Size const max_line_len(1024);
	char first_line[max_line_len];
	motif_info.getline( first_line, sizeof(first_line) );
	*/

	std::string first_line;
	std::getline(motif_info, first_line);
	std::istringstream line_in( first_line );

	line_in >> res1;
	line_in >> r1atom1;
	line_in >> r1atom2;
	line_in >> r1atom3;
	line_in >> res2;
	line_in >> r2atom1;
	line_in >> r2atom2;
	line_in >> r2atom3;

	//char array, old (keeping as backup)
	/*
	char rest_of_line[max_line_len];
	line_in.getline( rest_of_line, sizeof(rest_of_line) );
	*/

	std::string rest_of_line;
	std::getline(line_in, rest_of_line);

	std::string rest( rest_of_line );
	std::string::size_type index( rest.find("REMARK ") );
	std::string this_remark;
	bool has_remark( false );
	if ( index != std::string::npos ) {
		this_remark = rest.substr( index + 7 );
		has_remark = true;
	}

	motif_info >> motif_jump;

	SingleMotifOP helper( new SingleMotif( res1, r1atom1, r1atom2, r1atom3, r2atom1, r2atom2, r2atom3, motif_jump ) );

	retval = helper;


	//bool indicate whether to recursively call this function if a motif is encountered that fails the orthogonality check
	//bool call_next_motif = false;

	//mu_tr << "Stream state: " << motif_info.rdstate() << std::endl;
	//indicate to user that the motif failed the orthogonality test and what the motif is
	if ( motif_info.rdstate() != 0 ) {
		//if ( basic::options::option[ basic::options::OptionKeys::motifs::verbosity ]() ) {
		if ( true ) {
			mu_tr << "RT failed orthogonality check!" << std::endl;
			mu_tr << "Failed motif is:" << std::endl;
			mu_tr << res1 << "   " << r1atom1 << "   " << r1atom2 << "   " << r1atom3 << "   " << res2 << "   " << r2atom1 << "   " << r2atom2 << "   " << r2atom3 << "   " << rest << std::endl;
			mu_tr << motif_jump << std::endl;
		}
		//call_next_motif = true;
	}

	motif_info.clear();
	//set the state of stream to good
	motif_info.setstate( std::ios_base::goodbit );

	if ( has_remark ) {
		retval->store_remark( this_remark );
	}

	return motif_info;
}

// @brief writes a MotifLibrary to a .motifs file; converts the MotifLibrary to motifCOPs and then calls the overload that uses motifCOPs for simplicity
void
write_motifs_to_disk(MotifLibrary ml, std::string filename)
{
	write_motifs_to_disk(ml.library() , filename);
}

// @brief motifCOPS to a .motifs file
void
write_motifs_to_disk(MotifCOPs motifcops, std::string filename)
{
	//if there is not period in the file name, append ".motifs" to the end
	if ( filename.find('.') == std::string::npos ) {
		filename = filename + ".motifs";
	}

	utility::io::ozstream motif_output_file( filename );

	for ( protocols::motifs::MotifCOP & motifcop: motifcops ) {
		motif_output_file << *motifcop;
	}
}

// @brief function to hash out an input motif library into a standard map of motifCOPs
//inputs are initial motif library and map that is to be filled out
//map keys are tuples of 7 strings, which is the residue involved in the motif and then the names of the atoms involved (3 atoms on both sides of motif; we don't care about ligand name in key)
//
void hash_motif_library_into_map(protocols::motifs::MotifCOPs & input_library, std::map<motif_atoms, protocols::motifs::MotifCOPs> & mymap)
{
	//iterate over the input library
	for ( auto motifcop : input_library ) {
		//collect residue name from motifcop
		std::string motif_residue_name(motifcop->restype_name1());
		//declare tuple key
		motif_atoms key_tuple(motifcop->restype_name1(),motifcop->res1_atom1_name(),motifcop->res1_atom2_name(),motifcop->res1_atom3_name(),motifcop->res2_atom1_name(),motifcop->res2_atom2_name(),motifcop->res2_atom3_name());

		//push back motif to motifcops at key address
		mymap[key_tuple].push_back(motifcop);

	}
}


} // namespace motifs
} // namespace protocols
