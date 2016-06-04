// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/motifs/Motif.cc
/// @brief Implementation of interaction motifs
/// @author havranek, sthyme (sthyme@gmail.com)

// Unit Headers
#include <protocols/motifs/Motif.hh>

// Package Headers

// Project Headers
#include <protocols/dna/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/graph/Graph.hh>
#include <core/kinematics/Jump.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AtomType.hh> //Need this to prevent the compiling error: invalid use of incomplete type 'const struct core::chemical::AtomType
#include <core/chemical/AtomTypeSet.hh>

// Utility Headers
#include <utility/string_util.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/tools/make_map.hh>

#include <numeric/xyz.functions.hh>

// C++ Headers
#include <sstream>

#include <utility/vector1.hh>


namespace protocols {
namespace motifs {

static THREAD_LOCAL basic::Tracer mt( "protocols.motifs.Motif", basic::t_info );

Motif::Motif(
	std::string const & resname1,
	std::string const & res1_atom1,
	std::string const & res1_atom2,
	std::string const & res1_atom3,
	std::string const & resname2,
	std::string const & res2_atom1,
	std::string const & res2_atom2,
	std::string const & res2_atom3,
	core::kinematics::Jump const & orientation
) : restype_name1_( protocols::dna::dna_full_name3( resname1 ) ),
	res1_atom1_name_( res1_atom1 ),
	res1_atom2_name_( res1_atom2 ),
	res1_atom3_name_( res1_atom3 ),
	restype_name2_( protocols::dna::dna_full_name3( resname2 ) ),
	res2_atom1_name_( res2_atom1 ),
	res2_atom2_name_( res2_atom2 ),
	res2_atom3_name_( res2_atom3 ),
	forward_jump_( orientation ),
	backward_jump_( ( orientation.reversed() ) ),
	has_remark_( false ),
	has_path_( false )
{
	core::chemical::ResidueTypeSet & rsd_set( core::chemical::ChemicalManager::get_instance()->nonconst_residue_type_set( core::chemical::FA_STANDARD ) );
	core::chemical::ResidueType const & rsd_type( rsd_set.name_map( restype_name1_ ) );
	core::chemical::ResidueType const & rsd_type2( rsd_set.name_map( restype_name2_ ) );
	res1_atom1_index_ = rsd_type.atom_index( res1_atom1_name_ );
	res1_atom2_index_ = rsd_type.atom_index( res1_atom2_name_ );
	res1_atom3_index_ = rsd_type.atom_index( res1_atom3_name_ );
	res2_atom1_index_ = rsd_type2.atom_index( res2_atom1_name_ );
	res2_atom2_index_ = rsd_type2.atom_index( res2_atom2_name_ );
	res2_atom3_index_ = rsd_type2.atom_index( res2_atom3_name_ );
}

Motif::Motif(
	core::pose::Pose const & pose,
	Size const residue_position_1,
	char const chain1,
	std::string const & res1_atom1_name,
	std::string const & res1_atom2_name,
	std::string const & res1_atom3_name,
	Size const residue_position_2,
	char const chain2,
	std::string const & res2_atom1_name,
	std::string const & res2_atom2_name,
	std::string const & res2_atom3_name
) : restype_name1_( protocols::dna::dna_full_name3( pose.residue( pose.pdb_info()->pdb2pose( chain1, residue_position_1 ) ).name3() ) ),
	res1_atom1_name_( res1_atom1_name ),
	res1_atom2_name_( res1_atom2_name ),
	res1_atom3_name_( res1_atom3_name ),
	restype_name2_( protocols::dna::dna_full_name3( pose.residue( pose.pdb_info()->pdb2pose( chain2, residue_position_2 ) ).name3() ) ),
	res2_atom1_name_( res2_atom1_name ),
	res2_atom2_name_( res2_atom2_name ),
	res2_atom3_name_( res2_atom3_name ),
	forward_jump_( core::kinematics::Jump(
	core::kinematics::Stub(  pose.residue( pose.pdb_info()->pdb2pose( chain1, residue_position_1 ) ).atom( res1_atom2_name ).xyz(),
	pose.residue( pose.pdb_info()->pdb2pose( chain1, residue_position_1 ) ).atom( res1_atom1_name ).xyz(),
	pose.residue( pose.pdb_info()->pdb2pose( chain1, residue_position_1 ) ).atom( res1_atom2_name ).xyz(),
	pose.residue( pose.pdb_info()->pdb2pose( chain1, residue_position_1 ) ).atom( res1_atom3_name ).xyz() ),
	core::kinematics::Stub(  pose.residue( pose.pdb_info()->pdb2pose( chain2, residue_position_2 ) ).atom( res2_atom2_name ).xyz(),
	pose.residue( pose.pdb_info()->pdb2pose( chain2, residue_position_2 ) ).atom( res2_atom1_name ).xyz(),
	pose.residue( pose.pdb_info()->pdb2pose( chain2, residue_position_2 ) ).atom( res2_atom2_name ).xyz(),
	pose.residue( pose.pdb_info()->pdb2pose( chain2, residue_position_2 ) ).atom( res2_atom3_name ).xyz() ) ) ),
	backward_jump_( forward_jump_.reversed() ),
	has_remark_( false ),
	has_path_( false )
{
	core::chemical::ResidueTypeSet & rsd_set( core::chemical::ChemicalManager::get_instance()->nonconst_residue_type_set( core::chemical::FA_STANDARD ) );
	core::chemical::ResidueType const & rsd_type( rsd_set.name_map( restype_name1_ ) );
	core::chemical::ResidueType const & rsd_type2( rsd_set.name_map( restype_name2_ ) );
	res1_atom1_index_ = rsd_type.atom_index( res1_atom1_name_ );
	res1_atom2_index_ = rsd_type.atom_index( res1_atom2_name_ );
	res1_atom3_index_ = rsd_type.atom_index( res1_atom3_name_ );
	res2_atom1_index_ = rsd_type2.atom_index( res2_atom1_name_ );
	res2_atom2_index_ = rsd_type2.atom_index( res2_atom2_name_ );
	res2_atom3_index_ = rsd_type2.atom_index( res2_atom3_name_ );
}

Motif::Motif(
	core::pose::Pose const & pose,
	Size const residue_position_1,
	std::string const & res1_atom1_name,
	std::string const & res1_atom2_name,
	std::string const & res1_atom3_name,
	Size const residue_position_2,
	std::string const & res2_atom1_name,
	std::string const & res2_atom2_name,
	std::string const & res2_atom3_name
) : restype_name1_( protocols::dna::dna_full_name3( pose.residue( residue_position_1 ).name3() ) ),
	res1_atom1_name_( res1_atom1_name ),
	res1_atom2_name_( res1_atom2_name ),
	res1_atom3_name_( res1_atom3_name ),
	restype_name2_( protocols::dna::dna_full_name3( pose.residue( residue_position_2 ).name3() ) ),
	res2_atom1_name_( res2_atom1_name ),
	res2_atom2_name_( res2_atom2_name ),
	res2_atom3_name_( res2_atom3_name ),
	forward_jump_( core::kinematics::Jump(
	core::kinematics::Stub(  pose.residue( residue_position_1 ).atom( res1_atom2_name ).xyz(),
	pose.residue( residue_position_1 ).atom( res1_atom1_name ).xyz(),
	pose.residue( residue_position_1 ).atom( res1_atom2_name ).xyz(),
	pose.residue( residue_position_1 ).atom( res1_atom3_name ).xyz() ),
	core::kinematics::Stub(  pose.residue( residue_position_2 ).atom( res2_atom2_name ).xyz(),
	pose.residue( residue_position_2 ).atom( res2_atom1_name ).xyz(),
	pose.residue( residue_position_2 ).atom( res2_atom2_name ).xyz(),
	pose.residue( residue_position_2 ).atom( res2_atom3_name ).xyz() ) ) ),
	backward_jump_( forward_jump_.reversed() ),
	has_remark_( false ),
	has_path_( false )
{
	core::chemical::ResidueTypeSet & rsd_set( core::chemical::ChemicalManager::get_instance()->nonconst_residue_type_set( core::chemical::FA_STANDARD ) );
	core::chemical::ResidueType const & rsd_type( rsd_set.name_map( restype_name1_ ) );
	core::chemical::ResidueType const & rsd_type2( rsd_set.name_map( restype_name2_ ) );
	res1_atom1_index_ = rsd_type.atom_index( res1_atom1_name_ );
	res1_atom2_index_ = rsd_type.atom_index( res1_atom2_name_ );
	res1_atom3_index_ = rsd_type.atom_index( res1_atom3_name_ );
	res2_atom1_index_ = rsd_type2.atom_index( res2_atom1_name_ );
	res2_atom2_index_ = rsd_type2.atom_index( res2_atom2_name_ );
	res2_atom3_index_ = rsd_type2.atom_index( res2_atom3_name_ );
}

// Matt's new constructor: This constructor gets the jump for you - Rosetta numbering
// Residue 1 is protein, residue 2 is ligand. For residue 1, just get the residue. For residue 2, get the residue and the 3 atom numbers as a vector.
Motif::Motif(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	utility::vector1< Size > const & res2_atoms
):
	has_remark_( false ),
	has_path_( false )
{
	if ( (res1.is_protein() && res2.is_ligand() ) || (res1.is_protein() && res2.is_protein() ) ) {
		restype_name1_ = res1.name3();
		mt << "Res1 name is " <<  restype_name1_ << std::endl;
		res1_atom1_name_ = motifAtomIDs[restype_name1_][1];
		res1_atom2_name_ = motifAtomIDs[restype_name1_][2];
		res1_atom3_name_ = motifAtomIDs[restype_name1_][3];
		restype_name2_ = "LG1"; //We need to give the second residue (ligand residue) some random name
		core::chemical::AtomType res2_atom1_type = res2.atom_type(res2_atoms[1]);
		core::chemical::AtomType res2_atom2_type = res2.atom_type(res2_atoms[2]);
		core::chemical::AtomType res2_atom3_type = res2.atom_type(res2_atoms[3]);

		res2_atom1_name_ = res2_atom1_type.atom_type_name();
		mt << "atom1 name is " <<  res2_atom1_name_  << std::endl;
		res2_atom2_name_ = res2_atom2_type.atom_type_name();
		mt << "atom2 name is " <<   res2_atom2_name_  << std::endl;
		res2_atom3_name_ = res2_atom3_type.atom_type_name();
		mt << "atom3 name is " <<   res2_atom3_name_  << std::endl;
		//Can't use motif atom name for
		forward_jump_ = core::kinematics::Jump(
			core::kinematics::Stub(  res1.atom( res1_atom2_name_ ).xyz(),
			res1.atom( res1_atom1_name_ ).xyz(),
			res1.atom( res1_atom2_name_ ).xyz(),
			res1.atom( res1_atom3_name_ ).xyz() ),
			core::kinematics::Stub(  res2.atom( res2_atoms[2]  ).xyz(),
			res2.atom( res2_atoms[1] ).xyz(),
			res2.atom( res2_atoms[2] ).xyz(),
			res2.atom( res2_atoms[3] ).xyz() ) );
		backward_jump_ = forward_jump_.reversed();

	} else {
		mt << "Input motif residues do not match the expected combinations for this constructor; please use constructor that allows you to explicitly specify atoms involved." << std::endl;
	}
	core::chemical::ResidueTypeSet & rsd_set( core::chemical::ChemicalManager::get_instance()->nonconst_residue_type_set( core::chemical::FA_STANDARD ) );
	core::chemical::ResidueType const & rsd_type1( rsd_set.name_map( restype_name1_ ) );
	//core::chemical::ResidueType const & rsd_type2( rsd_set.name_map( restype_name2_ ) );
	res1_atom1_index_ = rsd_type1.atom_index( res1_atom1_name_ );
	res1_atom2_index_ = rsd_type1.atom_index( res1_atom2_name_ );
	res1_atom3_index_ = rsd_type1.atom_index( res1_atom3_name_ );
	res2_atom1_index_ = res2_atoms[1];
	res2_atom2_index_ = res2_atoms[2];
	res2_atom3_index_ = res2_atoms[3];
}

// Search constructor for reading out the motifs from the motif file for ligands
Motif::Motif(
	std::string const & resname1,
	std::string const & res1_atom1,
	std::string const & res1_atom2,
	std::string const & res1_atom3,
	std::string const & res2_atom1,
	std::string const & res2_atom2,
	std::string const & res2_atom3,
	core::kinematics::Jump const & orientation
):
	has_remark_( false ),
	has_path_( false )
{
	restype_name1_ = resname1;
	res1_atom1_name_ = res1_atom1;
	res1_atom2_name_ = res1_atom2;
	res1_atom3_name_ = res1_atom3;
	restype_name2_ = "LG1"; //We need to give the second residue (ligand residue) some random name
	res2_atom1_name_ = res2_atom1;
	res2_atom2_name_ = res2_atom2;
	res2_atom3_name_ = res2_atom3;

	forward_jump_ = orientation;
	backward_jump_ = forward_jump_.reversed();

	core::chemical::ResidueTypeSet & rsd_set( core::chemical::ChemicalManager::get_instance()->nonconst_residue_type_set( core::chemical::FA_STANDARD ) );
	core::chemical::ResidueType const & rsd_type1( rsd_set.name_map( restype_name1_ ) );

	res1_atom1_index_ = rsd_type1.atom_index( res1_atom1_name_ );
	res1_atom2_index_ = rsd_type1.atom_index( res1_atom2_name_ );
	res1_atom3_index_ = rsd_type1.atom_index( res1_atom3_name_ );
	// mt << "Res1: " <<  res1_atom1_name_ << res1_atom2_name_ << res1_atom3_name_ << " res2: " << res2_atom1_name_ << res2_atom2_name_ << res2_atom3_name_ << std::endl;
	//int res1_atom1_int_ = rsd_type.atom_type_index( res1_atom1_index_);
	//int res1_atom2_int_ = rsd_type.atom_type_index( res1_atom2_index_ );
	//int res1_atom3_int_ = rsd_type.atom_type_index( res1_atom3_index_ );
	core::chemical::AtomTypeSetCOP atset = core::chemical::ChemicalManager::get_instance()->atom_type_set( core::chemical::FA_STANDARD );
	res2_atom1_int_ = atset->atom_type_index(res2_atom1_name_);
	res2_atom2_int_ = atset->atom_type_index(res2_atom2_name_);
	res2_atom3_int_ = atset->atom_type_index(res2_atom3_name_);
	// mt << " res2: " << res2_atom1_int_ << ", " << res2_atom2_int_ << ", " << res2_atom3_int_ << std::endl;
}

Motif::Motif(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2
):
	has_remark_( false ),
	has_path_( false )
{
	if ( (res1.is_protein() && res2.is_protein() ) || (res1.is_protein() && res2.is_DNA() ) ) {
		restype_name1_ = protocols::dna::dna_full_name3( res1.name3() );
		res1_atom1_name_ = motifAtomIDs[restype_name1_][1];
		res1_atom2_name_ = motifAtomIDs[restype_name1_][2];
		res1_atom3_name_ = motifAtomIDs[restype_name1_][3];
		restype_name2_ = protocols::dna::dna_full_name3( res2.name3() );
		res2_atom1_name_ = motifAtomIDs[restype_name2_][1];
		res2_atom2_name_ = motifAtomIDs[restype_name2_][2];
		res2_atom3_name_ = motifAtomIDs[restype_name2_][3];
		forward_jump_ = core::kinematics::Jump(
			core::kinematics::Stub(  res1.atom( res1_atom2_name_ ).xyz(),
			res1.atom( res1_atom1_name_ ).xyz(),
			res1.atom( res1_atom2_name_ ).xyz(),
			res1.atom( res1_atom3_name_ ).xyz() ),
			core::kinematics::Stub(  res2.atom( res2_atom2_name_ ).xyz(),
			res2.atom( res2_atom1_name_ ).xyz(),
			res2.atom( res2_atom2_name_ ).xyz(),
			res2.atom( res2_atom3_name_ ).xyz() ) );
		backward_jump_ = forward_jump_.reversed();
	} else if ( res1.is_DNA() && res2.is_DNA() ) { // Can't use .build_motif_rotamers or .build_rotamers if both residues are DNA (at least right now)
		restype_name1_ = protocols::dna::dna_full_name3( res1.name3() );
		res1_atom1_name_ = basebaseAtomIDs[restype_name1_][1];
		res1_atom2_name_ = basebaseAtomIDs[restype_name1_][2];
		res1_atom3_name_ = basebaseAtomIDs[restype_name1_][3];
		restype_name2_ = protocols::dna::dna_full_name3( res2.name3() );
		res2_atom1_name_ = basebaseAtomIDs[restype_name2_][1];
		res2_atom2_name_ = basebaseAtomIDs[restype_name2_][2];
		res2_atom3_name_ = basebaseAtomIDs[restype_name2_][3];
		forward_jump_ = core::kinematics::Jump(
			core::kinematics::Stub(  res1.atom( res1_atom2_name_ ).xyz(),
			res1.atom( res1_atom1_name_ ).xyz(),
			res1.atom( res1_atom2_name_ ).xyz(),
			res1.atom( res1_atom3_name_ ).xyz() ),
			core::kinematics::Stub(  res2.atom( res2_atom2_name_ ).xyz(),
			res2.atom( res2_atom1_name_ ).xyz(),
			res2.atom( res2_atom2_name_ ).xyz(),
			res2.atom( res2_atom3_name_ ).xyz() ) );
		backward_jump_ = forward_jump_.reversed();
	} else {
		mt << "Input motif residues do not match the expected combinations for this constructor; please use constructor that allows you to explicitly specify atoms involved." << std::endl;
	}
	core::chemical::ResidueTypeSet & rsd_set( core::chemical::ChemicalManager::get_instance()->nonconst_residue_type_set( core::chemical::FA_STANDARD ) );
	core::chemical::ResidueType const & rsd_type( rsd_set.name_map( restype_name1_ ) );
	core::chemical::ResidueType const & rsd_type2( rsd_set.name_map( restype_name2_ ) );
	res1_atom1_index_ = rsd_type.atom_index( res1_atom1_name_ );
	res1_atom2_index_ = rsd_type.atom_index( res1_atom2_name_ );
	res1_atom3_index_ = rsd_type.atom_index( res1_atom3_name_ );
	res2_atom1_index_ = rsd_type2.atom_index( res2_atom1_name_ );
	res2_atom2_index_ = rsd_type2.atom_index( res2_atom2_name_ );
	res2_atom3_index_ = rsd_type2.atom_index( res2_atom3_name_ );
}

Motif::Motif(
	Motif const & src
) : utility::pointer::ReferenceCount(),
	restype_name1_( src.restype_name1_ ),
	res1_atom1_name_( src.res1_atom1_name_ ),
	res1_atom2_name_( src.res1_atom2_name_ ),
	res1_atom3_name_( src.res1_atom3_name_ ),
	res1_atom1_int_( src.res1_atom1_int_ ),
	res1_atom2_int_( src.res1_atom2_int_ ),
	res1_atom3_int_( src.res1_atom3_int_ ),
	res1_atom1_index_( src.res1_atom1_index_ ),
	res1_atom2_index_( src.res1_atom2_index_ ),
	res1_atom3_index_( src.res1_atom3_index_ ),
	restype_name2_( src.restype_name2_ ),
	res2_atom1_name_( src.res2_atom1_name_ ),
	res2_atom2_name_( src.res2_atom2_name_ ),
	res2_atom3_name_( src.res2_atom3_name_ ),
	res2_atom1_int_( src.res2_atom1_int_ ),
	res2_atom2_int_( src.res2_atom2_int_ ),
	res2_atom3_int_( src.res2_atom3_int_ ),
	res2_atom1_index_( src.res2_atom1_index_ ),
	res2_atom2_index_( src.res2_atom2_index_ ),
	res2_atom3_index_( src.res2_atom3_index_ ),
	forward_jump_( src.forward_jump_ ),
	backward_jump_( src.backward_jump_ ),
	remark_( src.remark_ ),
	path_( src.path_ ),
	has_remark_( src.has_remark_ ),
	has_path_( src.has_path_ )
{}

MotifOP
Motif::clone() const
{
	return MotifOP( new Motif( *this ) );
}

Motif::~Motif()
{}

bool
Motif::forward_check(
	core::conformation::Residue const & check_res
) const
{
	return( utility::trimmed_compare( check_res.name3(), restype_name1() ) );
}

bool
Motif::backward_check(
	core::conformation::Residue const & check_res
) const
{
	return( utility::trimmed_compare( check_res.name3(), restype_name2() ) );
}

bool
Motif::apply_check(
	core::pose::Pose const & pose,
	Size const pos
) const
{
	if ( pose.residue(pos).is_DNA() ) {
		return(
			utility::trimmed_compare(  protocols::dna::dna_full_name3(pose.residue( pos ).name3()), restype_name1_ ) ||
			utility::trimmed_compare(  protocols::dna::dna_full_name3(pose.residue( pos ).name3()), restype_name2_ )
		);
	}
	return( utility::trimmed_compare( pose.residue( pos ).name3(), restype_name1_ ) ||
		utility::trimmed_compare( pose.residue( pos ).name3(), restype_name2_ ) );
}
/*
// Search constructor for reading out the motifs from the motif file for ligands
Motif::Motif(
std::string const resname1,
std::string const res1_atom1,
std::string const res1_atom2,
std::string const res1_atom3,
std::string const res2_atom1,
std::string const res2_atom2,
std::string const res2_atom3,
core::kinematics::Jump const & orientation
)
{
restype_name1_ = resname1;
res1_atom1_name_ = res1_atom1;
res1_atom2_name_ = res1_atom2;
res1_atom3_name_ = res1_atom3;
restype_name2_ = "LG1"; //We need to give the second residue (ligand residue) some random name
res2_atom1_name_ = res2_atom1;
res2_atom2_name_ = res2_atom2;
res2_atom3_name_ = res2_atom3;
*/

/*
//Wasn't using this function, may want to put back in future.
void
Motif::generate_atom_ints(
)
{
core::chemical::AtomTypeSetCAP atset = core::chemical::ChemicalManager::get_instance()->atom_type_set( core::chemical::FA_STANDARD );
mt << "Res1: " <<  res1_atom1_name_ << res1_atom2_name_ << res1_atom3_name_ << " res2: " << res2_atom1_name_ << res2_atom2_name_ << res2_atom3_name_ << std::endl;

// core::chemical::ResidueTypeSet & rsd_set( core::chemical::ChemicalManager::get_instance()->nonconst_residue_type_set( core::chemical::FA_STANDARD ) ); // Unused variable causes warning
// core::chemical::ResidueType const & rsd_type( rsd_set.name_map( restype_name1_ ) ); // Unused variable causes warning
// Size res1_atom1_index = rsd_type.atom_index( res1_atom1_name_ ); // Unused variable causes warning
// Size res1_atom2_index = rsd_type.atom_index( res1_atom2_name_ ); // Unused variable causes warning
// Size res1_atom3_index = rsd_type.atom_index( res1_atom3_name_ ); // Unused variable causes warning

// FIXME: these local variables have the same name as data members!
// int res1_atom1_int_ = rsd_type.atom( res1_atom1_index).atom_type_index(); // Unused variable causes warning
// int res1_atom2_int_ = rsd_type.atom( res1_atom2_index ).atom_type_index(); // Unused variable causes warning
// int res1_atom3_int_ = rsd_type.atom( res1_atom3_index ).atom_type_index(); // Unused variable causes warning
int res2_atom1_int_ = atset->atom_type_index(res2_atom1_name_);
int res2_atom2_int_ = atset->atom_type_index(res2_atom2_name_);
int res2_atom3_int_ = atset->atom_type_index(res2_atom3_name_);
mt << " res2: " << res2_atom1_int_ << ", " << res2_atom2_int_ << ", " << res2_atom3_int_ << std::endl;
return;
}
*/

void
Motif::store_remark(
	std::string const & remark_in
)
{
	has_remark_ = true;
	remark_ = remark_in;
	return;
}

void
Motif::store_path(
	std::string const & path_in
)
{
	has_path_ = true;
	path_ = path_in;
	return;
}

core::pack::rotamer_set::RotamerSetOP
Motif::build_rotamers(
	core::pose::Pose & pose,
	Size const rotamer_build_position,
	Size const ex_,
	bool res2
) const
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
	task->set_bump_check( false );
	task->temporarily_fix_everything();
	task->temporarily_set_pack_residue( rotamer_build_position, true );
	utility::vector1< bool > aa_info( core::chemical::num_canonical_aas, false );

	// Here's the aa to build
	if ( res2 ) {
		aa_info[ aa_from_name( restype_name2() ) ] = true;
		//use_forward = true;  // set but never used ~Labonte
	} else {
		aa_info[ aa_from_name( restype_name1() ) ] = true;
		//use_forward = false;  // set but never used ~Labonte
	}
	// This functionality WILL fail if neither of the residues put into your constructor are amino acids
	// Should put a safety check here to make sure that only amino acids are being rotamerized

	task->nonconst_residue_task( rotamer_build_position ).restrict_absent_canonical_aas( aa_info );

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

	core::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, scorefxn, task );

	rotset->build_rotamers( pose, scorefxn, *task, packer_neighbor_graph, false );

	return rotset;
}

core::pack::rotamer_set::RotamerSetOP
Motif::build_inverted_rotamers(
	core::pose::Pose & pose,
	Size const motif_anchor_position,
	bool & use_forward,
	Size rotamer_build_position
) const
{
	using namespace core::pack::rotamer_set;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::pack;
	using namespace core::pack::task;

	Size extra_value( 0 );

	// If no build position was passed, it will be the default value ( 0 ),
	// and we need to find the closest bb position to the anchor residue
	if ( rotamer_build_position == Size( 0 ) ) {
		// Find closest
		// *********** TO DO ************
		mt << "Using nonsense build residue in build_motif_rotamer - default argument not yet supported" << std::endl;
		rotamer_build_position = 2;
	}

	RotamerSetOP rotset;

	// Check to make sure the motif applies to the anchor residue
	if ( !apply_check( pose, motif_anchor_position ) ) {
		mt << "Bailing from build_motif_rotamer - given anchor residue not in motif" << std::endl;
		return rotset;
	}

	use_forward = forward_check( pose.residue( motif_anchor_position ) );

	rotset = build_rotamers( pose, rotamer_build_position, extra_value, use_forward );

	// Invert the rotamer library as specified by the motif
	for ( Size ir = 1 , end_ir = rotset->num_rotamers() ; ir <= end_ir ; ++ir ) {
		place_residue( pose.residue( motif_anchor_position ), *(rotset->nonconst_rotamer( ir )) );
	}

	return rotset;
}

// place_atoms for ligands
void
Motif::place_atoms(
	core::conformation::Residue const & fixed,
	core::conformation::Residue & mobile,
	utility::vector1< Size > const & atoms,
	Size const & res2_atom1_index_in,
	Size const & res2_atom2_index_in,
	Size const & res2_atom3_index_in,
	bool one_three
) const
{
	std::string ligand_name( "LG1" );
	/* if( utility::trimmed_compare( protocols::dna::dna_full_name3( fixed.name3() ), restype_name1_ ) &&
	utility::trimmed_compare( protocols::dna::dna_full_name3( mobile.name3() ), restype_name2_ ) ) {
	return place_atoms_( fixed, mobile, true, atoms, one_three );

	} else if ( utility::trimmed_compare( protocols::dna::dna_full_name3( fixed.name3() ), restype_name2_ ) &&
	utility::trimmed_compare( protocols::dna::dna_full_name3( mobile.name3() ), restype_name1_ ) ) {
	return place_atoms_( fixed, mobile, false, atoms, one_three );

	} else*/ if ( utility::trimmed_compare( protocols::dna::dna_full_name3( fixed.name3() ), restype_name1_ ) &&
			utility::trimmed_compare( ligand_name , restype_name2_ ) ) {
		return place_atoms_( fixed, mobile, true, atoms, res2_atom1_index_in, res2_atom2_index_in, res2_atom3_index_in, one_three ); //This is Matt's new test for ligands
	} else {
		mt << "Bad Mojo a! call to Motif::place_atom() with wrong residue(s)!" << std::endl;
		mt << "Motif wants: " << restype_name1_ << " and " << restype_name2_ << std::endl;
		mt << "arguments are: " << fixed.name3() << " and " << mobile.name3() << std::endl;
		mt << "Neither order matches!" << std::endl;;
	}

	return;
}

// place_atom for ligands
void
Motif::place_atom(
	core::conformation::Residue const & fixed,
	core::conformation::Residue & mobile,
	core::conformation::Atom & atm,
	Size const & res2_atom1_index_in,
	Size const & res2_atom2_index_in,
	Size const & res2_atom3_index_in,
	//std::string const & atomtype,
	Size const & atomtype,
	bool one_three
) const
{
	return place_atom_( fixed, mobile, true, atm, res2_atom1_index_in, res2_atom2_index_in, res2_atom3_index_in, atomtype, one_three );
}

// place_residue for ligands
void
Motif::place_residue(
	core::conformation::Residue const & fixed,
	core::conformation::Residue & mobile,
	Size const & res2_atom1_index_in,
	Size const & res2_atom2_index_in,
	Size const & res2_atom3_index_in,
	bool one_three
) const
{
	//restype_name1_ is amino acid
	//restype_name2_ is ligand
	//fixed should be amino acid
	//mobile should be ligand (are there counterexamples?)
	//std::string ligand_name( "LG1" );
	if ( utility::trimmed_compare( protocols::dna::dna_full_name3( fixed.name3() ), restype_name1_ ) ) {
		return Motif::place_residue_( fixed, mobile, true,res2_atom1_index_in, res2_atom2_index_in, res2_atom3_index_in, one_three );

	} else if ( utility::trimmed_compare( protocols::dna::dna_full_name3( mobile.name3() ), restype_name1_ ) ) {
		return Motif::place_residue_( fixed, mobile, false, res2_atom1_index_in, res2_atom2_index_in, res2_atom3_index_in,  one_three );
	} else {
		mt << "Bad Mojo b! call to Motif::place_residue() with wrong residue(s)!" << std::endl;
		mt << "Motif wants: " << restype_name1_ << " and " << restype_name2_ << std::endl;
		mt << "arguments are: " << fixed.name3() << " and " << mobile.name3() << std::endl;
		mt << "Neither order matches!" << std::endl;
	}

	return;
}

//Matt's place_atom_ (this is where atoms are actually placed)
void
Motif::place_atom_(
	core::conformation::Residue const & fixed,
	core::conformation::Residue & mobile,
	bool forward,
	core::conformation::Atom & atm,
	Size const & res2_atom1_index_in,
	Size const & res2_atom2_index_in,
	Size const & res2_atom3_index_in,
	//std::string const & atomtype,
	Size const & atomtype,
	bool one_three
) const
{
	core::kinematics::Stub end_stub;
	core::kinematics::Stub mobile_stub;
	core::kinematics::Stub start_stub1;

	if ( forward ) {
		// Extract a stub from the fixed residue
		if ( one_three ) {
			core::kinematics::Stub start_stub(
				fixed.atom( res1_atom2_index_ ).xyz(),
				fixed.atom( res1_atom1_index_ ).xyz(),
				fixed.atom( res1_atom2_index_ ).xyz(),
				fixed.atom( res1_atom3_index_ ).xyz()
			);
			start_stub1 = start_stub;
		} else {
			core::kinematics::Stub start_stub(
				fixed.atom( res1_atom2_index_ ).xyz(),
				fixed.atom( res1_atom3_index_ ).xyz(),
				fixed.atom( res1_atom2_index_ ).xyz(),
				fixed.atom( res1_atom1_index_ ).xyz()
			);
			start_stub1 = start_stub;
		}

		mobile_stub.from_four_points(
			mobile.atom( res2_atom2_index_in ).xyz(),
			mobile.atom( res2_atom1_index_in ).xyz(),
			mobile.atom( res2_atom2_index_in ).xyz(),
			mobile.atom( res2_atom3_index_in ).xyz()
		);

		// Get the stub at the other end of the jump
		forward_jump_.make_jump( start_stub1, end_stub );
	} else {
		// Extract a stub from the fixed residue
		core::kinematics::Stub start_stub(
			fixed.atom( res2_atom2_index_in ).xyz(),
			fixed.atom( res2_atom1_index_in ).xyz(),
			fixed.atom( res2_atom2_index_in ).xyz(),
			fixed.atom( res2_atom3_index_in ).xyz()
		);

		if ( one_three ) {
			mobile_stub.from_four_points(
				mobile.atom( res1_atom2_index_ ).xyz(),
				mobile.atom( res1_atom1_index_ ).xyz(),
				mobile.atom( res1_atom2_index_ ).xyz(),
				mobile.atom( res1_atom3_index_ ).xyz()
			);
		} else {
			mobile_stub.from_four_points(
				mobile.atom( res1_atom2_index_ ).xyz(),
				mobile.atom( res1_atom3_index_ ).xyz(),
				mobile.atom( res1_atom2_index_ ).xyz(),
				mobile.atom( res1_atom1_index_ ).xyz()
			);
		}

		// Get the stub at the other end of the jump
		backward_jump_.make_jump( start_stub, end_stub );
	}

	// Apply to the mobile residue
	atm.xyz( end_stub.local2global( mobile_stub.global2local( mobile.xyz( atomtype ) ) ) );

	return;
}

void
Motif::place_atoms(
	core::conformation::Residue const & fixed,
	core::conformation::Residue & mobile,
	utility::vector1< std::string > const & atoms,
	bool one_three
) const
{
	if ( utility::trimmed_compare( protocols::dna::dna_full_name3( fixed.name3() ), restype_name1_ ) &&
			utility::trimmed_compare( protocols::dna::dna_full_name3( mobile.name3() ), restype_name2_ ) ) {
		return place_atoms_( fixed, mobile, true, atoms, one_three );

	} else if ( utility::trimmed_compare( protocols::dna::dna_full_name3( fixed.name3() ), restype_name2_ ) &&
			utility::trimmed_compare( protocols::dna::dna_full_name3( mobile.name3() ), restype_name1_ ) ) {
		return place_atoms_( fixed, mobile, false, atoms, one_three );

	} else {
		mt << "Bad Mojo! call to Motif::place_atom() with wrong residue(s)!" << std::endl;
		mt << "Motif wants: " << restype_name1_ << " and " << restype_name2_ << std::endl;
		mt << "arguments are: " << fixed.name3() << " and " << mobile.name3() << std::endl;
		mt << "Neither order matches!" << std::endl;
	}

	return;
}

void
Motif::place_atom(
	core::conformation::Residue const & fixed,
	core::conformation::Residue & mobile,
	core::conformation::Atom & atm,
	bool one_three,
	std::string const & atomtype
) const
{
	// PUT IN WARNING, ONLY CALL THIS FUNCTION WITH TRUE
	// IE, WITH ATM BEING AN ATOM FROM THE MOBILE RESIDUE, NOT THE FIXED
	return place_atom_( fixed, mobile, true, atm, one_three, atomtype );
	/*if( utility::trimmed_compare( protocols::dna::dna_full_name3( fixed.name3() ), restype_name1_ ) &&
	utility::trimmed_compare( protocols::dna::dna_full_name3( mobile.name3() ), restype_name2_ ) ) {
	return place_atom_( fixed, mobile, true, atm, one_three, atomtype );

	} else if ( utility::trimmed_compare( protocols::dna::dna_full_name3( fixed.name3() ), restype_name2_ ) &&
	utility::trimmed_compare( protocols::dna::dna_full_name3( mobile.name3() ), restype_name1_ ) ) {
	return place_atom_( fixed, mobile, false, atm, one_three, atomtype );

	} else {
	mt << "Bad Mojo! call to Motif::place_atom() with wrong residue(s)!\n";
	mt << "Motif wants: " << restype_name1_ << " and " << restype_name2_ << "\n";
	mt << "arguments are: " << fixed.name3() << " and " << mobile.name3() << "\n";
	mt << "Neither order matches!\n";
	}

	return;*/
}

void
Motif::place_residue(
	core::conformation::Residue const & fixed,
	core::conformation::Residue & mobile,
	bool one_three
) const
{
	if ( utility::trimmed_compare( protocols::dna::dna_full_name3( fixed.name3() ), restype_name1_ ) &&
			utility::trimmed_compare( protocols::dna::dna_full_name3( mobile.name3() ), restype_name2_ ) ) {
		return place_residue_( fixed, mobile, true, one_three );

	} else if ( utility::trimmed_compare( protocols::dna::dna_full_name3( fixed.name3() ), restype_name2_ ) &&
			utility::trimmed_compare( protocols::dna::dna_full_name3( mobile.name3() ), restype_name1_ ) ) {
		return place_residue_( fixed, mobile, false, one_three );

	} else {
		mt << "Bad Mojo! call to Motif::place_residue() with wrong residue(s)!" << std::endl;
		mt << "Motif wants: " << restype_name1_ << " and " << restype_name2_ << std::endl;
		mt << "arguments are: " << fixed.name3() << " and " << mobile.name3() << std::endl;
		mt << "Neither order matches!" << std::endl;
	}

	return;
}

void
Motif::place_residue_helper(
	core::kinematics::Stub & end_stub,
	core::kinematics::Stub & mobile_stub,
	core::conformation::Residue const & fixed,
	core::conformation::Residue & mobile,
	bool forward,
	bool one_three
) const {

	core::kinematics::Stub start_stub1;

	if ( forward ) {
		// Extract a stub from the fixed residue
		if ( one_three ) {
			core::kinematics::Stub start_stub(
				fixed.atom( res1_atom2_name_ ).xyz(),
				fixed.atom( res1_atom1_name_ ).xyz(),
				fixed.atom( res1_atom2_name_ ).xyz(),
				fixed.atom( res1_atom3_name_ ).xyz()
			);
			start_stub1 = start_stub;
		} else {
			core::kinematics::Stub start_stub(
				fixed.atom( res1_atom2_name_ ).xyz(),
				fixed.atom( res1_atom3_name_ ).xyz(),
				fixed.atom( res1_atom2_name_ ).xyz(),
				fixed.atom( res1_atom1_name_ ).xyz()
			);
			start_stub1 = start_stub;
		}

		mobile_stub.from_four_points(
			mobile.atom( res2_atom2_name_ ).xyz(),
			mobile.atom( res2_atom1_name_ ).xyz(),
			mobile.atom( res2_atom2_name_ ).xyz(),
			mobile.atom( res2_atom3_name_ ).xyz()
		);

		// Get the stub at the other end of the jump
		forward_jump_.make_jump( start_stub1, end_stub );

	} else {
		// Extract a stub from the fixed residue
		core::kinematics::Stub start_stub(
			fixed.atom( res2_atom2_name_ ).xyz(),
			fixed.atom( res2_atom1_name_ ).xyz(),
			fixed.atom( res2_atom2_name_ ).xyz(),
			fixed.atom( res2_atom3_name_ ).xyz()
		);

		// This one_three bool assumes that res1 is always the amino acid of a base-protein interaction
		// The base will never need automorphism
		if ( one_three ) {
			mobile_stub.from_four_points(
				mobile.atom( res1_atom2_name_ ).xyz(),
				mobile.atom( res1_atom1_name_ ).xyz(),
				mobile.atom( res1_atom2_name_ ).xyz(),
				mobile.atom( res1_atom3_name_ ).xyz()
			);
		} else {
			mobile_stub.from_four_points(
				mobile.atom( res1_atom2_name_ ).xyz(),
				mobile.atom( res1_atom3_name_ ).xyz(),
				mobile.atom( res1_atom2_name_ ).xyz(),
				mobile.atom( res1_atom1_name_ ).xyz()
			);
		}

		// Get the stub at the other end of the jump
		backward_jump_.make_jump( start_stub, end_stub );
	}
}

void
Motif::place_residue_(
	core::conformation::Residue const & fixed,
	core::conformation::Residue & mobile,
	bool forward,
	bool one_three
) const {
	core::kinematics::Stub end_stub;
	core::kinematics::Stub mobile_stub;
	//core::kinematics::Stub start_stub1;

	place_residue_helper( end_stub, mobile_stub, fixed, mobile, forward, one_three );

	// Apply to the mobile residue
	for ( Size i(1), end_i = mobile.natoms() ; i <= end_i ; ++i ) {
		mobile.set_xyz( i, end_stub.local2global( mobile_stub.global2local( mobile.xyz( i ) ) ) );
	}

	return;
}

void
Motif::place_residue_helper(
	core::kinematics::Stub & end_stub,
	core::kinematics::Stub & mobile_stub,
	core::conformation::Residue const & fixed,
	core::conformation::Residue & mobile,
	bool forward,
	Size const & res2_atom1_index_in,
	Size const & res2_atom2_index_in,
	Size const & res2_atom3_index_in,
	bool one_three
) const {
	core::kinematics::Stub start_stub1;
	if ( forward ) {
		// Extract a stub from the fixed residue
		if ( one_three ) {
			core::kinematics::Stub start_stub(
				fixed.atom( res1_atom2_name_ ).xyz(),
				fixed.atom( res1_atom1_name_ ).xyz(),
				fixed.atom( res1_atom2_name_ ).xyz(),
				fixed.atom( res1_atom3_name_ ).xyz()
			);
			start_stub1 = start_stub;
		} else {
			core::kinematics::Stub start_stub(
				fixed.atom( res1_atom2_name_ ).xyz(),
				fixed.atom( res1_atom3_name_ ).xyz(),
				fixed.atom( res1_atom2_name_ ).xyz(),
				fixed.atom( res1_atom1_name_ ).xyz()
			);
			start_stub1 = start_stub;
		}

		mobile_stub.from_four_points(
			mobile.atom( res2_atom2_index_in ).xyz(),
			mobile.atom( res2_atom1_index_in ).xyz(),
			mobile.atom( res2_atom2_index_in ).xyz(),
			mobile.atom( res2_atom3_index_in ).xyz()
		);
		// Get the stub at the other end of the jump
		forward_jump_.make_jump( start_stub1, end_stub );
	} else {
		// Extract a stub from the fixed residue
		core::kinematics::Stub start_stub(
			fixed.atom( res2_atom2_index_in ).xyz(),
			fixed.atom( res2_atom1_index_in ).xyz(),
			fixed.atom( res2_atom2_index_in ).xyz(),
			fixed.atom( res2_atom3_index_in ).xyz()
		);

		// This one_three bool assumes that res1 is always the amino acid of a base-protein interaction
		// The base will never need automorphism
		if ( one_three ) {
			mobile_stub.from_four_points(
				mobile.atom( res1_atom2_name_ ).xyz(),
				mobile.atom( res1_atom1_name_ ).xyz(),
				mobile.atom( res1_atom2_name_ ).xyz(),
				mobile.atom( res1_atom3_name_ ).xyz()
			);
		} else {
			mobile_stub.from_four_points(
				mobile.atom( res1_atom2_name_ ).xyz(),
				mobile.atom( res1_atom3_name_ ).xyz(),
				mobile.atom( res1_atom2_name_ ).xyz(),
				mobile.atom( res1_atom1_name_ ).xyz()
			);
		}
		// Get the stub at the other end of the jump
		backward_jump_.make_jump( start_stub, end_stub );
	}
}

// place_residue_ for ligand motifs
void
Motif::place_residue_(
	core::conformation::Residue const & fixed,
	core::conformation::Residue & mobile,
	bool forward,
	Size const & res2_atom1_index_in,
	Size const & res2_atom2_index_in,
	Size const & res2_atom3_index_in,
	bool one_three
) const
{
	core::kinematics::Stub end_stub;
	core::kinematics::Stub mobile_stub;
	// core::kinematics::Stub start_stub1;

	place_residue_helper( end_stub, mobile_stub, fixed, mobile, forward,
		res2_atom1_index_in, res2_atom2_index_in, res2_atom3_index_in, one_three );

	// Apply to the mobile residue
	for ( Size i(1), end_i = mobile.natoms() ; i <= end_i ; ++i ) {
		mobile.set_xyz( i, end_stub.local2global( mobile_stub.global2local( mobile.xyz( i ) ) ) );
	}
	return;
}

void
Motif::place_atom_(
	core::conformation::Residue const & fixed,
	core::conformation::Residue & mobile,
	bool forward,
	core::conformation::Atom & atm,
	bool one_three,
	std::string const & atomtype
) const
{
	core::kinematics::Stub end_stub;
	core::kinematics::Stub mobile_stub;
	core::kinematics::Stub start_stub1;

	if ( forward ) {
		// Extract a stub from the fixed residue
		if ( one_three ) {
			core::kinematics::Stub start_stub(
				fixed.atom( res1_atom2_index_ ).xyz(),
				fixed.atom( res1_atom1_index_ ).xyz(),
				fixed.atom( res1_atom2_index_ ).xyz(),
				fixed.atom( res1_atom3_index_ ).xyz()
			);
			start_stub1 = start_stub;
		} else {
			core::kinematics::Stub start_stub(
				fixed.atom( res1_atom2_index_ ).xyz(),
				fixed.atom( res1_atom3_index_ ).xyz(),
				fixed.atom( res1_atom2_index_ ).xyz(),
				fixed.atom( res1_atom1_index_ ).xyz()
			);
			start_stub1 = start_stub;
		}

		mobile_stub.from_four_points(
			mobile.atom( res2_atom2_index_ ).xyz(),
			mobile.atom( res2_atom1_index_ ).xyz(),
			mobile.atom( res2_atom2_index_ ).xyz(),
			mobile.atom( res2_atom3_index_ ).xyz()
		);

		// Get the stub at the other end of the jump
		forward_jump_.make_jump( start_stub1, end_stub );
	} else {
		// Extract a stub from the fixed residue
		core::kinematics::Stub start_stub(
			fixed.atom( res2_atom2_index_ ).xyz(),
			fixed.atom( res2_atom1_index_ ).xyz(),
			fixed.atom( res2_atom2_index_ ).xyz(),
			fixed.atom( res2_atom3_index_ ).xyz()
		);

		if ( one_three ) {
			mobile_stub.from_four_points(
				mobile.atom( res1_atom2_index_ ).xyz(),
				mobile.atom( res1_atom1_index_ ).xyz(),
				mobile.atom( res1_atom2_index_ ).xyz(),
				mobile.atom( res1_atom3_index_ ).xyz()
			);
		} else {
			mobile_stub.from_four_points(
				mobile.atom( res1_atom2_index_ ).xyz(),
				mobile.atom( res1_atom3_index_ ).xyz(),
				mobile.atom( res1_atom2_index_ ).xyz(),
				mobile.atom( res1_atom1_index_ ).xyz()
			);
		}

		// Get the stub at the other end of the jump
		backward_jump_.make_jump( start_stub, end_stub );
	}

	// Apply to the mobile residue
	atm.xyz( end_stub.local2global( mobile_stub.global2local( mobile.xyz( atomtype ) ) ) );

	return;
}

void
Motif::place_atoms_(
	core::conformation::Residue const & fixed,
	core::conformation::Residue & mobile,
	bool forward,
	utility::vector1< std::string > const & atoms,
	bool one_three
) const
{
	core::kinematics::Stub end_stub;
	core::kinematics::Stub mobile_stub;
	//core::kinematics::Stub start_stub1;

	place_residue_helper( end_stub, mobile_stub, fixed, mobile, forward, one_three );

	// Apply to the mobile residue
	for ( Size i(1), end_i = atoms.size() ; i <= end_i ; ++i ) {
		mobile.set_xyz( atoms[i], end_stub.local2global( mobile_stub.global2local( mobile.xyz( atoms[i] ) ) ) );
	}

	return;
}

// place_atoms_ for ligands
void
Motif::place_atoms_(
	core::conformation::Residue const & fixed,
	core::conformation::Residue & mobile,
	bool forward,
	utility::vector1< Size > const & atoms,
	Size const & res2_atom1_index_in,
	Size const & res2_atom2_index_in,
	Size const & res2_atom3_index_in,
	bool one_three
) const
{
	core::kinematics::Stub end_stub;
	core::kinematics::Stub mobile_stub;
	//core::kinematics::Stub start_stub1;

	place_residue_helper( end_stub, mobile_stub, fixed, mobile, forward,
		res2_atom1_index_in, res2_atom2_index_in, res2_atom3_index_in, one_three );

	// Apply to the mobile residue
	for ( Size i(1), end_i = atoms.size() ; i <= end_i ; ++i ) {
		mobile.set_xyz( atoms[i], end_stub.local2global( mobile_stub.global2local( mobile.xyz( atoms[i] ) ) ) );
	}

	return;
}

void
Motif::print( std::ostream & out ) const
{

	out << "SINGLE   " << restype_name1() << "   " << res1_atom1_name()
		<< "   " << res1_atom2_name()
		<< "   " << res1_atom3_name()
		<< "   " << restype_name2() <<  "   " << res2_atom1_name()
		<< "   " << res2_atom2_name()
		<< "   " << res2_atom3_name();// << '\n';

	if ( has_remark() ) {
		out << " ### REMARK " << remark();
	}
	out << '\n';

	out << forward_jump() << '\n';

}

std::string
Motif::print() const
{
	std::ostringstream os;
	print( os );
	return os.str();
}

std::ostream & operator <<(
	std::ostream & os,
	Motif const & motif
)
{
	os << motif.print();
	return os;
}

// Map of the atoms that are used in each motif, residues can be added to this file or modified if you want to use a different set of atoms (or you can use the alternative constructor that allows you set up the motif explicitly
std::map < std::string, utility::vector1< std::string > > Motif::motifAtomIDs(
	utility::tools::make_map(
	std::string("ADE"), utility::tools::make_vector1( std::string("N6"), std::string("C5"), std::string("N7") ),
	std::string("CYT"), utility::tools::make_vector1( std::string("N4"), std::string("C4"), std::string("C5") ),
	std::string("GUA"), utility::tools::make_vector1( std::string("O6"), std::string("C5"), std::string("N7") ),
	std::string("THY"), utility::tools::make_vector1( std::string("O4"), std::string("C5"), std::string("C7") ),
	std::string("ALA"), utility::tools::make_vector1( std::string("CB"), std::string("CA"), std::string("N") ),
	std::string("CYS"), utility::tools::make_vector1( std::string("SG"), std::string("CB"), std::string("CA") ),
	std::string("ASP"), utility::tools::make_vector1( std::string("OD1"), std::string("CG"), std::string("OD2") ),
	std::string("GLU"), utility::tools::make_vector1( std::string("OE1"), std::string("CD"), std::string("OE2") ),
	std::string("PHE"), utility::tools::make_vector1( std::string("CE1"), std::string("CZ"), std::string("CE2") ),
	std::string("HIS"), utility::tools::make_vector1( std::string("ND1"), std::string("CE1"), std::string("NE2") ),
	std::string("ILE"), utility::tools::make_vector1( std::string("CG1"), std::string("CB"), std::string("CG2") ),
	std::string("LYS"), utility::tools::make_vector1( std::string("NZ"), std::string("CE"), std::string("CD") ),
	std::string("LEU"), utility::tools::make_vector1( std::string("CD1"), std::string("CG"), std::string("CD2") ),
	std::string("MET"), utility::tools::make_vector1( std::string("CG"), std::string("SD"), std::string("CE") ),
	std::string("ASN"), utility::tools::make_vector1( std::string("OD1"), std::string("CG"), std::string("ND2") ),
	std::string("GLN"), utility::tools::make_vector1( std::string("OE1"), std::string("CD"), std::string("NE2") ),
	std::string("PRO"), utility::tools::make_vector1( std::string("CG"), std::string("CD"), std::string("NV") ),
	std::string("ARG"), utility::tools::make_vector1( std::string("NH1"), std::string("CZ"), std::string("NH2") ),
	std::string("SER"), utility::tools::make_vector1( std::string("OG"), std::string("CB"), std::string("CA") ),
	std::string("THR"), utility::tools::make_vector1( std::string("OG1"), std::string("CB"), std::string("CG2") ),
	std::string("VAL"), utility::tools::make_vector1( std::string("CG1"), std::string("CB"), std::string("CG2") ),
	std::string("TRP"), utility::tools::make_vector1( std::string("CG"), std::string("CD2"), std::string("CE2") ),
	std::string("TYR"), utility::tools::make_vector1( std::string("CE1"), std::string("CZ"), std::string("CE2") )
	) );

// Map of atoms used for base-base motifs
// Do not try to use rotamer building functions if the motifs are basebase
// This map isn't really used at this point
std::map < std::string, utility::vector1< std::string > > Motif::basebaseAtomIDs(
	utility::tools::make_map(
	std::string("ADE"), utility::tools::make_vector1( std::string("N6"), std::string("N1"), std::string("C2") ),
	std::string("CYT"), utility::tools::make_vector1( std::string("N4"), std::string("N3"), std::string("O2") ),
	std::string("GUA"), utility::tools::make_vector1( std::string("O6"), std::string("N1"), std::string("N2") ),
	std::string("THY"), utility::tools::make_vector1( std::string("O4"), std::string("N3"), std::string("O2") )
	) );

} // namespace motifs
} // namespace protocols
