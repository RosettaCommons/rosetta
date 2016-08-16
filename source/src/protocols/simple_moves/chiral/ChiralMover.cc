// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/chiral/ChiralMover.cc
/// @brief ChiralMover methods implemented
/// @author Kevin Drew, kdrew@nyu.edu

// Unit Headers
#include <protocols/simple_moves/chiral/ChiralMover.hh>
// Package Headers

// Project Headers
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Patch.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.functions.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
// Utility Headers
#include <numeric/xyz.functions.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
//#include <core/types.hh>

// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.chiral.ChiralMover" );

using namespace core;
using namespace conformation;
using namespace chemical;
using namespace core::id;

namespace protocols {
namespace simple_moves {
namespace chiral {

ResidueType const & get_chiral_residue_type( ResidueType const & rt, Chirality chirality )
{
	chemical::ResidueTypeSetCOP residue_type_set = rt.residue_type_set();
	//kdrew: first letters of a residuetype name (before '_p') are the letter code for the aa and is what is stored in the map
	//std::string base_name;
	//std::string patch_name;

	//kdrew: use for finding base_name instead of substr
	std::string const base_name( residue_type_base_name( rt ) );
	//std::string const patch_name( residue_type_all_patches_name( rt ) );
	TR << "base_name: " << base_name /*<< " patch_name: " << patch_name */<< std::endl;

	//kdrew: is residuetype patched?
	//Size base_end_pos = rt.name().find("_p");
	//TR << "base_end_pos: " << base_end_pos << std::endl;
	//if( base_end_pos != std::string::npos )
	//{
	// base_name = rt.name().substr( 0, base_end_pos );
	// //kdrew: the remaining porition of the string is the patch, reapply this porition to the end of the mapped string
	// patch_name = rt.name().substr( base_end_pos, rt.name().size() );
	//}
	////kdrew: if not patched just use the residuetype name
	//else
	//{
	// base_name = rt.name();
	// patch_name = "";
	//}

	TR << "restype: " << rt.name() << " " << rt.aa() << std::endl;
	//std::string chiral_name = rt.chiral_equivalent_name();
	utility::vector1< std::string > variant_types = rt.properties().get_list_of_variants();

	if ( !rt.is_l_aa() && !rt.is_d_aa() ) {
		TR << " possibly achiral (ex GLY)" <<  std::endl;
		return rt;
	}

	// Prepend or remove D depending on targeting.
	if ( chirality == L_CHIRALITY && rt.is_d_aa() ) {
		return residue_type_set->name_map( rt.name().substr( 1 ) );
	} else if ( chirality == D_CHIRALITY && rt.is_l_aa() ) {
		return residue_type_set->name_map( "D"+rt.name() );
	}

	// if all else fails, return original rt
	return rt;
}

/*
kdrew: the apply function changes a single residue's chirality
pose: pose to make change to
chiral_seq_pos: position of residue to change chirality
*/
void ChiralMover::apply( core::pose::Pose & pose )
{
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	//kdrew: assert for validity of parameters
	//runtime_assert ( chiral_seq_pos_ != 1 );

	Real phi_angle = pose.phi( chiral_seq_pos_ );
	Real psi_angle = pose.psi( chiral_seq_pos_ );


	TR << "phi_angle: " << phi_angle << " psi_angle: " << psi_angle << std::endl;

	//kdrew: get residue type
	ResidueType rtype = pose.residue_type( chiral_seq_pos_ );
	TR << "Current residue type: " << rtype.name()  << std::endl;
	TR << "Current residue type lower terminus: " << rtype.is_lower_terminus()  << std::endl;
	TR << "Current residue type upper terminus: " << rtype.is_upper_terminus()  << std::endl;

	AtomID const atom1( rtype.atom_index( "N" ), chiral_seq_pos_ );
	AtomID const atom2( rtype.atom_index( "CA" ), chiral_seq_pos_ );
	AtomID const atom3( rtype.atom_index( "C" ), chiral_seq_pos_ );
	AtomID const atom4( rtype.atom_index( "O" ), chiral_seq_pos_ );

	Real pseudo_psi = 0.0;
	if ( rtype.is_upper_terminus() ) {
		pseudo_psi = pose.conformation().torsion_angle( atom1, atom2, atom3, atom4);
		TR << "pseudo_psi: " << pseudo_psi << std::endl;
	}
	//kdrew: get chiral residue type
	ResidueType const & chiral_rtype = get_chiral_residue_type( rtype, chirality_ );
	TR << "Flipped residue type: " << chiral_rtype.name()  << " " << chiral_rtype.aa() << std::endl;
	if ( chiral_rtype.name() == rtype.name() ) {
		TR << " not making chiral change" << std::endl;
		return;
	}
	conformation::Residue res = pose.residue(chiral_seq_pos_);
	//kdrew: mutate to chiral residue type
	conformation::Residue replace_res ( chiral_rtype, true );

	if ( orient_functional_group_ ) {
		AtomIndices rtype_sidechain_atoms = rtype.all_sc_atoms();

		AtomIndices chiral_neighbor_ids = chiral_rtype.bonded_neighbor(chiral_rtype.first_sidechain_atom());
		AtomIndices neighbor_ids = rtype.bonded_neighbor(rtype.first_sidechain_atom());

		utility::vector1< std::pair< std::string, std::string > > atom_pairs;

		if ( !rtype.atom_type(rtype.first_sidechain_atom()).is_hydrogen() && !chiral_rtype.atom_type(chiral_rtype.first_sidechain_atom()).is_hydrogen()  ) {
			TR << rtype.name() << " " << rtype.atom_name(rtype.first_sidechain_atom()) << std::endl;
			atom_pairs.push_back( std::make_pair( rtype.atom_name(rtype.first_sidechain_atom()), chiral_rtype.atom_name(chiral_rtype.first_sidechain_atom() ) ));
		} else {
			TR << "First sidechain atoms are not heavy atoms " << std::endl;
		}

		//kdrew: loop through all neighboring atoms to first sidechain atom looking for the first two heavy atoms to use for alignment
		for ( core::Size j = 1; j <= neighbor_ids.size() && atom_pairs.size() < 3; ++j ) {
			//kdrew: if heavy atom
			if ( !rtype.atom_type(neighbor_ids[j]).is_hydrogen() ) {
				TR << rtype.name() << " " << rtype.atom_name(neighbor_ids[j]) << std::endl;
				//TR << chiral_rtype.name() << " " << chiral_rtype.atom_name(chiral_neighbor_ids[j]) << std::endl;

				//kdrew: assumes atom names are the same in both L and D residue types
				atom_pairs.push_back( std::make_pair( rtype.atom_name( neighbor_ids[j] ), rtype.atom_name( neighbor_ids[j] )));
			}
		}
		//if( !chiral_rtype.atom_type(chiral_neighbor_ids[j]).is_hydrogen() )

		//kdrew: if there are not enough side chain atoms (ex. ALA) choose an atom from the backbone (ex. N) for alignment
		if ( atom_pairs.size() < 3 ) {
			for ( core::Size j = 1; j <= neighbor_ids.size() && atom_pairs.size() < 3; ++j ) {
				TR << "looping through neighbors: " << rtype.name() << " " << rtype.atom_name(neighbor_ids[j]) << std::endl;

				//kdrew: if heavy atom
				if ( !rtype.atom_type(neighbor_ids[j]).is_hydrogen() ) {
					AtomIndices neighbor_neighbor_ids = rtype.bonded_neighbor( neighbor_ids[j]  );
					for ( core::Size jj = 1; jj <= neighbor_neighbor_ids.size() && atom_pairs.size() < 3; ++jj ) {
						if ( !rtype.atom_type(neighbor_neighbor_ids[jj]).is_hydrogen() ) {
							TR << "looping through neighbor neighbors: " << rtype.name() << " " << rtype.atom_name(neighbor_neighbor_ids[jj]) << std::endl;
							atom_pairs.push_back( std::make_pair( rtype.atom_name( neighbor_neighbor_ids[jj] ), rtype.atom_name( neighbor_neighbor_ids[jj] )));
						}
					}

				}
			}

		}

		//kdrew: assert that residue type has three side chain atoms
		runtime_assert_msg ( atom_pairs.size() == 3 , "not enough heavy atoms to align residues" );

		core::conformation::idealize_hydrogens( replace_res, pose.conformation() );
		replace_res.orient_onto_residue( res, atom_pairs );

		//kdrew: for all sidechain atoms in residue_type, copy xyz
		for ( core::Size j = 1; j <= rtype_sidechain_atoms.size(); ++j ) {
			replace_res.atom(rtype_sidechain_atoms[j]).xyz( res.atom(rtype_sidechain_atoms[j]).xyz() );
		}

		pose.replace_residue( chiral_seq_pos_ , replace_res, false );
	} else {

		pose.replace_residue( chiral_seq_pos_ , replace_res, true );
		//pose.dump_pdb( "rosetta_out_chiral_preidealized.pdb" );

		//kdrew: idealize alpha carbon hydrogen
		core::conformation::ResidueOP iires = pose.residue( chiral_seq_pos_ ).clone();
		core::conformation::idealize_hydrogens( *iires, pose.conformation() );
		pose.replace_residue( chiral_seq_pos_, *iires, false );
	}

	//pose.dump_pdb( "rosetta_out_chiral_postidealized.pdb" );
	pose.set_phi( chiral_seq_pos_, (-1.0 * phi_angle ) );
	//pose.dump_pdb( "rosetta_out_chiral_moved_phi.pdb" );
	pose.set_psi( chiral_seq_pos_, (-1.0 * psi_angle ) );
	//pose.dump_pdb( "rosetta_out_chiral_moved_phi_psi.pdb" );

	AtomID const atom1c( chiral_rtype.atom_index( "N" ), chiral_seq_pos_ );
	AtomID const atom2c( chiral_rtype.atom_index( "CA" ), chiral_seq_pos_ );
	AtomID const atom3c( chiral_rtype.atom_index( "C" ), chiral_seq_pos_ );
	AtomID const atom4c( chiral_rtype.atom_index( "O" ), chiral_seq_pos_ );

	if ( chiral_rtype.is_upper_terminus() ) {
		Real chiral_pseudo_psi = pose.conformation().torsion_angle( atom1c, atom2c, atom3c, atom4c);
		TR << "chiral pseudo_psi: " << chiral_pseudo_psi << std::endl;
		pose.conformation().set_torsion_angle( atom1c, atom2c, atom3c, atom4c, -1.0*pseudo_psi);
	}

	TR << "chiral phi_angle: " << pose.phi( chiral_seq_pos_ ) << " chiral psi_angle: " << pose.psi( chiral_seq_pos_ ) << std::endl;

	TR<< "exiting apply" << std::endl;
}

std::string
ChiralMover::get_name() const {
	return "ChiralMover";
}

/// @brief
ChiralMover::ChiralMover(
	core::Size chiral_seq_position
): Mover(), chiral_seq_pos_( chiral_seq_position ), chirality_( FLIP_CHIRALITY ), orient_functional_group_(false)
{
	Mover::type( "ChiralMover" );
}

ChiralMover::ChiralMover(
	core::Size chiral_seq_position,
	Chirality chirality
): Mover(), chiral_seq_pos_( chiral_seq_position ), chirality_( chirality ), orient_functional_group_(false)
{
	Mover::type( "ChiralMover" );
}

ChiralMover::ChiralMover(
	core::Size chiral_seq_position,
	bool orient_functional_group
): Mover(), chiral_seq_pos_( chiral_seq_position ), chirality_( FLIP_CHIRALITY ), orient_functional_group_(orient_functional_group)
{
	Mover::type( "ChiralMover" );
}

ChiralMover::ChiralMover(
	core::Size chiral_seq_position,
	Chirality chirality,
	bool orient_functional_group
): Mover(), chiral_seq_pos_( chiral_seq_position ), chirality_( chirality ), orient_functional_group_(orient_functional_group)
{
	Mover::type( "ChiralMover" );
}

ChiralMover::~ChiralMover(){}

}//chiral
}//simple_moves
}//protocols

