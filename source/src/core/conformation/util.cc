// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/util.cc
/// @author Phil Bradley


// Unit headers
#include <core/conformation/util.hh>

// Package headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/Stub.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/IdealBondLengthSet.hh>
#include <core/chemical/rings/util.hh>
#include <core/chemical/rna/util.hh> // for default root atom -- there is a choice encoded below for DNA vs. RNA vs. proteins
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/constants.hh>
#include <core/kinematics/util.hh>

// Numeric headers
#include <numeric/util.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/assert.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <algorithm> //std::min()

// Additional misc headers:
#include <utility/assert.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer TR( "core.conformation.util" );


namespace core {
namespace conformation {

void
orient_residue_for_ideal_bond(
	Residue & moving_rsd,
	chemical::ResidueConnection const & moving_connection,
	Residue const & fixed_rsd,
	chemical::ResidueConnection const & fixed_connection,
	Conformation const & conformation,
	bool lookup_bond_length
)
{

	// confirm that both connections are defined entirely by atoms internal to their residues,
	// so we don't need the conformation info (which might not be correct for the moving_rsd if we
	// are in the process of adding it to the conformation)
	//
	runtime_assert( moving_connection.icoor().is_internal() && fixed_connection.icoor().is_internal() );

	// a1 ---  a2  --- (a3)
	//
	//  (b2) ---  b3  --- b4
	//
	// a3 and b2 are the pseudo-atoms built from their respective residues using the ideal geometry in the
	// corresponding ResidueConnection objects
	//
	Vector
		a1(  fixed_connection.icoor().stub_atom2().xyz(  fixed_rsd, conformation ) ),
		b4( moving_connection.icoor().stub_atom2().xyz( moving_rsd, conformation ) ),

		a2(  fixed_rsd.xyz(  fixed_connection.atomno() ) ),
		b3( moving_rsd.xyz( moving_connection.atomno() ) ),

		a3(  fixed_connection.icoor().build(  fixed_rsd, conformation ) ),
		b2( moving_connection.icoor().build( moving_rsd, conformation ) );

	// we want to move b2 to align with a2 and b3 to align with a3. Torsion about the a2->a3 bond
	// (ie the inter-residue bond) determined by atoms a1 and b4 (torsion set to 0 by default).

	core::Size fixed_rsd_atom_type_index = fixed_rsd.atom_type_index(fixed_connection.atomno());
	core::Size moving_rsd_atom_type_index = moving_rsd.atom_type_index(moving_connection.atomno());
	if ( lookup_bond_length ) {
		if ( TR.visible() ) {
			TR << "moving atom_type_index " << moving_rsd_atom_type_index<< std::endl;
			TR << "fixed atom_type_index " << fixed_rsd_atom_type_index<< std::endl;
		}

		// Real bond_length= lookup_bond_length(fixed_rsd_atom_type_index, moving_rsd_atom_type_index);
		core::chemical::ChemicalManager *cm= core::chemical::ChemicalManager::get_instance();
		core::chemical::IdealBondLengthSetCOP ideal_bond_lengths( cm->ideal_bond_length_set( core::chemical::FA_STANDARD ) );

		Real bond_length= ideal_bond_lengths->get_bond_length( fixed_rsd_atom_type_index, moving_rsd_atom_type_index);
		Vector old_bond= a3-a2;
		Vector new_bond= old_bond;
		// new_bond.normalize(2);
		new_bond.normalize(bond_length);
		Vector offset = new_bond - old_bond;


		a3 += offset;
		//b2 -= offset;
	}


	kinematics::Stub A( a3, a2, a1 ), B( b3, b2, b4 );

	for ( Size i=1; i<= moving_rsd.natoms(); ++i ) {
		Vector global = A.local2global( B.global2local( moving_rsd.xyz(i)));
		//global.z(global.z() -1.1);
		moving_rsd.set_xyz(i, global );
	}

	Vector global_action_coord= A.local2global( B.global2local( moving_rsd.actcoord() ));
	//global_action_coord.z(global_action_coord.z() -1.1);
	//global_action_coord-=.47;
	moving_rsd.actcoord() = global_action_coord;

}


///////////////////////////////////////////////////////////////////////////////////////
void
insert_ideal_mainchain_bonds(
	Size const seqpos,
	Conformation & conformation
)
{
	using namespace id;

	Residue const & rsd( conformation.residue( seqpos ) );

	// create a mini-conformation with
	ConformationOP idl_op( new Conformation() );
	Conformation & idl = *idl_op;
	{
		ResidueOP idl_rsd( ResidueFactory::create_residue( rsd.type() ) );
		idl.append_residue_by_bond( *idl_rsd );
	}

	int idl_pos(1);

	// intra-residue mainchain bonds and angles
	Size const nbb( rsd.n_mainchain_atoms() );
	runtime_assert( nbb >= 2 ); // or logic gets a bit trickier
	for ( Size i=2; i<= nbb; ++i ) {
		AtomID
			bb1_idl( rsd.mainchain_atom(i-1), idl_pos ),
			bb1 ( rsd.mainchain_atom(i-1),  seqpos ),
			bb2_idl( rsd.mainchain_atom(  i), idl_pos ),
			bb2 ( rsd.mainchain_atom(  i),  seqpos );

		// bond length
		conformation.set_bond_length( bb1, bb2, idl.bond_length( bb1_idl, bb2_idl ) );

		if ( i<nbb ) {
			AtomID
				bb3_idl( rsd.mainchain_atom(i+1), idl_pos ),
				bb3 ( rsd.mainchain_atom(i+1),  seqpos );

			// bond angle
			conformation.set_bond_angle( bb1, bb2, bb3, idl.bond_angle( bb1_idl, bb2_idl, bb3_idl ) );
		}
	}

	if ( seqpos>1 && rsd.is_polymer() && !rsd.is_lower_terminus() ) {
		ResidueOP prev_rsd( ResidueFactory::create_residue( conformation.residue( seqpos-1 ).type() ) );
		idl.prepend_polymer_residue_before_seqpos( *prev_rsd, 1, true );
		idl_pos = 2;

		Size const nbb_prev( prev_rsd->n_mainchain_atoms() );
		AtomID
			bb1_idl( prev_rsd->mainchain_atom( nbb_prev ), idl_pos-1 ),
			bb1 ( prev_rsd->mainchain_atom( nbb_prev ),  seqpos-1 ),
			bb2_idl( rsd.mainchain_atom( 1 ), idl_pos ),
			bb2 ( rsd.mainchain_atom( 1 ),  seqpos ),
			bb3_idl( rsd.mainchain_atom( 2 ), idl_pos ),
			bb3 ( rsd.mainchain_atom( 2 ),  seqpos );
		conformation.set_bond_angle( bb1, bb2, bb3, idl.bond_angle( bb1_idl, bb2_idl, bb3_idl ) );
	}

	if ( seqpos < conformation.size() && rsd.is_polymer() && !rsd.is_upper_terminus() ) {
		ResidueOP next_rsd( ResidueFactory::create_residue( conformation.residue( seqpos+1 ).type() ) );
		idl.append_polymer_residue_after_seqpos( *next_rsd, idl_pos, true );

		AtomID
			bb1_idl( rsd.mainchain_atom( nbb-1 ), idl_pos ),
			bb1 ( rsd.mainchain_atom( nbb-1 ),  seqpos ),
			bb2_idl( rsd.mainchain_atom( nbb   ), idl_pos ),
			bb2 ( rsd.mainchain_atom( nbb   ),  seqpos ),
			bb3_idl( next_rsd->mainchain_atom( 1 ), idl_pos+1 ),
			bb3 ( next_rsd->mainchain_atom( 1 ),  seqpos+1 );

		// bond angle
		conformation.set_bond_angle( bb1, bb2, bb3, idl.bond_angle( bb1_idl, bb2_idl, bb3_idl ) );

		// bond length
		conformation.set_bond_length( bb2, bb3, idl.bond_length( bb2_idl, bb3_idl ) );
	}


}

/// @details  Just sets the two bond angles and the bond length between seqpos and seqpos+1

void
insert_ideal_bonds_at_polymer_junction(
	Size const seqpos,
	Conformation & conformation
)
{
	using namespace id;

	Residue const & rsd( conformation.residue( seqpos ) );
	runtime_assert( seqpos < conformation.size() && rsd.is_polymer() && !rsd.is_upper_terminus() );

	if ( TR.Warning.visible() && conformation.fold_tree().is_cutpoint( seqpos ) ) {
		TR.Warning << "insert_ideal_bonds_at_polymer_junction: seqpos (" << seqpos << ") is a foldtree cutpoint, " <<
			"returning!" << std::endl;
	}

	// create a mini-conformation with ideal bond lengths and angles
	ConformationOP idl_op( new Conformation() );
	Conformation & idl = *idl_op;
	{
		ResidueOP idl_rsd( ResidueFactory::create_residue( rsd.type() ) );
		idl.append_residue_by_bond( *idl_rsd );
	}

	int idl_pos(1);

	Size const nbb( rsd.n_mainchain_atoms() );
	runtime_assert( nbb >= 2 ); // or logic gets a bit trickier

	ResidueOP next_rsd( ResidueFactory::create_residue( conformation.residue( seqpos+1 ).type() ) );
	idl.append_polymer_residue_after_seqpos( *next_rsd, idl_pos, true /* append with ideal geometry */ );

	AtomID
		bb1_idl( rsd.mainchain_atom( nbb-1 ), idl_pos ),
		bb1 ( rsd.mainchain_atom( nbb-1 ),  seqpos ),
		bb2_idl( rsd.mainchain_atom( nbb   ), idl_pos ),
		bb2 ( rsd.mainchain_atom( nbb   ),  seqpos ),
		bb3_idl( next_rsd->mainchain_atom( 1 ), idl_pos+1 ),
		bb3 ( next_rsd->mainchain_atom( 1 ),  seqpos+1 ),
		bb4_idl( next_rsd->mainchain_atom( 2 ), idl_pos+1 ),
		bb4 ( next_rsd->mainchain_atom( 2 ),  seqpos+1 );

	// bond angle
	conformation.set_bond_angle( bb1, bb2, bb3, idl.bond_angle( bb1_idl, bb2_idl, bb3_idl ) );

	// bond length
	conformation.set_bond_length( bb2, bb3, idl.bond_length( bb2_idl, bb3_idl ) );

	// bond angle
	conformation.set_bond_angle( bb2, bb3, bb4, idl.bond_angle( bb2_idl, bb3_idl, bb4_idl ) );

	// fix up atoms that depend on the atoms across the bond for their torsion offsets
	conformation.rebuild_polymer_bond_dependent_atoms( seqpos );

}


///////////////////////////////////////////////////////////////////////////////////////
void
idealize_position(
	Size const seqpos,
	Conformation & conformation
)
{
	using namespace id;

	runtime_assert( seqpos > 0 );
	runtime_assert( conformation.size() >= seqpos );
	runtime_assert( conformation.size() > 0 );
	Residue const & rsd( conformation.residue( seqpos ) );

	//// create a mini-conformation with completely ideal residue ( and nbrs, if appropriate )

	ConformationOP idl_op( new Conformation() );
	Conformation & idl = *idl_op;
	{
		ResidueOP idl_rsd( ResidueFactory::create_residue( rsd.type() ) );
		idl.append_residue_by_bond( *idl_rsd );
	}

	int idl_pos(1);
	bool lower_connect( false ), upper_connect( false );

	if ( rsd.is_polymer() ) {
		// add polymer nbrs?
		if ( seqpos > 1 && !rsd.is_lower_terminus() && !conformation.fold_tree().is_cutpoint( seqpos-1 ) ) {
			lower_connect = true;
			ResidueOP prev_rsd( ResidueFactory::create_residue( conformation.residue( seqpos-1 ).type() ) );
			idl.prepend_polymer_residue_before_seqpos( *prev_rsd, 1, true );
			idl_pos = 2;
		}

		if ( seqpos < conformation.size() && !rsd.is_upper_terminus() && !conformation.fold_tree().is_cutpoint( seqpos ) ) {
			upper_connect = true;
			ResidueOP next_rsd( ResidueFactory::create_residue( conformation.residue( seqpos+1 ).type() ) );
			idl.append_polymer_residue_after_seqpos( *next_rsd, idl_pos, true );
		}
	}

	//// now set the torsion angles in the ideal conformation... This is to prepare for replacing rsd with
	//// the idealized version

	// chi angles
	for ( Size i=1; i<= rsd.nchi(); ++i ) {
		idl.set_torsion( TorsionID( idl_pos, CHI, i ), rsd.chi( i ) );
	}
	// mainchain torsions, if polymer residue
	if ( rsd.is_polymer() ) {
		Size const nbb( rsd.n_mainchain_atoms() );
		for ( Size i=1; i<= nbb; ++i ) {
			//if ( ( !lower_connect && i == 1 ) || ( !upper_connect && i >= nbb-1 ) ) continue;
			idl.set_torsion( TorsionID( idl_pos, BB, i ), rsd.mainchain_torsion( i ) );
		}
	}


	//// now we copy the mainchain bond geometry from ideal conf into the conf to be idealized
	if ( rsd.is_polymer() ) {

		// intra-residue mainchain bonds and angles
		Size const nbb( rsd.n_mainchain_atoms() );
		runtime_assert( nbb >= 2 ); // or logic gets a bit trickier
		for ( Size i=2; i<= nbb; ++i ) {
			AtomID
				bb1_idl( rsd.mainchain_atom(i-1), idl_pos ),
				bb1 ( rsd.mainchain_atom(i-1),  seqpos ),
				bb2_idl( rsd.mainchain_atom(  i), idl_pos ),
				bb2 ( rsd.mainchain_atom(  i),  seqpos );

			// bond length
			conformation.set_bond_length( bb1, bb2, idl.bond_length( bb1_idl, bb2_idl ) );

			if ( i<nbb ) {
				AtomID
					bb3_idl( rsd.mainchain_atom(i+1), idl_pos ),
					bb3 ( rsd.mainchain_atom(i+1),  seqpos );

				// bond angle
				conformation.set_bond_angle( bb1, bb2, bb3, idl.bond_angle( bb1_idl, bb2_idl, bb3_idl ) );
			}
		}

		if ( lower_connect ) {
			Residue const & prev_rsd( idl.residue( idl_pos-1 ) );

			Size const nbb_prev( prev_rsd.n_mainchain_atoms() );
			AtomID
				bb1_idl( prev_rsd.mainchain_atom( nbb_prev ), idl_pos-1 ),
				bb1 ( prev_rsd.mainchain_atom( nbb_prev ),  seqpos-1 ),
				bb2_idl(   rsd.mainchain_atom(  1 ), idl_pos   ),
				bb2 (   rsd.mainchain_atom(  1 ),  seqpos   ),
				bb3_idl(   rsd.mainchain_atom(  2 ), idl_pos   ),
				bb3 (   rsd.mainchain_atom(  2 ),  seqpos   );
			conformation.set_bond_angle( bb1, bb2, bb3, idl.bond_angle( bb1_idl, bb2_idl, bb3_idl ) );
		}

		if ( upper_connect ) {
			Residue const & next_rsd( idl.residue( idl_pos+1 ) );

			AtomID
				bb1_idl(   rsd.mainchain_atom( nbb-1 ), idl_pos   ),
				bb1 (   rsd.mainchain_atom( nbb-1 ),  seqpos   ),
				bb2_idl(   rsd.mainchain_atom( nbb   ), idl_pos   ),
				bb2 (   rsd.mainchain_atom( nbb   ),  seqpos   ),
				bb3_idl( next_rsd.mainchain_atom(  1 ), idl_pos+1 ),
				bb3 ( next_rsd.mainchain_atom(  1 ),  seqpos+1 );

			// bond angle
			conformation.set_bond_angle( bb1, bb2, bb3, idl.bond_angle( bb1_idl, bb2_idl, bb3_idl ) );

			// bond length
			conformation.set_bond_length( bb2, bb3, idl.bond_length( bb2_idl, bb3_idl ) );
		}
	} // rsd.is_polymer()


	//// now we orient the ideal residue onto the existing residue and replace it
	ResidueOP idl_rsd( idl.residue( idl_pos ).clone() );

	// EXTREMELY SUBTLE DANGEROUS THING ABOUT HAVING Residue const &'s lying around:
	// if we put "rsd" in place of conformation.residue( seqpos ) below, this will not behave properly.
	// that's because the xyz's in rsd don't get updated properly to reflect the changes to the conformation
	// a "refold" is only triggered by another call to conformation.residue
	// May want to consider making Residue smarter about this kind of thing, ie knowing its from a
	// Conformation... acting as an observer... could get pretty tricky though...
	//
	idl_rsd->orient_onto_residue( conformation.residue( seqpos ) ); // can't use rsd here!
	conformation.replace_residue( seqpos, *idl_rsd, false );

}


///////////////////////////////////////////////////////////////////////////////////////
// NOTE: conformation is needed for polymer neighbours
// @author Barak 07/2009
bool
is_ideal_position(
	Size const seqpos,
	Conformation const & conf,
	Real theta_epsilon,
	Real D_epsilon
)
{
	using namespace core::kinematics;
	runtime_assert( conf.size() >= seqpos  &&  seqpos >= 1 );

	Residue const rsd( conf.residue( seqpos ) );

	// I. Create mini-conformations for both idealized and original residue + nbrs, if appropriate )
	ConformationOP miniconf_op( new Conformation() );
	Conformation & miniconf = *miniconf_op;
	ConformationOP miniconf_idl_op( new Conformation() );
	Conformation & miniconf_idl = *miniconf_idl_op;
	{
		miniconf.append_residue_by_bond( rsd );
		ResidueOP prsd_idl( ResidueFactory::create_residue( rsd.type() ) );
		miniconf_idl.append_residue_by_bond( *prsd_idl );
	}
	// add polymer nbrs if needed
	int seqpos_miniconf(1);
	//bool lower_connect( false ), upper_connect( false );
	if ( rsd.is_polymer() ) {
		// prepending previous neighbour
		if ( seqpos > 1 && !rsd.is_lower_terminus() && !conf.fold_tree().is_cutpoint( seqpos-1 ) ) {
			//lower_connect = true;  // set but never used ~Labonte
			seqpos_miniconf = 2; // pushed forward by prependedq residue
			Residue const & prev_rsd = conf.residue( seqpos-1 );
			miniconf.safely_prepend_polymer_residue_before_seqpos // TODO: safely is probably unneeded here, verify this point
				( prev_rsd, 1, false /*build_ideal_geom*/);
			ResidueOP pprev_rsd_idl( ResidueFactory::create_residue( prev_rsd.type() ) );
			miniconf_idl.safely_prepend_polymer_residue_before_seqpos
				( *pprev_rsd_idl, 1, true /*build_ideal_geom*/);
		}
		// appending next neighbour
		if ( seqpos < conf.size() && !rsd.is_upper_terminus() && !conf.fold_tree().is_cutpoint( seqpos ) ) {
			//upper_connect = true;  // set but never used ~Labonte
			Residue const & next_rsd = conf.residue( seqpos+1 );
			miniconf.safely_append_polymer_residue_after_seqpos
				( next_rsd, seqpos_miniconf, false /*build_ideal_geom*/ );
			ResidueOP pnext_rsd_idl( ResidueFactory::create_residue( next_rsd.type() ) );
			miniconf_idl.safely_append_polymer_residue_after_seqpos
				( *pnext_rsd_idl, seqpos_miniconf, true /*build_ideal_geom*/ );
		}
	}


	// II. Compare angle and length of all residue bonded atoms in the new mini-conformations
	for ( Size atom = 1, atom_end = miniconf.residue( seqpos_miniconf ).natoms();
			atom <= atom_end;
			++atom ) {
		id::AtomID aid(atom, seqpos_miniconf);
		if ( miniconf.atom_tree().atom(aid).is_jump() ) { // skip non-bonded atoms
			continue;
		}
		id::DOF_ID dof_theta(aid, id::THETA); // bond angle
		id::DOF_ID dof_D(aid, id::D); // bond length
		if ( ! numeric::equal_by_epsilon(
				miniconf.dof( dof_theta ), miniconf_idl.dof( dof_theta ), theta_epsilon ) ) {
			if ( TR.visible() ) {
				TR << "Non-ideal residue detected: "
					<< " Residue #" << seqpos << " atom #" << atom << "( " << rsd.atom_name( atom ) << " ) " << ": "
					<< " Ideal theta=" << miniconf_idl.dof( dof_theta)
					<< ", Inspected theta=" << miniconf.dof( dof_theta )
					<< " (in Radians)"
					<< std::endl;
			}
			return false;
		}
		if ( ! numeric::equal_by_epsilon(
				miniconf.dof( dof_D), miniconf_idl.dof (dof_D), D_epsilon ) ) {
			if ( TR.visible() ) {
				TR << "Non-ideal residue detected: "
					<< " Residue #" << seqpos << " atom #" << atom << "( " << rsd.atom_name( atom ) << " ) " << ": "
					<< " Ideal D=" << miniconf_idl.dof( dof_D)
					<< ", Inspected D=" << miniconf.dof( dof_D )
					<< std::endl;
			}
			return false;
		}
	}

	return true;
}


/// @brief  Fills coords of target_rsd with coords from source_rsd of same atom_name, rebuilds others.
/// @details  If preserve_only_sidechain_dihedrals is true, then this function only copies mainchain coordinates,
/// and rebuilds all sidechain coordinates from scratch, setting side-chain dihedrals based on the source residue.
/// Otherwise, if false, it copies all the atoms that it can from the source residue, then rebuilds the rest.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
copy_residue_coordinates_and_rebuild_missing_atoms(
	Residue const & source_rsd,
	Residue & target_rsd,
	Conformation const & conformation,
	bool const preserve_only_sidechain_dihedrals
) {
	target_rsd.seqpos( source_rsd.seqpos() ); // in case fill_missing_atoms needs context info
	target_rsd.chain ( source_rsd.chain () );

	if ( preserve_only_sidechain_dihedrals ) { //If we're keeping only the sidechain dihedral information
		Size const natoms(target_rsd.natoms());
		Size const nmainchain( target_rsd.n_mainchain_atoms() );
		bool any_missing(false);
		utility::vector1< bool > missing( natoms, false );

		for ( Size ia=1; ia<=natoms; ++ia ) { //Loop through all atoms
			std::string const & atom_name( target_rsd.atom_name(ia) );
			if ( ia <= nmainchain && source_rsd.has( atom_name ) ) { //If this is a mainchain atom and it's present in both target and source
				target_rsd.atom(ia).xyz( source_rsd.atom(atom_name).xyz() );
			} else {
				if ( TR.Debug.visible() ) TR.Debug << "copy_residue_coordinates_and_rebuild_missing_atoms: missing atom " << target_rsd.name() << ' ' << atom_name << std::endl;
				any_missing=true;
				missing[ia]=true;
			}
		}
		if ( any_missing ) {
			target_rsd.fill_missing_atoms( missing, conformation );
		}

		//Set side-chain dihedrals
		for ( core::Size i=1, imax=std::min( target_rsd.nchi(), source_rsd.nchi() ); i<=imax; ++i ) { //Loop through all shared chi angles
			target_rsd.set_chi(i, source_rsd.chi(i)); //Set side-chain dihedral.
		}

	} else { //If we're trying to keep the sidechain atom position information
		copy_residue_coordinates_and_rebuild_missing_atoms( source_rsd, target_rsd, conformation );
	}
	return;
}


/// @details  For building variant residues, eg
/// @note  Need conformation for context in case we have to rebuild atoms, eg backbone H
void
copy_residue_coordinates_and_rebuild_missing_atoms(
	Residue const & source_rsd,
	Residue & target_rsd,
	Conformation const & conformation
)
{
	target_rsd.seqpos( source_rsd.seqpos() ); // in case fill_missing_atoms needs context info
	target_rsd.chain ( source_rsd.chain () );

	Size const natoms( target_rsd.natoms() );

	utility::vector1< bool > missing( natoms, false );
	bool any_missing( false );

	for ( Size i=1; i<= natoms; ++i ) {
		std::string const & atom_name( target_rsd.atom_name(i) );
		if ( source_rsd.has( atom_name ) ) {
			target_rsd.atom( i ).xyz( source_rsd.atom( atom_name ).xyz() );
		} else {
			if ( TR.Debug.visible() ) TR.Debug << "copy_residue_coordinates_and_rebuild_missing_atoms: missing atom " << target_rsd.name() << ' ' << atom_name << std::endl;
			any_missing = true;
			missing[i] = true;
		}
	}

	if ( any_missing ) {
		target_rsd.fill_missing_atoms( missing, conformation );
	}

}

/// @details  Helper function for below
std::ostream &
print_atom( id::AtomID const & id, Conformation const & conf, std::ostream & os )
{
	Residue const & rsd( conf.residue(id.rsd() ) );
	os << rsd.atom_name( id.atomno() ) << " (" << id.atomno() << ") " << rsd.name() << ' ' << rsd.seqpos();
	return os;
}


void
show_atom_tree(
	kinematics::tree::Atom const & atom,
	Conformation const & conf, std::ostream & os
) {
	os << "ATOM "; print_atom( atom.id(), conf, os ) << std::endl;
	os << "CHILDREN: ";
	for ( Size i=0; i< atom.n_children(); ++i ) {
		print_atom( atom.child(i)->id(), conf, os ) << ' ';
	}
	os << std::endl;
	for ( Size i=0; i< atom.n_children(); ++i ) {
		show_atom_tree( *atom.child(i), conf, os );
	}
}

/// helper function for residue replacement/residuetype switching
/// @note  Will call new_rsd->fill_missing_atoms if the new residue has atoms
/// that the old one doesn't

void
replace_conformation_residue_copying_existing_coordinates(
	conformation::Conformation & conformation,
	Size const seqpos,
	chemical::ResidueType const & new_rsd_type
)
{

	Residue const & old_rsd( conformation.residue( seqpos ) );
	ResidueOP new_rsd( ResidueFactory::create_residue( new_rsd_type ) );
	conformation::copy_residue_coordinates_and_rebuild_missing_atoms( old_rsd, *new_rsd, conformation );
	conformation.replace_residue( seqpos, *new_rsd, false );

}


/// @details E.g., make a terminus variant, and replace the original in pose.
/// @note This copies any atoms in common between old and new residues, rebuilding the others.
void
add_variant_type_to_conformation_residue(
	conformation::Conformation & conformation,
	chemical::VariantType const variant_type,
	Size const seqpos )
{
	Residue const & old_rsd( conformation.residue( seqpos ) );

	// the type of the desired variant residue
	chemical::ResidueTypeSet const & rsd_set( old_rsd.residue_type_set() );
	chemical::ResidueType const & new_rsd_type( rsd_set.get_residue_type_with_variant_added( old_rsd.type(), variant_type ) );

	core::conformation::replace_conformation_residue_copying_existing_coordinates( conformation, seqpos, new_rsd_type );
}


/// @details E.g., remove a terminus variant, and replace the original in pose.
/// @note This copies any atoms in common between old and new residues, rebuilding the others.
void
remove_variant_type_from_conformation_residue(
	conformation::Conformation & conformation,
	chemical::VariantType const variant_type,
	Size const seqpos )
{
	Residue const & old_rsd( conformation.residue( seqpos ) );

	// the type of the desired variant residue
	chemical::ResidueTypeSet const & rsd_set( old_rsd.residue_type_set() );
	chemical::ResidueType const & new_rsd_type( rsd_set.get_residue_type_with_variant_removed( old_rsd.type(), variant_type ) );

	core::conformation::replace_conformation_residue_copying_existing_coordinates( conformation, seqpos, new_rsd_type );
}


///////////////////////////////////////////////////////////////////////////////
void
add_lower_terminus_type_to_conformation_residue(
	conformation::Conformation & conformation,
	Size const seqpos
)
{
	core::conformation::add_variant_type_to_conformation_residue( conformation, chemical::LOWER_TERMINUS_VARIANT, seqpos );
}

///////////////////////////////////////////////////////////////////////////////
void
remove_lower_terminus_type_from_conformation_residue(
	conformation::Conformation & conformation,
	Size const seqpos
)
{
	remove_variant_type_from_conformation_residue( conformation, chemical::LOWER_TERMINUS_VARIANT, seqpos );
}

///////////////////////////////////////////////////////////////////////////////
void
add_upper_terminus_type_to_conformation_residue(
	conformation::Conformation & conformation,
	Size const seqpos
)
{
	core::conformation::add_variant_type_to_conformation_residue( conformation, chemical::UPPER_TERMINUS_VARIANT, seqpos );
}

///////////////////////////////////////////////////////////////////////////////
void
remove_upper_terminus_type_from_conformation_residue(
	conformation::Conformation & conformation,
	Size const seqpos
)
{
	remove_variant_type_from_conformation_residue( conformation, chemical::UPPER_TERMINUS_VARIANT, seqpos );
}


///////////////////////////////////////////////////////////////////////////////
/// @details  Build an atom-tree from a fold-tree and a set of residues
/// atoms in the tree are allocated with new (on the heap) and held in owning pointers in atom_pointer
/// @note  atom_pointer is cleared at the beginning and then filled with AtomOPs to all the atoms


void
build_tree(
	kinematics::FoldTree const & fold_tree,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomPointer2D & atom_pointer
)
{
	using conformation::Residue;

	// note -- we clear atom_pointer
	atom_pointer.clear();
	atom_pointer.resize( residues.size() );

	// build first atom. make it a "JumpAtom"
	{
		int const start_pos( fold_tree.root() );
		Residue const & rsd( *(residues[start_pos]) );

		build_residue_tree( residues, rsd, fold_tree, atom_pointer[ start_pos ] );

	}

	// traverse tree, build edges
	for ( kinematics::FoldTree::const_iterator it = fold_tree.begin(),
			it_end = fold_tree.end(); it != it_end; ++it ) {

		if ( it->is_jump() ) {
			build_jump_edge( *it, residues, atom_pointer );

		} else if ( it->label() == kinematics::Edge::PEPTIDE ) {
			// build a peptide edge
			build_polymer_edge( *it, residues, atom_pointer );

		} else if ( it->label() == kinematics::Edge::CHEMICAL ) {
			build_chemical_edge( *it, residues, atom_pointer );

		} else {
			if ( TR.Error.visible() ) {
				TR.Error << "Failed to identify kinematics::Edge label in core/kinematics/util.cc::build_tree()" << std::endl;
				TR.Error << "Label = " << it->label() << std::endl;
			}
			utility_exit();
		}
	}

	// now guarantee that jump stubs are residue-internal if desired
	for ( kinematics::FoldTree::const_iterator it = fold_tree.begin(),
			it_end = fold_tree.end(); it != it_end; ++it ) {
		if ( it->is_jump() && it->keep_stub_in_residue() ) {
			promote_sameresidue_child_of_jump_atom( *it, residues, atom_pointer );
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
/// @details root_atomno is the root for the sub atom-tree of this edge.
/// anchor_atomno is the entry point of this sub atom-tree into the main atom-tree.

void
build_jump_edge(
	kinematics::Edge const & edge,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomPointer2D & atom_pointer
)
{
	debug_assert( edge.is_jump() );

	int const estart  ( edge.start() );
	int const estop   ( edge.stop () );

	// these may have been set in the edge
	Size anchor_atomno, root_atomno;
	get_anchor_and_root_atoms( *residues[ estart ], *residues[ estop ], edge, anchor_atomno, root_atomno );

	// get the anchor atom
	kinematics::tree::AtomOP anchor_atom( atom_pointer( id::AtomID( anchor_atomno, estart ) ) );
	debug_assert( anchor_atom );

	// build the new residue
	build_residue_tree( root_atomno, *residues[ estop ], atom_pointer[ estop ], true /*Jump*/ );

	// now wire in the new residue connection
	anchor_atom->insert_atom( atom_pointer[ id::AtomID( root_atomno, estop ) ] );


	//  std::cout << "build_jump_edge: " << edge << ' ' <<  estop << ' ' << root_atomno << ' ' <<
	//   estart << ' ' << anchor_atomno << std::endl;

	//  std::cout << "TREE AFTER BEGIN: " << edge << std::endl;
	//  anchor_atom->show();
	// debug_assert( anchor_atom->id() == AtomID( anchor_atomno, estart ) );
	//  std::cout << "TREE AFTER END: " << edge << std::endl;

}


///////////////////////////////////////////////////////////////////////////////
/// @details assumes that the start residue of edge has already been built. Traverse
///the polymer edge residue by residue and after building sub atom-tree for this residue,
///attaches the edge's subtree to the anchor_atom in the previous residue.
///
void
build_polymer_edge(
	kinematics::Edge const & edge,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomPointer2D & atom_pointer
)
{
	int const start( edge.start() );
	int const stop ( edge.stop() );
	int const dir  ( edge.polymer_direction() );

	debug_assert( dir == 1 || dir == -1 );

	id::AtomID first_anchor;
	for ( int pos=start+dir; pos != stop + dir; pos += dir ) {
		runtime_assert( pos > 0 );
		conformation::Residue const & rsd( *residues[pos] );
		int const anchor_pos( pos-dir );
		Size anchor_atomno, root_atomno;
		get_anchor_and_root_atoms( *residues[ anchor_pos ], rsd, edge, anchor_atomno, root_atomno );

		// build the new residue tree, fill in the atom_pointer array
		build_residue_tree( root_atomno, rsd, atom_pointer[pos], false /*Jump*/ );

		kinematics::tree::AtomOP anchor_atom( atom_pointer( id::AtomID( anchor_atomno, anchor_pos ) ) );
		kinematics::tree::AtomOP   root_atom( atom_pointer( id::AtomID(   root_atomno,  pos ) ) );
		debug_assert( anchor_atom && root_atom );
		anchor_atom->insert_atom( root_atom );
		//std::cout << "build_polymer_edge: " << edge << ' ' <<  pos << ' ' << root_atomno << ' ' <<
		//   anchor_pos << ' ' << anchor_atomno << std::endl;
		//  if ( pos == start+dir ) first_anchor = AtomID( anchor_atomno, anchor_pos );
	}
	//  if ( start != stop ) {
	//   std::cout << "TREE AFTER BEGIN: " << edge << std::endl;
	//   atom_pointer[ first_anchor ]->show();
	//   std::cout << "TREE AFTER END: " << edge << std::endl;
	//  }
}

///////////////////////////////////////////////////////////////////////////////
/// @details assumes that the start residue of edge has already been built. Traverse
///the chemical edge residue by residue and after building sub atom-tree for this residue,
///attaches the edge's subtree to the anchor_atom in the previous residue.
///
void
build_chemical_edge(
	kinematics::Edge const & edge,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomPointer2D & atom_pointer
)
{

	int const estart  ( edge.start() );
	int const estop   ( edge.stop () );

	// these may have been set in the edge
	Size anchor_atomno, root_atomno;
	get_anchor_and_root_atoms( *residues[ estart ], *residues[ estop ], edge, anchor_atomno, root_atomno );

	build_residue_tree( root_atomno, *residues[ estop ], atom_pointer[ estop ], false /*Jump*/ );

	// get the anchor atom
	kinematics::tree::AtomOP anchor_atom( atom_pointer( id::AtomID( anchor_atomno, estart ) ) );
	kinematics::tree::AtomOP   root_atom( atom_pointer( id::AtomID(   root_atomno, estop  ) ) );
	debug_assert( anchor_atom && root_atom );

	anchor_atom->insert_atom( root_atom );
}

/////////////////////////////////////////////////////////////////////////////
///\brief get the root atom for building residue atom-tree given the folding direction "dir"
int
get_root_atomno(
	conformation::Residue const & rsd,
	int const dir // +1, -1, or "dir_jump"
)
{
	// if dir == 1 or -1, we are building a contiguous edge ("peptide edge")
	// in the fold_tree ( ==> residue should be a polymer type? )
	int const forward( 1 );
	int const backward( -1 );

	if ( dir == forward ) {
		// N for proteins, P for DNA
		if ( rsd.is_polymer() ) {
			if ( rsd.is_lower_terminus() && rsd.mainchain_atoms().size() > 0 ) {
				// this is a little strange to be folding through a terminal residue but kinematics doesn't
				// have to correlate perfectly with chemical
				return rsd.mainchain_atom(1);
			} else {
				return rsd.lower_connect_atom(); //mainchain_atoms()[1];
			}
		} else {
			return get_root_atomno( rsd, kinematics::dir_jump );
		}
	} else if ( dir == backward ) {
		if ( rsd.is_polymer() ) {
			if ( rsd.is_upper_terminus() && rsd.mainchain_atoms().size() > 0 ) {
				// this is a little strange to be folding through a terminal residue but kinematics doesn't
				// have to correlate perfectly with chemical
				return rsd.mainchain_atom( rsd.mainchain_atoms().size() );
			} else {
				// C for proteins, O3' for DNA
				return rsd.upper_connect_atom(); //mainchain_atoms()[ rsd.mainchain_atoms().size() ];
			}
		} else {
			return get_root_atomno( rsd, kinematics::dir_jump );
		}
	} else {
		debug_assert( dir == kinematics::dir_jump );
		// default for jumps, use the N-terminal attachment?
		if ( rsd.mainchain_atoms().empty() ) {
			// For non-polymeric ligands, match the ICOOR root.
			// (This should minimize issues with matching trees.)
			return rsd.type().atom_index( rsd.type().root_atom() );
		} else if ( rsd.type().has_variant_type( chemical::N_ACETYLATION ) ) {
			return 1; // awful hack! want N-CA-C triad to stay put when one of these residues is at root.
		} else {
			// PB 10/29/07 - changing to use an interior mainchain atom
			//
			// old choice of mainchain_atoms()[1] meant that changing omega of jump_pos-1 would propagate
			// C-terminal to jump_pos, which is bad for loop modeling where the usual assumption was that
			// we can use loop_begin-1 and loop_end+1 as jump points
			//
			Size const nbb( rsd.n_mainchain_atoms() ); //    1  2  3  4  5  6  -- nbb
			Size const mainchain_index( ( nbb-1 )/2 + 1 ); //   1  1  2  2  3  3  -- mainchain_index
			return rsd.mainchain_atoms()[ mainchain_index ];
		}
	}
}

/// @details  Determine which atom to use as the root of the root residue in the atomtree.
/// It is sometimes useful to be able to control the atom chosen as the root of the atomtree, eg in graphics.
/// The logic below uses atom_info stored in the foldtree for jump edges emanating from the root residue
/// to override the default atom (if the root residue is a jump_point and said atom_info exists)


Size
get_root_residue_root_atomno(
	conformation::Residue const & rsd,
	kinematics::FoldTree const & fold_tree
)
{
	Size const seqpos( rsd.seqpos() );
	debug_assert( seqpos == Size( fold_tree.root() ) ); // need to refactor foldtree to use Size instead of int

	Size root_atomno = get_root_atomno( rsd, kinematics::dir_jump ); // the default setting

	if ( fold_tree.is_jump_point( seqpos ) ) {
		for ( Size i=1; i<= fold_tree.num_jump(); ++i ) {
			kinematics::Edge const & edge( fold_tree.jump_edge( i ) );
			if ( seqpos == Size(edge.start()) && edge.has_atom_info() ) {
				root_atomno = rsd.atom_index( edge.upstream_atom() );
				break;
			}
		}
	}
	return root_atomno;
}

/// @details  A wrapper function to build an atom-tree for a residue. Uses information from the foldtree to
/// find the proper root atom if it is not defined, and determine the direction in which to build the residue.
/// The goal of this function is to allow the Conformation to rebuild parts of the AtomTree in a way that is
/// compatible with what would be built if we erased the entire atomtree and rebuilt it using the foldtree.
/// Also used to build the atomtree for the root residue when we are building the atomtree from scratch.

void
build_residue_tree(
	conformation::ResidueCOPs const & residues,
	conformation::Residue const & rsd,
	kinematics::FoldTree const & fold_tree,
	kinematics::AtomPointer1D & atom_ptr
)
{
	// determine root_atomno and whether root_atom is a jump_atom

	bool root_atom_is_jump_atom( false );
	int root_atomno(0);

	int const seqpos( rsd.seqpos() );

	if ( seqpos == fold_tree.root() ) {
		root_atom_is_jump_atom = true;

		root_atomno = get_root_residue_root_atomno( rsd, fold_tree );

	} else {
		// seqpos is not the root of the fold_tree
		kinematics::Edge const & edge( fold_tree.get_residue_edge( seqpos ) ); // the edge that builds seqpos

		if ( edge.is_peptide() ) {
			// peptide edge
			root_atomno = get_root_atomno( rsd, edge.polymer_direction() ); // the default setting
		} else if ( edge.is_jump() ) {
			// jump edge
			root_atom_is_jump_atom = true;
			if ( edge.has_atom_info() ) {
				root_atomno = rsd.atom_index( edge.downstream_atom() );
			} else {
				root_atomno = get_root_atomno( rsd, kinematics::dir_jump ); // the default setting
			}
		} else {
			//// chemical edge -- atomnumbers of connections should be programmed in
			//debug_assert( edge.has_atom_info() );
			//root_atomno = rsd.atom_index( edge.downstream_atom() );
			// Connection atoms are not always present in the edge, but might be;  else use defaults.
			int const estart  ( edge.start() );
			int const estop   ( edge.stop () );
			Size anchor_atno, root_atno;
			get_anchor_and_root_atoms( *residues[ estart ], *residues[ estop ], edge, anchor_atno, root_atno );
			root_atomno = (int) root_atno; // stupid Size/int confusion...
		}
	}
	debug_assert( root_atomno );

	// now call the main build_residue_tree function
	build_residue_tree( root_atomno, rsd, atom_ptr, root_atom_is_jump_atom );
}

///////////////////////////////////////////////////////////////////////////////
/// @brief Check if this atom neighbor has been black-listed ("CUT_BOND" in params file).
bool
check_good_neighbor( Size const & atom_index, utility::vector1< Size > const & cut_nbrs ) {
	for ( Size kk=1; kk<=cut_nbrs.size(); ++kk )  {
		if ( cut_nbrs[kk] == atom_index ) return false;
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////
/// @brief is atom2 the last atom of our chi angle ?
inline
bool
chi_continuation(
	Size const atom1, // current atom
	Size const atom2, // potential nbr
	utility::vector1< utility::vector1< Size > > const & chi_atoms
)
{
	int const nchi( chi_atoms.size() );

	for ( int i=1; i<= nchi; ++i ) {
		utility::vector1< Size > const & atoms( chi_atoms[i] );
		if ( atom1 == atoms[2] &&
				atom2 == atoms[1] ) return true;

		if ( atom1 == atoms[3] &&
				atom2 == atoms[4] ) return true;
	}
	return false;

}

///////////////////////////////////////////////////////////////////////////////
/// @brief would we be breaking into a chi angle by adding atom2 to the tree now?
inline
bool
chi_interruption(
	Size const atom1, // current atom
	Size const atom2, // potential nbr
	utility::vector1< utility::vector1< Size > > const & chi_atoms,
	utility::vector1< bool > const & is_done
)
{
	int const nchi( chi_atoms.size() );

	//not an interruption if atom1 and atom2 delineate a chi angle.
	for ( int i=1; i<= nchi; ++i ) {
		utility::vector1< Size > const & atoms( chi_atoms[i] );

		if ( atom2 == atoms[2] &&
				( atom1 == atoms[1] || atom1 == atoms[3] ) ) return false;

		if ( atom2 == atoms[3] &&
				( atom1 == atoms[2] || atom1 == atoms[4] ) ) return false;
	}


	for ( int i=1; i<= nchi; ++i ) {
		utility::vector1< Size > const & atoms( chi_atoms[i] );

		if ( atom2 == atoms[2] &&
				atom1 != atoms[1] &&
				atom1 != atoms[3] ) return true;

		if ( atom2 == atoms[3] &&
				atom1 != atoms[2] &&
				atom1 != atoms[4] ) return true;

		if ( (atom2 == atoms[1] && is_done[ atoms[4] ]) ||
				(atom2 == atoms[4] && is_done[ atoms[1] ]) ) return true;

	}
	return false;

}


///////////////////////////////////////////////////////////////////////////////
/// @brief simply fill the "links" by adding, for each atom, its bonded neighbors
void
setup_links_simple(
	conformation::Residue const & rsd,
	kinematics::Links & links
)
{
	int const natoms( rsd.natoms() );

	links.clear();
	links.resize( natoms );

	for ( int atomno1=1; atomno1<= natoms; ++atomno1 ) {
		utility::vector1< Size > const & nbrs( rsd.nbrs( atomno1 ) );
		utility::vector1< Size > const & cut_nbrs( rsd.cut_bond_neighbor( atomno1 ) );
		for ( Size jj=1; jj<= nbrs.size(); ++jj ) {
			if ( !check_good_neighbor( nbrs[jj], cut_nbrs ) ) continue;
			links[ atomno1 ].push_back( nbrs[jj] );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// called recursively

// Rule I -- dont enter a chi angle in the middle (atoms 2 or 3 ),
//     and don't enter a chi angle from one side if the other side's atoms have already been added to the tree
//
//
// Rule II -- if you're atom 2 and atom 1 hasn't been done, put that first
//   likewise for atoms 3 and 4
//
// Rule III -- mainchain atoms 1st, tiebreaking by rule II
//
// Rule IV -- heavyatoms before hydrogens, subject to other three rules
//
/// @brief set correct order for how atoms are linked to each other.
///
/// this function is called recursively.
/// @li atom1 is the root of links at the current level.
/// @li full_links store information about UNORDERED bonded neighbors for each atom in this residue.
/// @li is_done indicates which atoms have already been linked.
/// @li is_mainchain, is_chi, is_hydrogen and chi atoms are self explanatory
/// @li new_links store information about ORDERED bonded neighbors (links) for each atom in this residue.
void
setup_atom_links(
	int const atom1,
	kinematics::Links const & full_links,
	utility::vector1< bool > & is_done,
	utility::vector1< bool > const & is_mainchain,
	utility::vector1< bool > const & is_chi,
	utility::vector1< bool > const & is_hydrogen,
	utility::vector1< utility::vector1< Size > > const & chi_atoms,
	kinematics::Links & new_links
)
{
	is_done[atom1] = true;

	utility::vector1< Size > const & full_l( full_links[atom1] );
	utility::vector1< Size >    &  new_l(  new_links[atom1] );

	debug_assert( new_l.empty() );

	Size const n_nbrs( full_l.size() );

	if ( n_nbrs == 1 ) {
		Size const nbr( full_l[1] );
		if ( !is_done[nbr] ) {
			// special case -- we have no choice but to add this guy otherwise setup will fail
			new_l.push_back( nbr );
			setup_atom_links( nbr, full_links, is_done, is_mainchain, is_chi, is_hydrogen, chi_atoms, new_links );
		}
		return;
	}

	// top priority -- within a chi angle, and mainchain
	for ( Size i=1; i<= n_nbrs; ++i ) {
		int const atom2( full_l[i] );
		if ( is_done[ atom2 ] ) continue;

		if ( is_mainchain[ atom2 ] && is_chi[ atom1 ] && is_chi[ atom2 ] && chi_continuation( atom1, atom2, chi_atoms ) ) {
			new_l.push_back( atom2 );
			setup_atom_links( atom2, full_links, is_done, is_mainchain, is_chi,
				is_hydrogen, chi_atoms, new_links );
		}
	}

	// next priority -- mainchain
	for ( Size i=1; i<= n_nbrs; ++i ) {
		int const atom2( full_l[i] );
		if ( is_done[ atom2 ] ) continue;

		if ( is_mainchain[ atom2 ] ) {
			new_l.push_back( atom2 );
			setup_atom_links( atom2, full_links, is_done, is_mainchain, is_chi,
				is_hydrogen, chi_atoms, new_links );
		}
	}

	// next priority -- within a chi angle and heavy.
	for ( Size i=1; i<= n_nbrs; ++i ) {
		int const atom2( full_l[i] );
		if ( is_done[ atom2 ] ) continue;

		if ( is_chi[ atom1 ] && is_chi[ atom2 ] && chi_continuation( atom1, atom2, chi_atoms ) && !is_hydrogen[ atom2 ] ) {
			new_l.push_back( atom2 );
			setup_atom_links( atom2, full_links, is_done, is_mainchain, is_chi,
				is_hydrogen, chi_atoms, new_links );
		}
	}

	// next priority -- any heavy atom chi?
	for ( Size i=1; i<= n_nbrs; ++i ) {
		int const atom2( full_l[i] );
		if ( is_done[ atom2 ] ) continue;

		if ( is_chi[ atom2 ] && !is_hydrogen[ atom2 ] ) {
			new_l.push_back( atom2 );
			setup_atom_links( atom2, full_links, is_done, is_mainchain, is_chi,
				is_hydrogen, chi_atoms, new_links );
		}
	}


	// next priority -- any chi -- could be hydrogen
	for ( Size i=1; i<= n_nbrs; ++i ) {
		int const atom2( full_l[i] );
		if ( is_done[ atom2 ] ) continue;

		if ( is_chi[ atom2 ] ) {
			new_l.push_back( atom2 );
			setup_atom_links( atom2, full_links, is_done, is_mainchain, is_chi,
				is_hydrogen, chi_atoms, new_links );
		}
	}

	// next priority -- heavyatoms
	for ( Size i=1; i<= n_nbrs; ++i ) {
		int const atom2( full_l[i] );
		if ( is_done[ atom2 ] ) continue;
		if ( is_chi[ atom2 ] && chi_interruption( atom1, atom2, chi_atoms, is_done ) ) continue;
		if ( !is_hydrogen[ atom2 ] ) {
			new_l.push_back( atom2 );
			setup_atom_links( atom2, full_links, is_done, is_mainchain, is_chi,
				is_hydrogen, chi_atoms, new_links );
		}
	}


	// lowest priority -- hydrogens
	for ( Size i=1; i<= n_nbrs; ++i ) {
		int const atom2( full_l[i] );
		if ( is_done[ atom2 ] ) continue;
		if ( is_chi[ atom2 ] && chi_interruption( atom1, atom2, chi_atoms, is_done ) ) continue;
		new_l.push_back( atom2 );
		setup_atom_links( atom2, full_links, is_done, is_mainchain, is_chi,
			is_hydrogen, chi_atoms, new_links );
	}
}

///////////////////////////////////////////////////////////////////////////////
/// @brief given the root_atomno, set up rules for how other atoms are linked for this residue
///a wrapper function calling setup_atom_links recursively .
void
setup_links(
	conformation::Residue const & rsd,
	int const root_atomno,
	kinematics::Links & links
)
{
	typedef utility::vector1< bool > BVec;

	//////////////////////////////////
	// get the full list of atom nbrs:
	kinematics::Links full_links;

	setup_links_simple( rsd, full_links );

	// natoms
	Size const natoms( rsd.natoms() );

	////////////////
	// chi atoms
	BVec is_chi( natoms, false ); // vector bool
	utility::vector1< utility::vector1< Size > > const & chi_atoms
		( rsd.chi_atoms() );
	for ( Size i=1; i<= chi_atoms.size(); ++i ) {
		for ( int j=1; j<= 4; ++j ) {
			is_chi[ chi_atoms[i][j] ] = true;
		}
	}


	//////////////////
	// mainchain atoms
	BVec is_mainchain( natoms, false ); // vector bool
	utility::vector1< Size > const & mainchain( rsd.mainchain_atoms() );
	for ( Size i=1; i<= mainchain.size(); ++i ) {
		is_mainchain[ mainchain[i] ] = true;
	}
	if ( rsd.has_variant_type( chemical::CUTPOINT_LOWER ) ) {
		is_mainchain[ rsd.atom_index( "OVL1" ) ] = true;
		if ( rsd.has( "OVL2" ) ) {
			is_mainchain[ rsd.atom_index( "OVL2" ) ] = true;
		}
	}
	if ( rsd.has_variant_type( chemical::CUTPOINT_UPPER ) ) {
		is_mainchain[ rsd.atom_index( "OVU1" ) ] = true;
	}

	////////////
	// hydrogens
	BVec is_hydrogen( natoms, false );
	for ( Size i=1; i<= natoms; ++i ) {
		is_hydrogen[i] = rsd.atom_type(i).is_hydrogen();
	}


	links.clear();
	links.resize( natoms );

	BVec is_done( natoms, false ); // keep track of what's done
	setup_atom_links( root_atomno, full_links, is_done, is_mainchain,  is_chi,
		is_hydrogen, chi_atoms, links );


	// check for an atom that wasn't added
	for ( Size i=1; i<= natoms; ++i ) {
		if ( !is_done[ i ] ) {
			if ( TR.Error.visible() ) {
				TR.Error << "Unable to setup links for residue " << std::endl;
				for ( Size j=1; j<= natoms; ++j ) {
					TR.Error << rsd.atom_name(j) << ' ' << is_done[j] << ' ' << is_mainchain[j] << ' ' << is_chi[j] << ' ' <<
						is_hydrogen[j] << std::endl;
				}
			}
			utility_exit();
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
///\brief set up a local atom-tree for a residue from the defined root atom.
void
build_residue_tree(
	int const root_atomno,
	conformation::Residue const & rsd,
	kinematics::AtomPointer1D & atom_ptr,
	bool const root_is_jump_atom
)
{
	Size const natoms( rsd.natoms() );

	atom_ptr.clear();
	atom_ptr.resize( natoms ); // default C-TOR for AtomOP is 0 initialized

	// setup the atom nbrs so that the tree will match the desired torsion
	// angles
	kinematics::Links links;
	setup_links( rsd, root_atomno, links );

	// now build the residue atom-tree using the recursive function add_atom
	kinematics::add_atom( root_atomno, rsd.seqpos(), links, atom_ptr, root_is_jump_atom );

	// fill in the coordinates from the residue into the atomtree
	for ( Size i=1; i<= natoms; ++i ) {
		atom_ptr[i]->xyz( rsd.atom(i).xyz() );
	}
}

/// @details  Get the incoming connection and all outgoing connections from a residue

void
get_residue_connections(
	conformation::Residue const & new_rsd,
	kinematics::FoldTree const & fold_tree,
	conformation::ResidueCOPs const & residues,
	id::BondID & new_rsd_in,
	utility::vector1< id::BondID > & new_rsd_out
)
{

	Size const seqpos( new_rsd.seqpos() );

	// setup incoming connection
	if ( fold_tree.is_root( seqpos ) ) {
		new_rsd_in.atom1 = id::BOGUS_ATOM_ID;
		new_rsd_in.atom2 = id::AtomID( get_root_residue_root_atomno( new_rsd, fold_tree ), seqpos );
	} else {
		Size anchor_atomno, root_atomno, anchor_pos;
		kinematics::Edge const & edge( fold_tree.get_residue_edge( seqpos ) );
		if ( edge.is_polymer() ) anchor_pos = seqpos - edge.polymer_direction();
		else anchor_pos = edge.start();
		get_anchor_and_root_atoms( *residues[ anchor_pos ], new_rsd, edge, anchor_atomno, root_atomno );
		new_rsd_in.atom1 = id::AtomID( anchor_atomno, anchor_pos );
		new_rsd_in.atom2 = id::AtomID( root_atomno, seqpos );
	}


	// setup the outgoing connections
	new_rsd_out.clear();
	utility::vector1< kinematics::Edge > const outgoing_edges( fold_tree.get_outgoing_edges( seqpos ) );
	for ( utility::vector1< kinematics::Edge >::const_iterator it= outgoing_edges.begin(); it != outgoing_edges.end(); ++it ) {
		Size anchor_atomno, root_atomno;
		Size const root_pos( ( it->is_polymer() ) ? seqpos + it->polymer_direction() : it->stop() );
		get_anchor_and_root_atoms( new_rsd, *residues[ root_pos ], *it, anchor_atomno, root_atomno );
		new_rsd_out.push_back( id::BondID( id::AtomID( anchor_atomno, seqpos ), id::AtomID( root_atomno, root_pos ) ) );
	}
}


/// @details  Helper function for conformation routines.
/// Uses fold_tree to deduce the incoming/outgoing connections for the new residue and the old residue
/// We want it to be the case that the tree we get after this call is the same one that we would have gotten
/// by calling build_tree
void
replace_residue_in_atom_tree(
	conformation::Residue const & new_rsd,
	//conformation::Residue const & old_rsd,
	kinematics::FoldTree const & fold_tree,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomTree & atom_tree
)
{
	id::BondID new_rsd_in;
	utility::vector1< id::BondID > new_rsd_out;

	get_residue_connections( new_rsd, fold_tree, residues, new_rsd_in, new_rsd_out );

	// build the new atoms
	kinematics::AtomPointer1D new_atoms;
	build_residue_tree( residues, new_rsd, fold_tree, new_atoms );

	// now replace the atoms
	atom_tree.replace_residue_subtree( new_rsd_in, new_rsd_out, new_atoms );

	// preserve same-residue jump status if necessary
	if ( fold_tree.is_jump_point( new_rsd.seqpos() ) && !fold_tree.is_root( new_rsd.seqpos() ) ) {
		kinematics::Edge const & edge( fold_tree.get_residue_edge( new_rsd.seqpos() ) );
		if ( edge.is_jump() && edge.keep_stub_in_residue() ) {
			promote_sameresidue_child_of_jump_atom( edge, residues, atom_tree );
		}
	}
}

/// @brief  Inserts/ appends new residue subtree into an existing atomtree
/// @note  The foldtree must already have been changed to reflect the new residue
/// @note  The residues array should already have been inserted into
/// @note  The sequence position of the new residue is deduced from new_rsd.seqpos()
/// @note  This function handles renumbering of the atomtree if necessary
void
insert_residue_into_atom_tree(
	conformation::Residue const & new_rsd,
	kinematics::FoldTree const & fold_tree,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomTree & atom_tree
)
{

	Size const seqpos( new_rsd.seqpos() );
	Size const nres( fold_tree.nres() );
	Size const old_nres( atom_tree.size() );
	debug_assert( nres == residues.size() && seqpos <= nres && old_nres == nres-1 );

	/////////////////////////////////
	// setup for renumbering atomtree
	utility::vector1< int > old2new( old_nres, 0 );
	for ( Size i=1; i<= old_nres; ++i ) {
		if ( i< seqpos ) old2new[i] = i;
		else old2new[i] = i+1;
	}

	// renumber the atomtree, fold_tree
	atom_tree.update_sequence_numbering( nres, old2new );
	debug_assert( atom_tree.size() == nres );

	replace_residue_in_atom_tree( new_rsd, fold_tree, residues, atom_tree );

}


/// @details  Get the atom-index of the atom to which the residue at position seqpos should be anchored in
/// constructing the atomtree.

int
get_anchor_atomno( conformation::Residue const & anchor_rsd, Size const seqpos, kinematics::FoldTree const & fold_tree )
{
	int anchor_atomno(0);
	debug_assert( seqpos != (Size)fold_tree.root() );

	kinematics::Edge const & edge( fold_tree.get_residue_edge( seqpos ) );

	if ( edge.is_jump() ) {
		// jump edge
		debug_assert( (seqpos == (Size)edge.stop()) && ((Size)anchor_rsd.seqpos() == (Size)edge.start()) );
		if ( edge.has_atom_info() ) {
			anchor_atomno = anchor_rsd.atom_index( edge.upstream_atom() );
		} else {
			anchor_atomno = get_anchor_atomno( anchor_rsd, kinematics::dir_jump );
		}
	} else if ( edge.is_peptide() ) {
		// peptide edge
		int const dir( edge.polymer_direction() );
		debug_assert( anchor_rsd.seqpos() == seqpos-dir );
		anchor_atomno = get_anchor_atomno( anchor_rsd, dir );
	} else {
		// chemical edge
		debug_assert( seqpos == static_cast< Size >( edge.stop() ) && anchor_rsd.seqpos() == static_cast< Size >( edge.start() ) && edge.has_atom_info() );
		anchor_atomno = anchor_rsd.atom_index( edge.upstream_atom() );
	}
	debug_assert( anchor_atomno );
	return anchor_atomno;
}

///////////////////////////////////////////////////////////////////////////////
int
get_anchor_atomno(
	conformation::Residue const & rsd,
	int const dir // forward(1), backward(-1), or "dir_jump"
)
{
	// if dir == 1 or -1, we are building a contiguous edge ("peptide edge")
	// in the fold_tree ( ==> residue should be a polymer type? )
	int const forward( 1 ), backward( -1 );

	if ( dir == forward ) {
		// C for proteins, O3' for DNA
		if ( rsd.is_polymer() ) {
			if ( rsd.is_upper_terminus() ) {
				// this is a little strange to be folding through a terminal residue but kinematics doesn't
				// have to correlate perfectly with chemical
				return rsd.mainchain_atom( rsd.mainchain_atoms().size() );
			} else {
				return rsd.upper_connect_atom();
			}
		} else {
			return get_anchor_atomno( rsd, kinematics::dir_jump );
		}
	} else if ( dir == backward ) {
		// N for proteins, P for DNA
		if ( rsd.is_polymer() ) {
			if ( rsd.is_lower_terminus() ) {
				// this is a little strange to be folding through a terminal residue but kinematics doesn't
				// have to correlate perfectly with chemical
				return rsd.mainchain_atom(1);
			} else {
				return rsd.lower_connect_atom();
			}
		} else {
			return get_anchor_atomno( rsd, kinematics::dir_jump );
		}

	} else {
		debug_assert( dir == kinematics::dir_jump );
		// default for jumps, use the N-terminal attachment?
		if ( rsd.mainchain_atoms().empty() ) {
			// SHORT TERM HACK -- need to add some logic for root atomno
			return 1;
		} else {
			// PB 10/29/07 - changing to use an interior mainchain atom
			Size const nbb( rsd.n_mainchain_atoms() ); //    1  2  3  4  5  6  -- nbb
			Size const mainchain_index( ( nbb-1 )/2 + 1 ); //   1  1  2  2  3  3  -- mainchain_index
			return rsd.mainchain_atoms()[ mainchain_index ];
			//return rsd.mainchain_atoms()[1];
		}
	}
}

/// @details  Determines the anchor and root atom indices based on residue types and the foldtree edge.

void
get_anchor_and_root_atoms(
	conformation::Residue const & anchor_rsd,
	conformation::Residue const & root_rsd,
	kinematics::Edge const & edge,
	Size & anchor_atomno,
	Size & root_atomno
)
{
	ASSERT_ONLY(Size const anchor_pos( anchor_rsd.seqpos() ););
	ASSERT_ONLY(Size const root_pos( root_rsd.seqpos() ););

	if ( edge.is_polymer() ) {
		// POLYMER EDGE
		int const dir( edge.polymer_direction() );
		debug_assert( dir == 1 || dir == -1 );
		debug_assert( root_pos == anchor_pos + dir );
		anchor_atomno = get_anchor_atomno( anchor_rsd, dir );
		root_atomno   = get_root_atomno  (   root_rsd, dir );
	} else {
		debug_assert( anchor_pos == Size(edge.start()) && root_pos == Size(edge.stop()) );
		if ( edge.has_atom_info() ) {
			// JUMP OR CHEMICAL W/ ATOM INFO
			anchor_atomno = anchor_rsd.atom_index( edge.upstream_atom() );
			root_atomno   =   root_rsd.atom_index( edge.downstream_atom() );
		} else {
			if ( edge.is_jump() ) {
				// JUMP EDGE
				anchor_atomno = get_anchor_atomno( anchor_rsd, kinematics::dir_jump );
				root_atomno   =   get_root_atomno(   root_rsd, kinematics::dir_jump );
			} else if ( edge.is_chemical_bond() ) {
				// CHEMICAL EDGE
				get_chemical_root_and_anchor_atomnos( anchor_rsd, root_rsd, anchor_atomno, root_atomno );
			} else {
				utility_exit_with_message( "Unrecognized edge type!" );
			}
		}
	}
	return;
}


// this code assumed we were also being passed old_rsd:
// optionally debug some things before making atom_tree call
/**if ( debug ) {
BondID old_rsd_in;
vector1< BondID > old_rsd_out;
get_residue_connections( old_rsd, fold_tree, residues, old_rsd_in, old_rsd_out );
debug_assert( new_rsd_in.atom1 == old_rsd_in.atom1 && new_rsd_out.size() == old_rsd_out.size() );
for ( Size i=1; i<= new_rsd_out.size(); ++i ) {
debug_assert( new_rsd_out[i].atom2 == old_rsd_out[i].atom2 );
}
}**/


///////////////////////////////////////////////////////////////////////////////

void
promote_sameresidue_child_of_jump_atom(
	kinematics::Edge const & edge,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomTree & atom_tree
)
{
	debug_assert( edge.is_jump() );
	Size root_pos( edge.stop() ), anchor_atomno, root_atomno;
	get_anchor_and_root_atoms( *residues[ edge.start() ], *residues[ root_pos ], edge, anchor_atomno, root_atomno );
	atom_tree.promote_sameresidue_nonjump_child( id::AtomID( root_atomno, root_pos ) );
}


void
promote_sameresidue_child_of_jump_atom(
	kinematics::Edge const & edge,
	conformation::ResidueCOPs const & residues,
	kinematics::AtomPointer2D const & atom_pointer
)
{
	debug_assert( edge.is_jump() );
	Size root_pos( edge.stop() ), anchor_atomno, root_atomno;
	get_anchor_and_root_atoms( *residues[ edge.start() ], *residues[ root_pos ], edge, anchor_atomno, root_atomno );
	kinematics::tree::AtomOP root_atom( atom_pointer[ id::AtomID( root_atomno, root_pos ) ] );
	debug_assert( root_atom->is_jump() );
	kinematics::tree::AtomOP same_residue_child( 0 );
	for ( Size i=0; i< root_atom->n_nonjump_children(); ++i ) {
		kinematics::tree::AtomOP child( atom_pointer[ root_atom->get_nonjump_atom( i )->id() ] ); // want nonconst, use atom_pointer
		if ( Size(child->id().rsd()) == root_pos ) {
			same_residue_child = child;
			break;
		}
	}
	if ( same_residue_child ) {
		debug_assert( !same_residue_child->is_jump() );
		root_atom->delete_atom( same_residue_child );
		root_atom->insert_atom( same_residue_child );
	} else {
		if ( TR.Warning.visible() ) TR.Warning << "Unable to keep stub in residue: jump_atom has no non-jump, same-residue children!" << std::endl;
	}
}

void
get_chemical_root_and_anchor_atomnos(
	conformation::Residue const & rsd_anchor,
	conformation::Residue const & rsd_root,
	Size & anchor_atom_no,
	Size & root_atom_no
)
{
	debug_assert( rsd_anchor.is_bonded( rsd_root ) );

	Size first_connection_id = rsd_anchor.connections_to_residue( rsd_root )[ 1 ];
	chemical::ResConnID first_connection = rsd_anchor.connect_map( first_connection_id );
	anchor_atom_no = rsd_anchor.residue_connection( first_connection_id ).atomno();
	root_atom_no = rsd_root.residue_connection( first_connection.connid() ).atomno();
}

///////////////////////////////////////////////////////////////////////////////
//
// right now, this is used by the pose to setup a mapping which is
// passed in to the atomtree when replacing one residue with another
//
//
// maybe the behavior should really be to identify atoms with the same
// name in the two residues and map them...
//
// oh well, this will work for rewiring backbone connections properly;
// basically the atomtree will use this mapping when it is replacing
// the atoms of rsd1 with the subtree formed by the atoms of rsd2
// for this to work, all outgoing connections from the rsd1 atoms have
// to have their rsd1-atom represented in this map, so the atomtree
// knows where to attach the child when rsd2 is put in.


/// @details only map by atom number, not by identity currently.
void
setup_corresponding_atoms(
	id::AtomID_Map< id::AtomID > & atom_map,
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
)
{

	// PB -- this whole function is going to go away soon, so I'm not really worried about all the hacks

	int const seqpos1( rsd1.seqpos() );
	int const seqpos2( rsd2.seqpos() );

	utility::vector1< Size > const &
		mainchain1( rsd1.mainchain_atoms() ),
		mainchain2( rsd2.mainchain_atoms() );

	debug_assert( mainchain1.size() == mainchain2.size() );

	// special case for virtual residues -- fpd
	if ( mainchain1.size() == 0 &&
			rsd1.type().name3() == "XXX" && rsd2.type().name3() == "XXX" ) {
		debug_assert( rsd1.atoms().size() == rsd2.atoms().size() );

		for ( Size i=1, i_end = rsd1.atoms().size(); i<= i_end; ++i ) {
			atom_map.set( id::AtomID( i, seqpos1 ), id::AtomID( i, seqpos2 ) );
		}
	} else {
		for ( Size i=1, i_end = mainchain1.size(); i<= i_end; ++i ) {
			atom_map.set( id::AtomID( mainchain1[i], seqpos1 ),
				id::AtomID( mainchain2[i], seqpos2 ) );
		}
		if ( rsd1.is_DNA() && rsd2.is_DNA() ) {
			// special case for dna-dna jumps anchored at sidechain atoms
			for ( Size i=1; i<= 4; ++i ) {
				atom_map.set( id::AtomID( rsd1.chi_atoms(1)[i], seqpos1 ),
					id::AtomID( rsd2.chi_atoms(1)[i], seqpos2 ) );
			}
		}

		if ( rsd1.is_RNA() && rsd2.is_RNA() ) {
			// By the way, this is a total hack AND SHOULD NOT BE CHECKED IN
			// BEFORE CONSULTING WITH PHIL! -- rhiju
			// special case for rna-rna jumps anchored at sidechain atoms
			if ( rsd1.name1() == rsd2.name1() ) {
				for ( Size i=1; i<= rsd1.natoms(); ++i ) {
					atom_map.set( id::AtomID( i, seqpos1 ),
						id::AtomID( i, seqpos2 ) );
				}
			}
		}

	}

	/// Now include atoms involved in residue connections, since these might be anchor points for the atomtree
	for ( Size connid=1; connid<= rsd1.n_residue_connections() && connid <= rsd2.n_residue_connections(); ++connid ) {
		Size const atom1( rsd1.type().residue_connection( connid ).atomno() );
		Size const atom2( rsd2.type().residue_connection( connid ).atomno() );
		if ( rsd1.atom_name( atom1 ) == rsd2.atom_name( atom2 ) ) {
			atom_map.set( id::AtomID( atom1, seqpos1 ), id::AtomID( atom2, seqpos2 ) );
		}
	}
}

/// @brief Switch the disulfide state of a disulfide-forming residue (e.g. CYS->CYD or CYD->CYS or
/// DCYD->DCYS or DCYS->DCYD or whatnot).
/// @param[in] index Position of the residue to replace.
/// @param[in] cys_type_name3 The 3-letter name of the cys type to use: e.g. CYS
///  or CYD.  DEPRECATED and kept only for backward-compatibility.
/// @param[inout] conf The conformation to modify
/// @details Substitutes a residue with the given cys type, keeping as many of
///  the existing atom positions as possible.  If the original residue has a
///  disulfide variant it will be removed, otherwise a disulfide variant will
///  be added.  Should work with any ResidueTypeSet that has the proper
///  disulfide variants.  If the replacement fails for any reason a warning
///  will be printed.
/// @return true if the replacement was successful, false otherwise.
bool change_cys_state(
	Size const index,
	std::string const & cys_type_name3, //Deprecated.
	Conformation & conf
) {

	bool removing(false);

	// Cache information on old residue.
	Residue const & res( conf.residue( index ) );
	chemical::ResidueTypeSet const & residue_type_set = res.type().residue_type_set();

	// make sure we're working on a cys
	if ( ( ! res.type().is_sidechain_thiol() ) && ( ! res.type().is_disulfide_bonded() ) ) {
		if ( TR.Warning.visible() ) TR.Warning << "WARNING: change_cys_state() was called on non-cys-like residue " << index << ", skipping!" << std::endl;
		return false;
	}

	// Track the variant types of the old residue type.  We want the
	// new residue to have the same variant type as the old.
	utility::vector1< std::string > variant_types = res.type().properties().get_list_of_variants();

	// check and handle disulfide state
	if ( res.has_variant_type( chemical::DISULFIDE ) ) {
		// if the old residue has DISULFIDE variant type then we are removing a
		// disulfide, so remove the variant type from the list
		variant_types.erase( std::find( variant_types.begin(), variant_types.end(), "DISULFIDE" ) );
		removing = true;
		//TR << "We just erased the DISULFIDE variant type as desired, since we are going from rn " << rn << " to thiol type " << std::endl;

	} else /*If it does not have the DISULFIDE variant type*/ {
		//TR << "We just added the DISULFIDE variant type as desired, since we are going from rn " << rn << " to disulfide type " << std::endl;
		variant_types.push_back( "DISULFIDE" ); //chemical::DISULFIDE ); // creating a disulfide
		removing = false;
	}

	// The following is for backward-compatibility only, for functions that would request CYD when the residue was already CYD or whatnot.
	if ( cys_type_name3=="CYS" && !removing /*not removing if disulfide variant not found*/ ) return true; //If we're asked for a cys and we already have a cys, do nothing.
	else if ( cys_type_name3=="CYD" && removing /*removing if the disulfide variant was found*/ ) return true; //Similarly, if we're asked for a cyd and we already have a cyd, do nothing.

	// Get the residue type of the desired new residue type.
	chemical::ResidueTypeCOP replacement_type( residue_type_set.get_representative_type_name3( res.type().name3(), variant_types ) );
	if ( replacement_type ) {
		ResidueOP new_res = ResidueFactory::create_residue( *replacement_type, res, conf );
		copy_residue_coordinates_and_rebuild_missing_atoms( res, *new_res, conf );
		conf.replace_residue( index, *new_res, false );
		return true;
	}


	// If we are here then a residue type match wasn't found; issue error message.
	if ( TR.Error.visible() ) TR.Error << "ERROR: Couldn't find a " << (removing ? "disulfide-free" : "disulfide-bonded")  << " equivalent for residue " << res.name3() << index << "." <<std::endl;

	return false;
}

id::NamedAtomID
atom_id_to_named_atom_id(
	id::AtomID const & atom_id,
	conformation::Residue const & rsd
){
	std::string atom_name;
	Size rsd_id;
	if ( atom_id.atomno() <= rsd.atoms().size() ) {
		atom_name = rsd.atom_name( atom_id.atomno() );
		rsd_id = atom_id.rsd();
	} else {
		if ( TR.Error.visible() ) TR.Error << "[ ERROR ] Can't resolve atom_id " << atom_id << std::endl;
		rsd_id = 0;
		atom_name ="";
	}

	return id::NamedAtomID( atom_name, rsd_id );
}

id::AtomID
named_atom_id_to_atom_id(
	id::NamedAtomID const & named_atom_id,
	conformation::Residue const & rsd
){
	return id::AtomID( rsd.atom_index( named_atom_id.atom() ), named_atom_id.rsd() );
}

id::NamedStubID
stub_id_to_named_stub_id(
	id::StubID const & stub_id,
	conformation::Residue const & rsd
){
	using namespace core::id;
	if ( stub_id.center().valid() ) {
		return NamedStubID(
			NamedAtomID( atom_id_to_named_atom_id( stub_id.center() , rsd ) ),
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 1 ), rsd ) ),
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 2 ), rsd ) ),
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 3 ), rsd ) )
		);
	} else {
		return NamedStubID( // TODO does this make sense? if input stub is bad, what should be done?
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 1 ), rsd ) ),
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 2 ), rsd ) ),
			NamedAtomID( atom_id_to_named_atom_id( stub_id.atom( 3 ), rsd ) )
		);
	}
}

id::StubID
named_stub_id_to_stub_id(
	id::NamedStubID const & named_stub_id,
	conformation::Residue const & rsd
){
	using namespace core::id;
	if ( named_stub_id.center().valid() ) {
		return StubID(
			AtomID( named_atom_id_to_atom_id( named_stub_id.center() , rsd ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 1 ), rsd ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 2 ), rsd ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 3 ), rsd ) )
		);
	} else {
		return StubID( // TODO does this make sense? if input stub is bad, what should be done?
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 1 ), rsd ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 1 ), rsd ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 2 ), rsd ) ),
			AtomID( named_atom_id_to_atom_id( named_stub_id.atom( 3 ), rsd ) )
		);
	}
}


/// @brief Introduce cysteines at the specified location and define a
///  disulfide bond between them.
/// @details Does not do the repacking & minimization required to place the
///  disulfide correctly.
void
form_disulfide(
	Conformation & conformation,
	Size lower_res,
	Size upper_res,
	bool const preserve_d_residues,
	bool const force_d_residues
) {
	// Verify we're dealing with a FA conformation
	runtime_assert( conformation.is_fullatom() );

	chemical::ResidueTypeSetCOP restype_set =
		chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

	// Break existing disulfide bonds to lower
	if ( conformation.residue( lower_res ).has_variant_type( chemical::DISULFIDE ) ) { // full atom residue
		std::string lcatom = conformation.residue_type(lower_res).get_disulfide_atom_name();
		Size const connect_atom( conformation.residue( lower_res ).atom_index( lcatom ) );
		Size other_res( 0 );
		Size conn(0);
		for ( conn = conformation.residue( lower_res ).type().n_residue_connections(); conn >= 1; --conn ) {
			if ( Size( conformation.residue(lower_res).type().residue_connection(conn).atomno() ) == connect_atom ) {
				other_res = conformation.residue( lower_res ).connect_map( conn ).resid();
				break;
			}
		}
		if ( other_res == 0 ) {
			if ( TR.Error.visible() ) TR.Error << "Error: Residue " << lower_res << " was disulfide bonded but had no partner" << std::endl;
			utility_exit();
		}

		if ( other_res == upper_res ) {
			// Already a disulfide bond
			runtime_assert_msg(conformation.residue( upper_res ).connect_map( conn ).resid() == lower_res,
				"Error: Disulfide bond wasn't reciprocal");
			return;
		}

		// Break the disulfide bond to upper_res
		bool result = change_cys_state( other_res, "", conformation );
		runtime_assert_msg(result,"Error removing disulfide variant from "+conformation.residue(other_res).name3());
	} else {
		form_disulfide_helper(conformation, lower_res, restype_set, preserve_d_residues, force_d_residues); //This mutates the lower_res residue to a disulfide-forming variant, preserving other variant types in the process.
	}
	// Break existing disulfide bonds to upper
	if ( !upper_is_symm_equivalent_of_lower( conformation, lower_res, upper_res ) ) { //If the upper residue is the symmetric equivalent of the lower residue, then we've already altered its disulfide type, and can skip the next few steps.
		if ( conformation.residue( upper_res ).has_variant_type( chemical::DISULFIDE ) ) {
			std::string ucatom = conformation.residue_type(lower_res).get_disulfide_atom_name();
			Size const connect_atom( conformation.residue( upper_res ).atom_index( ucatom ) );
			Size other_res( 0 );
			Size conn(0);
			for ( conn = conformation.residue( upper_res ).type().n_residue_connections(); conn >= 1; --conn ) {
				if ( Size ( conformation.residue(upper_res).type().residue_connection(conn).atomno() ) == connect_atom ) {
					other_res = conformation.residue( upper_res ).connect_map( conn ).resid();
					break;
				}
			}
			if ( other_res == 0 ) {
				if ( TR.Error.visible() ) TR.Error << "Error: Residue " << upper_res << " was disulfide bonded but had no partner" << std::endl;
				utility_exit();
			}

			// Break the disulfide bond to lower_res
			bool result = change_cys_state( other_res, "", conformation );
			runtime_assert_msg(result,"Error removing disulfide variant from "+conformation.residue(other_res).name3());
		} else {
			form_disulfide_helper(conformation, upper_res, restype_set, preserve_d_residues, force_d_residues); //This mutates the upper_res residue to a disulfide-forming variant, preserving other variant types in the process.
		}
	} //If !upper_is_symm_equivalent_of_lower()

	// Both residues are now CYD
	runtime_assert( conformation.residue(lower_res).has_variant_type(chemical::DISULFIDE) );
	runtime_assert( conformation.residue(upper_res).has_variant_type(chemical::DISULFIDE) );

	//form the bond between the two residues
	// NOW we still need to get the lc and uc atoms, because we couldn't get them before
	// since this function doesn't necessarily start with a cys type
	// (damn this function!)
	std::string lcatom = conformation.residue_type(lower_res).get_disulfide_atom_name();
	std::string ucatom = conformation.residue_type(upper_res).get_disulfide_atom_name();
	conformation.declare_chemical_bond(lower_res,lcatom,upper_res,ucatom);

}

/// @brief Helper function for the form_disulfide function.
/// @details This function ensures that as a residue is mutated to a disulfide-bonding residue type,
/// all other variant types are preserved; it is used to avoid code duplication.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)
void form_disulfide_helper(
	core::conformation::Conformation &conformation,
	core::Size const res_index,
	core::chemical::ResidueTypeSetCOP restype_set,
	bool const preserve_d_residues,
	bool const force_d_residues
) {
	using namespace core::chemical;
	AA query_type( aa_cys );

	// amw: This was not explicit enough to me and I briefly altered this function's semantics!
	// If the residues aren't already a disulfide-forming type, we must make them one. Functions like
	// disulfide insertion mover assume we do this.
	if ( force_d_residues ) {
		query_type = aa_dcs;
	} else {
		if ( ! conformation.residue_type(res_index).is_disulfide_bonded() &&
				! conformation.residue_type(res_index).is_sidechain_thiol() &&
				conformation.residue_type(res_index).is_alpha_aa() ) {
			if ( !conformation.residue_type(res_index).is_d_aa() || !preserve_d_residues ) {
				query_type = aa_cys;
			} else {
				query_type = aa_dcs;
			}
		}
	}

	//Get the list of variant types for this residue, so that the replacement has the same variant types.
	utility::vector1< std::string > variant_types = conformation.residue(res_index).type().properties().get_list_of_variants();

	//Add a disulfide variant type to the list, since we're turning CYS into CYD and so forth.
	bool has_disulfide(false);
	for ( core::Size i = 1, imax = variant_types.size(); i <= imax; ++i ) {
		if ( variant_types[ i ] == "DISULFIDE" ) {
			has_disulfide = true;
			break;
		}
	}
	if ( !has_disulfide ) variant_types.push_back( "DISULFIDE" );

	//Get the suitable variant type, and create a residue with it:
	ResidueOP lower_cyd( ResidueFactory::create_residue( *(restype_set->get_representative_type_aa( query_type, variant_types )), conformation.residue( res_index ), conformation ) );

	copy_residue_coordinates_and_rebuild_missing_atoms( conformation.residue(res_index), *lower_cyd, conformation, true );
	conformation.replace_residue(res_index, *lower_cyd, false /*backbone already oriented*/); // doug

	return;
}

/// @brief Another helper function for the form_disulfide function.
/// @details Returns true if and only if the conformation is symmetric and upper_res is a symmetric copy of lower_res.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)
bool
upper_is_symm_equivalent_of_lower(
	core::conformation::Conformation const &conformation,
	core::Size const lower_res,
	core::Size const upper_res
) {
	core::conformation::symmetry::SymmetricConformationCOP sym_conf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const >(conformation.get_self_ptr()) );
	if ( !sym_conf ) return false;
	//From now on, can assume conformation IS symmetric.

	if ( sym_conf->Symmetry_Info()->bb_is_independent( lower_res ) ) {
		if ( sym_conf->Symmetry_Info()->bb_is_independent( upper_res ) ) return false; //If both residues are independent, they're not equivalent.
		if ( sym_conf->Symmetry_Info()->bb_follows( upper_res ) == lower_res ) return true;  //If the upper res follows the lower, then these ARE equivalent residues.
		return false; //Otherwise they're not, if the lower res is independent.
	}
	//From now on, assume lower_res depends on something
	if ( sym_conf->Symmetry_Info()->bb_is_independent( upper_res ) ) {
		if ( sym_conf->Symmetry_Info()->bb_follows(lower_res) == upper_res ) return true; //If the lower res follows the upper, then these ARE equivalent residues.
		return false; //Otherwise they're not, if the upper is independent.
	}

	//From now on, lower_res and upper_res can be assumed to depend on something.  If they depend on the same thing, they're equivalent.
	if ( sym_conf->Symmetry_Info()->bb_follows(lower_res) == sym_conf->Symmetry_Info()->bb_follows(upper_res) ) return true;

	return false;
}

/// @brief Find whether there is a disulfide defined between two residues
///
/// @details We define a disulfide to exist between a pair of residues iff
///  -# They are both cysteines
///  -# They are bonded by their sidechains
bool
is_disulfide_bond( conformation::Conformation const& conformation, Size residueA_pos, Size residueB_pos)
{
	Residue const& A = conformation.residue(residueA_pos);
	Residue const& B = conformation.residue(residueB_pos);

	if ( !A.is_protein() || !B.is_protein() ) {
		return false;
	}

	//both Cys or CysD
	//if( A.type().name1() != 'C' || B.type().name1() != 'C' )
	if ( ( ! A.type().is_sidechain_thiol() && ! A.type().is_disulfide_bonded() )
			|| ( ! B.type().is_sidechain_thiol() && ! B.type().is_disulfide_bonded() ) ) {
		return false;
	}

	//std::string n1 = A.type().name3(); std::string n2 = B.type().name3();

	//bonded
	runtime_assert( A.type().has( A.type().get_disulfide_atom_name() ) ); //should be fa or centroid
	Size a_connect_atom = A.atom_index( A.type().get_disulfide_atom_name() );

	for ( Size connection = A.type().n_residue_connections(); connection >= 1; --connection ) {
		//check if A bonded to B
		if ( (Size) A.type().residue_connection( connection ).atomno() == a_connect_atom && //bond to sg, not the backbone
				A.connect_map( connection ).resid() == residueB_pos ) { //bonded to B
			return true;
		}
	}

	return false;
}

/// @brief Generate a list of all disulfide bonds in the conformation
void
disulfide_bonds( conformation::Conformation const& conformation, utility::vector1< std::pair<Size,Size> > & disulfides )
{
	for ( Size i = 1; i<= conformation.size(); ++i ) {
		Residue const& res(conformation.residue(i));

		// Skip things besides CYD
		//if( !(res.aa() == chemical::aa_cys && res.has_variant_type(chemical::DISULFIDE) ))
		if ( !(res.type().is_disulfide_bonded() && res.has_variant_type(chemical::DISULFIDE) ) ) {
			continue;
		}
		Size connect_atom( 0);
		if ( res.type().get_disulfide_atom_name() == "NONE" ) {
			if ( TR.Warning.visible() ) TR.Warning << "Warning: unable to establish which atom to use for the disulfide to residue " << i << std::endl;
			continue;
		} else {
			connect_atom = res.atom_index( res.type().get_disulfide_atom_name() ) ;
		}

		Size other_res(0);
		Size conn;
		for ( conn = conformation.residue( i ).type().n_residue_connections(); conn >= 1; --conn ) {
			if ( Size( conformation.residue( i ).type().residue_connection(conn).atomno() ) == connect_atom ) {
				other_res = conformation.residue( i ).connect_map( conn ).resid();
				break;
			}
		}
		if ( other_res == 0 ) {
			TR.Error << "Error: Residue " << i << " was disulfide bonded but had no partner" << std::endl;
			utility_exit();
		}

		// Output the pair once
		if ( i < other_res ) {
			disulfides.push_back( std::make_pair(i, other_res) );
		}
	}
}


// Is the query atom in this residue axial or equatorial to the given ring or neither?
/// @details This function calculates an average plane and determines whether the coordinates of a given atom are
/// axial or equatorial to it (or neither).  The attachment atom is auto-detected.
/// @param   <residue>:    The Residue containing the atoms in question.
/// @param   <query_atom>: The index of the atom in question.
/// @param   <ring_atoms>: A list of indices for the atoms of a monocyclic ring system in sequence.
/// @return  An AxEqDesignation enum type value: AXIAL, EQUATORIAL, or NEITHER
/// @author  Labonte <JWLabonte@jhu.edu>
chemical::rings::AxEqDesignation
is_atom_axial_or_equatorial_to_ring(
	Residue const & residue,
	uint query_atom,
	utility::vector1< uint > const & ring_atoms )
{
	using namespace utility;
	using namespace numeric;

	vector1< uint > const bonded_heavy_atoms( residue.get_adjacent_heavy_atoms( query_atom ) );
	Size const n_bonded_heavy_atoms( bonded_heavy_atoms.size() );
	for ( uint i( 1 ); i <= n_bonded_heavy_atoms; ++i ) {
		if ( ring_atoms.contains( bonded_heavy_atoms[ i ] ) ) {
			// We found an attachment point.
			// Now we need to get the coordinates of all the important atoms and call the function that does the actual
			// math.
			xyzVector< Distance > const query_atom_coords( residue.xyz( query_atom ) );
			xyzVector< Distance > const attachment_atom_coords( residue.xyz( bonded_heavy_atoms[ i ] ) );

			Size const n_ring_atoms( ring_atoms.size() );
			vector1< xyzVector< Distance > > ring_atom_coords( n_ring_atoms );
			for ( uint j( 1 ); j <= n_ring_atoms; ++j ) {
				ring_atom_coords[ j ] = residue.xyz( ring_atoms[ j ] );
			}
			return chemical::rings::is_atom_axial_or_equatorial_to_ring(
				query_atom_coords, attachment_atom_coords, ring_atom_coords );
		}
	}
	TR.Warning << "The attachment point for the query atom is not found in the ring; ";
	TR.Warning << "an axial/equatorial designation is meaningless." << std::endl;
	return chemical::rings::NEITHER;
}

// Is the query atom in this residue axial or equatorial or neither?
/// @details This function calculates an average plane and determines whether the coordinates of a given atom are
/// axial or equatorial to it (or neither).  The ring is requested from the Residue.
/// @param   <residue>:    The Residue containing the atoms in question.
/// @param   <query_atom>: The index of the atom in question.
/// @return  An AxEqDesignation enum type value: AXIAL, EQUATORIAL, or NEITHER
/// @author  Labonte <JWLabonte@jhu.edu>
/*chemical::rings::AxEqDesignation
is_atom_axial_or_equatorial( Residue const & residue, uint query_atom )
{
chemical::ResidueType const & rsd_type( residue.type() );
if ( ! rsd_type.is_cyclic() ) {
TR.Warning << "Queried atom must be on a cyclic residue." << std::endl;
return chemical::rings::NEITHER;
}



return is_atom_axial_or_equatorial_to_ring( residue, query_atom, rsd_type.ring_atoms() );
}*/

} // namespace conformation
} // namespace core
