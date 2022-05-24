// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/util.cc
/// @brief Utility functions for setting up thioether cyclization.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/thioether_util.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/id/DOF_ID.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/CharmmPeriodicFunc.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>

// Protocols headers
#include <protocols/simple_moves/DeclareBond.hh>

// Numeric headers
#include <numeric/constants.hh>

// Basic headers
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.cyclic_peptide.crosslinker.util" );


namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

/// @brief Given a pose and two residues, set up the thioether variant types.
/// @details Sidechainres gets SIDECHAIN_CONJUGATION; ntermres gets ACETYLATED_NTERMINUS_CONNECTION_VARIANT.
/// The pose is modified by this operation.
void
set_up_thioether_variants(
	core::pose::Pose & pose,
	core::Size const ntermres,
	core::Size const sidechainres
) {
	std::string const errmsg( "Error in protocols::cyclic_peptide::crosslinker::set_up_thioether_variants(): " );
	runtime_assert_string_msg(
		pose.residue_type(sidechainres).is_sidechain_thiol(),
		errmsg + "Position " + std::to_string( sidechainres )
		+ " was selected, but this is not a residue with a thiol-containing sidechain."
	);

	runtime_assert_string_msg(
		pose.residue(ntermres).connected_residue_at_lower() == 0,
		errmsg + "Position " + std::to_string( ntermres ) + " was selected as the N-terminal residue for a "
		"thioether connection, but this residue has another residue connected at its lower terminus."
	);

	if ( pose.residue_type(ntermres).is_lower_terminus() ) {
		core::pose::remove_lower_terminus_type_from_pose_residue( pose, ntermres );
		TR << "Removed lower terminus type from residue " << pose.residue_type(ntermres).base_name() << ntermres << "." << std::endl;
	}

	core::pose::add_variant_type_to_pose_residue( pose, core::chemical::SIDECHAIN_CONJUGATION, sidechainres );
	TR << "Added sidechain conjugation type to residue " << pose.residue_type(sidechainres).base_name() << sidechainres << "." << std::endl;
	core::pose::add_variant_type_to_pose_residue( pose, core::chemical::ACETYLATED_NTERMINUS_CONNECTION_VARIANT, ntermres );
	TR << "Added acetylated N-terminus connection type to residue " << pose.residue_type(ntermres).base_name() << ntermres << "." << std::endl;

	pose.update_residue_neighbors();
	// This connects to ii+1 by default, but that's appropriate -- it just got this connection harmed
	// by adding the patch anyway.
	pose.conformation().update_polymeric_connection( ntermres );

	// Force rebuild of the "H":
	if ( pose.residue_type(ntermres).has("H") ) {
		pose.set_xyz( core::id::NamedAtomID( "H", ntermres ), pose.residue_type( ntermres ).icoor( pose.residue_type( ntermres ).atom_index("H") ).build( pose.residue( ntermres ), pose.conformation() ) );
	}

	pose.update_residue_neighbors();
}

/// @brief Set up the mover that creates thioether lariat bonds.
void
set_up_thioether_bond_mover(
	protocols::simple_moves::DeclareBond & termini,
	core::pose::Pose const & pose,
	core::Size const ntermres,
	core::Size const sidechainres
) {
	//Get the name of the first sidechain connection atom:
	core::chemical::ResidueType const & restype( pose.residue_type( ntermres ) );
	core::Size const ntermres_sc_connection_id( restype.n_possible_residue_connections() ); //We will assume that the highest-numbered residue connection is the sidechain connection ID.
	core::Size const ntermres_sc_connection_atom_index( restype.residue_connect_atom_index( ntermres_sc_connection_id ) );
	std::string const ntermres_sc_connection_atom( restype.atom_name( ntermres_sc_connection_atom_index ) );

	//Get the name of the second sidechain connection atom:
	core::chemical::ResidueType const & restype2( pose.residue_type( sidechainres ) );
	core::Size const sidechainres_sc_connection_id( restype2.n_possible_residue_connections() ); //We will assume that the highest-numbered residue connection is the sidechain connection ID.
	core::Size const sidechainres_sc_connection_atom_index( restype2.residue_connect_atom_index( sidechainres_sc_connection_id ) );
	std::string const sidechainres_sc_connection_atom( restype2.atom_name( sidechainres_sc_connection_atom_index ) );

	TR << "Setting up thioether lariat covalent bond from " << restype.base_name() << ntermres << ", atom " << ntermres_sc_connection_atom << " to residue " << restype2.base_name() << sidechainres << ", atom " << sidechainres_sc_connection_atom << "." << std::endl;

	termini.set( sidechainres, sidechainres_sc_connection_atom, ntermres, ntermres_sc_connection_atom, false );
}

/// @brief Correct the bond angles and bond lenghts for virtual atoms at thioether bonds.
void
correct_thioether_virtuals(
	core::pose::Pose & pose,
	core::Size const ntermres,
	core::Size const sidechainres
) {
	{ //Update ideal bondlength and bond angle on sidechain.
		core::id::AtomID const atom_virt_sidechainres(
			pose.residue_type(sidechainres).atom_index("V1"),
			sidechainres
		);
		core::id::AtomID const atom_virt_sidechainres_parent(
			pose.residue_type(sidechainres).icoor( atom_virt_sidechainres.atomno() ).stub_atom1().atomno(),
			sidechainres
		);
		core::id::AtomID const atom_virt_sidechainres_grandparent(
			pose.residue_type(sidechainres).icoor( atom_virt_sidechainres.atomno() ).stub_atom2().atomno(),
			sidechainres
		);

		pose.conformation().set_bond_length( atom_virt_sidechainres_parent, atom_virt_sidechainres, THIOETHER_UTIL_THIOETHER_BOND_LENGTH );
		pose.conformation().set_bond_angle( atom_virt_sidechainres_grandparent, atom_virt_sidechainres_parent, atom_virt_sidechainres, THIOETHER_UTIL_THIOETHER_BOND_C_ANGLE );
	}

	{ //Update ideal bondlength and bond angle on N-terminus.
		core::id::AtomID const atom_virt_ntermres(
			pose.residue_type(ntermres).atom_index("VTH"),
			ntermres
		);
		core::id::AtomID const atom_virt_ntermres_parent(
			pose.residue_type(ntermres).icoor( atom_virt_ntermres.atomno() ).stub_atom1().atomno(),
			ntermres
		);
		core::id::AtomID const atom_virt_ntermres_grandparent(
			pose.residue_type(ntermres).icoor( atom_virt_ntermres.atomno() ).stub_atom2().atomno(),
			ntermres
		);

		pose.conformation().set_bond_length( atom_virt_ntermres_parent, atom_virt_ntermres, THIOETHER_UTIL_THIOETHER_BOND_LENGTH );
		pose.conformation().set_bond_angle( atom_virt_ntermres_grandparent, atom_virt_ntermres_parent, atom_virt_ntermres, THIOETHER_UTIL_THIOETHER_BOND_N_ANGLE );
	}

	pose.update_residue_neighbors();
}

/// @brief Given a pose and two residues to constrain (the first being the residue with the modified N-terminus,
/// and the second being the one with the thiol-containing sidechain), add constraints for a thioether linkage.
void
set_up_thioether_constraints(
	core::pose::Pose & pose,
	core::Size const ntermres,
	core::Size const sidechainres
) {
	// The central bond that's not part of the scoring function is from the
	// sidechainres cysteine SG to the ntermres acetyl CP2.
	using namespace core::pose;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;
	using namespace core::id;

	TR << "Setting up thioether lariat constraints." << std::endl;

	//The residue types of the N and C residues:
	core::chemical::ResidueType const & ntype( pose.residue_type(ntermres) );
	core::chemical::ResidueType const & ctype( pose.residue_type(sidechainres) );

	//The four atoms defining the thioether bond:
	AtomID const atom_a(
		ctype.icoor(ctype.residue_connect_atom_index(ctype.n_possible_residue_connections())).stub_atom1().atomno(),
		sidechainres );
	AtomID const atom_b(
		ctype.residue_connect_atom_index( ctype.n_possible_residue_connections() ),
		sidechainres );
	AtomID const atom_c(
		ntype.residue_connect_atom_index( ntype.n_possible_residue_connections() ),
		ntermres );
	AtomID const atom_d(
		ntype.icoor( ntype.residue_connect_atom_index( ntype.n_possible_residue_connections() ) ).stub_atom1().atomno(),
		ntermres );

	AtomID const atom_virt_sidechainres(
		pose.residue_type(sidechainres).atom_index("V1"),
		sidechainres
	);
	AtomID const atom_virt_ntermres(
		pose.residue_type(ntermres).atom_index("VTH"),
		ntermres
	);

	TR << "The following four atoms define the thioether bond:" << std::endl;
	TR << "1.\tRes=" << atom_a.rsd() << "\tAtom=" << pose.residue(atom_a.rsd()).atom_name(atom_a.atomno()) << std::endl;
	TR << "2.\tRes=" << atom_b.rsd() << "\tAtom=" << pose.residue(atom_b.rsd()).atom_name(atom_b.atomno()) << std::endl;
	TR << "3.\tRes=" << atom_c.rsd() << "\tAtom=" << pose.residue(atom_c.rsd()).atom_name(atom_c.atomno()) << std::endl;
	TR << "4.\tRes=" << atom_d.rsd() << "\tAtom=" << pose.residue(atom_d.rsd()).atom_name(atom_d.atomno()) << std::endl;

	{//Thioether bond length constraint:
		FuncOP harmfunc1( utility::pointer::make_shared< HarmonicFunc >( THIOETHER_UTIL_THIOETHER_BOND_LENGTH, 0.01) );
		ConstraintCOP distconst1( utility::pointer::make_shared< AtomPairConstraint >( atom_b, atom_c, harmfunc1 ) );
		pose.add_constraint (distconst1);
	}

	{//Thioether virtual atom constraints:
		FuncOP harmfunc2( utility::pointer::make_shared< HarmonicFunc >( 0.0, 0.01) );
		FuncOP harmfunc3( utility::pointer::make_shared< HarmonicFunc >( 0.0, 0.01) );
		ConstraintCOP distconst2( utility::pointer::make_shared< AtomPairConstraint >( atom_virt_ntermres, atom_b, harmfunc2 ) );
		ConstraintCOP distconst3( utility::pointer::make_shared< AtomPairConstraint >( atom_virt_sidechainres, atom_c, harmfunc3 ) );
		pose.add_constraint (distconst2);
		pose.add_constraint (distconst3);
	}

	{ //Thioether dihedral angle constraints, from methionine
		// (TODO -- change these if we sample a trans-proline.)
		FuncOP circharmfunc1( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi, 0.24, 1 ) );
		ConstraintCOP dihedconst1( utility::pointer::make_shared< DihedralConstraint >( atom_a, atom_b, atom_c, atom_d, circharmfunc1) );
		pose.add_constraint (dihedconst1);
		FuncOP circharmfunc2( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 0.37, 3 ) );
		ConstraintCOP dihedconst2( utility::pointer::make_shared< DihedralConstraint >( atom_a, atom_b, atom_c, atom_d, circharmfunc2) );
		pose.add_constraint (dihedconst2);
	}

	{ //Thioether bond angle constraints:
		FuncOP circharmfunc2a( utility::pointer::make_shared< CircularHarmonicFunc >( THIOETHER_UTIL_THIOETHER_BOND_C_ANGLE, 0.02) );
		FuncOP circharmfunc2b( utility::pointer::make_shared< CircularHarmonicFunc >( THIOETHER_UTIL_THIOETHER_BOND_N_ANGLE, 0.02) );
		ConstraintCOP angleconst1( utility::pointer::make_shared< AngleConstraint >( atom_a, atom_b, atom_c, circharmfunc2a) );
		ConstraintCOP angleconst2( utility::pointer::make_shared< AngleConstraint >( atom_b, atom_c, atom_d, circharmfunc2b) );
		pose.add_constraint (angleconst1);
		pose.add_constraint (angleconst2);
	}

	TR << "Finished setting up constraints." << std::endl;
}

} //crosslinker
} //cyclic_peptide
} //protocols

