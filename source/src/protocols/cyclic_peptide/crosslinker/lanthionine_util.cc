// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/util.cc
/// @brief Utility functions for setting up lanthionine cyclization.
/// @author Clay Tydings (claiborne.w.tydings@vanderbilt.edu)

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/lanthionine_util.hh>

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

/// @brief Given a pose and two residues, set up the lanthionine variant types.
/// @details Sidechainres gets SIDECHAIN_CONJUGATION; dalares gets ACETYLATED_NTERMINUS_CONNECTION_VARIANT.
/// The pose is modified by this operation.
void
set_up_lanthionine_variants(
	core::pose::Pose & pose,
	core::Size const dalares,
	core::Size const cysres
) {
	std::string const errmsg( "Error in protocols::cyclic_peptide::crosslinker::set_up_lanthionine_variants(): " );
	runtime_assert_string_msg(
		pose.residue_type(cysres).is_sidechain_thiol(),
		errmsg + "Position " + std::to_string( cysres )
		+ " was selected, but this is not a residue with a thiol-containing sidechain."
	);

	runtime_assert_string_msg(
		//pose.residue(dalares).connected_residue_at_lower() == 0,
		//pose.residue_type(dalares).is_d_aa(),
		pose.residue_type( dalares ).base_name() == "DBS"  || pose.residue_type( dalares ).base_name() == "DBR"
		|| pose.residue_type( dalares ).base_name() == "DDBS"  || pose.residue_type( dalares ).base_name() == "DDBR"
		|| pose.residue_type(dalares).base_name() == "DALA" || pose.residue_type(dalares).base_name() == "ALA",
		errmsg + "Position " + std::to_string( dalares ) + " was selected as the non-cys residue for a "
		"lanthionine connection, but this residue is not a D-ALA or ALA or Butyrine Variant."
	);

	//if ( pose.residue_type(dalares).is_lower_terminus() ) {
	// core::pose::remove_lower_terminus_type_from_pose_residue( pose, dalares );
	// TR << "Removed lower terminus type from residue " << pose.residue_type(dalares).base_name() << dalares << "." << std::endl;
	//}

	core::pose::add_variant_type_to_pose_residue( pose, core::chemical::SIDECHAIN_CONJUGATION, cysres );
	TR << "Added sidechain conjugation type to residue " << pose.residue_type(cysres).base_name() << cysres << "." << std::endl;

	core::pose::add_variant_type_to_pose_residue( pose, core::chemical::SC_BRANCH_POINT, dalares );
	TR << "Added acetylated SC_BRANCH_POINT connection type to residue " << pose.residue_type(dalares).base_name() << dalares << "." << std::endl;

	pose.update_residue_neighbors();
	// This connects to ii+1 by default, but that's appropriate -- it just got this connection harmed
	// by adding the patch anyway.
	//pose.conformation().update_polymeric_connection( dalares );

	// Force rebuild of the "H":
	//if ( pose.residue_type(dalares).has("H") ) {
	// pose.set_xyz( core::id::NamedAtomID( "H", dalares ), pose.residue_type( dalares ).icoor( pose.residue_type( dalares ).atom_index("H") ).build( pose.residue( dalares ), pose.conformation() ) );
	//}

	//pose.update_residue_neighbors();
}

/// @brief Set up the mover that creates lanthionine lariat bonds.
void
set_up_lanthionine_bond_mover(
	protocols::simple_moves::DeclareBond & termini,
	core::pose::Pose const & pose,
	core::Size const dalares,
	core::Size const cysres
) {
	//Get the name of the first sidechain connection atom:
	core::chemical::ResidueType const & restype( pose.residue_type( dalares ) );
	core::Size const dalares_sc_connection_id( restype.n_possible_residue_connections() ); //We will assume that the highest-numbered residue connection is the sidechain connection ID.
	core::Size const dalares_sc_connection_atom_index( restype.residue_connect_atom_index( dalares_sc_connection_id ) );
	std::string const dalares_sc_connection_atom( restype.atom_name( dalares_sc_connection_atom_index ) );

	//Get the name of the second sidechain connection atom:
	core::chemical::ResidueType const & restype2( pose.residue_type( cysres ) );
	core::Size const cysres_sc_connection_id( restype2.n_possible_residue_connections() ); //We will assume that the highest-numbered residue connection is the sidechain connection ID.
	core::Size const cysres_sc_connection_atom_index( restype2.residue_connect_atom_index( cysres_sc_connection_id ) );
	std::string const cysres_sc_connection_atom( restype2.atom_name( cysres_sc_connection_atom_index ) );

	TR << "Setting up lanthionine linker covalent bond from " << restype.base_name() << dalares << ", atom " << dalares_sc_connection_atom << " to residue " << restype2.base_name() << cysres << ", atom " << cysres_sc_connection_atom << "." << std::endl;

	termini.set( cysres, cysres_sc_connection_atom, dalares, dalares_sc_connection_atom, false );
}

/// @brief Correct the bond angles and bond lengths for virtual atoms at lanthionine bonds.
void
correct_lanthionine_virtuals(
	core::pose::Pose & pose,
	core::Size const dalares,
	core::Size const cysres
) {
	{ //Update ideal bondlength and bond angle on sidechain.
		core::id::AtomID const atom_virt_cysres(
			pose.residue_type(cysres).atom_index("V1"),
			cysres
		);
		core::id::AtomID const atom_virt_cysres_parent(
			pose.residue_type(cysres).icoor( atom_virt_cysres.atomno() ).stub_atom1().atomno(),
			cysres
		);
		core::id::AtomID const atom_virt_cysres_grandparent(
			pose.residue_type(cysres).icoor( atom_virt_cysres.atomno() ).stub_atom2().atomno(),
			cysres
		);

		pose.conformation().set_bond_length( atom_virt_cysres_parent, atom_virt_cysres, LANTHIONINE_UTIL_LANTHIONINE_BOND_LENGTH );
		pose.conformation().set_bond_angle( atom_virt_cysres_grandparent, atom_virt_cysres_parent, atom_virt_cysres, LANTHIONINE_UTIL_LANTHIONINE_BOND_S_ANGLE );
	}

	{ //Update ideal bondlength and bond angle on N-terminus.
		core::id::AtomID const atom_virt_dalares(
			pose.residue_type(dalares).atom_index("V1"),
			dalares
		);
		core::id::AtomID const atom_virt_dalares_parent(
			pose.residue_type(dalares).icoor( atom_virt_dalares.atomno() ).stub_atom1().atomno(),
			dalares
		);
		core::id::AtomID const atom_virt_dalares_grandparent(
			pose.residue_type(dalares).icoor( atom_virt_dalares.atomno() ).stub_atom2().atomno(),
			dalares
		);

		pose.conformation().set_bond_length( atom_virt_dalares_parent, atom_virt_dalares, LANTHIONINE_UTIL_LANTHIONINE_BOND_LENGTH );
		pose.conformation().set_bond_angle( atom_virt_dalares_grandparent, atom_virt_dalares_parent, atom_virt_dalares, LANTHIONINE_UTIL_LANTHIONINE_BOND_C_ANGLE );
	}

	pose.update_residue_neighbors();
}

/// @brief Given a pose and two residues to constrain (the first being the residue with the modified N-terminus,
/// and the second being the one with the thiol-containing sidechain), add constraints for a lanthionine linkage.
void
set_up_lanthionine_constraints(
	core::pose::Pose & pose,
	core::Size const dalares,
	core::Size const cysres
) {
	// The central bond that's not part of the scoring function is from the
	// cysres cysteine SG to the dalares acetyl CP2.
	using namespace core::pose;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;
	using namespace core::id;

	TR << "Setting up lanthionine lariat linker." << std::endl;

	//The residue types of the D-ALA and Cys residues:
	core::chemical::ResidueType const & dalatype( pose.residue_type(dalares) );
	core::chemical::ResidueType const & cystype( pose.residue_type(cysres) );

	//The four atoms defining the lanthionine bond:
	//TODO change this to six atoms
	AtomID const atom_ao(
		cystype.icoor(cystype.icoor(cystype.residue_connect_atom_index(cystype.n_possible_residue_connections())).stub_atom1().atomno()).stub_atom1().atomno(),
		cysres );
	AtomID const atom_a(
		cystype.icoor(cystype.residue_connect_atom_index(cystype.n_possible_residue_connections())).stub_atom1().atomno(),
		cysres );
	AtomID const atom_b(
		cystype.residue_connect_atom_index( cystype.n_possible_residue_connections() ),
		cysres );
	AtomID const atom_c(
		dalatype.residue_connect_atom_index( dalatype.n_possible_residue_connections() ),
		dalares );
	AtomID const atom_d(
		dalatype.icoor( dalatype.residue_connect_atom_index( dalatype.n_possible_residue_connections() ) ).stub_atom1().atomno(),
		dalares );
	AtomID const atom_do(
		dalatype.icoor(dalatype.icoor(dalatype.residue_connect_atom_index(dalatype.n_possible_residue_connections())).stub_atom1().atomno()).stub_atom1().atomno(),
		dalares );

	AtomID const atom_virt_cysres(
		pose.residue_type(cysres).atom_index("V1"),
		cysres
	);
	AtomID const atom_virt_dalares(
		pose.residue_type(dalares).atom_index("V1"),
		dalares
	);

	TR << "The following six atoms define the lanthionine bond:" << std::endl;
	TR << "0.\tRes=" << atom_ao.rsd() << "\tAtom=" << pose.residue(atom_ao.rsd()).atom_name(atom_ao.atomno()) << std::endl;
	TR << "1.\tRes=" << atom_a.rsd() << "\tAtom=" << pose.residue(atom_a.rsd()).atom_name(atom_a.atomno()) << std::endl;
	TR << "2.\tRes=" << atom_b.rsd() << "\tAtom=" << pose.residue(atom_b.rsd()).atom_name(atom_b.atomno()) << std::endl;
	TR << "3.\tRes=" << atom_c.rsd() << "\tAtom=" << pose.residue(atom_c.rsd()).atom_name(atom_c.atomno()) << std::endl;
	TR << "4.\tRes=" << atom_d.rsd() << "\tAtom=" << pose.residue(atom_d.rsd()).atom_name(atom_d.atomno()) << std::endl;
	TR << "5.\tRes=" << atom_do.rsd() << "\tAtom=" << pose.residue(atom_do.rsd()).atom_name(atom_do.atomno()) << std::endl;

	{//Lanthionine bond length constraint:
		FuncOP harmfunc1( utility::pointer::make_shared< HarmonicFunc >( LANTHIONINE_UTIL_LANTHIONINE_BOND_LENGTH, 0.16) );
		ConstraintCOP distconst1( utility::pointer::make_shared< AtomPairConstraint >( atom_b, atom_c, harmfunc1 ) );
		pose.add_constraint (distconst1);
	}

	{//Lanthionine virtual atom constraints:
		FuncOP harmfunc2( utility::pointer::make_shared< HarmonicFunc >( 0.0, 0.01) );
		FuncOP harmfunc3( utility::pointer::make_shared< HarmonicFunc >( 0.0, 0.01) );
		ConstraintCOP distconst2( utility::pointer::make_shared< AtomPairConstraint >( atom_virt_dalares, atom_b, harmfunc2 ) );
		ConstraintCOP distconst3( utility::pointer::make_shared< AtomPairConstraint >( atom_virt_cysres, atom_c, harmfunc3 ) );
		pose.add_constraint (distconst2);
		pose.add_constraint (distconst3);
	}

	if ( pose.residue_type(dalares).base_name() == "DALA" || pose.residue_type(dalares).base_name() == "ALA" ) {
		if ( pose.residue_type(
				dalares).is_d_aa() ) { //Llanthionine dihedral angle (N-Ca-Cb-S, from DALA) constraints for D_AA steriochem
			FuncOP n_ca_cb_s_harm0(utility::pointer::make_shared<CharmmPeriodicFunc>(0, 9.73, 0));
			ConstraintCOP n_ca_cb_s_cnst0(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm0));
			pose.add_constraint(n_ca_cb_s_cnst0);
			FuncOP n_ca_cb_s_harm1(utility::pointer::make_shared<CharmmPeriodicFunc>(0, 0.602, 1));
			ConstraintCOP n_ca_cb_s_cnst1(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm1));
			pose.add_constraint(n_ca_cb_s_cnst1);
			FuncOP n_ca_cb_s_harm2(utility::pointer::make_shared<CharmmPeriodicFunc>(0, -1.412, 2));
			ConstraintCOP n_ca_cb_s_cnst2(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm2));
			pose.add_constraint(n_ca_cb_s_cnst2);
			FuncOP n_ca_cb_s_harm3(utility::pointer::make_shared<CharmmPeriodicFunc>(0, -3.249, 3));
			ConstraintCOP n_ca_cb_s_cnst3(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm3));
			pose.add_constraint(n_ca_cb_s_cnst3);
			FuncOP n_ca_cb_s_harm4(utility::pointer::make_shared<CharmmPeriodicFunc>(0, -0.023, 4));
			ConstraintCOP n_ca_cb_s_cnst4(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm4));
			pose.add_constraint(n_ca_cb_s_cnst4);
			FuncOP n_ca_cb_s_harm5(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, 0.232, 1));
			ConstraintCOP n_ca_cb_s_cnst5(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm5));
			pose.add_constraint(n_ca_cb_s_cnst5);
			FuncOP n_ca_cb_s_harm6(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, -0.82, 2));
			ConstraintCOP n_ca_cb_s_cnst6(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm6));
			pose.add_constraint(n_ca_cb_s_cnst6);
			FuncOP n_ca_cb_s_harm7(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, -0.633, 3));
			ConstraintCOP n_ca_cb_s_cnst7(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm7));
			pose.add_constraint(n_ca_cb_s_cnst7);
			FuncOP n_ca_cb_s_harm8(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, -0.023, 4));
			ConstraintCOP n_ca_cb_s_cnst8(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm8));
			pose.add_constraint(n_ca_cb_s_cnst8);
		} else if ( pose.residue_type(dalares).is_l_aa() ) {
			//.base_name() == "ALA"
			//Lanthionine dihedral angle (N-Ca-Cb-S, from ALA) constraints for L_AA steriochem
			FuncOP n_ca_cb_s_harm0(utility::pointer::make_shared<CharmmPeriodicFunc>(0, 6.304, 0));
			ConstraintCOP n_ca_cb_s_cnst0(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm0));
			pose.add_constraint(n_ca_cb_s_cnst0);
			FuncOP n_ca_cb_s_harm1(utility::pointer::make_shared<CharmmPeriodicFunc>(0, 0.838, 1));
			ConstraintCOP n_ca_cb_s_cnst1(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm1));
			pose.add_constraint(n_ca_cb_s_cnst1);
			FuncOP n_ca_cb_s_harm2(utility::pointer::make_shared<CharmmPeriodicFunc>(0, 0.439, 2));
			ConstraintCOP n_ca_cb_s_cnst2(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm2));
			pose.add_constraint(n_ca_cb_s_cnst2);
			FuncOP n_ca_cb_s_harm3(utility::pointer::make_shared<CharmmPeriodicFunc>(0, -3.417, 3));
			ConstraintCOP n_ca_cb_s_cnst3(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm3));
			pose.add_constraint(n_ca_cb_s_cnst3);
			FuncOP n_ca_cb_s_harm4(utility::pointer::make_shared<CharmmPeriodicFunc>(0, 0.116, 4));
			ConstraintCOP n_ca_cb_s_cnst4(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm4));
			pose.add_constraint(n_ca_cb_s_cnst4);
			FuncOP n_ca_cb_s_harm5(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, -0.562, 1));
			ConstraintCOP n_ca_cb_s_cnst5(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm5));
			pose.add_constraint(n_ca_cb_s_cnst5);
			FuncOP n_ca_cb_s_harm6(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, 1.234, 2));
			ConstraintCOP n_ca_cb_s_cnst6(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm6));
			pose.add_constraint(n_ca_cb_s_cnst6);
			FuncOP n_ca_cb_s_harm7(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, 0.144, 3));
			ConstraintCOP n_ca_cb_s_cnst7(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm7));
			pose.add_constraint(n_ca_cb_s_cnst7);
			FuncOP n_ca_cb_s_harm8(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, 0.116, 4));
			ConstraintCOP n_ca_cb_s_cnst8(
				utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
				n_ca_cb_s_harm8));
			pose.add_constraint(n_ca_cb_s_cnst8);
		}

		{ //Lanthionine dihedral angle (Ca-Cb-S-Cb, from ALA) constraint
			FuncOP ala_ca_cb_s_cb_harm0(utility::pointer::make_shared<CharmmPeriodicFunc>(0, 4.892, 0));
			ConstraintCOP ala_ca_cb_s_cb_cnst0(
				utility::pointer::make_shared<DihedralConstraint>(atom_a, atom_b, atom_c, atom_d,
				ala_ca_cb_s_cb_harm0));
			pose.add_constraint(ala_ca_cb_s_cb_cnst0);
			FuncOP ala_ca_cb_s_cb_harm1(utility::pointer::make_shared<CharmmPeriodicFunc>(0, -2.008, 1));
			ConstraintCOP ala_ca_cb_s_cb_cnst1(
				utility::pointer::make_shared<DihedralConstraint>(atom_a, atom_b, atom_c, atom_d,
				ala_ca_cb_s_cb_harm1));
			pose.add_constraint(ala_ca_cb_s_cb_cnst1);
			FuncOP ala_ca_cb_s_cb_harm2(utility::pointer::make_shared<CharmmPeriodicFunc>(0, -0.383, 2));
			ConstraintCOP ala_ca_cb_s_cb_cnst2(
				utility::pointer::make_shared<DihedralConstraint>(atom_a, atom_b, atom_c, atom_d,
				ala_ca_cb_s_cb_harm2));
			pose.add_constraint(ala_ca_cb_s_cb_cnst2);
			FuncOP ala_ca_cb_s_cb_harm3(utility::pointer::make_shared<CharmmPeriodicFunc>(0, -1.45, 3));
			ConstraintCOP ala_ca_cb_s_cb_cnst3(
				utility::pointer::make_shared<DihedralConstraint>(atom_a, atom_b, atom_c, atom_d,
				ala_ca_cb_s_cb_harm3));
			pose.add_constraint(ala_ca_cb_s_cb_cnst3);
			FuncOP ala_ca_cb_s_cb_harm4(utility::pointer::make_shared<CharmmPeriodicFunc>(0, 0.37, 4));
			ConstraintCOP ala_ca_cb_s_cb_cnst4(
				utility::pointer::make_shared<DihedralConstraint>(atom_a, atom_b, atom_c, atom_d,
				ala_ca_cb_s_cb_harm4));
			pose.add_constraint(ala_ca_cb_s_cb_cnst4);
			FuncOP ala_ca_cb_s_cb_harm5(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, -0.156, 1));
			ConstraintCOP ala_ca_cb_s_cb_cnst5(
				utility::pointer::make_shared<DihedralConstraint>(atom_a, atom_b, atom_c, atom_d,
				ala_ca_cb_s_cb_harm5));
			pose.add_constraint(ala_ca_cb_s_cb_cnst5);
			FuncOP ala_ca_cb_s_cb_harm6(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, 1.384, 2));
			ConstraintCOP ala_ca_cb_s_cb_cnst6(
				utility::pointer::make_shared<DihedralConstraint>(atom_a, atom_b, atom_c, atom_d,
				ala_ca_cb_s_cb_harm6));
			pose.add_constraint(ala_ca_cb_s_cb_cnst6);
			FuncOP ala_ca_cb_s_cb_harm7(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, -0.045, 3));
			ConstraintCOP ala_ca_cb_s_cb_cnst7(
				utility::pointer::make_shared<DihedralConstraint>(atom_a, atom_b, atom_c, atom_d,
				ala_ca_cb_s_cb_harm7));
			pose.add_constraint(ala_ca_cb_s_cb_cnst7);
			FuncOP ala_ca_cb_s_cb_harm8(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, 0.37, 4));
			ConstraintCOP ala_ca_cb_s_cb_cnst8(
				utility::pointer::make_shared<DihedralConstraint>(atom_a, atom_b, atom_c, atom_d,
				ala_ca_cb_s_cb_harm8));
			pose.add_constraint(ala_ca_cb_s_cb_cnst8);

			//Lanthionine dihedral angle (Cb-S-Cb-Ca, from CYS) constraints
			FuncOP cys_ca_cb_s_cb_harm0(utility::pointer::make_shared<CharmmPeriodicFunc>(0, 4.592, 0));
			ConstraintCOP cys_ca_cb_s_cb_cnst0(
				utility::pointer::make_shared<DihedralConstraint>(atom_ao, atom_a, atom_b, atom_c,
				cys_ca_cb_s_cb_harm0));
			pose.add_constraint(cys_ca_cb_s_cb_cnst0);
			FuncOP cys_ca_cb_s_cb_harm1(utility::pointer::make_shared<CharmmPeriodicFunc>(0, -1.923, 1));
			ConstraintCOP cys_ca_cb_s_cb_cnst1(
				utility::pointer::make_shared<DihedralConstraint>(atom_ao, atom_a, atom_b, atom_c,
				cys_ca_cb_s_cb_harm1));
			pose.add_constraint(cys_ca_cb_s_cb_cnst1);
			FuncOP cys_ca_cb_s_cb_harm2(utility::pointer::make_shared<CharmmPeriodicFunc>(0, -0.463, 2));
			ConstraintCOP cys_ca_cb_s_cb_cnst2(
				utility::pointer::make_shared<DihedralConstraint>(atom_ao, atom_a, atom_b, atom_c,
				cys_ca_cb_s_cb_harm2));
			pose.add_constraint(cys_ca_cb_s_cb_cnst2);
			FuncOP cys_ca_cb_s_cb_harm3(utility::pointer::make_shared<CharmmPeriodicFunc>(0, -1.378, 3));
			ConstraintCOP cys_ca_cb_s_cb_cnst3(
				utility::pointer::make_shared<DihedralConstraint>(atom_ao, atom_a, atom_b, atom_c,
				cys_ca_cb_s_cb_harm3));
			pose.add_constraint(cys_ca_cb_s_cb_cnst3);
			FuncOP cys_ca_cb_s_cb_harm4(utility::pointer::make_shared<CharmmPeriodicFunc>(0, 0.384, 4));
			ConstraintCOP cys_ca_cb_s_cb_cnst4(
				utility::pointer::make_shared<DihedralConstraint>(atom_ao, atom_a, atom_b, atom_c,
				cys_ca_cb_s_cb_harm4));
			pose.add_constraint(cys_ca_cb_s_cb_cnst4);
			FuncOP cys_ca_cb_s_cb_harm5(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, -0.171, 1));
			ConstraintCOP cys_ca_cb_s_cb_cnst5(
				utility::pointer::make_shared<DihedralConstraint>(atom_ao, atom_a, atom_b, atom_c,
				cys_ca_cb_s_cb_harm5));
			pose.add_constraint(cys_ca_cb_s_cb_cnst5);
			FuncOP cys_ca_cb_s_cb_harm6(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, 1.255, 2));
			ConstraintCOP cys_ca_cb_s_cb_cnst6(
				utility::pointer::make_shared<DihedralConstraint>(atom_ao, atom_a, atom_b, atom_c,
				cys_ca_cb_s_cb_harm6));
			pose.add_constraint(cys_ca_cb_s_cb_cnst6);
			FuncOP cys_ca_cb_s_cb_harm7(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, -0.002, 3));
			ConstraintCOP cys_ca_cb_s_cb_cnst7(
				utility::pointer::make_shared<DihedralConstraint>(atom_ao, atom_a, atom_b, atom_c,
				cys_ca_cb_s_cb_harm7));
			pose.add_constraint(cys_ca_cb_s_cb_cnst7);
			FuncOP cys_ca_cb_s_cb_harm8(
				utility::pointer::make_shared<CharmmPeriodicFunc>(numeric::constants::d::pi_over_2, 0.384, 4));
			ConstraintCOP cys_ca_cb_s_cb_cnst8(
				utility::pointer::make_shared<DihedralConstraint>(atom_ao, atom_a, atom_b, atom_c,
				cys_ca_cb_s_cb_harm8));
			pose.add_constraint(cys_ca_cb_s_cb_cnst8);
		}
	} else {
		//There are four possible steriochem configurations for methyllanthionine linkers
		// the CA can be D(S)/L(R) and the CB can be R/S
		if ( pose.residue_type(dalares).is_d_aa() ) {
			//Methyllanthionine dihedral angle (N-Ca-Cb-S, from DBB/AAB) constraints for 2S (D_AA) steriochem
			FuncOP s2_ca_cb_s_cb_harm0( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 4.707, 0 ) );
			ConstraintCOP s2_ca_cb_s_cb_cnst0( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, s2_ca_cb_s_cb_harm0) );
			pose.add_constraint (s2_ca_cb_s_cb_cnst0);
			FuncOP s2_ca_cb_s_cb_harm1( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 0.073, 1 ) );
			ConstraintCOP s2_ca_cb_s_cb_cnst1( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, s2_ca_cb_s_cb_harm1) );
			pose.add_constraint (s2_ca_cb_s_cb_cnst1);
			FuncOP s2_ca_cb_s_cb_harm2( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 0.917, 2 ) );
			ConstraintCOP s2_ca_cb_s_cb_cnst2( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, s2_ca_cb_s_cb_harm2) );
			pose.add_constraint (s2_ca_cb_s_cb_cnst2);
			FuncOP s2_ca_cb_s_cb_harm3( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, -2.833, 3 ) );
			ConstraintCOP s2_ca_cb_s_cb_cnst3( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, s2_ca_cb_s_cb_harm3) );
			pose.add_constraint (s2_ca_cb_s_cb_cnst3);
			FuncOP s2_ca_cb_s_cb_harm4( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 0.011, 4 ) );
			ConstraintCOP s2_ca_cb_s_cb_cnst4( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, s2_ca_cb_s_cb_harm4) );
			pose.add_constraint (s2_ca_cb_s_cb_cnst4);
			FuncOP s2_ca_cb_s_cb_harm5( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, 0.43, 1 ) );
			ConstraintCOP s2_ca_cb_s_cb_cnst5( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, s2_ca_cb_s_cb_harm5) );
			pose.add_constraint (s2_ca_cb_s_cb_cnst5);
			FuncOP s2_ca_cb_s_cb_harm6( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, 1.276, 2 ) );
			ConstraintCOP s2_ca_cb_s_cb_cnst6( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, s2_ca_cb_s_cb_harm6) );
			pose.add_constraint (s2_ca_cb_s_cb_cnst6);
			FuncOP s2_ca_cb_s_cb_harm7( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, -0.58, 3 ) );
			ConstraintCOP s2_ca_cb_s_cb_cnst7( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, s2_ca_cb_s_cb_harm7) );
			pose.add_constraint (s2_ca_cb_s_cb_cnst7);
			FuncOP s2_ca_cb_s_cb_harm8( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, 0.011, 4 ) );
			ConstraintCOP s2_ca_cb_s_cb_cnst8( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, s2_ca_cb_s_cb_harm8) );
			pose.add_constraint (s2_ca_cb_s_cb_cnst8);
		} else if ( pose.residue_type(dalares).is_l_aa() ) {
			//Methyllanthionine dihedral angle (N-Ca-Cb-S, from DBB/AAB) constraints for 2R (L_AA) steriochem
			FuncOP r2_ca_cb_s_cb_harm0( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 4.626, 0 ) );
			ConstraintCOP r2_ca_cb_s_cb_cnst0( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, r2_ca_cb_s_cb_harm0) );
			pose.add_constraint (r2_ca_cb_s_cb_cnst0);
			FuncOP r2_ca_cb_s_cb_harm1( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 0.524, 1 ) );
			ConstraintCOP r2_ca_cb_s_cb_cnst1( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, r2_ca_cb_s_cb_harm1) );
			pose.add_constraint (r2_ca_cb_s_cb_cnst1);
			FuncOP r2_ca_cb_s_cb_harm2( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 1.18, 2 ) );
			ConstraintCOP r2_ca_cb_s_cb_cnst2( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, r2_ca_cb_s_cb_harm2) );
			pose.add_constraint (r2_ca_cb_s_cb_cnst2);
			FuncOP r2_ca_cb_s_cb_harm3( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, -3.111, 3 ) );
			ConstraintCOP r2_ca_cb_s_cb_cnst3( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, r2_ca_cb_s_cb_harm3) );
			pose.add_constraint (r2_ca_cb_s_cb_cnst3);
			FuncOP r2_ca_cb_s_cb_harm4( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, -0.074, 4 ) );
			ConstraintCOP r2_ca_cb_s_cb_cnst4( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, r2_ca_cb_s_cb_harm4) );
			pose.add_constraint (r2_ca_cb_s_cb_cnst4);
			FuncOP r2_ca_cb_s_cb_harm5( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, -0.458, 1 ) );
			ConstraintCOP r2_ca_cb_s_cb_cnst5( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, r2_ca_cb_s_cb_harm5) );
			pose.add_constraint (r2_ca_cb_s_cb_cnst5);
			FuncOP r2_ca_cb_s_cb_harm6( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, 1.396, 2 ) );
			ConstraintCOP r2_ca_cb_s_cb_cnst6( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, r2_ca_cb_s_cb_harm6) );
			pose.add_constraint (r2_ca_cb_s_cb_cnst6);
			FuncOP r2_ca_cb_s_cb_harm7( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, 0.368, 3 ) );
			ConstraintCOP r2_ca_cb_s_cb_cnst7( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, r2_ca_cb_s_cb_harm7) );
			pose.add_constraint (r2_ca_cb_s_cb_cnst7);
			FuncOP r2_ca_cb_s_cb_harm8( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, -0.074, 4 ) );
			ConstraintCOP r2_ca_cb_s_cb_cnst8( utility::pointer::make_shared< DihedralConstraint >( atom_b, atom_c, atom_d, atom_do, r2_ca_cb_s_cb_harm8) );
			pose.add_constraint (r2_ca_cb_s_cb_cnst8);
		}

		if ( pose.residue_type(dalares).base_name() == "DBS" || pose.residue_type(dalares).base_name() == "DDBS" ) {
			//Methyllanthionine dihedral angle (Ca-Cb-S-Cb, from CYS) constraints for 3S steriochem
			FuncOP s3_ca_cb_s_cb_harm0( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 3.8, 0 ) );
			ConstraintCOP s3_ca_cb_s_cb_cnst0( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, s3_ca_cb_s_cb_harm0) );
			pose.add_constraint (s3_ca_cb_s_cb_cnst0);
			FuncOP s3_ca_cb_s_cb_harm1( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, -1.478, 1 ) );
			ConstraintCOP s3_ca_cb_s_cb_cnst1( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, s3_ca_cb_s_cb_harm1) );
			pose.add_constraint (s3_ca_cb_s_cb_cnst1);
			FuncOP s3_ca_cb_s_cb_harm2( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 0.217, 2 ) );
			ConstraintCOP s3_ca_cb_s_cb_cnst2( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, s3_ca_cb_s_cb_harm2) );
			pose.add_constraint (s3_ca_cb_s_cb_cnst2);
			FuncOP s3_ca_cb_s_cb_harm3( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, -1.391, 3 ) );
			ConstraintCOP s3_ca_cb_s_cb_cnst3( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, s3_ca_cb_s_cb_harm3) );
			pose.add_constraint (s3_ca_cb_s_cb_cnst3);
			FuncOP s3_ca_cb_s_cb_harm4( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 0.246, 4 ) );
			ConstraintCOP s3_ca_cb_s_cb_cnst4( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, s3_ca_cb_s_cb_harm4) );
			pose.add_constraint (s3_ca_cb_s_cb_cnst4);
			FuncOP s3_ca_cb_s_cb_harm5( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, -0.336, 1 ) );
			ConstraintCOP s3_ca_cb_s_cb_cnst5( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, s3_ca_cb_s_cb_harm5) );
			pose.add_constraint (s3_ca_cb_s_cb_cnst5);
			FuncOP s3_ca_cb_s_cb_harm6( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, 1.164, 2 ) );
			ConstraintCOP s3_ca_cb_s_cb_cnst6( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, s3_ca_cb_s_cb_harm6) );
			pose.add_constraint (s3_ca_cb_s_cb_cnst6);
			FuncOP s3_ca_cb_s_cb_harm7( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, 0.114, 3 ) );
			ConstraintCOP s3_ca_cb_s_cb_cnst7( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, s3_ca_cb_s_cb_harm7) );
			pose.add_constraint (s3_ca_cb_s_cb_cnst7);
			FuncOP s3_ca_cb_s_cb_harm8( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, 0.246, 4 ) );
			ConstraintCOP s3_ca_cb_s_cb_cnst8( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, s3_ca_cb_s_cb_harm8) );
			pose.add_constraint (s3_ca_cb_s_cb_cnst8);
		} else if ( pose.residue_type(dalares).base_name() == "DBR" || pose.residue_type(dalares).base_name() == "DDBR" ) {
			//Methyllanthionine dihedral angle (Ca-Cb-S-Cb, from CYS) constraints for 3R steriochem
			FuncOP r3_ca_cb_s_cb_harm0( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 2.318, 0 ) );
			ConstraintCOP r3_ca_cb_s_cb_cnst0( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, r3_ca_cb_s_cb_harm0) );
			pose.add_constraint (r3_ca_cb_s_cb_cnst0);
			FuncOP r3_ca_cb_s_cb_harm1( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, -0.946, 1 ) );
			ConstraintCOP r3_ca_cb_s_cb_cnst1( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, r3_ca_cb_s_cb_harm1) );
			pose.add_constraint (r3_ca_cb_s_cb_cnst1);
			FuncOP r3_ca_cb_s_cb_harm2( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 0.183, 2 ) );
			ConstraintCOP r3_ca_cb_s_cb_cnst2( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, r3_ca_cb_s_cb_harm2) );
			pose.add_constraint (r3_ca_cb_s_cb_cnst2);
			FuncOP r3_ca_cb_s_cb_harm3( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, -1.323, 3 ) );
			ConstraintCOP r3_ca_cb_s_cb_cnst3( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, r3_ca_cb_s_cb_harm3) );
			pose.add_constraint (r3_ca_cb_s_cb_cnst3);
			FuncOP r3_ca_cb_s_cb_harm4( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 0.15, 4 ) );
			ConstraintCOP r3_ca_cb_s_cb_cnst4( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, r3_ca_cb_s_cb_harm4) );
			pose.add_constraint (r3_ca_cb_s_cb_cnst4);
			FuncOP r3_ca_cb_s_cb_harm5( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, 0.471, 1 ) );
			ConstraintCOP r3_ca_cb_s_cb_cnst5( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, r3_ca_cb_s_cb_harm5) );
			pose.add_constraint (r3_ca_cb_s_cb_cnst5);
			FuncOP r3_ca_cb_s_cb_harm6( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, 1.078, 2 ) );
			ConstraintCOP r3_ca_cb_s_cb_cnst6( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, r3_ca_cb_s_cb_harm6) );
			pose.add_constraint (r3_ca_cb_s_cb_cnst6);
			FuncOP r3_ca_cb_s_cb_harm7( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, -0.211, 3 ) );
			ConstraintCOP r3_ca_cb_s_cb_cnst7( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, r3_ca_cb_s_cb_harm7) );
			pose.add_constraint (r3_ca_cb_s_cb_cnst7);
			FuncOP r3_ca_cb_s_cb_harm8( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, 0.15, 4 ) );
			ConstraintCOP r3_ca_cb_s_cb_cnst8( utility::pointer::make_shared< DihedralConstraint >( atom_ao, atom_a, atom_b, atom_c, r3_ca_cb_s_cb_harm8) );
			pose.add_constraint (r3_ca_cb_s_cb_cnst8);
		}

		//Non-sterio dependent constraints
		{//Methyllanthionine dihedral angle (Ca-Cb-S-Cb, from DBB/AAB) constraint
			FuncOP cys_ca_cb_s_cb_harm0( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 4.86, 0 ) );
			ConstraintCOP cys_ca_cb_s_cb_cnst0( utility::pointer::make_shared< DihedralConstraint >( atom_a, atom_b, atom_c, atom_d, cys_ca_cb_s_cb_harm0) );
			pose.add_constraint (cys_ca_cb_s_cb_cnst0);
			FuncOP cys_ca_cb_s_cb_harm1( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, -2.577, 1 ) );
			ConstraintCOP cys_ca_cb_s_cb_cnst1( utility::pointer::make_shared< DihedralConstraint >( atom_a, atom_b, atom_c, atom_d, cys_ca_cb_s_cb_harm1) );
			pose.add_constraint (cys_ca_cb_s_cb_cnst1);
			FuncOP cys_ca_cb_s_cb_harm2( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, -0.294, 2 ) );
			ConstraintCOP cys_ca_cb_s_cb_cnst2( utility::pointer::make_shared< DihedralConstraint >( atom_a, atom_b, atom_c, atom_d, cys_ca_cb_s_cb_harm2) );
			pose.add_constraint (cys_ca_cb_s_cb_cnst2);
			FuncOP cys_ca_cb_s_cb_harm3( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, -1.159, 3 ) );
			ConstraintCOP cys_ca_cb_s_cb_cnst3( utility::pointer::make_shared< DihedralConstraint >( atom_a, atom_b, atom_c, atom_d, cys_ca_cb_s_cb_harm3) );
			pose.add_constraint (cys_ca_cb_s_cb_cnst3);
			FuncOP cys_ca_cb_s_cb_harm4( utility::pointer::make_shared< CharmmPeriodicFunc >( 0, 0.392, 4 ) );
			ConstraintCOP cys_ca_cb_s_cb_cnst4( utility::pointer::make_shared< DihedralConstraint >( atom_a, atom_b, atom_c, atom_d, cys_ca_cb_s_cb_harm4) );
			pose.add_constraint (cys_ca_cb_s_cb_cnst4);
			FuncOP cys_ca_cb_s_cb_harm5( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, -0.236, 1 ) );
			ConstraintCOP cys_ca_cb_s_cb_cnst5( utility::pointer::make_shared< DihedralConstraint >( atom_a, atom_b, atom_c, atom_d, cys_ca_cb_s_cb_harm5) );
			pose.add_constraint (cys_ca_cb_s_cb_cnst5);
			FuncOP cys_ca_cb_s_cb_harm6( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, 1.566, 2 ) );
			ConstraintCOP cys_ca_cb_s_cb_cnst6( utility::pointer::make_shared< DihedralConstraint >( atom_a, atom_b, atom_c, atom_d, cys_ca_cb_s_cb_harm6) );
			pose.add_constraint (cys_ca_cb_s_cb_cnst6);
			FuncOP cys_ca_cb_s_cb_harm7( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, -0.103, 3 ) );
			ConstraintCOP cys_ca_cb_s_cb_cnst7( utility::pointer::make_shared< DihedralConstraint >( atom_a, atom_b, atom_c, atom_d, cys_ca_cb_s_cb_harm7) );
			pose.add_constraint (cys_ca_cb_s_cb_cnst7);
			FuncOP cys_ca_cb_s_cb_harm8( utility::pointer::make_shared< CharmmPeriodicFunc >( numeric::constants::d::pi_over_2, 0.392, 4 ) );
			ConstraintCOP cys_ca_cb_s_cb_cnst8( utility::pointer::make_shared< DihedralConstraint >( atom_a, atom_b, atom_c, atom_d, cys_ca_cb_s_cb_harm8) );
			pose.add_constraint (cys_ca_cb_s_cb_cnst8);
		}
	}

	{ //Lanthionine bond angle constraints sd from methionine:
		FuncOP circharmfunc2a( utility::pointer::make_shared< CircularHarmonicFunc >( LANTHIONINE_UTIL_LANTHIONINE_BOND_S_ANGLE, 0.18) );
		FuncOP circharmfunc2b( utility::pointer::make_shared< CircularHarmonicFunc >( LANTHIONINE_UTIL_LANTHIONINE_BOND_C_ANGLE, 0.14) );
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

