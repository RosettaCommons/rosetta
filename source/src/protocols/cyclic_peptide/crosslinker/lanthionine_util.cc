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
//#include <core/scoring/func/CharmmPeriodicFunc.hh>
#include <core/scoring/func/AmberPeriodicFunc.hh>
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

	//The six atoms defining the lanthionine bond:
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
	AtomID const atom_bb_c(
		pose.residue_type(dalares).atom_index("C"),
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

	//Using the CYS N/C-CA-Cb-S and Met C-C-S-C parameters from the Amber ff14SB forcefield
	//Starting with N to S and C to S
	FuncOP n_ca_cb_s_harm1(utility::pointer::make_shared<AmberPeriodicFunc>(0, 0.602, 1));
	ConstraintCOP n_ca_cb_s_cnst1(
		utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
		n_ca_cb_s_harm1));
	pose.add_constraint(n_ca_cb_s_cnst1);
	FuncOP n_ca_cb_s_harm2(utility::pointer::make_shared<AmberPeriodicFunc>(numeric::constants::d::pi_over_2, 0.398, 2));
	ConstraintCOP n_ca_cb_s_cnst2(
		utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
		n_ca_cb_s_harm2));
	pose.add_constraint(n_ca_cb_s_cnst2);
	FuncOP n_ca_cb_s_harm3(utility::pointer::make_shared<AmberPeriodicFunc>(0, 0.323, 3));
	ConstraintCOP n_ca_cb_s_cnst3(
		utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
		n_ca_cb_s_harm3));
	pose.add_constraint(n_ca_cb_s_cnst3);
	FuncOP n_ca_cb_s_harm4(utility::pointer::make_shared<AmberPeriodicFunc>(0, 0.278, 4));
	ConstraintCOP n_ca_cb_s_cnst4(
		utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_do,
		n_ca_cb_s_harm4));
	pose.add_constraint(n_ca_cb_s_cnst4);
	FuncOP c_ca_cb_s_harm1(utility::pointer::make_shared<AmberPeriodicFunc>(0, 0.469, 1));
	ConstraintCOP c_ca_cb_s_cnst1(
		utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_bb_c,
		c_ca_cb_s_harm1));
	pose.add_constraint(c_ca_cb_s_cnst1);
	FuncOP c_ca_cb_s_harm2(utility::pointer::make_shared<AmberPeriodicFunc>(numeric::constants::d::pi_over_2, 0.021, 2));
	ConstraintCOP c_ca_cb_s_cnst2(
		utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_bb_c,
		c_ca_cb_s_harm2));
	pose.add_constraint(c_ca_cb_s_cnst2);
	FuncOP c_ca_cb_s_harm3(utility::pointer::make_shared<AmberPeriodicFunc>(0, 0.323, 3));
	ConstraintCOP c_ca_cb_s_cnst3(
		utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_bb_c,
		c_ca_cb_s_harm3));
	pose.add_constraint(c_ca_cb_s_cnst3);
	FuncOP c_ca_cb_s_harm4(utility::pointer::make_shared<AmberPeriodicFunc>(0, 0.064, 4));
	ConstraintCOP c_ca_cb_s_cnst4(
		utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_d, atom_bb_c,
		c_ca_cb_s_harm4));
	pose.add_constraint(c_ca_cb_s_cnst4);

	//Now adding C-C-S-C - starting with ALA side
	FuncOP ala_ca_cb_s_cb_harm1(utility::pointer::make_shared<AmberPeriodicFunc>(0, 0.057, 1));
	ConstraintCOP ala_ca_cb_s_cb_cnst1(
		utility::pointer::make_shared<DihedralConstraint>(atom_a, atom_b, atom_c, atom_d,
		ala_ca_cb_s_cb_harm1));
	pose.add_constraint(ala_ca_cb_s_cb_cnst1);
	FuncOP ala_ca_cb_s_cb_harm2(utility::pointer::make_shared<AmberPeriodicFunc>(0, 0.414, 1));
	ConstraintCOP ala_ca_cb_s_cb_cnst2(
		utility::pointer::make_shared<DihedralConstraint>(atom_a, atom_b, atom_c, atom_d,
		ala_ca_cb_s_cb_harm2));
	pose.add_constraint(ala_ca_cb_s_cb_cnst2);
	FuncOP ala_ca_cb_s_cb_harm3(utility::pointer::make_shared<AmberPeriodicFunc>(0, 0.442, 1));
	ConstraintCOP ala_ca_cb_s_cb_cnst3(
		utility::pointer::make_shared<DihedralConstraint>(atom_a, atom_b, atom_c, atom_d,
		ala_ca_cb_s_cb_harm3));
	pose.add_constraint(ala_ca_cb_s_cb_cnst3);
	FuncOP ala_ca_cb_s_cb_harm4(utility::pointer::make_shared<AmberPeriodicFunc>(0, 0.247, 1));
	ConstraintCOP ala_ca_cb_s_cb_cnst4(
		utility::pointer::make_shared<DihedralConstraint>(atom_a, atom_b, atom_c, atom_d,
		ala_ca_cb_s_cb_harm4));
	pose.add_constraint(ala_ca_cb_s_cb_cnst4);

	//Now for the CYS angle
	FuncOP cys_ca_cb_s_cb_harm1(utility::pointer::make_shared<AmberPeriodicFunc>(0, 0.057, 1));
	ConstraintCOP cys_ca_cb_s_cb_cnst1(
		utility::pointer::make_shared<DihedralConstraint>(atom_ao, atom_a, atom_b, atom_c,
		cys_ca_cb_s_cb_harm1));
	pose.add_constraint(cys_ca_cb_s_cb_cnst1);
	FuncOP cys_ca_cb_s_cb_harm2(utility::pointer::make_shared<AmberPeriodicFunc>(0, 0.414, 1));
	ConstraintCOP cys_ca_cb_s_cb_cnst2(
		utility::pointer::make_shared<DihedralConstraint>(atom_ao, atom_a, atom_b, atom_c,
		cys_ca_cb_s_cb_harm2));
	pose.add_constraint(cys_ca_cb_s_cb_cnst2);
	FuncOP cys_ca_cb_s_cb_harm3(utility::pointer::make_shared<AmberPeriodicFunc>(0, 0.442, 1));
	ConstraintCOP cys_ca_cb_s_cb_cnst3(
		utility::pointer::make_shared<DihedralConstraint>(atom_ao, atom_a, atom_b, atom_c,
		cys_ca_cb_s_cb_harm3));
	pose.add_constraint(cys_ca_cb_s_cb_cnst3);
	FuncOP cys_ca_cb_s_cb_harm4(utility::pointer::make_shared<AmberPeriodicFunc>(0, 0.247, 1));
	ConstraintCOP cys_ca_cb_s_cb_cnst4(
		utility::pointer::make_shared<DihedralConstraint>(atom_ao, atom_a, atom_b, atom_c,
		cys_ca_cb_s_cb_harm4));
	pose.add_constraint(cys_ca_cb_s_cb_cnst4);

	{ //Lanthionine bond angle constraints sd from methionine:
		FuncOP circharmfunc2a( utility::pointer::make_shared< CircularHarmonicFunc >( LANTHIONINE_UTIL_LANTHIONINE_BOND_S_ANGLE, 0.18) );
		FuncOP circharmfunc2b( utility::pointer::make_shared< CircularHarmonicFunc >( LANTHIONINE_UTIL_LANTHIONINE_BOND_C_ANGLE, 0.14) );
		ConstraintCOP angleconst1( utility::pointer::make_shared< AngleConstraint >( atom_a, atom_b, atom_c, circharmfunc2a) );
		ConstraintCOP angleconst2( utility::pointer::make_shared< AngleConstraint >( atom_b, atom_c, atom_d, circharmfunc2b) );
		pose.add_constraint (angleconst1);
		//This is adding the SCC angle for lanthionine. methyllanthionine is slightly different
		if ( pose.residue_type(dalares).base_name() == "DALA" || pose.residue_type(dalares).base_name() == "ALA" ) {
			pose.add_constraint (angleconst2);
		}
		//if methylanthionine, add additional constraints to ensure good geometry around the methyl extension
		if ( pose.residue_type(dalares).base_name() == "DBR" || pose.residue_type(dalares).base_name() == "DDBR"  ||
				pose.residue_type(dalares).base_name() == "DBS" || pose.residue_type(dalares).base_name() == "DDBS" ) {
			//First, use slightly different angle for methyllanthionine S-CB-CA
			FuncOP circharmfunc2c( utility::pointer::make_shared< CircularHarmonicFunc >( LANTHIONINE_UTIL_METHYLLANTHIONINE_BOND_CA_ANGLE, 0.14) );
			ConstraintCOP angleconst3( utility::pointer::make_shared< AngleConstraint >( atom_b, atom_c, atom_d, circharmfunc2c) );
			pose.add_constraint (angleconst3);
			//Now, prep for the S-CB-CG angle
			AtomID const atom_cg(
				pose.residue_type(dalares).atom_index("CG"),
				dalares
			);
			//Adding an angle constraint for the CG-CB-S connection abt 115 degrees
			//distance is ~2.76
			TR << "Adding additional angle constraint for the lanthionine bond:" << std::endl;
			TR << "1.\tRes=" << atom_cg.rsd() << "\tAtom=" << pose.residue(atom_cg.rsd()).atom_name(atom_cg.atomno()) << std::endl;
			TR << "2.\tRes=" << atom_c.rsd() << "\tAtom=" << pose.residue(atom_c.rsd()).atom_name(atom_c.atomno()) << std::endl;
			TR << "3.\tRes=" << atom_b.rsd() << "\tAtom=" << pose.residue(atom_b.rsd()).atom_name(atom_b.atomno()) << std::endl;
			FuncOP circharmfunc2d( utility::pointer::make_shared< CircularHarmonicFunc >( LANTHIONINE_UTIL_METHYLLANTHIONINE_BOND_CG_ANGLE, 0.1) ); //from 0.18
			ConstraintCOP angleconst4( utility::pointer::make_shared< AngleConstraint >( atom_cg, atom_c, atom_b, circharmfunc2d) );
			pose.add_constraint (angleconst4);
			//Finally strengthen the CA-CB-CG angle term for doing cartesian min
			FuncOP circharmfunc2e( utility::pointer::make_shared< CircularHarmonicFunc >( LANTHIONINE_UTIL_METHYLLANTHIONINE_BOND_CB_ANGLE, 0.18) );
			ConstraintCOP angleconst5( utility::pointer::make_shared< AngleConstraint >( atom_cg, atom_c, atom_d, circharmfunc2e) );
			pose.add_constraint (angleconst5);
			//Improper dihedral param as well
			//FuncOP harmdistfunc( utility::pointer::make_shared< HarmonicFunc >( 2.76, 0.07));
			//ConstraintCOP distcst ( utility::pointer::make_shared< AtomPairConstraint >( atom_cg, atom_b, harmdistfunc) );
			//pose.add_constraint( distcst );
			//The S-CB-CB-CA dihedral is 122.6 degrees with wb97xd/6-311+g(d,p) geom opt
			if ( pose.residue_type(dalares).base_name() == "DBS" || pose.residue_type(dalares).base_name() == "DDBS" ) {
				FuncOP circharmfuncdih(utility::pointer::make_shared<CircularHarmonicFunc>(2.139, 0.14));
				ConstraintCOP dihconst(
					utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_cg, atom_d,
					circharmfuncdih));
				pose.add_constraint(dihconst);
			} else { //must be DBR
				FuncOP circharmfuncdih(utility::pointer::make_shared<CircularHarmonicFunc>(-2.139, 0.14));
				ConstraintCOP dihconst(
					utility::pointer::make_shared<DihedralConstraint>(atom_b, atom_c, atom_cg, atom_d,
					circharmfuncdih));
				pose.add_constraint(dihconst);
			}
		}
	}

	TR << "Finished setting up constraints." << std::endl;
}

} //crosslinker
} //cyclic_peptide
} //protocols

