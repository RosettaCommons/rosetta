// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/CyclizationMover.hh
/// @brief Implimentation file for CyclizationMover
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

// unit headers
#include <protocols/simple_moves/CyclizationMover.hh>

// protocols headers
#include <protocols/minimization_packing/MinMover.hh>

// core headers
#include <core/types.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueConnection.hh>

#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>

#include <core/scoring/func/CharmmPeriodicFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>

#include <core/kinematics/MoveMap.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// numeric
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>

// c++ headers
#include <string>
#include <algorithm>
#include <utility>

namespace protocols {
namespace simple_moves {

static basic::Tracer TR( "protocols.simple_moves.CyclizationMover" );

/// @brief Default constructor
CyclizationMover::CyclizationMover( core::Size chain_to_cyclize, bool add_constraints = true, bool minimize = true, core::Size minimization_rebuild_rounds = 3 ) :
	protocols::moves::Mover("CyclizationMover"),
	chain_to_cyclize_( chain_to_cyclize ),
	nterm_rsd_num_( 0 ),
	cterm_rsd_num_( 0 ),
	add_constraints_( add_constraints ),
	minimize_( minimize ),
	minimization_rebuild_rounds_( minimization_rebuild_rounds ),
	score_fxn_( /* 0 */ ),
	move_map_( /* 0 */ ),
	mm_torsion_library_(core::scoring::ScoringManager::get_instance()->get_MMTorsionLibrary() )
{}

CyclizationMover::CyclizationMover( core::Size chain_to_cyclize, bool add_constraints, bool minimize, core::Size minimization_rebuild_rounds, core::scoring::ScoreFunctionOP score_fxn, core::kinematics::MoveMapOP move_map ) :
	chain_to_cyclize_( chain_to_cyclize ),
	nterm_rsd_num_( 0 ),
	cterm_rsd_num_( 0 ),
	add_constraints_( add_constraints ),
	minimize_( minimize ),
	minimization_rebuild_rounds_( minimization_rebuild_rounds ),
	score_fxn_(std::move( score_fxn )),
	move_map_(std::move( move_map )),
	mm_torsion_library_( core::scoring::ScoringManager::get_instance()->get_MMTorsionLibrary() )
{}

/// @brief Sets up inter residue cyclic connections and potentially adds constraints, and minimizes the pose
void
CyclizationMover::apply( core::pose::Pose & pose )
{
	using namespace core;

	// check to see if specified chain exists
	runtime_assert( chain_to_cyclize_ <= pose.conformation().num_chains() );

	// get residue numbers of chain N-terminus and C-terminus
	nterm_rsd_num_ = pose.conformation().chain_begin( chain_to_cyclize_ );
	cterm_rsd_num_ = pose.conformation().chain_end( chain_to_cyclize_ );

	// setup residue connections and constraints
	setup_connections( pose );

	// add constraints to maintain the cyclization to the pose
	if ( add_constraints_ ) {
		setup_constraints ( pose );
	}

	// minimize the pose to bring
	if ( minimize_ && minimization_rebuild_rounds_ > 0 ) {
		setup_scorefunction();
		setup_minimizer( pose );
		for ( Size i( 1 ); i <= minimization_rebuild_rounds_; ++i ) {
			minimize_rebuild( pose );
		}
	}
}

/// @brief Modifes terminal ResidueTypes to cyclized variants and then uses connformation::detect_bonds() to create a residue connection
void
CyclizationMover::setup_connections( core::pose::Pose & pose )
{

	using namespace core;
	using namespace chemical;

	TR << "Setting up residue connections..." << std::endl;

	// make sure the seqence positions have been initialized to meaningful values
	runtime_assert( nterm_rsd_num_ != 0 && cterm_rsd_num_ != 0 );

	// make sure the two residuetypes are peptide or peptoid as that is all the NtermConnect and CtermConnect patches apply to
	runtime_assert( pose.residue( nterm_rsd_num_ ).type().is_protein() || pose.residue( nterm_rsd_num_ ).type().is_peptoid() );
	runtime_assert( pose.residue( cterm_rsd_num_ ).type().is_protein() || pose.residue( cterm_rsd_num_ ).type().is_peptoid() );

	// AMW: as of now, we don't need to add special patches for cyclization!
	// The actual upper and lower connects can be used.
	
	// get types name for N-terminus and C-terminus (manipulating strings like this is a little hacky)
	std::string nterm_connect_type_name(
		pose.residue_type( nterm_rsd_num_ ).is_peptoid() ?
		pose.residue( nterm_rsd_num_ ).type().name3() + ":peptoid_cutpoint_upper" :
		pose.residue( nterm_rsd_num_ ).type().name3() + ":protein_cutpoint_upper"
	);
	std::string cterm_connect_type_name(
		pose.residue( cterm_rsd_num_ ).type().name3() + ":protein_cutpoint_lower" //Both peptoids and proteins use the protein_cutpoint_lower patch.
	);

	// remove spaces if the name3 really only has 2 letters, Damn it Tim!
	//nterm_connect_type_name.erase( std::remove( nterm_connect_type_name.begin(), nterm_connect_type_name.end(), ' ' ), nterm_connect_type_name.end() );
	//cterm_connect_type_name.erase( std::remove( cterm_connect_type_name.begin(), cterm_connect_type_name.end(), ' ' ), cterm_connect_type_name.end() );

	// get CtermConnect and NtermConnect variant types
	ResidueTypeSetCOP rsd_type_set( pose.residue_type_set_for_pose( FULL_ATOM_t ) );
	ResidueType const & nterm_connect_type( rsd_type_set->name_map( nterm_connect_type_name ) );
	ResidueType const & cterm_connect_type( rsd_type_set->name_map( cterm_connect_type_name ) );

	// replace N-terminus and C-terminus with connect varients aligning new residue to old ones
	pose::replace_pose_residue_copying_existing_coordinates( pose, nterm_rsd_num_, nterm_connect_type );
	pose::replace_pose_residue_copying_existing_coordinates( pose, cterm_rsd_num_, cterm_connect_type );

	// tell rosetta to detect the new bonds
	pose.conformation().detect_bonds();
	pose.conformation().rebuild_residue_connection_dependent_atoms( nterm_rsd_num_, 2 );
	pose.conformation().rebuild_residue_connection_dependent_atoms( cterm_rsd_num_, 2 );

	// TODO: replace these calls with something more specfic
	//pose.conformation().show_residue_connections();
	//pose.fold_tree();
}

/// @brief Creates constraints to maintain the proper conformation and adds it to the pose.
void
CyclizationMover::setup_constraints( core::pose::Pose & pose )
{
	// TODO: make the constraints funcs based on mm_strech, mm_bend, mm_twist

	using namespace core;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring::constraints;
	using namespace scoring::func;

	TR << "Setting up constraints to maintain cycle based on polymeric base types..." << std::endl;

	// get connect variants of ResidueTypes of N-terminus and C-terminus
	ResidueType const & nterm_connect_type( pose.residue( nterm_rsd_num_ ).type() );
	ResidueType const & cterm_connect_type( pose.residue( cterm_rsd_num_ ).type() );

	// get base variants of ResidueTypes of N-terminus and C-terminus
	ResidueTypeSetCOP rsd_type_set( pose.residue_type_set_for_pose( FULL_ATOM_t ) );
	// remove spaces if the name3 really only has 2 letters, Damn it Tim!
	std::string nterm_base_type_name( nterm_connect_type.name3() );
	std::string cterm_base_type_name( cterm_connect_type.name3() );
	nterm_base_type_name.erase( std::remove( nterm_base_type_name.begin(), nterm_base_type_name.end(), ' ' ), nterm_base_type_name.end() );
	cterm_base_type_name.erase( std::remove( cterm_base_type_name.begin(), cterm_base_type_name.end(), ' ' ), cterm_base_type_name.end() );

	ResidueType const & nterm_base_type( rsd_type_set->name_map( nterm_base_type_name ) );
	ResidueType const & cterm_base_type( rsd_type_set->name_map( cterm_base_type_name ) );

	// make sure the two residuetypes are peptide or peptoid as that is all the NtermConnect and CtermConnect patches apply to
	runtime_assert( nterm_connect_type.is_protein() || nterm_connect_type.is_peptoid() );
	runtime_assert( cterm_connect_type.is_protein() || cterm_connect_type.is_peptoid() );

	// get ideal internal coordinates of nterm lower connect and cterm upper connect
	AtomICoor nterm_base_lc_icoor( nterm_base_type.lower_connect().icoor() );
	AtomICoor cterm_base_uc_icoor( cterm_base_type.upper_connect().icoor() );

	TR << "DEBUG NTERM LOWER ICOOR:\t" << nterm_base_type.lower_connect_id() << "\t"
		<< nterm_base_lc_icoor.d()     << "\t" << nterm_base_lc_icoor.stub_atom1().atomno()  << "\t"
		<< nterm_base_lc_icoor.theta() << "\t" << nterm_base_lc_icoor.stub_atom2().atomno()  << "\t"
		<< nterm_base_lc_icoor.phi()   << "\t" << nterm_base_lc_icoor.stub_atom3().atomno()  << "\t"
		<< std::endl;

	TR << "DEBUG CTERM UPPER ICOOR:\t" << cterm_base_type.upper_connect_id() << "\t"
		<< cterm_base_uc_icoor.d()     << "\t" << cterm_base_uc_icoor.stub_atom1().atomno()  << "\t"
		<< cterm_base_uc_icoor.theta() << "\t" << cterm_base_uc_icoor.stub_atom2().atomno()  << "\t"
		<< cterm_base_uc_icoor.phi()   << "\t" << cterm_base_uc_icoor.stub_atom3().atomno()  << "\t"
		<< std::endl;

	// create AtomIDs for the two central atoms and the atoms one bond away from them
	id::AtomID nterm_n(  nterm_base_lc_icoor.stub_atom1().atomno(), nterm_rsd_num_ );
	id::AtomID nterm_ca( nterm_base_lc_icoor.stub_atom2().atomno(), nterm_rsd_num_ );
	id::AtomID cterm_o(  cterm_base_type.atom_index( "O" ),         cterm_rsd_num_ );
	id::AtomID cterm_c(  cterm_base_uc_icoor.stub_atom1().atomno(), cterm_rsd_num_ );
	id::AtomID cterm_ca( cterm_base_uc_icoor.stub_atom2().atomno(), cterm_rsd_num_ );

	// create an AtomPairConstraint between the two central atoms
	TR << "Nterm atom " << nterm_base_lc_icoor.stub_atom1().atomno() << " is " << nterm_base_lc_icoor.d() << " from its lower connect" << std::endl;
	TR << "Cterm atom " << cterm_base_uc_icoor.stub_atom1().atomno() << " is " << cterm_base_uc_icoor.d() << " from its upper connect" << std::endl;
	Real ap_cst_length( ( nterm_base_lc_icoor.d() + cterm_base_uc_icoor.d() ) / 2 );
	TR << "Adding AtomPairConstraint of length:" << ap_cst_length << std::endl;
	ConstraintOP b1( new AtomPairConstraint( nterm_n, cterm_c, core::scoring::func::FuncOP( new HarmonicFunc( ap_cst_length, 0.01 ) ) ) );
	pose.add_constraint( b1 );

	// create three AngleConstraints
	TR << "Nterm atom " << nterm_base_lc_icoor.stub_atom2().atomno() << " and " << nterm_base_lc_icoor.stub_atom1().atomno() << " make and angle of " << nterm_base_lc_icoor.theta() << "with its lower connect" << std::endl;
	Real a1_cst_radian( numeric::constants::r::pi - nterm_base_lc_icoor.theta() );
	TR << "Adding AngleConstraint of radian: " << a1_cst_radian << std::endl;
	ConstraintOP a1( new AngleConstraint( nterm_ca, nterm_n, cterm_c, core::scoring::func::FuncOP( new CircularHarmonicFunc( a1_cst_radian, 0.1 ) ) ) );
	pose.add_constraint( a1 );

	TR << "Cterm atom " << cterm_base_uc_icoor.stub_atom2().atomno() << " and " << cterm_base_uc_icoor.stub_atom1().atomno() << " make and angle of " << cterm_base_uc_icoor.theta() << "with its upper connect" << std::endl;
	Real a2_cst_radian( numeric::constants::r::pi - cterm_base_uc_icoor.theta() );
	TR << "Adding AngleConstraint of radian: " << a2_cst_radian << std::endl;
	ConstraintOP a2( new AngleConstraint( nterm_n, cterm_c, cterm_ca, core::scoring::func::FuncOP( new CircularHarmonicFunc( a2_cst_radian, 0.1 ) ) ) );
	pose.add_constraint( a2 );

	TR << "Nterm atom N with Cterm atom C and O makes an angle of ~123.4 degrees" << std::endl;
	Real a3_cst_radian( numeric::constants::r::pi *123.4/180.0 );
	TR << "Adding AngleConstraint of radian: " << a3_cst_radian << std::endl;
	ConstraintOP a3( new AngleConstraint( nterm_n, cterm_c, cterm_o, core::scoring::func::FuncOP( new CircularHarmonicFunc( a3_cst_radian, 0.05 ) ) ) );
	pose.add_constraint( a3 );

	// create improper torsion constraint for CO planarity
	TR << "Nterm atom N with Cterm atom C and O and CA makes an improper torsion of ~180 degrees" << std::endl;
	Real d4_cst_radian( numeric::constants::r::pi );
	TR << "Adding AngleConstraint of radian: " << d4_cst_radian << std::endl;
	ConstraintOP d4( new DihedralConstraint( nterm_n, cterm_ca, cterm_c, cterm_o, core::scoring::func::FuncOP( new CircularHarmonicFunc( d4_cst_radian, 0.05 ) ) ) );
	pose.add_constraint( d4 );


	// create dihedral constraints based on MMTorsion term

	/*
	std::string nterm_n_name(  nterm_connect_type.atom_name( nterm_base_lc_icoor.stub_atom1().atomno() ) );
	std::string nterm_ca_name( nterm_connect_type.atom_name( nterm_base_lc_icoor.stub_atom2().atomno() ) );
	std::string cterm_c_name(  cterm_connect_type.atom_name( cterm_base_uc_icoor.stub_atom1().atomno() ) );
	std::string cterm_ca_name( cterm_connect_type.atom_name( cterm_base_uc_icoor.stub_atom2().atomno() ) );
	std::string nterm_3_name(  nterm_connect_type.atom_name( nterm_base_lc_icoor.stub_atom3().atomno() ));
	std::string cterm_3_name(  cterm_connect_type.atom_name( cterm_base_uc_icoor.stub_atom3().atomno() ));

	TR << "Adding dihedral constraint between " << nterm_rsd_num_ << ":" << nterm_n_name << " " << nterm_rsd_num_ << ":" << nterm_ca_name << " " << cterm_rsd_num_ << ":" << cterm_c_name << " " << cterm_rsd_num_ << ":" << cterm_ca_name << std::endl;
	TR << "Atom 3 names: " << nterm_rsd_num_ << ":" << nterm_3_name << " " << cterm_rsd_num_ << ":" << cterm_3_name << std::endl;
	*/


	// get residue types
	conformation::Residue const & rsd1 = pose.residue( nterm_rsd_num_ );
	conformation::Residue const & rsd2 = pose.residue( cterm_rsd_num_ );
	chemical::ResidueType const & rsd1_type = pose.residue( nterm_rsd_num_ ).type();
	chemical::ResidueType const & rsd2_type = pose.residue( cterm_rsd_num_ ).type();

	TR << "residue_pair_energy: processing residues " << rsd1.seqpos() << "." << rsd1_type.name() << "-" << rsd2.seqpos() << "." << rsd2_type.name() << std::endl;

	utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );

	for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {

		Size const resconn_id1( r1_resconn_ids[ii] );
		Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

		Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
		Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

		TR << "Found residue connection id " << resconn_id1 << "-" << resconn_id2 << ": " << rsd1.atom_name( resconn_atomno1 ) << "-" << rsd2.atom_name( resconn_atomno2 ) << std::endl;

		Size const resconn_mmat1 = rsd1_type.atom( resconn_atomno1 ).mm_atom_type_index();
		Size const resconn_mmat2 = rsd2_type.atom( resconn_atomno2 ).mm_atom_type_index();

		/// Iterate across all atom-quadrouples that define dihedral angles spanning the interface.
		/// 1. iterate over all pairs of pairs within 1 bond of either residue connection atom.
		/// 2. iterate over all triples on residue 1 within 2 bonds of resconn_atomno1.
		/// 3. iterate over all triples on residue 2 within 2 bonds of resconn_atomno2.


		{ // Scope Section 1.
			Size const mmat2 = resconn_mmat1;
			Size const mmat3 = resconn_mmat2;

			utility::vector1< chemical::two_atom_set > const & rsd1_atoms_wi1_bond_of_ii(
				rsd1_type.atoms_within_one_bond_of_a_residue_connection( resconn_id1 ));

			utility::vector1< chemical::two_atom_set > const & rsd2_atoms_wi1_bond_of_ii(
				rsd2_type.atoms_within_one_bond_of_a_residue_connection( resconn_id2 ));

			for ( Size jj = 1; jj <= rsd1_atoms_wi1_bond_of_ii.size(); ++jj ) {
				debug_assert( rsd1_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno1 );
				Size const jj_term_atomno = rsd1_atoms_wi1_bond_of_ii[ jj ].key2();
				Size const mmat1 = rsd1_type.atom( jj_term_atomno ).mm_atom_type_index();

				for ( Size kk = 1; kk <= rsd2_atoms_wi1_bond_of_ii.size(); ++kk ) {
					debug_assert( rsd2_atoms_wi1_bond_of_ii[ kk ].key1() == resconn_atomno2 );
					Size const kk_term_atomno = rsd2_atoms_wi1_bond_of_ii[ kk ].key2();
					Size const mmat4 = rsd2_type.atom( kk_term_atomno ).mm_atom_type_index();

					//Real const angle = numeric::dihedral_radians(
					// rsd1.xyz( jj_term_atomno ),
					// rsd1.xyz( resconn_atomno1 ),
					// rsd2.xyz( resconn_atomno2 ),
					// rsd2.xyz( kk_term_atomno ) );

					TR << "Section 1:"
						<< "r1 " << jj_term_atomno  << " " << rsd1.atom_name( jj_term_atomno )  << "(" << rsd1_type.atom( jj_term_atomno ).mm_name()  << ") - "
						<< "r1 " << resconn_atomno1 << " " << rsd1.atom_name( resconn_atomno1 ) << "(" << rsd1_type.atom( resconn_atomno1 ).mm_name() << ") - "
						<< "r2 " << resconn_atomno2 << " " << rsd2.atom_name( resconn_atomno2 ) << "(" << rsd2_type.atom( resconn_atomno2 ).mm_name() << ") - "
						<< "r2 " << kk_term_atomno  << " " << rsd2.atom_name( kk_term_atomno )  << "(" << rsd2_type.atom( kk_term_atomno ).mm_name()  << ")" << std::endl;

					// get mm params and add constraints for each set of params
					scoring::mm::mm_torsion_library_citer_pair pair = mm_torsion_library_.lookup( mmat1, mmat2, mmat3, mmat4 );

					id::AtomID a1( jj_term_atomno, nterm_rsd_num_ );
					id::AtomID a2( resconn_atomno1, nterm_rsd_num_ );
					id::AtomID a3( resconn_atomno2, cterm_rsd_num_ );
					id::AtomID a4( kk_term_atomno, cterm_rsd_num_ );

					for ( auto i = pair.first, e = pair.second; i != e; ++i ) {
						scoring::func::FuncOP cpf_temp( new CharmmPeriodicFunc( (i->second).key3() /*x0*/, (i->second).key1() /*n*/, (i->second).key2() /*k*/ ) );
						ConstraintOP d_temp( new DihedralConstraint( a1, a2, a3, a4, cpf_temp ) );
						pose.add_constraint( d_temp );
						TR << "Added Charmm constraint x0/n/k: " <<  (i->second).key3() << "/" << (i->second).key1() << "/" << (i->second).key2() << std::endl;
					}

				}

			}
		} // end Scope section 1.

		{ // Scope section 2.
			Size const mmat3 = resconn_mmat1;
			Size const mmat4 = resconn_mmat2;

			utility::vector1< chemical::three_atom_set > const & rsd1_atoms_wi2_bonds_of_ii(
				rsd1_type.atoms_within_two_bonds_of_a_residue_connection( resconn_id1 ));

			for ( Size jj = 1; jj <= rsd1_atoms_wi2_bonds_of_ii.size(); ++jj ) {
				debug_assert( rsd1_atoms_wi2_bonds_of_ii[ jj ].key1() == resconn_atomno1 );

				Size const jj_atom2 = rsd1_atoms_wi2_bonds_of_ii[ jj ].key2();
				Size const mmat2 = rsd1_type.atom( jj_atom2 ).mm_atom_type_index();

				Size const jj_atom1 = rsd1_atoms_wi2_bonds_of_ii[ jj ].key3();
				Size const mmat1 = rsd1_type.atom( jj_atom1 ).mm_atom_type_index();

				//Real const angle = numeric::dihedral_radians(
				// rsd1.xyz( jj_atom1 ),
				// rsd1.xyz( jj_atom2 ),
				// rsd1.xyz( resconn_atomno1 ),
				// rsd2.xyz( resconn_atomno2 ) );

				TR << "Section 2: "
					<< "r1 " << jj_atom1        << " " << rsd1.atom_name( jj_atom1 )        << "(" << rsd1_type.atom( jj_atom1 ).mm_name()  << ") - "
					<< "r1 " << jj_atom2        << " " << rsd1.atom_name( jj_atom2 )        << "(" << rsd1_type.atom( jj_atom2 ).mm_name()  << ") - "
					<< "r1 " << resconn_atomno1 << " " << rsd1.atom_name( resconn_atomno1 ) << "(" << rsd1_type.atom( resconn_atomno1 ).mm_name() << ") - "
					<< "r2 " << resconn_atomno2 << " " << rsd2.atom_name( resconn_atomno2 ) << "(" << rsd2_type.atom( resconn_atomno2 ).mm_name() << ")" << std::endl;

				// get mm params and add constraints for each set of params
				scoring::mm::mm_torsion_library_citer_pair pair = mm_torsion_library_.lookup( mmat1, mmat2, mmat3, mmat4 );

				id::AtomID a1( jj_atom1, nterm_rsd_num_ );
				id::AtomID a2( jj_atom2, nterm_rsd_num_ );
				id::AtomID a3( resconn_atomno1, nterm_rsd_num_ );
				id::AtomID a4( resconn_atomno2, cterm_rsd_num_ );

				for ( auto i = pair.first, e = pair.second; i != e; ++i ) {
					scoring::func::FuncOP cpf_temp( new CharmmPeriodicFunc( (i->second).key3() /*x0*/, (i->second).key1() /*n*/, (i->second).key2() /*k*/ ) );
					ConstraintOP d_temp( new DihedralConstraint( a1, a2, a3, a4, cpf_temp ) );
					pose.add_constraint( d_temp );
					TR << "Added Charmm constraint x0/n/k: " <<  (i->second).key3() << "/" << (i->second).key1() << "/" << (i->second).key2() << std::endl;
				}
			}
		} // end Scope section 2.


		{ // Scope section 3.
			Size const mmat1 = resconn_mmat1;
			Size const mmat2 = resconn_mmat2;

			utility::vector1< chemical::three_atom_set > const & rsd2_atoms_wi2_bonds_of_ii(
				rsd2_type.atoms_within_two_bonds_of_a_residue_connection( resconn_id2 ));

			for ( Size jj = 1; jj <= rsd2_atoms_wi2_bonds_of_ii.size(); ++jj ) {
				debug_assert( rsd2_atoms_wi2_bonds_of_ii[ jj ].key1() == resconn_atomno2 );

				Size const jj_atom3 = rsd2_atoms_wi2_bonds_of_ii[ jj ].key2();
				Size const mmat3 = rsd2_type.atom( jj_atom3 ).mm_atom_type_index();

				Size const jj_atom4 = rsd2_atoms_wi2_bonds_of_ii[ jj ].key3();
				Size const mmat4 = rsd2_type.atom( jj_atom4 ).mm_atom_type_index();

				//Real const angle = numeric::dihedral_radians(
				// rsd1.xyz( resconn_atomno1 ),
				// rsd2.xyz( resconn_atomno2 ),
				// rsd2.xyz( jj_atom3 ),
				// rsd2.xyz( jj_atom4 ) );

				TR << "Section 3: "
					<< "r1 " << resconn_atomno1 << " " << rsd1.atom_name( resconn_atomno1 ) << "(" << rsd1_type.atom( resconn_atomno1 ).mm_name() << ") - "
					<< "r2 " << resconn_atomno2 << " " << rsd2.atom_name( resconn_atomno2 ) << "(" << rsd2_type.atom( resconn_atomno2 ).mm_name() << ") - "
					<< "r2 " << jj_atom3        << " " << rsd1.atom_name( jj_atom3 )        << "(" << rsd1_type.atom( jj_atom3 ).mm_name()  << ") - "
					<< "r2 " << jj_atom4        << " " << rsd1.atom_name( jj_atom4 )        << "(" << rsd1_type.atom( jj_atom4 ).mm_name()  << ")" << std::endl;

				// get mm params and add constraints for each set of params
				scoring::mm::mm_torsion_library_citer_pair pair = mm_torsion_library_.lookup( mmat1, mmat2, mmat3, mmat4 );

				id::AtomID a1( resconn_atomno1, nterm_rsd_num_ );
				id::AtomID a2( resconn_atomno2, cterm_rsd_num_ );
				id::AtomID a3( jj_atom3, cterm_rsd_num_ );
				id::AtomID a4( jj_atom4, cterm_rsd_num_ );

				for ( auto i = pair.first, e = pair.second; i != e; ++i ) {
					scoring::func::FuncOP cpf_temp( new CharmmPeriodicFunc( (i->second).key3() /*x0*/, (i->second).key1() /*n*/, (i->second).key2() /*k*/ ) );
					ConstraintOP d_temp( new DihedralConstraint( a1, a2, a3, a4, cpf_temp ) );
					pose.add_constraint( d_temp );
					TR << "Added Charmm constraint x0/n/k: " <<  (i->second).key3() << "/" << (i->second).key1() << "/" << (i->second).key2() << std::endl;
				}
			}
		} // end Scope section 3.


	}
}


/// @brief Setup the score function with the appropriate weights on the constraints.
/// If user has provided a score function the and set the weights to a specific value,
/// that value is kept. If weight is not set, prints a warning.
void
CyclizationMover::setup_scorefunction()
{
	using namespace core;
	using namespace scoring;

	if ( score_fxn_ == nullptr ) {
		TR << "Creating score function and setting geometric constraint weights to 1" << std::endl;
		score_fxn_ = get_score_function();
		score_fxn_->set_weight( atom_pair_constraint, 1 );
		score_fxn_->set_weight( angle_constraint, 1 );
		score_fxn_->set_weight( dihedral_constraint, 10 );
	} else if ( score_fxn_->get_weight( atom_pair_constraint ) == 0 ) {
		TR.Warning << "atom_pair_constraint weight set to zero. Cyclization constraints will not work properly." << std::endl;
	} else if ( score_fxn_->get_weight( angle_constraint ) == 0 ) {
		TR.Warning << "angle_constraint weight set to zero. Cyclization constraints will not work properly" << std::endl;
	} else if ( score_fxn_->get_weight( dihedral_constraint ) == 0 ) {
		TR.Warning << "dihedral_constraint weight set to zero. Cyclization constraints will not work properly" << std::endl;
	}
}


/// @brief
void
CyclizationMover::setup_minimizer( core::pose::Pose & pose )
{
	using namespace core;
	using namespace kinematics;

	if ( move_map_ == nullptr ) {
		TR << "Creating move map and setting all backbone and sidechain DOFs to movable" << std::endl;
		move_map_ = core::kinematics::MoveMapOP( new MoveMap() );
		move_map_->set_bb_true_range( pose.conformation().chain_begin( chain_to_cyclize_ ), pose.conformation().chain_end( chain_to_cyclize_ ) );
		move_map_->set_chi_true_range( pose.conformation().chain_begin( chain_to_cyclize_ ), pose.conformation().chain_end( chain_to_cyclize_ ) );
	}
}

/// @brief Minimize the pose using the score function, and constraints we setup
void
CyclizationMover::minimize_rebuild( core::pose::Pose & pose )
{
	using namespace core;
	using namespace protocols;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	minimization_packing::MinMoverOP min_mover( new minimization_packing::MinMover( move_map_, score_fxn_, option[ run::min_type ].value(), 0.01, true ) );

	min_mover->apply( pose );

	// rebuild the conection dependant atoms
	pose.conformation().rebuild_residue_connection_dependent_atoms( nterm_rsd_num_, 2 );
	pose.conformation().rebuild_residue_connection_dependent_atoms( cterm_rsd_num_, 2 );
}

}//simple_moves
}//protocols
