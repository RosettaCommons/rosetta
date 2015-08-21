// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunction.hh
/// @brief  molecular mechanics torsion energy
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit headers
#include <core/scoring/methods/MMTorsionEnergy.hh>
#include <core/scoring/methods/MMTorsionEnergyCreator.hh>
#include <core/scoring/mm/MMTorsionScore.hh>

// Project headers
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
//#include <core/scoring/ScoringManager.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>
#include <numeric/deriv/dihedral_deriv.hh>

// C++ headers
#include <iostream>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the MMTorsionEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
MMTorsionEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MMTorsionEnergy );
}

ScoreTypes
MMTorsionEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( mm_twist );
	return sts;
}


static thread_local basic::Tracer TR( "core.scoring.methods.MMTorsionEnergy" );


typedef std::pair< mm::mm_torsion_atom_quad, core::Real > mm_torsion_atom_quad_angle_pair;
typedef utility::vector1< mm_torsion_atom_quad_angle_pair >::const_iterator mmtaqap_iter;

MMTorsionEnergy::MMTorsionEnergy() :
	parent( methods::EnergyMethodCreatorOP( new MMTorsionEnergyCreator ) )
{}

/// clone
EnergyMethodOP
MMTorsionEnergy::clone() const
{
	return EnergyMethodOP( new MMTorsionEnergy );
}


void
MMTorsionEnergy::setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const
{
	pose.update_residue_neighbors();
}


void
MMTorsionEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}


void
MMTorsionEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}


bool
MMTorsionEnergy::defines_intrares_energy( EnergyMap const & ) const
{
	return true;
}

void
MMTorsionEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap & emap
) const
{
	Real energy = 0;
	/*
	// if not sequential
	if( rsd2.seqpos() != rsd1.seqpos() + 1 ) {
	// add energy to emap
	emap[ mm_twist ] += energy;
	} else {
	// get residue types
	chemical::ResidueType const & rsd1_type = rsd1.type();
	chemical::ResidueType const & rsd2_type = rsd2.type();

	// get residue 1 upper connect (r1uc) and residue 2 lower connect (r2lc)
	Size r1uc = rsd1_type.upper_connect_atom();
	Size r2lc = rsd2_type.lower_connect_atom();

	// generate list of new dihedrals about the bond formed between the r1uc and the r2lc
	// as well as anything with a path distance of 1 from either the r1uc (r1ca) and r2lc (r2ca)
	utility::vector1< mm_torsion_atom_quad_angle_pair > mm_torsion_atom_quad_angle_pairs_;

	// get all the atoms that are one away from the r1uc
	utility::vector1< Size > r1ucd1;
	for ( Size i = 1; i <= rsd1_type.natoms(); ++i ) {
	if ( ( rsd1_type.path_distance( r1uc, i ) == 1 ) ) {
	r1ucd1.push_back( i );
	}
	}

	// get all the atoms that are one away from the r2lc
	utility::vector1< Size > r2lcd1;
	for ( Size j = 1; j <= rsd2_type.natoms(); ++j ) {
	if ( rsd2_type.path_distance( r2lc, j ) == 1  ) {
	r2lcd1.push_back( j );
	}
	}

	// torsions about the bonds imediatly preceding r1uc and r2lc bond
	for ( utility::vector1< Size >::iterator r1ca = r1ucd1.begin(); r1ca != r1ucd1.end(); ++r1ca ) {
	for ( Size r1ta = 1; r1ta <= rsd1_type.natoms(); ++r1ta ) {
	if ( ( rsd1_type.path_distance( *r1ca, r1ta ) == 1 ) && r1ta != r1uc) {
	// get mm atom type index
	Size mmat1 = rsd1_type.mm_atom_type_index( r1ta );
	Size mmat2 = rsd1_type.mm_atom_type_index(*r1ca );
	Size mmat3 = rsd1_type.mm_atom_type_index( r1uc );
	Size mmat4 = rsd2_type.mm_atom_type_index( r2lc );

	// get angle
	Real angle = numeric::dihedral_radians
	( rsd1.atom( r1ta ).xyz(),
	rsd1.atom(*r1ca ).xyz(),
	rsd1.atom( r1uc ).xyz(),
	rsd2.atom( r2lc ).xyz() );

	// make pair
	mm_torsion_atom_quad_angle_pair temp( mm::mm_torsion_atom_quad( mmat1, mmat2, mmat3, mmat4 ), angle );

	// add to list
	mm_torsion_atom_quad_angle_pairs_.push_back( temp );
	}
	}
	}

	// torsions about the r1uc and r2lc bond
	for ( utility::vector1< Size >::iterator terminal_atom1 = r1ucd1.begin();
	terminal_atom1 != r1ucd1.end(); ++terminal_atom1 ) {
	for ( utility::vector1< Size >::iterator terminal_atom2 = r2lcd1.begin();
	terminal_atom2 != r2lcd1.end(); ++terminal_atom2 ) {
	// get mm atom type index
	Size mmat1 = rsd1_type.mm_atom_type_index( *terminal_atom1 );
	Size mmat2 = rsd1_type.mm_atom_type_index( r1uc );
	Size mmat3 = rsd2_type.mm_atom_type_index( r2lc );
	Size mmat4 = rsd2_type.mm_atom_type_index( *terminal_atom2 );

	// get angle
	Real angle = numeric::dihedral_radians
	( rsd1.atom( *terminal_atom1 ).xyz(),
	rsd1.atom( r1uc ).xyz(),
	rsd2.atom( r2lc ).xyz(),
	rsd2.atom( *terminal_atom2 ).xyz() );

	// make pair
	mm_torsion_atom_quad_angle_pair temp( mm::mm_torsion_atom_quad( mmat1, mmat2, mmat3, mmat4 ), angle );

	// add to list
	mm_torsion_atom_quad_angle_pairs_.push_back( temp );
	}
	}

	// torsions about the bonds imediatly following r1uc and r2lc bond
	for ( utility::vector1< Size >::iterator r2ca = r2lcd1.begin(); r2ca != r2lcd1.end(); ++r2ca ) {
	for ( Size r2ta = 1; r2ta <= rsd2_type.natoms(); ++r2ta ) {
	if ( ( rsd2_type.path_distance( *r2ca, r2ta ) == 1 ) && r2ta != r2lc ) {
	// get mm atom type index
	Size mmat1 = rsd1_type.mm_atom_type_index( r1uc );
	Size mmat2 = rsd2_type.mm_atom_type_index( r2lc );
	Size mmat3 = rsd2_type.mm_atom_type_index(*r2ca );
	Size mmat4 = rsd2_type.mm_atom_type_index( r2ta );

	// get angle
	Real angle = numeric::dihedral_radians
	( rsd1.atom( r1uc ).xyz(),
	rsd2.atom( r2lc ).xyz(),
	rsd2.atom(*r2ca ).xyz(),
	rsd2.atom( r2ta ).xyz() );

	// make pair
	mm_torsion_atom_quad_angle_pair temp( mm::mm_torsion_atom_quad( mmat1, mmat2, mmat3, mmat4 ), angle );

	// add to list
	mm_torsion_atom_quad_angle_pairs_.push_back( temp );
	}
	}
	}

	// score each pair
	for ( Size i = 1; i <= mm_torsion_atom_quad_angle_pairs_.size(); ++i )
	{
	// score dihedral
	energy += potential_.core::scoring::mm::MMTorsionScore::score
	( mm_torsion_atom_quad_angle_pairs_[ i ].first, mm_torsion_atom_quad_angle_pairs_[ i ].second );
	}

	// add energy to emap
	emap[ mm_twist ] += energy;
	}
	*/

	if ( ! rsd1.is_bonded( rsd2 ) ) return;

	// get residue types
	chemical::ResidueType const & rsd1_type = rsd1.type();
	chemical::ResidueType const & rsd2_type = rsd2.type();

	TR(basic::t_trace) << "residue_pair_energy: processing residues "
		<< rsd1.seqpos() << "." << rsd1_type.name() << "-"
		<< rsd2.seqpos() << "." << rsd2_type.name() << std::endl;

	utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );

	for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {

		Size const resconn_id1( r1_resconn_ids[ii] );
		Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

		Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
		Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

		TR(basic::t_trace) << "Found residue connection id " << resconn_id1 << "-" << resconn_id2 << ": "
			<< rsd1.atom_name( resconn_atomno1 ) << "-" << rsd2.atom_name( resconn_atomno2 ) << std::endl;

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

					Real const angle = numeric::dihedral_radians(
						rsd1.xyz( jj_term_atomno ),
						rsd1.xyz( resconn_atomno1 ),
						rsd2.xyz( resconn_atomno2 ),
						rsd2.xyz( kk_term_atomno ) );

					TR(basic::t_trace)
						<< "r1 " << jj_term_atomno << " " << rsd1.atom_name( jj_term_atomno ) << "("
						<< rsd1_type.atom( jj_term_atomno ).mm_name() << ") - "
						<< "r1 " << resconn_atomno1 << " " << rsd1.atom_name( resconn_atomno1 ) << "("
						<< rsd1_type.atom( resconn_atomno1 ).mm_name() << ") - "
						<< "r2 " << resconn_atomno2 << " " << rsd2.atom_name( resconn_atomno2 ) << "("
						<< rsd2_type.atom( resconn_atomno2 ).mm_name() << ") - "
						<< "r2 " << kk_term_atomno << " " << rsd2.atom_name( kk_term_atomno ) << "("
						<< rsd2_type.atom( kk_term_atomno ).mm_name() << ")" << std::endl;

					energy += potential_.core::scoring::mm::MMTorsionScore::score(
						mm::mm_torsion_atom_quad( mmat1, mmat2, mmat3, mmat4 ), angle );

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

				Real const angle = numeric::dihedral_radians(
					rsd1.xyz( jj_atom1 ),
					rsd1.xyz( jj_atom2 ),
					rsd1.xyz( resconn_atomno1 ),
					rsd2.xyz( resconn_atomno2 ) );

				energy += potential_.core::scoring::mm::MMTorsionScore::score(
					mm::mm_torsion_atom_quad( mmat1, mmat2, mmat3, mmat4 ), angle );


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

				Real const angle = numeric::dihedral_radians(
					rsd1.xyz( resconn_atomno1 ),
					rsd2.xyz( resconn_atomno2 ),
					rsd2.xyz( jj_atom3 ),
					rsd2.xyz( jj_atom4 ) );

				energy += potential_.core::scoring::mm::MMTorsionScore::score(
					mm::mm_torsion_atom_quad( mmat1, mmat2, mmat3, mmat4 ), angle );


			}
		} // end Scope section 3.
	}

	emap[ mm_twist ] += energy;

}

void
MMTorsionEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	Real energy = 0;

	// get residue type
	chemical::ResidueType const & rsd_type = rsd.type();

	TR(basic::t_trace) << "Intrares Processing residue " << rsd.seqpos() << " " << rsd_type.name()
		<< std::endl;

	// for each dihedral angle in the residue type
	for ( Size dihe = 1; dihe <= rsd_type.ndihe(); ++dihe ) {
		// get ResidueType ints
		int rt1 = ( rsd_type.dihedral( dihe ) ).key1();
		int rt2 = ( rsd_type.dihedral( dihe ) ).key2();
		int rt3 = ( rsd_type.dihedral( dihe ) ).key3();
		int rt4 = ( rsd_type.dihedral( dihe ) ).key4();

		// get mm atom type index
		int mmat1 = rsd_type.atom( rt1 ).mm_atom_type_index();
		int mmat2 = rsd_type.atom( rt2 ).mm_atom_type_index();
		int mmat3 = rsd_type.atom( rt3 ).mm_atom_type_index();
		int mmat4 = rsd_type.atom( rt4 ).mm_atom_type_index();

		// get angle
		Real angle = numeric::dihedral_radians
			( rsd.atom( rt1 ).xyz(), rsd.atom( rt2 ).xyz(),
			rsd.atom( rt3 ).xyz(), rsd.atom( rt4 ).xyz() );

		TR(basic::t_trace) << rt1 << " " << rsd.atom_name( rt1 ) << "(" << rsd_type.atom( rt1 ).mm_name()
			<< ") - " << rt2 << " " << rsd.atom_name( rt2 ) << "(" << rsd_type.atom( rt2 ).mm_name()
			<< ") - " << rt3 << " " << rsd.atom_name( rt3 ) << "(" << rsd_type.atom( rt3 ).mm_name()
			<< ") - " << rt4 << " " << rsd.atom_name( rt4 ) << "(" << rsd_type.atom( rt4 ).mm_name()
			<< ") angle " << angle << std::endl;

		// score dihedral
		energy += potential_.core::scoring::mm::MMTorsionScore::score
			(  mm::mm_torsion_atom_quad( mmat1, mmat2, mmat3, mmat4 ), angle );
	}

	// add energy to emap
	emap[ mm_twist ] += energy;
}

void
MMTorsionEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	//TR << "Evaluating derivatives for rsd# " << id.rsd() << " atom# " << id.atomno() << std::endl;

	Vector LF1( 0.0 ); // Local F1 -- for debuggin'
	Vector LF2( 0.0 ); // Local F2 -- for debuggin'

	//Vector f1(0.0), f2(0.0);
	/// 1. Evaluate intraresidue torsion angle derivatives.
	///    for torsions where this atom is scored
	/// 2. Evaluate interresidue torsion angle derivatives.
	///    KNOWN BUG: This does not evaluate derivatives for torsion angles
	///    spanning three or more residues.
	///    a.  Calculate derivatives when exactly three of the atoms come from id.rsd()
	///    b.  Calculate derivatives when exactly two of the atoms come from id.rsd()
	///    c.  Calculate derivatives when exactly one of the atoms comes from id.rsd()

	core::chemical::ResidueType const & restype( pose.residue_type( id.rsd() ) );
	core::conformation::Residue const & res( pose.residue( id.rsd() ) );
	Size const atomno( id.atomno());

	/// 1.
	utility::vector1< Size > const & diheds( restype.dihedrals_for_atom( id.atomno() ) );

	for ( Size ii = 1, ii_end = diheds.size(); ii <= ii_end; ++ii ) {
		chemical::dihedral_atom_set const & ii_dihed( restype.dihedral( diheds[ ii ] ) );
		debug_assert( ii_dihed.key1() == atomno || ii_dihed.key2() == atomno || ii_dihed.key3() == atomno  || ii_dihed.key4() == atomno );

		//TR << "  1. Deriv for angle# " << angs[ ii ] << " between " << ii_bangle.key1() << " " << ii_bangle.key2()  << " " << ii_bangle.key3() << std::endl;

		Size const mmat1 = restype.atom( ii_dihed.key1()).mm_atom_type_index();
		Size const mmat2 = restype.atom( ii_dihed.key2()).mm_atom_type_index();
		Size const mmat3 = restype.atom( ii_dihed.key3()).mm_atom_type_index();
		Size const mmat4 = restype.atom( ii_dihed.key4()).mm_atom_type_index();

		Vector f1(0.0), f2(0.0);
		Real theta(0.0);
		if ( ii_dihed.key1() == atomno ) {
			numeric::deriv::dihedral_p1_cosine_deriv(
				res.xyz( ii_dihed.key1() ),
				res.xyz( ii_dihed.key2() ),
				res.xyz( ii_dihed.key3() ),
				res.xyz( ii_dihed.key4() ),
				theta, f1, f2 );
		} else if ( ii_dihed.key2() == atomno ) {
			numeric::deriv::dihedral_p2_cosine_deriv(
				res.xyz( ii_dihed.key1() ),
				res.xyz( ii_dihed.key2() ),
				res.xyz( ii_dihed.key3() ),
				res.xyz( ii_dihed.key4() ),
				theta, f1, f2 );
		} else if ( ii_dihed.key3() == atomno ) {
			numeric::deriv::dihedral_p2_cosine_deriv(
				res.xyz( ii_dihed.key4() ),
				res.xyz( ii_dihed.key3() ),
				res.xyz( ii_dihed.key2() ),
				res.xyz( ii_dihed.key1() ),
				theta, f1, f2 );
		} else {
			numeric::deriv::dihedral_p1_cosine_deriv(
				res.xyz( ii_dihed.key4() ),
				res.xyz( ii_dihed.key3() ),
				res.xyz( ii_dihed.key2() ),
				res.xyz( ii_dihed.key1() ),
				theta, f1, f2 );
		}

		Real const dE_dtheta = weights[ mm_twist ] * potential_.core::scoring::mm::MMTorsionScore::dscore(
			mm::mm_torsion_atom_quad( mmat1, mmat2, mmat3, mmat4 ), theta );
		//TR << "    theta " << numeric::conversions::degrees( theta ) << " dE_dtheta " << dE_dtheta << " f1: " << f1.x() << " " << f1.y() << " " << f1.z() << " f2: " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;


		LF1 += dE_dtheta * f1;
		LF2 += dE_dtheta * f2;

	}

	{ // Scope 2a.  Bond angles involving three atoms on this residue.
		utility::vector1< std::pair< Size, Size > > const & interres_wi2_for_this_atom(
			restype.within2bonds_sets_for_atom( atomno ));
		for ( Size ii = 1; ii <= interres_wi2_for_this_atom.size(); ++ii ) {
			Size const ii_resconn   = interres_wi2_for_this_atom[ ii ].first;
			Size const ii_whichpair = interres_wi2_for_this_atom[ ii ].second;

			chemical::three_atom_set const & ii_triple = restype.
				atoms_within_two_bonds_of_a_residue_connection( ii_resconn )[ ii_whichpair ];
			debug_assert( ii_triple.key1() == atomno || ii_triple.key2() == atomno || ii_triple.key3() == atomno );

			/// Find the neighbor residue and atom
			Size const ii_neighb      = res.residue_connection_partner( ii_resconn );
			Size const neighb_resconn = res.residue_connection_conn_id( ii_resconn );
			conformation::Residue const & neighb_res( pose.residue( ii_neighb ));
			chemical::ResidueType const & neighb_restype( pose.residue_type( ii_neighb ) );
			Size const neighb_atom    = neighb_restype.residue_connection( neighb_resconn ).atomno();


			//TR << "  2a. Deriv for interres angle between (" << ii_neighb << "," << neighb_atom << ") " << ii_pair.key1()  << " " << ii_pair.key2() << std::endl;

			Size const mmat1 = neighb_restype.atom( neighb_atom ).mm_atom_type_index();
			Size const mmat2 = restype.atom( ii_triple.key1()).mm_atom_type_index();
			Size const mmat3 = restype.atom( ii_triple.key2()).mm_atom_type_index();
			Size const mmat4 = restype.atom( ii_triple.key3()).mm_atom_type_index();


			Vector f1(0.0), f2(0.0);
			Real theta(0.0);
			if ( ii_triple.key1() == atomno ) {
				numeric::deriv::dihedral_p2_cosine_deriv(
					neighb_res.xyz( neighb_atom ),
					res.xyz( ii_triple.key1() ),
					res.xyz( ii_triple.key2() ),
					res.xyz( ii_triple.key3() ),
					theta, f1, f2 );
			} else if ( ii_triple.key2() == atomno ) {
				numeric::deriv::dihedral_p2_cosine_deriv(
					res.xyz( ii_triple.key3() ),
					res.xyz( ii_triple.key2() ),
					res.xyz( ii_triple.key1() ),
					neighb_res.xyz( neighb_atom ),
					theta, f1, f2 );
			} else {
				numeric::deriv::dihedral_p1_cosine_deriv(
					res.xyz( ii_triple.key3() ),
					res.xyz( ii_triple.key2() ),
					res.xyz( ii_triple.key1() ),
					neighb_res.xyz( neighb_atom ),
					theta, f1, f2 );
			}
			Real const dE_dtheta = weights[ mm_twist ] * potential_.core::scoring::mm::MMTorsionScore::dscore(
				mm::mm_torsion_atom_quad( mmat1, mmat2, mmat3, mmat4 ), theta );

			//TR << "    theta " << numeric::conversions::degrees( theta ) << " dE_dtheta " << dE_dtheta << " f1: " << f1.x() << " " << f1.y() << " " << f1.z() << " f2: " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;

			LF1 += dE_dtheta * f1;
			LF2 += dE_dtheta * f2;
		}
	} // end scope 2a

	/// b. Bond angles involving two atoms on this residue.
	/// Iterate across residue connections for this atom -- for each residue connection,
	/// find the residue connection partner, and the list of within-1-bond atom pairs
	/// for the atom forming the residue connection.  Then iterate across this list.
	{ // Scope 2b.
		utility::vector1< std::pair< Size, Size > > const & interres_wi1_for_this_atom(
			restype.within1bonds_sets_for_atom( atomno ));
		for ( Size ii = 1; ii <= interres_wi1_for_this_atom.size(); ++ii ) {
			Size const ii_resconn   = interres_wi1_for_this_atom[ ii ].first;
			Size const ii_whichpair = interres_wi1_for_this_atom[ ii ].second;

			chemical::two_atom_set const & ii_pair = restype.
				atoms_within_one_bond_of_a_residue_connection( ii_resconn )[ ii_whichpair ];
			debug_assert( ii_pair.key1() == atomno || ii_pair.key2() == atomno );


			Size const mmat1 = restype.atom( ii_pair.key2()).mm_atom_type_index();
			Size const mmat2 = restype.atom( ii_pair.key1()).mm_atom_type_index();

			/// Find the neighbor residue and connection atom
			Size const ii_neighb      = res.residue_connection_partner( ii_resconn );
			Size const neighb_resconn = res.residue_connection_conn_id( ii_resconn );
			conformation::Residue const & neighb_res( pose.residue( ii_neighb ));
			chemical::ResidueType const & neighb_restype( pose.residue_type( ii_neighb ) );
			Size const neighb_atom1   = neighb_restype.residue_connection( neighb_resconn ).atomno();
			Size const mmat3          = neighb_restype.atom( neighb_atom1 ).mm_atom_type_index();

			utility::vector1< chemical::two_atom_set > const & neighb_atoms_wi1_bond_of_ii(
				neighb_restype.atoms_within_one_bond_of_a_residue_connection( neighb_resconn ));

			for ( Size jj = 1; jj <= neighb_atoms_wi1_bond_of_ii.size(); ++jj ) {
				chemical::two_atom_set neighb_pair( neighb_atoms_wi1_bond_of_ii[ jj ] );
				debug_assert( neighb_pair.key1() == neighb_atom1 );

				Size const mmat4 = neighb_restype.atom( neighb_pair.key2() ).mm_atom_type_index();

				Vector f1(0.0), f2(0.0);
				Real theta(0.0);
				if ( ii_pair.key1() == atomno ) {
					numeric::deriv::dihedral_p2_cosine_deriv(
						res.xyz( ii_pair.key2() ),
						res.xyz( ii_pair.key1() ),
						neighb_res.xyz( neighb_atom1 ),
						neighb_res.xyz( neighb_pair.key2() ),
						theta, f1, f2 );
				} else {
					numeric::deriv::dihedral_p1_cosine_deriv(
						res.xyz( ii_pair.key2() ),
						res.xyz( ii_pair.key1() ),
						neighb_res.xyz( neighb_atom1 ),
						neighb_res.xyz( neighb_pair.key2() ),
						theta, f1, f2 );
				}

				Real const dE_dtheta = weights[ mm_twist ] * potential_.core::scoring::mm::MMTorsionScore::dscore(
					mm::mm_torsion_atom_quad( mmat1, mmat2, mmat3, mmat4 ), theta );

				//TR << "    theta " << numeric::conversions::degrees( theta ) << " dE_dtheta " << dE_dtheta << " f1: " << f1.x() << " " << f1.y() << " " << f1.z() << " f2: " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;

				LF1 += dE_dtheta * f1;
				LF2 += dE_dtheta * f2;
			}
		}
	} /// end scope 2b.

	{ // Scope 2c.  All the inter-residue connections involving this atom and all within-2-bond triples
		/// of the bonded residues.
		utility::vector1< Size > const & connections( restype.residue_connections_for_atom( atomno ) );
		Size const mmat1 = restype.atom( atomno ).mm_atom_type_index();

		for ( Size ii = 1; ii <= connections.size(); ++ii ) {
			Size const ii_resconn   = connections[ ii ];

			/// Find the neighbor residue and connection atom
			Size const ii_neighb      = res.residue_connection_partner( ii_resconn );
			Size const neighb_resconn = res.residue_connection_conn_id( ii_resconn );
			conformation::Residue const & neighb_res( pose.residue( ii_neighb ));
			chemical::ResidueType const & neighb_restype( pose.residue_type( ii_neighb ) );
			Size const neighb_atom1    = neighb_restype.residue_connection( neighb_resconn ).atomno();
			Size const mmat2          = neighb_restype.atom( neighb_atom1 ).mm_atom_type_index();

			utility::vector1< chemical::three_atom_set > const & neighb_atoms_wi2(
				neighb_restype.atoms_within_two_bonds_of_a_residue_connection( neighb_resconn ));

			for ( Size jj = 1; jj <= neighb_atoms_wi2.size(); ++jj ) {

				chemical::three_atom_set const & neighb_triple = neighb_atoms_wi2[ jj ];
				debug_assert( neighb_triple.key1() == neighb_atom1 );

				Size const neighb_atom2 = neighb_triple.key2();
				Size const mmat3 = neighb_restype.atom( neighb_atom2 ).mm_atom_type_index();

				Size const neighb_atom3 = neighb_triple.key3();
				Size const mmat4 = neighb_restype.atom( neighb_atom3 ).mm_atom_type_index();

				//TR << "  2c. Deriv for interres angle between (" << ii_neighb << "," << neighb_atom2 << ") ("<< ii_neighb << "," << neighb_atom1 << ") " << atomno << std::endl;


				Vector f1(0.0), f2(0.0);
				Real theta(0.0);

				numeric::deriv::dihedral_p1_cosine_deriv(
					res.xyz( atomno ),
					neighb_res.xyz( neighb_atom1 ),
					neighb_res.xyz( neighb_atom2 ),
					neighb_res.xyz( neighb_atom3 ),
					theta, f1, f2 );

				Real const dE_dtheta = weights[ mm_twist ] * potential_.core::scoring::mm::MMTorsionScore::dscore(
					mm::mm_torsion_atom_quad( mmat1, mmat2, mmat3, mmat4 ), theta );

				//TR << "    theta " << numeric::conversions::degrees( theta ) << " dE_dtheta " << dE_dtheta << " f1: " << f1.x() << " " << f1.y() << " " << f1.z() << " f2: " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;

				LF1 += dE_dtheta * f1;
				LF2 += dE_dtheta * f2;

			}
		}
	} // end scope 2c.

	//TR << "Finished dihedral derivatives for rsd# " << id.rsd() << " atom# " << id.atomno() <<  " F1: " << LF1.x() << " " << LF1.y() << " " << LF1.z() << " F2: " << LF2.x() << " " << LF2.y() << " " << LF2.z() << std::endl;
	F1 += LF1;
	F2 += LF2;


}

/// @brief MMTorsionEnergy does not have an atomic interation threshold
Distance
MMTorsionEnergy::atomic_interaction_cutoff() const
{
	return 0.0;
}

/// @brief MMTorsionEnergy is context independent; indicates that no context graphs are required
void
MMTorsionEnergy::indicate_required_context_graphs(utility::vector1< bool > & ) const
{}
core::Size
MMTorsionEnergy::version() const
{
	return 1; // Initial versioning
}

} // namespace methods
} // namespace scoring
} // namespace core
