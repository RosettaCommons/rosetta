// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/MMBondLengthAnergy.cc
/// @brief  Molecular mechanics bond length score class
/// @author Frank DiMaio (based on Colin Smith's MMBondAngle potential)

// Unit headers
#include <core/scoring/methods/MMBondLengthEnergy.hh>
#include <core/scoring/methods/MMBondLengthEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Utility headers
// AUTO-REMOVED #include <utility/string_util.hh>

// Numeric headers
// AUTO-REMOVED #include <numeric/xyz.functions.hh>
#include <numeric/deriv/distance_deriv.hh>

// C++ headers
#include <iostream>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the MMBondLengthEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
MMBondLengthEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return new MMBondLengthEnergy( options );
}

ScoreTypes
MMBondLengthEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( mm_stretch );
	return sts;
}


static thread_local basic::Tracer TR( "core.mm.MMBondLengthEnergy" );

MMBondLengthEnergy::MMBondLengthEnergy( methods::EnergyMethodOptions const & /*options*/):
	parent( methods::EnergyMethodCreatorOP( new MMBondLengthEnergyCreator ) )
{}

MMBondLengthEnergy::MMBondLengthEnergy( MMBondLengthEnergy const & src ):
	parent( src )
{}

MMBondLengthEnergy::~MMBondLengthEnergy() {}


/// clone
EnergyMethodOP
MMBondLengthEnergy::clone() const
{
	return new MMBondLengthEnergy( *this );
}

///
void
MMBondLengthEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const & ) const
{
	pose.update_residue_neighbors();
	pose.update_actcoords();
}

///
void
MMBondLengthEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
	pose.update_actcoords();
}

///
void
MMBondLengthEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
	pose.update_actcoords();
}

///
bool
MMBondLengthEnergy::defines_intrares_energy( EnergyMap const & ) const
{
	return true;
}

void
MMBondLengthEnergy::residue_pair_energy(
   conformation::Residue const & rsd1,
	 conformation::Residue const & rsd2,
	 pose::Pose const & /*pose*/,
	 ScoreFunction const & ,
	 EnergyMap & emap
) const
{
	// bail out if the residues aren't bonded
	if (!rsd1.is_bonded(rsd2)) return;

	Real energy = 0;

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

		// check for vrt
		if ( rsd1_type.atom_type(resconn_atomno1).is_virtual() || rsd2_type.atom_type(resconn_atomno2).is_virtual() )
			continue;

		TR(basic::t_trace) << "Found residue connection id " << resconn_id1 << "-" << resconn_id2 << ": "
			<< rsd1.atom_name( resconn_atomno1 ) << "-" << rsd2.atom_name( resconn_atomno2 ) << std::endl;

		Size const resconn_mmat1 = rsd1_type.atom( resconn_atomno1 ).mm_atom_type_index();
		Size const resconn_mmat2 = rsd2_type.atom( resconn_atomno2 ).mm_atom_type_index();

		Real const d = 
			( rsd1.atom( resconn_atomno1 ).xyz()-rsd2.atom( resconn_atomno2 ).xyz() ).length();

		TR(basic::t_trace)
			<< "r1 " << resconn_atomno1 << " " << rsd1.atom_name( resconn_atomno1 ) << "("
			<< rsd1_type.atom( resconn_atomno1 ).mm_name() << ") - "
			<< "r2 " << resconn_atomno2 << " " << rsd2.atom_name( resconn_atomno2 ) << "("
			<< rsd2_type.atom( resconn_atomno2 ).mm_name() << ")" << std::endl;
 
		energy += potential_.mm::MMBondLengthScore::score
			(  mm::mm_bondlength_atom_pair( resconn_mmat1, resconn_mmat2 ), d );
	}
	emap[ mm_stretch ] += energy;
}

void
MMBondLengthEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap & emap
) const {
	Real energy = 0;

	// get residue type
	chemical::ResidueType const & rsd_type = rsd.type();

	TR(basic::t_trace) << "MMIntrares Processing residue " << rsd.seqpos() << " " << rsd_type.name()
		<< std::endl;

	// for each bonded atom
	for (Size atm_i=1; atm_i<=rsd_type.natoms(); ++atm_i) {
		chemical::AtomIndices atm_nbrs = rsd_type.nbrs( atm_i );
		for (Size j=1; j<=atm_nbrs.size(); ++j) {
			Size atm_j = atm_nbrs[j];
			if ( atm_i<atm_j ) { // only score each bond once -- use restype index to define ordering
				int mmat1 = rsd_type.atom( atm_i ).mm_atom_type_index();
				int mmat2 = rsd_type.atom( atm_j ).mm_atom_type_index();

				// check for vrt
				if ( rsd_type.atom_type(atm_i).is_virtual() || rsd_type.atom_type(atm_j).is_virtual() )
					continue;

				Real const d = ( rsd.atom( atm_i ).xyz()-rsd.atom( atm_j ).xyz() ).length();
	
				TR(basic::t_trace)
					<< "r1 " << atm_i << " " << rsd.atom_name( atm_i ) << "("
					<< rsd_type.atom( atm_i ).mm_name() << ") - "
					<< "r2 " << atm_j << " " << rsd.atom_name( atm_j ) << "("
					<< rsd_type.atom( atm_j ).mm_name() << ")" << std::endl;

				// score bond angle
				energy += potential_.mm::MMBondLengthScore::score
					(  mm::mm_bondlength_atom_pair( mmat1, mmat2 ), d );
			}
		}
	}

	// add energy to emap
	emap[ mm_stretch ] += energy;
}

void
MMBondLengthEnergy::eval_atom_derivative(
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

	Vector LF1( 0.0 );
	Vector LF2( 0.0 );

	core::chemical::ResidueType const & rsd_type( pose.residue_type( id.rsd() ) );
	core::conformation::Residue const & rsd( pose.residue( id.rsd() ) );
	Size const atomno( id.atomno());

	// 1 intrares bonds
	chemical::AtomIndices atm_nbrs = rsd_type.nbrs( atomno );
	for (Size j=1; j<=atm_nbrs.size(); ++j) {
		Size atm2 = atm_nbrs[j];

		int mmat1 = rsd_type.atom( atomno ).mm_atom_type_index();
		int mmat2 = rsd_type.atom( atm2 ).mm_atom_type_index();

		// check for vrt
		if ( rsd_type.atom_type(atomno).is_virtual() || rsd_type.atom_type(atm2).is_virtual() )
			continue;

		Vector f1(0.0), f2(0.0);
		Real d=0;
		numeric::deriv::distance_f1_f2_deriv( rsd.xyz( atomno ), rsd.xyz( atm2 ), d, f1, f2 );

		Real dE_dd = weights[ mm_stretch ] * potential_.mm::MMBondLengthScore::dscore(
			mm::mm_bondlength_atom_pair( mmat1, mmat2 ), d );

		LF1 += dE_dd * f1;
		LF2 += dE_dd * f2;
	}

	// 2 interres bonds
	utility::vector1< Size > const & connections( rsd_type.residue_connections_for_atom( atomno ) );
	for ( Size ii = 1; ii <= connections.size(); ++ii ) {
		Size const ii_resconn   = connections[ ii ];

		Size const ii_neighb      = rsd.residue_connection_partner( ii_resconn );
		Size const neighb_resconn = rsd.residue_connection_conn_id( ii_resconn );
		conformation::Residue const & neighb_res( pose.residue( ii_neighb ));
		chemical::ResidueType const & neighb_restype( pose.residue_type( ii_neighb ) );
		Size const neighb_atom1    = neighb_restype.residue_connection( neighb_resconn ).atomno();

		Size const mmat1 = rsd_type.atom( atomno ).mm_atom_type_index();
		Size const mmat2 = neighb_restype.atom( neighb_atom1 ).mm_atom_type_index();

		// check for vrt
		if ( rsd_type.atom_type(atomno).is_virtual() || neighb_restype.atom_type(neighb_atom1).is_virtual() )
			continue;

		Vector f1(0.0), f2(0.0);
		Real d=0;
		numeric::deriv::distance_f1_f2_deriv( rsd.xyz( atomno ), neighb_res.xyz( neighb_atom1 ), d, f1, f2 );

		Real dE_dd = weights[ mm_stretch ] * potential_.mm::MMBondLengthScore::dscore(
			mm::mm_bondlength_atom_pair( mmat1, mmat2 ), d );

		LF1 += dE_dd * f1;
		LF2 += dE_dd * f2;
	}

	F1 += LF1;
	F2 += LF2;
}


/// @brief MMBondLengthEnergy does not have an atomic interation threshold
Distance
MMBondLengthEnergy::atomic_interaction_cutoff() const
{
	return 0.0;
}

/// @brief MMBondLengthEnergy is context independent; indicates that no context graphs are required
void
MMBondLengthEnergy::indicate_required_context_graphs(utility::vector1< bool > & ) const
{}

core::Size
MMBondLengthEnergy::version() const
{
	return 1; // Initial versioning
}

} // namespace methods
} // namespace scoring
} // namespace core
