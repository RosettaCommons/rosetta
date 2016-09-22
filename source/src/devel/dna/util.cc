// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file util.cc
/// @details Miscellaneous 'public' functions related to DNA that aren't really protocols.  That is, functions here should be those that are expected to be included outside of this directory.  The rest should just be in util.  (The existence of two separate utility files is intended to minimize dependencies, and to make general use functions more obvious.)
/// @author ashworth

#include <devel/dna/util.hh>

#include <protocols/dna/util.hh>
using namespace protocols::dna;

#include <core/types.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <basic/Tracer.hh>


// Utility Headers
using utility::vector1;

// Numeric Headers
#include <numeric/xyzVector.hh>
typedef numeric::xyzVector< core::Real > xyzVec;

// ObjexxFCL Headers

// C++ headers
#include <iostream>
#include <string>
#include <list>

#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <numeric/random/random.fwd.hh>
#include <ObjexxFCL/format.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pack;
using namespace rotamer_set;
using namespace ObjexxFCL::format;

namespace devel {
namespace dna {

static THREAD_LOCAL basic::Tracer TR( "devel.dna.util" );

///////////////////////////////////////////////////////////////////////////////
void
detect_interface_residues(
	pose::Pose const & pose,
	utility::vector1< Size > const & pos_list,
	Real const contact_threshold,
	utility::vector1< bool > & interface
)
{
	Size const nres( pose.size() );

	interface.resize( nres, false );

	for ( Size j=1; j<= nres; ++j ) {
		interface[j] = false;

		conformation::Residue const & rsd2( pose.residue( j ) );
		if ( !rsd2.is_protein() ) continue;

		for ( Size n=1; n<= pos_list.size(); ++n ) {
			conformation::Residue const & rsd1( pose.residue( pos_list[n] ) );
			assert( rsd1.is_DNA() );

			if ( rsd1.xyz( rsd1.nbr_atom() ).distance( rsd2.xyz( rsd2.nbr_atom() ) ) <= contact_threshold ) {
				interface[ j ] = true;
				break;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
void
detect_interface_by_nbrs(
	pose::Pose const & pose,
	utility::vector1< bool > & interface
)
{
	Size const nres( pose.size() );

	interface.clear();
	interface.resize( nres, false );

	for ( Size j=1; j<= nres; ++j ) {
		interface[j] = false;

		conformation::Residue const & rsd2( pose.residue( j ) );
		if ( !rsd2.is_protein() ) continue;

		for ( Size i=1; i<= nres; ++i ) {
			conformation::Residue const & rsd1( pose.residue( i ) );
			if ( !rsd1.is_DNA() ) continue;

			if ( rsd1.nbr_atom_xyz().distance( rsd2.nbr_atom_xyz() ) <= rsd1.nbr_radius() + rsd2.nbr_radius() + 5.5 ) {
				interface[ j ] = true;
				break;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
void
detect_allatom_interface_residues(
	pose::Pose const & pose,
	utility::vector1< Size > const & pos_list,
	Real const contact_threshold,
	utility::vector1< bool > & interface
)
{
	Size const nres( pose.size() );

	interface.resize( nres, false );

	for ( Size j=1; j<= nres; ++j ) {
		interface[j] = false;

		conformation::Residue const & rsd2( pose.residue( j ) );
		if ( !rsd2.is_protein() ) continue;

		for ( Size n=1; n<= pos_list.size(); ++n ) {
			conformation::Residue const & rsd1( pose.residue( pos_list[n] ) );
			assert( rsd1.is_DNA() );

			for ( Size ii=1; ii<= rsd1.natoms(); ++ii ) {
				for ( Size jj=1; jj<= rsd2.natoms(); ++jj ) {

					if ( rsd1.xyz(ii).distance( rsd2.xyz(jj) ) <= contact_threshold ) {
						interface[ j ] = true;
						break;
					}
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
void
calc_protein_DNA_rmsd(
	pose::Pose const & pose,
	pose::Pose const & reference_pose,
	Real & ca_rmsd,
	Real & interface_ca_rmsd,
	Real & interface_allatom_rmsd
)
{
	Size const nres( pose.size() );

	utility::vector1< int > pos_list;
	utility::vector1< bool > interface;
	for ( Size i=1; i<= nres; ++i ) {
		if ( pose.residue(i).is_DNA() ) pos_list.push_back(i);
	}

	detect_allatom_interface_residues( reference_pose, pos_list, 4.5, interface );

	ca_rmsd = 0.0;
	interface_ca_rmsd = 0.0;
	interface_allatom_rmsd = 0.0;

	int protres(0), intres(0), intatoms(0);
	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) ), ref_rsd( reference_pose.residue(i) );
		if ( !rsd.is_protein() ) continue;

		Real const ca_dis2( rsd.xyz("CA").distance_squared( ref_rsd.xyz("CA") ) );
		++protres;
		ca_rmsd += ca_dis2;
		if ( interface[i] ) {
			++intres;
			interface_ca_rmsd += ca_dis2;
			for ( Size j=1; j<= rsd.nheavyatoms(); ++j ) {
				++intatoms;
				interface_allatom_rmsd += rsd.xyz(j).distance( ref_rsd.xyz( rsd.atom_name(j) ) );
			}
		}
	}

	ca_rmsd = std::sqrt( ca_rmsd/ protres );
	interface_ca_rmsd = std::sqrt( interface_ca_rmsd / intres );
	interface_allatom_rmsd = std::sqrt( interface_allatom_rmsd / intatoms );
}


void
calc_DNA_bb_rmsd(
	pose::Pose const & pose,
	pose::Pose const & reference_pose,
	Real & dna_bb_rmsd
)
{
	Size const nres( pose.size() );

	dna_bb_rmsd = 0.0;

	int dna_bb_atoms(0);
	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) ), ref_rsd( reference_pose.residue(i) );
		if ( !rsd.is_DNA() ) continue;

		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			if ( rsd.atom_is_backbone( j ) ) {
				dna_bb_rmsd += rsd.xyz(j).distance( ref_rsd.xyz( rsd.atom_name(j) ) );
				++dna_bb_atoms;
			}
		}
	}

	dna_bb_rmsd = std::sqrt( dna_bb_rmsd/ dna_bb_atoms );
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// @details
/// uses sasa calculation to look at packing quality and unsatisfied hbonds
///
void
analyze_interface_sasa(
	pose::Pose const & complex_pose, // since we need to update nbr info
	int const split_jump,
	Real & bsasa14,
	Real & bsasa5,
	utility::vector1< core::id::AtomID > & buried_unsatisfied_donors,
	utility::vector1< core::id::AtomID > & buried_unsatisfied_acceptors
)
{

	pose::Pose partner1, partner2;
	pose::partition_pose_by_jump( complex_pose, split_jump, partner1, partner2 );

	analyze_interface_sasa( complex_pose, partner1, partner2, bsasa14, bsasa5,
		buried_unsatisfied_donors, buried_unsatisfied_acceptors );
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// @details
/// uses sasa calculation to look at packing quality and unsatisfied hbonds


void
analyze_interface_sasa(
	pose::Pose const & complex_pose,
	pose::Pose const & partner1,
	pose::Pose const & partner2,
	Real & bsasa14,
	Real & bsasa5,
	utility::vector1< core::id::AtomID > & buried_unsatisfied_donors,
	utility::vector1< core::id::AtomID > & buried_unsatisfied_acceptors
)
{
	using scoring::calc_total_sasa;


	Real const hbond_energy_threshold( -0.01 );
	Real const burial_threshold( 0.01 );
	Real const probe_radius( 1.0 ); // less than 1.4 --> fewer things "buried" --> more conservative

	//////////////
	// buried sasa
	Real const total_sasa_complex14( calc_total_sasa( complex_pose, 1.4 ) );
	Real const total_sasa_complex5 ( calc_total_sasa( complex_pose, 0.5 ) );

	Real const total_sasa_partner1_14( calc_total_sasa( partner1, 1.4 ) );
	Real const total_sasa_partner1_5 ( calc_total_sasa( partner1, 0.5 ) );

	Real const total_sasa_partner2_14( calc_total_sasa( partner2, 1.4 ) );
	Real const total_sasa_partner2_5 ( calc_total_sasa( partner2, 1.4 ) );

	bsasa14 = ( total_sasa_partner1_14 + total_sasa_partner2_14 - total_sasa_complex14 )/2;
	bsasa5  = ( total_sasa_partner1_5  + total_sasa_partner2_5  - total_sasa_complex5  )/2;

	////////////////////////////////////
	// now for buried unsatisfied hbonds
	pose::Pose pose;
	pose = complex_pose; // since we have to perform a non-const operation, updating the nbrs

	id::AtomID_Map< Real > atom_sasa_complex, total_hbond_energy;
	utility::vector1< Real > rsd_sasa_complex;
	scoring::calc_per_atom_sasa( pose, atom_sasa_complex, rsd_sasa_complex, probe_radius, true );

	// HACK!!!!!!!!!!!!!! TO TRIGGER 10A nbr graph update
	scoring::hbonds::HBondSet hbond_set;
	{
		using namespace scoring;
		ScoreFunction sf;
		sf.set_weight( hbond_sc, 1.0 );
		sf(pose);
		pose.update_residue_neighbors();
		scoring::hbonds::fill_hbond_set( pose, false, hbond_set );
	}

	core::pose::initialize_atomid_map( total_hbond_energy, pose, 0.0 );

	for ( Size i=1; i<= Size(hbond_set.nhbonds()); ++i ) {
		if ( hbond_set.allow_hbond(i) ) {
			scoring::hbonds::HBond const & hb( hbond_set.hbond(i) );
			id::AtomID const hatm( hb.don_hatm(), hb.don_res() );
			id::AtomID const aatm( hb.acc_atm() , hb.acc_res() );
			total_hbond_energy[ hatm ] += hb.energy(); // unweighted
			total_hbond_energy[ aatm ] += hb.energy(); // unweighted
		}
	}

	// now find unsat+buried
	buried_unsatisfied_donors.clear();
	buried_unsatisfied_acceptors.clear();

	for ( Size i=1; i<= pose.size(); ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );

		// donors
		for ( auto
				hnum  = rsd.Hpos_polar().begin(),
				hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
			id::AtomID const hatm( *hnum, i );
			if ( total_hbond_energy[ hatm ] >= hbond_energy_threshold &&
					atom_sasa_complex[ hatm ] <= burial_threshold ) {
				TR << "UNSAT_DONOR " << I(4,i) << ' ' << rsd.name3() << ' ' << rsd.atom_name(*hnum) << std::endl;
				buried_unsatisfied_donors.push_back( hatm );
				//++buried_unsatisfied_donors;
			}
		}

		// acceptors
		for ( auto
				anum  = rsd.accpt_pos().begin(),
				anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {
			id::AtomID const aatm( *anum, i );
			if ( total_hbond_energy[ aatm ] >= hbond_energy_threshold &&
					atom_sasa_complex[ aatm ] <= burial_threshold ) {
				TR << "UNSAT_ACCEPTOR " << I(4,i) << ' ' << rsd.name3() << ' ' << rsd.atom_name(*anum) << std::endl;
				buried_unsatisfied_acceptors.push_back( aatm );
				//++buried_unsatisfied_acceptors;
			}
		}
	} // i

}

//////////////////////////////////////////////////////////////////////////////////////////////
void
make_base_pair_mutation(
	pose::Pose & pose,
	Size const seqpos,
	chemical::AA const & na
)
{
	using namespace chemical;
	using namespace conformation;
	using namespace scoring::dna;

	ResidueTypeSetCOP residue_set( pose.residue(1).residue_type_set() );
	BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );

	for ( int r=1; r<= 2; ++r ) {
		Size const pos( r == 1 ? seqpos : partner[seqpos] );
		if ( pos == 0 ) continue; // unpaired
		AA const aa( r == 1 ? na : dna_base_partner( na ) );

		Residue const & existing_residue( pose.residue( pos ) );
		assert( existing_residue.is_DNA() );

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


//////////////////////////////////////////////////////////////////////////////////////////////
void
randomize_motif_sequence(
	pose::Pose & pose,
	Size const motif_begin,
	Size const motif_size
)
{
	using namespace chemical;

	for ( Size i=motif_begin; i< motif_begin+motif_size; ++i ) {
		AA const random_na( AA( first_DNA_aa + static_cast< int >( numeric::random::uniform() * 4 ) ) );
		assert( random_na >= first_DNA_aa && random_na <= last_DNA_aa );

		make_base_pair_mutation( pose, i, random_na );
	}
}


void
check_residue_proximity_to_dna(
	core::Size const ppos,
	Positions const & dna_design_positions,
	core::pose::Pose const & pose,
	bool & close,  // output
	bool & contact // output
)
{
	Real const close_threshold2( 10 * 10 ), contact_threshold2( 3.4 * 3.4 ), z_cutoff( 3.0 ); // z_cutoff should be higher?

	check_residue_proximity_to_dna( ppos, dna_design_positions, pose, close_threshold2, contact_threshold2, z_cutoff,
		close, contact );
}


void
check_residue_proximity_to_dna(
	core::Size const ppos,
	Positions const & dna_design_positions,
	core::pose::Pose const & pose,
	core::Real const close_threshold2,    // A**2
	core::Real const contact_threshold2,  // A**2
	core::Real const z_cutoff,             // A
	bool & close,  // output
	bool & contact // output
)
{
	close = false;
	contact = false;

	Residue const & res( pose.residue( ppos ) );
	assert( res.is_protein() );

	Real shortest_arg_dis2(10000), dis2;

	for ( core::Size dna_design_position : dna_design_positions ) {
		Residue const & dres( pose.residue(dna_design_position) );
		close = close || close_to_dna( res, dres, close_threshold2, true );
		if ( close ) {
			if ( z_axis_dist( res, dres ) < z_cutoff ) {
				dis2 = argrot_dna_dis2( pose, ppos, res, dres, contact_threshold2, true );
				if ( dis2 < shortest_arg_dis2 ) shortest_arg_dis2 = dis2;
				if ( shortest_arg_dis2 < contact_threshold2 ) contact = true;
			}
		}
	}
}

}
}
