// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_VDW_Energy.cc
/// @brief  Statistically derived clash check for RNA.
/// @author Rhiju Das


// Unit headers
#include <core/scoring/rna/RNA_VDW_Energy.hh>
#include <core/scoring/rna/RNA_VDW_EnergyCreator.hh>

// Package headers
#include <core/scoring/rna/RNA_AtomVDW.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <basic/Tracer.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>


#include <core/id/AtomID.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <ObjexxFCL/format.hh>

using ObjexxFCL::format::I;

static THREAD_LOCAL basic::Tracer tr( "core.scoring.rna.RNA_VDW_Energy" );

namespace core {
namespace scoring {
namespace rna {


/// @details This must return a fresh instance of the RNA_VDW_Energy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_VDW_EnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new RNA_VDW_Energy );
}

ScoreTypes
RNA_VDW_EnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rna_vdw );
	return sts;
}


RNA_VDW_Energy::RNA_VDW_Energy() :
	parent( methods::EnergyMethodCreatorOP( new RNA_VDW_EnergyCreator ) ),
	rna_atom_vdw_( ScoringManager::get_instance()->get_RNA_AtomVDW() ), // need to make the table choice configurable
	vdw_scale_factor_( 0.8 ) // hack from rosetta++
{}


/// clone
methods::EnergyMethodOP
RNA_VDW_Energy::clone() const
{
	return methods::EnergyMethodOP( new RNA_VDW_Energy );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_VDW_Energy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
	setup_atom_numbers_for_vdw_calculation( pose );
}


//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_VDW_Energy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
	setup_atom_numbers_for_vdw_calculation( pose );
}


//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_VDW_Energy::setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const &  ) const
{
	pose.update_residue_neighbors();
	setup_atom_numbers_for_vdw_calculation( pose );
}


//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_VDW_Energy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	//std::cout << "Checking out VDW: " << rsd1.seqpos() << " " << rsd2.seqpos() << std::endl;

	rna::RNA_ScoringInfo const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	utility::vector1< bool > const & is_magnesium = rna_scoring_info.is_magnesium();

	Size const pos1 = rsd1.seqpos();
	Size const pos2 = rsd2.seqpos();

	if ( !rsd1.is_RNA() && !is_magnesium[ pos1 ] ) return;
	if ( !rsd2.is_RNA() && !is_magnesium[ pos2 ] ) return;

	char const which_nucleotide1 = rsd1.name1(); //a,c,g,u
	char const which_nucleotide2 = rsd2.name1(); //a,c,g,u

	Real score( 0.0 );

	//Precomputed list of atom numbers... cached inside the pose.
	utility::vector1< utility::vector1< Size > > const &
		atom_numbers_for_vdw_calculation( rna_scoring_info.atom_numbers_for_vdw_calculation() );

	utility::vector1< Size > const & atom_numbers1 ( atom_numbers_for_vdw_calculation[ pos1 ]  );

	Size const num_vdw_atoms1( atom_numbers1.size() );

	// no countpair for RNA!
	// Don't need to loop over all atoms, just some representatives...
	for ( Size m = 1; m <= num_vdw_atoms1; ++m ) {

		Size const i = atom_numbers1[ m ];
		if ( rsd1.is_virtual( i ) ) continue;

		Vector const & i_xyz( rsd1.xyz( i ) );

		utility::vector1< Size > const & atom_numbers2 ( atom_numbers_for_vdw_calculation[ pos2 ]  );
		Size const num_vdw_atoms2( atom_numbers2.size() );

		for ( Size n = 1; n <= num_vdw_atoms2; ++n ) {

			Size const j = atom_numbers2[ n ];
			if ( rsd2.is_virtual( j ) ) continue;

			Real const bump_dsq( rna_atom_vdw_.bump_parameter( m, n, which_nucleotide1, which_nucleotide2 ) );
			Real const clash( bump_dsq - i_xyz.distance_squared( rsd2.xyz( j ) ) );

			if ( clash > 0.0 ) {
				score += ( clash * clash ) / bump_dsq;
				//     tr << "BUMP " << rsd1.name3() << I(4,rsd1.seqpos() ) << ' ' << rsd1.atom_name(i) << " --- " << rsd2.name3() << I(4,rsd2.seqpos() ) << ' ' << rsd2.atom_name(j) << " Penalty: " << ( clash * clash ) / bump_dsq << " sqrt(atomvdw) " << sqrt( bump_dsq ) <<  " dist "  << i_xyz.distance( rsd2.xyz(j)) << "  index: "  << m << " " << n << std::endl;
			}
		}
	}

	emap[ rna_vdw ] += score * vdw_scale_factor_; // vdw prefactor!
}


/////////////////////////////////
Size
RNA_VDW_Energy::get_vdw_atom_number(
	utility::vector1< utility::vector1< Size > > const & atom_numbers_for_vdw_calculation,
	Size const & pos1,
	Size const & i
) const {
	Size m( 0 );
	bool is_vdw_atom( false );

	for ( m = 1; m <= atom_numbers_for_vdw_calculation.size(); ++m ) {
		if ( atom_numbers_for_vdw_calculation[ pos1 ][ m ] == i ) {
			is_vdw_atom = true; break;
		}
	}

	if ( is_vdw_atom ) return m;
	return 0;
}


/////////////////////////////////
void
RNA_VDW_Energy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {

	Size const pos1( atom_id.rsd() );
	Size const i   ( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( pos1 ) );

	rna::RNA_ScoringInfo const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	utility::vector1< bool > const & is_magnesium = rna_scoring_info.is_magnesium();

	if ( !rsd1.is_RNA() && !is_magnesium[ pos1 ]  ) return;

	char const which_nucleotide1 = rsd1.name1(); //a,c,g,u

	int const pos1_map( domain_map( pos1 ) );
	bool const pos1_fixed( pos1_map != 0 );

	Vector const & i_xyz( rsd1.xyz( i ) );

	utility::vector1< utility::vector1< Size > > const &
		atom_numbers_for_vdw_calculation( rna_scoring_info.atom_numbers_for_vdw_calculation() );

	Size m = get_vdw_atom_number( atom_numbers_for_vdw_calculation, rsd1.seqpos(), i );
	if ( m == 0 )  return;

	// cached energies object
	Energies const & energies( pose.energies() );

	// the neighbor/energy links
	EnergyGraph const & energy_graph( energies.energy_graph() );

	// loop over *all* nbrs of rsd1 (not just upper or lower)
	for ( graph::Graph::EdgeListConstIter
			iru  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
			irue = energy_graph.get_node( pos1 )->const_edge_list_end();
			iru != irue; ++iru ) {
		Size const pos2( ( *iru )->get_other_ind( pos1 ) );

		if ( pos1_fixed && pos1_map == domain_map( pos2 ) ) continue; // fixed wrt one another

		conformation::Residue const & rsd2( pose.residue( pos2 ) );
		char const which_nucleotide2 = rsd2.name1(); //a,c,g,u

		if ( !rsd2.is_RNA()  && !is_magnesium[ pos2 ] ) continue;

		runtime_assert( pos2 != pos1 );

		utility::vector1< Size > const & atom_numbers2 ( atom_numbers_for_vdw_calculation[ pos2 ]  );

		for ( Size n = 1; n <= atom_numbers2.size(); ++n ) {

			Size const j = atom_numbers2[ n ];

			Vector const & j_xyz( rsd2.xyz( j ) );
			Vector const f2( i_xyz - j_xyz );
			Real const dis2( f2.length_squared() );
			Real const bump_dsq( rna_atom_vdw_.bump_parameter( m, n, which_nucleotide1, which_nucleotide2 ) );

			if ( dis2 < bump_dsq ) {
				Real const dE_dr_over_r = vdw_scale_factor_ * weights[ rna_vdw ] * 4.0 * ( dis2 - bump_dsq ) / bump_dsq;
				Vector const f1( i_xyz.cross( j_xyz ) );
				F1 += dE_dr_over_r * f1;
				F2 += dE_dr_over_r * f2;
			}
		}
	} // loop over nbrs of rsd1
}


/// @brief RNA_VDW_Energy distance cutoff
Distance
RNA_VDW_Energy::atomic_interaction_cutoff() const
{
	return 5.0; /// now subtracted off 3.0 from cutoffs in centroid params files
	//return 0.0; /// since all the cutoffs for centroid mode are rolled into the cendist check
}

/// @brief RNA_VDW_Energy
void
RNA_VDW_Energy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{
}

//////////////////////////////////////////////////////////////////////////////////////////
//And make atom_numbers_for_vdw_calculation a vector of vectors, and
//   have a copy available in RNA_Scoring_info for setting (in this funciton) and getting in
//    the score & derivative functions above.

void
RNA_VDW_Energy::setup_atom_numbers_for_vdw_calculation( pose::Pose & pose ) const
{
	//We don't know a priori which atom numbers correspond to which
	// atom names (e.g., O2' on an adenosine could be different depending
	// on whether its at a chainbreak, terminus, etc.)
	//Better to do a quick setup every time to pinpoint atoms that require
	//  monitoring for VDW clashes.

	rna::RNA_ScoringInfo & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );

	if ( rna_scoring_info.vdw_calculation_annotated_sequence() == pose.annotated_sequence() ) return; // should be up to date
	rna_scoring_info.set_vdw_calculation_annotated_sequence( pose.annotated_sequence() );

	utility::vector1< utility::vector1< Size > > &
		atom_numbers_for_vdw_calculation( rna_scoring_info.nonconst_atom_numbers_for_vdw_calculation() );
	utility::vector1< bool > & is_magnesium = rna_scoring_info.nonconst_is_magnesium();

	Size const total_residue( pose.total_residue() );

	atom_numbers_for_vdw_calculation.resize( total_residue );
	is_magnesium.resize( total_residue );

	for ( Size i = 1; i <= total_residue; i++ ) {
		conformation::Residue const & rsd( pose.residue( i ) );

		is_magnesium[ i ] = ( rsd.name3() == " MG" );
		atom_numbers_for_vdw_calculation[ i ].clear();

		if ( !rsd.is_RNA() && !is_magnesium[ i ] ) continue;
		
		//a,c,g, or u?
		char const which_nucleotide = rsd.name1();
		//What atom names to look at?
		utility::vector1< std::string > const vdw_atom_list = rna_atom_vdw_.vdw_atom_list( which_nucleotide );
		for ( Size m = 1; m <= vdw_atom_list.size(); m++ ) {
			Size const vdw_atom_index =  rsd.atom_index( vdw_atom_list[ m ] );
			atom_numbers_for_vdw_calculation[ i ].push_back( vdw_atom_index );
			//std::cout << i << " " << rsd.name1() << " " << m << " of " << vdw_atom_list.size() << ": " << vdw_atom_list[ m ] << " " << vdw_atom_index << std::endl;
		}
	}
}

core::Size
RNA_VDW_Energy::version() const
{
	return 1; // Initial versioning
}


} //rna
} //scoring
} //core
