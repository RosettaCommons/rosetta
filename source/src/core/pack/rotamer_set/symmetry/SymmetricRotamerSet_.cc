// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.cc
/// @brief  amino acid rotamer set class implementation for symmetric packing
/// @author Ingemar Andre

// Unit Headers
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>


// Project Headers
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
#include <basic/Tracer.hh>

#include <core/graph/Graph.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/Methods.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2.hh>

// C++ headers
#include <string>
#include <iostream>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <ObjexxFCL/FArray2D.hh>

namespace core {
namespace pack {
namespace rotamer_set {
namespace symmetry {

static thread_local basic::Tracer tt( "core.pack.rotamer_set.symmetry.SymmetricRotamerSet_", basic::t_info );

SymmetricRotamerSet_::SymmetricRotamerSet_()
: RotamerSet_()
{}

SymmetricRotamerSet_::~SymmetricRotamerSet_() {}

// @ details The packer's 1-body is a combination of the rotamer internal energies (the
// context dependent and independent one body energies), the intra-residue
// energies defined by the two body energies, and the sum of the
// two body energies with the background residues in this repacking.  The
// PackerTask tells the RotamerSet_ what residues are part of the
// background and which are being repacked.
// This is the symmetrical version of the class. It adds the score for the clone residues into
// the energy for the master residue
void
SymmetricRotamerSet_::compute_one_body_energies(
	pose::Pose const & pose,
	scoring::ScoreFunction const & sf,
	task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph,
	utility::vector1< core::PackerEnergy > & energies
) const
{
	using namespace conformation;
	using namespace scoring;

	//fpd  make sure the pose is symmetric
	if ( !core::pose::symmetry::is_symmetric( pose ) ) {
		RotamerSet_::compute_one_body_energies( pose, sf, task, packer_neighbor_graph, energies );
		return;
	}

	// find SymmInfo
  SymmetricConformation const & SymmConf (
    dynamic_cast<SymmetricConformation const &> ( pose.conformation() ) );
  SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	std::fill( energies.begin(), energies.end(), core::PackerEnergy( 0.0 ) );
	utility::vector1< core::PackerEnergy > temp_energies = energies;

	int const nrotamers = num_rotamers(); // does not change in this function
	Size const theresid = resid();

	for ( int ii = 1; ii <= nrotamers; ++ii ) {
		EnergyMap emap;
		sf.eval_ci_1b( *rotamer( ii ), pose, emap );
		sf.eval_cd_1b( *rotamer( ii ), pose, emap );
		energies[ ii ] += static_cast< core::PackerEnergy > (sf.weights().dot( emap ))*symm_info->score_multiply_factor(); // precision loss here. Multiply with the number of total subunits to get the total energy
	}


	sf.evaluate_rotamer_intrares_energies( *this, pose, temp_energies );
	// multiply the rotamer_intrares energy with the number of subunits in the system
	PackerEnergyMultiply( temp_energies, symm_info->score_multiply_factor() );
    PackerEnergyAdd( energies, temp_energies );

	// define a factory
	RotamerSetFactory rsf;

	for ( graph::Graph::EdgeListConstIter
			ir  = packer_neighbor_graph->get_node( theresid )->const_edge_list_begin(),
			ire = packer_neighbor_graph->get_node( theresid )->const_edge_list_end();
			ir != ire; ++ir ) {
		int const neighbor_id( (*ir)->get_other_ind( theresid ) );
		if ( task.pack_residue( neighbor_id ) ) continue;
		Residue const & neighbor( pose.residue( neighbor_id ) );
		uint const symm_clone ( symm_info->chi_follows( neighbor_id ) );
		if ( symm_clone != 0 && task.pack_residue( symm_clone )
					&& symm_clone != theresid )  continue;
		// Residue is not interating with itself
		if ( symm_clone != theresid ) {
			// evaluate the energy for one rotamer pair interaction and
			// multiply it with the number of subunits in the system
			std::fill( temp_energies.begin(), temp_energies.end(), core::PackerEnergy( 0.0 ) );
			sf.evaluate_rotamer_background_energies( *this, neighbor, pose, temp_energies );
			PackerEnergyMultiply( temp_energies, symm_info->score_multiply_factor() );
			PackerEnergyAdd( energies, temp_energies );
//					PackerEnergyAdd( energies, temp_energies );
		} else {
			// We have a self interaction. We have to calculate the rotamer-rotamer pair interaction
			// where both rotamers change at the same time to the same state. We go through all rotamers
			// and calculated interaction energies with translated copies of itself
			for ( int jj = 1; jj <= nrotamers; ++jj ) {
				// make a new rotamer set that is going to be translated to the neighbor interation residue
				conformation::ResidueOP sym_rsd( this->rotamer( jj )->clone() );
				RotamerSetOP one_rotamer_set = rsf.create_rotamer_set( *sym_rsd );
			  one_rotamer_set->set_resid( theresid );
			  one_rotamer_set->add_rotamer( *sym_rsd );
				// place rotamer set at neighbor position
				RotamerSetOP sym_rotset(
            orient_rotamer_set_to_symmetric_partner( pose, sym_rsd, neighbor_id, one_rotamer_set ) );
				sf.prepare_rotamers_for_packing( pose, *one_rotamer_set );
				sf.prepare_rotamers_for_packing( pose, *sym_rotset );
				// make a temporary core::PackerEnergy object with on rotamer in it
		    ObjexxFCL::FArray2D< core::PackerEnergy > temp_energies( 1, 1, 0.0 );
				// evaluate the energy for this rotamer-rotamer interaction
		    sf.evaluate_rotamer_pair_energies(
              *one_rotamer_set, *sym_rotset, pose, temp_energies );
				// add the energy of this interaction. Mulitply with the number of subunits in
				// the system
				energies[ jj ] += temp_energies[ 0 ]*symm_info->score_multiply( theresid, neighbor_id );
			}
		}
	}

	// long-range energy interactions with background
	// Iterate across the long range energy functions and use the iterators generated
	// by the LRnergy container object
	for ( ScoreFunction::LR_2B_MethodIterator
			lr_iter = sf.long_range_energies_begin(),
			lr_end  = sf.long_range_energies_end();
			lr_iter != lr_end; ++lr_iter ) {
		LREnergyContainerCOP lrec = pose.energies().long_range_container( (*lr_iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue; // only score non-emtpy energies.

		// Potentially O(N) operation leading to O(N^2) behavior
		for ( ResidueNeighborConstIteratorOP
						rni = lrec->const_neighbor_iterator_begin( theresid ),
						rniend = lrec->const_neighbor_iterator_end( theresid );
					(*rni) != (*rniend); ++(*rni) ) {
			Size const neighbor_id = rni->neighbor_id();
		debug_assert( neighbor_id != theresid );
			if ( task.pack_residue( neighbor_id ) ) continue;

			// Residue is not interating with itself
			if ( symm_info->chi_follows( theresid ) != neighbor_id ) {
			// evaluate the energy for one rotamer pair interaction and
			// multiply it with the number of subunits in the system
				std::fill( temp_energies.begin(), temp_energies.end(), core::PackerEnergy( 0.0 ) );
				(*lr_iter)->evaluate_rotamer_background_energies(
					*this, pose.residue( neighbor_id ), pose, sf,
					sf.weights(), temp_energies );
					PackerEnergyMultiply( temp_energies, symm_info->score_multiply_factor() );
					PackerEnergyAdd( energies, temp_energies );
			} else {
				// We have a self interaction. We have to calculate the rotamer-rotamer pair interaction
				// where both rotamers change at the same time to the same state. We go through all rotamers
				// and calculated interaction energies with translated copies of itself
				//Residue const & neighbor( pose.residue( neighbor_id ) );
				for ( int jj = 1; jj <= nrotamers; ++jj ) {
					// make a new rotamer set that is going to be translated to the neighbor interation residue
					conformation::ResidueOP sym_rsd( this->rotamer( jj )->clone() );
					RotamerSetOP one_rotamer_set = rsf.create_rotamer_set( *sym_rsd );
					one_rotamer_set->set_resid( theresid );
					one_rotamer_set->add_rotamer( *sym_rsd );
					// place rotamer set at neighbor position
					RotamerSetOP sym_rotset(
						orient_rotamer_set_to_symmetric_partner( pose, sym_rsd, neighbor_id, one_rotamer_set ) );
					sf.prepare_rotamers_for_packing( pose, *one_rotamer_set );
					sf.prepare_rotamers_for_packing( pose, *sym_rotset );
					ObjexxFCL::FArray2D< core::PackerEnergy > temp_energies( 1, 1, 0.0 );
					// evaluate the energy for this rotamer-rotamer interaction
					(*lr_iter)->evaluate_rotamer_pair_energies(
               *one_rotamer_set, *sym_rotset, pose, sf, sf.weights(), temp_energies );
					// add the energy of this interaction. Mulitply with the number of subunits in
					// the system
					 energies[ jj ] += temp_energies[ 0 ]*symm_info->score_multiply( theresid, neighbor_id );
				}
			}
		} // (potentially) long-range neighbors of theresid [our resid()]
	} // long-range energy functions
}
// @details function for mulitplying core::PackerEnergy objects. Should make a * operator in core::PackerEnergy
// instead
void
SymmetricRotamerSet_::PackerEnergyMultiply( utility::vector1< core::PackerEnergy > & energies,
											Size factor ) const
{
		int const nrotamers = num_rotamers();
		for ( int ii = 1; ii <= nrotamers; ++ii ){
			energies[ii] = energies[ii]*factor;
	}
}

// @details function for adding core::PackerEnergy objects. Should make a + operator in core::PackerEnergy
// instead. The objects must have the same size
void
SymmetricRotamerSet_::PackerEnergyAdd( utility::vector1< core::PackerEnergy > & energies,
								 utility::vector1< core::PackerEnergy > const & add ) const
{
debug_assert( energies.size() == add.size() );
	int const nrotamers = num_rotamers();
	for ( int ii = 1; ii <= nrotamers; ++ii ){
		energies[ii] = energies[ii] + add[ii];
	}
}

// @details function for subtract core::PackerEnergy objects. Should make a - operator in core::PackerEnergy
// instead. The objects must have the same size
void
SymmetricRotamerSet_::PackerEnergySubtract( utility::vector1< core::PackerEnergy > & energies,
                 utility::vector1< core::PackerEnergy > const & subtract ) const
{
 debug_assert( energies.size() == subtract.size() );
	int const nrotamers = num_rotamers();
  for ( int ii = 1; ii <= nrotamers; ++ii ){
    energies[ii] = energies[ii] - subtract[ii];
  }
}

// @details translate a rotamer_set to a different sequence position
RotamerSetOP
SymmetricRotamerSet_::orient_rotamer_set_to_symmetric_partner(
  pose::Pose const & pose,
  conformation::ResidueOP residue_in,
  int const & sympos,
	RotamerSetOP rotset_in
) const
{

  RotamerSetFactory rsf;
  RotamerSetOP sym_rotamer_set = rsf.create_rotamer_set( *residue_in );
  for (Rotamers::const_iterator
    rot     = rotset_in->begin(),
    rot_end = rotset_in->end();
    rot != rot_end; ++rot) {
      conformation::Residue target_rsd( *residue_in );

			// peptoids have a different orientation function due to the placement of the first side chain atom
			if ( target_rsd.type().is_peptoid() ) {
				target_rsd.orient_onto_residue_peptoid( pose.residue( sympos ), pose.conformation() );
			} else {
				target_rsd.orient_onto_residue( pose.residue( sympos ) );
			}

			target_rsd.copy_residue_connections_from( pose.residue( sympos ) );
			sym_rotamer_set->set_resid( sympos );
      sym_rotamer_set->add_rotamer( target_rsd );
  }
	return sym_rotamer_set;
}

}	// symmetry
} // rotamer_set
} // pack
} // core
