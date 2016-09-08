// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file cst_util.hh
/// @brief set of functions for adding constraints used in relax
/// @author James Thompson

#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/types.hh>

#include <protocols/relax/cst_util.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.relax.cst_util" );

namespace protocols {
namespace relax {

void coordinate_constrain_selection(
	core::pose::Pose & pose,
	core::sequence::SequenceAlignment aln,
	core::Real coord_sdev
) {

	if ( aln.size() == 0 ) {
		tr.Debug << "Empty alignment => no constraints!\n";
		tr.Debug.flush();
	}
	// add coordinate constraints over aligned residues
	using core::Size;
	using core::Real;
	using core::id::AtomID;
	using core::conformation::ResidueFactory;
	using namespace core::scoring::constraints;
	if ( pose.residue(pose.size()).name() != "VRT" ) {
		pose.append_residue_by_jump(
			*ResidueFactory::create_residue(
			pose.residue(1).residue_type_set()->name_map( "VRT" )
			),
			static_cast< Size > (pose.size() / 2)
		);
	}
	Size nres = pose.size();

	core::id::SequenceMapping map = aln.sequence_mapping(1,2);

	utility::vector1< core::Real > coord_sdevs;
	for ( Size idx = 1; idx <= nres - 1; ++idx ) {
		if ( map[idx] != 0 ) {
			coord_sdevs.push_back(coord_sdev);
		}
	}

	ConstraintSetOP cst_set = generate_bb_coordinate_constraints(
		pose, coord_sdevs
	);

	//for ( Size idx = 1; idx <= nres - 1; ++idx ) {
	// if ( map[idx] != 0 ) {
	//  using namespace core::conformation;
	//  Residue const & rsd( pose.residue(idx) );
	//  for ( Size ii = 1; ii <= rsd.last_backbone_atom(); ++ii ) {
	//   using namespace core::scoring::constraints;
	//   pose.add_constraint(
	//    new CoordinateConstraint(
	//     AtomID(ii,idx), AtomID(1,nres), rsd.xyz(ii),
	//     new HarmonicFunc(0.0,coord_sdev)
	//    )
	//   );
	//  } // ii
	// } // if residue is aligned
	//} // nres
} // coordinate_constraint_selection

core::scoring::constraints::ConstraintSetOP
generate_bb_coordinate_constraints(
	core::pose::Pose & pose,
	utility::vector1< core::Real > const & coord_sdevs
) {
	using core::Size;
	using core::Real;
	using core::id::AtomID;
	using namespace core::conformation;
	using namespace core::scoring::constraints;

	ConstraintSetOP cst_set( new ConstraintSet );
	for ( Size idx = 1; idx <= pose.size(); ++idx ) {
		if ( coord_sdevs.size() >= idx ) {
			Residue const & rsd( pose.residue(idx) );
			core::Real const coord_sdev( coord_sdevs[idx] );
			if ( coord_sdev > 0 ) {
				for ( Size ii = 1; ii <= rsd.last_backbone_atom(); ++ii ) {
					core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc(0.0,coord_sdev) );
					cst_set->add_constraint(
						ConstraintCOP( ConstraintOP( new CoordinateConstraint(
						AtomID(ii,idx), AtomID(1,pose.size()), rsd.xyz(ii),
						fx
						) ) )
					);
				}
			} // coord_sdev > 0
		} // have a coord_sdev for idx
	} // idx

	return cst_set;
} // generate_bb_coordinate_constraints

void delete_virtual_residues(
	core::pose::Pose & pose
) {
	using core::Size;
	// remove virtual residues
	for ( Size idx = 1; idx <= pose.size(); ++idx ) {
		if ( pose.residue_type(idx).name() == "VRT" ) {
			pose.conformation().delete_residue_slow(idx);
		}
	}
}

utility::vector1< core::Real >
get_per_residue_scores(
	core::pose::Pose & pose,
	core::scoring::ScoreType scoretype
) {
	using core::Size;
	using core::Real;
	using utility::vector1;
	using namespace core::scoring;
	vector1< Real > scores;
	for ( Size jj = 1; jj <= pose.size(); ++jj ) {
		EnergyMap rsd_energies( pose.energies().residue_total_energies(jj) );
		Real const per_residue_score( rsd_energies[ scoretype ] );
		scores.push_back( per_residue_score );
	}
	return scores;
}

void add_virtual_residue_to_cterm(
	core::pose::Pose & pose
) {
	using core::Size;
	using core::conformation::ResidueFactory;
	if ( pose.residue(pose.size()).name() != "VRT" ) {
		pose.append_residue_by_jump(
			*ResidueFactory::create_residue(
			pose.residue(1).residue_type_set()->name_map( "VRT" )
			),
			static_cast< Size > (pose.size() / 2)
		);
	}
}

void derive_sc_sc_restraints(
	core::pose::Pose & pose,
	core::Real const upper_dist_cutoff
) {
	using core::Size;
	using core::Real;
	using core::id::AtomID;
	using namespace core::scoring::constraints;

	Real const cst_sdev ( 2.0 );
	//Real const cst_width( 2.0 );

	for ( Size ii = 1, end_ii = pose.size(); ii <= end_ii; ++ii ) {
		for ( Size jj = ii+1, end_jj = pose.size(); jj <= end_jj; ++jj ) {
			const core::conformation::Residue& resi = pose.residue(ii);
			const core::conformation::Residue& resj = pose.residue(jj);

			for ( Size atm_ii = resi.first_sidechain_atom(); atm_ii <= resi.natoms(); ++atm_ii ) {
				for ( Size atm_jj = resj.first_sidechain_atom(); atm_jj <= resj.natoms(); ++atm_jj ) {
					// determine whether Euclidean distance is below threshold
					Real const distance(
						resi.xyz(atm_ii).distance( resj.xyz(atm_jj) ));
					if ( distance <= upper_dist_cutoff ) {
						core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( distance, cst_sdev ) );
						pose.add_constraint(
							core::scoring::constraints::ConstraintCOP( core::scoring::constraints::ConstraintOP( new AtomPairConstraint(
							AtomID(atm_ii,ii), AtomID(atm_jj,jj),
							fx
							) ) ) );
						tr << "adding restraint from "
							<< "AtomID(" << atm_ii << "," << ii
							<< "AtomID(" << atm_jj << "," << jj
							<< ") with distance = " << distance << std::endl;
					}
				} // atm_jj
			} // atm_ii
		} // res jj
	} // res ii
	tr.flush_all_channels();
} // derive_sc_sc_restraints

} // relax
} // protocols
