// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/drug_design/util.cc
/// @brief Utilities for DrugDesign
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/drug_design/util.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomRefMapping.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResidueLevelTask.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>

#include <utility/graph/Graph.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

#include <numeric/model_quality/rms.hh>
#include <numeric/xyz.functions.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.drug_design.util" );


namespace protocols {
namespace drug_design {

void align_residues(
	core::conformation::Residue & residue,
	core::conformation::Residue const & target,
	core::chemical::IndexIndexMapping const & map // from target to residue
){
	utility::vector1<core::PointPosition> coords, target_coords;
	utility::vector1<core::Size> vectmap;
	for ( core::Size ii(1); ii <= target.natoms(); ++ii ) {
		if ( map[ii] != map.invalid_entry() ) {
			target_coords.push_back( target.xyz( ii ) );
			coords.push_back( residue.xyz( map[ii] ) );
			vectmap.push_back( map[ii] );
		} else {
			vectmap.push_back(0);
		}
	}

	if ( coords.size() == 0 ) {
		return; // No atoms to overlay, therefore no alignment can be done.
	}

	core::Vector orig_center( numeric::center_of_mass(coords) ), new_center( numeric::center_of_mass(target_coords) );
	numeric::xyzMatrix<core::Real> rot_matrix;

	if ( coords.size() >= 3 ) {
		utility::vector1<core::Real> weights(target_coords.size(),1.0);
		core::Real sigma3 = 0.0; //unused but findUU() needs it

		// Find the rotation matrix which minimizes rmsd
		// (after translating both sets to around the origin)
		numeric::model_quality::findUU(target_coords, coords, weights, target_coords.size(), rot_matrix, sigma3);
	} else if ( coords.size() == 2 ) {
		// Because findUU is being an issue with two points
		rot_matrix = numeric::alignVectorSets(coords[1]-orig_center, coords[2]-orig_center, target_coords[1]-new_center, target_coords[2]-new_center);
	} else if ( coords.size() == 1 ) {
		rot_matrix.to_identity(); // No rotation needed for single point alignment
	} else {
		return; // No atoms to overlay, therefore no alignment can be done. (Redundant with above to make clang_tidy happy.)
	}

	// As the rotation found is for the center of mass being around the origin,
	// We need to translate the starting residue to the origin before applying the rotation
	for ( core::Size ii(1); ii <= residue.natoms(); ++ii ) {
		residue.set_xyz( ii, rot_matrix*( residue.xyz(ii) - orig_center ) + new_center );
	}

	if ( TR.Debug.visible() ) {
		TR << "align_residues() Super/No Super Delta: " << core::scoring::residue_rmsd_nosuper(target, residue, vectmap) - core::scoring::residue_rmsd_super(target, residue, vectmap) << std::endl;
	}
}

void
place_new_restype(
	core::pose::Pose & pose,
	core::Size position,
	core::chemical::ResidueType const & new_restype,
	core::chemical::IndexIndexMapping const & map, // from original residuetype in pose to new_restype
	std::string const & method_mc
) {
	std::string method = utility::lowercased( method_mc );
	if ( method.empty() || method == "default" ) {
		// Current "best" approach is the rotamer align method (mostly from lack of choice).
		place_new_restype_rotamer_align( pose, position, new_restype, map );
	} else if ( method == "rotamer_align" ) {
		place_new_restype_rotamer_align( pose, position, new_restype, map );
	} else if ( method == "none" || method == "no_align" ) {
		place_new_restype_no_align( pose, position, new_restype );
	}
}

void
place_new_restype_no_align(
	core::pose::Pose & pose,
	core::Size position,
	core::chemical::ResidueType const & new_restype
) {
	core::conformation::Residue new_res( new_restype, true );

	pose.replace_residue(position, new_res, false /* don't align */);
}

void
place_new_restype_rotamer_align(
	core::pose::Pose & pose,
	core::Size position,
	core::chemical::ResidueType const & new_restype,
	core::chemical::IndexIndexMapping const & map // from original residuetype in pose to new_restype
) {

	if ( map.empty() ) {
		TR.Warning << "When placing residue, the atom-atom mapping is empty." << std::endl;
	}
	std::clock_t begin = std::clock();

	// I used to use RTMin for this, but stopped because it was taking too much time.
	// Instead, we pick the rotamer that has the best superimposable rmsd, and then
	// minimize that in context.

	// Setup a minimal scorefunction which just has the overlay terms
	// and minimal clash avoidance terms.
	core::scoring::ScoreFunction scorefxn;
	// The weights here are a total guess
	scorefxn.set_weight(core::scoring::fa_rep, 0.5);
	scorefxn.set_weight(core::scoring::fa_intra_rep_xover4, 0.5);
	scorefxn.set_weight(core::scoring::coordinate_constraint, 10.0);

	scorefxn(pose); // Need an initial scoring such that the neighbor graph exists

	core::conformation::Residue const & current_residue( pose.residue( position ) );
	utility::vector1< core::conformation::ResidueOP > possible_rotamers;

	core::pack::rotamers::SingleResidueRotamerLibraryCOP rotlib( core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( new_restype ) );

	if ( rotlib ) {
		core::pack::task::PackerTaskOP help_task = core::pack::task::TaskFactory::create_packer_task( pose );
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			if ( i == position ) help_task->nonconst_residue_task( i ).restrict_to_repacking();
			else help_task->nonconst_residue_task( i ).prevent_repacking();
		}
		utility::graph::GraphOP neighbor_graph( new utility::graph::Graph( pose.energies().energy_graph() ) );
		//core::graph::GraphOP neighbor_graph( new core::graph::Graph( pose.energies().energy_graph() ) );
		utility::vector1< utility::vector1< core::Real > > extra_chi_steps( new_restype.nchi() );

		utility::vector1< core::conformation::ResidueOP > all_rotamers;
		rotlib->fill_rotamer_vector( pose, scorefxn, *help_task, neighbor_graph, new_restype.get_self_ptr(), current_residue, extra_chi_steps, true, all_rotamers );

		TR << "Number of Total Rotamers generated: " << all_rotamers.size() << std::endl;

		utility::vector1<core::Size> atom_map(current_residue.natoms(), 0);
		for ( core::Size aa(1); aa <= current_residue.natoms(); ++aa ) {
			atom_map[aa] = map[aa];
		}

		core::Real min_rmsd = 999999999;
		//core::Size min_rotamer = 0;
		utility::vector1< core::Real > rmsds;
		for ( core::Size rr(1); rr <= all_rotamers.size(); ++rr ) {
			core::Real current_rmsd( core::scoring::residue_rmsd_super( current_residue, *all_rotamers[rr], atom_map, /*skip_hydro=*/ false ) );
			if ( utility::isnan(current_rmsd) ) {
				TR.Warning << "NaN encountered during alignment!" << std::endl;
			}

			rmsds.push_back( current_rmsd );
			if ( ! utility::isnan(current_rmsd) && ( min_rmsd == 999999999 || current_rmsd < min_rmsd ) ) {
				//min_rotamer = rr;
				min_rmsd = current_rmsd;
			}

		}
		TR << "When placing ligand, RMSD for best rotamer is " << min_rmsd << std::endl;
		for ( core::Size rr(1); rr <= rmsds.size(); ++rr ) {
			if ( ! utility::isnan(rmsds[rr]) && (rmsds[rr] < min_rmsd + 0.5) ) {
				align_residues( *all_rotamers[rr], current_residue, map );
				possible_rotamers.push_back( all_rotamers[rr] );
			}
		}
	}

	if ( possible_rotamers.empty() ) {
		if ( rotlib ) {
			TR.Warning << "Alignment for every rotamer resulted in NaN." << std::endl;
		}
		core::conformation::ResidueOP only_rot( new core::conformation::Residue(new_restype, true) );
		align_residues( *only_rot, current_residue, map );
		possible_rotamers.push_back( only_rot );
	}

	TR << "Number of Reasonable Rotamers generated: " << possible_rotamers.size() << std::endl;

	runtime_assert( ! possible_rotamers.empty() );

	// We need to remove any constraints on the off chance they're not fully compatible with the ligand
	// (e.g. if an atom doesn't exist - that can be checked with the AtomExistsFilter,
	// but only after we have placed the ligand in the pose.)
	core::scoring::constraints::ConstraintSetOP stored_constraints;
	if ( pose.constraint_set() ) {
		stored_constraints = pose.constraint_set()->clone();
	}
	pose.remove_constraints();

	// Setup constraints for the matching atoms
	core::scoring::constraints::ConstraintCOPs atom_constraints;
	{ // Scope to limit residue, as we'll be replacing it later in the function
		core::Size ft_root( pose.fold_tree().root() );
		core::conformation::Residue const & residue( pose.residue( position ) );
		for ( core::Size ii(1); ii <= residue.natoms(); ++ii ) {
			if ( map[ii] != map.invalid_entry() ) {
				using namespace core::scoring;
				func::FuncOP function( new func::HarmonicFunc( 0.0, 0.5 ) );
				constraints::ConstraintCOP constraint( new constraints::CoordinateConstraint(
					core::id::AtomID(map[ii], position),
					core::id::AtomID(1, ft_root),
					residue.xyz(ii),
					function ));
				atom_constraints.push_back( constraint );
			}
		}
	}
	pose.add_constraints( atom_constraints );

	// Figure out which reasonable rotamer has the lowest energy (constraint and repulsives)
	core::Real min_score(999999999);
	core::Size min_rotamer(0);
	for ( core::Size rr(1); rr <= possible_rotamers.size(); ++rr ) {
		pose.replace_residue(position, *possible_rotamers[rr], false /* don't re-align */);
		core::Real current_score( scorefxn(pose) );
		if ( min_rotamer == 0  || current_score < min_score ) {
			min_rotamer = rr;
			min_score = current_score;
		}
	}

	pose.replace_residue(position, *possible_rotamers[min_rotamer], false /* don't realign */);

	// We probably don't need to do any sort of minimization here
	// If someone wants to minimize, they should do it themselves later in the protocol

	// Replace the constraint set with the original one.
	pose.constraint_set( stored_constraints );

	std::clock_t end = std::clock();
	TR << "Placement took: " << double(end-begin)/ CLOCKS_PER_SEC << "s total" << std::endl;
}




} //protocols
} //drug_design


