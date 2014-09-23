// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/relax/AtomCoordinateCstMover.cc
/// @brief Add coordinate constraints to a pose, as for the relax protocol
/// @author Rocco Moretti (rmoretti@u.washington.edu)

#include <protocols/relax/AtomCoordinateCstMover.hh>
#include <protocols/relax/AtomCoordinateCstMoverCreator.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/SequenceMapping.hh>
#include <core/pack/task/TaskFactory.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/datacache/DataMap.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/sequence/util.hh>
#include <core/pose/util.hh>

#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/BoundConstraint.hh>

#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

static thread_local basic::Tracer TR( "protocols.relax.AtomCoordinateCstMover" );

namespace protocols {
namespace relax {

std::string
AtomCoordinateCstMoverCreator::keyname() const
{
	return AtomCoordinateCstMoverCreator::mover_name();
}

protocols::moves::MoverOP
AtomCoordinateCstMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AtomCoordinateCstMover() );
}

std::string
AtomCoordinateCstMoverCreator::mover_name()
{
	return "AtomCoordinateCstMover";
}


////////////////////////////////////////////////////////////////////////////////////////////////////

AtomCoordinateCstMover::AtomCoordinateCstMover() :
		cst_sd_( 0.5 ),
		bounded_( false ),
		cst_width_( 0 ),
		cst_sidechain_( false ),
		amb_hnq_( false )
{}

AtomCoordinateCstMover::AtomCoordinateCstMover( AtomCoordinateCstMover const & other ) : protocols::moves::Mover(other),
		refpose_( other.refpose_ ),
		cst_sd_( other.cst_sd_ ),
		bounded_( other.bounded_ ),
		cst_width_( other.cst_width_ ),
		cst_sidechain_( other.cst_sidechain_ ),
		amb_hnq_( other.amb_hnq_ ),
		loop_segments_( other.loop_segments_ ),
		task_segments_( other.task_segments_ )
{}

AtomCoordinateCstMover::~AtomCoordinateCstMover() {}

protocols::moves::MoverOP AtomCoordinateCstMover::fresh_instance() const {
	return protocols::moves::MoverOP( new AtomCoordinateCstMover );
}
protocols::moves::MoverOP AtomCoordinateCstMover::clone() const {
	return protocols::moves::MoverOP( new AtomCoordinateCstMover( *this ) );
}

void
AtomCoordinateCstMover::set_loop_segments( protocols::loops::LoopsCOP loops) {
	loop_segments_ = loops;
}

void
AtomCoordinateCstMover::set_task_segments( core::pack::task::TaskFactoryCOP taskfactory) {
	task_segments_ = taskfactory;
}

void
AtomCoordinateCstMover::set_refstruct( core::pose::PoseCOP ref ) {
	refpose_ = ref;
}

void AtomCoordinateCstMover::apply( core::pose::Pose & pose) {

	using namespace core::scoring::constraints;

	core::pose::Pose constraint_target_pose = pose;
	core::id::SequenceMapping seq_map; // A mapping of pose -> constraint_target_pose numbering

	if( refpose_ ) {
		constraint_target_pose = *refpose_;

		if (  pose.total_residue() == constraint_target_pose.total_residue() &&
				( ! cst_sidechain_ || pose.sequence() != constraint_target_pose.sequence() ) ) {
			// We match in size and (for sidechains) sequence - we're looking at the traditional 1:1 mapping.
			seq_map = core::id::SequenceMapping::identity( pose.total_residue() );
		} else {
			// Try to match on a PDB-identity basis, or a sequence alignment basis if that fails.
			TR << "Length " << (cst_sidechain_?"and/or identities ":"") <<
					"of input structure and native don't match - aligning on PDB identity or sequence." << std::endl;
			seq_map = core::pose::sequence_map_from_pdbinfo( pose, constraint_target_pose );
		}
		// Align the native pose to the input pose to avoid rotation/translation based
		//  errors.
		//fpd  (Only if not already rooted on a VRT to avoid problems with density/symmetry)
		if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
			core::id::SequenceMapping rev_seq_map( seq_map ); // constraint_target_pose -> pose mapping
			rev_seq_map.reverse();
			core::sequence::calpha_superimpose_with_mapping(constraint_target_pose, pose, rev_seq_map);
		}
	} else {
		// Aligning to input - mapping is 1:1
		seq_map = core::id::SequenceMapping::identity( pose.total_residue() );
	}

	core::Size nres = pose.total_residue();

	utility::vector1< bool > constrain_residues( nres, false );
	core::pack::task::PackerTaskOP task;
	if( task_segments_ ) {
		task = task_segments_->create_task_and_apply_taskoperations(pose);
	}
	for( core::Size ii = 1; ii <= nres; ++ii ) {
		// Constrain residue if it's in a loop, or in the task, or if we don't have task/loop limitations.
		if ( (loop_segments_ && loop_segments_->is_loop_residue( ii )) ||
				( task && task->being_packed(ii) ) ||
				( !loop_segments_ && !task ) ) {
			constrain_residues[ ii ] = true;
		}
	}

	// Warn about not having a virtual root (but go ahead with constraints).
	if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
		TR.Warning << "WARNING: Adding coordinate constraints to a pose without a virtual root - results may not be as expected." << std::endl;
	}

	for ( core::Size i = 1; i<= nres; ++i ) {
		if ( pose.residue( i ).aa() == core::chemical::aa_vrt ) continue; // Skip virtual residues.

		if ( constrain_residues[i] ) {
			core::Size j(seq_map[i]);
			if( j == 0 ) continue;
			assert( j <= constraint_target_pose.total_residue() ); // Should be, if map was set up properly.

			core::conformation::Residue const & pose_i_rsd( pose.residue(i) );
			core::conformation::Residue const & targ_j_rsd( constraint_target_pose.residue(j) );
			core::Size last_atom( pose_i_rsd.last_backbone_atom() );
			core::Size last_targ_atom( targ_j_rsd.last_backbone_atom() );
			bool use_atom_names(false);
			if ( cst_sidechain_ ) {
				last_atom = pose_i_rsd.nheavyatoms();
				last_targ_atom = targ_j_rsd.nheavyatoms();
				use_atom_names = pose_i_rsd.name() != targ_j_rsd.name(); // Don't bother with lookup if they're the same residue type.
			}
			if ( !use_atom_names && last_atom != last_targ_atom ) {
				TR.Warning << "Warning: Coordinate constraint reference residue has different number of " << (cst_sidechain_?"heavy":"backbone") << " atoms: ref. "
					<< targ_j_rsd.name() << " (res " << j << ") versus  " << pose_i_rsd.name() << " (res " << i << "). - skipping." << std::endl;
				continue;
			}
			for ( core::Size ii = 1; ii<= last_atom; ++ii ) {
				core::Size jj(ii);
				if ( use_atom_names ) {
					std::string atomname( pose_i_rsd.atom_name(ii) );
					if ( ! targ_j_rsd.has(atomname) ) {
						TR.Debug << "Skip adding coordinate constraints for atom " << atomname << " of residue " << i << " (" << pose_i_rsd.name() <<
							") - not found in residue " << j << " (" << targ_j_rsd.name() << ") of reference structure." << std::endl;
						continue;
					}
					jj = targ_j_rsd.atom_index( atomname );
				}
				core::scoring::func::FuncOP function;
				if( bounded_ ) {
					function = core::scoring::func::FuncOP( new BoundFunc( 0, cst_width_, cst_sd_, "xyz" ) );
				} else {
					function = core::scoring::func::FuncOP( new core::scoring::func::HarmonicFunc( 0.0, cst_sd_ ) );
				}

				// Rely on shortcutting evaluation to speed things up - get to else clause as soon as possible.
				if( amb_hnq_ && cst_sidechain_ &&
						( (pose_i_rsd.aa() == core::chemical::aa_asn && targ_j_rsd.aa() == core::chemical::aa_asn &&
								( pose_i_rsd.atom_name(ii) == " OD1" || pose_i_rsd.atom_name(ii) == " ND2" )) ||
						(pose_i_rsd.aa() == core::chemical::aa_gln && targ_j_rsd.aa() == core::chemical::aa_gln &&
								( pose_i_rsd.atom_name(ii) == " OE1" || pose_i_rsd.atom_name(ii) == " NE2" )) ||
						(pose_i_rsd.aa() == core::chemical::aa_his && targ_j_rsd.aa() == core::chemical::aa_his &&
								(pose_i_rsd.atom_name(ii) == " ND1" || pose_i_rsd.atom_name(ii) == " NE2" ||
										pose_i_rsd.atom_name(ii) == " CD2" || pose_i_rsd.atom_name(ii) == " CE1")) ) ) {
					std::string atom1, atom2;
					if ( pose_i_rsd.aa() == core::chemical::aa_asn ) {
						atom1 = " OD1";
						atom2 = " ND2";
					} else if ( pose_i_rsd.aa() == core::chemical::aa_gln ) {
						atom1 = " OE1";
						atom2 = " NE2";
					} else if ( pose_i_rsd.aa() == core::chemical::aa_his &&
							(pose_i_rsd.atom_name(ii) == " ND1" || pose_i_rsd.atom_name(ii) == " CD2") ) {
						atom1 = " ND1";
						atom2 = " CD2";
					} else if ( pose_i_rsd.aa() == core::chemical::aa_his &&
							(pose_i_rsd.atom_name(ii) == " NE2" || pose_i_rsd.atom_name(ii) == " CE1") ) {
						atom1 = " NE2";
						atom2 = " CE1";
					} else {
						utility_exit_with_message("Logic error in AtomCoordinateConstraints");
					}
					core::scoring::constraints::AmbiguousConstraintOP amb_constr( new core::scoring::constraints::AmbiguousConstraint );
					amb_constr->add_individual_constraint( ConstraintCOP( new CoordinateConstraint(  core::id::AtomID(ii,i),
							core::id::AtomID(1,pose.fold_tree().root()), targ_j_rsd.xyz( atom1 ), function ) ) );
					amb_constr->add_individual_constraint( ConstraintCOP( new CoordinateConstraint(  core::id::AtomID(ii,i),
							core::id::AtomID(1,pose.fold_tree().root()), targ_j_rsd.xyz( atom2 ), function ) ) );
					pose.add_constraint( amb_constr );
				} else {
					pose.add_constraint( core::scoring::constraints::ConstraintCOP( new CoordinateConstraint( core::id::AtomID(ii,i),
							core::id::AtomID(1,pose.fold_tree().root()), targ_j_rsd.xyz( jj ), function ) ) );
				}
			} // for atom
		} // if(loop)
	} // for residue
}

void
AtomCoordinateCstMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	cst_sd( tag->getOption< core::Real >( "coord_dev", 0.5 ) );
	bounded( tag->getOption< bool >( "bounded", false ) );
	cst_width( tag->getOption< core::Real >( "bound_width", 0 ) );
	cst_sidechain( tag->getOption< bool >( "sidechain", false ) );
	ambiguous_hnq( tag->getOption< bool >( "flip_hnq", false ) );

	if ( tag->hasOption("task_operations") ) { // Only do task operations if set
		set_task_segments( protocols::rosetta_scripts::parse_task_operations(tag, data) );
	}

	if ( tag->getOption< bool >( "native", false ) ) {
		refpose_ = get_native_pose();
		if ( ! refpose_ ) {
			throw utility::excn::EXCN_RosettaScriptsOption("Use native for AtomCoordinateCstMover specified, but not native pose is availible.");
		}
	}

}

} // namespace relax
} // namespace protocols
