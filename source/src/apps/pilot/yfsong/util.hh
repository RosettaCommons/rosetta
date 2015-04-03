// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief pilot app utilities
/// @details
/// @author Yifan Song

#ifndef apps_pilot_yfsong_util_HH
#define apps_pilot_yfsong_util_HH

#include <apps/pilot/yfsong/AlignChunkMover.hh>

#include <ObjexxFCL/FArray2D.hh>

#include <protocols/simple_moves/MutateResidue.hh>

#include <core/kinematics/Edge.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/VariantType.hh>

#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <utility/pointer/access_ptr.hh>

#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>

namespace challenge {
	basic::options::IntegerOptionKey	max_registry_shift("challenge:max_registry_shift");
	basic::options::BooleanOptionKey	virtual_loops("challenge:virtual_loops");
	basic::options::BooleanOptionKey	revert_real_loops("challenge:revert_real_loops");
}


namespace apps {
namespace pilot {

using namespace core;
using namespace core::chemical;
using namespace core::kinematics;
using namespace ObjexxFCL;
using namespace protocols::moves;
using namespace protocols::loops;
using namespace basic::options;
using namespace basic::options::OptionKeys;

class SampleSecondaryStructureAlignmentMover: public protocols::moves::Mover
{
public:

SampleSecondaryStructureAlignmentMover(numeric::random::RandomGenerator & RG,
									   utility::vector1 < utility::file::FileName > const & template_filenames)
	: RG_(RG) {
		//initialize template structures
		read_template_structures(template_poses_, template_filenames);
	}

	void read_template_structures(utility::vector1 <core::pose::PoseOP> & template_structures, utility::vector1 < utility::file::FileName > const & template_filenames)
	{
		template_structures.clear();

		template_structures.resize(template_filenames.size());
		for (core::Size i_ref=1; i_ref<= template_filenames.size(); ++i_ref) {
			template_structures[i_ref] = new core::pose::Pose();
			core::import_pose::pose_from_pdb( *(template_structures[i_ref]), template_filenames[i_ref] );

			protocols::jumping::Dssp dssp_obj( *template_structures[i_ref] );
			dssp_obj.insert_ss_into_pose( *template_structures[i_ref] );
		}
	}

void
set_loops_to_virt_ala(core::pose::Pose & pose, Loops loops)
{
	chemical::ResidueTypeSet const& restype_set( pose.residue(1).residue_type_set() );

	for (Size iloop=1; iloop<=loops.num_loop(); ++iloop) {
		/*
		if ( loops[iloop].start() == 1 ) {
			if (loops[iloop].stop() <=2) continue;
		}
		else if (loops[iloop].stop() - loops[iloop].start() <= 3) {
			continue;
		}
		else {
			Size nres(pose.total_residue());
			while (pose.residue_type(nres).is_virtual_residue()) {
				--nres;
			}
			if ( loops[iloop].stop() == nres ) {
				if (loops[iloop].start() >=nres-1) continue;
			}
		}
		 */

		for (Size ires=loops[iloop].start(); ires<=loops[iloop].stop(); ++ires) {

			// Create the new residue and replace it
			conformation::ResidueOP new_res = conformation::ResidueFactory::create_residue(
																						   restype_set.name_map("VBB"), pose.residue(ires),
																						   pose.conformation());
			// Make sure we retain as much info from the previous res as possible
			conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue(ires),
																			 *new_res, pose.conformation() );
			pose.replace_residue(ires, *new_res, false );
			//core::pose::add_variant_type_to_pose_residue(pose, "VIRTUAL_BB", ires);
		}
	}
}

void
revert_loops_to_original(core::pose::Pose & pose, Loops loops)
{
	std::string sequence = core::sequence::read_fasta_file( option[ in::file::fasta ]()[1] )[1]->sequence();
	chemical::ResidueTypeSet const& restype_set( pose.residue(1).residue_type_set() );

	for (Size iloop=1; iloop<=loops.num_loop(); ++iloop) {
		/*
		if ( loops[iloop].start() == 1 ) {
			if (loops[iloop].stop() <=2) continue;
		}
		else if (loops[iloop].stop() - loops[iloop].start() <= 3) {
			continue;
		}
		else {
			Size nres(pose.total_residue());
			while (pose.residue_type(nres).is_virtual_residue()) {
				--nres;
			}
			if ( loops[iloop].stop() == nres ) {
				if (loops[iloop].start() >=nres-1) continue;
			}
		}
		*/
		for (Size ires=loops[iloop].start(); ires<=loops[iloop].stop(); ++ires) {
			utility::vector1< VariantType > variant_types = pose.residue_type(ires).variant_types();
			MutateResidue mutate_mover(ires, sequence[ires-1]);
			mutate_mover.apply(pose);
			for (Size i_var = 1; i_var <=variant_types.size(); ++i_var) {
				core::pose::add_variant_type_to_pose_residue(pose, variant_types[i_var], ires);
			}
		}
	}
}

void add_gap_constraints_to_pose(core::pose::Pose & pose, Size seq_distance_to_gap, Loops const & chunks, Real stdev=1.) {
	basic::Tracer TR( "pilot.yfsong.util" );
	using namespace ObjexxFCL::format;

	// add constraints
	Real boundary(999.);
	if (seq_distance_to_gap==1) {
		boundary=11.;
	}
	else if (seq_distance_to_gap==2) {
		boundary=18.;
	}
	else if (seq_distance_to_gap==3) {
		boundary=24.5;
	}
	else if (seq_distance_to_gap==4) {
		boundary=31.;
	}
	else {
		return;
	}

	for (Size i=1; i<chunks.num_loop(); ++i) {
		int gap_start = chunks[i].stop() - seq_distance_to_gap;
		int gap_end   = chunks[i+1].start() + seq_distance_to_gap;

		if ( gap_start >= 1 && gap_end <= pose.total_residue() ) {
			if (!pose.residue_type(gap_start).is_protein()) continue;
			if (!pose.residue_type(gap_end).is_protein()) continue;
			Size iatom = pose.residue_type(gap_start).atom_index("CA");
			Size jatom = pose.residue_type(gap_end).atom_index("CA");

			TR << "Add constraint to residue " << I(4,gap_start) << " and residue " << I(4,gap_end) << std::endl;
			pose.add_constraint(
								new core::scoring::constraints::AtomPairConstraint(
														   core::id::AtomID(iatom,gap_start),
														   core::id::AtomID(jatom,gap_end),
														   new core::scoring::constraints::BoundFunc( 0, boundary, stdev, "gap" ) )
								);
		}
	}
}

void
setup_startree(core::pose::Pose & pose) {
	std::string cut_point_decision("middle");

	basic::Tracer TR( "pilot.yfsong.util" );

	if (!option[ OptionKeys::in::file::psipred_ss2 ].user() ) {
		utility_exit_with_message("Error in reading psipred_ss2 file, is the -in:file:psipred_ss2 flag set correctly?");
	}

	bool check_psipred = set_secstruct_from_psipred_ss2(pose);
	assert (check_psipred);

	// Build the star fold tree, identify jumps
	if (option[challenge::ss].user()) {
		ss_chunks_pose_ = extract_secondary_structure_chunks( pose, option[challenge::ss]() );
	}
	else {
		ss_chunks_pose_ = extract_secondary_structure_chunks( pose );
	}
	ss_chunks_pose_.sequential_order();
	TR.Debug << "Target secondary chunks:" << std::endl;
	TR.Debug << ss_chunks_pose_ << std::endl;
	loops_pose_ = ss_chunks_pose_.invert(pose.total_residue());
	TR.Debug << "Target loops: " << pose.total_residue() << std::endl;
	TR.Debug << loops_pose_ << std::endl;
	TR.flush();

	if (option[challenge::virtual_loops]()) {
		set_loops_to_virt_ala(pose, loops_pose_);
	}

	// complete the chunks to cover the whole protein and customize cutpoints
	// cutpoints in the middle of the loop
	Loops chunks(ss_chunks_pose_);
	for (Size i=1; i<=chunks.num_loop(); ++i) {
		if ( cut_point_decision == "middle") {
			if (i==1) {
				chunks[i].set_start(1);
			}
			else {
				Size new_start = (ss_chunks_pose_[i-1].stop() + ss_chunks_pose_[i].start() + 1) / 2;
				//Size new_start = (ss_chunks_pose_[i-1].stop() + 1);
				chunks[i].set_start( new_start );
			}

			if (i==chunks.num_loop()) {
				chunks[i].set_stop(pose.total_residue());
			}
			else {
				Size new_stop = (ss_chunks_pose_[i].stop() + ss_chunks_pose_[i+1].start() - 1) / 2;
				chunks[i].set_stop( new_stop );
			}
		}
		else if ( cut_point_decision == "beginning" ) {
			if (i==1) {
				chunks[i].set_start(1);
			}
			else {
				Size new_start = (ss_chunks_pose_[i-1].stop() + 1);
				chunks[i].set_start( new_start );
			}

			if (i==chunks.num_loop()) {
				chunks[i].set_stop(pose.total_residue());
			}
		}
	}
	//add_gap_constraints_to_pose(pose, 4, chunks);

	TR.Debug << "Chunks: " << pose.total_residue() << std::endl;
	TR.Debug << chunks << std::endl;

	StarTreeBuilder builder;
	TR.Debug << pose.fold_tree() << std::endl;
	if (chunks.num_loop() > 0) {
		builder.set_up(chunks, &pose);
	}
	TR.Debug << pose.fold_tree() << std::endl;
	protocols::nonlocal::add_cutpoint_variants(&pose);
}

numeric::xyzVector<Real>
center_of_mass(core::pose::Pose const & pose) {
	int  nres( pose.total_residue() ), nAtms = 0;
	numeric::xyzVector<core::Real> massSum(0.,0.,0.), CoM;
	for ( int i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		if (rsd.aa() == core::chemical::aa_vrt) continue;

		for ( Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );
			massSum += atom.xyz();
			nAtms++;
		}
	}
	CoM = massSum / (core::Real)nAtms;
	return CoM;
}

void
translate_virt_to_CoM(core::pose::Pose & pose) {
	numeric::xyzVector<Real> CoM;
	CoM = center_of_mass(pose);
	numeric::xyzVector<Real> curr_pos = pose.residue(pose.total_residue()).xyz(1);
	numeric::xyzVector<Real> translation = CoM - curr_pos;

	basic::Tracer TR( "pilot.yfsong.util" );
	using namespace ObjexxFCL::format;
	TR.Debug << F(8,3,translation.x()) << F(8,3,translation.y()) << F(8,3,translation.z()) << std::endl;

	// apply transformation
	utility::vector1< core::id::AtomID > ids;
	utility::vector1< numeric::xyzVector<core::Real> > positions;
	for (Size iatom = 1; iatom <= pose.residue_type(pose.total_residue()).natoms(); ++iatom) {
		numeric::xyzVector<core::Real> atom_xyz = pose.xyz( core::id::AtomID(iatom,pose.total_residue()) );
		ids.push_back(core::id::AtomID(iatom,pose.total_residue()));
		positions.push_back(atom_xyz + translation);
	}
	pose.batch_set_xyz(ids,positions);
}

Loops loops() {
	return loops_pose_;
}

void
apply(core::pose::Pose & pose) {
	setup_startree(pose);

	// Initialize the structure
	MultiTemplateAlignChunkMover multi_align_mover(RG_, template_poses_, ss_chunks_pose_, all_chunks);
	multi_align_mover.apply(pose);
	translate_virt_to_CoM(pose);

	Size max_registry_shift = option[challenge::max_registry_shift]();
	MultiTemplateAlignChunkMoverOP random_align_mover(
					  new MultiTemplateAlignChunkMover(RG_, template_poses_, ss_chunks_pose_, random_chunk, max_registry_shift) );

	CustomFragmentMoverOP fragment_insertion_mover(
												new CustomFragmentMover(RG_, ss_chunks_pose_)
												   );

	for (Size i=1;i<=5;++i) {
		RandomMoverOP random_mover( new RandomMover() );
		Real weight = 0.1 * (Real)i;
		random_mover->add_mover(random_align_mover, weight);
		random_mover->add_mover(fragment_insertion_mover, 1.-weight);
		//random_mover->add_mover(new HelixMover(RG_));

		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		MoverOP sampling_mover = new RationalMonteCarlo(
														random_mover,
														scorefxn,
														200,
														2.0,
														true );

		sampling_mover->apply(pose);
	}

	if (option[challenge::revert_real_loops]()) {
		revert_loops_to_original(pose, loops_pose_);
	}

	basic::Tracer TR( "pilot.yfsong.util" );
	for (Size ires=1; ires<=pose.total_residue(); ++ires) {
		using namespace ObjexxFCL::format;
		TR.Debug << "Trial counter:" << I(4,ires) << I(8, random_align_mover->trial_counter(ires)) << std::endl;
	}
}

std::string
get_name() const {
	return "SampleSecondaryStructureAlignmentMover";
}

private:
numeric::random::RandomGenerator & RG_;
utility::vector1 < core::pose::PoseOP > template_poses_;

Loops ss_chunks_pose_;
Loops loops_pose_;

};

} // pilot
} // apps

#endif
