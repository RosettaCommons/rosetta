// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Align a random jump to template
/// @detailed
/// @author Yifan Song

#include <protocols/comparative_modeling/hybridize/FoldTreeHybridize.hh>
#include <protocols/comparative_modeling/hybridize/ChunkTrialMover.hh>
#include <protocols/comparative_modeling/hybridize/WeightedFragmentTrialMover.hh>
#include <protocols/comparative_modeling/hybridize/HybridizeFoldtreeDynamic.hh>

#include <core/import_pose/import_pose.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>

#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/kinematics_util.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <protocols/simple_moves/rational_mc/RationalMonteCarlo.hh>
#include <protocols/simple_moves/MutateResidue.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <numeric/random/DistributionSampler.hh>
#include <numeric/util.hh>

static numeric::random::RandomGenerator RG(42136);
static basic::Tracer TR( "protocols.comparative_modeling.hybridize.FoldTreeHybridize" );

namespace protocols {
namespace comparative_modeling {
namespace hybridize {

using namespace core;
using namespace chemical;
using namespace core::kinematics;
using namespace ObjexxFCL;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::loops;
using namespace numeric::model_quality;
using namespace id;
using namespace basic::options;
using namespace basic::options::OptionKeys;
	

FoldTreeHybridize::FoldTreeHybridize (
		core::Size const initial_template_index,
		utility::vector1 < core::pose::PoseCOP > const & template_poses,
		utility::vector1 < core::Real > const & template_wts,
		utility::vector1 < protocols::loops::Loops > const & template_chunks,
		utility::vector1 < protocols::loops::Loops > const & template_contigs,
		utility::vector1 < core::fragment::FragSetOP > & frag_libs)
{
	//initialize template structures
	initial_template_index_ = initial_template_index;
	template_poses_ = template_poses;
	template_wts_ = template_wts;
	template_chunks_ = template_chunks;
	template_contigs_ = template_contigs;

	// normalize weights
	core::Real weight_sum = 0.0;
	for (int i=1; i<=template_poses_.size(); ++i) weight_sum += template_wts_[i];
	for (int i=1; i<=template_poses_.size(); ++i) template_wts_[i] /= weight_sum;

	// abinitio frag9,frag3 flags
	frag_libs_ = frag_libs;
}
	
void
FoldTreeHybridize::set_loops_to_virt_ala(core::pose::Pose & pose, Loops loops)
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
FoldTreeHybridize::revert_loops_to_original(core::pose::Pose & pose, Loops loops)
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


Real
FoldTreeHybridize::gap_distance(Size Seq_gap)
{
	core::Real gap_torr_0( 4.0);
	core::Real gap_torr_1( 7.5);
	core::Real gap_torr_2(11.0);
	core::Real gap_torr_3(14.5);
	core::Real gap_torr_4(18.0);
	core::Real gap_torr_5(21.0);
	core::Real gap_torr_6(24.5);
	core::Real gap_torr_7(27.5);
	core::Real gap_torr_8(31.0);
	
	switch (Seq_gap) {
		case 0:
			return gap_torr_0;
			break;
		case 1:
			return gap_torr_1;
			break;
		case 2:
			return gap_torr_2;
			break;
		case 3:
			return gap_torr_3;
			case 4:
			return gap_torr_4;
			break;
		case 5:
			return gap_torr_5;
			break;
		case 6:
			return gap_torr_6;
			break;
		case 7:
			return gap_torr_7;
			break;
		case 8:
			return gap_torr_8;
			break;
		default:
			return 9999.;
	}
	return 9999.;
}

void FoldTreeHybridize::add_gap_constraints_to_pose(core::pose::Pose & pose, Loops const & chunks, int gap_edge_shift, Real stdev) {
	basic::Tracer TR("pilot.yfsong.util");
	using namespace ObjexxFCL::fmt;
	for (Size i=1; i<chunks.num_loop(); ++i) {
		int gap_start = chunks[i].stop()	+ gap_edge_shift;
		int gap_stop  = chunks[i+1].start() - gap_edge_shift;
		int gap_size = gap_stop - gap_start - 1;
		if (gap_size < 0) continue;
		if (gap_size > 8) continue;
		if (!pose.residue_type(gap_start).is_protein()) continue;
		if (!pose.residue_type(gap_stop ).is_protein()) continue;
		Size iatom = pose.residue_type(gap_start).atom_index("CA");
		Size jatom = pose.residue_type(gap_stop ).atom_index("CA");
		
		TR << "Add constraint to residue " << I(4,gap_start) << " and residue " << I(4,gap_stop) << std::endl;
		pose.add_constraint(
							new core::scoring::constraints::AtomPairConstraint(
																			   core::id::AtomID(iatom,gap_start),
																			   core::id::AtomID(jatom,gap_stop),
																			   new core::scoring::constraints::BoundFunc( 0, gap_distance(gap_size), stdev, "gap" ) )
							);
		
	}
}

/*
 void add_gap_constraints_to_pose_extend_gap(core::pose::Pose & pose, Size seq_distance_to_gap, Loops const & chunks, Real stdev=0.1) {
 basic::Tracer TR("pilot.yfsong.util");
 using namespace ObjexxFCL::fmt;
 
 // add constraints
 Real boundary(gap_distance(seq_distance_to_gap*2));
 
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
 */
	
protocols::loops::Loops FoldTreeHybridize::renumber_template_chunks(
												 protocols::loops::Loops & template_chunk,
												 core::pose::PoseCOP template_pose
												 )
{
	protocols::loops::Loops renumbered_template_chunks(template_chunk);
	for (core::Size ichunk = 1; ichunk<=template_chunk.num_loop(); ++ichunk) {
		Size seqpos_start_templ = template_chunk[ichunk].start();
		Size seqpos_start_target = template_pose->pdb_info()->number(seqpos_start_templ);
		renumbered_template_chunks[ichunk].set_start( seqpos_start_target );
		
		Size seqpos_stop_templ = template_chunk[ichunk].stop();
		Size seqpos_stop_target = template_pose->pdb_info()->number(seqpos_stop_templ);
		renumbered_template_chunks[ichunk].set_stop( seqpos_stop_target );
	}
	return renumbered_template_chunks;
}

/*
void
FoldTreeHybridize::setup_foldtree(core::pose::Pose & pose) {
	std::string cut_point_decision("middle");
	
	if (!option[ OptionKeys::in::file::psipred_ss2 ].user() ) {
		utility_exit_with_message("Error in reading psipred_ss2 file, is the -in:file:psipred_ss2 flag set correctly?");
	}
	
	bool check_psipred = set_secstruct_from_psipred_ss2(pose);
	assert (check_psipred);
	
	core::Size gap_size = RG.random_range(0,2);
	core::Size minimum_length_of_chunk_helix = 4;
	core::Size minimum_length_of_chunk_strand = 2;
	
	// extract secondary structure chunks from target psipred
	ss_chunks_pose_ = extract_secondary_structure_chunks( pose, option[cm::hybridize::ss](), gap_size, minimum_length_of_chunk_helix,minimum_length_of_chunk_strand );
	ss_chunks_pose_.sequential_order();
	
	// combine chunk definition from secondary structure and template
	utility::vector1 < protocols::loops::Loops > renumbered_template_contigs(template_contigs_);
	for (core::Size icontigs = 1; icontigs<=template_contigs_.size(); ++icontigs) {
		for (core::Size ichunk = 1; ichunk<=template_contigs_[icontigs].num_loop(); ++ichunk) {
			Size seqpos_start_templ = template_contigs_[icontigs][ichunk].start();
			Size seqpos_start_target = template_poses_[icontigs]->pdb_info()->number(seqpos_start_templ);
			renumbered_template_contigs[icontigs][ichunk].set_start( seqpos_start_target );

			Size seqpos_stop_templ = template_contigs_[icontigs][ichunk].stop();
			Size seqpos_stop_target = template_poses_[icontigs]->pdb_info()->number(seqpos_stop_templ);
			renumbered_template_contigs[icontigs][ichunk].set_stop( seqpos_stop_target );
		}
	}
	HybridizeFoldtreeDynamic foldtree_mover(renumbered_template_contigs);
	foldtree_mover.initialize(pose, ss_chunks_pose_); // combine psipred and template contigs for chunk initialization
	TR.Debug << pose.fold_tree() << std::endl;
}
*/

void
FoldTreeHybridize::setup_foldtree(core::pose::Pose & pose) {
    if (!option[ OptionKeys::in::file::psipred_ss2 ].user() ) {
		utility_exit_with_message("Error in reading psipred_ss2 file, is the -in:file:psipred_ss2 flag set correctly?");
	}
	bool check_psipred = set_secstruct_from_psipred_ss2(pose);
	assert (check_psipred);

	// combine:
	// (a) contigs in the current template
	utility::vector1< bool > template_mask( pose.total_residue(), false );
	protocols::loops::Loops my_chunks(template_chunks_[initial_template_index_]);

	for (core::Size icontig = 1; icontig<=template_chunks_[initial_template_index_].num_loop(); ++icontig) {
		Size seqpos_start_target = template_poses_[initial_template_index_]->pdb_info()->number(
			template_chunks_[initial_template_index_][icontig].start());
		my_chunks[icontig].set_start( seqpos_start_target );
		Size seqpos_stop_target = template_poses_[initial_template_index_]->pdb_info()->number(
			template_chunks_[initial_template_index_][icontig].stop());
		my_chunks[icontig].set_stop( seqpos_stop_target );

		for (Size j=seqpos_start_target; j<=seqpos_stop_target; ++j) template_mask[j] = true;
	}
	TR.Debug << "Chunks of initial template: " << initial_template_index_ << std::endl;
	TR.Debug << template_chunks_[initial_template_index_] << std::endl;
	TR.Debug << "Chunks from initial template: " << std::endl;
	TR.Debug << my_chunks << std::endl;

	if ( option[cm::hybridize::add_non_init_chunks]() ) {
	// (b) probabilistically sampled chunks from all other templates _outside_ these residues
	utility::vector1< std::pair< core::Real, protocols::loops::Loop > >  wted_insertions_to_consider;
	for (core::Size itempl = 1; itempl<=template_chunks_.size(); ++itempl) {
		if (itempl == initial_template_index_)
			continue;
		for (core::Size icontig = 1; icontig<=template_chunks_[itempl].num_loop(); ++icontig) {
			// remap
			Size seqpos_start_target = template_poses_[itempl]->pdb_info()->number(template_chunks_[itempl][icontig].start());
			Size seqpos_stop_target = template_poses_[itempl]->pdb_info()->number(template_chunks_[itempl][icontig].stop());

			bool uncovered = true;
			for (Size j=seqpos_start_target; j<=seqpos_stop_target && uncovered ; ++j) 
				uncovered &= !template_mask[j];

			if (uncovered) {
				protocols::loops::Loop new_loop = template_chunks_[itempl][icontig];
				new_loop.set_start( seqpos_start_target );
				new_loop.set_stop( seqpos_stop_target );
				wted_insertions_to_consider.push_back( std::pair< core::Real, protocols::loops::Loop >( template_wts_[itempl] , new_loop ) );
			}
		}
	}

	// (c) randomly shuffle, then add each with given prob
	TR.Debug << "Chunks from all template: " << std::endl;
	TR.Debug << my_chunks << std::endl;
	std::random_shuffle ( wted_insertions_to_consider.begin(), wted_insertions_to_consider.end() );
	for (int i=1; i<=wted_insertions_to_consider.size(); ++i) {
		// ensure the insert is still valid
		bool uncovered = true;
		for (Size j=wted_insertions_to_consider[i].second.start(); j<=wted_insertions_to_consider[i].second.stop() && uncovered ; ++j) 
			uncovered &= !template_mask[j];

		if (!uncovered) continue;

		core::Real selector = numeric::random::uniform();
		TR << "Consider " << wted_insertions_to_consider[i].second.start() << "," << wted_insertions_to_consider[i].second.stop() << std::endl;
		if (selector <= wted_insertions_to_consider[i].first) {
			TR << " ====> taken!" << std::endl;
			my_chunks.add_loop( wted_insertions_to_consider[i].second );
			for (Size j=wted_insertions_to_consider[i].second.start(); j<=wted_insertions_to_consider[i].second.stop(); ++j) template_mask[j] = true;
		}
	}
	}
	
	my_chunks.sequential_order();
	TR << "Chunks used for foldtree setup: " << std::endl;
	TR << my_chunks << std::endl;

	HybridizeFoldtreeDynamic foldtree_mover;
	foldtree_mover.initialize(pose, my_chunks);
	TR << pose.fold_tree() << std::endl;
}


numeric::xyzVector<Real>
FoldTreeHybridize::center_of_mass(core::pose::Pose const & pose) {
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
FoldTreeHybridize::translate_virt_to_CoM(core::pose::Pose & pose) {
	numeric::xyzVector<Real> CoM;
	CoM = center_of_mass(pose);
	numeric::xyzVector<Real> curr_pos = pose.residue(pose.total_residue()).xyz(1);
	numeric::xyzVector<Real> translation = CoM - curr_pos;
	
	basic::Tracer TR("pilot.yfsong.util");
	using namespace ObjexxFCL::fmt;
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

/*
Loops FoldTreeHybridize::loops() {
	return loops_pose_;	
}
*/

utility::vector1< core::Real > FoldTreeHybridize::get_residue_weights_from_loops(core::pose::Pose & pose) {
	utility::vector1< core::Real > residue_weights(pose.total_residue());
	protocols::loops::Loops renumbered_template_chunks = renumber_template_chunks(
		template_contigs_[initial_template_index_], template_poses_[initial_template_index_]);
	TR.Debug << "Insert fragment for these residues:";
	for ( Size ires=1; ires<= pose.total_residue(); ++ires ) {
		if (! renumbered_template_chunks.has(ires) ) {
			residue_weights[ires] = 1.0;
			TR.Debug << " " << ires;
		}
		else {
			residue_weights[ires] = option[cm::hybridize::frag_weight_aligned]();
		}
	}
	TR.Debug << std::endl;
	return residue_weights;
}

void FoldTreeHybridize::backup_original_foldtree(core::pose::Pose const & pose) {
	orig_ft_ = pose.conformation().fold_tree();
	orig_n_residue_ = pose.total_residue();
}

void FoldTreeHybridize::restore_original_foldtree(core::pose::Pose & pose) {
	if (pose.total_residue() > orig_n_residue_) {
		pose.conformation().delete_residue_range_slow(orig_n_residue_+1, pose.total_residue());
	}
	pose.conformation().fold_tree( orig_ft_ );
}

void
FoldTreeHybridize::apply(core::pose::Pose & pose) {
	backup_original_foldtree(pose);
	setup_foldtree(pose);
	
	// Initialize the structure
	bool use_random_template = false;
	ChunkTrialMover initialize_chunk_mover(template_poses_, template_chunks_, ss_chunks_pose_, use_random_template, all_chunks);
	initialize_chunk_mover.set_template(initial_template_index_);
	initialize_chunk_mover.apply(pose);
	translate_virt_to_CoM(pose);
	
	//pose.dump_pdb("after_init.pdb");
	
	use_random_template = true;
	Size max_registry_shift = option[cm::hybridize::max_registry_shift]();
	ChunkTrialMoverOP random_sample_chunk_mover(
		new ChunkTrialMover(template_poses_, template_chunks_, ss_chunks_pose_, use_random_template, random_chunk, max_registry_shift)
	);
	
	utility::vector1< core::Real > residue_weights( get_residue_weights_from_loops(pose) );
	WeightedFragmentTrialMoverOP fragment_trial_mover(
		new WeightedFragmentTrialMover(frag_libs_, residue_weights)
	);
	
	// just fragment mover
	(*scorefxn_)(pose);
	protocols::moves::MonteCarloOP mc_frag = new protocols::moves::MonteCarlo( pose, *scorefxn_, 2.0 );
	for (Size i = 1; i<= 2000; ++i) {
		fragment_trial_mover->apply(pose);
		(*scorefxn_)(pose);
		mc_frag->boltzmann(pose);
	}
	mc_frag->recover_low(pose);
	
	//utility_exit_with_message("debugging");

	for (Size i=1;i<=4;++i) {
		if (i==3) scorefxn_->set_weight( core::scoring::linear_chainbreak, 0.5 );
		if (i==4) scorefxn_->set_weight( core::scoring::linear_chainbreak, 2.0 );

		RandomMoverOP random_mover( new RandomMover() );
		Real weight = 0.05 * (Real)i;
		random_mover->add_mover(random_sample_chunk_mover, 1. - weight);
		random_mover->add_mover(fragment_trial_mover, weight);
		//random_mover->add_mover(new HelixMover(RG));

		(*scorefxn_)(pose);
		protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( pose, *scorefxn_, 2.0 );

		for (Size i = 1; i<= 500; ++i) {
			random_mover->apply(pose);
			(*scorefxn_)(pose);
			mc->boltzmann(pose);
		}
		mc->recover_low(pose);
	}

	restore_original_foldtree(pose);

	basic::Tracer TR("pilot.yfsong.util");
	for (Size ires=1; ires<=pose.total_residue(); ++ires) {
		using namespace ObjexxFCL::fmt;
		TR.Debug << "Trial counter:" << I(4,ires) << I(8, random_sample_chunk_mover->trial_counter(ires)) << std::endl;
	}
}
	
std::string FoldTreeHybridize::get_name() const
{
	return "FoldTreeHybridize";
}

} // hybridize 
} // comparative_modeling 
} // protocols

