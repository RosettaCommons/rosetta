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
#include <protocols/comparative_modeling/hybridize/util.hh>

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
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/kinematics_util.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

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
	

FoldTreeHybridize::FoldTreeHybridize() :
foldtree_mover_()
{
	init();
}

FoldTreeHybridize::FoldTreeHybridize (
		core::Size const initial_template_index,
		utility::vector1 < core::pose::PoseCOP > const & template_poses,
		utility::vector1 < core::Real > const & template_wts,
		utility::vector1 < protocols::loops::Loops > const & template_chunks,
		utility::vector1 < protocols::loops::Loops > const & template_contigs,
		core::fragment::FragSetOP fragments3_in,
		core::fragment::FragSetOP fragments9_in )
{
	init();

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
	frag_libs3_.push_back(fragments3_in);
	frag_libs9_.push_back(fragments9_in);
	frag_libs_.push_back(fragments3_in);
	frag_libs_.push_back(fragments9_in);
}
	
void
FoldTreeHybridize::init() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	increase_cycles_ = option[cm::hybridize::stage1_increase_cycles]();
	add_non_init_chunks_ = option[cm::hybridize::add_non_init_chunks]();
	frag_weight_aligned_ = option[cm::hybridize::frag_weight_aligned]();
	max_registry_shift_ = option[cm::hybridize::max_registry_shift]();
	frag_insertion_weight_ = 0.5;
	domain_assembly_ = false;
	// default scorefunction
	set_scorefunction ( core::scoring::ScoreFunctionFactory::create_score_function( "score3" ) );
}

void
FoldTreeHybridize::set_loops_to_virt_ala(core::pose::Pose & pose, Loops loops)
{
	chemical::ResidueTypeSet const& restype_set( pose.residue(1).residue_type_set() );
	
	for (Size iloop=1; iloop<=loops.num_loop(); ++iloop) {
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
			return gap_torr_0; break;
		case 1:
			return gap_torr_1; break;
		case 2:
			return gap_torr_2; break;
		case 3:
			return gap_torr_3; break;
		case 4:
			return gap_torr_4; break;
		case 5:
			return gap_torr_5; break;
		case 6:
			return gap_torr_6; break;
		case 7:
			return gap_torr_7; break;
		case 8:
			return gap_torr_8; break;
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


void
FoldTreeHybridize::setup_foldtree(core::pose::Pose & pose) {
	//if (!option[ OptionKeys::in::file::psipred_ss2 ].user() ) {
	//	utility_exit_with_message("Error in reading psipred_ss2 file, is the -in:file:psipred_ss2 flag set correctly?");
	//}
	//bool check_psipred = set_secstruct_from_psipred_ss2(pose);
	//assert (check_psipred);

	// combine:
	// (a) contigs in the current template

	core::Size nres = pose.total_residue();

	//symmetry
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
		nres = symm_info->num_independent_residues();
	}

	utility::vector1< bool > template_mask( nres, false );
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

	if ( add_non_init_chunks_ || domain_assembly_ ) {
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

	foldtree_mover_.initialize(pose, my_chunks);
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

utility::vector1< core::Real > FoldTreeHybridize::get_residue_weights_from_loops(core::pose::Pose & pose) {
	utility::vector1< core::Real > residue_weights(pose.total_residue());
	TR.Debug << "Insert fragment for these residues:";
	for ( Size ires=1; ires<= pose.total_residue(); ++ires ) {
		if (domain_assembly_) {
			bool residue_in_template = false;
			for (Size i_template=1; i_template<=template_poses_.size(); ++i_template) {
				protocols::loops::Loops renumbered_template_chunks = renumber_template_chunks(
									  template_contigs_[i_template], template_poses_[i_template]);
				if (renumbered_template_chunks.has(ires)) {
					residue_in_template = true;
					break;
				}
			}
			if (! residue_in_template ) {
				residue_weights[ires] = 1.0;
				TR.Debug << " " << ires;
			}
			else {
				residue_weights[ires] = frag_weight_aligned_;
			}
		}
		else {
			protocols::loops::Loops renumbered_template_chunks
			= renumber_template_chunks(
									   template_contigs_[initial_template_index_], template_poses_[initial_template_index_]);
			
			if (! renumbered_template_chunks.has(ires) ) {
				residue_weights[ires] = 1.0;
				TR.Debug << " " << ires;
			}
			else {
				residue_weights[ires] = frag_weight_aligned_;
			}
		}
	}
	// reset linker fragment insertion weights
	if (domain_assembly_) {
	for (Size i_template=1; i_template<=template_poses_.size(); ++i_template) {
		int coverage_start = template_poses_[i_template]->pdb_info()->number(1);
		int coverage_end   = template_poses_[i_template]->pdb_info()->number(1);
		for (Size ires = 1; ires <= template_poses_[i_template]->total_residue(); ++ires) {
			if (template_poses_[i_template]->pdb_info()->number(ires) < coverage_start) {
				coverage_start = template_poses_[i_template]->pdb_info()->number(ires);
			}
			if (template_poses_[i_template]->pdb_info()->number(ires) > coverage_end) {
				coverage_end = template_poses_[i_template]->pdb_info()->number(ires);
			}
		}
		
		for (int shift = -5; shift<=5; ++shift) {
			int ires = coverage_start + shift;
			if (ires >= 1 && ires <= pose.total_residue()) {
				residue_weights[ires] = 1.;
				TR.Debug << " " << ires;
			}

			ires = coverage_end + shift;
			if (ires >= 1 && ires <= pose.total_residue()) {
				residue_weights[ires] = 1.;
				TR.Debug << " " << ires;
			}
		}
	}
	}
	TR.Debug << std::endl;
	return residue_weights;
}

void FoldTreeHybridize::backup_original_foldtree(core::pose::Pose const & /*pose*/ ) {
	//orig_ft_ = pose.conformation().fold_tree();
	//orig_n_residue_ = pose.total_residue();
}

void FoldTreeHybridize::restore_original_foldtree(core::pose::Pose & pose) {
	foldtree_mover_.reset(pose);
	/*
	if (pose.total_residue() > orig_n_residue_) {
		pose.conformation().delete_residue_range_slow(orig_n_residue_+1, pose.total_residue());
	}
	protocols::loops::remove_cutpoint_variants( pose );
	pose.conformation().fold_tree( orig_ft_ );
	*/
}

void
FoldTreeHybridize::setup_scorefunctions( 
		core::scoring::ScoreFunctionOP score0,
		core::scoring::ScoreFunctionOP score1,
		core::scoring::ScoreFunctionOP score2,
		core::scoring::ScoreFunctionOP score5,
		core::scoring::ScoreFunctionOP score3) {
	core::Real lincb_orig = scorefxn_->get_weight( core::scoring::linear_chainbreak );
	core::Real cb_orig = scorefxn_->get_weight( core::scoring::chainbreak );
	core::Real cst_orig = scorefxn_->get_weight( core::scoring::atom_pair_constraint );
	
	score0->reset();
	score0->set_weight( core::scoring::vdw, 0.1*scorefxn_->get_weight( core::scoring::vdw ) );

	score1->reset();
	score2->set_weight( core::scoring::linear_chainbreak, 0.1*lincb_orig );
	score1->set_weight( core::scoring::atom_pair_constraint, 0.1*cst_orig );
	score1->set_weight( core::scoring::vdw, scorefxn_->get_weight( core::scoring::vdw ) );
	score1->set_weight( core::scoring::env, scorefxn_->get_weight( core::scoring::env ) );
	score1->set_weight( core::scoring::cen_env_smooth, scorefxn_->get_weight( core::scoring::cen_env_smooth ) );
	score1->set_weight( core::scoring::pair, scorefxn_->get_weight( core::scoring::pair ) );
	score1->set_weight( core::scoring::cen_pair_smooth, scorefxn_->get_weight( core::scoring::cen_pair_smooth ) );
	score1->set_weight( core::scoring::hs_pair, scorefxn_->get_weight( core::scoring::hs_pair ) );
	score1->set_weight( core::scoring::ss_pair, 0.3*scorefxn_->get_weight( core::scoring::ss_pair ) );
	score1->set_weight( core::scoring::sheet, scorefxn_->get_weight( core::scoring::sheet ) );
	//STRAND_STRAND_WEIGHTS 1 11
	core::scoring::methods::EnergyMethodOptions score1_options(score1->energy_method_options());
	score1_options.set_strand_strand_weights(1,11);
	score1->set_energy_method_options(score1_options);
	

	score2->reset();
	score2->set_weight( core::scoring::linear_chainbreak, 0.25*lincb_orig );
	score2->set_weight( core::scoring::atom_pair_constraint, 0.25*cst_orig );
	score2->set_weight( core::scoring::vdw, scorefxn_->get_weight( core::scoring::vdw ) );
	score2->set_weight( core::scoring::env, scorefxn_->get_weight( core::scoring::env ) );
	score2->set_weight( core::scoring::cen_env_smooth, scorefxn_->get_weight( core::scoring::cen_env_smooth ) );
	score2->set_weight( core::scoring::cbeta, 0.25*scorefxn_->get_weight( core::scoring::cbeta ) );
	score2->set_weight( core::scoring::cbeta_smooth, 0.25*scorefxn_->get_weight( core::scoring::cbeta_smooth ) );
	score2->set_weight( core::scoring::cenpack, 0.5*scorefxn_->get_weight( core::scoring::cenpack ) );
	score2->set_weight( core::scoring::cenpack_smooth, 0.5*scorefxn_->get_weight( core::scoring::cenpack_smooth ) );
	score2->set_weight( core::scoring::pair, scorefxn_->get_weight( core::scoring::pair ) );
	score2->set_weight( core::scoring::cen_pair_smooth, scorefxn_->get_weight( core::scoring::cen_pair_smooth ) );
	score2->set_weight( core::scoring::hs_pair, scorefxn_->get_weight( core::scoring::hs_pair ) );
	score2->set_weight( core::scoring::ss_pair, 0.3*scorefxn_->get_weight( core::scoring::ss_pair ) );
	score2->set_weight( core::scoring::sheet, scorefxn_->get_weight( core::scoring::sheet ) );
	//STRAND_STRAND_WEIGHTS 1 6
	core::scoring::methods::EnergyMethodOptions score2_options(score1->energy_method_options());
	score2_options.set_strand_strand_weights(1,6);
	score2->set_energy_method_options(score2_options);

	score5->reset();
	score5->set_weight( core::scoring::linear_chainbreak, 0.25*lincb_orig );
	score5->set_weight( core::scoring::atom_pair_constraint, 0.25*cst_orig );
	score5->set_weight( core::scoring::vdw, scorefxn_->get_weight( core::scoring::vdw ) );
	score5->set_weight( core::scoring::env, scorefxn_->get_weight( core::scoring::env ) );
	score5->set_weight( core::scoring::cen_env_smooth, scorefxn_->get_weight( core::scoring::cen_env_smooth ) );
	score5->set_weight( core::scoring::cbeta, 0.25*scorefxn_->get_weight( core::scoring::cbeta ) );
	score5->set_weight( core::scoring::cbeta_smooth, 0.25*scorefxn_->get_weight( core::scoring::cbeta_smooth ) );
	score5->set_weight( core::scoring::cenpack, 0.5*scorefxn_->get_weight( core::scoring::cenpack ) );
	score5->set_weight( core::scoring::cenpack_smooth, 0.5*scorefxn_->get_weight( core::scoring::cenpack_smooth ) );
	score5->set_weight( core::scoring::pair, scorefxn_->get_weight( core::scoring::pair ) );
	score5->set_weight( core::scoring::cen_pair_smooth, scorefxn_->get_weight( core::scoring::cen_pair_smooth ) );
	score5->set_weight( core::scoring::hs_pair, scorefxn_->get_weight( core::scoring::hs_pair ) );
	score5->set_weight( core::scoring::ss_pair, 0.3*scorefxn_->get_weight( core::scoring::ss_pair ) );
	score5->set_weight( core::scoring::sheet, scorefxn_->get_weight( core::scoring::sheet ) );
	//STRAND_STRAND_WEIGHTS 1 11
	core::scoring::methods::EnergyMethodOptions score5_options(score1->energy_method_options());
	score5_options.set_strand_strand_weights(1,11);
	score5->set_energy_method_options(score5_options);

	score3 = scorefxn_->clone();
}


void
FoldTreeHybridize::apply(core::pose::Pose & pose) {
	//backup_original_foldtree(pose);
	setup_foldtree(pose);

	setup_centroid_constraints( pose, template_poses_, template_wts_, cst_file_ );
	
	// Initialize the structure
	bool use_random_template = false;
	ChunkTrialMover initialize_chunk_mover(template_poses_, template_chunks_, ss_chunks_pose_, use_random_template, all_chunks);
	initialize_chunk_mover.set_template(initial_template_index_);
	initialize_chunk_mover.apply(pose);
	translate_virt_to_CoM(pose);
	
	use_random_template = true;
	Size max_registry_shift = max_registry_shift_;
	ChunkTrialMoverOP random_sample_chunk_mover(
		new ChunkTrialMover(template_poses_, template_chunks_, ss_chunks_pose_, use_random_template, random_chunk, max_registry_shift) );

	utility::vector1< core::Real > residue_weights( get_residue_weights_from_loops(pose) );
	utility::vector1< core::Size > jump_anchors = foldtree_mover_.get_anchors();

	// ab initio ramping up of weights
	// set up scorefunctions
	core::scoring::ScoreFunctionOP score0=scorefxn_->clone(),
	                               score1=scorefxn_->clone(),
	                               score2=scorefxn_->clone(),
	                               score5=scorefxn_->clone(),
	                               score3=scorefxn_->clone();
	setup_scorefunctions( score0, score1, score2, score5, score3 );

	// set up movers
	RandomMoverOP stage0mover( new RandomMover() );
	if ( frag_insertion_weight_ < 1. ) {
		stage0mover->add_mover(random_sample_chunk_mover, 1.-frag_insertion_weight_);
	}
	if ( frag_insertion_weight_ > 0. ) {
		core::Real sum_weight = 0.;
		for ( Size ires=1; ires<= residue_weights.size(); ++ires ) sum_weight += residue_weights[ires];
		if (sum_weight > 1e-6) {
			WeightedFragmentTrialMoverOP fragment_trial_mover(
															  new WeightedFragmentTrialMover(frag_libs_, residue_weights, jump_anchors) );
			stage0mover->add_mover(fragment_trial_mover, frag_insertion_weight_);
		}
	}
	
	// stage 1
 	protocols::moves::MonteCarloOP mc1 = new protocols::moves::MonteCarlo( pose, *score0, 2.0 );
 	for (int i=1; i<=2000*increase_cycles_; ++i) {
 		stage0mover->apply(pose);
 	 	(*score0)(pose); mc1->boltzmann(pose);
 	}
 	mc1->show_scores();
 	mc1->show_counters();
 	mc1->recover_low(pose);

	// stage 2
	protocols::moves::MonteCarloOP mc2 = new protocols::moves::MonteCarlo( pose, *score1, 2.0 );
	for (int i=1; i<=2000*increase_cycles_; ++i) {
		stage0mover->apply(pose);
	 	(*score1)(pose); mc2->boltzmann(pose);
	}
	mc2->show_scores();
	mc2->show_counters();
	mc2->recover_low(pose);

	// stage 3
	for (int nmacro=1; nmacro<=10; ++nmacro) {
		protocols::moves::MonteCarloOP mc3a = new protocols::moves::MonteCarlo( pose, *score2, 2.0 );
		for (int i=1; i<=200*increase_cycles_; ++i) {
			stage0mover->apply(pose);
			(*score2)(pose); mc3a->boltzmann(pose);
		}
		mc3a->show_scores();
		mc3a->show_counters();
		mc3a->recover_low(pose);

		protocols::moves::MonteCarloOP mc3b = new protocols::moves::MonteCarlo( pose, *score5, 2.0 );
		for (int i=1; i<=200*increase_cycles_; ++i) {
			stage0mover->apply(pose);
			(*score5)(pose); mc3b->boltzmann(pose);
		}
		mc3b->show_scores();
		mc3b->show_counters();
		mc3b->recover_low(pose);
	}

	// stage 4 -- ramp up chainbreak
	core::Real lincb_orig = scorefxn_->get_weight( core::scoring::linear_chainbreak );
	if (lincb_orig == 0) lincb_orig = 1.0;  // force chainbreak
	score3->set_weight( core::scoring::linear_chainbreak, 0.0 );
	for (int nmacro=1; nmacro<=4; ++nmacro) {
		score3->set_weight( core::scoring::linear_chainbreak, 0.25 * nmacro * lincb_orig );

		protocols::moves::MonteCarloOP mc5 = new protocols::moves::MonteCarlo( pose, *score3, 2.0 );
		for (int i=1; i<=(core::Size)(500*increase_cycles_); ++i) {
			stage0mover->apply( pose );
			(*score3)(pose); mc5->boltzmann(pose);
		}
		mc5->show_scores();
		mc5->show_counters();
		mc5->recover_low(pose);
	}

	pose.remove_constraints();
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

