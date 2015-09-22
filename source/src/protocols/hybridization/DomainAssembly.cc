// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Domain Assembly
/// @details
/// @author Yifan Song

#include <protocols/hybridization/DomainAssembly.hh>
#include <protocols/hybridization/TMalign.hh>
#include <protocols/hybridization/util.hh>
#include <protocols/docking/DockingLowRes.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/types.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/Energies.hh>

//numeric
#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>

// options
#include <basic/options/option.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>

//utility
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.hybridization.DomainAssembly" );

namespace protocols {
//namespace comparative_modeling {
namespace hybridization {

bool TMalign_poses(core::pose::Pose & aligned_pose,
	core::pose::Pose const & ref_pose,
	std::list <Size> const & residue_list,
	std::list <Size> const & ref_residue_list) {
	TMalign tm_align;
	std::string seq_pose, seq_ref, aligned;
	core::id::AtomID_Map< core::id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, aligned_pose, core::id::BOGUS_ATOM_ID );
	core::Size n_mapped_residues=0;
	core::Size normalize_length = aligned_pose.total_residue() < ref_pose.total_residue() ? aligned_pose.total_residue() : ref_pose.total_residue();

	tm_align.apply(aligned_pose, ref_pose, residue_list, ref_residue_list);
	tm_align.alignment2AtomMap(aligned_pose, ref_pose, residue_list, ref_residue_list, n_mapped_residues, atom_map);
	tm_align.alignment2strings(seq_pose, seq_ref, aligned);
	core::Real TMscore = tm_align.TMscore(normalize_length);

	using namespace ObjexxFCL::format;
	TR << "Align domain with TMscore of " << F(8,3,TMscore) << std::endl;
	TR << seq_pose << std::endl;
	TR << aligned << std::endl;
	TR << seq_ref << std::endl;

	std::list <Size> full_residue_list;
	for ( Size ires=1; ires<=aligned_pose.total_residue(); ++ires ) {
		full_residue_list.push_back(ires);
	}

	if ( n_mapped_residues >= 6 ) {
		utility::vector1< core::Real > aln_cutoffs;
		aln_cutoffs.push_back(6);
		aln_cutoffs.push_back(4);
		aln_cutoffs.push_back(3);
		aln_cutoffs.push_back(2);
		aln_cutoffs.push_back(1.5);
		aln_cutoffs.push_back(1);
		core::Real min_coverage = 0.2;
		partial_align(aligned_pose, ref_pose, atom_map, full_residue_list, true, aln_cutoffs, min_coverage); // iterate_convergence = true
		return true;
	}

	return false;
}

DomainAssembly::DomainAssembly(
	utility::vector1 < core::pose::PoseOP > & poses,
	utility::vector1 < core::Real > & domain_assembly_weights
)
{
	// initialize the order of assembly
	utility::vector1 <core::Size> pose_index_to_add;
	utility::vector1 <core::Size> pose_index_added;
	for ( Size i=1; i<=poses.size(); ++i ) {
		pose_index_to_add.push_back(i);
	}
	while ( pose_index_to_add.size() > 0 ) {
		core::Real sum_weight(0.0);
		utility::vector1 < core::Real > weights_work;
		for ( Size i=1; i<=pose_index_to_add.size(); ++i ) {
			Size k_pose = pose_index_to_add[i];
			sum_weight += domain_assembly_weights[k_pose];
			weights_work.push_back(domain_assembly_weights[k_pose]);
		}
		Size ipose;
		if ( sum_weight > 1e-6 ) {
			numeric::random::WeightedSampler weighted_sampler_;
			weighted_sampler_.weights(weights_work);
			ipose = pose_index_to_add[weighted_sampler_.random_sample(numeric::random::rg())];
		} else {
			ipose = pose_index_to_add[numeric::random::rg().random_range(1,pose_index_to_add.size())];
		}

		poses_.push_back(poses[ipose]);
		pose_index_added.push_back(ipose);

		pose_index_to_add.clear();
		for ( Size ipose=1; ipose<=poses.size(); ++ipose ) {
			bool added = false;
			for ( Size j=1; j<=pose_index_added.size(); ++j ) {
				if ( ipose == pose_index_added[j] ) {
					added = true;
					break;
				}
			}
			if ( !added ) {
				pose_index_to_add.push_back(ipose);
			}
		}
	}
}

void remove_residues(core::pose::Pose & pose,
	utility::vector1 <int> const resnum_list,
	utility::vector1<int> & remaining_resnum) {

	remaining_resnum.clear();
	utility::vector1<Size> delete_list;
	for ( Size ires=1; ires <= pose.total_residue(); ++ires ) {
		bool deleting = false;
		for ( Size jres=1; jres<=resnum_list.size(); ++jres ) {
			if ( resnum_list[jres] == pose.pdb_info()->number(ires) ) {
				deleting = true;
				delete_list.push_back(ires);
				break;
			}
		}
		if ( !deleting ) {
			remaining_resnum.push_back(pose.pdb_info()->number(ires));
		}
	}

	for ( Size i=delete_list.size(); i>=1; --i ) {
		pose.conformation().delete_residue_slow(delete_list[i]);
	}
}

utility::vector1<Size> find_uncovered_residues (core::pose::Pose const & pose,
	utility::vector1 <Size> const covered_residues) {
	utility::vector1<Size> uncovered_residues;
	for ( Size ires = 1; ires <= pose.total_residue(); ++ires ) {
		if ( pose.pdb_info() != 0 ) {
			bool covered(false);
			for ( Size i=1; i<=covered_residues.size(); ++i ) {
				if ( pose.pdb_info()->number(ires) == (int)covered_residues[i] ) {
					covered = true;
					break;
				}
			}
			if ( !covered ) {
				uncovered_residues.push_back(ires);
			}
		}
	}
	return uncovered_residues;
}

core::Real
gap_distance(Size Seq_gap)
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
	case 0 :
		return gap_torr_0; break;
	case 1 :
		return gap_torr_1; break;
	case 2 :
		return gap_torr_2; break;
	case 3 :
		return gap_torr_3; break;
	case 4 :
		return gap_torr_4; break;
	case 5 :
		return gap_torr_5; break;
	case 6 :
		return gap_torr_6; break;
	case 7 :
		return gap_torr_7; break;
	case 8 :
		return gap_torr_8; break;
	default :
		return 9999.;
	}
	return 9999.;
}

void
DomainAssembly::run()
{
	if ( poses_.size() < 2 ) return;

	utility::vector1<int> coverage_start;
	utility::vector1<int> coverage_end;

	int range_start = poses_[1]->pdb_info()->number(1);
	int range_end   = poses_[1]->pdb_info()->number(1);
	for ( Size ires = 1; ires <= poses_[1]->total_residue(); ++ires ) {
		if ( poses_[1]->pdb_info()->number(ires) < range_start ) {
			range_start = poses_[1]->pdb_info()->number(ires);
		}
		if ( poses_[1]->pdb_info()->number(ires) > range_end ) {
			range_end = poses_[1]->pdb_info()->number(ires);
		}
	}
	coverage_start.push_back(range_start);
	coverage_end.push_back(range_end);

	for ( Size ipose = 2; ipose <= poses_.size(); ++ipose ) {
		// check if the current pose is overlapped with any previous poses
		range_start = poses_[ipose]->pdb_info()->number(1);
		range_end   = poses_[ipose]->pdb_info()->number(1);
		for ( Size ires = 1; ires <= poses_[ipose]->total_residue(); ++ires ) {
			if ( poses_[ipose]->pdb_info()->number(ires) < range_start ) {
				range_start = poses_[ipose]->pdb_info()->number(ires);
			}
			if ( poses_[ipose]->pdb_info()->number(ires) > range_end ) {
				range_end = poses_[ipose]->pdb_info()->number(ires);
			}
		}
		coverage_start.push_back(range_start);
		coverage_end.push_back(range_end);

		bool align_success(false);
		for ( Size jpose = 1; jpose < ipose; ++jpose ) {
			int covered_size = coverage_end[jpose] - coverage_start[jpose] + 1;
			int overlap_start = range_start > coverage_start[jpose] ? range_start:coverage_start[jpose];
			int overlap_end   = range_end   < coverage_end[jpose]   ? range_end:coverage_end[jpose];
			int overlap = overlap_end - overlap_start; // end could be smaller than start
			Size normalize_length = covered_size < (int)poses_[ipose]->total_residue() ? covered_size:poses_[ipose]->total_residue();

			// if overlap, use TMalign to orient poses_[ipose]
			if ( overlap > 0.3 * normalize_length ) {
				// collect residues in ipose
				std::list <Size> i_residue_list; // residue numbers in the overlapped region
				std::list <Size> full_residue_list; // all residue numbers in ipose, used for transformation after alignment
				for ( Size ires = 1; ires <= poses_[ipose]->total_residue(); ++ires ) {
					full_residue_list.push_back(ires);
					if ( poses_[ipose]->pdb_info()->number(ires) >= overlap_start &&
							poses_[ipose]->pdb_info()->number(ires) <= overlap_end ) {
						i_residue_list.push_back(ires);
					}

				}

				// collect residues in the pose to be aligned to
				std::list <Size> j_residue_list;
				for ( Size jres = 1; jres <= poses_[jpose]->total_residue(); ++jres ) {
					if ( poses_[jpose]->pdb_info()->number(jres) >= overlap_start &&
							poses_[jpose]->pdb_info()->number(jres) <= overlap_end ) {
						j_residue_list.push_back(jres);
					}
				}
				if ( j_residue_list.size() < 0.3 * normalize_length ) continue;

				TMalign tm_align;
				std::string seq_pose, seq_ref, aligned;
				core::id::AtomID_Map< core::id::AtomID > atom_map;
				core::pose::initialize_atomid_map( atom_map, *poses_[ipose], core::id::BOGUS_ATOM_ID );
				core::Size n_mapped_residues=0;

				tm_align.apply(*poses_[ipose], *poses_[jpose], i_residue_list, j_residue_list);
				tm_align.alignment2AtomMap(*poses_[ipose], *poses_[jpose], i_residue_list, j_residue_list, n_mapped_residues, atom_map);
				tm_align.alignment2strings(seq_pose, seq_ref, aligned);
				core::Real TMscore = tm_align.TMscore(normalize_length);

				using namespace ObjexxFCL::format;
				TR << "Align domain with TMscore of " << F(8,3,TMscore) << std::endl;
				TR << seq_pose << std::endl;
				TR << aligned << std::endl;
				TR << seq_ref << std::endl;

				if ( TMscore > 0.35 ) {
					if ( n_mapped_residues >= 6 ) {
						utility::vector1< core::Real > aln_cutoffs;
						aln_cutoffs.push_back(6);
						aln_cutoffs.push_back(4);
						aln_cutoffs.push_back(3);
						aln_cutoffs.push_back(2);
						aln_cutoffs.push_back(1.5);
						aln_cutoffs.push_back(1);
						core::Real min_coverage = 0.2;
						partial_align(*poses_[ipose], *poses_[jpose], atom_map, full_residue_list, true, aln_cutoffs, min_coverage); // iterate_convergence = true
						align_success = true;
					}
				}

				if ( align_success ) break;
			}
		}

		// if not overlap, use docking to align to the largest non-overlapped region
		if ( !align_success ) {
			using namespace ObjexxFCL::format;
			// construct a pose for domain assembly
			core::pose::PoseOP full_length_pose;
			//Size first_domain_end;
			utility::vector1<int> covered_resnum;
			utility::vector1<int> resnum;
			for ( Size ires=1; ires<=poses_[ipose]->total_residue(); ++ires ) {
				covered_resnum.push_back(poses_[ipose]->pdb_info()->number(ires));
			}
			for ( Size jpose = 1; jpose < ipose; ++jpose ) {
				utility::vector1<Size> uncovered_residues = find_uncovered_residues(*poses_[jpose], covered_resnum);
				if ( uncovered_residues.size() == 0 ) continue;

				core::pose::Pose inserted_pose(*poses_[jpose]);
				utility::vector1<int> remaining_resnum;
				remove_residues(inserted_pose, covered_resnum, remaining_resnum);
				for ( Size i=1; i<=remaining_resnum.size(); ++i ) resnum.push_back(remaining_resnum[i]);
				for ( Size i=1; i<=remaining_resnum.size(); ++i ) covered_resnum.push_back(remaining_resnum[i]);

				if ( jpose == 1 ) {
					full_length_pose = core::pose::PoseOP( new core::pose::Pose(inserted_pose) );
					core::pose::PDBInfoOP pdb_info;
					full_length_pose->pdb_info(pdb_info);
				} else {
					Size njump = full_length_pose->fold_tree().num_jump();
					full_length_pose->conformation().insert_conformation_by_jump( inserted_pose.conformation(),
						full_length_pose->total_residue() + 1, njump+1,
						full_length_pose->total_residue() );
				}
				//full_length_pose->dump_pdb("full_length_"+I(1,ipose)+"_"+I(1,jpose)+".pdb");
			}

			Size nres_domain1 = full_length_pose->total_residue();

			for ( Size ires=1; ires<=poses_[ipose]->total_residue(); ++ires ) resnum.push_back(poses_[ipose]->pdb_info()->number(ires));

			Size jump_num = full_length_pose->fold_tree().num_jump()+1;
			full_length_pose->conformation().insert_conformation_by_jump( poses_[ipose]->conformation(),
				full_length_pose->total_residue() + 1, jump_num,
				full_length_pose->total_residue() );
			//TR << full_length_pose->fold_tree() << std::endl;
			//full_length_pose->dump_pdb("full_length_"+I(1,ipose)+"before_docking.pdb");

			for ( Size ires=1; ires <= nres_domain1; ++ires ) {
				for ( Size jres=nres_domain1+1; jres <= full_length_pose->total_residue(); ++jres ) {
					Size seq_sep = std::abs(resnum[ires] - resnum[jres]);
					if ( seq_sep>=6 && seq_sep <=8 ) {
						core::Real gd = gap_distance(seq_sep);
						Size iatom = full_length_pose->residue_type(ires).atom_index("CA");
						Size jatom = full_length_pose->residue_type(jres).atom_index("CA");

						TR.Debug << "Adding constraints to residue " << ires << " and " << jres << std::endl;
						core::scoring::func::FuncOP fx( new core::scoring::constraints::BoundFunc( 0, gd, 5., "gap" ) );
						full_length_pose->add_constraint( core::scoring::constraints::ConstraintCOP( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AtomPairConstraint(
							core::id::AtomID(iatom, ires),
							core::id::AtomID(jatom, jres),
							fx
							) ) ) );
					}
				}
			}

			/*
			if ( poses_[1]->pdb_info()->number(1) < poses_[ipose]->pdb_info()->number(1) ) {
			full_length_pose = new core::pose::Pose(*poses_[1]);

			core::pose::PDBInfoOP pdb_info;
			full_length_pose->pdb_info(pdb_info);

			// remove overlap residues
			first_domain_end = full_length_pose->total_residue();
			while (poses_[1]->pdb_info()->number(first_domain_end) >= range_start) {
			--first_domain_end;
			}
			if (first_domain_end < full_length_pose->total_residue()) {
			full_length_pose->conformation().delete_residue_range_slow(first_domain_end+1, full_length_pose->total_residue());
			}

			full_length_pose->conformation().insert_conformation_by_jump( poses_[ipose]->conformation(),
			full_length_pose->total_residue() + 1,
			full_length_pose->total_residue(),
			1);

			utility::vector1 < utility::vector1 <Size > > covered_range;
			for (Size jpose = 1; jpose < ipose; ++jpose) {
			Size resi_start(0);
			Size resi_end(0);
			for (Size jres = 1; jres <= poses_[jpose]->total_residue(); ++jres) {
			if ( poses_[jpose]->pdb_info()->number(jres) > range_end ) {
			if (resi_start == 0) {
			resi_start = jres;
			resi_end = jres;
			}
			else {
			if ( poses_[jpose]->pdb_info()->number(jres) < poses_[jpose]->pdb_info()->number(resi_start) ) {
			// add residue
			resi_start = jres;
			}
			if ( poses_[jpose]->pdb_info()->number(jres) > poses_[jpose]->pdb_info()->number(resi_end) ) {
			// add residue
			resi_end = jres;
			}
			}
			}
			}
			if (resi_start != 0) {
			for (Size irange=1; irange <= covered_range.size(); ++irange) {

			}

			core::pose::Pose inserted_pose(*poses_[jpose], resi_start, resi_end);
			full_length_pose->conformation().insert_conformation_by_jump( poses_[ipose]->conformation(),
			full_length_pose->total_residue() + 1,
			full_length_pose->total_residue(),
			1);
			}

			}
			else {
			full_length_pose = new core::pose::Pose(*poses_[ipose]);

			core::pose::PDBInfoOP pdb_info;
			full_length_pose->pdb_info(pdb_info);

			first_domain_end = full_length_pose->total_residue();
			full_length_pose->conformation().insert_conformation_by_jump( poses_[1]->conformation(),
			full_length_pose->total_residue() + 1,
			full_length_pose->total_residue(),
			1);
			}
			*/

			// docking
			/*
			Size iatom = full_length_pose->residue_type(first_domain_end-2).atom_index("CA");
			Size jatom = full_length_pose->residue_type(first_domain_end+3).atom_index("CA");
			TR << "Adding constraints to residue " << first_domain_end-2 << " and " << first_domain_end+3 << std::endl;
			full_length_pose->add_constraint(
			new core::scoring::constraints::AtomPairConstraint(
			core::id::AtomID(iatom, first_domain_end-2),
			core::id::AtomID(jatom, first_domain_end+3),
			new core::scoring::constraints::BoundFunc( 0, 24., 5., "gap" ) )
			);
			*/

			using namespace basic::options;
			using namespace basic::options::OptionKeys;
			scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function(
				"score4_smooth", "" );
			scorefxn_->set_weight( core::scoring::atom_pair_constraint, 1.0 );

			// randomize orientation
			protocols::rigid::RigidBodyRandomizeMoverOP rb_mover( new protocols::rigid::RigidBodyRandomizeMover(*full_length_pose, jump_num) );
			rb_mover->apply(*full_length_pose);

			// put the two domains apart
			core::scoring::ScoreFunctionOP simple_scorefxn = core::scoring::ScoreFunctionFactory::create_score_function(
				"score0", "" );

			protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( *full_length_pose, jump_num) );
			translate->step_size( 1000. );
			core::Vector translation_axis(1.0,0.0,0.0);
			translate->trans_axis(translation_axis);
			translate->apply( *full_length_pose );

			core::Real unbound_score = (*simple_scorefxn)(*full_length_pose);

			translation_axis = core::Vector(-1.0,0.0,0.0);
			translate->trans_axis(translation_axis);
			translate->apply( *full_length_pose );

			core::Real max_step_size(5.0);
			core::Size counter(0);
			while ( (*simple_scorefxn)(*full_length_pose) > unbound_score+10. ) {
				//TR << "Random translation " << (*simple_scorefxn)(*full_length_pose) << std::endl;
				translate->step_size( max_step_size*numeric::random::rg().uniform() );
				translation_axis = core::Vector(numeric::random::rg().uniform(),numeric::random::rg().uniform(),numeric::random::rg().uniform());
				translate->trans_axis(translation_axis);
				translate->apply( *full_length_pose );
				counter++;
				if ( counter > 100 ) {
					max_step_size *= 1.5;
					counter = 0;
				}
			}
			//full_length_pose->dump_pdb("full_length_"+I(1,ipose)+"before_docking2.pdb");

			using namespace protocols::docking;
			DockingLowResOP docking_mover( new DockingLowRes(scorefxn_, jump_num) );
			docking_mover->set_outer_cycles(20);
			docking_mover->apply(*full_length_pose);
			//scorefxn_->show(*full_length_pose);
			//full_length_pose->dump_pdb("test_full_length.pdb");

			// align full_length_pose to pose1_
			std::list <Size> residue_list1;
			std::list <Size> residue_list2;
			for ( Size ires = 1; ires <= poses_[1]->total_residue(); ++ires ) {
				residue_list1.push_back(ires);
				residue_list2.push_back(ires);
			}
			TMalign_poses(*full_length_pose, *poses_[1], residue_list1, residue_list2);
			//full_length_pose->dump_pdb("full_length_after_align.pdb");

			// align ipose to full_length_pose
			residue_list1.clear();
			residue_list2.clear();
			for ( Size ires = 1; ires <= poses_[ipose]->total_residue(); ++ires ) {
				residue_list1.push_back(ires);
				residue_list2.push_back(nres_domain1+ires);
			}

			TMalign_poses(*poses_[ipose], *full_length_pose, residue_list1, residue_list2);

			//full_length_pose->dump_pdb("full_length_"+I(1,ipose)+"after_docking.pdb");
			//poses_[ipose]->dump_pdb("aligned_pose.pdb");
		}
	}
}

void DomainAssembly::apply(
	core::pose::Pose & pose
)
{
	using namespace protocols::docking;

	//core::pose::PoseOP lower_pose(pose1_->pdb_info()->number(1) < pose2_->pdb_info()->number(1) ? pose1_ : pose2_);
	//core::pose::PoseOP higher_pose(pose1_->pdb_info()->number(1) < pose2_->pdb_info()->number(1) ? pose2_ : pose1_);
	core::pose::PoseOP lower_pose(pose1_);
	core::pose::PoseOP higher_pose(pose2_);

	int jump_pos1 = numeric::random::rg().random_range(1, lower_pose->total_residue());
	int jump_pos2 = numeric::random::rg().random_range(lower_pose->total_residue()+1, lower_pose->total_residue()+higher_pose->total_residue());

	core::kinematics::FoldTree ft = pose.fold_tree();
	int jump_num = ft.new_jump(jump_pos1, jump_pos2, (int)lower_pose->total_residue()+1);
	pose.fold_tree( ft );

	pose.copy_segment(lower_pose->total_residue(), *lower_pose, 1, 1);
	pose.copy_segment(higher_pose->total_residue(), *higher_pose, lower_pose->total_residue()+1, 1);


	Size gap_start;
	if ( lower_pose->total_residue() > 3 ) {
		gap_start = lower_pose->total_residue() - 3;
	} else {
		gap_start = 1;
	}

	Size gap_stop;
	if ( higher_pose->total_residue() > 3 ) {
		gap_stop = 3;
	} else {
		gap_stop = higher_pose->total_residue();
	}
	//int seq_sep = higher_pose->pdb_info()->number(gap_stop) - higher_pose->pdb_info()->number(gap_start);

	//if (seq_sep <= 8) {
	if ( !lower_pose->residue_type(gap_start).is_protein() ) utility_exit_with_message("Error! not an amino acid!");
	if ( !higher_pose->residue_type(gap_stop ).is_protein() ) utility_exit_with_message("Error! not an amino acid!");
	Size iatom = pose.residue_type(gap_start).atom_index("CA");
	Size jatom = pose.residue_type(lower_pose->total_residue()+gap_stop ).atom_index("CA");
	TR << "Adding constraints to residue " << gap_start << " and " << lower_pose->total_residue()+gap_stop << std::endl;
	core::scoring::func::FuncOP fx( new core::scoring::constraints::BoundFunc( 0, 24., 5., "gap" ) );
	pose.add_constraint( core::scoring::constraints::ConstraintCOP( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AtomPairConstraint(
		core::id::AtomID(iatom, gap_start),
		core::id::AtomID(jatom, lower_pose->total_residue()+gap_stop),
		fx
		) ) ) );
	//}

	DockingLowResOP docking_mover( new DockingLowRes(scorefxn_, jump_num) );
	docking_mover->apply(pose);
}

} // hybridize
//} // comparative_modeling
} // protocols
