// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Align a random jump to template
/// @details
/// @author Yifan Song

#ifndef apps_pilot_yfsong_AlignChunkMover_HH
#define apps_pilot_yfsong_AlignChunkMover_HH

#include <apps/pilot/yfsong/AlignChunkMover.fwd.hh>

namespace challenge {
basic::options::StringOptionKey  ss("challenge:ss");
basic::options::BooleanOptionKey    aligned("challenge:aligned");
basic::options::IntegerVectorOptionKey    chunk_mapping("challenge:chunk_mapping");
}

namespace apps {
namespace pilot {

using namespace core;
using namespace core::kinematics;
using namespace ObjexxFCL;
using namespace protocols::moves;
using namespace protocols::loops;
using namespace protocols::nonlocal;
using namespace numeric::model_quality;
using namespace id;
using namespace basic::options;
using namespace basic::options::OptionKeys;

bool discontinued_upper(core::pose::Pose const & pose, Size const seqpos) {
	core::Real N_C_cutoff(2.0);

	if ( seqpos == 1 ) return true;
	if ( !pose.residue_type(seqpos).is_polymer() ) return true;
	if ( !pose.residue_type(seqpos-1).is_polymer() ) return true;
	if ( pose.residue_type(seqpos).is_protein() && pose.residue_type(seqpos-1).is_protein() ) {
		if ( pose.residue(seqpos).xyz("N").distance(pose.residue(seqpos-1).xyz("C")) > N_C_cutoff ) {
			return true;
		}
	}
	return false;
}

bool discontinued_lower(core::pose::Pose const & pose, Size const seqpos) {
	core::Real N_C_cutoff(2.0);

	if ( seqpos == pose.size() ) return true;
	if ( !pose.residue_type(seqpos).is_polymer() ) return true;
	if ( !pose.residue_type(seqpos+1).is_polymer() ) return true;
	if ( pose.residue_type(seqpos).is_protein() && pose.residue_type(seqpos+1).is_protein() ) {
		if ( pose.residue(seqpos).xyz("C").distance(pose.residue(seqpos+1).xyz("N")) > N_C_cutoff ) {
			return true;
		}
	}
	return false;
}

enum AlignOption { all_chunks, random_chunk };

std::list < Size >
downstream_residues_from_jump(core::pose::Pose const & pose, Size const jump_number) {
	std::list < Size > residue_list;
	utility::vector1< Edge > edges = pose.fold_tree().get_outgoing_edges(pose.fold_tree().jump_edge(jump_number).stop());

	for ( Size i_edge = 1; i_edge <= edges.size(); ++i_edge ) {
		if ( !edges[i_edge].is_polymer() ) continue;
		Size start = edges[i_edge].start() <= edges[i_edge].stop() ? edges[i_edge].start() : edges[i_edge].stop();
		Size stop  = edges[i_edge].start() <= edges[i_edge].stop() ? edges[i_edge].stop()  : edges[i_edge].start();
		for ( Size ires = start; ires <= stop; ++ires ) {
			residue_list.push_back(ires);
		}

	}
	residue_list.sort();
	residue_list.unique();
	return residue_list;
}

class AlignChunkMover: public protocols::moves::Mover
{
public:

	AlignChunkMover(numeric::random::RandomGenerator & RG) :
		RG_(RG), registry_shift_(0), reset_torsion_unaligned_(true), align_to_ss_only_(true), copy_ss_torsion_only_(false), secstruct_('L')
	{
		align_trial_counter_.clear();
	}

	// atom_map: from mod_pose to ref_pose
	void
	get_superposition_transformation(
		pose::Pose const & mod_pose,
		pose::Pose const & ref_pose,
		id::AtomID_Map< id::AtomID > const & atom_map,
		numeric::xyzMatrix< core::Real > &R,
		numeric::xyzVector< core::Real > &preT,
		numeric::xyzVector< core::Real > &postT
	)
	{
		// count number of atoms for the array
		Size total_mapped_atoms(0);
		for ( Size ires=1; ires<= mod_pose.size(); ++ires ) {
			for ( Size iatom=1; iatom<= mod_pose.residue(ires).natoms(); ++iatom ) {
				AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
				if ( !aid.valid() ) continue;

				++total_mapped_atoms;
			}
		}

		preT = postT = numeric::xyzVector< core::Real >(0,0,0);
		if ( total_mapped_atoms <= 2 ) {
			R.xx() = R.yy() = R.zz() = 1;
			R.xy() = R.yx() = R.zx() = R.zy() = R.yz() = R.xz() = 0;
			return;
		}

		ObjexxFCL::FArray2D< core::Real > final_coords( 3, total_mapped_atoms );
		ObjexxFCL::FArray2D< core::Real > init_coords( 3, total_mapped_atoms );
		preT = postT = numeric::xyzVector< core::Real >(0,0,0);
		Size atomno(0);
		for ( Size ires=1; ires<= mod_pose.size(); ++ires ) {
			for ( Size iatom=1; iatom<= mod_pose.residue(ires).natoms(); ++iatom ) {
				AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
				if ( !aid.valid() ) continue;
				++atomno;

				numeric::xyzVector< core::Real > x_i = mod_pose.residue(ires).atom(iatom).xyz();
				preT += x_i;
				numeric::xyzVector< core::Real > y_i = ref_pose.xyz( aid );
				postT += y_i;

				for ( int j=0; j<3; ++j ) {
					init_coords(j+1,atomno) = x_i[j];
					final_coords(j+1,atomno) = y_i[j];
				}
			}
		}

		preT /= (float) total_mapped_atoms;
		postT /= (float) total_mapped_atoms;
		for ( int i=1; i<=(int)total_mapped_atoms; ++i ) {
			for ( int j=0; j<3; ++j ) {
				init_coords(j+1,i) -= preT[j];
				final_coords(j+1,i) -= postT[j];
			}
		}

		// get optimal superposition
		// rotate >init< to >final<
		ObjexxFCL::FArray1D< numeric::Real > ww( total_mapped_atoms, 1.0 );
		ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
		numeric::Real ctx;

		numeric::model_quality::findUU( init_coords, final_coords, ww, total_mapped_atoms, uu, ctx );
		R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
		R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
		R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );
	}

	void
	apply_transform(
		pose::Pose & mod_pose,
		std::list <Size> const & residue_list,
		numeric::xyzMatrix< core::Real > const & R, numeric::xyzVector< core::Real > const & preT, numeric::xyzVector< core::Real > const & postT
	)
	{ // translate xx2 by COM and fill in the new ref_pose coordinates
		utility::vector1< core::id::AtomID > ids;
		utility::vector1< numeric::xyzVector<core::Real> > positions;

		Vector x2;
		FArray2D_double xx2;
		FArray1D_double COM(3);
		for ( std::list<Size>::const_iterator it = residue_list.begin();
				it != residue_list.end();
				++it ) {
			Size ires = *it;
			for ( Size iatom=1; iatom<= mod_pose.residue_type(ires).natoms(); ++iatom ) { // use residue_type to prevent internal coord update
				ids.push_back(core::id::AtomID(iatom,ires));
				positions.push_back(postT + (R*( mod_pose.xyz(core::id::AtomID(iatom,ires)) - preT )));
			}
		}
		mod_pose.batch_set_xyz(ids,positions);
	}

	void align_chunk(core::pose::Pose & pose) {
		std::list <Size> residue_list;
		for ( Size ires_pose=seqpos_start_; ires_pose<=seqpos_stop_; ++ires_pose ) {
			residue_list.push_back(ires_pose);
		}

		numeric::xyzMatrix< core::Real > R;
		numeric::xyzVector< core::Real > preT;
		numeric::xyzVector< core::Real > postT;

		get_superposition_transformation( pose, *template_pose_, atom_map_, R, preT, postT );
		apply_transform( pose, residue_list, R, preT, postT );
	}

	void set_template(core::pose::PoseCOP template_pose,
		std::map <core::Size, core::Size> const & sequence_alignment ) {
		template_pose_ = template_pose;
		sequence_alignment_ = sequence_alignment;
	}

	void set_aligned_chunk(core::pose::Pose const & pose, Size const jump_number) {
		jump_number_ = jump_number;

		std::list < Size > downstream_residues = downstream_residues_from_jump(pose, jump_number_);
		seqpos_start_ = downstream_residues.front();
		seqpos_stop_ = downstream_residues.back();

		// make sure it is continuous, may not be necessary if the function gets expanded to handle more than 1 chunk
		assert(downstream_residues.size() == (seqpos_stop_ - seqpos_start_ + 1));
	}

	void set_reset_torsion_unaligned(bool reset_torsion_unaligned) {
		reset_torsion_unaligned_ = reset_torsion_unaligned;
	}

	void steal_torsion_from_template(core::pose::Pose & pose) {
		basic::Tracer TR( "pilot.yfsong.util" );
		using namespace ObjexxFCL::format;
		for ( Size ires_pose=seqpos_start_; ires_pose<=seqpos_stop_; ++ires_pose ) {
			if ( reset_torsion_unaligned_ ) {
				pose.set_omega(ires_pose, 180);

				if ( secstruct_ == 'H' ) {
					pose.set_phi(ires_pose,  -60);
					pose.set_psi(ires_pose,  -45);
				} else {
					pose.set_phi(ires_pose,  -110);
					pose.set_psi(ires_pose,   130);
				}
			}

			if ( sequence_alignment_local_.find(ires_pose) != sequence_alignment_local_.end() ) {
				core::Size jres_template = sequence_alignment_local_.find(ires_pose)->second;
				if ( !discontinued_upper(*template_pose_,jres_template) ) {
					TR.Debug << "template phi: " << I(4,jres_template) << F(8,2, template_pose_->phi(jres_template)) << std::endl;
					pose.set_phi(ires_pose, template_pose_->phi(jres_template));
				}
				if ( !discontinued_lower(*template_pose_,jres_template) ) {
					TR.Debug << "template psi: " << I(4,jres_template) << F(8,2, template_pose_->psi(jres_template)) << std::endl;
					pose.set_psi(ires_pose, template_pose_->psi(jres_template));
				}
				pose.set_omega(ires_pose, template_pose_->omega(jres_template));

				while ( ires_pose > align_trial_counter_.size() ) {
					align_trial_counter_.push_back(0);
				}
				++align_trial_counter_[ires_pose];
			}
			TR.Debug << "torsion: " << I(4,ires_pose) << F(8,3, pose.phi(ires_pose)) << F(8,3, pose.psi(ires_pose)) << std::endl;
		}
	}

	bool get_local_sequence_mapping(core::pose::Pose const & pose,
		int registry_shift = 0,
		Size MAX_TRIAL = 10 )
	{
		core::Size counter = 0;
		while ( counter < MAX_TRIAL ) {
			++counter;
			sequence_alignment_local_.clear();
			core::pose::initialize_atomid_map( atom_map_, pose, core::id::AtomID::BOGUS_ATOM_ID() );

			core::Size seqpos_pose = RG_.random_range(seqpos_start_, seqpos_stop_);

			if ( sequence_alignment_.find(seqpos_pose+registry_shift) == sequence_alignment_.end() ) continue;
			core::Size seqpos_template = sequence_alignment_.find(seqpos_pose+registry_shift)->second;

			if ( align_to_ss_only_ ) {
				if ( template_pose_->secstruct(seqpos_template) == 'L' ) continue;
			}
			if ( template_pose_->secstruct(seqpos_template) != 'L' ) {
				secstruct_ = template_pose_->secstruct(seqpos_template);
			}


			core::Size atom_map_count = 0;
			for ( Size ires_pose=seqpos_pose; ires_pose>=seqpos_start_; --ires_pose ) {
				int jres_template = ires_pose + seqpos_template - seqpos_pose;
				if ( discontinued_upper(*template_pose_,jres_template) ) break;

				if ( jres_template <= 0 || jres_template > template_pose_->total_residue() ) continue;
				if ( !template_pose_->residue_type(jres_template).is_protein() ) continue;
				if ( copy_ss_torsion_only_ ) {
					if ( template_pose_->secstruct(jres_template) == 'L' ) continue;
				}

				sequence_alignment_local_[ires_pose] = jres_template;
				core::id::AtomID const id1( pose.residue_type(ires_pose).atom_index("CA"), ires_pose );
				core::id::AtomID const id2( template_pose_->residue_type(jres_template).atom_index("CA"), jres_template );
				atom_map_[ id1 ] = id2;
				++atom_map_count;
			}
			for ( Size ires_pose=seqpos_pose+1; ires_pose<=seqpos_stop_; ++ires_pose ) {
				int jres_template = ires_pose + seqpos_template - seqpos_pose;
				if ( discontinued_lower(*template_pose_,jres_template) ) break;

				if ( jres_template <= 0 || jres_template > template_pose_->total_residue() ) continue;
				if ( !template_pose_->residue_type(jres_template).is_protein() ) continue;
				if ( copy_ss_torsion_only_ ) {
					if ( template_pose_->secstruct(jres_template) == 'L' ) continue;
				}

				sequence_alignment_local_[ires_pose] = jres_template;
				core::id::AtomID const id1( pose.residue_type(ires_pose).atom_index("CA"), ires_pose );
				core::id::AtomID const id2( template_pose_->residue_type(jres_template).atom_index("CA"), jres_template );
				atom_map_[ id1 ] = id2;
				++atom_map_count;
			}

			if ( atom_map_count >=3 ) {
				return true;
			}
		}
		return false;
	}

	void set_registry_shift(int registry_shift) {
		registry_shift_ = registry_shift;
	}

	Size trial_counter(Size ires) {
		if ( ires <= align_trial_counter_.size() ) {
			return align_trial_counter_[ires];
		}
		return 0;
	}

	void
	apply(core::pose::Pose & pose) {
		// apply alignment
		bool success = get_local_sequence_mapping(pose, registry_shift_);
		if ( !success ) return;

		steal_torsion_from_template(pose);
		align_chunk(pose);
	}

	std::string
	get_name() const {
		return "AlignChunkMover";
	}

private:
	numeric::random::RandomGenerator & RG_;
	core::pose::PoseCOP template_pose_;
	std::map <core::Size, core::Size> sequence_alignment_;

	Size jump_number_; // the jump to be realigned
	Size seqpos_start_; // start and end seqpose of the chunk
	Size seqpos_stop_;
	int registry_shift_;
	char secstruct_;

	// parameters of the protocol
	bool reset_torsion_unaligned_; // reset torsion of unaligned region to the default value for the given secondary structure
	bool align_to_ss_only_; // only use the secondary structure portion to align to the template
	bool copy_ss_torsion_only_; // only copy the secondary structure information from the template

	std::map <core::Size, core::Size> sequence_alignment_local_;
	core::id::AtomID_Map< core::id::AtomID > atom_map_; // atom map for superposition
	utility::vector1 <Size> align_trial_counter_;
};

class MultiTemplateAlignChunkMover: public protocols::moves::Mover
{

public:
	MultiTemplateAlignChunkMover(numeric::random::RandomGenerator & RG,
		utility::vector1 < core::pose::PoseOP > const & template_poses,
		Loops ss_chunks_pose,
		AlignOption align_option = all_chunks,
		Size max_registry_shift = 0) :
		RG_(RG),
		template_poses_(template_poses),
		align_option_(align_option),
		align_chunk_(RG_),
		max_registry_shift_input_(max_registry_shift)
	{
		bool alignment_from_template = option[challenge::aligned]();

		// set up secstruct chunks
		template_ss_chunks_.clear();
		template_ss_chunks_.resize(template_poses_.size());
		for ( core::Size i_template=1; i_template<=template_poses_.size(); ++i_template ) {
			// find ss chunks in template
			template_ss_chunks_[i_template] = extract_secondary_structure_chunks(*template_poses_[i_template]);
		}

		Size count = 0;
		for ( core::Size i_template=1; i_template<=template_poses_.size(); ++i_template ) {
			if ( template_ss_chunks_[i_template].size() != 0 ) ++count;
		}
		if ( count == 0 ) {
			utility_exit_with_message("Template structures need at least one secondary structure for this protocol");
		}

		sequence_alignments_.clear();
		for ( core::Size i_template=1; i_template<=template_poses_.size(); ++i_template ) {
			std::map <core::Size, core::Size> sequence_alignment;
			if ( alignment_from_template ) {
				get_alignment_from_template(template_poses_[i_template], sequence_alignment);
			} else {
				std::map <core::Size, core::Size> chunk_mapping;
				if ( option[challenge::chunk_mapping].user() ) {
					for ( Size i=1; i<=option[challenge::chunk_mapping]().size(); ++i ) {
						chunk_mapping[i] = option[challenge::chunk_mapping]()[i];
					}
				}

				get_alignment_from_chunk_mapping(chunk_mapping, template_ss_chunks_[i_template], ss_chunks_pose, sequence_alignment);
			}
			sequence_alignments_.push_back(sequence_alignment);
		}
	}

	void
	get_alignment_from_template(core::pose::PoseCOP const template_pose, std::map <core::Size, core::Size> & seqpos_alignment) {
		// specific to this case, alignment comes from residue number
		for ( core::Size ires=1; ires<=template_pose->total_residue(); ++ires ) {
			seqpos_alignment[template_pose->pdb_info()->number(ires)] = ires;
		}
	}


	void
	get_alignment_from_chunk_mapping(std::map <core::Size, core::Size> const & chunk_mapping,
		Loops const template_ss_chunks,
		Loops const target_ss_chunks,
		std::map <core::Size, core::Size> & sequence_alignment)
	{
		max_registry_shift_.resize(target_ss_chunks.size());
		for ( Size i_chunk_pose = 1; i_chunk_pose <= target_ss_chunks.size(); ++i_chunk_pose ) {
			if ( chunk_mapping.find(i_chunk_pose) == chunk_mapping.end() ) continue;
			core::Size j_chunk_template = chunk_mapping.find(i_chunk_pose)->second;

			Size respos_mid_pose = (target_ss_chunks[i_chunk_pose].start() + target_ss_chunks[i_chunk_pose].stop()) / 2;
			Size respos_mid_template = (template_ss_chunks[j_chunk_template].start() + template_ss_chunks[j_chunk_template].stop()) / 2;
			int offset = respos_mid_template - respos_mid_pose;

			using namespace ObjexxFCL::format;
			if ( target_ss_chunks[i_chunk_pose].length() <= template_ss_chunks[j_chunk_template].length() ) {
				max_registry_shift_[i_chunk_pose] = max_registry_shift_input_ + template_ss_chunks[j_chunk_template].length() - target_ss_chunks[i_chunk_pose].length();
				for ( Size ires=target_ss_chunks[i_chunk_pose].start(); ires<=target_ss_chunks[i_chunk_pose].stop(); ++ires ) {
					sequence_alignment[ires] = ires+offset;
					//std::cout << I(4, ires) << I(4, ires+offset) << std::endl;
				}
			} else {
				max_registry_shift_[i_chunk_pose] = max_registry_shift_input_ + target_ss_chunks[i_chunk_pose].length() - template_ss_chunks[j_chunk_template].length();
				for ( Size ires_templ=template_ss_chunks[j_chunk_template].start(); ires_templ<=template_ss_chunks[j_chunk_template].stop(); ++ires_templ ) {
					sequence_alignment[ires_templ-offset] = ires_templ;
					//std::cout << I(4, ires_templ) << I(4, ires_templ-offset) << std::endl;
				}
			}
		}
	}

	void pick_random_template() {
		assert(template_poses_.size() != 0);

		template_number_ = 0;
		while ( !template_number_ ) {
			template_number_ = RG_.random_range(1, template_poses_.size());
			if ( template_ss_chunks_[template_number_].size() == 0 ) template_number_ = 0;
		}
	}

	void pick_random_chunk(core::pose::Pose & pose) {
		jump_number_ = RG_.random_range(1, pose.num_jump());
	}

	Size trial_counter(Size ires) {
		return align_chunk_.trial_counter(ires);
	}

	void
	apply(core::pose::Pose & pose) {
		max_registry_shift_.resize(pose.num_jump(), max_registry_shift_input_);

		// pick a random template
		pick_random_template();
		align_chunk_.set_template(template_poses_[template_number_], sequence_alignments_[template_number_]);

		// random chunk or loop all chunks
		if ( align_option_ == random_chunk ) {
			// pick a random jump
			pick_random_chunk(pose);
			align_chunk_.set_aligned_chunk(pose, jump_number_);
			align_chunk_.set_reset_torsion_unaligned(false);

			// apply alignment
			int registry_shift = RG_.random_range(-max_registry_shift_[jump_number_], max_registry_shift_[jump_number_]);
			align_chunk_.set_registry_shift(registry_shift);
			align_chunk_.apply(pose);
		} else {
			// loop over all jumps
			for ( core::Size jump_number=1; jump_number<=pose.num_jump(); ++jump_number ) {
				align_chunk_.set_aligned_chunk(pose, jump_number);

				// apply alignment
				int registry_shift = RG_.random_range(-max_registry_shift_[jump_number], max_registry_shift_[jump_number]);
				align_chunk_.set_registry_shift(registry_shift);
				align_chunk_.apply(pose);
			}
		}
	}

	std::string
	get_name() const {
		return "MultiTemplateAlignChunkMover";
	}

private:
	AlignChunkMover align_chunk_;
	numeric::random::RandomGenerator & RG_;
	AlignOption align_option_;

	utility::vector1 < core::pose::PoseCOP > template_poses_;
	utility::vector1 < Loops > template_ss_chunks_;
	utility::vector1 < std::map <core::Size, core::Size> > sequence_alignments_;
	Size max_registry_shift_input_;
	utility::vector1 < Size > max_registry_shift_;

	Size template_number_; // the jump to be realigned
	Size jump_number_; // the jump to be realigned
}; //class MultiTemplateAlignChunkMover

class CustomStarTreeMover: public protocols::moves::Mover
{
	void
	apply(core::pose::Pose & pose)
	{
		basic::Tracer TR( "pilot.yfsong.util" );
		if ( ! option[ OptionKeys::in::file::psipred_ss2 ].user() ) {
			utility_exit_with_message("Error in reading psipred_ss2 file, is the -in:file:psipred_ss2 flag set correctly?");
		}
		bool check_psipred = set_secstruct_from_psipred_ss2(pose);
		assert (check_psipred);

		// Build the star fold tree, identify jumps
		if ( option[challenge::ss].user() ) {
			ss_chunks_pose_ = extract_secondary_structure_chunks( pose, option[challenge::ss]() );
		} else {
			ss_chunks_pose_ = extract_secondary_structure_chunks( pose );
		}
		ss_chunks_pose_.sequential_order();
		TR.Debug << "Target secondary chunks:" << std::endl;
		TR.Debug << ss_chunks_pose_ << std::endl;
		loops_pose_ = ss_chunks_pose_.invert(pose.size());
		TR.Debug << "Target loops: " << pose.size() << std::endl;
		TR.Debug << loops_pose_ << std::endl;
		TR.flush();

		/*
		if (option[challenge::virtual_loop]()) {
		set_loops_to_virt_ala(pose, loops_pose_);
		}
		*/
		// complete the chunks to cover the whole protein and customize cutpoints
		Loops chunks(ss_chunks_pose_);
		for ( Size i=1; i<=chunks.num_loop(); ++i ) {
			if ( i==1 ) {
				chunks[i].set_start(1);
			} else {
				Size new_start = (ss_chunks_pose_[i-1].stop() + ss_chunks_pose_[i].start() + 1) / 2;
				//Size new_start = (ss_chunks_pose_[i-1].stop() + 1);
				chunks[i].set_start( new_start );
			}

			if ( i==chunks.num_loop() ) {
				chunks[i].set_stop(pose.size());
			} else {
				Size new_stop = (ss_chunks_pose_[i].stop() + ss_chunks_pose_[i+1].start() - 1) / 2;
				chunks[i].set_stop( new_stop );
			}
		}
		TR.Debug << "Chunks: " << pose.size() << std::endl;
		TR.Debug << chunks << std::endl;

		StarTreeBuilder builder;
		if ( chunks.num_loop() > 0 ) {
			builder.set_up(chunks, &pose);
		}
		TR.Debug << pose.fold_tree() << std::endl;
		protocols::nonlocal::add_cutpoint_variants(&pose);
	}

	std::string
	get_name() const {
		return "CustomStarTreeMover";
	}
private:
	Loops ss_chunks_pose_;
	Loops loops_pose_;

}; // class CustomStarTreeMover

/// @brief A Mover that insert fragments to the loop region
class CustomFragmentMover: public protocols::moves::Mover
{
public:
	CustomFragmentMover(numeric::random::RandomGenerator & RG,
		Loops const & ss_chunks) : RG_(RG)
	{
		protocols::loops::read_loop_fragments(frag_libs_);

		frag_insert_pos_weights_.resize(frag_libs_.size());
		for ( Size i_frag_set = 1; i_frag_set<=frag_libs_.size(); ++i_frag_set ) {
			for ( Size i_frame = 1; i_frame <= frag_libs_[i_frag_set]->nr_frames(); ++i_frame ) {
				core::fragment::FrameIterator frame_it = frag_libs_[i_frag_set]->begin();
				advance(frame_it, i_frame-1);
				core::Size seqpos_start = (*frame_it)->start();
				core::Size seqpos_end   = (*frame_it)->end();
				for ( Size seqpos = seqpos_start; seqpos <= seqpos_end; ++seqpos ) {
					if ( ! ss_chunks.has(seqpos) ) {
						frag_insert_pos_weights_[i_frag_set].push_back(i_frame);
					}
				}
			}
		}
		// code shamelessly stolen from nonlocal/SingleFragmentMover.cc
		/*
		for (Size i_frag_set = 1; i_frag_set<=frag_libs.size(); ++i_frag_set) {
		for (core::fragment::FrameIterator i = frag_libs[i_frag_set]->begin(); i != frag_libs[i_frag_set]->end(); ++i) {
		core::Size position = (*i)->start();
		core::fragment::Frame frame = **i;
		//library_[position] = **i;
		}
		}
		*/
	}
	void
	apply(core::pose::Pose & pose) {
		Size i_frag_set = RG_.random_range(1, frag_libs_.size());

		// pick insertion position
		//Size insert_pos = RG_.random_range(1, frag_libs_[i_frag_set]->nr_frames());
		Size insert_pos = frag_insert_pos_weights_[i_frag_set][RG_.random_range(1, frag_insert_pos_weights_[i_frag_set].size())];
		core::fragment::FrameIterator frame_it = frag_libs_[i_frag_set]->begin();
		advance(frame_it, insert_pos-1);
		Size i_frag = RG_.random_range(1, frame_it->nr_frags());

		frame_it->apply( i_frag, pose );
	}

	std::string
	get_name() const {
		return "CustomFragmentMover";
	}
private:
	numeric::random::RandomGenerator & RG_;
	utility::vector1< core::fragment::FragSetOP > frag_libs_;
	utility::vector1< utility::vector1< Size > > frag_insert_pos_weights_;
};

class HelixMover: public protocols::moves::Mover
{
public:
	HelixMover(numeric::random::RandomGenerator & RG) : RG_(RG)
	{
	}

	numeric::xyzVector<core::Real>
	get_helix_center(core::pose::Pose const & pose, Size const ires) {
		assert(ires > 1 && ires < pose.size());
		for ( Size i = ires-1; i<=ires+1; ++i ) {
			assert(pose.residue_type(i).is_protein());
		}

		numeric::xyzVector<core::Real> center_vect = pose.residue(ires-1).xyz("CA") + pose.residue(ires+1).xyz("CA") - 2.*pose.residue(ires).xyz("CA");
		center_vect.normalize();
		numeric::xyzVector<core::Real> helix_center = pose.residue(ires).xyz("CA") + center_vect * 2.3;
		return helix_center;

		//using namespace ObjexxFCL::format;
		//std::cout << "ATOM   9999  X   DUM A 999    " << F(8,3,helix_center.x())<< F(8,3,helix_center.y())<< F(8,3,helix_center.z()) << std::endl;
	}

	bool
	get_helix_centers(core::pose::Pose const & pose, Size const seqpos_start, Size const seqpos_stop,
		utility::vector1< numeric::xyzVector<core::Real> > & helix_centers)
	{
		assert(helix_centers.size() == pose.size());
		utility::vector1< bool > updated(pose.size(), false);
		Size counts(0);
		for ( Size ires=seqpos_start+1; ires<seqpos_stop; ++ires ) {
			if ( discontinued_upper(pose, ires) || discontinued_lower(pose, ires) ) continue;
			if ( pose.secstruct(ires-1) == 'H' &&
					pose.secstruct(ires)   == 'H' &&
					pose.secstruct(ires+1) == 'H' ) {

				helix_centers[ires] = get_helix_center(pose, ires);
				updated[ires] = true;
				++counts;
			}
		}
		if ( counts < 4 ) return false;

		// extrapolate helix center from the edge of helices, does not follow the trace of protein anymore
		bool updating = true;
		while ( updating ) {
			updating = false;
			for ( Size ires=seqpos_start; ires<=seqpos_stop; ++ires ) {
				if ( updated[ires] ) continue;
				if ( ires-1>=1 && ires+1 <= pose.size() ) {
					if ( !discontinued_upper(pose, ires) && !discontinued_lower(pose, ires) ) { // check if there is chain breaks
						if ( updated[ires-1] && updated[ires+1] ) {
							helix_centers[ires] = (helix_centers[ires-1] + helix_centers[ires+1]) / 2.;
							updated[ires] = true;
							updating = true;
							continue;
						}
					}
				}
				if ( ires-2>=1 ) {
					if ( updated[ires-2] && updated[ires-1] ) {
						if ( !discontinued_upper(pose, ires) && !discontinued_upper(pose, ires-1) ) {
							helix_centers[ires] = 2. * helix_centers[ires-1] - helix_centers[ires-2];
							updated[ires] = true;
							updating = true;
							continue;
						}
					}
				}
				if ( ires+2 <= pose.size() ) {
					if ( updated[ires+2] && updated[ires+1] ) {
						if ( !discontinued_lower(pose, ires) && !discontinued_lower(pose, ires+1) ) {
							helix_centers[ires] = 2. * helix_centers[ires+1] - helix_centers[ires+2];
							updated[ires] = true;
							updating = true;
							continue;
						}
					}
				}
			}
		}
		return true;
	}

	bool
	get_helix_vectors(core::pose::Pose const & pose,
		Size const seqpos_start, Size const seqpos_stop,
		utility::vector1< numeric::xyzVector<core::Real> > & helix_vectors) {
		utility::vector1< numeric::xyzVector<core::Real> > helix_centers(pose.size());
		bool success = get_helix_centers(pose,seqpos_start,seqpos_stop,helix_centers);
		if ( !success ) return false;

		for ( Size ires=seqpos_start; ires<seqpos_stop; ++ires ) {
			helix_vectors[ires] = helix_centers[ires+1] - helix_centers[ires];
		}
		helix_vectors[seqpos_stop] = helix_vectors[seqpos_stop-1];

		//smooth
		//smoothed_helix_vectors.resize(helix_length);
		/*
		for (Size ires=seqpos_start+1; ires<seqpos_stop; ++ires) {
		if (ires == 1) {
		smoothed_helix_vectors[ires] = helix_vectors[ires] + helix_vectors[ires+1];
		}
		else if (ires == helix_vectors.size()) {
		smoothed_helix_vectors[ires] = helix_vectors[ires-1] + helix_vectors[ires];
		}
		else {
		smoothed_helix_vectors[ires] = helix_vectors[ires-1] + 2. * helix_vectors[ires] + helix_vectors[ires+1];
		}
		smoothed_helix_vectors[ires].normalize();
		}
		*/
		return true;
	}

	void
	slide_non_ideal_helix(core::pose::Pose & pose, Size const seqpos_start, Size const seqpos_stop, core::Real const distance)
	{
		utility::vector1< numeric::xyzVector<core::Real> > helix_vectors(pose.size());
		bool success = get_helix_vectors(pose, seqpos_start, seqpos_stop, helix_vectors);
		if ( !success ) return;

		if ( distance > 0 ) {
			for ( Size ires = seqpos_start; ires < seqpos_stop; ++ires ) {
				helix_vectors[ires] += helix_vectors[ires+1];
				helix_vectors[ires].normalize();
			}
		} else {
			for ( Size ires = seqpos_stop; ires > seqpos_start; ++ires ) {
				helix_vectors[ires] += helix_vectors[ires-1];
				helix_vectors[ires].normalize();
			}
		}

		// apply transformation to each residue in the chain
		utility::vector1< core::id::AtomID > ids;
		utility::vector1< numeric::xyzVector<core::Real> > positions;
		for ( Size ires = seqpos_start; ires <= seqpos_stop; ++ires ) {
			for ( Size iatom = 1; iatom <= pose.residue_type(ires).natoms(); ++iatom ) {
				numeric::xyzVector<core::Real> atom_xyz = pose.xyz( core::id::AtomID(iatom,ires) );
				ids.push_back(core::id::AtomID(iatom,ires));
				positions.push_back(atom_xyz + helix_vectors[ires]*distance);
			}
		}
		pose.batch_set_xyz(ids,positions);
	}

	void
	apply(core::pose::Pose & pose)
	{
		Size jump_number = RG_.random_range(1, pose.num_jump());
		std::list < Size > downstream_residues = downstream_residues_from_jump(pose, jump_number);
		Size seqpos_start = downstream_residues.front();
		Size seqpos_stop = downstream_residues.back();

		Real distance = 2. * RG_.uniform();
		slide_non_ideal_helix(pose, seqpos_start, seqpos_stop, distance);
	}
	std::string
	get_name() const {
		return "HelixMover";
	}

private:
	numeric::random::RandomGenerator & RG_;

}; //class HelixMover

} // pilot
} // apps

#endif
