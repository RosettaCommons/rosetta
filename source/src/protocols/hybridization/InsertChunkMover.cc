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
/// @details
/// @author Yifan Song

#include <protocols/hybridization/InsertChunkMover.hh>
#include <protocols/hybridization/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/util/kinematics_util.hh>

//#include <core/kinematics/FoldTree.hh>

// history
#include <protocols/hybridization/TemplateHistory.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>


#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

#include <protocols/moves/Mover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>
#include <utility/stream_util.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.hybridization.InsertChunkMover" );

using utility::operator <<;

namespace protocols {
namespace hybridization {

using namespace core;
using namespace id;
using namespace ObjexxFCL;

InsertChunkMover::InsertChunkMover() {
	moves::Mover::type( "InsertChunkMover" );
	init();
}

void
InsertChunkMover::init() {
	registry_shift_ = 0;
	anchor_insert_only_ = false;
	align_to_ss_only_ = false;
	copy_ss_torsion_only_ = false;
	secstruct_ = 'L';

	template_id_ = 0;

	jump_number_ = 1;
	seqpos_start_ = 0;
	seqpos_stop_ = 0;
	seqpos_aligned_start_ = 0;
	seqpos_aligned_stop_ = 0;

	success_ = false;

	align_trial_counter_.clear();
}

InsertChunkMover::~InsertChunkMover() {}

void InsertChunkMover::set_template(core::pose::PoseCOP template_pose, core::Size template_id,
	std::map <core::Size, core::Size> const & sequence_alignment) {
	template_pose_ = template_pose;
	template_id_ = template_id;
	sequence_alignment_ = sequence_alignment;
}

void InsertChunkMover::set_aligned_chunk(core::pose::Pose const & pose, Size const jump_number, bool anchor_insert_only_in) {
	jump_number_ = jump_number;
	anchor_insert_only_ = anchor_insert_only_in;

	std::list < Size > downstream_residues = downstream_residues_from_jump(pose, jump_number_);
	seqpos_start_ = downstream_residues.front();
	seqpos_stop_ = downstream_residues.back();

	TR.Debug << "Setting Chunk Insertion from residue " << seqpos_aligned_start_ << " to " << seqpos_stop_ << " for jump " << jump_number_ << std::endl;

	// make sure it is continuous, may not be necessary if the function gets expanded to handle more than 1 chunk
	assert(downstream_residues.size() == (seqpos_stop_ - seqpos_start_ + 1));
}

bool InsertChunkMover::get_local_sequence_mapping(core::pose::Pose & pose,
	int registry_shift,
	Size MAX_TRIAL) {
	core::Size counter = 0;
	//TR.Debug << sequence_alignment_ << std::endl;
	while ( counter < MAX_TRIAL ) {
		++counter;
		sequence_alignment_local_.clear();
		core::pose::initialize_atomid_map( atom_map_, pose, core::id::BOGUS_ATOM_ID );

		//fpd pick a random downstream residue and steal it's position from a template
		//fpd if anchor_insert_only_ is set, use the jump anchor position
		core::Size seqpos_pose = numeric::random::rg().random_range(seqpos_start_, seqpos_stop_);
		if ( anchor_insert_only_ ) {
			seqpos_pose = pose.fold_tree().downstream_jump_residue( jump_number_ );
		}

		//TR.Debug << "Align Seqpos: " << seqpos_pose << std::endl;
		if ( sequence_alignment_.find(seqpos_pose+registry_shift) == sequence_alignment_.end() ) continue;
		core::Size seqpos_template = sequence_alignment_.find(seqpos_pose+registry_shift)->second;
		//TR.Debug << "Found Seqpos: " << seqpos_pose+registry_shift << " -> " << seqpos_template << std::endl;

		if ( align_to_ss_only_ && template_pose_->secstruct(seqpos_template) == 'L' ) continue;

		//TR.Debug << "Passed SS" << std::endl;

		if ( template_pose_->secstruct(seqpos_template) != 'L' ) {
			secstruct_ = template_pose_->secstruct(seqpos_template);
		}

		// collect local alignment for stealing torsion
		seqpos_aligned_start_ = seqpos_pose;
		seqpos_aligned_stop_ = seqpos_pose;
		for ( Size ires_pose=seqpos_pose; ires_pose>=seqpos_start_; --ires_pose ) {
			if ( sequence_alignment_.find(ires_pose+registry_shift) == sequence_alignment_.end() ) {
				break;
			}

			core::Size jres_template = sequence_alignment_.find(ires_pose+registry_shift)->second;
			if ( jres_template <= 0 || jres_template > template_pose_->total_residue() ) continue;

			if ( !template_pose_->residue_type(jres_template).is_protein() ) continue;
			if ( copy_ss_torsion_only_ && template_pose_->secstruct(jres_template) == 'L' ) continue;

			sequence_alignment_local_[ires_pose] = jres_template;
			seqpos_aligned_start_ = ires_pose;

			if ( discontinued_upper(*template_pose_,jres_template) ) {
				//TR.Debug << "Disconnect upper: " << ires_pose << " "  << jres_template << std::endl;
				break;
			}
		}

		for ( Size ires_pose=seqpos_pose+1; ires_pose<=seqpos_stop_; ++ires_pose ) {
			if ( sequence_alignment_.find(ires_pose+registry_shift) == sequence_alignment_.end() ) break;
			core::Size jres_template = sequence_alignment_.find(ires_pose+registry_shift)->second;

			if ( jres_template <= 0 || jres_template > template_pose_->total_residue() ) continue;

			if ( !template_pose_->residue_type(jres_template).is_protein() ) continue;
			if ( copy_ss_torsion_only_ && template_pose_->secstruct(jres_template) == 'L' ) continue;

			sequence_alignment_local_[ires_pose] = jres_template;
			seqpos_aligned_stop_ = ires_pose;

			if ( discontinued_lower(*template_pose_,jres_template) ) {
				//TR.Debug << "Disconnect lower: " << ires_pose << " "  << jres_template << std::endl;
				break;
			}
		}

		// collect atom_map for superposition
		core::Size atom_map_count = 0;
		for ( Size ires_pose=seqpos_pose; ires_pose>=seqpos_start_; --ires_pose ) {
			if ( sequence_alignment_.find(ires_pose+registry_shift) == sequence_alignment_.end() ) break;
			if ( !pose.residue_type(ires_pose).is_protein() ) continue;

			core::Size jres_template = sequence_alignment_.find(ires_pose+registry_shift)->second;
			if ( jres_template <= 0 || jres_template > template_pose_->total_residue() ) continue;

			if ( !template_pose_->residue_type(jres_template).is_protein() ) continue;
			if ( copy_ss_torsion_only_ && template_pose_->secstruct(jres_template) == 'L' ) continue;

			core::id::AtomID const id1( pose.residue_type(ires_pose).atom_index("CA"), ires_pose );
			core::id::AtomID const id2( template_pose_->residue_type(jres_template).atom_index("CA"), jres_template );
			atom_map_[ id1 ] = id2;
			++atom_map_count;

			if ( discontinued_upper(*template_pose_,jres_template) ) break;
		}
		for ( Size ires_pose=seqpos_pose+1; ires_pose<=seqpos_stop_; ++ires_pose ) {
			if ( sequence_alignment_.find(ires_pose+registry_shift) == sequence_alignment_.end() ) break;
			if ( !pose.residue_type(ires_pose).is_protein() ) continue;

			core::Size jres_template = sequence_alignment_.find(ires_pose+registry_shift)->second;
			if ( jres_template <= 0 || jres_template > template_pose_->total_residue() ) continue;

			if ( !template_pose_->residue_type(jres_template).is_protein() ) continue;
			if ( copy_ss_torsion_only_ ) {
				if ( template_pose_->secstruct(jres_template) == 'L' ) continue;
			}

			core::id::AtomID const id1( pose.residue_type(ires_pose).atom_index("CA"), ires_pose );
			core::id::AtomID const id2( template_pose_->residue_type(jres_template).atom_index("CA"), jres_template );
			atom_map_[ id1 ] = id2;
			++atom_map_count;
			if ( discontinued_lower(*template_pose_,jres_template) ) break;
		}

		// fpd we need at least 3 residues aligned
		if ( atom_map_count >=3 ) {
			//TR.Debug << sequence_alignment_local_ << std::endl;
			return true;
		}
	}

	//TR << "Failing to get aligned: " << sequence_alignment_local_ << std::endl;
	sequence_alignment_local_.clear();
	return false;
}


void InsertChunkMover::set_registry_shift(int registry_shift) {
	registry_shift_ = registry_shift;
}

Size InsertChunkMover::trial_counter(Size ires) {
	if ( ires <= align_trial_counter_.size() ) {
		return align_trial_counter_[ires];
	}
	return 0;
}

void
InsertChunkMover::apply(core::pose::Pose & pose) {
	// apply alignment
	success_ = get_local_sequence_mapping(pose, registry_shift_, (seqpos_stop_-seqpos_start_+1));
	if ( !success_ ) return;

	set_bb_xyz_aligned(pose);
	check_overlap(pose);
}

void InsertChunkMover::set_bb_xyz_aligned(core::pose::Pose & pose) {
	utility::vector1< core::id::AtomID > ids;
	utility::vector1< numeric::xyzVector<core::Real> > positions;
	utility::vector1< core::id::AtomID > sch_ids;
	utility::vector1< numeric::xyzVector<core::Real> > sch_positions;

	Size jump_residue_pose = pose.fold_tree().downstream_jump_residue(jump_number_);
	TR.Debug << "Jump residue: " << jump_residue_pose << std::endl;

	// copy xyz of the backbone
	for ( Size ires_pose=seqpos_aligned_start_; ires_pose<=seqpos_aligned_stop_; ++ires_pose ) {
		if ( sequence_alignment_local_.find(ires_pose) != sequence_alignment_local_.end() ) {
			core::Size jres_template = sequence_alignment_local_.find(ires_pose)->second;
			TR.Debug << "Copy xyz of residue " << ires_pose << std::endl;
			for ( Size iatom=1; iatom <= pose.residue_type(ires_pose).last_backbone_atom(); ++iatom ) { // use residue_type to prevent internal coord update
				std::string atom_name(pose.residue_type(ires_pose).atom_name(iatom));
				if ( template_pose_->residue_type(jres_template).has(atom_name) ) {
					Size jatom = template_pose_->residue_type(jres_template).atom_index(atom_name);
					ids.push_back(core::id::AtomID(iatom,ires_pose));
					positions.push_back(template_pose_->xyz(core::id::AtomID(jatom,jres_template)));
				} else {
					sch_ids.push_back(core::id::AtomID(iatom,ires_pose));
				}
			}

			for ( Size iatom=pose.residue_type(ires_pose).last_backbone_atom()+1;
					iatom<= pose.residue_type(ires_pose).natoms(); ++iatom ) { // use residue_type to prevent internal coord update
				sch_ids.push_back(core::id::AtomID(iatom,ires_pose));
			}

			while ( ires_pose > align_trial_counter_.size() ) {
				align_trial_counter_.push_back(0);
			}
			++align_trial_counter_[ires_pose];
		}
	}

	pose.batch_set_xyz(ids,positions);

	// idealize sidechains
	for ( Size iatom = 1; iatom <= sch_ids.size(); ++iatom ) {
		sch_positions.push_back(
			pose.residue(sch_ids[iatom].rsd()).build_atom_ideal( sch_ids[iatom].atomno(), pose.conformation() )
		);
	}
	pose.batch_set_xyz(sch_ids,sch_positions);

	// idealize the connection between copied and uncopied region
	if ( seqpos_aligned_start_ > seqpos_start_ ) {
		core::conformation::idealize_position(seqpos_aligned_start_, pose.conformation());
		core::conformation::idealize_position(seqpos_aligned_start_-1, pose.conformation());
	}
	if ( seqpos_aligned_stop_ < seqpos_stop_ ) {
		core::conformation::idealize_position(seqpos_aligned_stop_, pose.conformation());
		core::conformation::idealize_position(seqpos_aligned_stop_+1, pose.conformation());
	}

	runtime_assert( pose.data().has( core::pose::datacache::CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY ) );
	TemplateHistory &history =
		*( utility::pointer::static_pointer_cast< protocols::hybridization::TemplateHistory > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY ) ));
	history.set( seqpos_start_, seqpos_stop_, template_id_ );
}

void InsertChunkMover::check_overlap(core::pose::Pose & pose) {
	bool overlapped = false;
	for ( Size ires=seqpos_start_; ires<= seqpos_stop_; ++ires ) {
		if ( !pose.residue_type(ires).has("CA") ) continue;
		for ( Size jres=1; jres<= pose.total_residue(); ++jres ) {
			if ( jres >=seqpos_start_ && jres<= seqpos_stop_ ) continue;
			if ( !pose.residue_type(jres).has("CA") ) continue;
			numeric::xyzVector < core::Real > xyz_iatom (pose.residue(ires).xyz("CA"));
			numeric::xyzVector < core::Real > xyz_jatom (pose.residue(jres).xyz("CA"));
			if ( xyz_iatom.distance_squared(xyz_jatom) < 1e-4 ) {
				overlapped = true;
				break;
			}
		}
		if ( overlapped ) break;
	}
	if ( overlapped ) {
		utility::vector1< core::id::AtomID > ids;
		utility::vector1< numeric::xyzVector<core::Real> > positions;
		numeric::xyzVector<core::Real> trans(
			2.*numeric::random::rg().uniform()-1.,
			2.*numeric::random::rg().uniform()-1.,
			2.*numeric::random::rg().uniform()-1.);

		for ( Size ires=seqpos_start_; ires<= seqpos_stop_; ++ires ) {
			for ( Size iatom=1; iatom<= pose.residue_type(ires).natoms(); ++iatom ) { // use residue_type to prevent internal coord update
				ids.push_back(core::id::AtomID(iatom,ires));
				positions.push_back( pose.residue(ires).xyz(iatom) + trans);
			}
		}
		pose.batch_set_xyz(ids,positions);
	}
}

std::string
InsertChunkMover::get_name() const {
	return "InsertChunkMover";
}


} // hybridization
} // protocols
