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

#include <protocols/comparative_modeling/hybridize/InsertChunkMover.hh>
#include <protocols/comparative_modeling/hybridize/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/util/kinematics_util.hh>

#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>

#include <numeric/xyz.functions.hh>

#include <protocols/moves/Mover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>

#include <basic/Tracer.hh>

static numeric::random::RandomGenerator RG(1183103);
static basic::Tracer TR( "protocols.comparative_modeling.hybridize.InsertChunkMover" );

namespace protocols {
namespace comparative_modeling {
namespace hybridize {

using namespace core;
using namespace id;
using namespace ObjexxFCL;

InsertChunkMover::InsertChunkMover() : 
registry_shift_(0), reset_torsion_unaligned_(false), align_to_ss_only_(false), copy_ss_torsion_only_(false), secstruct_('L')
{
	align_trial_counter_.clear();
}

InsertChunkMover::~InsertChunkMover(){}

// atom_map: from mod_pose to ref_pose
void
InsertChunkMover::get_superposition_transformation(
                                                   pose::Pose const & mod_pose,
                                                   pose::Pose const & ref_pose,
                                                   id::AtomID_Map< id::AtomID > const & atom_map,
                                                   numeric::xyzMatrix< core::Real > &R, numeric::xyzVector< core::Real > &preT, numeric::xyzVector< core::Real > &postT )
{
	// count number of atoms for the array
	Size total_mapped_atoms(0);
	for ( Size ires=1; ires<= mod_pose.total_residue(); ++ires ) {
		for ( Size iatom=1; iatom<= mod_pose.residue(ires).natoms(); ++iatom ) {
			AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
			if (!aid.valid()) continue;
			
			++total_mapped_atoms;
		}
	}
	
	preT = postT = numeric::xyzVector< core::Real >(0,0,0);
	if (total_mapped_atoms <= 2) {
		R.xx() = R.yy() = R.zz() = 1;
		R.xy() = R.yx() = R.zx() = R.zy() = R.yz() = R.xz() = 0;
		return;
	}
	
	ObjexxFCL::FArray2D< core::Real > final_coords( 3, total_mapped_atoms );
	ObjexxFCL::FArray2D< core::Real > init_coords( 3, total_mapped_atoms );
	preT = postT = numeric::xyzVector< core::Real >(0,0,0);
	Size atomno(0);
	for ( Size ires=1; ires<= mod_pose.total_residue(); ++ires ) {
		for ( Size iatom=1; iatom<= mod_pose.residue(ires).natoms(); ++iatom ) {
			AtomID const & aid( atom_map[ id::AtomID( iatom,ires ) ] );
			if (!aid.valid()) continue;
			++atomno;
			
			numeric::xyzVector< core::Real > x_i = mod_pose.residue(ires).atom(iatom).xyz();
			preT += x_i;
			numeric::xyzVector< core::Real > y_i = ref_pose.xyz( aid );
			postT += y_i;
			
			for (int j=0; j<3; ++j) {
				init_coords(j+1,atomno) = x_i[j];
				final_coords(j+1,atomno) = y_i[j];
			}
		}
	}
	
	preT /= (float) total_mapped_atoms;
	postT /= (float) total_mapped_atoms;
	for (int i=1; i<=(int)total_mapped_atoms; ++i) {
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
InsertChunkMover::apply_transform(
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
	for (std::list<Size>::const_iterator it = residue_list.begin();
		 it != residue_list.end();
		 ++it) {
		Size ires = *it;
		for ( Size iatom=1; iatom<= mod_pose.residue_type(ires).natoms(); ++iatom ) { // use residue_type to prevent internal coord update
			ids.push_back(core::id::AtomID(iatom,ires));
			positions.push_back(postT + (R*( mod_pose.xyz(core::id::AtomID(iatom,ires)) - preT )));
		}
	}
	mod_pose.batch_set_xyz(ids,positions);
}
	
void InsertChunkMover::align_chunk(core::pose::Pose & pose) {
	std::list <Size> residue_list;
	for (Size ires_pose=seqpos_start_; ires_pose<=seqpos_stop_; ++ires_pose) {
		residue_list.push_back(ires_pose);
	}
	
	numeric::xyzMatrix< core::Real > R;
	numeric::xyzVector< core::Real > preT;
	numeric::xyzVector< core::Real > postT;
	
	get_superposition_transformation( pose, *template_pose_, atom_map_, R, preT, postT );
	apply_transform( pose, residue_list, R, preT, postT );
}
	
void InsertChunkMover::set_template(core::pose::PoseCOP template_pose,
				  std::map <core::Size, core::Size> const & sequence_alignment ) {
	template_pose_ = template_pose;
	sequence_alignment_ = sequence_alignment;
}

void InsertChunkMover::set_aligned_chunk(core::pose::Pose const & pose, Size const jump_number) {
	jump_number_ = jump_number;
	
	std::list < Size > downstream_residues = downstream_residues_from_jump(pose, jump_number_);
	seqpos_start_ = downstream_residues.front();
	seqpos_stop_ = downstream_residues.back();
	
	// make sure it is continuous, may not be necessary if the function gets expanded to handle more than 1 chunk
	assert(downstream_residues.size() == (seqpos_stop_ - seqpos_start_ + 1));
}

void InsertChunkMover::set_reset_torsion_unaligned(bool reset_torsion_unaligned) {
	reset_torsion_unaligned_ = reset_torsion_unaligned;
}
	
void InsertChunkMover::steal_torsion_from_template(core::pose::Pose & pose) {
	using namespace ObjexxFCL::fmt;
	for (Size ires_pose=seqpos_start_; ires_pose<=seqpos_stop_; ++ires_pose) {
		if (reset_torsion_unaligned_) {
			pose.set_omega(ires_pose, 180);
			
			if (secstruct_ == 'H') {
				pose.set_phi(ires_pose,	 -60);
				pose.set_psi(ires_pose,	 -45);
			}
			else {
				pose.set_phi(ires_pose,	 -110);
				pose.set_psi(ires_pose,	  130);
			}
		}
		
		if (sequence_alignment_local_.find(ires_pose) != sequence_alignment_local_.end()) {
			core::Size jres_template = sequence_alignment_local_.find(ires_pose)->second;
			if ( !discontinued_upper(*template_pose_,jres_template) ) {
				TR.Debug << "template phi: " << I(4,jres_template) << F(8,2, template_pose_->phi(jres_template)) << std::endl;
				pose.set_phi(ires_pose,	template_pose_->phi(jres_template));
			}
			if ( !discontinued_lower(*template_pose_,jres_template) ) {
				TR.Debug << "template psi: " << I(4,jres_template) << F(8,2, template_pose_->psi(jres_template)) << std::endl;
				pose.set_psi(ires_pose,	template_pose_->psi(jres_template));
			}
			pose.set_omega(ires_pose,	template_pose_->omega(jres_template));
			
			while (ires_pose > align_trial_counter_.size()) {
				align_trial_counter_.push_back(0);
			}
			++align_trial_counter_[ires_pose];
		}
		TR.Debug << "torsion: " << I(4,ires_pose) << F(8,3, pose.phi(ires_pose)) << F(8,3, pose.psi(ires_pose)) << std::endl;
	}
}
	
bool InsertChunkMover::get_local_sequence_mapping(core::pose::Pose const & pose,
								int registry_shift,
								Size MAX_TRIAL)
{
	core::Size counter = 0;
	while (counter < MAX_TRIAL) {
		++counter;
		sequence_alignment_local_.clear();
		core::pose::initialize_atomid_map( atom_map_, pose, core::id::BOGUS_ATOM_ID );
		
		core::Size seqpos_pose = RG.random_range(seqpos_start_, seqpos_stop_);

		if (sequence_alignment_.find(seqpos_pose+registry_shift) == sequence_alignment_.end()) continue;
		core::Size seqpos_template = sequence_alignment_.find(seqpos_pose+registry_shift)->second;
		
		if(align_to_ss_only_) {
			if (template_pose_->secstruct(seqpos_template) == 'L') continue;
		}
		if (template_pose_->secstruct(seqpos_template) != 'L') {
			secstruct_ = template_pose_->secstruct(seqpos_template);
		}

		
		core::Size atom_map_count = 0;
		for (Size ires_pose=seqpos_pose; ires_pose>=seqpos_start_; --ires_pose) {
			int jres_template = ires_pose + seqpos_template - seqpos_pose;
			if ( jres_template <= 0 || jres_template > template_pose_->total_residue() ) continue;
			if (discontinued_upper(*template_pose_,jres_template)) break;

			if ( !template_pose_->residue_type(jres_template).is_protein() ) continue;
			if(copy_ss_torsion_only_) {
				if (template_pose_->secstruct(jres_template) == 'L') continue;
			}
			
			sequence_alignment_local_[ires_pose] = jres_template;
			core::id::AtomID const id1( pose.residue_type(ires_pose).atom_index("CA"), ires_pose );
			core::id::AtomID const id2( template_pose_->residue_type(jres_template).atom_index("CA"), jres_template );
			atom_map_[ id1 ] = id2;
			++atom_map_count;
		}
		for (Size ires_pose=seqpos_pose+1; ires_pose<=seqpos_stop_; ++ires_pose) {
			int jres_template = ires_pose + seqpos_template - seqpos_pose;
			if ( jres_template <= 0 || jres_template > template_pose_->total_residue() ) continue;
			if (discontinued_lower(*template_pose_,jres_template)) break;
			
			if ( !template_pose_->residue_type(jres_template).is_protein() ) continue;
			if(copy_ss_torsion_only_) {
				if (template_pose_->secstruct(jres_template) == 'L') continue;
			}
			
			sequence_alignment_local_[ires_pose] = jres_template;
			core::id::AtomID const id1( pose.residue_type(ires_pose).atom_index("CA"), ires_pose );
			core::id::AtomID const id2( template_pose_->residue_type(jres_template).atom_index("CA"), jres_template );
			atom_map_[ id1 ] = id2;
			++atom_map_count;
		}
		
		if (atom_map_count >=3) {
			return true;
		}
	}
	return false;	
}
	
void InsertChunkMover::set_registry_shift(int registry_shift) {
	registry_shift_ = registry_shift;
}
	
Size InsertChunkMover::trial_counter(Size ires) {
	if (ires <= align_trial_counter_.size()) {
		return align_trial_counter_[ires];	
	}
	return 0;
}
	
void
InsertChunkMover::apply(core::pose::Pose & pose) {
	// apply alignment
	success_ = get_local_sequence_mapping(pose, registry_shift_);
	if (!success_) return;
	
	steal_torsion_from_template(pose);
	align_chunk(pose);
}
	
std::string
InsertChunkMover::get_name() const {
	return "InsertChunkMover";
}
	
	
} // hybridize 
} // comparative_modeling 
} // protocols
