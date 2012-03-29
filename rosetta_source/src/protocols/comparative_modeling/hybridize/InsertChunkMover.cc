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
#include <core/conformation/util.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/util/kinematics_util.hh>

//#include <core/kinematics/FoldTree.hh>

// history
#include <protocols/comparative_modeling/hybridize/TemplateHistory.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

//
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

void InsertChunkMover::set_template(core::pose::PoseCOP template_pose, core::Size template_id,
				  std::map <core::Size, core::Size> const & sequence_alignment ) {
	template_pose_ = template_pose;
	template_id_ = template_id;
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


void InsertChunkMover::steal_torsion_and_bonds_from_template(core::pose::Pose & pose) {
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

			//fpd  bondlengths and angles
			core::conformation::Residue const &res = template_pose_->residue(jres_template);
			for ( unsigned int iatom = 1; iatom <= res.natoms(); ++iatom ) {
				/*
				std::string atom_name(pose.residue_type(ires_pose).atom_name(iatom));
                Size jatom;
				if (template_pose_->residue_type(jres_template).has(atom_name)) {
                    jatom = template_pose_->residue_type(jres_template).atom_index(atom_name);
				}
				else {
					continue;
				}
				 */
				core::id::AtomID atm_ij     ( iatom, ires_pose);
				core::id::AtomID atm_ijtempl( iatom, jres_template );
				core::id::DOF_ID d_ij		( atm_ij, core::id::D );
				core::id::DOF_ID theta_ij   ( atm_ij, core::id::THETA );
				core::id::DOF_ID phi_ij     ( atm_ij, core::id::PHI );
				
				if ( pose.has_dof(d_ij) ) {
					//core::id::AtomID stub_atom1 = pose.conformation().atom_tree().get_stub_atom1_id(atm_ij);
					//core::id::AtomID stub_atom2 = pose.conformation().atom_tree().get_stub_atom2_id(atm_ij);
					//core::id::AtomID stub_atom3 = pose.conformation().atom_tree().get_stub_atom3_id(atm_ij);
					//core::id::AtomID stub2 
				}
				
				core::id::DOF_ID theta_ijtempl(  atm_ijtempl, core::id::THETA );
				core::id::DOF_ID d_ijtempl(  atm_ijtempl, core::id::D );

				if (discontinued_upper(*template_pose_,jres_template) && iatom == 1) continue; // don't steal bond distance/angle across jump
				if ( template_pose_->has_dof(d_ijtempl) && pose.has_dof(d_ij) ) {
					pose.set_dof( d_ij, template_pose_->dof(d_ijtempl) );
				}

				if (discontinued_upper(*template_pose_,jres_template) && iatom == 2) continue; // don't steal bond distance/angle across jump
				if (discontinued_lower(*template_pose_,jres_template) && iatom == 3) continue; // don't steal bond distance/angle across jump
				if ( template_pose_->has_dof(theta_ijtempl) && pose.has_dof(theta_ij) ) {
					pose.set_dof( theta_ij, template_pose_->dof(theta_ijtempl) );
				TR.Debug << "template dof: " << I(4,jres_template) << " " << template_pose_->residue(jres_template).name3() << " " 
				<< template_pose_->residue_type(jres_template).atom_name(iatom) << F(8,3,template_pose_->dof(theta_ijtempl)) << F(8,3,template_pose_->dof(d_ijtempl)) << std::endl;
				TR.Debug << "target dof:   " << I(4,ires_pose) << " " << pose.residue(ires_pose).name3() << " "
				<< pose.residue_type(ires_pose).atom_name(iatom) << F(8,3,pose.dof(theta_ij)) << F(8,3,pose.dof(d_ij)) << std::endl;
				}

				//fpd make sure we're not inserting improper bonds
				if (template_pose_->has_dof(d_ijtempl) && template_pose_->dof(d_ijtempl) > 2.5 && iatom != 6) {
					TR << "WARNING! Inserting long bond "
					   << template_pose_->dof(d_ijtempl) << " at atom " << iatom << " res " << ires_pose << " from templ res " << jres_template << std::endl;
				}
			}

			while (ires_pose > align_trial_counter_.size()) {
				align_trial_counter_.push_back(0);
			}
			++align_trial_counter_[ires_pose];
		}
		//TR.Debug << "torsion: " << I(4,ires_pose) << F(8,3, pose.phi(ires_pose)) << F(8,3, pose.psi(ires_pose)) << std::endl;
	}
    TR.Debug << "End " << std::endl;
}


bool InsertChunkMover::get_local_sequence_mapping(core::pose::Pose & pose,
								int registry_shift,
								Size MAX_TRIAL)
{
	core::Size counter = 0;
	TR.Debug << sequence_alignment_ << std::endl;
	while (counter < MAX_TRIAL) {
		++counter;
		sequence_alignment_local_.clear();
		core::pose::initialize_atomid_map( atom_map_, pose, core::id::BOGUS_ATOM_ID );
		
		core::Size seqpos_pose = RG.random_range(seqpos_start_, seqpos_stop_);
		TR.Debug << "Align Seqpos: " << seqpos_pose << std::endl;
		if (sequence_alignment_.find(seqpos_pose+registry_shift) == sequence_alignment_.end()) continue;
		core::Size seqpos_template = sequence_alignment_.find(seqpos_pose+registry_shift)->second;

        //TR << "local_align: " << seqpos_pose << " " << seqpos_template << " " << registry_shift << " " << seqpos_start_ << " " << seqpos_stop_ << std::endl;

		if(align_to_ss_only_) {
			if (template_pose_->secstruct(seqpos_template) == 'L') continue;
		}
		if (template_pose_->secstruct(seqpos_template) != 'L') {
			secstruct_ = template_pose_->secstruct(seqpos_template);
		}

        // collect local alignment for stealing torsion
		seqpos_aligned_start_ = seqpos_pose;
		seqpos_aligned_stop_ = seqpos_pose;
		for (Size ires_pose=seqpos_pose; ires_pose>=seqpos_start_; --ires_pose) {
    		if (sequence_alignment_.find(ires_pose+registry_shift) == sequence_alignment_.end()) break;
    		core::Size jres_template = sequence_alignment_.find(ires_pose+registry_shift)->second;
			if ( jres_template <= 0 || jres_template > template_pose_->total_residue() ) continue;
            
			if ( !template_pose_->residue_type(jres_template).is_protein() ) continue;
			if(copy_ss_torsion_only_) {
				if (template_pose_->secstruct(jres_template) == 'L') continue;
			}
			
			sequence_alignment_local_[ires_pose] = jres_template;
			seqpos_aligned_start_ = ires_pose;
			
			if (discontinued_upper(*template_pose_,jres_template)) {
				TR.Debug << "Disconnect upper: " << ires_pose << " "  << jres_template << std::endl;
				break;
			}
		}
		for (Size ires_pose=seqpos_pose+1; ires_pose<=seqpos_stop_; ++ires_pose) {
    		if (sequence_alignment_.find(ires_pose+registry_shift) == sequence_alignment_.end()) break;
    		core::Size jres_template = sequence_alignment_.find(ires_pose+registry_shift)->second;
			if ( jres_template <= 0 || jres_template > template_pose_->total_residue() ) continue;
			
			if ( !template_pose_->residue_type(jres_template).is_protein() ) continue;
			if(copy_ss_torsion_only_) {
				if (template_pose_->secstruct(jres_template) == 'L') continue;
			}
			
			sequence_alignment_local_[ires_pose] = jres_template;
			seqpos_aligned_stop_ = ires_pose;
			
			if (discontinued_lower(*template_pose_,jres_template)) {
				TR.Debug << "Disconnect lower: " << ires_pose << " "  << jres_template << std::endl;
				break;
			}
		}
        
        // collect atom_map for superposition
		core::Size atom_map_count = 0;
		for (Size ires_pose=seqpos_pose; ires_pose>=seqpos_start_; --ires_pose) {
    		if (sequence_alignment_.find(ires_pose+registry_shift) == sequence_alignment_.end()) break;
    		core::Size jres_template = sequence_alignment_.find(ires_pose+registry_shift)->second;
			if ( jres_template <= 0 || jres_template > template_pose_->total_residue() ) continue;

			if ( !template_pose_->residue_type(jres_template).is_protein() ) continue;
			if(copy_ss_torsion_only_) {
				if (template_pose_->secstruct(jres_template) == 'L') continue;
			}
			
            core::id::AtomID const id1( pose.residue_type(ires_pose).atom_index("CA"), ires_pose );
			core::id::AtomID const id2( template_pose_->residue_type(jres_template).atom_index("CA"), jres_template );
			atom_map_[ id1 ] = id2;
			++atom_map_count;

			if (discontinued_upper(*template_pose_,jres_template)) break;
		}
		for (Size ires_pose=seqpos_pose+1; ires_pose<=seqpos_stop_; ++ires_pose) {
    		if (sequence_alignment_.find(ires_pose+registry_shift) == sequence_alignment_.end()) break;
    		core::Size jres_template = sequence_alignment_.find(ires_pose+registry_shift)->second;
			if ( jres_template <= 0 || jres_template > template_pose_->total_residue() ) continue;
			
			if ( !template_pose_->residue_type(jres_template).is_protein() ) continue;
			if(copy_ss_torsion_only_) {
				if (template_pose_->secstruct(jres_template) == 'L') continue;
			}
			
			core::id::AtomID const id1( pose.residue_type(ires_pose).atom_index("CA"), ires_pose );
			core::id::AtomID const id2( template_pose_->residue_type(jres_template).atom_index("CA"), jres_template );
			atom_map_[ id1 ] = id2;
			++atom_map_count;
			if (discontinued_lower(*template_pose_,jres_template)) break;
		}
		
		if (atom_map_count >=3) {
            TR.Debug << sequence_alignment_local_ << std::endl;
			using namespace basic::options;
			using namespace basic::options::OptionKeys;
			if ( option[cm::hybridize::move_anchor]() ) {
				change_anchor(pose, jump_number_, seqpos_pose);
			}
			return true;
		}
	}
    TR.Debug << "Failing to get aligned: " << sequence_alignment_local_ << std::endl;
    sequence_alignment_local_.clear();
	return false;	
}

void InsertChunkMover::change_anchor(core::pose::Pose & pose,
				   Size jump_number,
				   Size new_anchor_seqpos) {
	
	// Update the fold tree with the new jump points
	kinematics::FoldTree f ( pose.conformation().fold_tree() );
	TR.Debug << "starting tree:" << f << std::endl;
	
	// Setup the lists of jumps and cuts
	Size num_jumps( f.num_jump() );
	Size num_cuts( f.num_cutpoint() );
	
	utility::vector1< int > cuts_vector( f.cutpoints() );
	ObjexxFCL::FArray1D_int cuts( num_cuts );
	ObjexxFCL::FArray2D_int jumps( 2, num_jumps );
	
	// Initialize jumps
	for ( Size i = 1; i<= num_jumps; ++i ) {
		int down ( f.downstream_jump_residue(i) );
		int up ( f.upstream_jump_residue(i) );
		if ( down < up ) {
			jumps(1,i) = down;
			jumps(2,i) = up;
		} else {
			jumps(1,i) = up;
			jumps(2,i) = down;
		}
	}
	
	for ( Size i = 1; i<= num_cuts; ++i ) {
		cuts(i) = cuts_vector[i];
	}
	
	// This is the basejump
	int root ( f.root() );
	int residue_that_builds_anchor( f.upstream_jump_residue( jump_number ) );
	
	jumps(1, jump_number ) = new_anchor_seqpos;
	jumps(2, jump_number ) = residue_that_builds_anchor;
	
	f.tree_from_jumps_and_cuts( pose.conformation().size(), num_jumps, jumps, cuts );
	f.reorder( root );

	TR.Debug << "final tree:" << f << std::endl;

	pose.conformation().fold_tree( f );
	
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
	success_ = get_local_sequence_mapping(pose, registry_shift_, (seqpos_stop_-seqpos_start_+1));
	if (!success_) return;
	
	//steal_torsion_from_template(pose);
	//align_chunk(pose);
	set_bb_xyz_aligned(pose);
	//steal_torsion_and_bonds_from_template(pose);
	//align_chunk(pose);
	check_overlap(pose);
}
    
void InsertChunkMover::set_bb_xyz_aligned(core::pose::Pose & pose) {
    utility::vector1< core::id::AtomID > ids;
	utility::vector1< numeric::xyzVector<core::Real> > positions;
    utility::vector1< core::id::AtomID > sch_ids;
	utility::vector1< numeric::xyzVector<core::Real> > sch_positions;

    utility::vector1< core::Size > aligned_residues;
    utility::vector1< core::Size > non_aligned_residues;
    bool aligned;
    Size jump_residue_pose = pose.fold_tree().downstream_jump_residue(jump_number_);
	TR.Debug << "Jump residue: " << jump_residue_pose << std::endl;

	/*
    if (sequence_alignment_local_.find(jump_residue_pose) != sequence_alignment_local_.end()) {
        aligned = true;
    }
    else {
        return;
    }
    
    bool all_aligned_lower = true;
	Size non_aligned_lower_start(seqpos_stop_+1);
    for (Size ires_pose=jump_residue_pose; ires_pose<=seqpos_stop_; ++ires_pose) {
		if (sequence_alignment_local_.find(ires_pose) == sequence_alignment_local_.end()) {
            all_aligned_lower = false;
            non_aligned_lower_start = ires_pose;
            break;
        }
    }
    bool all_aligned_upper = true;
    Size non_aligned_upper_end(seqpos_start_-1);
    for (Size ires_pose=jump_residue_pose; ires_pose>=seqpos_start_; --ires_pose) {
		if (sequence_alignment_local_.find(ires_pose) == sequence_alignment_local_.end()) {
            all_aligned_upper = false;
            non_aligned_upper_end = ires_pose;
            break;
        }
    }
    */
	
	// copy xyz of the backbone
    for (Size ires_pose=seqpos_aligned_start_; ires_pose<=seqpos_aligned_stop_; ++ires_pose) {
		if (sequence_alignment_local_.find(ires_pose) != sequence_alignment_local_.end()) {
			core::Size jres_template = sequence_alignment_local_.find(ires_pose)->second;
			TR.Debug << "Copy xyz of residue " << ires_pose << std::endl;
            for ( Size iatom=1; iatom <= pose.residue_type(ires_pose).last_backbone_atom(); ++iatom ) { // use residue_type to prevent internal coord update
                std::string atom_name(pose.residue_type(ires_pose).atom_name(iatom));
                if (template_pose_->residue_type(jres_template).has(atom_name)) {
                    Size jatom = template_pose_->residue_type(jres_template).atom_index(atom_name);
                    ids.push_back(core::id::AtomID(iatom,ires_pose));
                    positions.push_back(template_pose_->xyz(core::id::AtomID(jatom,jres_template)));
                }
                else {
                    sch_ids.push_back(core::id::AtomID(iatom,ires_pose));
                }
            }

            for ( Size iatom=pose.residue_type(ires_pose).last_backbone_atom()+1;
                 iatom<= pose.residue_type(ires_pose).natoms(); ++iatom ) { // use residue_type to prevent internal coord update
                sch_ids.push_back(core::id::AtomID(iatom,ires_pose));
            }
            
			while (ires_pose > align_trial_counter_.size()) {
				align_trial_counter_.push_back(0);
			}
			++align_trial_counter_[ires_pose];
		}
	}
	pose.batch_set_xyz(ids,positions);

	// idealize sidechains
    for (Size iatom = 1; iatom <= sch_ids.size(); ++iatom) {
        sch_positions.push_back(
                                pose.residue(sch_ids[iatom].rsd()).build_atom_ideal(
                                                                                    sch_ids[iatom].atomno(),
                                                                                    pose.conformation()
                                                                                    )
                                );
    }
    pose.batch_set_xyz(sch_ids,sch_positions);
    
	// idealize the connection between copied and uncopied region
    if (seqpos_aligned_start_ > seqpos_start_) {
        core::conformation::idealize_position(seqpos_aligned_start_, pose.conformation());
        core::conformation::idealize_position(seqpos_aligned_start_-1, pose.conformation());
    }
    if (seqpos_aligned_stop_ < seqpos_stop_) {
        core::conformation::idealize_position(seqpos_aligned_stop_, pose.conformation());
        core::conformation::idealize_position(seqpos_aligned_stop_+1, pose.conformation());
    }
	
	runtime_assert( pose.data().has( core::pose::datacache::CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY ) );
	TemplateHistory &history = 
	*( static_cast< TemplateHistory* >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY )() ));
	history.set( seqpos_start_, seqpos_stop_, template_id_ );
}
	
void InsertChunkMover::check_overlap(core::pose::Pose & pose) {
	bool overlapped = false;
	for ( Size ires=seqpos_start_; ires<= seqpos_stop_; ++ires ) {
		if (!pose.residue_type(ires).has("CA")) continue;
		for ( Size jres=1; jres<= pose.total_residue(); ++jres ) {
			if (jres >=seqpos_start_ && jres<= seqpos_stop_) continue;
			if (!pose.residue_type(jres).has("CA")) continue;
			numeric::xyzVector < core::Real > xyz_iatom (pose.residue(ires).xyz("CA"));
			numeric::xyzVector < core::Real > xyz_jatom (pose.residue(jres).xyz("CA"));
			if (xyz_iatom.distance_squared(xyz_jatom) < 1e-4) {
				overlapped = true;
				break;
			}
		}
		if (overlapped) break;
	}
	if (overlapped) {
		utility::vector1< core::id::AtomID > ids;
		utility::vector1< numeric::xyzVector<core::Real> > positions;
		numeric::xyzVector<core::Real> trans(2.*RG.uniform()-1.,
											 2.*RG.uniform()-1.,
											 2.*RG.uniform()-1.);
		
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
	
	
} // hybridize 
} // comparative_modeling 
} // protocols
