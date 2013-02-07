// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file LoophashLoopInserter.cc
///
/// @brief Use loophash to find a fragment from loop_anchor to loop_anchor+1. Build ideal residues
/// using the loophash sequence. Don't idealize the last loop connections (this results in a broken loop)
/// @author Tim Jacobs

//Unit
#include <protocols/loophash/LoophashLoopInserter.hh>
#include <protocols/loophash/LoophashLoopInserterCreator.hh>

//Numeric
#include <numeric/random/random.hh>
#include <numeric/geometry/hashing/SixDHasher.fwd.hh>

//Basic
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>

//Core
#include <core/pose/Pose.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

//protocols
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loops/Loop.hh>

//utility
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

namespace protocols {
namespace loophash {

static basic::Tracer TR( "protocols.loophash.LoophashLoopInserter" );

LoophashLoopInserterCreator::LoophashLoopInserterCreator() {}
LoophashLoopInserterCreator::~LoophashLoopInserterCreator() {}
protocols::loops::loop_creation::LoopInserterOP LoophashLoopInserterCreator::create_loop_inserter() const {
	return new LoophashLoopInserter;
}

std::string LoophashLoopInserterCreator::inserter_name() const {
	return "LoophashLoopInserter";
}
	
LoophashLoopInserter::LoophashLoopInserter():
num_flanking_residues_to_match_(3),
modify_flanking_regions_(false)
{
	using namespace protocols::loophash;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	
//	min_torsion_rms_=option[lh::min_bbrms];
//	max_torsion_rms_=option[lh::max_bbrms];
//	max_lh_radius_=option[lh::max_radius];
//	
//	loop_sizes_ = option[lh::loopsizes]();

	// initialize lhlibrary
	utility::vector1<core::Size> actual_lh_fragment_sizes(loop_sizes_.size());
	for(core::Size i=1; i<=loop_sizes_.size(); ++i)
	{
		actual_lh_fragment_sizes[i]=loop_sizes_[i]+(2*num_flanking_residues_to_match_);
	}
	lh_library_ = new LoopHashLibrary( actual_lh_fragment_sizes );
	lh_library_->load_mergeddb();
}

protocols::loops::Loop
LoophashLoopInserter::insert_loop(
	core::pose::Pose & pose,
	core::Size loop_anchor
){
	using namespace core;
	using namespace protocols::loophash;
	
	core::pose::Pose centroid_pose = pose;
	core::util::switch_to_residue_type_set( centroid_pose, core::chemical::CENTROID);
	
	core::Size lh_fragment_begin = loop_anchor-num_flanking_residues_to_match_+1;
	core::Size lh_fragment_end = loop_anchor+num_flanking_residues_to_match_;
	
	//Collect backbone segments from the specified number of residues before and after the loop
	//to be built
	BackboneSegment pose_bs_1;
	BackboneSegment pose_bs_2;
	pose_bs_1.read_from_pose( centroid_pose, lh_fragment_begin, num_flanking_residues_to_match_-1 );
	pose_bs_2.read_from_pose( centroid_pose, loop_anchor+1, num_flanking_residues_to_match_-1 );
	
	numeric::geometry::hashing::Real6 loop_transform;
	TR << "Getting transform from residues " << lh_fragment_begin << " and " << lh_fragment_end+1 << std::endl;
	if(!get_rt_over_leap_without_foldtree_bs( centroid_pose, lh_fragment_begin, lh_fragment_end+1, loop_transform )){
		utility_exit_with_message("Unable to find rigid body transform over jump");
	}
	
	core::Size num_filtered_fragments = 0;
	
	//Gradually increase radius size until the number of fragments meeting torsion RMSD requirments
	//is at least the minimum (for all loop sizes combined), or until the maximum radius is reached
	std::map< core::Size, std::vector<core::Size> > hash_buckets;
//	for( Size radius = 0; radius <= max_lh_radius_; radius++ )
//	{
		num_filtered_fragments = 0;
		for( core::Size i = 0; i < lh_library_->hash_sizes().size(); i++ )
		{
			core::Size loop_size = lh_library_->hash_sizes()[ i ];
			LoopHashMap &hashmap = lh_library_->gethash( loop_size );
			
			std::vector<core::Size> leap_index_bucket;
//			hashmap.radial_lookup( radius, loop_transform, leap_index_bucket);
			hashmap.radial_lookup( max_lh_radius_, loop_transform, leap_index_bucket);
			
			TR.Debug << "radius, loop_size, lookup_size = " << max_lh_radius_ << ", " <<
				loop_size << "," << leap_index_bucket.size() << std::endl;
			
			for(  std::vector<core::Size>::const_iterator it = leap_index_bucket.begin();
				it != leap_index_bucket.end();
				++it )
			{
				// Get the actual strucure index (not just the bin index)
				core::Size retrieve_index = (core::Size) (*it);
				LeapIndex cp = hashmap.get_peptide( retrieve_index );
				
				// Retrieve the backbone structure for the pre and post-loop segments
				BackboneSegment new_bs_1;
				BackboneSegment new_bs_2;
				
				//offset is angle-based, not residue based, hence the x3
				core::Size bs_2_offset = cp.offset+(((loop_size-num_flanking_residues_to_match_)+1)*3);
				
				//subtract 1 from num_flanking_residues_to_match_ because the last residue will have undefined torsions
				lh_library_->backbone_database().get_backbone_segment( cp.index, cp.offset,
					num_flanking_residues_to_match_-1 , new_bs_1 );
				lh_library_->backbone_database().get_backbone_segment( cp.index, bs_2_offset,
					num_flanking_residues_to_match_-1, new_bs_2 );
				
				// Check the values against against any RMS limitations
				core::Real bb_rms_1 = get_rmsd( pose_bs_1, new_bs_1 );
				core::Real bb_rms_2 = get_rmsd( pose_bs_2, new_bs_2 );
				
//				TR.Debug << "BB rmsds for fragment: " << bb_rms_1 << " " << bb_rms_2 << std::endl;
				
				if( (bb_rms_1 > min_torsion_rms_) && (bb_rms_1 < max_torsion_rms_ ) &&
					(bb_rms_2 > min_torsion_rms_) && (bb_rms_2 < max_torsion_rms_ ) )
				{
					hash_buckets[loop_size].push_back( *it );
					num_filtered_fragments++;
//					TR.Debug << "Accepted pre-jump rms : " << bb_rms_1 << std::endl;
//					TR.Debug << "Accepted post-jump rms : " << bb_rms_2 << std::endl;
				}
			}
						
			TR.Debug << num_filtered_fragments << " low-rmsd(" << min_torsion_rms_ <<
				", " << max_torsion_rms_ << ") fragments after loop size " << loop_size << std::endl;
//		}
//				
//		//not enough low-rmsd fragments, so increase radius if possible
//		if( num_filtered_fragments < min_filtered_fragments && radius != max_lh_radius_) {
//			hash_buckets.clear();
//			num_filtered_fragments=0;
//		}
//		//have enough fragments
//		else{
//			break;
//		}
	}
	if(num_filtered_fragments==0)
	{
		std::stringstream err;
		err << "No loophash fragments found for transform between residues " << loop_anchor << " and " << loop_anchor+1
			<< " consider loosening torsion rmsd or increasing max loophash radius size." << std::endl;
		utility_exit_with_message(err.str());
	}
	
	TR << "Found " << num_filtered_fragments << " fragments that meet torsion requirements. Picking one randomly to insert" << std::endl;

	//pick a random sized loop length
	core::Size hash_index = numeric::random::random_range(0, lh_library_->hash_sizes().size()-1);
	core::Size lh_fragment_size = lh_library_->hash_sizes()[ hash_index ];
	core::Size loop_size = lh_fragment_size - (2*num_flanking_residues_to_match_);
	TR << "Inserting a " << loop_size << " element loop" << std::endl;

	//pick a random fragment from the filtered fragments
	std::vector<Size> filtered_leap_index_bucket = hash_buckets[lh_fragment_size];
	core::Size bucket_index = numeric::random::random_range(0, filtered_leap_index_bucket.size()-1);
	core::Size retrieve_index = filtered_leap_index_bucket[bucket_index];
	
	LoopHashMap & hashmap = lh_library_->gethash( lh_fragment_size );
	LeapIndex cp = hashmap.get_peptide( retrieve_index );
	
	BackboneSegment lh_fragment_bs;
	lh_library_->backbone_database().get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , lh_fragment_bs );
	
	lh_fragment_bs.print();
	
	//get the extra data from the backbone database so we have a sequence to build
	BBData bb;
	BBExtraData extra_data;
	lh_library_->backbone_database().get_protein( cp.index, bb );
	lh_library_->backbone_database().get_extra_data(bb.extra_key, extra_data);
	std::string sequence = extra_data.sequence;
	TR << "PDB: " << extra_data.pdb_id << std::endl;
	TR << "Sequence: " << sequence << std::endl;
	
	core::Size seq_offset = (cp.offset/3)+(num_flanking_residues_to_match_*3);//should I be subtracting 1 here?
	
	core::Size n_term_append_size=numeric::random::random_range(1, loop_size-1);
	core::Size c_term_append_size=loop_size-n_term_append_size;
	TR.Debug << "Attaching " << n_term_append_size << " residues to N-termini and "
		<< c_term_append_size << " residues to C-termini" << std::endl;
	
	//append idealized residues to the pose
	//TODO: use the rotamer id from the backbone db extra data
	chemical::ResidueTypeSetCAP restype_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	for(core::Size i=0; i<n_term_append_size; ++i){
		
		core::Size append_seqpos = loop_anchor+i;
	
		char aa_char = sequence[seq_offset+i];
		chemical::AA aa = chemical::aa_from_oneletter_code(aa_char);
		chemical::ResidueTypeCOPs const & rsd_type_list( restype_set->aa_map( aa ) );
		if(rsd_type_list.empty()){
			utility_exit_with_message("Could not find residue type for AA: " + aa_char);
		}
		chemical::ResidueType const & rsd_type( *(rsd_type_list[ 1 ]) );
		conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( rsd_type ) );
		
		TR.Debug << "Attaching residue " << aa_char << " after seqpos " << append_seqpos << std::endl;
		pose.append_polymer_residue_after_seqpos(*new_rsd, append_seqpos, true);
		++lh_fragment_end;
	}
	core::Size prepend_seqpos = loop_anchor+loop_size-c_term_append_size+1;
	for(core::Size i=loop_size-1; i>=loop_size-c_term_append_size; --i){
		
		char aa_char = sequence[seq_offset+i];
		chemical::AA aa = chemical::aa_from_oneletter_code(aa_char);
		chemical::ResidueTypeCOPs const & rsd_type_list( restype_set->aa_map( aa ) );
		if(rsd_type_list.empty()){
			utility_exit_with_message("Could not find residue type for AA: " + aa_char);
		}
		chemical::ResidueType const & rsd_type( *(rsd_type_list[ 1 ]) );
		conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( rsd_type ) );
		
		TR.Debug << "Attaching residue " << aa_char << " before seqpos " << prepend_seqpos << std::endl;
		pose.prepend_polymer_residue_before_seqpos(*new_rsd, prepend_seqpos, true);
		++lh_fragment_end;
	}
	
	protocols::loops::Loop loop(
		loop_anchor+1,
		loop_anchor+loop_size,
		loop_anchor+n_term_append_size
	);
	
	//Apply the torsions of the loop residue to the pose
	if(modify_flanking_regions_){
		kinematics::FoldTree ft;
		
		ft.add_edge(1, lh_fragment_begin-1, kinematics::Edge::PEPTIDE);
		
		ft.add_edge(1, lh_fragment_begin, 1);
		ft.add_edge(lh_fragment_begin, loop_anchor+n_term_append_size, kinematics::Edge::PEPTIDE);
		
		ft.add_edge(1, lh_fragment_end, 2);
		ft.add_edge(lh_fragment_end, loop_anchor+n_term_append_size+1, kinematics::Edge::PEPTIDE);
		
		ft.add_edge(1, lh_fragment_end+1, 3);
		ft.add_edge(lh_fragment_end+1, pose.total_residue(), kinematics::Edge::PEPTIDE);
		
		if(!ft.check_fold_tree())
		{
			utility_exit_with_message("LoophashLoopInserter made a bad fold tree. File a bug!");
		}
		pose.fold_tree(ft);
		
		pose.dump_pdb("before_torsion_apply.pdb");
		TR.Debug << "Foldtree for torsion apply " << ft << std::endl;
		TR.Debug << "Applying " << lh_fragment_bs.length() << " lh fragment torsions starting at residue " << lh_fragment_begin << std::endl;
		lh_fragment_bs.apply_to_pose(pose, lh_fragment_begin, false);
	}
	else{
		std::vector<Real> full_phi = lh_fragment_bs.phi();
		std::vector<Real> full_psi = lh_fragment_bs.psi();
		std::vector<Real> full_omega = lh_fragment_bs.omega();
		
		std::vector<Real> loop_phi(full_phi.begin()+num_flanking_residues_to_match_, full_phi.end()-num_flanking_residues_to_match_);
		std::vector<Real> loop_psi(full_psi.begin()+num_flanking_residues_to_match_, full_psi.end()-num_flanking_residues_to_match_);
		std::vector<Real> loop_omega(full_omega.begin()+num_flanking_residues_to_match_, full_omega.end()-num_flanking_residues_to_match_);
		
		BackboneSegment loop_bs(loop_phi, loop_psi, loop_omega);
		TR.Debug << "Applying " << loop_bs.length() << " lh fragment torsions starting at residue " << loop_anchor+1 << std::endl;
		loop_bs.apply_to_pose(pose, loop_anchor+1, false);
	}
	
	return loop;
}

void
LoophashLoopInserter::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	
			
	/***********REQUIRED OPTIONS************/
	//Maximum RMSD of torsion angles to flanking residues
	if(tag->hasOption("max_torsion_rms") || option[lh::max_bbrms].user()){
		max_torsion_rms_ =
			tag->getOption<core::Real>("max_torsion_rms", option[lh::max_bbrms].value());
	}
	else{
		utility_exit_with_message("You must specify the maximum torsion rmsd using the max_torsion_rms tag or the lh::max_bbrms option");
	}
	
	//Minimum RMSD of torsion angles to flanking residues
	if(tag->hasOption("min_torsion_rms") || option[lh::min_bbrms].user()){
		min_torsion_rms_ =
			tag->getOption<core::Real>("min_torsion_rms", option[lh::min_bbrms].value());
	}
	else{
		utility_exit_with_message("You must specify the minimum torsion rmsd using the min_torsion_rms tag or the lh::min_bbrms option");
	}
	
	if(tag->hasOption("loop_sizes") || option[lh::loopsizes].user()){
		if(option[lh::loopsizes].user())
		{
			loop_sizes_=option[lh::loopsizes].value();
		}
		else
		{
			loop_sizes_.clear();
			utility::vector1<std::string> loop_sizes_strings =
				utility::string_split(tag->getOption< std::string >("loop_sizes"), ',');
			for(core::Size i=1; i<=loop_sizes_strings.size(); ++i)
			{
				loop_sizes_.push_back(utility::string2int(loop_sizes_strings[i]));
			}
		}
	}
	else{
		utility_exit_with_message("You must specify desired loop sizes through either the loop_sizes tag or the lh::loopsizes option");
	}
	
	/***********OPTIONAL OPTIONS************/
	modify_flanking_regions_ =
		tag->getOption<bool>("modify_flanking_regions", false);
		
	num_flanking_residues_to_match_ =
		tag->getOption<core::Size>("num_flanking_residues_to_match", 3);
	
	//Radius to use for loophash fragment lookup
	max_lh_radius_ =
		tag->getOption<core::Real>("max_lh_radius", option[lh::max_radius]);
		
	//INITIALIZE DB, PROBABLY SHOULDN'T DO THAT HERE
	utility::vector1<core::Size> actual_lh_fragment_sizes(loop_sizes_.size());
	for(core::Size i=1; i<=loop_sizes_.size(); ++i)
	{
		actual_lh_fragment_sizes[i]=loop_sizes_[i]+(2*num_flanking_residues_to_match_);
	}
	lh_library_ = new LoopHashLibrary( actual_lh_fragment_sizes );
	lh_library_->load_mergeddb();
}

} //protocols
} //loophash
