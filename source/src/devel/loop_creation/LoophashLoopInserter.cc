// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LoophashLoopInserter.cc
///
/// @brief Use loophash to find a fragment from loop_anchor to loop_anchor+1. Build ideal residues
/// using the loophash sequence. Don't idealize the last loop connections (this results in a broken loop)
/// @author Tim Jacobs

//Unit
//#include <protocols/loophash/LoophashLoopInserter.hh>
//#include <protocols/loophash/LoophashLoopInserterCreator.hh>
#include <devel/loop_creation/LoophashLoopInserter.hh>
#include <devel/loop_creation/LoophashLoopInserterCreator.hh>

//Numeric
#include <numeric/random/random.hh>
#include <numeric/geometry/hashing/SixDHasher.fwd.hh>

//Basic
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>

//Core
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
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
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_closure/ccd/ccd_closure.hh>

//utility
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

//numeric
#include <numeric/random/random_permutation.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif


//namespace protocols {
//namespace loophash {
namespace devel {
namespace loop_creation {

static THREAD_LOCAL basic::Tracer TR( "protocols.loophash.LoophashLoopInserter" );

//****CREATOR METHODS****//
std::string
LoophashLoopInserterCreator::keyname() const
{
	return LoophashLoopInserterCreator::mover_name();
}

protocols::moves::MoverOP
LoophashLoopInserterCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoophashLoopInserter );
}

std::string
LoophashLoopInserterCreator::mover_name()
{
	return "LoophashLoopInserter";
}
//****END CREATOR METHODS****//

LoophashLoopInserter::LoophashLoopInserter():
	max_closure_deviation_(5.0),
	num_flanking_residues_to_match_(3),
	modify_flanking_regions_(true),
	lh_initialized_(false)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	min_torsion_rms_=option[lh::min_bbrms].value();
	max_torsion_rms_=option[lh::max_bbrms].value();
	max_lh_radius_=option[lh::max_radius].value();
}

protocols::moves::MoverOP
LoophashLoopInserter::clone() const {
	return( protocols::moves::MoverOP( new LoophashLoopInserter( *this ) ) );
}
protocols::moves::MoverOP
LoophashLoopInserter::fresh_instance() const {
	return protocols::moves::MoverOP( new LoophashLoopInserter );
}

std::string
LoophashLoopInserter::get_name() const {
	return "LoophashLoopInserter";
}

void
LoophashLoopInserter::init(
	core::pose::Pose & pose
){
	using namespace core;
	using namespace protocols::loophash;
	using namespace basic::resource_manager;
	using namespace basic::options;

	if ( !lh_initialized_ ) {
		if ( ResourceManager::get_instance()->has_resource("LoopHashLibrary") ) {
			TR << "Retrieving lh library from resource manager." << std::endl;
			lh_library_ = get_resource<LoopHashLibrary>( "LoopHashLibrary" );
		} else {
			TR << "Initializing lh library from command line" << std::endl;
			utility::vector1<core::Size> actual_lh_fragment_sizes(loop_sizes_.size());
			for ( core::Size i=1; i<=loop_sizes_.size(); ++i ) {
				actual_lh_fragment_sizes[i]=loop_sizes_[i]+(2*num_flanking_residues_to_match_);
			}
			lh_library_ = protocols::loophash::LoopHashLibraryOP( new LoopHashLibrary( actual_lh_fragment_sizes ) );
			lh_library_->load_mergeddb();
		}
		lh_initialized_=true;
	}

	//Sanity checks
	if ( loop_anchor()<=0 || loop_anchor()>pose.size()-1 ) {
		std::stringstream err;
		err << "Loop anchor " << loop_anchor() << " is invalid" << std::endl;
		utility_exit_with_message(err.str());
	}
	if ( lh_library_->hash_sizes().size()<=0 ) {
		std::stringstream err;
		err << "No loop sizes for loophash" << std::endl;
		utility_exit_with_message(err.str());
	}

	TR << "Attempting build new loop residues between " << loop_anchor() << " and " << loop_anchor()+1 << std::endl;

	//Remove variants from anchor positions
	pose::remove_upper_terminus_type_from_pose_residue(pose, loop_anchor());
	pose::remove_lower_terminus_type_from_pose_residue(pose, loop_anchor()+1);
}

void
LoophashLoopInserter::apply(
	core::pose::Pose & pose
){
	using namespace core;
	using namespace protocols::loophash;

	init(pose);

	core::Size lh_fragment_begin = loop_anchor()-num_flanking_residues_to_match_+1;
	core::Size lh_fragment_end = loop_anchor()+num_flanking_residues_to_match_;

	HashBuckets hash_buckets = find_fragments(pose, lh_fragment_begin, lh_fragment_end);
	if ( hash_buckets.size()==0 ) {
		std::stringstream err;
		err << "No loophash fragments found for transform between residues " << loop_anchor() << " and " << loop_anchor()+1
			<< " consider loosening torsion rmsd or increasing max loophash radius size." << std::endl;
		utility_exit_with_message(err.str());
	}

	std::pair<core::Size, core::Size> random_fragment =
		get_random_fragment(hash_buckets);

	clock_t start_time = clock();
	std::pair<core::Real,core::Real> deviations =
		build_loop(pose, lh_fragment_begin, lh_fragment_end, random_fragment.first, random_fragment.second);
	clock_t build_time = clock() - start_time;
	TR.Debug << "Clocks - Build fragment in: " << build_time << std::endl;

	TR << "Deviations after initial loop insert: " << deviations.first << " " << deviations.second << std::endl;
}

/// @brief get a random fragment length and fragment retrieval index from the given hash bucket.
std::pair<core::Size,core::Size>
LoophashLoopInserter::get_random_fragment(
	HashBuckets hash_buckets
){
	//pick a random key index from the hash buckets
	core::Size random_key_index = numeric::random::random_range(0, hash_buckets.size()-1);

	//iterate this random number of times
	HashBuckets::const_iterator it=hash_buckets.begin();
	for ( core::Size i=1; i<=random_key_index; ++i ) {
		++it;
	}
	//Get the vector of retreival indexes corresponding to this random key
	core::Size lh_fragment_size = it->first;
	std::vector<Size> filtered_leap_index_bucket = it->second;

	//pick a random retrieval index (fragment) from the filtered list
	core::Size bucket_index = numeric::random::random_range(0, filtered_leap_index_bucket.size()-1);
	core::Size retrieve_index = filtered_leap_index_bucket[bucket_index];
	TR.Debug << "Random retrieval index: " << retrieve_index << std::endl;

	return std::make_pair(lh_fragment_size, retrieve_index);
}

LoophashLoopInserter::HashBuckets
LoophashLoopInserter::find_fragments(
	core::pose::Pose const & pose,
	core::Size lh_fragment_begin,
	core::Size lh_fragment_end
){
	core::Size max_size=0;
	core::Size min_size=100000;
	for ( core::Size const & i : lh_library_->hash_sizes() ) {
		max_size = std::max(max_size, i);
		min_size = std::min(min_size, i);
	}
	return find_fragments(pose, lh_fragment_begin, lh_fragment_end, min_size, max_size);
}

LoophashLoopInserter::HashBuckets
LoophashLoopInserter::find_fragments(
	core::pose::Pose const & pose,
	core::Size lh_fragment_begin,
	core::Size lh_fragment_end,
	core::Size min_fragment_size,
	core::Size max_fragment_size
){
	using namespace protocols::loophash;
	clock_t start_time = clock();

	core::pose::Pose centroid_pose = pose;
	core::util::switch_to_residue_type_set( centroid_pose, core::chemical::CENTROID_t );

	//Collect backbone segments from the specified number of residues before and after the loop
	//to be built. Don't use the last residue before and after the jump because it will have an undefined
	//psi (n-terminal pre-jump residue) or phi (c-terminal post-jump residue)
	BackboneSegment pose_bs_1;
	pose_bs_1.read_from_pose( centroid_pose, lh_fragment_begin, num_flanking_residues_to_match_-1 );
	TR << "Attemping to find fragments matching " << num_flanking_residues_to_match_-1
		<< " pre-jump residues, starting at " << lh_fragment_begin << std::endl;

	BackboneSegment pose_bs_2;
	pose_bs_2.read_from_pose( centroid_pose, loop_anchor()+2, num_flanking_residues_to_match_-1 );
	TR << "Attemping to find fragments matching " << num_flanking_residues_to_match_-1
		<< " pre-jump residues, starting at " << loop_anchor()+2 << std::endl;

	numeric::geometry::hashing::Real6 loop_transform;
	TR << "Getting transform from residues " << lh_fragment_begin << " and " << lh_fragment_end+1 << std::endl;
	if ( !get_rt_over_leap_without_foldtree_bs( centroid_pose, lh_fragment_begin, lh_fragment_end+1, loop_transform ) ) {
		utility_exit_with_message("Unable to find rigid body transform over jump");
	}
	TR << "Max lh fragment size: " << max_fragment_size << std::endl;

	HashBuckets hash_buckets;
	core::Size num_filtered_fragments = 0;
	for ( core::Size i = 0; i < lh_library_->hash_sizes().size(); i++ ) {
		core::Size loop_size = lh_library_->hash_sizes()[ i ];
		if ( loop_size > max_fragment_size ) { continue; }
		if ( loop_size < min_fragment_size ) { continue; }

		LoopHashMap &hashmap = lh_library_->gethash( loop_size );

		std::vector<core::Size> leap_index_bucket;
		hashmap.radial_lookup( core::Size(max_lh_radius_), loop_transform, leap_index_bucket);

		TR.Debug << "radius, loop_size, lookup_size = " << max_lh_radius_ << ", " <<
			loop_size << "," << leap_index_bucket.size() << std::endl;

		for (  std::vector<core::Size>::const_iterator it = leap_index_bucket.begin();
				it != leap_index_bucket.end();
				++it ) {
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

			if ( (bb_rms_1 > min_torsion_rms_) && (bb_rms_1 < max_torsion_rms_ ) &&
					(bb_rms_2 > min_torsion_rms_) && (bb_rms_2 < max_torsion_rms_ ) ) {
				hash_buckets[loop_size].push_back( *it );
				num_filtered_fragments++;
			}
		}

		TR.Debug << num_filtered_fragments << " low-rmsd(" << min_torsion_rms_ <<
			", " << max_torsion_rms_ << ") fragments after loop size " << loop_size << std::endl;
	}

	clock_t fragment_time = clock() - start_time;
	TR.Debug << "Clocks - Found " << num_filtered_fragments << " filtered fragments (" <<
		hash_buckets.size() << " loop sizes) for jump from " << lh_fragment_begin << " to " << lh_fragment_end
		<< " - time: " << fragment_time << std::endl;

	return hash_buckets;
}

std::pair<core::Real,core::Real>
LoophashLoopInserter::build_loop(
	core::pose::Pose & pose,
	core::Size lh_fragment_begin,
	core::Size lh_fragment_end,
	core::Size lh_fragment_size,
	core::Size retrieve_index
){
	using namespace core;
	using namespace protocols::loophash;

	core::Size loop_size = lh_fragment_size - (2*num_flanking_residues_to_match_);
	TR << "Inserting a " << loop_size << " element loop" << std::endl;

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

	TR.Debug << "Loop size: " << loop_size << ", attaching " << n_term_append_size << " residues to N-termini and "
		<< c_term_append_size << " residues to C-termini" << std::endl;

	//append idealized residues to the pose
	//TODO: use the rotamer id from the backbone db extra data
	chemical::ResidueTypeSetCOP restype_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	for ( core::Size i=0; i<n_term_append_size; ++i ) {

		core::Size append_seqpos = loop_anchor()+i;

		char aa_char = sequence[seq_offset+i];
		chemical::AA aa = chemical::aa_from_oneletter_code(aa_char);
		chemical::ResidueTypeCOP rsd_type( restype_set->get_representative_type_aa( aa ) );
		if ( ! rsd_type ) {
			std::string err = "Could not find residue type for AA: ";
			err+=aa_char;
			utility_exit_with_message(err);
		}
		conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( *rsd_type ) );

		TR.Debug << "Attaching residue " << aa_char << " after seqpos " << append_seqpos << std::endl;
		pose.conformation().safely_append_polymer_residue_after_seqpos(*new_rsd, append_seqpos, true);
		++lh_fragment_end;
	}

	//Prevent Size wrap around
	core::Size prepend_seqpos = loop_anchor()+loop_size-c_term_append_size+1;
	if ( loop_size-c_term_append_size > 0 ) {
		for ( core::Size i=loop_size-1; i>=loop_size-c_term_append_size; --i ) {

			char aa_char = sequence[seq_offset+i];
			chemical::AA aa = chemical::aa_from_oneletter_code(aa_char);
			chemical::ResidueTypeCOP rsd_type( restype_set->get_representative_type_aa( aa ) );
			if ( ! rsd_type ) {
				std::string err = "Could not find residue type for AA: ";
				err+=aa_char;
				utility_exit_with_message(err);
			}
			conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( *rsd_type ) );

			TR.Debug << "Attaching residue " << aa_char << " before seqpos " << prepend_seqpos << std::endl;
			pose.conformation().safely_prepend_polymer_residue_before_seqpos(*new_rsd, prepend_seqpos, true);
			++lh_fragment_end;
		}
	}

	protocols::loops::Loop loop(
		loop_anchor()+1,
		loop_anchor()+loop_size,
		loop_anchor()+n_term_append_size
	);

	//Apply the torsions of the loop residue to the pose
	if ( modify_flanking_regions_ ) {
		modified_range(lh_fragment_begin, lh_fragment_end);
		protocols::loops::Loop extended_loop(lh_fragment_begin, lh_fragment_end, loop.cut());
		protocols::loops::set_single_loop_fold_tree(pose, extended_loop);


		TR.Debug << "Applying " << lh_fragment_bs.length() << " lh fragment torsions starting at residue " << lh_fragment_begin << std::endl;
		lh_fragment_bs.apply_to_pose(pose, lh_fragment_begin, false);
	} else {
		modified_range(loop.start(), loop.stop());
		protocols::loops::set_single_loop_fold_tree(pose, loop);

		std::vector<Real> full_phi = lh_fragment_bs.phi();
		std::vector<Real> full_psi = lh_fragment_bs.psi();
		std::vector<Real> full_omega = lh_fragment_bs.omega();

		std::vector<Real> loop_phi(full_phi.begin()+num_flanking_residues_to_match_, full_phi.end()-num_flanking_residues_to_match_);
		std::vector<Real> loop_psi(full_psi.begin()+num_flanking_residues_to_match_, full_psi.end()-num_flanking_residues_to_match_);
		std::vector<Real> loop_omega(full_omega.begin()+num_flanking_residues_to_match_, full_omega.end()-num_flanking_residues_to_match_);

		BackboneSegment loop_bs(loop_phi, loop_psi, loop_omega);
		TR.Debug << "Applying " << loop_bs.length() << " lh fragment torsions starting at residue " << loop_anchor()+1 << std::endl;
		loop_bs.apply_to_pose(pose, loop_anchor()+1, false);
	}

	created_loop_=loop;

	return protocols::loops::loop_closure::ccd::get_deviation(pose, loop.cut());
}

void
LoophashLoopInserter::parse_my_tag(
	utility::tag::TagCOP tag, basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::loophash;

	parse_loop_anchor(tag);

	//Maximum RMSD of torsion angles to flanking residues
	if ( tag->hasOption("max_torsion_rms") || option[lh::max_bbrms].user() ) {
		max_torsion_rms_ =
			tag->getOption<core::Real>("max_torsion_rms", option[lh::max_bbrms].value());
	} else {
		utility_exit_with_message("You must specify the maximum torsion rmsd using the max_torsion_rms tag or the lh::max_bbrms option");
	}

	//Minimum RMSD of torsion angles to flanking residues
	if ( tag->hasOption("min_torsion_rms") || option[lh::min_bbrms].user() ) {
		min_torsion_rms_ =
			tag->getOption<core::Real>("min_torsion_rms", option[lh::min_bbrms].value());
	} else {
		utility_exit_with_message("You must specify the minimum torsion rmsd using the min_torsion_rms tag or the lh::min_bbrms option");
	}

	max_closure_deviation_ =
		tag->getOption<core::Real>("max_closure_deviation", 1);

	if ( tag->hasOption("loop_sizes") || option[lh::loopsizes].user() ) {
		if ( tag->hasOption("loop_sizes") ) {
			loop_sizes_.clear();
			utility::vector1<std::string> loop_sizes_strings =
				utility::string_split(tag->getOption< std::string >("loop_sizes"), ',');
			for ( core::Size i=1; i<=loop_sizes_strings.size(); ++i ) {
				loop_sizes_.push_back(utility::string2int(loop_sizes_strings[i]));
			}
		} else {
			loop_sizes_=option[lh::loopsizes].value();
		}
	} else {
		utility_exit_with_message("You must specify desired loop sizes through either the loop_sizes tag or the lh::loopsizes option");
	}

	if ( tag->hasOption("max_lh_radius") || option[lh::max_radius].user() ) {
		max_lh_radius_ =
			tag->getOption<core::Size>("max_lh_radius", option[lh::max_radius].value());
	} else {
		utility_exit_with_message("You must specify max radius for loophash through either the max_radius tag or the lh::max_radius option");
	}

	if ( tag->hasOption("modify_flanking_regions") ) {
		modify_flanking_regions_ =
			tag->getOption<bool>("modify_flanking_regions");
	}

	if ( tag->hasOption("num_flanking_residues_to_match") ) {
		num_flanking_residues_to_match_ =
			tag->getOption<core::Size>("num_flanking_residues_to_match");
	}

}

} //loop creation
} //devel
//} //protocols
//} //loophash
