// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopHashSampler.cc
/// @brief
/// @author Mike Tyka

#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <protocols/loophash/BackboneDB.hh>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <utility/string_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>

#include <utility/vector1.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>


#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif


namespace protocols {
namespace loophash {

static thread_local basic::Tracer TR( "LocalHashSampler" );

LoopHashSampler::LoopHashSampler(
	LoopHashLibraryOP library,
	LocalInserterOP inserter
):
	library_(library),
	inserter_(inserter),
	start_res_ ( 2 ),
	stop_res_  ( 0 ),
	min_bbrms_ ( 0.0 ),
	max_bbrms_ ( 100000.0 ),
	min_rms_   ( 0.0 ),
	max_rms_   ( 100.0 ),
	max_struct_ (10),
	max_struct_per_radius_ (10),
	nonideal_ ( false ),
	nprefilter_ ( 0 ) // OBSOLETE?
{
	set_defaults();
}

LoopHashSampler::~LoopHashSampler() {}

void
LoopHashSampler::set_defaults(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	set_max_radius(  option[ lh::max_radius ]() );

	set_min_bbrms(  option[ lh::min_bbrms ]() );
	set_max_bbrms(  option[ lh::max_bbrms ] ()  );
	set_min_rms(    option[ lh::min_rms ]() );
	set_max_rms(    option[ lh::max_rms ]() );
	set_max_struct(    option[ lh::max_struct ]() );
	set_max_struct_per_radius(    option[ lh::max_struct_per_radius ]() );
	set_max_nstruct( 10000000 ); // OBSOLETE?

	filter_by_phipsi_ = option[ lh::filter_by_phipsi ]();
}

bool cmp( core::pose::Pose a, core::pose::Pose b) {
	std::string scoreterm = "censcore";
	core::Real as, bs;
	core::pose::getPoseExtraScore( a, scoreterm, as ); core::pose::getPoseExtraScore( b, scoreterm, bs);
	return as < bs;
}

// returns a vector of real numbers, one per residue of pose, giving the
// sampling weight. weight of 1.0 corresponds to "normal", i.e. unmodified sampling weight.
utility::vector1 < core::Real > extract_sample_weights( const core::pose::Pose &pose ){
	std::string sample_weight_str;
	core::pose::get_comment(pose, "sample_weight", sample_weight_str);

	utility::vector1 < std::string > sample_weight_input_parameters;
	sample_weight_input_parameters = utility::split(sample_weight_str);

	utility::vector1 < core::Real > sample_weight;
	for ( core::Size res_count = 1; res_count <= pose.total_residue(); ++res_count ) {
		core::Real new_sample_weight = 1.0;
		if ( res_count < sample_weight_input_parameters.size() ) {
			new_sample_weight = utility::string2float( sample_weight_input_parameters[res_count] );
		}
		sample_weight.push_back( new_sample_weight );
	}

	return sample_weight;
}


bool is_valid_backbone(
	const std::string &sequence,
	const core::Size &ir,  // sequence offset
	const std::vector< core::Real > &phi,
	const std::vector< core::Real > &psi,
	bool &filter_pro,
	bool &filter_beta,
	bool &filter_gly
){
	runtime_assert( phi.size() == psi.size() )

		// Check phi/psi angles against the sequence
		// Pose counts residues starting from one, so offset that
		filter_pro = false;
	filter_beta = false;
	filter_gly = false;

	// now check every residue
	for ( core::Size bs_position = 0; bs_position < phi.size() ; ++bs_position ) {
		int sequence_position = ir - 1 + bs_position;

		// Proline
		if ( sequence[sequence_position] == 'P' ) {
			if ( phi[bs_position] < -103 || phi[bs_position] > -33 ) filter_pro = true;
		}
		// Beta branched residues
		if ( sequence[sequence_position] == 'I' || sequence[sequence_position] == 'V' || sequence[sequence_position] == 'T' ) {
			if ( phi[bs_position] > -40 ) filter_beta = true;
		}
		// Non glycine residues are confined to only part of the positive phi region
		// populated by glycine residues
		if ( sequence[sequence_position] != 'G' ) {
			if ( phi[bs_position] > 70 ) filter_gly = true;
		}
		if ( sequence[sequence_position] != 'G' ) {
			if ( psi[bs_position] < -75 && psi[bs_position] > -170 ) filter_gly = true;
		}
	}

	// were any of the filters triggered ? only return true if all filters are false!
	return !( filter_pro || filter_beta || filter_gly );
}

// Just a handy datastructure to carry over some statistics together with the actual retrieve index
struct FilterBucket {
	FilterBucket():
		retrieve_index(0),
		BBrms(0),
		filter_pro(false),
		filter_beta(false),
		filter_gly(false)
	{}
	core::Size retrieve_index;
	core::Real BBrms;
	bool filter_pro;
	bool filter_beta;
	bool filter_gly;
};

// @brief create a set of structures for a the given range of residues and other parameters
void
LoopHashSampler::build_structures(
	const core::pose::Pose& start_pose,
	std::vector< core::io::silent::SilentStructOP > &lib_structs
)
{
	using namespace core;
	using namespace core::pose;
	using namespace conformation;
	using namespace kinematics;
	using namespace numeric::geometry::hashing;
	using namespace optimization;
	using namespace id;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	runtime_assert( library_ != 0 );

	long starttime = time(NULL);

	// Statistics counters
	Size count_filter_rejects = 0;
	Size count_total_loops = 0;
	Size count_loop_builds = 0;
	Size count_filter_pro = 0;
	Size count_filter_beta = 0;
	Size count_filter_gly = 0;
	Size count_rejected_carms = 0;
	Size count_rejected_bbrms = 0;
	Size count_max_rad = 0;

	// Parameters
	core::Size models_build_this_loopsize_max =         std::max( Size(1), Size( max_struct_ / library_->hash_sizes().size()) );
	core::Size models_build_this_loopsize_per_rad_max = std::max( Size(1), Size( models_build_this_loopsize_max * 2 / max_radius_ ));
	core::Size fragments_tried_this_loopsize_max =      models_build_this_loopsize_max * 200;
	TR <<  "LoopHashSampler limits: " << max_struct_ << " " << models_build_this_loopsize_max << "  " << models_build_this_loopsize_per_rad_max << "  " << fragments_tried_this_loopsize_max << std::endl;

	std::string sequence = start_pose.sequence();

	core::pose::Pose original_pose = start_pose;
	core::pose::Pose edit_pose = start_pose;
	core::optimization::MinimizerOptions options( "lbfgs_armijo", 0.2, true , false );
	core::optimization::MinimizerOptions options2( "lbfgs_armijo", 0.02,true , false );

	kinematics::MoveMap final_mm;
	final_mm.set_bb(true);

	Size nres = start_pose.total_residue();
	Size ir, jr;

	core::Size start_res = start_res_;
	core::Size stop_res = stop_res_;

	// figure out start and stop residues
	if ( stop_res == 0 ) stop_res = nres;  // to do the whole protein just set stop_res to 0
	start_res = std::max( start_res, (core::Size)2 ); // dont start before 2  - WHY ? << cos you need a stub of at least 2 residues to calculate a proper takeoff point. Why ? I dont know. Rosettavoodoo. This knowledge has been lost in history. Historians have struggled for centuries to recover it.

	if ( start_res > stop_res ) stop_res = start_res;

	TR << "Running: Start:" << start_res << "  End: " << stop_res << std::endl;
	for ( ir = start_res; ir <= stop_res; ir ++ ) {

		// Loop over loopsizes in library
		for ( core::Size k = 0; k < library_->hash_sizes().size(); k ++ ) {
			core::Size loop_size = library_->hash_sizes()[ k ];

			jr = ir + loop_size;
			if ( ir > nres ) continue;
			if ( jr > nres ) continue;

			// get the rigid body transform for the current segment
			BackboneSegment pose_bs;
			pose_bs.read_from_pose( start_pose, ir, loop_size );
			Real6 loop_transform;
			if ( !get_rt_over_leap( original_pose, ir, jr, loop_transform ) ) continue;

			LoopHashMap &hashmap = library_->gethash( loop_size );


			// Now we compute the per residue sample weight averaged over the segment

			// we want models_build_this_loopsize_max models no matter what
			// but if there is a brrms constraint, we might never reach x models
			// some breakpoint, like x bins checked or x frags checked
			// or radius check
			core::Size fragments_tried_this_loopsize = 0;
			core::Size models_build_this_loopsize = 0;

			for ( Size radius = 0; radius <= max_radius_; radius++ ) {
				count_max_rad = std::max( count_max_rad, radius );
				core::Size models_build_this_loopsize_this_rad = 0;
				std::vector < core::Size > leap_index_bucket;
				std::vector < FilterBucket > filter_leap_index_bucket;

				hashmap.radial_lookup( radius, loop_transform, leap_index_bucket );     // grab list of fragments using radial lookup out from our loop transform
				TR.Debug << "Rad: " << radius << "  " << leap_index_bucket.size() << std::endl;
				if ( leap_index_bucket.size() == 0 )  continue;                           // no fragments found

				// Now for every hit, get the internal coordinates and make a short list of replacement loops
				// according to the RMS criteria
				for (  std::vector < core::Size >::const_iterator it = leap_index_bucket.begin();
						it != leap_index_bucket.end();
						++it ) {

					// Get the actual strucure index (not just the bin index)
					core::Size retrieve_index = (core::Size) (*it);
					LeapIndex cp = hashmap.get_peptide( retrieve_index );

					// Retrieve the actual backbone structure
					BackboneSegment new_bs;
					library_->backbone_database().get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );

					// Check the values against against any RMS limitations
					// if violated then skip rest of loop
					core::Real BBrms = get_rmsd( pose_bs, new_bs );
					if ( ( BBrms < min_bbrms_) || ( BBrms > max_bbrms_ ) ) {
						count_rejected_bbrms ++;
						continue;
					}


					FilterBucket bucket;
					bucket.retrieve_index = *it;  // save the bucket index for the next step later
					bucket.BBrms = BBrms;         // also save the back bone RMS for later analysis & stats

					bool is_valid =
						is_valid_backbone( sequence, ir, new_bs.phi(), new_bs.psi(),    // input is sequence, current position in sequence, and the phi/psi's of the proposed angles.
						bucket.filter_pro, bucket.filter_beta, bucket.filter_gly );  // output is a bunch of booleans giving information about any clashes.

					// count rejection stats
					if ( bucket.filter_pro )  count_filter_pro ++;
					if ( bucket.filter_beta ) count_filter_beta ++;
					if ( bucket.filter_gly )  count_filter_gly ++;

					if ( (!get_filter_by_phipsi()) || is_valid ) { // should we filter at all and if so is it valid.
						filter_leap_index_bucket.push_back( bucket ); // add to our short list of good fragments
					} else {
						count_filter_rejects++;            // or increment reject counter
					}

					count_total_loops++;
					fragments_tried_this_loopsize++;
					if ( fragments_tried_this_loopsize > fragments_tried_this_loopsize_max ) break; // continue with however many are in the bucket now, and break at end
				}

				// treat the fragments in a random order so shuffle them up
				//std::random__shuffle( filter_leap_index_bucket.begin(), filter_leap_index_bucket.end());
				numeric::random::random_permutation(filter_leap_index_bucket.begin(), filter_leap_index_bucket.end(), numeric::random::rg());

				// Now create models and check rms after insertion
				for (  std::vector < FilterBucket >::const_iterator it = filter_leap_index_bucket.begin();
						it != filter_leap_index_bucket.end();
						++it ) {

					clock_t starttime = clock();

					core::Size retrieve_index = it->retrieve_index;
					LeapIndex cp = hashmap.get_peptide( retrieve_index );

					BackboneSegment new_bs;
					library_->backbone_database().get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );

					core::pose::Pose newpose( start_pose );
					//transfer_phi_psi( start_pose, newpose );   //fpd necessary??

					core::Real final_rms = inserter_->make_local_bb_change( newpose, original_pose, new_bs, ir );
					count_loop_builds++;

					bool isok = false;
					if ( ( final_rms < max_rms_ ) && ( final_rms > min_rms_) ) {

						core::pose::Pose mynewpose( start_pose );

						transfer_phi_psi( newpose, mynewpose );
						transfer_jumps( newpose, mynewpose );

						core::io::silent::SilentStructOP new_struct = nonideal_ ?
							core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary") :
							core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
						new_struct->fill_struct( mynewpose );    // make the silent struct from the copy pose
						new_struct->energies_from_pose( newpose ); // take energies from the modified pose, not the copy pose
						new_struct->add_energy( "lh_carms", final_rms );
						new_struct->add_energy( "lh_bbrms", it->BBrms );
						new_struct->add_energy( "lh_radius", radius );
						new_struct->add_energy( "lh_loopsize", loop_size );
						new_struct->add_energy( "lh_filter_pro", it->filter_pro );
						new_struct->add_energy( "lh_filter_beta", it->filter_beta );
						new_struct->add_energy( "lh_filter_gly", it->filter_gly );

						//TR << "SAMPLER: " << new_struct->get_energy("censcore") << std::endl;
						// Add donor history for this round of loophash only

						// Assume extra data is loade, because we need it!
						BBData bb;
						BBExtraData bbextra;
						library_->backbone_database().get_protein( cp.index, bb );

						std::string donorhistory = new_struct->get_comment("donorhistory");
						if ( library_->backbone_database().extra_size() <= bb.extra_key ) {
							std::cerr << "ERROR: No extra data ?: " << library_->backbone_database().extra_size() << " < " << bb.extra_key << std::endl;

							donorhistory = donorhistory
								+       utility::to_string( loop_size )
								+ "/" + utility::to_string( library_->loopdb_range().first + cp.index )
								+ "/" + utility::to_string( cp.offset/3 )
								+ "/" + utility::to_string( ir-1 ) + ";";
						} else {

							library_->backbone_database().get_extra_data( bb.extra_key, bbextra );

							donorhistory = donorhistory
								+       utility::to_string( loop_size )
								+ "/" + utility::to_string( bbextra.pdb_id)
								//+ "/" + utility::to_string( bbextra.sequence.substr(cp.offset/3, 5) )
								+ "/" + utility::to_string( library_->loopdb_range().first + cp.index )
								+ "/" + utility::to_string( cp.offset/3 )
								+ "/" + utility::to_string( ir-1 ) + ";";
						}
						new_struct->erase_comment( "donorhistory" );
						new_struct->add_comment( "donorhistory", donorhistory );
						lib_structs.push_back( new_struct );

						models_build_this_loopsize++;
						models_build_this_loopsize_this_rad++;


						isok = true;
					} else {
						count_rejected_carms ++;
					}

					//if ( lib_structs.size() > 2  ) return;

					clock_t endtime = clock();

					TR.Debug << "Clocks: " << endtime - starttime << "  " << final_rms << (isok ? " OK" : " Reject") << std::endl;

					if ( models_build_this_loopsize >= models_build_this_loopsize_max ) break;
					if ( models_build_this_loopsize_this_rad >= models_build_this_loopsize_per_rad_max ) break;

				}
				// To break out of the outer for loop when these conditions are met
				if ( models_build_this_loopsize >= models_build_this_loopsize_max ) break;
				if ( fragments_tried_this_loopsize >= fragments_tried_this_loopsize_max ) break;
			}

			TR.Debug << " IR: " << ir << " LS: " << loop_size
				<< " Frag: " << fragments_tried_this_loopsize << " ( " << fragments_tried_this_loopsize_max << " ) "
				<< " Modls: " << models_build_this_loopsize   << " ( " << models_build_this_loopsize_max << " ) "
				<< std::endl;


		} // Loop over fragment sizes
	} // Loop iver residue window


	// Now just print some final statistics
	long endtime = time(NULL);
	TR.Info << "LHS: " << start_res << "-" << stop_res << ":  "
		<< " struc " << lib_structs.size()
		<< " (max) " << max_struct_
		<< " secs " << endtime - starttime << " secs "
		<< " Total: "   << count_total_loops
		<< " RjTor: "  << count_filter_rejects
		<< " RjPro: "  << count_filter_pro
		<< " RjBeta: " << count_filter_beta
		<< " RjGly: "  << count_filter_gly
		<< " RjCA("    << min_rms_ << "-" << max_rms_ << "): "   << count_rejected_carms
		<< " RjBB: "   << count_rejected_bbrms
		<< " Built: "   << count_loop_builds
		<< " MaxRad: "  << count_max_rad

		<< std::endl;

	for ( std::vector< core::io::silent::SilentStructOP >::iterator it=lib_structs.begin();
			it != lib_structs.end(); ++it ) {
		TR.Debug << "Samples: " << (*it)->get_energy("censcore") << std::endl;
	}


}

// closes gaps.  focuses on decreasing variability, not for use with loophash mpi
// very very similar to build_structures(), but too lazy to make separate functions to reduce redundancy
void
LoopHashSampler::close_gaps(
	const core::pose::Pose& /*start_pose*/,
	std::vector< core::pose::Pose> &/*lib_structs*/,
	core::Size /*loop_size*/
)
{
}


} // namespace loops
} // namespace protocols


