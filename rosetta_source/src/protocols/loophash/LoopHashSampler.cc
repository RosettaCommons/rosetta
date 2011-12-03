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
// AUTO-REMOVED #include <utility/sort_predicates.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>

#include <utility/vector1.hh>


#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif






namespace protocols {
namespace loophash {

  static basic::Tracer TR("LocalHashSampler");

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
	nprefilter_ ( 0 ), // OBSOLETE?
	nonideal_ ( false )
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
	core::pose::getPoseExtraScores( a, scoreterm, as ); core::pose::getPoseExtraScores( b, scoreterm, bs);
	return as < bs;
}

  // @brief create a set of structures for a the given range of residues and other parameters
 void
 LoopHashSampler::build_structures(
		const core::pose::Pose& start_pose,
    std::vector< core::io::silent::SilentStructOP > &lib_structs,
		core::Size round
	)
  {
    using namespace core;
    using namespace core::pose;
    using namespace conformation;
    using namespace kinematics;
    using namespace numeric::geometry::hashing;
    using namespace optimization;
    using namespace id;

		runtime_assert( library_ );
		TR.Trace << "Testing libstructs: " << "NStruct: " << lib_structs.size() << std::endl;


		long starttime = time(NULL);

		Size impossible_torsion_rejects = 0;
		Size impossible_torsion_candidates = 0;
    std::string sequence = start_pose.sequence();

    core::pose::Pose original_pose = start_pose;
    core::pose::Pose edit_pose = start_pose;
	 	TR.Debug  << "Setting options.." << std::endl;
	  //core::optimization::MinimizerOptions options( "dfpmin", 0.2, true , false );
	  //core::optimization::MinimizerOptions options2( "dfpmin", 0.02,true , false );
	  core::optimization::MinimizerOptions options( "lbfgs_armijo", 0.2, true , false );
	  core::optimization::MinimizerOptions options2( "lbfgs_armijo", 0.02,true , false );

    kinematics::MoveMap final_mm;
    final_mm.set_bb(true);
   // setup movemap & minimisation

  //				for ( Size ii=ir; ii<= jr; ++ii ) {
  //					final_mm.set_bb( ii, true );
  //					if ( newpose->residue(ii).aa() == chemical::aa_pro ) final_mm.set( TorsionID( phi_torsion, BB, ii ), false );
  //				}

    Size nres = start_pose.total_residue();
    Size ir, jr;
    //Size newpep_index = 0;

    //core::Size backbone_offset;
    //bbdb_.add_pose( pose, backbone_offset );

		// Get the sample weight array, split it and cast it as int
		// rosetta string functions suck
		std::string sample_weight_str;
		core::pose::get_comment(start_pose, "sample_weight", sample_weight_str);
		std::list < std::string > t;
		t = utility::split_to_list(sample_weight_str);
		utility::vector1 < core::Size > sample_weight;
		for (std::list<std::string>::const_iterator iterator = t.begin(), end = t.end(); iterator != end; ++iterator) {
				    sample_weight.push_back( utility::string2int( *iterator ) );
		}


    int runcount=0;
    runcount++;

    core::Size start_res = start_res_;
    core::Size stop_res = stop_res_;

    // figure out start and stop residues
    if ( stop_res == 0 ) stop_res = nres;  // to do the whole protein just set stop_res to 0
    start_res = std::max( start_res, (core::Size)2 ); // dont start before 2  - WHY ?
		if( start_res > stop_res ) stop_res = start_res;

	 	TR << "Running: Start:" << start_res << "  End: " << stop_res << std::endl;
    for( ir = start_res; ir <= stop_res; ir ++ ){
      for( core::Size k = 0; k < library_->hash_sizes().size(); k ++ ){
        core::Size loop_size = library_->hash_sizes()[ k ];

        jr = ir + loop_size;
        if ( ir > nres ) continue;
        if ( jr > nres ) continue;

        // get the rigid body transform for the current segment
        BackboneSegment pose_bs;
        pose_bs.read_from_pose( start_pose, ir, loop_size );
        Real6 loop_transform;
        if(!get_rt_over_leap( original_pose, ir, jr, loop_transform )) continue;

        LoopHashMap &hashmap = library_->gethash( loop_size );


				// Now we compute the per residue sample weight averaged over the segment
				core::Real avg_sw = 0;
				for( core::Size m = ir; m <= jr; m++ ) {
						avg_sw += sample_weight[m];
				}
				avg_sw = avg_sw/loop_size;
				if (avg_sw == 0 ) {
					//TR << "Sample Weight is 0" << std::endl;
					break;
				}

				// make up some function that uses sample weight to generate model number cutoffs
				// and one that gives a max_bbrms and min_bbrms
				// Make it slightly dependent on round
				//core::Size sw_nmodels = (int)(avg_sw/5*(1.0+(round-1)/60));
				core::Size sw_nmodels = (int)(avg_sw*max_struct_/50*(1.0+(round-1)/60));
				core::Size sw_nfrags = (int)avg_sw*400;
				
				// Limit how many structures chosen from a given radius
				// Make it dependent on round, less for higher rounds
				//core::Size sw_nmodels_per_rad = (int)(avg_sw/50*3);
				core::Size sw_nmodels_per_rad = (int)(avg_sw*max_struct_per_radius_/50);

				/* dont change these based on avg_sw for now
				set_min_rms(  avg_sw/10   );
				set_max_rms(  avg_sw/8  );
				set_min_bbrms( avg_sw/20  );  
				set_max_bbrms(  avg_sw/12  );
				*/


				// we want sw_nmodels models no matter what
				// but if there is a brrms constraint, we might never reach x models
				// some breakpoint, like x bins checked or x frags checked
				// or radius check
				core::Size fragments_seen = 0;
				core::Size models_seen = 0;

				for( Size radius = 0; radius <= max_radius_; radius++ ) {
					core::Size models_seen_this_rad = 0;
					std::vector < core::Size > leap_index_bucket;
					std::vector < core::Size > filter_leap_index_bucket;
					hashmap.radial_lookup( radius, loop_transform, leap_index_bucket );
					//TR << "radius, lookup_size = " << radius << ", " << leap_index_bucket.size() << std::endl;

					if( leap_index_bucket.size() == 0) {
						continue;
					}

					// Now for every hit, get the internal coordinates and make a short list of replacement loops
					// according to the RMS criteria
					for(  std::vector < core::Size >::const_iterator it = leap_index_bucket.begin();
							it != leap_index_bucket.end();
							++it ){

						// Get the actual strucure index (not just the bin index)
						core::Size retrieve_index = (core::Size) (*it);
						LeapIndex cp = hashmap.get_peptide( retrieve_index );

						// Retrieve the actual backbone structure
						BackboneSegment new_bs;
						library_->backbone_database().get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );

						// Check the values against against any RMS limitations
						core::Real BBrms = get_rmsd( pose_bs, new_bs );
						if( ( BBrms > min_bbrms_) && ( BBrms < max_bbrms_ ) ){
							// Next part is ported from score_fragment.f so these are magical cutoffs
							const std::vector< core::Real > phi = new_bs.phi();
							const std::vector< core::Real > psi = new_bs.psi();
							// Check phi/psi angles against the sequence
							// Pose counts residues starting from one, so offset that
							bool offlimits=false;
							if( get_filter_by_phipsi() ){
								for( core::Size bs_position = 0; bs_position < phi.size() ; ++bs_position ){
									int sequence_position = ir - 1 + bs_position;
									offlimits = true;
									// Proline
									if( sequence[sequence_position] == 'P' ) {
										 if( phi[bs_position] < -103 || phi[bs_position] > -33 ) break;
									}
									// Beta branched residues
									if( sequence[sequence_position] == 'I' || sequence[sequence_position] == 'V' || sequence[sequence_position] == 'T' ) {
											if( phi[bs_position] > -40 ) break; 
									}
									// Non glycine residues are confined to only part of the positive phi region
									// populated by glycine residues
									if( sequence[sequence_position] != 'I' || sequence[sequence_position] != 'V' || sequence[sequence_position] != 'T' || 
											sequence[sequence_position] != 'P'||  sequence[sequence_position] != 'G' ) {
											if( phi[bs_position] > 70 ) break; 
									}
									if( sequence[sequence_position] != 'G' ) {
										 if( psi[bs_position] < -75 && psi[bs_position] > -170 ) break;
									}
									// Residues other than glycine preceding prolines are quite restricted
//									if( sequence[sequence_position] == 'P' ) {
//										 if( phi[bs_position] < 40 && sequence[sequence_position] != 'G' ) {
//													if( psi[bs_position] > -25 || phi[bs_position] < -90 ) break;
//										 }
//									}
									offlimits = false;
								}
							}
							
							impossible_torsion_candidates++;
							if( offlimits ) {
									impossible_torsion_rejects++;
							} else {
									filter_leap_index_bucket.push_back( *it );
							}
							fragments_seen++;
							if( fragments_seen > sw_nfrags ) break; // continue with however many are in the bucket now, and break at end
						}
					}

					std::random_shuffle( filter_leap_index_bucket.begin(), filter_leap_index_bucket.end());

					// Now create models and check rms after insertion
					for(  std::vector < core::Size >::const_iterator it = filter_leap_index_bucket.begin();
							it != filter_leap_index_bucket.end();
							++it ){

						clock_t starttime = clock();

						core::Size retrieve_index = *it;
						LeapIndex cp = hashmap.get_peptide( retrieve_index );

						BackboneSegment new_bs;
						library_->backbone_database().get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );

						core::pose::Pose newpose( start_pose );
						//transfer_phi_psi( start_pose, newpose );			//fpd necessary??

						core::Real final_rms = inserter_->make_local_bb_change( newpose, original_pose, new_bs, ir );

						bool isok = false;
						if ( ( final_rms < max_rms_ ) && ( final_rms > min_rms_) ){

							core::pose::Pose mynewpose( start_pose );

							transfer_phi_psi( newpose, mynewpose );
							transfer_jumps( newpose, mynewpose );

							core::io::silent::SilentStructOP new_struct = nonideal_ ?
								core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary") :
								core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
							new_struct->fill_struct( mynewpose );    // make the silent struct from the copy pose
							new_struct->energies_from_pose( newpose ); // take energies from the modified pose, not the copy pose
							//TR << "SAMPLER: " << new_struct->get_energy("censcore") << std::endl;
							// Add donor history for this round of loophash only

							// Assume extra data is loade, because we need it!
							BBData bb;
							BBExtraData bbextra;
							library_->backbone_database().get_protein( cp.index, bb );
							
							std::string donorhistory = new_struct->get_comment("donorhistory");
							if( library_->backbone_database().extra_size() <= bb.extra_key ){ 
								std::cerr << "ERROR: No extra data ?: " << library_->backbone_database().extra_size() << " < " << bb.extra_key << std::endl;
								
								donorhistory = donorhistory
												+       utility::to_string( loop_size )
												+ "/" + utility::to_string( library_->loopdb_range().first + cp.index )
												+ "/" + utility::to_string( cp.offset/3 )
												+ "/" + utility::to_string( ir-1 ) + ";";
							}else{

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
							models_seen++;
							models_seen_this_rad++;
							if( models_seen > sw_nmodels ) break;
							if( models_seen_this_rad > sw_nmodels_per_rad ) break;
							isok = true;
						}

						//if ( lib_structs.size() > 2  ) return;

						clock_t endtime = clock();

						TR.Debug << "Clocks: " << endtime - starttime << "  " << final_rms << (isok ? " OK" : " Reject") << std::endl;

				}
				// To break out of the outer for loop when these conditions are met
				if( models_seen > sw_nmodels ) break;
				if( fragments_seen > sw_nfrags) break;
      }
    }
	}


		long endtime = time(NULL);


		TR.Info << "LHS: " << start_res << "-" << stop_res << ":  " <<lib_structs.size() << "struc " << endtime - starttime << " secs, rej: " << impossible_torsion_rejects << "/" << impossible_torsion_candidates << std::endl;

		for( std::vector< core::io::silent::SilentStructOP >::iterator it=lib_structs.begin();
				it != lib_structs.end(); ++it ){
			TR.Debug << "SAMPLER2" << (*it)->get_energy("censcore") << std::endl;
		}



}

// closes gaps.  focuses on decreasing variability, not for use with loophash mpi
// very very similar to build_structures(), but too lazy to make separate functions to reduce redundancy
 void
 LoopHashSampler::close_gaps(
		const core::pose::Pose& start_pose,
    std::vector< core::pose::Pose> &lib_structs,
		core::Size loop_size,
		core::Size round
	)
  {
    using namespace core;
    using namespace core::pose;
    using namespace conformation;
    using namespace kinematics;
    using namespace numeric::geometry::hashing;
    using namespace optimization;
    using namespace id;

		runtime_assert( library_ );
		TR.Trace << "Testing libstructs: " << "NStruct: " << lib_structs.size() << std::endl;


		long starttime = time(NULL);


    core::pose::Pose original_pose = start_pose;
    core::pose::Pose edit_pose = start_pose;


	 	TR.Debug  << "Setting options.." << std::endl;
	  //core::optimization::MinimizerOptions options( "dfpmin", 0.2, true , false );
	  //core::optimization::MinimizerOptions options2( "dfpmin", 0.02,true , false );
	  core::optimization::MinimizerOptions options( "lbfgs_armijo", 0.2, true , false );
	  core::optimization::MinimizerOptions options2( "lbfgs_armijo", 0.02,true , false );

    kinematics::MoveMap final_mm;
    final_mm.set_bb(true);
   // setup movemap & minimisation

  //				for ( Size ii=ir; ii<= jr; ++ii ) {
  //					final_mm.set_bb( ii, true );
  //					if ( newpose->residue(ii).aa() == chemical::aa_pro ) final_mm.set( TorsionID( phi_torsion, BB, ii ), false );
  //				}

    Size nres = start_pose.total_residue();
    Size ir, jr;
    //Size newpep_index = 0;

    //core::Size backbone_offset;
    //bbdb_.add_pose( pose, backbone_offset );


    int runcount=0;
    runcount++;

    core::Size start_res = start_res_;
    core::Size stop_res = stop_res_;

    // figure out start and stop residues
    if ( stop_res == 0 ) stop_res = nres;  // to do the whole protein just set stop_res to 0
    start_res = std::max( start_res, (core::Size)2 ); // dont start before 2  - WHY ?
		if( start_res > stop_res ) stop_res = start_res;

    for( ir = start_res; ir <= stop_res; ir ++ ){
		
			jr = ir + loop_size;
			if ( ir > nres ) continue;
			if ( jr > nres ) continue;

			// get the rigid body transform for the current segment
			BackboneSegment pose_bs;
			pose_bs.read_from_pose( start_pose, ir, loop_size );
			Real6 loop_transform;
			if(!get_rt_over_leap( original_pose, ir, jr, loop_transform )) continue;

			LoopHashMap &hashmap = library_->gethash( loop_size );


			core::Size sw_nmodels = 10 ;
			core::Size sw_nfrags = 100;
			core::Size sw_nmodels_per_rad = 5;

			/* dont change these based on avg_sw for now
			set_min_rms(  avg_sw/10   );
			set_max_rms(  avg_sw/8  );
			set_min_bbrms( avg_sw/20  );  
			set_max_bbrms(  avg_sw/12  );
			*/


			// we want sw_nmodels models no matter what
			// but if there is a brrms constraint, we might never reach x models
			// some breakpoint, like x bins checked or x frags checked
			// or radius check
			core::Size fragments_seen = 0;
			core::Size models_seen = 0;

		
			for( Size radius = 0; radius <= 1; radius++ ) {
				core::Size models_seen_this_rad = 0;
				std::vector < core::Size > leap_index_bucket;
				std::vector < core::Size > filter_leap_index_bucket;
				hashmap.radial_lookup( radius, loop_transform, leap_index_bucket );
				//TR << "radius, lookup_size = " << radius << ", " << leap_index_bucket.size() << std::endl;

				if( leap_index_bucket.size() == 0) {
					continue;
				}

				// Now for every hit, get the internal coordinates and make a short list of replacement loops
				// according to the RMS criteria
				for(  std::vector < core::Size >::const_iterator it = leap_index_bucket.begin();
						it != leap_index_bucket.end();
						++it ){

					// Get the actual strucure index (not just the bin index)
					core::Size retrieve_index = (core::Size) (*it);
					LeapIndex cp = hashmap.get_peptide( retrieve_index );

					// Retrieve the actual backbone structure
					BackboneSegment new_bs;
					library_->backbone_database().get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );

					// Check the values against against any RMS limitations
					core::Real BBrms = get_rmsd( pose_bs, new_bs );
					if( ( BBrms > min_bbrms_) && ( BBrms < max_bbrms_ ) ){
						filter_leap_index_bucket.push_back( *it );
						fragments_seen++;
						if( fragments_seen > sw_nfrags ) break; // continue with however many are in the bucket now, and break at end
					}
				}

				std::random_shuffle( filter_leap_index_bucket.begin(), filter_leap_index_bucket.end());

				// Now create models and check rms after insertion
				for(  std::vector < core::Size >::const_iterator it = filter_leap_index_bucket.begin();
						it != filter_leap_index_bucket.end();
						++it ){

					clock_t starttime = clock();

					core::Size retrieve_index = *it;
					LeapIndex cp = hashmap.get_peptide( retrieve_index );

					BackboneSegment new_bs;
					library_->backbone_database().get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );

					core::pose::Pose newpose( start_pose );
					//transfer_phi_psi( start_pose, newpose );			//fpd necessary??

					core::Real final_rms = inserter_->make_local_bb_change_close_gaps( newpose, original_pose, new_bs, ir );

					bool isok = false;
					if ( ( final_rms < max_rms_ ) && ( final_rms > min_rms_) ){
						lib_structs.push_back( newpose );
						models_seen++;
						models_seen_this_rad++;
						if( models_seen > sw_nmodels ) break;
						if( models_seen_this_rad > sw_nmodels_per_rad ) break;
						isok = true;
					}

					//if ( lib_structs.size() > 2  ) return;

					clock_t endtime = clock();

					TR.Debug << "Clocks: " << endtime - starttime << "  " << final_rms << (isok ? " OK" : " Reject") << std::endl;

			}
			// To break out of the outer for loop when these conditions are met
			if( models_seen > sw_nmodels ) break;
			if( fragments_seen > sw_nfrags) break;
		}
	}
}


} // namespace loops
} // namespace protocols




