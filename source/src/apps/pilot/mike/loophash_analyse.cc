// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Mike Tyka
/// @brief


// libRosetta headers
#include <protocols/jd2/JobDistributor.hh>

#include <boost/cstdint.hpp>
#include <boost/unordered_map.hpp>
#include <core/fragment/picking/VallChunk.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/fragment/picking/VallProvider.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.hh>
#include <basic/options/option.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <core/scoring/rms_util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <devel/init.hh>
#include <numeric/HomogeneousTransform.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/match/Hit.fwd.hh>
#include <protocols/match/Hit.hh>
#include <numeric/geometry/hashing/SixDHasher.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <protocols/loophash/BackboneDB.hh>

// C++ headers
//#include <cstdlib>

#include <iostream>
#include <string>
#include <cstdio>
#include <algorithm>

// option key includes
#include <basic/options/keys/broker.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/util/SwitchResidueTypeSet.hh>


static THREAD_LOCAL basic::Tracer TR( "main" );

using namespace protocols::moves;
using namespace core::scoring;
using namespace core;
using namespace core::pose;
using namespace conformation;
using namespace kinematics;
using namespace protocols::match;
using namespace core::fragment::picking;


using namespace protocols::loophash;


class LoopHash_Analyze;
typedef utility::pointer::owning_ptr< LoopHash_Analyze > LoopHash_AnalyzeOP;
typedef utility::pointer::owning_ptr< LoopHash_Analyze const > LoopHash_AnalyzeCOP;

class LoopHash_Analyze: public protocols::moves::Mover {
public:

  LoopHash_Analyze(
    LoopHashLibraryOP library
  ):
   library_(library)

  {
  }

	virtual void apply( core::pose::Pose& pose );

  virtual protocols::moves::MoverOP clone() const {
		return new LoopHash_Analyze( *this );
	}


	virtual std::string get_name() const {
		return "LoopHash_Analyze";
	}

	virtual	protocols::moves::MoverOP	fresh_instance() const {
		return new LoopHash_Analyze( library_ );
	}

private:
  LoopHashLibraryOP library_;

};

void
LoopHash_Analyze::apply( core::pose::Pose& pose )
{
  if( !library_ ) return;

  using namespace basic::options;
  using namespace basic::options::OptionKeys;

	core::pose::Pose native;
	core::import_pose::pose_from_file( native, option[ basic::options::OptionKeys::in::file::native ] , core::import_pose::PDB_file);

	core::Size nres = pose.size();

	core::Size start_res = 1 ;
	core::Size stop_res = nres;
	std::string sequence = pose.sequence();

	for( core::Size ir = start_res; ir <= stop_res; ir ++ ){
		for( core::Size k = 0; k < library_->hash_sizes().size(); k ++ ){
			std::cout << "-----------------------------------------------------" << std::endl;

			core::Size loop_size = library_->hash_sizes()[ k ];
			core::Real displace =  loop_size ;
			core::Real displace_angle =  15.0 * loop_size / 6.0 ;

			core::Size jr = ir + loop_size;
			if ( ir > nres ) continue;
			if ( jr > nres ) continue;

			// get the backbone segment from the pose
			BackboneSegment pose_bs;
			pose_bs.read_from_pose( pose, ir, loop_size );
			Real6 loop_transform;
			if(!get_rt_over_leap( pose, ir, jr, loop_transform )) continue;

			BackboneSegment native_bs;
			native_bs.read_from_pose( native, ir, loop_size );

      LoopHashMap &hashmap = library_->gethash( loop_size );


			// Get the fragment bucket

			for( core::Real expand = 0; expand <= 1; expand += 1.0 ){
      std::vector < core::Size > leap_index_bucket;
      std::vector < core::Size > leap_index_bucket_unique;

			core::Real x=0;
    	core::Real y=0;
    	core::Real z=0;
			core::Real a=0;
    	core::Real b=0;
    	core::Real c=0;

			for( x = -displace*expand; x <= displace*expand; x += displace ){
			for( y = -displace*expand; y <= displace*expand; y += displace ){
			for( z = -displace*expand; z <= displace*expand; z += displace ){

				for( a = -displace_angle*expand; a <= displace_angle*expand; a += displace_angle ){
				for( b = -displace_angle*expand; b <= displace_angle*expand; b += displace_angle ){
				for( c = -displace_angle*expand; c <= displace_angle*expand; c += displace_angle ){


						Real6 loop_transform_disp = loop_transform;
						loop_transform_disp[1] += x;
						loop_transform_disp[2] += y;
						loop_transform_disp[3] += z;
						loop_transform_disp[4] += a;
						loop_transform_disp[5] += b;
						loop_transform_disp[6] += c;

						hashmap.lookup( loop_transform_disp, leap_index_bucket );
				}
				}
				}

			}
		  }
			}

			std::sort(   leap_index_bucket.begin(), leap_index_bucket.end() );

			std::vector < core::Size >::const_iterator last;
      for(  std::vector < core::Size >::const_iterator it = leap_index_bucket.begin();
            it != leap_index_bucket.end(); ++it )
			{
				if( it !=  leap_index_bucket.begin() ){
					if( *last != *it ) leap_index_bucket_unique.push_back( *it );
				}
			  last = it;
			}


			{
			core::Size retrieve_index;
			core::Size count_full = 0;
			core::Size count_rama_filter = 0;
			core::Real lowest_rms = 100000.0;
			core::Size lowest_rms_index = 0;
      for(  std::vector < core::Size >::const_iterator it = leap_index_bucket_unique.begin();
            it != leap_index_bucket_unique.end(); ++it ){

				// Get the actual strucure index (not just the bin index)
				core::Size retrieve_index = (core::Size) (*it);
				LeapIndex cp = hashmap.get_peptide( retrieve_index );

				// Retrieve the actual backbone structure
				BackboneSegment new_bs;
				library_->backbone_database().get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );

				core::Real BBrms = get_rmsd( native_bs, new_bs );
				count_full ++;

				if( BBrms < lowest_rms ){
					lowest_rms = BBrms;
					lowest_rms_index = retrieve_index;
				}
      }

			LeapIndex cp = hashmap.get_peptide( lowest_rms_index );
			// Retrieve the actual backbone structure
			BackboneSegment new_bs;
			library_->backbone_database().get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );

			LocalInserter_SimpleMinOP simple_inserter( new LocalInserter_SimpleMin() );
			core::pose::Pose newpose( pose );
			core::util::switch_to_residue_type_set( newpose, core::chemical::CENTROID_t );
			transfer_phi_psi( pose, newpose );
			core::Real final_rms = simple_inserter->make_local_bb_change( newpose, pose, new_bs, ir );

			newpose.dump_pdb( "repseg_" +  utility::to_string( loop_size ) + "_" + utility::to_string( ir ) + ".pdb" );

			core::Real total_rms = core::scoring::native_CA_rmsd( native, newpose );

			std::cout << "LAbin " << loop_size  << "  " << expand << " " << ir << " " << count_full << " " << lowest_rms_index << " " << lowest_rms << " "  << count_rama_filter << "  " << final_rms << "  " <<  total_rms << std::endl;

			}


			}  // expansion loop


			// go through ALL the fragments and find the best.
			{
			core::Size retrieve_index;
			core::Size count_full = 0;
			core::Size count_rama_filter = 0;
			core::Real lowest_rms = 100000.0;
			core::Size lowest_rms_index = 0;
			for( retrieve_index = 0; retrieve_index < hashmap.n_loops(); retrieve_index++ ){

				LeapIndex cp = hashmap.get_peptide( retrieve_index );

				// Retrieve the actual backbone structure
				BackboneSegment new_bs;
				library_->backbone_database().get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );

				// Check the values against against any RMS limitations imposed by the caller
				core::Real BBrms = get_rmsd( native_bs, new_bs );
				//std::cout << loop_size << "_" << retrieve_index << " " << BBrms << std::endl;

				// Next part is ported from score_fragment.f so these are magical cutoffs
				const std::vector< core::Real > phi = new_bs.phi();
				const std::vector< core::Real > psi = new_bs.psi();
				bool offlimits = false; //true;


				// Check phi/psi angles against the sequence
				// Pose counts residues starting from one, so offset that
//				for( core::Size g = 0; g < phi.size(); g++ ){
//						core::Size m = g + ir;
//						// Proline
//						if( sequence[m] == 'P' ) {
//							 if( phi[g] < -103 || phi[g] > -33 ) break;
//						}
//						// Beta branched residues
//						if( sequence[m] == 'I' || sequence[m] == 'V' || sequence[m] == 'T' ) {
//							if( phi[g] > -40 ) break;
//						}
//						// Non glycine residues are confined to only part of the positive phi region
//						// populated by glycine residues
//						if( sequence[m] != 'I' || sequence[m] != 'V' || sequence[m] != 'T' ||
//								sequence[m] != 'P'||  sequence[m] != 'G' ) {
//								if( phi[g] > 70 ) break;
//						}
//						if( sequence[m] != 'G' ) {
//							 if( psi[g] < -75 && psi[g] > -170 ) break;
//						}
//						// Residues other than glycine preceding prolines are quite restricted
////						if( sequence[m] == 'P' ) {
////							 if( m > 1 ){
////								 if( phi[m-1] < 40 && sequence[m-1] != 'G' ) {
////											if( psi[m-1] > -25 || phi[m-1] < -90 ) break;
////								 }
////							 }
////						}
//						offlimits = false;
//				}
				count_full++;
				if( offlimits ) {
					count_rama_filter++;
					continue;
				}


				if( BBrms < lowest_rms ){
					lowest_rms = BBrms;
					lowest_rms_index = retrieve_index;
				}
			}


			LeapIndex cp = hashmap.get_peptide( lowest_rms_index );
			// Retrieve the actual backbone structure
			BackboneSegment new_bs;
			library_->backbone_database().get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );

			LocalInserter_SimpleMinOP simple_inserter( new LocalInserter_SimpleMin() );
			core::pose::Pose newpose( pose );
			core::util::switch_to_residue_type_set( newpose, core::chemical::CENTROID_t );
			transfer_phi_psi( pose, newpose );
			core::Real final_rms = simple_inserter->make_local_bb_change( newpose, pose, new_bs, ir );

			newpose.dump_pdb( "repseg_" +  utility::to_string( loop_size ) + "_" + utility::to_string( ir ) + ".pdb" );

			core::Real total_rms = core::scoring::native_CA_rmsd( native, newpose );


			// print statistics
			std::cout << "LAall " << loop_size  << "  " << ir << " all " << count_full << " " << lowest_rms_index << " " << lowest_rms << " " << count_rama_filter << "  " << final_rms << "  " <<  total_rms << std::endl;


			// Calculate the RT of the original fragment by inserting it raw and remeasuring.
			core::pose::Pose fragpose( pose );
			core::util::switch_to_residue_type_set( fragpose, core::chemical::CENTROID_t );
			transfer_phi_psi( pose, fragpose );
			new_bs.apply_to_pose( fragpose, ir );
			Real6 loop_transform_frag;
			if(!get_rt_over_leap( fragpose, ir, jr, loop_transform_frag )) continue;

			std::cout << "GD:   " << loop_size  << "  " << ir << "  "  << loop_transform[1]      << "  " << loop_transform[2]      << "  " << loop_transform[3]      << "  "
			                                                << loop_transform[4]      << "  " << loop_transform[5]      << "  " << loop_transform[6]      << std::endl;
      std::cout << "GF:   " << loop_size  << "  " << ir << "  "  << loop_transform_frag[1] << "  " << loop_transform_frag[2] << "  " << loop_transform_frag[3] << "  "
			                                                << loop_transform_frag[4] << "  " << loop_transform_frag[5] << "  " << loop_transform_frag[6] << std::endl;
			std::cout << "GDGF: " << loop_size  << "  " << ir << "  "
												 << (loop_transform[1] - loop_transform_frag[1]) / displace << " "
												 << (loop_transform[2] - loop_transform_frag[2]) / displace << " "
												 << (loop_transform[3] - loop_transform_frag[3]) / displace << " "
			                   << (loop_transform[4] - loop_transform_frag[4]) / displace_angle << " "
												 << (loop_transform[5] - loop_transform_frag[5]) / displace_angle << " "
												 << (loop_transform[6] - loop_transform_frag[6]) / displace_angle << " "
												 << std::endl;

			}


		}
	}


}


int
main( int argc, char * argv [] )
{
    try {
	using namespace protocols;
	using namespace protocols::jd2;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;


	// initialize core
	devel::init(argc, argv);

#ifdef USEMPI
	int mpi_rank_, mpi_npes_;
	MPI_Comm_rank( MPI_COMM_WORLD, ( int* )( &mpi_rank_ ) );
	MPI_Comm_size( MPI_COMM_WORLD, ( int* )( &mpi_npes_ ) );

	// unless you are rank one - go into infinite sleep loop
	if( mpi_rank_ != 0 ){
		TR << "NOT RANK 0: Sleeping .. " << std::endl;
		while(true){
			sleep( 10 );
		}
	}

#endif


	utility::vector1 < core::Size > loop_sizes = option[lh::loopsizes]();
	LoopHashLibraryOP loop_hash_library = new LoopHashLibrary( loop_sizes );

  loop_hash_library->load_db();
  LoopHash_AnalyzeOP lh_analyze = new LoopHash_Analyze( loop_hash_library );

  // Normal mode with external loophash library
  try{
    protocols::jd2::JobDistributor::get_instance()->go( lh_analyze );
  } catch ( utility::excn::EXCN_Base& excn ) {
    std::cerr << "Exception: " << std::endl;
    excn.show( std::cerr );
    std::cout << "Exception: " << std::endl;
    excn.show( std::cout ); //so its also seen in a >LOG file
  }

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}


