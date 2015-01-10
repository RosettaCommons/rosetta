// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

#include <core/io/pdb/pose_io.hh>
#include <core/pose/util.hh>
#include <protocols/wum/WorkUnitList.hh>
#include <protocols/wum/WorkUnitManager.hh>
#include <protocols/wum/MPI_WorkUnitManager_Slave.hh>
#include <protocols/relax/WorkUnit_BatchRelax.hh>
#include <protocols/mpi_refinement/MPI_Refinement.hh>
#include <protocols/mpi_refinement/MPI_Refine_Master.hh>
#include <protocols/mpi_refinement/MPI_Refine_Emperor.hh>
#include <protocols/mpi_refinement/WorkUnit_Aggressive.hh>
#include <protocols/mpi_refinement/WorkUnit_Loop.hh>
#include <protocols/mpi_refinement/WorkUnit_Relax.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

#include <cstdio>

#include <devel/init.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/excn/Exceptions.hh>


class MPI_Refinement_Launcher {
 public:

  void run(){

    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    using namespace protocols::wum;
    using namespace protocols::relax;
    using namespace protocols::mpi_refinement;
    
    core::pose::PoseOP native_pose_;
    if ( option[ in::file::native ].user() ) {
      native_pose_ = core::pose::PoseOP( new core::pose::Pose() );
      core::import_pose::pose_from_pdb( *native_pose_, option[ in::file::native ]() );
      core::pose::set_ss_from_phipsi( *native_pose_ );
      core::util::switch_to_residue_type_set( *native_pose_, core::chemical::CENTROID);
    }

    WorkUnitList wulist;

    // Create a workunit for waiting
    WorkUnit_WaitOP wait_wu( new WorkUnit_Wait( 1 ) );
    wulist.register_work_unit( "waitwu", wait_wu );

    // Create a workunit for tranfer of library structures
    WorkUnit_SilentStructStoreOP resultstore( new WorkUnit_SilentStructStore( ) );
    wulist.register_work_unit( "resultstore", resultstore );
    wulist.register_work_unit( "resultfeedback", resultstore->clone() );

    // Create a workunit for tranfer of termination signal
    WorkUnit_SilentStructStoreOP terminate( new WorkUnit_SilentStructStore( ) );
    wulist.register_work_unit( "terminate", terminate );
    wulist.register_work_unit( "terminated", terminate->clone() );

    // Create a workunit for requesting a new library structure
    WorkUnit_SilentStructStoreOP getnewstruct( new WorkUnit_SilentStructStore( ) );
    wulist.register_work_unit( "getnewstruct", getnewstruct );

    // Create a workunit for movers: (is there a smarter way than this?)
		// DONT USE CAPITAL for WU names to make it less confusing
		// also keep sync b/w WU name in code & schedule cmd file
    // LoopHash
		/*
    WorkUnit_SamplerOP sampler1 = new WorkUnit_LoopHash( );
    sampler1->init_from_cmd( (core::Size) mpi_rank() );
    wulist.register_work_unit( "global_loophasher", sampler1 );
    wulist.register_work_unit( "local_loophasher", sampler1->clone() );
		*/

    // bbGauss
    WorkUnit_SamplerOP sampler2( new WorkUnit_bbGauss( ) );
    sampler2->init_from_cmd( (core::Size) mpi_rank() );
    wulist.register_work_unit( "bbgauss", sampler2 );

    // CombinePose
    WorkUnit_SamplerOP sampler3( new WorkUnit_CombinePose( ) );
    sampler3->init_from_cmd( (core::Size) mpi_rank() );
    wulist.register_work_unit( "combine", sampler3 );

    // CartNormalMode
    WorkUnit_SamplerOP sampler4( new WorkUnit_NormalMode( ) );
    sampler4->init_from_cmd( (core::Size) mpi_rank() );
    wulist.register_work_unit( "nm", sampler4 );
    wulist.register_work_unit( "nmcen", sampler4->clone() );

		// MD
    WorkUnit_SamplerOP sampler5( new WorkUnit_MD( ) );
    sampler5->init_from_cmd( (core::Size) mpi_rank() );
    wulist.register_work_unit( "md", sampler5 );

		// Relax
    WorkUnit_SamplerOP sampler6( new WorkUnit_Relax( ) );
    sampler6->init_from_cmd( (core::Size) mpi_rank() );
    wulist.register_work_unit( "relax", sampler6 );
    wulist.register_work_unit( "rerelax", sampler6->clone() );

    // FragInsert
    WorkUnit_SamplerOP sampler7( new WorkUnit_FragInsert( ) );
    sampler7->init_from_cmd( (core::Size) mpi_rank() );
    wulist.register_work_unit( "fraginsert", sampler7 );

    // KIC
    WorkUnit_SamplerOP sampler8( new WorkUnit_KicCloser( ) );
    sampler8->init_from_cmd( (core::Size) mpi_rank() );
    wulist.register_work_unit( "kiccloser", sampler8 );

    // RamaPert
    WorkUnit_SamplerOP sampler9( new WorkUnit_RamaPerturber( ) );
    sampler9->init_from_cmd( (core::Size) mpi_rank() );
    wulist.register_work_unit( "ramapert", sampler9 );

    // RamaPert
    WorkUnit_SamplerOP sampler10( new WorkUnit_PartialAbinitio( ) );
    sampler10->init_from_cmd( (core::Size) mpi_rank() );
    wulist.register_work_unit( "partialabinitio", sampler10 );

    // batch relaxing
    WorkUnit_BatchRelaxOP batchrelax( new WorkUnit_BatchRelax_and_PostRescore( ) );
    batchrelax->set_native_pose( native_pose_ );
    wulist.register_work_unit( "batchrelax", batchrelax );

    WorkUnitManagerOP wu_manager;

    // Now define the structure of the algorithm, i.e. assign roles to individual nodes

    const core::Size emperor_node = 0;
    core::Size n_masters = 1;

    if( option[ OptionKeys::wum::n_masters ].user() ){
      n_masters = option[ OptionKeys::wum::n_masters ]();
    } else {
      //n_masters = mpi_npes() / option[ OptionKeys::wum::n_slaves_per_master ]();
			std::cerr << "wum::n_masters should be specified!" << std::endl;
			return;
    }

    if ( mpi_rank() == 0 ){
      wu_manager = WorkUnitManagerOP( new MPI_Refine_Emperor( ) );

    } else if ( (int(mpi_rank()) > int(0)) && int(mpi_rank()) <= int(n_masters) ) {
      core::Size master_rank =  mpi_rank() - 1; // master rank, just a serial number identifying the master
      wu_manager = WorkUnitManagerOP( new MPI_Refine_Master( emperor_node , master_rank ) );

    } else {
			core::Size slave_master;

			if( option[ OptionKeys::lh::mpi_master_cpu_weight ].user() ){
				utility::vector1< core::Real > slave_frac =	option[ OptionKeys::lh::mpi_master_cpu_weight ]();
				core::Real fsum( 0.0 );
				for( core::Size i = 1; i <= slave_frac.size(); ++i ) fsum += slave_frac[i];
				runtime_assert( fsum > 0.000001 );
				runtime_assert( slave_frac.size() == n_masters );

				for( core::Size i = 1; i <= slave_frac.size(); ++i ) slave_frac[i] /= fsum;

				int ncores = mpi_npes(); //MPI::COMM_WORLD.Get_size();
				core::Size nslave = (core::Size)(ncores) - n_masters - 1;

				core::Real fslave( (core::Real)(mpi_rank())/nslave );

				fsum = 0.0;
				core::Size i_master( 1 );
				for( i_master = 1; i_master <= n_masters; ++i_master ){
					fsum += slave_frac[i_master];
					if( fslave <= fsum ) break;
				}
				slave_master = i_master;

			} else {
				slave_master = (mpi_rank() %n_masters)+1;  // which master does this slave belong to ?
			}

      wu_manager = WorkUnitManagerOP( new MPI_WorkUnitManager_Slave( slave_master ) );
    }

    // make sure all the necessary work unit have been properly initiated.
    wu_manager->register_work_units( wulist );

    // now launch this node
    wu_manager->go();
  }
};


int
main( int argc, char * argv [] )
{

  try {
    // initialize core
    devel::init(argc, argv);
    
    MPI_Refinement_Launcher launch_mpi;
    launch_mpi.run();
    
#ifdef USEMPI
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
#endif
  } catch ( utility::excn::EXCN_Base const & e ) {
    std::cerr << "caught exception " << e.msg() << std::endl;
  }
  return 0;
}




