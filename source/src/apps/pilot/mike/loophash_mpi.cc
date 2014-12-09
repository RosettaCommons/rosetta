// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Mike Tyka
/// @brief




#include <core/io/pdb/pose_io.hh>
#include <core/pose/util.hh>
#include <protocols/wum/WorkUnitList.hh>
#include <protocols/wum/WorkUnitManager.hh>
#include <protocols/wum/MPI_WorkUnitManager_Slave.hh>
#include <protocols/relax/WorkUnit_BatchRelax.hh>
#include <protocols/loophash/MPI_LoopHashRefine.hh>
#include <protocols/loophash/MPI_LoopHashRefine_Master.hh>
#include <protocols/loophash/MPI_LoopHashRefine_Emperor.hh>
#include <protocols/loophash/WorkUnit_LoopHash.hh>
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


class MPI_LoopHash_Launcher {
 public:


  void run(){
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols::wum;
		using namespace protocols::relax;
		using namespace protocols::loophash;

		core::pose::PoseOP native_pose_;
		if ( option[ in::file::native ].user() ) {
			native_pose_ = new core::pose::Pose();
			core::import_pose::pose_from_pdb( *native_pose_, option[ in::file::native ]() );
			core::pose::set_ss_from_phipsi( *native_pose_ );
			core::util::switch_to_residue_type_set( *native_pose_, core::chemical::CENTROID);
		}

		WorkUnitList wulist;

		// Create a workunit for waiting
		WorkUnit_WaitOP wait_wu = new WorkUnit_Wait( 2 );
		wulist.register_work_unit( "waitwu", wait_wu );

		// Create a workunit for tranfer of library structures
		WorkUnit_SilentStructStoreOP resultpack = new WorkUnit_SilentStructStore( );
		wulist.register_work_unit( "resultpack", resultpack );

		// Create a workunit for requesting a new library structure
		WorkUnit_SilentStructStoreOP getnewstruct = new WorkUnit_SilentStructStore( );
		wulist.register_work_unit( "getnewstruct", getnewstruct );

		// Create a workunit for loop hashing
		WorkUnit_LoopHashOP loophasher = new WorkUnit_LoopHash( );
		loophasher->init_from_cmd( (core::Size) mpi_rank() );
		wulist.register_work_unit( "loophasher", loophasher );

		// Create a workunit for batch relaxing
		WorkUnit_BatchRelaxOP batchrelax = new WorkUnit_BatchRelax_and_PostRescore( );
		batchrelax->set_native_pose( native_pose_ );
		wulist.register_work_unit( "batchrelax", batchrelax );






		WorkUnitManagerOP wu_manager;

		// Now define the structure of the algorithm, i.e. assign roles to individual nodes

		const core::Size emperor_node = 0;
		core::Size n_masters = 1;

		if( option[ OptionKeys::wum::n_masters ].user() ){
			n_masters = option[ OptionKeys::wum::n_masters ]();
		}else{
			n_masters = mpi_npes() / option[ OptionKeys::wum::n_slaves_per_master ]();
		}


		if ( mpi_rank() == 0 ){
			wu_manager = new MPI_LoopHashRefine_Emperor( );
		}
		else if ( (int(mpi_rank()) > int(0)) && int(mpi_rank()) <= int(n_masters) )
		{
			core::Size master_rank =  mpi_rank() - 1; // master rank, just a serial number identifying the master
			wu_manager = new MPI_LoopHashRefine_Master( emperor_node , master_rank );
		}
		else
		{
			core::Size slave_master = (mpi_rank() %n_masters)+1;  // which master does this slave belong to ?
			wu_manager = new MPI_WorkUnitManager_Slave( slave_master );
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

	MPI_LoopHash_Launcher launch_mpi_loophash;
	launch_mpi_loophash.run();

	#ifdef USEMPI
		MPI_Barrier( MPI_COMM_WORLD );
		MPI_Finalize();
	#endif
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}




