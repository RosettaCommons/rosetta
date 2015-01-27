// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/wum/MPI_WorkUnitManager_Slave.hh
/// @brief
/// @author Mike Tyka

#ifndef INCLUDED_protocols_wum_MPI_WorkUnitManager_Slave_hh
#define INCLUDED_protocols_wum_MPI_WorkUnitManager_Slave_hh

#ifdef USEMPI
#include <mpi.h> //keep this first
#else
#define MPI_ANY_SOURCE 0
#endif

#include <protocols/wum/WorkUnitManager.hh>
#include <protocols/wum/MPI_WorkUnitManager.hh>
//  you cannot #include yourself #include <protocols/wum/MPI_WorkUnitManager_Slave.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <string>
#include <vector>

#include <utility/vector1.hh>





namespace protocols {
namespace wum {


class MPI_WorkUnitManager_Slave: public MPI_WorkUnitManager {
  public:
    MPI_WorkUnitManager_Slave( core::Size my_master );

    virtual ~MPI_WorkUnitManager_Slave(){}

		virtual void go();

	protected:
		void init(){};

		virtual void process_inbound_wus();

		virtual void process_outbound_wus();

		/// @brief Slave: call a master to ask for more work
		virtual void request_new_jobs();

		core::Size get_my_master(){ return my_master_; };

	private:
		const core::Size my_master_;
  	bool terminate_;
};


}
}

#endif

