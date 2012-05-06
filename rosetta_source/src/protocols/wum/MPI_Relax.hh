// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/wum/MPI_Relax.hh
/// @brief
/// @author Mike Tyka

#ifndef INCLUDED_protocols_wum_MPI_Relax_hh
#define INCLUDED_protocols_wum_MPI_Relax_hh

#include <protocols/wum/MPI_WorkUnitManager.hh>
#include <protocols/wum/WorkUnitManager.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/moves/Mover.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <string>
#include <vector>

namespace protocols {
namespace wum {



class MPI_Relax: public MPI_WorkUnitManager {
  public:
    MPI_Relax():
			MPI_WorkUnitManager(),
			max_out_queue_size_(100)
		{
			register_work_units();
		};

    virtual ~MPI_Relax(){}

		virtual void register_work_units();
	protected:
		void virtual init_master();

		void virtual init_slave();

		void virtual process_inbound_wus_master();

		void virtual process_outbound_wus_master();

		bool fill_outbound_queue();

	private:

		core::import_pose::pose_stream::MetaPoseInputStream input_;

		core::Size max_out_queue_size_;
};



}
}

#endif

