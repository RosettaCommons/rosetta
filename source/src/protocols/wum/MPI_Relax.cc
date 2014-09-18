// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopHashMap.cc
/// @brief
/// @author Mike Tyka


// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

#include <protocols/wum/MPI_Relax.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/util.hh>
#include <protocols/wum/SilentStructStore.hh>



namespace protocols {
namespace wum {


static thread_local basic::Tracer TR( "MPI_Relax" );

void
MPI_Relax::init_master(){
	// Open input stream
  TR << "Opening input streams! " << std::endl;
	input_ = core::import_pose::pose_stream::streams_from_cmd_line();

}


bool
MPI_Relax::fill_outbound_queue()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	TR << "Reading in structures..." << std::endl;

	core::chemical::ResidueTypeSetCAP rsd_set;
	if ( option[ in::file::fullatom ]() ) {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	} else {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	}

	core::Size count = 0;
	while( input_.has_another_pose() && ( outbound().size() < max_out_queue_size_ ) ) {
		TR << "Reading in pose: " << count << std::endl;
		core::pose::Pose pose;
		input_.fill_pose( pose, *rsd_set );

		TR << "Adding pose: " << count << std::endl;

		WorkUnitBaseCOP base_cop = work_unit_list().get_work_unit("relax");
		WorkUnit_MoverWrapperOP new_wu = dynamic_cast<  WorkUnit_MoverWrapper* >( &(*(base_cop->clone())));
		if( !new_wu ){
			TR << "ERROR ERROR ERROR" << std::endl;
			continue;
		}
		new_wu->decoys().add( pose );
		outbound().add( new_wu );
		count ++;
	}
	TR << "Added " << count << " workunits to outbound queue" << std::endl;

	if( count == 0) return false;

	return true;
}


void MPI_Relax::register_work_units(){
	TR << "Register movers..." << std::endl;
	protocols::moves::MoverOP relax_pose = protocols::relax::generate_relax_from_cmd();
	WorkUnit_MoverWrapperOP wu_relax = new WorkUnit_MoverWrapper( relax_pose );
	work_unit_list().register_work_unit( "relax", wu_relax );
}


void
MPI_Relax::init_slave(){
	TR << "MPI_Relax slave. " << std::endl;
}

void
MPI_Relax::process_outbound_wus_master(){
	// just make sure the queue stays full
	fill_outbound_queue();
}


void
MPI_Relax::process_inbound_wus_master(){
 	using namespace basic::options;
  using namespace basic::options::OptionKeys;

	while( inbound().size() > 0 )
	{
		WorkUnitBaseOP  next_wu =  inbound().pop_next();

		WorkUnit_SilentStructStore* structure_wu = dynamic_cast<  WorkUnit_SilentStructStore * > ( (WorkUnitBase*) (&(*next_wu)) );

		if ( structure_wu == NULL ){
			TR << "Cannot save structural data for WU: " << std::endl;
			next_wu->print( std::cout );
		} else {
			TR << "Saving decoy store.. " << std::endl;
			SilentStructStore &decoys = structure_wu->decoys();
			if( decoys.size() == 0 ){
				TR << "ERROR: WU did not contain any structures. " << std::endl;
			}else{
				core::io::silent::SilentFileData sfd;
				std::string filename = option[ OptionKeys::out::file::silent ]();

				for( SilentStructStore::const_iterator it = decoys.begin();
						 it != decoys.end();
						 ++it )
				{
					sfd.write_silent_struct( (*(*it)), filename );
				}
			}

		}


	}



}




}
}

