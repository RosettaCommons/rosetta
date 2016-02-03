// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/mpi_refinement/Serial_Refine.cc
/// @brief
/// @author Hahnbeom Park

#include <protocols/mpi_refinement/Serial_Refine.hh>
#include <protocols/mpi_refinement/util.hh>
#include <protocols/mpi_refinement/WorkUnit_Aggressive.hh>
#include <protocols/mpi_refinement/WorkUnit_Relax.hh>
#include <protocols/mpi_refinement/WorkUnit_Loop.hh>
#include <protocols/mpi_refinement/StructAvrgMover.hh>
#include <protocols/mpi_refinement/Scheduler.hh>
#include <protocols/mpi_refinement/MultiObjective.hh>
//#include <protocols/mpi_refinement/MPI_Refinement.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/wum/SilentStructStore.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <protocols/relax/AtomCoordinateCstMover.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>

/// ObjexxFCL headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <numeric/random/random.hh>
#include <utility/string_util.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

//Auto Headers
#include <utility/vector1.hh>
#include <numeric/random/random.hh>

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

namespace protocols {
namespace mpi_refinement {

using namespace protocols::wum;

static basic::Tracer TR("protocols.mpi_refinement.SerialRefine");


Serial_Refine::Serial_Refine()
{
	set_defaults();
	init();
}

Serial_Refine::~Serial_Refine(){}

void
Serial_Refine::set_defaults(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
}

void
Serial_Refine::init(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ OptionKeys::lh::mpi_master_schfile ].user() ) {
		TR << "Setup scheduler as given mpi_master_schfile: ";
		TR << option[ OptionKeys::lh::mpi_master_schfile ].user() << std::endl;
		scheduler_.set_random( false );
		scheduler_.prepare_search_stage( 0 );
	} else {
		utility_exit_with_message( "Scheduler needs to be specified." );
	}

	native_given_ = false;
	if ( option[ in::file::native ].user() ) {
		native_given_ = true;
		core::chemical::ResidueTypeSetCOP rsd_set
			= core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		core::import_pose::pose_from_file( native_pose_, *rsd_set, option[ in::file::native ]() , core::import_pose::PDB_file);
	}

	fobj_ = MultiObjectiveOP( new MultiObjective() );
}

void
Serial_Refine::load_structures_from_cmdline_into_library(
	core::pose::Pose const & pose,
	protocols::wum::SilentStructStore &library )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::pose::Pose pose_work( pose );

	// Minimize pose prior to store
	if ( option[ OptionKeys::lh::mpi_packmin_init ]() ) {
		core::scoring::ScoreFunctionOP scorefxn = fobj_->get_scorefxn( 1 ); // this should be standard energy...
		scorefxn->set_weight( core::scoring::coordinate_constraint, 10.0 ); //10.0 will allow about 0.1 ang diff...

		protocols::relax::AtomCoordinateCstMover coord_cst_mover;
		coord_cst_mover.cst_sd( 1.0 );
		coord_cst_mover.apply( pose_work );

		ramp_minpack_pose( pose_work, scorefxn ); // this is cartmin!
	}

	core::pose::set_ss_from_phipsi( pose_work );
	core::io::silent::SilentStructOP ss =
		core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary");
	ss->fill_struct( pose_work );
	ss->add_energy( "state", 0 );     // state: 0 init, 1 perturbed, 2 relaxed

	if ( native_given_ ) add_poseinfo_to_ss( *ss, native_pose_, "" );

	library.add( ss );

	fobj_->add_objective_function_info( library );
	TR << "Added " << library.size() << " structures to library " << std::endl;

}

core::Real
Serial_Refine::apply( core::pose::Pose &pose,
	utility::vector1< core::Size > fixres
)
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::wum;

	long starttime = time(NULL);

	// load structure
	load_structures_from_cmdline_into_library( pose, library_ref_ );

	// set movemap from input fixres definition
	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_bb( true );
	mm->set_chi( true );
	mm->set_jump( true );
	mm->set( core::id::THETA, true );
	mm->set( core::id::D, true );

	// fix mm for fixres
	for ( core::Size ires = 1; ires <= pose.total_residue(); ++ires ) {
		if ( fixres.contains( ires ) ) {
			mm->set_bb ( ires, false );
			mm->set_chi( ires, false );
			for ( core::Size j=1; j<=pose.residue_type(ires).natoms(); ++j ) {
				mm->set( core::id::DOF_ID(core::id::AtomID(j,ires), core::id::THETA ), false );
				mm->set( core::id::DOF_ID(core::id::AtomID(j,ires), core::id::D ), false );
			}
		}
	}

	// Get next queue from scheduler
	MethodParams params = scheduler_.get_params();
	TR << "Running mover: " << params.movertype << ", ngen ";
	TR << library_ref_.size() << std::endl;
	core::Size imovetype( 1 );

	while ( true ) {
		SilentStructStore method_decoys;

		core::Size igen( 0 );
		// Run pert->relax->store
		for ( core::Size irun = 1; irun <= params.nrun; ++irun ) {
			core::Size rerelax_type = params.rerelax_type;

			for ( SilentStructStore::iterator it = library_ref_.begin(),
					end =  library_ref_.end(); it != end; ++it ) {
				TR << "Create on structure: " << (*it)->decoy_tag() << std::endl;

				// perturb
				SilentStructStore pert_structs = perturb( params, *it );
				fobj_->add_objective_function_info( pert_structs );
				for ( SilentStructStore::const_iterator it2 = pert_structs.begin();
						it2 != pert_structs.end(); ++it2, ++igen ) {
					std::stringstream decoytag;
					decoytag << params.movertype << "." << igen;
					(*it2)->set_decoy_tag( decoytag.str() );
				}

				dump_structures( pert_structs, false, "pert" );
				// rerelax
				WorkUnit_SamplerOP relax_wu( new WorkUnit_Relax( rerelax_type, 0, 1, 0.0 ) );
				relax_wu->decoys().add( pert_structs );
				relax_wu->run();

				// store
				SilentStructStore rerelax_structs = relax_wu->decoys();
				//for( SilentStructStore::const_iterator it = rerelax_structs.begin();
				//   it != rerelax_structs.end(); ++it )

				fobj_->add_objective_function_info( rerelax_structs );

				dump_structures( rerelax_structs, false, "rerelax" );
				method_decoys.add( rerelax_structs );
			}
		}

		// shave fraction of structs before storing
		method_decoys.all_add_energy( "samplemethod", imovetype );
		method_decoys.sort_by( "goap" );
		core::Size npop( 0 );
		core::Size n_to_pop = core::Size( method_decoys.size()*params.fshave1 );
		while ( npop < n_to_pop ) {
			method_decoys.store().pop_back();
			npop++;
		}

		// store
		//fobj_->add_objective_function_info( method_decoys );
		library_central_.add( method_decoys );

		scheduler_.proceed();
		params = scheduler_.get_params();
		if ( scheduler_.roundtype().compare("done") == 0 ) break;

		TR << "Running mover: " << params.movertype << ", ngen ";
		TR << library_central_.size() << std::endl;
		imovetype ++;
	}

	// Here methodpicker can come
	utility::vector1< core::Size > methods_picked = scheduler_.methods_picked();
	TR << "methods picked: ";
	for ( core::Size i = 1; i <= methods_picked.size(); ++i ) TR << " " << methods_picked[i];
	TR << std::endl;

	// finally, average
	pose::Pose avrg_pose =
		get_average_structure( library_central_, methods_picked,
		"samplemethod", true );

	long endtime = time(NULL);
	// report time
	core::Real ds = (core::Real)(endtime - starttime);
	core::Size elapsemin = (core::Size)(ds/60.0);
	core::Size elapsesec = (core::Size)(ds - 60.0*elapsemin);
	TR << "Total time: " << elapsemin << " min " << elapsesec << " sec." << std::endl;

	pose = avrg_pose;
	// TODO: fobj_ to accept pose and return a vector of scores
	return 0.0;
}

SilentStructStore
Serial_Refine::perturb( MethodParams const &params,
	const core::io::silent::SilentStructOP &start_struct )
{
	WorkUnit_SamplerOP new_wu;
	std::string movername( params.movertype );

	if ( movername.compare("relax") == 0 ) {
		new_wu = WorkUnit_SamplerOP( new WorkUnit_Relax( params.relax_type, params.score_type, params.nperrun, params.cstw ) );

	} else if ( movername.compare("md") == 0 ) {
		new_wu = WorkUnit_SamplerOP( new WorkUnit_MD( params.relax_type, params.score_type, params.nperrun, params.cstw ) );

	} else if ( movername.compare("partialabinitio") == 0 ) {
		core::Size res1, res2;
		bool is_terminus;
		get_loop_info( start_struct, res1, res2, is_terminus );
		TR << "assigned loop: " << res1 << " " << res2 << std::endl;

		new_wu = WorkUnit_SamplerOP( new WorkUnit_PartialAbinitio( params.nperrun, true ) );
	}

	// Execute mover
	new_wu->decoys().add( start_struct );
	new_wu->run();

	return new_wu->decoys();
}

void
Serial_Refine::dump_structures( protocols::wum::SilentStructStore const &new_structs,
	bool score_only,
	std::string prefix ) const
{
	core::io::silent::SilentFileData sfd;
	std::string filename = prefix + ".out";

	core::Size istr( 0 );
	for ( SilentStructStore::const_iterator it = new_structs.begin();
			it != new_structs.end(); ++it ) {
		istr++;
		sfd.write_silent_struct( *(*it), filename, score_only );
	}
}

core::pose::Pose
Serial_Refine::get_average_structure( SilentStructStore &decoys,
	utility::vector1< core::Size > const touse,
	std::string const columnname,
	bool const minimize
) const
{
	assert( decoys.size() > 0 );

	// averager: construct via SilentStructStore
	pose::Pose pose;

	// Filter to be used for averaging
	SilentStructStore decoys_touse;
	for ( SilentStructStore::const_iterator it = decoys.begin();
			it != decoys.end(); ++it ) {
		core::Size const i = (*it)->get_energy( columnname );
		if ( touse.contains( i ) ) decoys_touse.add( *it );
	}

	io::silent::SilentStructOP ss = decoys.get_struct( 0 );
	ss->fill_pose( pose ); //default is FA_STANDARD

	// pose being used as reference
	StructAvrgMover averager( pose, decoys_touse, minimize );
	// pose being used as output
	averager.apply( pose );

	TR << "Structure averaging done." << std::endl;

	return pose;
}


} // namespace mpi_refinement
} // namespace protocols


