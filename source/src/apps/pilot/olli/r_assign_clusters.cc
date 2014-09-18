// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange

#include <protocols/moves/Mover.hh>
#include <devel/init.hh>

#include <protocols/jd2/NoOutputJobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/SilentFileJobInputter.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
//#include <core/io/silent/ProteinSilentStruct.hh>

#include <protocols/moves/mc_convergence_checks/Pool_ConvergenceCheck.hh>

// Utility headers
#include <basic/options/option_macros.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// option key includes
#include <ObjexxFCL/format.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/exit.hh>

//Auto Headers
#include <utility/excn/EXCN_Base.hh>
#include <numeric/random/random.hh>

static thread_local basic::Tracer tr( "main" );

using namespace core;
using namespace protocols;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::moves::mc_convergence_checks;

OPT_2GRP_KEY( File, in, file, cluster )

void register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  OPT( in::file::s );
  OPT( in::file::silent );
  OPT( out::file::silent );
  NEW_OPT( in::file::cluster, "read cluster information from this silent-file", "cluster.out" );
}

// Forward
class AssignClusterToolMover : public moves::Mover {
public:
  AssignClusterToolMover() :
    pool_( new Pool_RMSD( option[ in::file::cluster ]() ) ) {
  }
  virtual void apply( core::pose::Pose& );
  virtual void apply( core::io::silent::SilentStruct const& );

	///@brief returns true if this mover does not use the pose that is passed into the apply function
	/// this allows the jobdistributor to safe time, because creating a pose can be bery slow
	virtual bool accesses_jobinputter_in_apply() const {
		return true;
	}

	virtual moves::MoverOP fresh_instance() const {
		return new AssignClusterToolMover( *this );
	}
	std::string get_name() const { return "AssignClusterToolMover"; }

private:
  Pool_RMSD_OP pool_;
};

void AssignClusterToolMover::apply( core::pose::Pose &pose ) {
  protocols::jd2::JobDistributor* jd = protocols::jd2::JobDistributor::get_instance();
  protocols::jd2::SilentFileJobInputterOP sfd_inputter(
     dynamic_cast< protocols::jd2::SilentFileJobInputter* >( jd->job_inputter().get() )
  );

  using namespace core::io::silent;
  if ( sfd_inputter ) {
    SilentStruct const& current_struct( sfd_inputter->struct_from_job( jd->current_job() ) );
    apply( current_struct );
  } else { // okay no idea how to get the silent-struct directly use conventional path-way
    jd->job_inputter()->pose_from_job( pose, jd->current_job() );
    SilentStructOP current_struct = SilentStructFactory::get_instance()->get_silent_struct_out(); //new ProteinSilentStruct;
    //unfortunately in general it might be difficult to find out what is appropriate
    // should one use get_silent_struct_out ?
    current_struct->fill_struct( pose );
    apply( *current_struct );
  }
}

void AssignClusterToolMover::apply( core::io::silent::SilentStruct const& pss ) {
  std::string decoy;
  core::Real value;
	core::Size const id( pool_->evaluate( pss, decoy, value ) );

  protocols::jd2::JobDistributor* jd = protocols::jd2::JobDistributor::get_instance();
  jd->current_job()->add_string_string_pair( "cluster_tag", decoy );
  jd->current_job()->add_string_real_pair( "cluster_rmsd", value );
	jd->current_job()->add_string_string_pair( "cluster_id", ObjexxFCL::lead_zero_string_of( id, 8 ) );
}


void run() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  if ( !basic::options::option[ basic::options::OptionKeys::in::path::database ].user() ) {
    basic::options::option[ basic::options::OptionKeys::in::path::database ].def( "/work/olange/minirosetta_database");
  }
  moves::MoverOP mover = new AssignClusterToolMover;
  jd2::SilentFileJobOutputterOP jd_out = new jd2::SilentFileJobOutputter;
  jd_out->set_write_no_structures( true ); //only score file
  protocols::jd2::JobDistributor::get_instance()->go( mover, jd_out );

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// =============================== MAIN ============================================================
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
  register_options();
  devel::init( argc, argv );

  try{
    run();
  } catch ( utility::excn::EXCN_Base& excn ) {
    excn.show( std::cerr );
  }
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
  return 0;
}


