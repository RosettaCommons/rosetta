// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file xml_task_printer.cc
/// @brief takes in a task, outputs it
/// @author Ben Stranges (stranges@unc.edu)

// Unit headers
#include <devel/init.hh>
#include <protocols/moves/Mover.hh>
//project Headers

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility Headers
#include <basic/Tracer.hh>

#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
#include <protocols/toolbox/task_operations/RestrictToInterfaceVectorOperation.hh>

// Option keys
#include <basic/options/keys/optE.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// Job distributor
#include <protocols/jd2/JobDistributor.hh>

static basic::Tracer TR("XMLprinter");

using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using basic::T;
using basic::Error;
using basic::Warning;

// application specific options
namespace XMLprinter {
  basic::options::BooleanOptionKey const make_individual_files( "make_individual_files" );
	basic::options::BooleanOptionKey const print_pymol_selection( "print_pymol_selection" );
}
// mover deffinition
class XMLprinterMover : public protocols::moves::Mover {
public:

  XMLprinterMover();

  virtual void apply( core::pose::Pose& pose );

  virtual protocols::moves::MoverOP clone() const {
    return new XMLprinterMover( *this );
  }

  virtual
  std::string
  get_name() const {
    return "XMLprinterMover";
  }

  virtual protocols::moves::MoverOP fresh_instance() const {
    return new XMLprinterMover;
  }
  core::pack::task::TaskFactoryOP setup_tf( core::pack::task::TaskFactoryOP task_factory );
  std::set< core::Size > fill_designable_set( core::pose::Pose & pose, core::pack::task::TaskFactoryOP & tf );
  std::set< core::Size > fill_packable_set( core::pose::Pose & pose, core::pack::task::TaskFactoryOP & tf );

  void write_output_file(core::pose::Pose& pose, std::string filename, std::set<core::Size> design_set, std::set<core::Size> packable_set);
private:
  //std::string tagfile_name_;
  scoring::ScoreFunctionOP scorefxn_;
  bool output_files_;
	bool pymol_selection_;
};

//constructor
XMLprinterMover::XMLprinterMover(){
  //tagfile_name_ =  option[ optE::parse_tagfile ]() ;
  scorefxn_ = core::scoring::getScoreFunction();
  output_files_ = option[ XMLprinter::make_individual_files ];
	pymol_selection_ = option[ XMLprinter::print_pymol_selection ];
}

///@brief load custom TaskOperations according to an xml-like utility::tag file
core::pack::task::TaskFactoryOP
XMLprinterMover::setup_tf( core::pack::task::TaskFactoryOP task_factory ) {
 using namespace core::pack::task::operation;
 if ( option[ optE::parse_tagfile ].user() ) {
   std::string tagfile_name( option[ optE::parse_tagfile ]() );
   TaskOperationFactory::TaskOperationOPs tops;
   TaskOperationFactory::get_instance()->newTaskOperations( tops, tagfile_name );
   for ( TaskOperationFactory::TaskOperationOPs::iterator it( tops.begin() ), itend( tops.end() ); it != itend; ++it ) {
     task_factory->push_back( *it );
   }
 } else {
   task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );
 }

 return task_factory;
}

///@brief return the set of residues that are designable based given pose
std::set< Size > XMLprinterMover::fill_designable_set( core::pose::Pose & pose, core::pack::task::TaskFactoryOP & tf ) {

  std::set< Size > designable_set;
  core::pack::task::PackerTaskOP design_task( tf->create_task_and_apply_taskoperations( pose ) );

  //#ifndef NDEBUG
  TR<< "Task for " << pose.pdb_info()->name() << " is: \n" << *(design_task)  << std::endl;
  //#endif

    // iterate over all residues
  for ( Size ii = 1; ii<= design_task->total_residue(); ++ii ) {
    if( design_task->being_designed( ii ) )
      designable_set.insert( ii );
  }

  return designable_set;

} // end fill_designable_set

///@brief return the set of residues that are packable based given pose
std::set< Size > XMLprinterMover::fill_packable_set( core::pose::Pose & pose, core::pack::task::TaskFactoryOP & tf ) {

  std::set< Size > packable_set;
  core::pack::task::PackerTaskOP packable_task( tf->create_task_and_apply_taskoperations( pose ) );

  //#ifndef NDEBUG
  TR<< "Task for " << pose.pdb_info()->name() << " is: \n" << *(packable_task)  << std::endl;
  //#endif

    // iterate over all residues
  for ( Size ii = 1; ii<= packable_task->total_residue(); ++ii ) {
    if( packable_task->being_packed( ii ) )
      packable_set.insert( ii );
  }

  return packable_set;

} // end fill_packable_set

///@brief, write individual files for each input structure
void XMLprinterMover::write_output_file(pose::Pose& pose, std::string filename, std::set<core::Size> design_set, std::set<core::Size> packable_set ){
  utility::io::ozstream output_stream;
  utility::io::ozstream packable_output_stream;
  std::string const output_name = filename + "_designdef.txt";
	std::string const packable_output_name = filename + "_packabledef.txt";
  output_stream.open(output_name);
  //itterate and write
/// measure seq recov
 TR << "Number of residues being designed: " << design_set.size() << std::endl;
  for ( std::set< core::Size >::const_iterator it = design_set.begin(), end = design_set.end(); it != end; ++it ) {
    //first print the chain, then the residue number from the pdb
    output_stream << pose.pdb_info()->chain(*it) << " " << pose.pdb_info()->number(*it) <<std::endl;
  }
 output_stream.close();
 packable_output_stream.open(packable_output_name);
 TR << "Number of residues being packed: " << packable_set.size() << std::endl;
  for ( std::set< core::Size >::const_iterator it = packable_set.begin(), end = packable_set.end(); it != end; ++it ) {
    //first print the chain, then the residue number from the pdb
    packable_output_stream << pose.pdb_info()->chain(*it) << " " << pose.pdb_info()->number(*it) <<std::endl;
  }
	packable_output_stream.close();

}//end write_output_file

//begin mover apply
void XMLprinterMover::apply (pose::Pose& pose ) {
  utility::file::FileName filename = pose.pdb_info()->name();
  std::string pdbname( filename.base() );


  // score the input pose for kicks and giggles
  (*scorefxn_)( pose );

  // figure out the task & neighbor info
  core::pack::task::TaskFactoryOP task_factory( new core::pack::task::TaskFactory );
  // setup what residues we are going to look at...
  setup_tf( task_factory );
  //debugging
	TR<< pose.pdb_info()->name() << " foldtree: "<< pose.fold_tree() << std::endl;
  //task_factory->push_back( new protocols::toolbox::task_operations::RestrictToInterfaceVectorOperation);
  std::set< Size > design_set;
  design_set = fill_designable_set( pose, task_factory );
  std::set< Size > packable_set;
	packable_set = fill_packable_set( pose, task_factory );

  //now make output
  if(output_files_){
    TR<< "Writing output for:" << pdbname << std::endl;
    write_output_file(pose, pdbname, design_set, packable_set);
  }
	if(pymol_selection_){
		TR << "pymol selection not implemented." << std::endl;
}

}//end apply

// run protocol
int
main( int argc, char * argv [] )
{

	try {

  //using namespace protocols;
  using namespace protocols::jd2;
  using namespace protocols::moves;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  option.add( XMLprinter::make_individual_files, "Make indivdiual files for output").def(false);
	option.add( XMLprinter::print_pymol_selection, "Tracer output of pymol selection").def(true);
  // initialize core
  devel::init(argc, argv);
  protocols::jd2::JobDistributor::get_instance()->go( new XMLprinterMover );
  std::cout << "Done! -------------------------------\n";

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}//end main
