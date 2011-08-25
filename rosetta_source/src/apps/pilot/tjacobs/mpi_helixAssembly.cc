// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file mpi_helixAssembly.cc
///
/// @brief

/// @author Tim jacobs

// Unit headers
#include <devel/init.hh>
#include <apps/pilot/tjacobs/HelixAssemblyMover.cc>

/// Core headers
#include <core/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>

//Utility headers
#include <utility/mpi_util.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.functions.hh>

// option key includes
#include <basic/options/util.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//#ifdef USEMPI
#include <mpi.h>
//#endif

#include <stdio.h>


struct job{
  std::string pdbToSearch;
  std::string queryPdb;
  int round;
  int frag1_start;
  int frag1_end;
  int frag2_start;
  int frag2_end;
};

// run protocol
int
main( int argc, char * argv [] )
{
  using namespace std;
  using namespace utility;
  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  option.add( HelixAssembly::query_structure_path, "" );
  option.add( HelixAssembly::frag1_start, "" );
  option.add( HelixAssembly::frag1_end, "" );
  option.add( HelixAssembly::frag2_start, "" );
  option.add( HelixAssembly::frag2_end, "" );
  option.add( HelixAssembly::single_helix_rmsd_cutoff, "");
  option.add( HelixAssembly::helix_pair_rmsd_cutoff, "");
  option.add( HelixAssembly::helix_cap_distance_cutoff, "");
  option.add( HelixAssembly::helix_contact_distance_cutoff, "");
  option.add( HelixAssembly::minimum_helix_contacts, "");
  option.add( HelixAssembly::helices_to_add, "");
//  numHelicesToAdd

  option.add( HelixAssembly::foobar, "");

  // initialize core
  devel::init(argc, argv);

//  MPI_Init(&argc,&argv);

  if(mpi_rank()==0){

      int outputFilesCounter = 0;
      //Read in pdb library (should be handled by -l later)
      std::string pdbLibraryFilename( option[ HelixAssembly::foobar ].value() );
      utility::vector1<std::string> pdb_library;

      utility::io::izstream libraryData( pdbLibraryFilename.c_str() );
      if ( !libraryData.good() ) {
          utility_exit_with_message("unable to open input file file: "+pdbLibraryFilename+"\n");
      }
      while ( libraryData.good() ) {
          std::string name;
          libraryData.getline(name);
          if ( libraryData.good() ) pdb_library.push_back( name );
      }

      //Read original query pdb into a string
      std::string originalQueryPdbPath( option[ HelixAssembly::query_structure_path].value());
      utility::io::izstream originalQueryPdbData( originalQueryPdbPath.c_str() );
      std::string originalQueryPdbString;
      utility::slurp(originalQueryPdbData, originalQueryPdbString);

      vector1<job> jobQueue;
      vector1<int> freeProcessors;
      Size helices_to_add = option[ HelixAssembly::helices_to_add ];

      //Add all of the first-round jobs to the job queue
      for(Size i=1; i <= pdb_library.size(); ++i){
          job newJob;
          //newJob.queryPdb = /*read from input here*/;

          //read pdb into string
          utility::io::izstream data( pdb_library[i].c_str() );
          std::string pdbString;
          utility::slurp(data, pdbString);

          newJob.queryPdb = originalQueryPdbString;
          newJob.pdbToSearch = pdbString;
          newJob.round = 1;
          jobQueue.push_back(newJob);
      }

      cout << "Done populating the job queue with " << jobQueue.size() << " jobs." << endl;

      //distribute initial jobs to all processors
      for(Size i=1; i < mpi_nprocs(); ++i){
          if(jobQueue.size() > 0){
              job tempJob = jobQueue[jobQueue.size()];
              jobQueue.pop_back();
              send_string_to_node(i, tempJob.queryPdb);
              send_string_to_node(i, tempJob.pdbToSearch);
              send_integer_to_node(i, tempJob.round);
              cout << "Job send to node: " << i << endl;
          }
          else{
              freeProcessors.push_back(i);
          }
      }

      //wait for a job to finish, queue all jobs resulting from that job finishing, and start a new job.
      //Do this until all jobs are finished.
      while(true){
          //receive the signal from any node, that this specified node is finished
          int completedNode = receive_integer_from_node(MPI_ANY_SOURCE);
          cout << "Job completed on node: " << completedNode << endl;

          //the completed job first tells us how many pdbs were generated
          int numberOfPdbsToReturn = receive_integer_from_node(completedNode);
          cout << "The job created " << numberOfPdbsToReturn << " pdbs. Good job, job!" << endl;

          //for each generated pdb, add to a vector
          vector1<string> returnedPdbs;
          for(Size j=1; j <= numberOfPdbsToReturn; ++j){
              returnedPdbs.push_back(receive_string_from_node(completedNode));
          }

          //What round was this job in? if less than the last round, add jobs accordingly
          int roundCompleted = receive_integer_from_node(completedNode);
          cout << "The job that just finished was in round " << roundCompleted << " of " << helices_to_add << endl;
          if(roundCompleted < helices_to_add){
              for(Size i=1; i <= returnedPdbs.size(); ++i){
                  for(Size j=1; j <= pdb_library.size(); ++j){
                      job tempJob;
                      tempJob.queryPdb=returnedPdbs[i];
                      tempJob.pdbToSearch=pdb_library[j];
                      tempJob.round=++roundCompleted;
                      jobQueue.push_back(tempJob);
                  }
              }
          }
          //If this is the last round, dump the final output
          else{
              //These PDBs have completed their last round, print them to disk
              cout << "dumping " << returnedPdbs.size() << " output files." << endl;
              for(Size i=1; i <= returnedPdbs.size(); ++i){
                  outputFilesCounter++;
                  utility::io::ozstream outputStream;
                  outputStream.open("outputFile_" + utility::to_string(outputFilesCounter) + ".pdb");
                  outputStream << returnedPdbs[i];
                  outputStream.close();
              }
          }

          //Use the newly free processor for the next job
          if(jobQueue.size() > 0){
              //LOOP THROUGH FREE PROCESSORS HERE

              job tempJob = jobQueue[jobQueue.size()];
              jobQueue.pop_back();
              send_string_to_node(completedNode, tempJob.queryPdb);
              send_string_to_node(completedNode, tempJob.pdbToSearch);
              send_integer_to_node(completedNode, tempJob.round);
          }
          else{
              break;
          }
      }
  }
  else{ //ALL OTHER NODES

      //while(receive_bool_from_node(0)){
      while(true){
        //receive the query pdb, pdbToSearch, and round we are in
        //std::string queryPdb = receive_string_from_node(0);
        std::string queryPdb = receive_string_from_node(0);
        std::string pdbToSearch = receive_string_from_node(0);
        int round = receive_integer_from_node(0);

        cout << "Node " << mpi_rank() << " received a job, ready for action!" << endl;

        HelixAssemblyMover helixAssembler;
        helixAssembler.set_query_structure_string(queryPdb);

        core::pose::Pose poseToSearch;
  //      PDBInfoOP pdb_info(new PDBInfo());
  //      poseToSearch.pdb_info(pdb_info);
        core::import_pose::pose_from_pdbstring(poseToSearch, pdbToSearch);

        cout << "Node " << mpi_rank() << " successfully created a pose with " << poseToSearch.total_residue() << " residues" << endl;
        vector1<string> returnedPdbs = helixAssembler.apply(poseToSearch);
        cout << "Finished apply! returning " << returnedPdbs.size() << " pdb(s)." << endl;

        //tell head node that this node is done
        send_integer_to_node(0, mpi_rank());

        //tell head node how many pdbs we are sending
        send_integer_to_node(0, returnedPdbs.size());
        for(Size i=1; i <= returnedPdbs.size(); ++i){
            send_string_to_node(0, returnedPdbs[i]);
        }

        //let the head node what round we just finished
        send_integer_to_node(0, round);
      }
  }
  MPI_Finalize();
}
