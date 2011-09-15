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
#include <devel/helixAssembly/HelixAssemblyMover.hh>
#include <devel/helixAssembly/HelixAssemblyJob.hh>

// Devel headers
#include <devel/init.hh>

//Utility headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.functions.hh>

// option key includes
#include <basic/options/util.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/HelixAssembly.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Boost MPI includes
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>

// C++ headers
#include <stdio.h>

BOOST_IS_MPI_DATATYPE(HelixAssemblyJob)
static basic::Tracer TR("mpi_helixAssembly");

// run protocol
int
main( int argc, char * argv [] )
{
  using namespace std;
  using namespace utility;
  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  // initialize devel
  devel::init(argc, argv);

  boost::mpi::communicator world;

  //Should I check to make sure MPI is initialized?
  if(world.rank()==0){

      int job_id_counter = 1;
      int num_completed_jobs = 0;
      int currently_running_jobs = 0;

      vector1<file::FileName> input_lists( option[ in::file::l ]() );
      utility::vector1<file::FileName> pdb_library;

      for(Size p = 1; p <= input_lists.size(); p++){
        utility::io::izstream current_input_list( input_lists[p] );
        if ( !current_input_list.good() ) {
            utility_exit_with_message("unable to open input file file: "+input_lists[p].name()+"\n");
        }
        while ( current_input_list.good() ) {
            std::string name;
            current_input_list.getline(name);
            if ( current_input_list.good() ) pdb_library.push_back( file::FileName(name) );
        }
      }

      //Read original query pdb into a string
      std::string originalQueryPdbPath( option[ HelixAssembly::query_structure_path].value());
      utility::io::izstream originalQueryPdbData( originalQueryPdbPath.c_str() );
      std::string originalQueryPdbString;
      utility::slurp(originalQueryPdbData, originalQueryPdbString);

      vector1<HelixAssemblyJob> job_queue;
      vector1<int> freeProcessors;
      Size helices_to_add = option[ HelixAssembly::helices_to_add ];

      TR << "library size: " << pdb_library.size() << endl;

      //Add all of the first-round jobs to the job queue
      for(Size i=1; i <= pdb_library.size(); ++i){
          //create a HelixAssemblyJob
          HelixAssemblyJob temp_job;
          temp_job.set_job_name(pdb_library[i].base());
          temp_job.set_query_structure(originalQueryPdbString);
          //THIS is your problem.
          temp_job.set_search_index(i);
          temp_job.set_frag1_start(option[ HelixAssembly::frag1_start]);
          temp_job.set_frag1_end(option[ HelixAssembly::frag1_end]);
          temp_job.set_frag2_start(option[ HelixAssembly::frag2_start]);
          temp_job.set_frag2_end(option[ HelixAssembly::frag2_end]);
          temp_job.set_round(1);
          job_queue.push_back(temp_job);
      }
      TR << "Done populating the job queue with " << job_queue.size() << " jobs." << endl;

      //distribute initial jobs to all processors
      for(Size i=1; i < world.size(); ++i){
          if(job_queue.size() > 0){
              HelixAssemblyJob temp_job = job_queue[job_queue.size()];
              job_queue.pop_back();

              //lazily load the search string
              utility::io::izstream data( pdb_library[temp_job.get_search_index()].name().c_str() );
              std::string pdbString;
              utility::slurp(data, pdbString);
              temp_job.set_search_structure(pdbString);

              //send the job to the node
              world.send(i, 0, true);
              world.send(i, 0, boost::mpi::skeleton(temp_job));
              world.send(i, 0, boost::mpi::get_content(temp_job));
              ++currently_running_jobs;
              TR << "Job sent to node: " << i << endl;
          }
          //if we have more processors than we have jobs, keep track of which processors are free (for use by jobs created later)
          else{
              freeProcessors.push_back(i);
          }
      }

      //wait for a job to finish, queue all jobs resulting from that job finishing, and start a new job.
      //Do this until all jobs are finished.
      while(true){
          int completed_node;
          world.recv(boost::mpi::any_source,0,completed_node);

          TR << "Node " << completed_node << " finished a job" << endl;
          ++num_completed_jobs;
          --currently_running_jobs;

          TR << "Number of completed jobs: " << num_completed_jobs << endl;
          cout << "Number of completed jobs: " << num_completed_jobs << endl;
          cout << "Number of currently running jobs: " << num_completed_jobs << endl;

          int num_new_jobs;
          world.recv(completed_node,0,num_new_jobs);
          TR << "Node " << completed_node << " wants to add " << num_new_jobs << " new job(s)" << endl;
          vector1<HelixAssemblyJob> new_jobs;
          //receive all new jobs from completed node
          for(Size i=1; i<=num_new_jobs; ++i){
              HelixAssemblyJob temp_job;
              world.recv(completed_node,0,boost::mpi::skeleton(temp_job));
              world.recv(completed_node,0,boost::mpi::get_content(temp_job));
              new_jobs.push_back(temp_job);
          }

          TR << "Added all the new jobs" << endl;

          for(Size i=1; i <= new_jobs.size(); ++i){
              if(new_jobs[i].get_round() <= option[HelixAssembly::helices_to_add]){
                  TR << "job " << new_jobs[i].get_job_name() << " is on round " << new_jobs[i].get_round() <<
                      " of " << option[HelixAssembly::helices_to_add] << endl;
                  for(Size j=1; j <= pdb_library.size(); ++j){
                      //copy incomplete job and add a pose to search, then add to the queue
                      HelixAssemblyJob temp_job(new_jobs[i]);

                      //only thing we should need to change about the job is the search structure, which only the head node knows
                      //we use an index here so that we don't need to hang onto the string forever (just-in-time loading).
                      temp_job.set_search_index(j);
                      temp_job.set_job_name(temp_job.get_job_name() + "_" + pdb_library[j].base());
                      job_queue.push_back(temp_job);
                  }
                  //Incomplete job output for debugging purposes
                  TR << "Outputting incomplete job " << new_jobs[i].get_job_name() << endl;
                  utility::io::ozstream outputStream;
                  outputStream.open(new_jobs[i].get_job_name() + "_" + to_string(job_id_counter) + "_bundle.pdb");
                  ++job_id_counter;
                  outputStream << new_jobs[i].get_query_structure();
                  outputStream.close();
              }
              else{
                  TR << "Outputting job " << new_jobs[i].get_job_name() << endl;
                  //done. output PDBs
                  utility::io::ozstream outputStream;
                  outputStream.open(new_jobs[i].get_job_name() + "_" + to_string(job_id_counter) + "_bundle.pdb");
                  ++job_id_counter;
                  outputStream << new_jobs[i].get_query_structure();
                  outputStream.close();
              }
          }

          TR << "Jobs left: " << job_queue.size() << endl;
          cout << "Jobs left: " << job_queue.size() << endl;

          //use newly free processor
          if(job_queue.size() > 0){
              //LOOP THROUGH FREE PROCESSORS HERE

              HelixAssemblyJob temp_job = job_queue[job_queue.size()];
              job_queue.pop_back();
              TR << "Sending job " << temp_job.get_job_name() << " (" << temp_job.get_round() <<
                                    ":" << option[HelixAssembly::helices_to_add] << ")" << endl;

              //lazily load the search string
              utility::io::izstream data( pdb_library[temp_job.get_search_index()].name().c_str() );
              std::string pdbString;
              utility::slurp(data, pdbString);
              temp_job.set_search_structure(pdbString);

              world.send(completed_node, 0, true);
              world.send(completed_node, 0, boost::mpi::skeleton(temp_job));
              world.send(completed_node, 0, boost::mpi::get_content(temp_job));
              ++currently_running_jobs;
          }
          //we're all done, wait for all the currently running jobs to finish then tell all the other processors we're done.
          else if(currently_running_jobs == 0){
              //tell all processors that we don't have any more jobs
              for(Size i=1; i < world.size(); ++i){
                  world.send(i, 0, false);
              }
              break;
          }
      }
  }
  else{ //ALL OTHER NODES

      //while(receive_bool_from_node(0)){
      bool moreJobs;
      while(true){
          world.recv(0, 0, moreJobs);
          if(moreJobs){
              HelixAssemblyJob received_job;
              world.recv(0, 0, boost::mpi::skeleton(received_job));
              world.recv(0, 0, boost::mpi::get_content(received_job));
              TR << "Node " << world.rank() << " received a job" << endl;

              HelixAssemblyMover helixAssembler;
              vector1<HelixAssemblyJob> returnedJobs = helixAssembler.apply(received_job);
              TR << "Finished apply! returning " << returnedJobs.size() << " job(s)." << endl;

              //tell head node which job finished
              world.send(0,0,world.rank());

              //send how many jobs to expect, then send each job
              int num_jobs_to_send = returnedJobs.size();
              world.send(0,0,num_jobs_to_send);
              for(Size i=1; i<=returnedJobs.size(); ++i){
                  world.send(0,0,boost::mpi::skeleton(returnedJobs[i]));
                  world.send(0,0,boost::mpi::get_content(returnedJobs[i]));
              }
          }
          else{
              break;
          }
      }
  }
  MPI_Finalize();
  TR << "------------DONE!------------" << endl;
}
