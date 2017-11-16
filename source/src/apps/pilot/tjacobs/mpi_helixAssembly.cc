// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file mpi_sewing.cc
///
/// @brief

/// @author Tim jacobs

// Unit headers
#include <protocols/sewing/HelixAssemblyMover.hh>
#include <protocols/sewing/HelixAssemblyJob.hh>
#include <protocols/features/sewing/HelicalFragment.hh>
#include <protocols/sewing/NativeAtom.hh>
#include <protocols/sewing/NativeResidue.hh>

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
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#ifdef USEBOOSTMPI

//Boost MPI includes
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
BOOST_IS_MPI_DATATYPE(HelixAssemblyJob)
#endif

// C++ headers
#include <stdio.h>

static basic::Tracer TR( "mpi_sewing" );

// run protocol
int
main( int argc, char * argv [] )
{

	try {


#ifdef USEBOOSTMPI

  using namespace std;
  using namespace utility;
  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  // initialize devel
  devel::init(argc, argv);

  boost::mpi::communicator world;

  ///*****DEBUG*****/////

//  if(world.rank()==0){
//      std::map<int, std::vector<std::string> > test;
//      test[1].push_back("FOO");
//      test[50].push_back("BAR");
//      world.send(1, 1, boost::mpi::skeleton(test));
//      world.send(1, 1, boost::mpi::get_content(test));
//      cout << "sent" << endl;
//  }
//  else{
//      std::map<int, std::vector<std::string> > test;
//      world.recv(0,1,boost::mpi::skeleton(test));
//      world.recv(0,1,boost::mpi::get_content(test));
//      cout << test[1][0] << endl;
//      cout << test[50][0] << endl;
//      exit(1);
//  }

  ////***END DEBUG***////

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
      std::string originalQueryPdbPath( option[ sewing::query_structure_path].value());
      utility::io::izstream originalQueryPdbData( originalQueryPdbPath.c_str() );
      std::string originalQueryPdbString;
      utility::slurp(originalQueryPdbData, originalQueryPdbString);

      vector1<HelixAssemblyJob> job_queue;
      vector1<int> freeProcessors;
      Size helices_to_add = option[ sewing::helices_to_add ];

      TR << "library size: " << pdb_library.size() << endl;

      //Add all of the first-round jobs to the job queue
      for(Size i=1; i <= pdb_library.size(); ++i){
          //create a HelixAssemblyJob
          HelixAssemblyJob temp_job;
          temp_job.set_name(pdb_library[i].base());
          temp_job.set_query_structure(originalQueryPdbString);
          temp_job.set_search_index(i);
          temp_job.set_remaining_rounds(option[sewing::helices_to_add]);
          temp_job.set_direction_needed(true); //doesn't matter for first round
          temp_job.set_first_round(true);

          core::Size frag1_start = option[ sewing::frag1_start];
          core::Size frag1_end = option[ sewing::frag1_end];
          core::Size frag2_start = option[ sewing::frag2_start];
          core::Size frag2_end = option[ sewing::frag2_end];

          HelicalFragment fragment1(frag1_start, frag1_end, true);
          fragment1.set_pdb_source(pdb_library[i].name());

          HelicalFragment fragment2(frag2_start, frag2_end, false);
          fragment2.set_pdb_source(pdb_library[i].name());

          std::vector<HelicalFragment> fragments;
          fragments.push_back(fragment1);
          temp_job.set_query_frag_1_index(0);
          fragments.push_back(fragment2);
          temp_job.set_query_frag_2_index(1);
          temp_job.set_fragments(fragments);

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

          TR << "Number of currently running jobs: " << currently_running_jobs << endl;
          cout << "Number of currently running jobs: " << currently_running_jobs << endl;

          //receive each new job individually from the completed node
          std::vector<HelixAssemblyJob> new_jobs;

          core::Size num_new_jobs;
          world.recv(completed_node,0,num_new_jobs);

          cout << "Number of jobs to receive " << num_new_jobs << endl;

          for(Size i=0; i<num_new_jobs; i++){
              HelixAssemblyJob new_job;
              world.recv(completed_node,0,boost::mpi::skeleton(new_job));
              world.recv(completed_node,0,boost::mpi::get_content(new_job));

              new_jobs.push_back(new_job);
          }

          //for some reason this causes a segfault
//          world.recv(completed_node,0,boost::mpi::skeleton(new_jobs));
//          world.recv(completed_node,0,boost::mpi::get_content(new_jobs));

          TR << "Node " << completed_node << " added " << new_jobs.size() << " new job(s)" << endl;

          for(Size i=0; i < new_jobs.size(); ++i){
              if(new_jobs[i].get_remaining_rounds() > 0){
                  TR << "job " << new_jobs[i].get_name() << " has " << new_jobs[i].get_remaining_rounds() <<
                      " remaining" << endl;
                  for(Size j=1; j <= pdb_library.size(); ++j){
                      //copy incomplete job and add a pose to search, then add to the queue
                      HelixAssemblyJob temp_job(new_jobs[i]);

                      //only thing we should need to change about the job is the search structure, which only the head node knows
                      //we use an index here so that we don't need to hang onto the string forever (just-in-time loading).
                      temp_job.set_search_index(j);
                      temp_job.set_name(temp_job.get_name() + "_" + pdb_library[j].base());

                      //Add a partially incomplete Fragment Residue info to the end of the list of fragment residue infos. This
                      //will be finished in the HelixAssemblyMover(which is when we know the resnums of the new fragmnet)
//                      HelicalFragment newFragmentInfo;
//                      newFragmentInfo.set_pdb_source_file(pdb_library[i].name());
//                      temp_job.get_fragment_info().push_back(newFragmentInfo);

                      job_queue.push_back(temp_job);
                  }
                  //Incomplete job output for debugging purposes
//                  TR << "Outputting incomplete job " << new_jobs[i].get_name() << endl;
//                  utility::io::ozstream outputStream;
//                  outputStream.open(new_jobs[i].get_name() + "_" + to_string(job_id_counter) + "_bundle.pdb");
//                  ++job_id_counter;
//                  outputStream << new_jobs[i].get_query_structure();
//                  outputStream.close();
              }
              else{
                  TR << "Outputting job " << new_jobs[i].get_name() << endl;
                  //done. output PDBs
                  utility::io::ozstream pdb_output_stream;
                  pdb_output_stream.open(new_jobs[i].get_name()+ "_bundle.pdb");
                  ++job_id_counter;
                  pdb_output_stream << new_jobs[i].get_query_structure();
                  pdb_output_stream.close();

                  utility::io::ozstream residue_output_stream;
                  residue_output_stream.open(new_jobs[i].get_name()+ "_residues.txt");
                  residue_output_stream << new_jobs[i].printBundleResidues();
                  residue_output_stream.close();
              }
          }

          TR << "Jobs left in queue: " << job_queue.size() << endl;
          cout << "Jobs left in queue: " << job_queue.size() << endl;

          size_t job_queue_mem_size(0);
          for(core::Size test=1; test<=job_queue.size(); ++test){
        	  job_queue_mem_size += sizeof(string)+ job_queue[test].get_query_structure().capacity() * sizeof(char);
        	  job_queue_mem_size += sizeof(string)+ job_queue[test].get_search_structure().capacity() * sizeof(char);
          }

          cout << "Job queue memory size: " << job_queue_mem_size << endl;

          //use newly free processor
          if(job_queue.size() > 0){
              //LOOP THROUGH FREE PROCESSORS HERE

              HelixAssemblyJob temp_job = job_queue[job_queue.size()];
              job_queue.pop_back();
              TR << "Sending job " << temp_job.get_name() << " to node " << completed_node << endl;

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
              int pid = getpid();
              cout << "Node " << world.rank() << " with PID " << pid << " received job: " << received_job.get_name() << endl;

              HelixAssemblyMover helixAssembler(received_job.get_query_frag_1(), received_job.get_query_frag_2());
              std::vector<HelixAssemblyJob> returned_jobs = helixAssembler.apply(received_job);
//              std::vector<HelixAssemblyJob> returned_jobs;
//              while(true){
//            	  int counter = 1;
//            	  std::vector<HelixAssemblyJob> temp = helixAssembler.apply(received_job);
//            	  returned_jobs.insert(returned_jobs.end(), temp.begin(), temp.end());
//            	  cout << "Round " << counter << "finished" << endl;
//            	  counter++;
//              }
//              helixAssembler.apply(received_job);
//              cout << "Node " << world.rank() << " finished " << received_job.get_name() << ", returning " <<
//                  returned_jobs.size() << " new job(s)." << endl;

              ///MEMORY TESTS
//              std::vector<HelixAssemblyJob> returned_jobs;
//              HelixAssemblyJob return_job1(received_job);
//              HelixAssemblyJob return_job2(received_job);
//              returned_jobs.push_back(return_job1);
//              returned_jobs.push_back(return_job2);

              //tell head node which job finished
              world.send(0,0,world.rank());

              //For whatever reason this causes a segfault
//              world.send(0,0,boost::mpi::skeleton(returned_jobs));
//              world.send(0,0,boost::mpi::get_content(returned_jobs));

              world.send(0,0,returned_jobs.size());

              for(Size i=0; i < returned_jobs.size(); ++i){
                  world.send(0,0,boost::mpi::skeleton(returned_jobs[i]));
                  world.send(0,0,boost::mpi::get_content(returned_jobs[i]));
              }
          }
          else{
              break;
          }
      }
  }
  MPI_Finalize();
  TR << "------------DONE!------------" << endl;

#endif

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
