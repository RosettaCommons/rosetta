// C++ headers
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>

#include <boost/algorithm/string.hpp>

// Rosetta headers
#include <utility/pointer/owning_ptr.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>

#include <devel/init.hh>
#include <protocols/relax/FastRelax.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>

#include <protocols/evaluation/EvaluatorFactory.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/batch_relax.OptionKeys.gen.hh>

#ifdef USEMPI
#include <mpi.h>
#endif

// defines
#define TAG_BATCH_ASSIGN  0
#define TAG_WRITE_REQUEST  1
#define TAG_WRITE_APPROVE  2
#define TAG_WRITE_SUCCESS  3
#define TAG_EXIT_ACK 4
#define TAG_SLAVE_ERROR 5

#define BATCH_STRING_SIZE 2000

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

// new options
OPT_1GRP_KEY(String, batchrelax_mpi, jobfile)

using namespace std;

void wait_for_gdb()
{
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == i)
        sleep(5);
}

string split_string( string str, int index ) {
	vector<string> split_string;
	boost::split( split_string, str, boost::is_any_of("\t ") );
	return split_string[index];
}

bool isUnique( utility::vector1<string> & str )
{
	//brute force
	for ( utility::vector1<string>::iterator it = str.begin(); it != str.end(); ++it ) {
		for ( utility::vector1<string>::iterator it2 = it + 1; it2 != str.end(); ++it2 ) {
			if ( *it == *it2 )
				return false;
		}
	}
	return true;
}

static basic::Tracer TR("main");

int main(int argc, char *argv[]) {

	try {

	using namespace core;
	using namespace protocols;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::io::silent;
	using io::silent::SilentStructFactory;
	using io::silent::SilentStructOP;

	// new options
	NEW_OPT(batchrelax_mpi::jobfile, "jobfile - format is <silentfile native.pdb> on each line", "");

	// initialize core
	devel::init(argc, argv);

	// initialize evaluators
	evaluation::PoseEvaluatorsOP evaluators_( new protocols::evaluation::PoseEvaluators() );
	evaluation::EvaluatorFactory::get_instance()->add_all_evaluators(*evaluators_);

	// handle options
	core::Size nstruct = option[ OptionKeys::out::nstruct ];
	core::Size batch_size = option[ OptionKeys::batch_relax::batch_size ];

	// init residue set
	core::chemical::ResidueTypeSetCAP rsd_set;
	if ( option[ in::file::fullatom ]() )
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	else
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );

	// init sfxn
	core::scoring::ScoreFunctionOP scorefxn;
	scorefxn = core::scoring::getScoreFunction();

	//mpi init
	int numprocs = 1, rank = 1;
	int EXIT_BATCH_ID = -1;
	int BLANK;

	#ifdef USEMPI
	MPI_Status status;
	//MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

// finish init, start batchrelaxing

	//handle master here
	if (rank == 0 && numprocs > 1) {
		// create job list
		if ( ! option[ OptionKeys::batchrelax_mpi::jobfile ].user() )
			utility_exit_with_message( "supply a jobfile with flag -batchrelax_mpi::jobfile in format <silentfile native.pdb>" );
		string jobs_file_str = option[ OptionKeys::batchrelax_mpi::jobfile ];
		ifstream jobs_file( jobs_file_str.c_str() );
		vector<string> jobs;
		string line;
		map<string, bool> filelock;			//only used on master node

		if (jobs_file.is_open()) {
			while ( jobs_file.good() ) {
				getline (jobs_file,line);
				if (line.size() > 0) {
					jobs.push_back( line );
					//cout << line << endl;
				}
			}
			jobs_file.close();
		}

		// parse jobs list and split into batches
		vector<string> batches;				// yes, im mixing vector0 and vector1 in this app - sue me
		SilentFileData sfd;
		string silent_file_name;
		int infile_size;				// number of incoming structures
		int num_batches_in_file;
		utility::vector1<string> tags;
		for ( vector<string>::iterator it = jobs.begin(); it != jobs.end(); ++it ) {	// for each job in list, get structures and divide into batches
			silent_file_name = split_string( *it, 0 );				// sf path is first word on job line (native path is second word)
			tags = sfd.read_tags_fast(silent_file_name);
			infile_size = tags.size();
			if (!isUnique( tags ))
				utility_exit_with_message( silent_file_name + " does not have unique tags. this breaks the batching algorithm because SilentFileData is too stupid to handle tag collision and doesn't support indexing by number" );
			num_batches_in_file = infile_size / batch_size + 1;
			int start, end = 1;
			bool breakLoop = false;
			// create batches by computing tag indexes and pushing start/end of range with silent file name
			for ( start = 1; ; start += batch_size ) {
				end = start + batch_size - 1;
				if ( end >= infile_size ) {
					end = infile_size;							// don't run over end of file, recompute index of last structure
					breakLoop = true;							// we're done with this file, signal exit from loop
				}
				for ( int nstruct_i = 0; nstruct_i < nstruct; ++nstruct_i )			// duplicate batch for each nstruct
					batches.push_back( *it + " " + SSTR(start) + " " + SSTR(end) );
				if (breakLoop)
					break;
			}
		}

		// if batch size is less than 1/4 of batch_size, move it to the end - "stable sort" preserves input order
		// do this later - it's not THAT important, but a good idea nonetheless

		//for ( vector<string>::iterator it = batches.begin(); it != batches.end(); ++it )
			//cout << *it << endl;

		#ifdef USEMPI
		// broadcast batch list - send size, then each batch line by line
		int batch_size = batches.size();
		char* batch_string = new char[BATCH_STRING_SIZE];

		MPI_Bcast( &batch_size, 1, MPI_INT, 0, MPI_COMM_WORLD );
		for ( vector<string>::iterator it = batches.begin(); it != batches.end(); ++it ) {
			if ( it->size() > BATCH_STRING_SIZE )
				utility_exit_with_message("batch string size exceeds BATCH_STRING_SIZE (" + SSTR(BATCH_STRING_SIZE) + ")");
			std::strcpy( batch_string, it->c_str() );
			MPI_Bcast( batch_string, BATCH_STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD );
		}

		// initial batch assignments - iterate thru nodes and send
		int next_batch;
		for (next_batch = 0; next_batch < batches.size(); ++next_batch) {
			if (next_batch > numprocs - 2)	// - 2 to count head node and 0-based indexing
				break;			// all slaves are full, move on
			MPI_Send( &next_batch, 1, MPI_INT, next_batch+1, TAG_BATCH_ASSIGN, MPI_COMM_WORLD );
		}

		// if we have more procs than batches, kill off the remaining processors
		if ( next_batch <= numprocs - 2 )
		{
			for ( int ii = next_batch; ii <= numprocs - 2; ++ii ) {
                                MPI_Send( &EXIT_BATCH_ID, 1, MPI_INT, ii+1, TAG_BATCH_ASSIGN, MPI_COMM_WORLD );
			}
		}

		// go thru rest of batches - block until a node frees up and then assign that node a new batch
		int batch_id;
		vector< pair<int,int> > write_approve_queue;	// pair( batch_id, node_id )
		bool all_batches_sent = false;
		int num_exit_acks = 1;				//include self by starting with 1
		while ( num_exit_acks < numprocs ) {

			//try to clear the write approve queue
			if ( write_approve_queue.size() ) {
				TR << "write approve queue size: " << write_approve_queue.size() << endl;
				for ( vector< pair<int,int> >::iterator it = write_approve_queue.begin(); it != write_approve_queue.end(); ) {
					if ( ! filelock[ split_string(batches[it->first],0) ] ) {	// if unlocked, send write approval and erase from queue
						MPI_Send( &batch_id, 1, MPI_INT, it->second, TAG_WRITE_APPROVE, MPI_COMM_WORLD );
						filelock[ split_string(batches[batch_id],0) ] = true;
						it = write_approve_queue.erase( it );
					}
					else
						++it;
				}
			}

			//get next message
			// TAG_WRITE_REQUEST - when slave is done, it will ask master to check for filelock on output
			// 			if locked, push it on the write_approve_queue with batch_id/node_id
			// 			else, master will send write approval to this slave only
			// 			and master will lock this file until TAG_WRITE_SUCCESS
			// TAG_WRITE_SUCCESS - when slave completes a write, it will ask master to release file lock
			// 			master will release lock
			// 			if there are batches left, assign them - else send EXIT_BATCH_ID to indicate exit
			// TAG_EXIT_ACK      - when slave receives batch_id = EXIT_BATCH_ID, it will acknowledge and terminate
			// 			master will increment exit_ack count
			MPI_Recv( &batch_id, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

			if (status.MPI_TAG == TAG_WRITE_REQUEST) {
				TR << "master Recv loop: processing TAG_WRITE_REQUEST from " << status.MPI_SOURCE << endl;
				if ( filelock[ split_string(batches[batch_id],0) ] )
					write_approve_queue.push_back( make_pair( batch_id, status.MPI_SOURCE ) );

				else {
					MPI_Send( &batch_id, 1, MPI_INT, status.MPI_SOURCE, TAG_WRITE_APPROVE, MPI_COMM_WORLD );
					filelock[ split_string(batches[batch_id],0) ] = true;
				}
			}
			else if (status.MPI_TAG == TAG_WRITE_SUCCESS) {
				TR << "master Recv loop: processing TAG_WRITE_SUCCESS from " << status.MPI_SOURCE << endl;
				filelock[ split_string(batches[batch_id],0) ] = false;
				if ( next_batch < batches.size() ) {
					MPI_Send( &next_batch, 1, MPI_INT, status.MPI_SOURCE, TAG_BATCH_ASSIGN, MPI_COMM_WORLD );
					++next_batch;
				}
				else {
					MPI_Send( &EXIT_BATCH_ID, 1, MPI_INT, status.MPI_SOURCE, TAG_BATCH_ASSIGN, MPI_COMM_WORLD );
				}
			}
			else if (status.MPI_TAG == TAG_EXIT_ACK) {
				TR << "master Recv loop: processing TAG_EXIT_ACK from " << status.MPI_SOURCE << endl;
				++num_exit_acks;
			}
			else if (status.MPI_TAG == TAG_SLAVE_ERROR) {
				TR << "master Recv loop: processing TAG_SLAVE_ERROR from " << status.MPI_SOURCE << endl;
				++num_exit_acks;
			}
			//else

			// insert timeout code here

		} //end recv loop
		#endif
	} //end master


	//handle slaves here
	else if (numprocs > 1) {
		vector<string> batches;
		int batch_id;
		SilentFileDataOP sfd_in;
		core::pose::PoseOP native_pose;
		string silent_filename = "";
		string silent_out_filename = "";

		#ifdef USEMPI

		// receive batch list here
		char batch_string[BATCH_STRING_SIZE];
		int expected_batch_size;

		// get batch size and then get all batch strings
		MPI_Bcast( &expected_batch_size, 1, MPI_INT, 0, MPI_COMM_WORLD );
		for( int i = 0; i < expected_batch_size; ++i )
		{
			MPI_Bcast( &batch_string, BATCH_STRING_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD );
			batches.push_back( batch_string );
		}
		TR << "expected batch_size: " << expected_batch_size << " actual batch_size: " << batches.size() << endl;

		// block for messages until receive batch_id: EXIT_BATCH_ID, then quit loop
		while ( true )
		{
			MPI_Recv( &batch_id, 1, MPI_INT, 0, TAG_BATCH_ASSIGN, MPI_COMM_WORLD, &status );
			if ( batch_id == EXIT_BATCH_ID ) {
				MPI_Send( &BLANK, 1, MPI_INT, 0, TAG_EXIT_ACK, MPI_COMM_WORLD );
				break;
			}
			TR << "slave Recv loop: received batch_id " << batch_id << " corresponding to " << batches[batch_id] << endl;

			// run batch
		//	cout << "running " << batches[batch_id] << endl;
		//	wait_for_gdb();

			SilentFileDataOP sfd_out = new SilentFileData();

			if ( silent_filename != split_string(batches[batch_id],0) ) {
				silent_filename = split_string(batches[batch_id],0);
				sfd_in = new SilentFileData();
				sfd_in->read_file( silent_filename );

				string native_pdb_file = split_string(batches[batch_id],1);
				native_pose = core::import_pose::pose_from_pdb( native_pdb_file );
			}
			int start = atof(split_string(batches[batch_id],2).c_str());
			int end = atof(split_string(batches[batch_id],3).c_str());

			// grab the input structures that we want to work on
			utility::vector1<string> tags = sfd_in->read_tags_fast(silent_filename);
			std::vector<SilentStructOP> input_structs;
			for ( int tag_index = start; tag_index <= end; ++tag_index) {
				//SilentStructOP new_struct = sfd_in.get_structure( tags[tag_index] ).clone();
				input_structs.push_back( sfd_in->get_structure( tags[tag_index] ).clone() );
			}

			nstruct = 1; 	// use batch list to distribute nstruct rather than forcing a sequential loop
					// get rid of deep copy when it becomes obvious i'll never use a sequential loop

			for ( core::Size j=0;j<nstruct;j++ ) {
				// make a deep copy of the input_structs list so we can reuse input_structs
				std::vector<SilentStructOP> relax_structs;
				for( std::vector < SilentStructOP >::const_iterator it = input_structs.begin();
					it != input_structs.end();
					++ it )
				{
					SilentStructOP new_struct;
					new_struct = (*it)->clone();
					relax_structs.push_back( new_struct );
				}

				protocols::relax::FastRelax relax( scorefxn, option[ OptionKeys::relax::sequence_file ]() );
				relax.set_native_pose( native_pose );

				TR << "BATCHSIZE: " <<  relax_structs.size() << endl;
				long starttime = time(NULL);
				relax.batch_apply( relax_structs );
				long endtime = time(NULL);
				TR << "TIME: " << endtime - starttime << " seconds" << endl;
				TR << "RESULTING DECOYS SIZE: " << relax_structs.size() << endl;

				// Now save the resulting decoys
				for( std::vector < SilentStructOP >::const_iterator it = relax_structs.begin();
					it != relax_structs.end();
					++ it )
				{
					/*if( evaluators_->size()){
						core::pose::Pose cpose;
						input.fill_pose( cpose, *rsd_set );
						evaluators_->apply( cpose, "tag" , *(*it) );
					}*/
					sfd_out->add_structure( *(*it) );
				}
			} //end nstruct loop

			// request write approval
			MPI_Send( &batch_id, 1, MPI_INT, 0, TAG_WRITE_REQUEST, MPI_COMM_WORLD );
			// block for write approval
			MPI_Recv( &batch_id, 1, MPI_INT, 0, TAG_WRITE_APPROVE, MPI_COMM_WORLD, &status );
			// write
			// todo: use JobOutputter::affixed_number_name convention for output names
			silent_out_filename = option[ OptionKeys::out::prefix ]();
			silent_out_filename += utility::file_basename( silent_filename ) + ".batchrelax";
			silent_out_filename += option[ OptionKeys::out::suffix ]();
			TR << "Saving to " << silent_out_filename << endl;
			sfd_out->write_all( silent_out_filename );
			//cout << "writing to " << split_string( batches[batch_id], 0 ) << endl;
			MPI_Send( &batch_id, 1, MPI_INT, 0, TAG_WRITE_SUCCESS, MPI_COMM_WORLD );
		}
		#endif
	}
	//if no slaves, run batches on master
	else {
		TR << "only 1 core? use src/apps/pilot/mtyka/batchrelax.cc - dont want to refactor for this use case" << endl;
	}

	#ifdef USEMPI
	MPI_Finalize();
	#endif

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << endl;
		return -1;
	}

	return 0;
}
