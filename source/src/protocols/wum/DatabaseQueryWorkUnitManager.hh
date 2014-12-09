// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DatabaseQueryWorkUnitManager.hh
///
/// @brief An MPI work unit manager that runs a given database query and creates a DatabaseEntryWorkUnit for each of the result rows. This class is templated such
///        that the specific type of DatabaseEntryWorkUnit created is the template class.

/// @author Tim Jacobs

#ifndef INCLUDED_protocols_wum_DatabaseQueryWorkUnitManager_hh
#define INCLUDED_protocols_wum_DatabaseQueryWorkUnitManager_hh

//Unit
#include <protocols/wum/DatabaseEntryWorkUnit.hh>

#include <protocols/wum/MPI_WorkUnitManager.hh>
#include <protocols/wum/WorkUnitBase.hh>

//Core
#include <core/types.hh>

//Utility and basic
#include <basic/database/sql_utils.hh>

#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>

//Basic
#include <basic/Tracer.hh>

//C++
#include <string>
#include <map>

static thread_local basic::Tracer TR( "protocols.wum.DatabaseQueryWorkUnitManager" );

namespace protocols{
namespace wum{

template <class T>
class DatabaseQueryWorkUnitManager : public protocols::wum::MPI_WorkUnitManager {
public:

    DatabaseQueryWorkUnitManager(core::Size master_rank, utility::sql_database::sessionOP db_session, std::string query_string, std::string wu_type);

    //void set_defaults();

    virtual ~DatabaseQueryWorkUnitManager(){}

	virtual void go();

protected: // overloaded functions

    virtual void init(){/*no resume functionality yet*/}

	virtual void process_inbound_wus();

	virtual void process_outbound_wus();

protected: // added functions

    void create_work_units_from_query(std::string wu_type);

protected: // Accesors

	core::Size master_rank(){ return master_rank_; }

	core::Size my_emperor(){ return my_emperor_; }

private:

	// static settings
	const core::Size my_emperor_;
	const core::Size master_rank_;

    // database specific
	utility::sql_database::sessionOP db_session_;
    std::string query_string_;
};

template <class T>
DatabaseQueryWorkUnitManager<T>::DatabaseQueryWorkUnitManager(core::Size master_rank, utility::sql_database::sessionOP db_session, std::string query_string, std::string wu_type):
MPI_WorkUnitManager( 'M' ),
my_emperor_(0),
master_rank_( master_rank ),
db_session_( db_session ),
query_string_(query_string)
{
    create_work_units_from_query(wu_type);
}

template <class T>
void
DatabaseQueryWorkUnitManager<T>::create_work_units_from_query(std::string wu_type){
    using cppdb::statement;
    using cppdb::result;
    using namespace std;

    //ensure that template class is derived from database entry work unit
    //    BOOST_CONCEPT_ASSERT((DatabaseEntryWorkUnit<T>));

	TR << "Adding DB WU for each row in query result" << std::endl;

    statement query_statement(basic::database::safely_prepare_statement(query_string_,db_session_));

    result res(basic::database::safely_read_from_database(query_statement));

    core::Size wu_counter = 0;
    while(res.next()){
        ++wu_counter;

        map<string,string> row_map;

        for(int i=0; i<res.cols(); ++i){
            string key(res.name(i));

            string value("");
            res.fetch(i,value);

            row_map[key]=value;
        }

        //create the speicific
        //        T new_wu = new T(row_map);
        //        new_wu->set_wu_type(wu_type);

        DatabaseEntryWorkUnitOP new_wu = new T(row_map);
        new_wu->set_wu_type(wu_type);

        outbound().add( new_wu );
    }
    TR << "Added " << wu_counter << " database WUs to queue" << endl;
}

template <class T>
void
DatabaseQueryWorkUnitManager<T>::go(){
	// initialize master (this is a virtual function call and this function is overloaded by the children of this class)
	TR << "Init Master: " << protocols::wum::mpi_rank() << std::endl;
	init();

	TR << "Master Node: Waiting for job requests..." << std::endl;
	while(true){

        //check the outbound queue and send MPI messages accordingly - If nothing outbound then prepare to recieve MPI messages.
        //The tag of the MPI message dictates the next action (ie. sending a work unit/recieving a work unit)
		TR << "Master: processing msgs.." << std::endl;
		process_incoming_msgs();

        //Collect all the incoming work units, cast them to the proper type and execute their jobs
		TR << "Master: process incoming" << std::endl;
		process_inbound_wus();

		TR << "Master: process outbound" << std::endl;
		process_outbound_wus();

		// ok, we've done all our work, now wait until we hear from our slaves
		process_incoming_msgs( true );

		print_stats_auto();
	}
}

/// @brief Process the inbound WUs (these will be the same DB work units send during the constructor - but now they have the data we want).
template <class T>
void
DatabaseQueryWorkUnitManager<T>::process_inbound_wus(){

    if( inbound().size() > 0 ){
		TR << "Processing inbound WUs on master.." << std::endl;
	}

    //Create a database transaction for all database queries created by the slaves
    db_session_->begin();
	while( inbound().size() > 0 )
	{
		protocols::wum::WorkUnitBaseOP  next_wu =  inbound().pop_next();
		runtime_assert( next_wu );

		// skip returning waiting WUs
		if ( next_wu->get_wu_type() == "waitwu" ) continue;

		// Upcast to a DatabaseEntryWorkUnit WU
		DatabaseEntryWorkUnit* db_wu = dynamic_cast<  DatabaseEntryWorkUnit * > ( (protocols::wum::WorkUnitBase*) (&(*next_wu)) );

		// If upcast was unsuccessful - warn and ignore.
		if ( db_wu == NULL ){
			TR << "WARNING: Master recieved a non-db WU" << std::endl;
			next_wu->print( TR );
			continue;
		}

        cppdb::statement query_statement(basic::database::safely_prepare_statement(db_wu->result_query_string(),db_session_));
        basic::database::safely_write_to_database(query_statement);

	}
    //Commit the transaction
    db_session_->commit();
}

/// @brief create the database WUs and add them to the outgoing queue
template <class T>
void
DatabaseQueryWorkUnitManager<T>::process_outbound_wus(){
    //All outbound wus are taken care of in the constructor
}

} //namespace wum
} //namespace protocols

#endif
