// -*- mode:c++;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Javier Castellanos javier@cyrusbio.com
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/relax/FastRelax.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

#include "boost/date_time/posix_time/posix_time.hpp"
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <iostream>
#include <sstream>
#include <map>

#include <apps/pilot/cyrus/micro/micro_utils.hh>
#include <apps/pilot/cyrus/micro/micro_callbacks.hh>

#include <SimpleAmqpClient/SimpleAmqpClient.h>
#include <SimpleAmqpClient/AmqpLibraryException.h>
#include <jsoncpp/json/json.h>

const std::string UUID = boost::lexical_cast<std::string>(boost::uuids::random_generator()());

static THREAD_LOCAL basic::Tracer TR( "micro:" + UUID);

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace AmqpClient;


// OPTION REGISTRATION
OPT_KEY( String, broker )
OPT_KEY( Integer, port )
OPT_KEY( String, exchange )
OPT_KEY( String, work_manager )
OPT_KEY( String, reply_to )
OPT_KEY( String, username )
OPT_KEY( String, password )
OPT_KEY( String, vhost )
OPT_KEY( Integer, n_messages )
OPT_KEY( Integer, connection_attempts )
OPT_KEY( Integer, trial_timeout )
OPT_KEY( Integer, lifespan )

OPT_KEY( Boolean, exchange_passive )
OPT_KEY( Boolean, exchange_durable )
OPT_KEY( Boolean, exchange_auto_delete )

OPT_KEY( Boolean, queue_passive )
OPT_KEY( Boolean, queue_durable )
OPT_KEY( Boolean, queue_exclusive )
OPT_KEY( Boolean, queue_auto_delete )

OPT_KEY( Boolean, reply_queue_passive )
OPT_KEY( Boolean, reply_queue_durable )
OPT_KEY( Boolean, reply_queue_exclusive )
OPT_KEY( Boolean, reply_queue_auto_delete )

OPT_KEY( Boolean, consumer_no_local )
OPT_KEY( Boolean, consumer_no_ack )
OPT_KEY( Boolean, consumer_exclusive )

typedef boost::function<Json::Value (const Json::Value&)>  TaskCallback;
typedef std::map< std::string, TaskCallback> Callbacks;

std::string connect_queue( Channel::ptr_t channel, std::string const & queue, std::string const & consumer_tag = "")
{
	if(consumer_tag != "") {
		TR << "Cancelling consumer " << consumer_tag << std::endl;
		channel->BasicCancel(consumer_tag);
	}
	TR << "Connecting worker to queue " << queue << " consumer_tag = " << consumer_tag<< std::endl;
	channel->DeclareExchange(option[exchange], Channel::EXCHANGE_TYPE_DIRECT, option[exchange_passive], option[exchange_durable], option[exchange_auto_delete]);
	channel->DeclareQueue(queue,option[queue_passive], option[queue_durable], option[queue_exclusive], option[queue_auto_delete]);
	return channel->BasicConsume(queue, "", option[consumer_no_local], option[consumer_no_ack], option[consumer_exclusive], 1);
}


class Application
{
	public:
		void register_options()
		{
			using namespace basic::options;
			using namespace basic::options::OptionKeys;
			NEW_OPT(broker, "address of the amqp broker","127.0.0.1");
			NEW_OPT(port, "amqp broker port", 5672);
			NEW_OPT(exchange, "exchange from were to consume messages", "cad");
			NEW_OPT(work_manager, "Address of the work manager","");
			NEW_OPT(reply_to, "reply to queue", "");
			NEW_OPT(username, "broker username", "guest");
			NEW_OPT(password, "broker password", "guest");
			NEW_OPT(vhost, "vhost", "/");
			NEW_OPT(n_messages, "How many messages to consume, by default infinite(-1)", -1);
			NEW_OPT(connection_attempts, "How many time to try to establish a connection to the broker", 10);
			NEW_OPT(trial_timeout, "For many milliseconds to sleep before connections attempts.", 200);
			NEW_OPT(lifespan, "After a job is done exit if the executable has been running for longer that lifespan(seconds)", 600);

			NEW_OPT(exchange_passive, "", false);
			NEW_OPT(exchange_durable, "", true);
			NEW_OPT(exchange_auto_delete, "", true);

			NEW_OPT(queue_passive, "", false);
			NEW_OPT(queue_durable, "", true);
			NEW_OPT(queue_exclusive, "", false);
			NEW_OPT(queue_auto_delete, "", true);

			NEW_OPT(reply_queue_passive, "", false);
			NEW_OPT(reply_queue_durable, "", true);
			NEW_OPT(reply_queue_exclusive, "", false);
			NEW_OPT(reply_queue_auto_delete, "", true);

			NEW_OPT(consumer_no_local, "", false);
			NEW_OPT(consumer_no_ack, "", false);
			NEW_OPT(consumer_exclusive, "", false);
		}

		void register_callbacks()
		{
			using namespace cyrus::micro::callbacks;
			callback_map.insert(std::pair<std::string, TaskCallback>("pdb2pose",  pdb_to_pose));
			callback_map.insert(std::pair<std::string, TaskCallback>("pose2pdb",  pose_to_pdb));
			callback_map.insert(std::pair<std::string, TaskCallback>("fast_relax",  fast_relax));
			callback_map.insert(std::pair<std::string, TaskCallback>("repack",  repack));
			callback_map.insert(std::pair<std::string, TaskCallback>("design",  design));
			callback_map.insert(std::pair<std::string, TaskCallback>("minimize",  minimize));
		}

		Json::Value execute_task( Json::Value & doc)
		{
			std::string task = doc["job_info"].get("task", "empty").asString();
			TR << "task: " << task << std::endl;

			Json::Value result;
			if(callback_map.count(task) == 0)
			{
				TR << "Can't find callback for " << task << std::endl;
				doc["job_info"]["status"] = "TASK_NOT_FOUND";
			}
			else
			{
				try
				{
					result =  callback_map[task](doc);
					doc["job_info"]["status"] = "EXECUTED";
				}
				catch (utility::excn::EXCN_Base const & e)
				{
					doc["job_info"]["status"] = "FAILED";
				}
			}
			return result;
		}

	private:
		Callbacks callback_map;
};

int main( int argc, char** argv )
{

	Application App;
	App.register_options();
	App.register_callbacks();
	Json::FastWriter fastWriter;

	boost::posix_time::ptime init_time = boost::posix_time::second_clock::universal_time();

	try
	{

		devel::init(argc, argv);

		Channel::ptr_t channel = NULL;
		int connection_attempt = 1;
		while(channel == NULL)
		{
			try
			{
				channel = Channel::Create(option[broker], option[port], option[username], option[password], option[vhost]);
			}
			catch ( AmqpClient::AmqpLibraryException const & e )
			{
				if( connection_attempt == option[connection_attempts])
				{
					TR.Error << "Unable to start connection with broker" << std::endl;
					exit(-1);
				}
				else
				{
					TR.Warning << "Unable to establish connection, attempt " << connection_attempt << " of " << option[connection_attempts] << std::endl;
					usleep(option[trial_timeout]*1000);
					connection_attempt +=1;
				}
			}
		}

		Json::Value work = cyrus::micro::utils::get_work_queue(option[work_manager]);
		std::string consumer_tag = connect_queue(channel, work["queue"].asString());
		TR << "broker assign consumer tag " << consumer_tag << std::endl;

		std::string reply_to_queue(option[reply_to]);

		if( reply_to_queue != "" )
		{
			channel->DeclareQueue(reply_to_queue, option[reply_queue_passive], option[reply_queue_durable], option[reply_queue_exclusive], option[reply_queue_auto_delete]);
			channel->BindQueue(reply_to_queue, "amq.direct", "");
		}

		TR << "Worker " << UUID << " is ready to consume messages." <<  std::endl;
		int n_message = 0;
		int no_queue_trial = 0;
		while(n_message != option[n_messages])
		{
			boost::posix_time::ptime now = boost::posix_time::second_clock::universal_time();
			if((now - init_time).total_seconds() > option[lifespan])
			{
				TR.Info << "lifespan exceded! turning off worker " << UUID <<std::endl;
				exit(0);
			}

			if(option[n_messages] > 0 && n_message == option[n_messages]) {
				TR.Info << "Max number of messages consumed! turning off worker " << UUID <<std::endl;
				exit(0);
			}
				
			Envelope::ptr_t input_envelope;
			bool new_msg = channel->BasicConsumeMessage(consumer_tag, input_envelope, 100);
			if( !new_msg)
			{
				Json::Value new_work = cyrus::micro::utils::get_work_queue(option[work_manager]);
				if(new_work["queue"].asString() != work["queue"].asString()) {
					TR << "Changing queue to " << new_work["queue"].asString() << std::endl;	
					work = new_work;
					if(work["queue"].asString() == "none") {
						no_queue_trial++;
						boost::posix_time::ptime now = boost::posix_time::second_clock::universal_time();
						if(no_queue_trial == 12 && (now - init_time).total_seconds() > 10.0) {
							TR << "Unable to find an available queue" << std::endl;
							break;
						} else {
							TR << "No available queue, sleeping for 10 seconds." << std::endl;
							sleep(5);
							continue;
						}
					} else {
						no_queue_trial = 0;
						try {
							connect_queue(channel, work["queue"].asString(), consumer_tag);
						} catch(AmqpClient::ConsumerTagNotFoundException const & e) {
							sleep(5);
							continue;
						}
					}
				}
			}
			else
			{
				TR << "Processing message number "<< n_message  << std::endl;
				BasicMessage::ptr_t msg_in = input_envelope->Message();

				Json::Value task;
				Json::Reader reader;
				bool parsingSuccessful = reader.parse( msg_in->Body(), task, false );
				TR << "Message received, job = " << fastWriter.write(task["job_info"]) << std::endl;

				if(!parsingSuccessful)
				{
					TR.Error << "Failed to parse message" << std::flush << std::endl;
					continue;
				}

				TR << "Executing task...";
				Json::Value task_result;
				try
				{
					task_result = App.execute_task(task);
					TR << " Done, acknowledging input for job = " << fastWriter.write(task["job_info"]) << std::endl;
					channel->BasicAck(input_envelope);
				}
				catch ( utility::excn::EXCN_Base const & e )
				{
					TR.Error << "job failed! requeing message for job = " << fastWriter.write(task["job_info"]) << std::endl;
					// reject with requeue
					channel->BasicReject(input_envelope, true);
				}

				// Add job info
				task_result["job_info"] = task["job_info"];

				Json::FastWriter writer;
				std::string json_output = writer.write(task_result);

				TR.Debug << "payload: "<< json_output << std::endl;

				BasicMessage::ptr_t result_msg = BasicMessage::Create();
				result_msg->CorrelationId(msg_in->CorrelationId());
				result_msg->Body(json_output);
				if( reply_to_queue == "" )
				{
					TR.Debug << "sending message to queue "<< msg_in->ReplyTo() << std::endl;
					channel->BasicPublish("", msg_in->ReplyTo(), result_msg);
				}
				else
				{
					channel->BasicPublish("", reply_to_queue , result_msg);

				}

				n_message +=1;
			}
		}

	}
	catch ( utility::excn::EXCN_Base const & e )
	{
		std::cerr << "Caught exception " << e.msg() << std::flush << std::endl;
		return -1;
	}
	return 0;
}
