// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Ken Jung
/// @brief testing elscripts

#include <devel/init.hh>

#include <protocols/elscripts/SingleNode.hh>
#ifdef USEBOOSTMPI
#include <protocols/elscripts/MPI_Master.hh>
#include <protocols/elscripts/MPI_Slave.hh>
#endif

#include <utility/lua/LuaObject.hh>
#include <utility/lua/LuaIterator.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/els.OptionKeys.gen.hh>

#include <basic/Tracer.hh>
#include <boost/algorithm/string.hpp>
#include <alloca.h>

static basic::Tracer TR( "elscripts" );

#ifdef USEBOOSTMPI
namespace mpi = boost::mpi;
#endif

int
main( int argc, char * argv [] )
{
	try {

		// command line flag preprocessing
		// also known as c string hell
		std::string elscriptsscript;
		for ( int i = 1; i < argc; i++ ) {
			if ( strcmp(argv[i], "-els:script") == 0 ) {
				if ( i + 1 < argc ) {
					elscriptsscript = std::string(argv[i+1]);
					break;
				}
			}
		}

		int result_argc;
		char ** result_argv;

		if ( ! elscriptsscript.empty() ) {
			lua_State * tmpstate = luaL_newstate();
			// really shouldn't need lua libraries for loading flags, but some people might want print or something
			luaL_openlibs(tmpstate);

			// load lua script
			int err = luaL_dofile(tmpstate, elscriptsscript.c_str());
			if ( err == 1 ) {
				std::cout << "Loading lua script '" << elscriptsscript << "' failed. Error is:" << std::endl;
				std::cout << lua_tostring(tmpstate, -1) << std::endl;
			}
			std::vector< char * > mod_argv;
			int mod_argc = 0;
			utility::lua::LuaObject flags(luabind::globals(tmpstate)["flags"]);
			if ( flags ) {
				for ( utility::lua::LuaIterator i=flags.begin(), end; i != end; ++i ) {
					std::string flag = (*i).to<std::string>();
					std::vector<std::string> strs;
					boost::split(strs, flag, boost::is_any_of(" "));
					for ( std::vector<std::string>::iterator itr = strs.begin(); itr != strs.end(); itr++ ) {
						// this is where one goes to read about the dangers of alloca
						char * cflag = (char * ) alloca( itr->size() + 1 );
						std::copy(itr->begin(), itr->end(), cflag);
						cflag[itr->size()] ='\0';
						mod_argv.push_back( cflag );
						mod_argc++;
					}
				}

				// put elscripts parsed flags BEFORE command line flags
				// that way, command line can be used to override elscripts parsed flags
				// I can't think of a case where you would want it the other way around

				result_argc = argc + mod_argc;
				result_argv = ( char ** ) alloca( (argc + mod_argc) * sizeof( char * ) );

				result_argv[0] = argv[0];

				int idx = 0;
				for ( ; idx < mod_argv.size(); idx++ ) {
					result_argv[idx+1] = mod_argv[idx];
				}

				for ( int j = 1; j < argc; j++ ) {
					result_argv[idx+j] = argv[j];
				}
			}
		} else {
			std::cout << "Need to specify -els:script ! " << std::endl;
			std::exit(9);
		}


#ifdef USEBOOSTMPI
	mpi::environment env(result_argc, result_argv);
	mpi::communicator world;
#endif

		devel::init(result_argc, result_argv);

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		bool pool = option[OptionKeys::els::pool]();
		int num_masters =
			(option[OptionKeys::els::num_traj]() / option[OptionKeys::els::traj_per_master]() ) +
			(!( option[OptionKeys::els::num_traj]() % option[OptionKeys::els::traj_per_master]() == 0 )) ;

		using namespace protocols::elscripts;
#ifdef USEBOOSTMPI

	// rank 0 is the pool if enabled, otherwise masters start at 0
	// rank 0 to num_masters-1 are all masters
	// slaves are everything after that, % num_masters to give assignments

	if( world.size() == 1 || option[OptionKeys::els::singlenode]() ) {
		// single node version
		SingleNode role;
		role.go();
	} else {

		if( pool && world.rank() == 0 ) {
			// i am the pool!
			// but i dont have pool, so placeholder slave for now
			// Pool role( world );
			MPI_Slave role( world, 0 );
			role.go();
		} else if ( world.rank() <= num_masters - !pool ) {
			// i am a master
			std::vector< int > slaves;
			for( int i = num_masters + pool; i < world.size(); i++ ) {
				if( (i - pool ) % num_masters == world.rank() + pool )
					slaves.push_back(i);
			}
			MPI_Master role( world, slaves, option[OptionKeys::els::num_traj]() );
			role.go();
		} else {
			// i am a slave
			int master = ( (world.rank() - pool) % num_masters ) + pool;
			MPI_Slave role( world, master );
			role.go();
		}

	}
#else
		// single node version
		SingleNode role;
		role.go();
#endif

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

