// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ken Jung
/// @brief lua test
/// just playing around with lua exported cpp objects and the wrapper

#include <devel/init.hh>

// this includes lua.h and luabind.h
#include <utility/lua/LuaObject.hh>
#include <utility/lua/LuaIterator.hh>
#include <utility/excn/Exceptions.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <basic/Tracer.hh>
#include <iostream>

static THREAD_LOCAL basic::Tracer trmain( "test" );

OPT_1GRP_KEY(File, m, file)

int
main( int argc, char * argv [] )
{
    try {
	using namespace utility::lua;
	NEW_OPT(m::file, "file", "");
	devel::init(argc, argv);

	std::cout << "I'm in cpp now, going to execute some lua now" << std::endl;
	lua_State * lstate_ = luaL_newstate();
	luaL_openlibs(lstate_);
	std::string somelua =
		"a = 5; "
		"b = 5; "
		"print( 'Printing a+b from lua, '..(a + b) );";
	std::cout << " Will be executing the following lua string" << std::endl;
	std::cout << somelua << std::endl;
	int err = luaL_dostring ( lstate_, somelua.c_str() );
	if( err == 1) {
		trmain << "Lua interpreting of string failed. Error is:" << std::endl;
		trmain << lua_tostring(lstate_, -1) << std::endl;
	}

	LuaObject a( luabind::globals(lstate_)["a"]);
	std::cout << "Printing value of lua var 'a' from cpp, it is " << a.to<int>() << std::endl;

	// test lua table iteration
	somelua =
		"c = {1,2,3,4,"
		"{0,-1,-2},5,6,7};";
	std::cout << " Will be executing the following lua string" << std::endl;
	std::cout << somelua << std::endl;

	err = luaL_dostring ( lstate_, somelua.c_str() );
	if( err == 1) {
		trmain << "Lua interpreting of string failed. Error is:" << std::endl;
		trmain << lua_tostring(lstate_, -1) << std::endl;
	}

	LuaObject c( luabind::globals(lstate_)["c"]);
	std::cout << "Iterating through lua var 'c' from cpp" << std::endl;
	for (LuaIterator i=c.begin(), end; i != end; ++i) {
		if( luabind::type( (*i).raw() ) == LUA_TTABLE ) {
			std::cout << "Entering nested table" << std::endl;
			for (LuaIterator j=(*i).begin(), end; j != end; ++j) {
				std::cout << (*j).to<int>() << std::endl;
			}
			std::cout << "Leaving nested table" << std::endl;
		} else {
			std::cout << (*i).to<int>() << std::endl;
		}
	}

	// trying out pose
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;
	luabind::open(lstate_);
	core::pose::lregister_Pose(lstate_);
	core::pose::lregister_util(lstate_);
	core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	PoseOP pose_sp(new Pose());
	core::import_pose::pose_from_file( *pose_sp, *residue_set, option[ m::file]().name() , core::import_pose::PDB_file);
	luabind::globals(lstate_)["p"] = pose_sp;
	somelua =
		"print('Printing pose info from lua');"
		"print(p:total_residue());"
		"core.pose.setExtraScore(p, 'test', 10.0);"
		"print(core.pose.getExtraScore(p, 'test'));";
	err = luaL_dostring ( lstate_, somelua.c_str() );
	if( err == 1) {
		trmain << "Lua interpreting of string failed. Error is:" << std::endl;
		trmain << lua_tostring(lstate_, -1) << std::endl;
	}

    } catch ( utility::excn::EXCN_Base const & e ) {
			std::cerr << "caught exception " << e.msg() << std::endl;
			return -1;
    }
    return 0;
}
