// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/elscripts/BaseRole.cc
/// @brief  BaseRole, which handles a lot of common functions between all the elscripts roles
/// 				including a lot of lua stuff
/// @author Ken Jung


#include <fstream>
#include <streambuf>

#ifdef USELUA
#include <protocols/elscripts/BaseRole.hh>
#include <protocols/elscripts/modules/MonteCarlo.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/els.OptionKeys.gen.hh>

#include <utility/lua/LuaIterator.hh>
#include "luabind/class_info.hpp"


#include <utility/Factory.hh>
#include <protocols/inputter/InputterStream.hh>
#include <protocols/inputter/Inputter.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/outputter/Outputter.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/filters/FilterFactory.hh>

#include <core/io/serialization/Pipe.fwd.hh>
#include <core/pose/util.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>

// starting includes for initialize_scorefxns()
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.fwd.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
//ending

// stupid fucking calculatorfactory
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/InterfaceSasaDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <core/pose/metrics/simple_calculators/InterfaceDeltaEnergeticsCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
// end stupid fucking calculatorfactory

#include <basic/Tracer.hh>

// elscripts master
namespace protocols {
namespace elscripts {

static thread_local basic::Tracer TR( "protocols.elscripts.BaseRole" );

void lregister_BaseRole( lua_State * lstate ) {
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("elscripts")
		[
			luabind::class_<BaseRole>("BaseRole")
					.def("reparse_def", &BaseRole::reparse_def)
		]
	];
}
void BaseRole::lua_init(){
  using namespace basic::options;
    //load lua
    lstate_ = luaL_newstate();
    luabind::open(lstate_);
    luabind::bind_class_info(lstate_);
    luaL_openlibs(lstate_);

    // start everything into a elscripts namespace
    std::string action = "DELIM(\n"
      "els = {}\n"
      "setmetatable(els, {__index = _G })\n"
      "do\n"
      "  local _ENV = els\n"
      "  dtasks = {}\n"
      "  dscorefxns = {}\n"
      "  dmovers = {}\n"
      "  dfilters = {}\n"
      "  dinputters = {}\n"
      "  dinputterstream = {}\n"
      "  doutputters = {}\n"
      "  dworkunits = {}\n"
      "end\n"
      ")DELIM";
    int err = luaL_dostring ( lstate_, action.c_str() );
    if( err == 1) {
      TR << "Creating els namespace failed. Error is:" << std::endl;
      TR << lua_tostring(lstate_, -1) << std::endl;
      std::exit(9);
    }

		// register some classes with lua
		core::io::serialization::lregister_PipeMap(lstate_);
		core::io::serialization::lregister_Pipe(lstate_);
		core::pose::lregister_Pose(lstate_);
		core::pose::lregister_util(lstate_);
		protocols::moves::lregister_Mover(lstate_);
		protocols::filters::lregister_Filter(lstate_);
		core::scoring::lregister_ScoreFunction(lstate_);
		protocols::moves::lregister_SerializableState(lstate_);
		protocols::outputter::lregister_Outputter(lstate_);
		protocols::inputter::lregister_Inputter(lstate_);
		core::pack::task::operation::lregister_TaskOperation(lstate_);

    // load elscripts vars
    if( option[OptionKeys::els::vars].user() ) {
      err = luaL_dostring ( lstate_, option[OptionKeys::els::vars]().c_str() );
      if( err == 1) {
        TR << "Lua interpreting of elscripts::vars failed. Error is:" << std::endl;
        TR << lua_tostring(lstate_, -1) << std::endl;
      }
    }

		// load lua script
    if( ! option[OptionKeys::els::script].user() ){
      TR << "Must specify a elscript to run" << std::endl;
      std::exit(9);
    } else {
      std::string scriptname = option[OptionKeys::els::script]().name();
      // err = luaL_doString(lstate_, scriptname.c_str());
      // We need to load the script inside a closure, to share ENV among all the user-defined fxns
      // I'm sure there's a beautiful way to do this in pure lua with upvalues, but I'm ???
      // Instead, we'll do it the ugly way and prepend lines to the actual script
      // Because I can't figure out how to combine closures, I'll predefine loop_every()
      // With this setup, I also can't, after the fact, fix missing run() or proceed().  Oh well
      std::ifstream script( scriptname );
      std::string script_contents((std::istreambuf_iterator<char>(script)),
                      std::istreambuf_iterator<char>());
      std::string closure_script = "DELIM(\n"
        "do\n"
        "  local _ENV = _ENV\n"
        "  function els_setenv(t) _ENV = t end\n"
        "  function loop_every()\n"
        "    master:make_wu_until_limit(\"default\", 1)\n"
        "  end\n"
        "  )DELIM"
        // I guess manually load in modules here
        + modules::MonteCarlo
        + script_contents
        + "(\n"
        "end\n"
        ")";
      err = luaL_dostring ( lstate_, closure_script.c_str() );
      if( err == 1) {
        TR << "Loading lua script '" << scriptname << "' failed. Error is:" << std::endl;
        TR << lua_tostring(lstate_, -1) << std::endl;
        std::exit(9);
      }

      // run definitions() in lua script
      err = luaL_dostring ( lstate_, "definitions()" );
      if( err == 1) {
        TR << "Calling lua function definitions() failed. Error is:" << std::endl;
        TR << lua_tostring(lstate_, -1) << std::endl;
        std::exit(9);
      }

    }

}
void BaseRole::instantiate_workunits() {
  // run workunits() in lua script
  // this "instantiates" the wus on the lua sid
  // must be called after everything else has already been instantiated
  // so that nothing is pointing at nil
  int err = luaL_dostring ( lstate_, "define_workunits()" );
  if( err == 1) {
    TR << "Calling lua function define_workunits() failed. Error is:" << std::endl;
    TR << lua_tostring(lstate_, -1) << std::endl;
    std::exit(9);
  }
}

void BaseRole::instantiate_scorefxns() {
	using namespace core::scoring;
	using namespace core::scoring::symmetry;
	using namespace utility::lua;

	TR << "----------Instantiating Score Functions----------" << std::endl;
	// i might move this somewhere else, it is quite involved
	// ported straight from mini/src/protocols/jd2/parser/ScoreFunctionLoader.cc

  luabind::globals(lstate_)["els"]["scorefxns"] = luabind::newtable( lstate_ );
  scorefxns_.raw( luabind::globals(lstate_)["els"]["scorefxns"] );

	LuaObject dscorefxns(luabind::globals(lstate_)["els"]["dscorefxns"]);
	for (LuaIterator i=dscorefxns.begin(), end; i != end; ++i) {
		TR << "Instantiating score function named " << i.skey() << std::endl;
		// scorefxn
		ScoreFunctionOP scorefxn;
		if( (*i)["weights"] && (*i)["patch"]  ) {
			scorefxn = ScoreFunctionFactory::create_score_function(
					(*i)["weights"].to<std::string>(),
					(*i)["patch"].to<std::string>()
			);
		} else if ( (*i)["weights"] ) {
			scorefxn = ScoreFunctionFactory::create_score_function( (*i)["weights"].to<std::string>() );
		} else {
			scorefxn = new ScoreFunction;
			scorefxn->reset();
		}

		// reweight
		if( (*i)["reweight"]) {
			for (LuaIterator j=(*i)["reweight"].begin(), end; j != end; ++j) {
				ScoreType const type = score_type_from_name( j.skey() );
				scorefxn->set_weight( type, (*j).to<core::Real>() );
			}
		}

		// energy method options
		if( (*i)["emoptions"] ) {

			core::scoring::methods::EnergyMethodOptions emoptions( scorefxn->energy_method_options() );
			core::scoring::hbonds::HBondOptionsOP hboptions( emoptions.hbond_options() );

			for (LuaIterator j=(*i)["emoptions"].begin(), end; j != end; ++j) {
				if( j.skey() == "softrep_etable" && (*j).to<bool>() ) {
						emoptions.etable_type( core::scoring::FA_STANDARD_SOFT );
				} else if( j.skey() == "exclude_protein_protein_fa_elec" ) {
					emoptions.exclude_protein_protein_fa_elec( (*j).to<bool>() );
				} else if( j.skey() == "exclude_DNA_DNA" ) {
					emoptions.exclude_DNA_DNA( (*j).to<bool>() );
				} else if( j.skey() == "use_hb_env_dep_DNA" ) {
					hboptions->use_hb_env_dep_DNA( (*j).to<bool>() );
				} else if( j.skey() == "use_hb_env_dep" ) {
					hboptions->use_hb_env_dep( (*j).to<bool>() );
				} else if( j.skey() == "smooth_hb_env_dep" ) {
					hboptions->smooth_hb_env_dep( (*j).to<bool>() );
				} else if( j.skey() == "decompose_bb_hb_into_pair_energies" ) {
					hboptions->decompose_bb_hb_into_pair_energies( (*j).to<bool>() );
				}

				scorefxn->set_energy_method_options( emoptions );
			}
		}

		// hotspot hashing
		// why is this more special than a normal reweight?
		if( (*i)["hs_hash"] ) {
			scorefxn->set_weight( backbone_stub_constraint, (*i)["hs_hash"].to<core::Real>() );
		}

		//symmetry
		bool const scorefxn_symm = (*i)["symmetric"] ? (*i)["symmetric"].to<bool>() : 0 ;
		if (scorefxn_symm) {
			scorefxn = ScoreFunctionOP( new SymmetricScoreFunction( scorefxn ) );
		}

    // owning ptrs suck
    ScoreFunctionSP tmpsp( scorefxn.get() );
    scorefxn.relinquish_ownership();
    luabind::globals(lstate_)["els"]["scorefxns"][ i.skey() ] = tmpsp;
	}
	// add default score12 if score12 isn't yet taken
  if( luabind::type( scorefxns_["score12"].raw() ) == 0 ) {
    ScoreFunctionOP scorefxn = get_score_function();
    ScoreFunctionSP tmpsp( scorefxn.get() );
    scorefxn.relinquish_ownership();
    luabind::globals(lstate_)["els"]["scorefxns"]["score12"] = tmpsp;
  }
	TR << "----------Finished Instantiating Score Functions----------" << std::endl;
}

void BaseRole::instantiate_tasks() {
	using namespace core::pack::task::operation;
	using namespace utility::lua;

	TR << "----------Instantiating TaskOperations----------" << std::endl;

  luabind::globals(lstate_)["els"]["tasks"] = luabind::newtable( lstate_ );
  tasks_.raw( luabind::globals(lstate_)["els"]["tasks"] );

	LuaObject dtasks( luabind::globals(lstate_)["els"]["dtasks"]);
	//Factory<MMover> * mmoverfactory = Factory<MMover>::get_instance();
	for (LuaIterator i=dtasks.begin(), end; i != end; ++i) {
		TR << "Instantiating task operation " << (*i)["class"].to<std::string>() << " named " << i.skey() << std::endl;
		TaskOperationOP task( TaskOperationFactory::get_instance()->newTaskOperation( (*i)["class"].to<std::string>(), NULL) );
		//MMoverSP mmover ( mmoverfactory->from_string( (*i)["class"].to<std::string>() ) );
		task->parse_def( (*i) );
    TaskOperationSP tmpsp( task.get() );
    task.relinquish_ownership();
		tasks_.raw()[ i.skey() ] = tmpsp;
	}
	TR << "----------Finished Instantiating TaskOperations----------" << std::endl;
}

void BaseRole::instantiate_movers() {
	using namespace utility;
	using namespace utility::lua;
	using namespace protocols::moves;

	TR << "----------Instantiating Movers----------" << std::endl;

  luabind::globals(lstate_)["els"]["movers"] = luabind::newtable( lstate_ );
  movers_.raw( luabind::globals(lstate_)["els"]["movers"] );

	LuaObject dmovers( luabind::globals(lstate_)["els"]["dmovers"]);
	Factory<Mover> * moverfactory = Factory<Mover>::get_instance();
	for (LuaIterator i=dmovers.begin(), end; i != end; ++i) {
    TR << "Instantiating mover " << (*i)["class"].to<std::string>() << " named " << i.skey() << std::endl;
    if( moverfactory->has_string( (*i)["class"].to<std::string>() ) ) {
      MoverSP mover = moverfactory->from_string( (*i)["class"].to<std::string>() );
      mover->parse_def( (*i), scorefxns_, tasks_, mover_cache_ );
      luabind::globals(lstate_)["els"]["movers"][ i.skey() ] = mover;
    } else {
      MoverOP mover( MoverFactory::get_instance()->newMover( (*i)["class"].to<std::string>() ) );
      mover->parse_def( (*i), scorefxns_, tasks_, mover_cache_ );
      MoverSP tmpsp( mover.get() );
      mover.relinquish_ownership();
      luabind::globals(lstate_)["els"]["movers"][ i.skey() ] = tmpsp;
    }
	}
	TR << "----------Finished Instantiating Movers----------" << std::endl;
}

void BaseRole::instantiate_filters() {
	using namespace utility::lua;
	using namespace protocols::filters;

	TR << "----------Instantiating Filters----------" << std::endl;

  luabind::globals(lstate_)["els"]["filters"] = luabind::newtable( lstate_ );
  filters_.raw( luabind::globals(lstate_)["els"]["filters"] );

	LuaObject dfilters( luabind::globals(lstate_)["els"]["dfilters"]);
	for (LuaIterator i=dfilters.begin(), end; i != end; ++i) {
		TR << "Instantiating filter " << (*i)["class"].to<std::string>() << " named " << i.skey() << std::endl;
		FilterOP filter( FilterFactory::get_instance()->newFilter( (*i)["class"].to<std::string>() ) );
		filter->parse_def( (*i), scorefxns_, tasks_);
    FilterSP tmpsp( filter.get() );
    filter.relinquish_ownership();
    luabind::globals(lstate_)["els"]["filters"][ i.skey() ] = tmpsp;
	}
	TR << "----------Finished Instantiating Filters----------" << std::endl;
}
void BaseRole::instantiate_output() {
	using namespace utility::lua;
	using namespace protocols::outputter;

	TR << "----------Instantiating Outputters----------" << std::endl;

  luabind::globals(lstate_)["els"]["outputters"] = luabind::newtable( lstate_ );
  outputters_.raw( luabind::globals(lstate_)["els"]["outputters"] );

	LuaObject doutputters( luabind::globals(lstate_)["els"]["doutputters"] );
  utility::Factory<Outputter> * outputterfactory = utility::Factory<Outputter>::get_instance();
	for (LuaIterator i=doutputters.begin(), end; i != end; ++i) {
		TR << "Instantiating outputter " << (*i)["class"].to<std::string>() << " named " << i.skey() << std::endl;
		OutputterSP outputter ( outputterfactory->from_string( (*i)["class"].to<std::string>() ) );
    outputter->lregister( lstate_ );
		outputter->parse_def( (*i), tasks_ );
    luabind::globals(lstate_)["els"]["outputters"][ i.skey() ] = outputter;
	}
	TR << "----------Finished Instantiating Outputters----------" << std::endl;
}
void BaseRole::instantiate_inputters() {
	using namespace utility::lua;
	using namespace protocols::inputter;

	TR << "----------Instantiating Inputters----------" << std::endl;

  luabind::globals(lstate_)["els"]["inputters"] = luabind::newtable( lstate_ );
  inputters_.raw( luabind::globals(lstate_)["els"]["inputters"] );

	LuaObject dinputters( luabind::globals(lstate_)["els"]["dinputters"] );
  utility::Factory<Inputter> * inputterfactory = utility::Factory<Inputter>::get_instance();
	for (LuaIterator i=dinputters.begin(), end; i != end; ++i) {
		TR << "Instantiating inputter " << (*i)["class"].to<std::string>() << " named " << i.skey() << std::endl;
		InputterSP inputter ( inputterfactory->from_string( (*i)["class"].to<std::string>() ) );
    inputter->lregister( lstate_ );
		inputter->parse_def( (*i), tasks_, inputters_ );
    luabind::globals(lstate_)["els"]["inputters"][ i.skey() ] = inputter;
	}

	TR << "----------Finished Instantiating Inputters----------" << std::endl;
}

void BaseRole::instantiate_inputterstream() {
	using namespace utility::lua;
	using namespace protocols::inputter;
	TR << "----------Adding Inputters to InputterStream----------" << std::endl;
	LuaObject dinputterstream( luabind::globals(lstate_)["els"]["dinputterstream"]);
	for (LuaIterator i=dinputterstream.begin(), end; i != end; ++i) {
		inputterstream_->add_inputter( inputters_[ (*i).to<std::string>() ].to<InputterSP>() );
	}
  if( inputterstream_->size() == 0 ) {
    LuaObject dinputters( luabind::globals(lstate_)["els"]["dinputters"] );
    // inputterstream wasn't defined, so add all inputters to it
    for (LuaIterator i=dinputters.begin(), end; i != end; ++i) {
      inputterstream_->add_inputter( inputters_[ i.skey() ].to<InputterSP>() );
    }
  }
	TR << "----------Finished adding Inputters to InputterStream----------" << std::endl;
}

void BaseRole::update_mover_cache_mem() {
  if( mover_cache_length_ != mover_cache_->size() ) {
    // cache_mem_ is out of date
    // actually could be out of data even if size is same, but w.e, someone else can write a wrapper map class that keeps track of real memory usage while maintaining performance
#ifdef SERIALIZATION
/*
      std::stringstream s;
      core::io::serialization::toBinary(s, mover_cache_);
      mover_cache_mem_ =  s.str().length();
      mover_cache_length_ = mover_cache_->size();
*/
      mover_cache_mem_ =  0;
      mover_cache_length_ = mover_cache_->size();
#else
      TR << "Memory usage tracked only if compiled with serialization support" << std::endl;
#endif
  }
}
void BaseRole::register_calculators() {
  // this is so stupid, why is this not like all the other factories and register at init
  // and furthermore wtf is a calculator and how is it different from a filter
  // should totally be reimplemented
	core::Size const chain1( 1 ), chain2( 2 );
	using namespace core::pose::metrics;

	if( !CalculatorFactory::Instance().check_calculator_exists( "sasa_interface" ) ){
		PoseMetricCalculatorOP int_sasa_calculator = new core::pose::metrics::simple_calculators::InterfaceSasaDefinitionCalculator( chain1, chain2 );
		CalculatorFactory::Instance().register_calculator( "sasa_interface", int_sasa_calculator );
	}

	if( !CalculatorFactory::Instance().check_calculator_exists( "sasa" ) ){
		PoseMetricCalculatorOP sasa_calculator = new core::pose::metrics::simple_calculators::SasaCalculatorLegacy();
		CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );
	}
	if( !CalculatorFactory::Instance().check_calculator_exists( "ligneigh" ) ){
		PoseMetricCalculatorOP lig_neighbor_calc = new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator( chain1, chain2 );
  	CalculatorFactory::Instance().register_calculator( "ligneigh", lig_neighbor_calc );
	}
	if( !CalculatorFactory::Instance().check_calculator_exists( "liginterfE" ) ){
  	PoseMetricCalculatorOP lig_interf_E_calc = new core::pose::metrics::simple_calculators::InterfaceDeltaEnergeticsCalculator( "ligneigh" );
  	CalculatorFactory::Instance().register_calculator( "liginterfE", lig_interf_E_calc );
	}
}

void BaseRole::reparse_def( std::string const & type, std::string const & name ) {
	using namespace utility::lua;
  if( type == "tasks" ) {
    using namespace core::pack::task::operation;
    LuaObject dtasks( luabind::globals(lstate_)["els"]["dtasks"]);
    for (LuaIterator i=dtasks.begin(), end; i != end; ++i) {
      if( i.skey() == name ) {
        tasks_[ i.skey() ].to<TaskOperationSP>()->parse_def( (*i) );
        break;
      }
    }
  }
  if( type == "movers" ) {
    using namespace protocols::moves;
    LuaObject dmovers( luabind::globals(lstate_)["els"]["dmovers"]);
    for (LuaIterator i=dmovers.begin(), end; i != end; ++i) {
      if( i.skey() == name ) {
        movers_[ i.skey() ].to<MoverSP>()->parse_def( (*i), scorefxns_, tasks_, mover_cache_ );
        break;
      }
    }
  }
}

} //elscripts
} //protocols
#endif
