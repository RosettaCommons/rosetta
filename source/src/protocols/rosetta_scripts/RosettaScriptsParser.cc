// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rosetta_scripts/RosettaScriptsParser.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Interface base class Parser
/// @author Sarel Fleishman sarelf@u.washington.edu

#include <protocols/jd2/Job.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>

// Package headers
#include <protocols/jd2/parser/DataLoader.hh>
#include <protocols/jd2/parser/DataLoaderFactory.hh>

// Project headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>

#include <ObjexxFCL/string.functions.hh>
#include <basic/options/option.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/filters/FilterFactory.hh>
#include <protocols/filters/BasicFilters.hh>

#include <protocols/environment/Environment.hh>
#include <protocols/environment/ClaimingMover.hh>
#include <protocols/environment/EnvBeginMover.hh>
#include <protocols/environment/EnvEndMover.hh>

#include <utility/tag/Tag.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/metrics/CalculatorFactory.hh>

// movers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/moves/NullMover.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH


//scoring grids
//Pose Metric Calculators for filters
#include <core/pose/metrics/simple_calculators/InterfaceSasaDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <core/pose/metrics/simple_calculators/InterfaceDeltaEnergeticsCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>

// constraint types

// Utility headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <set>

// option key includes

#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/datacache/DataMap.hh>

#include <utility/vector0.hh>
#include <basic/Tracer.hh>


namespace protocols {
namespace rosetta_scripts {

using namespace core;
  using namespace basic::options;
  using namespace scoring;
using namespace moves;
using namespace jd2;

static basic::Tracer TR( "protocols.rosetta_scripts.RosettaScriptsParser" );

RosettaScriptsParser::RosettaScriptsParser() :
	Parser()
	//,ddfilter_factory_( new DockDesignFilterFactory )
{
  register_factory_prototypes();
}

RosettaScriptsParser::~RosettaScriptsParser(){}

typedef utility::tag::TagCOP TagCOP;
typedef utility::vector0< TagCOP > TagCOPs;

/// @details Uses the Tag interface to the xml reader library in boost to parse an xml file that contains design
/// protocol information. A sample protocol file can be found in src/pilot/apps/sarel/dock_design.protocol.
/// SCOREFXNS provides a way to define scorefunctions as they are defined in the rosetta database, using the
/// weights/patch convenctions. Several default scorefunctions are preset and can be used without defining them
/// explicitly.
/// FILTERS defines a set of filters that can be used together with the dockdesign movers to prune out poses that
/// don't meet certain criteria
/// MOVERS defines the movers that will be used
/// PROTOCOLS is the only order-sensitive section where subsequent movers and filters are expected to be defined.
/// These movers and filters were defined in the previous two sections. The order in which the protocols are
/// specified by the user will be maintained by the DockDesign mover.
/// APPLY_TO_POSE This section allows for certain movers to be applied to the pose prior to starting the DockDesign
/// protocol. For instance, these movers could set constraints, such as favor_native_residue. In this case, for
/// example, the weights of res_type_constraint in all of the scorefunctions that are defined in SCOREFXNS or by
/// default are set to 1.0, but can be changed by the user globally (in the definition of the weight on the constraint),
/// or in particular for each of the scorefunctions by changing the relevant term (which is set by default to the
/// global value).
/// OUTPUT is a section which allows the XML control of how things are output
/// Notice that the order of the sections by which the protocol is written doesn't matter, BUT the order of the
/// mover-filter pairs in PROTOCOLS section does matter.
bool
RosettaScriptsParser::generate_mover_from_pose( JobCOP, Pose & pose, MoverOP & in_mover, bool new_input, std::string const xml_fname ){

  bool modified_pose( false );

  if( !new_input ) return modified_pose;

  std::string const dock_design_filename( xml_fname == "" ? option[ OptionKeys::parser::protocol ] : xml_fname );
  TR << "dock_design_filename=" << dock_design_filename << std::endl;
  utility::io::izstream fin;
  fin.open(dock_design_filename.c_str() );
  if( !fin.good() ){
    utility_exit_with_message("Unable to open Rosetta Scripts XML file: '" + dock_design_filename + "'.");
  }

  TagCOP tag;
  if( option[ OptionKeys::parser::script_vars ].user() ) {
    std::stringstream fin_sub;
    substitute_variables_in_stream(fin, option[ OptionKeys::parser::script_vars ], fin_sub);
    tag = utility::tag::Tag::create( fin_sub );
    }
  else {
    tag = utility::tag::Tag::create(fin);
  }

  fin.close();
  TR << "Parsed script:" << "\n";
  TR << tag;
  TR.flush();
  runtime_assert( tag->getName() == "dock_design" || tag->getName() == "ROSETTASCRIPTS" );

  in_mover = generate_mover_for_protocol(pose, modified_pose, tag);

  return modified_pose;
}

MoverOP RosettaScriptsParser::parse_protocol_tag(Pose & pose, utility::tag::TagCOP protocol_tag)
{
  bool modified_pose = false;

  MoverOP mover =  generate_mover_for_protocol(pose, modified_pose, protocol_tag);

  if (modified_pose)
  {
    utility_exit_with_message("parse_protocol_tag resulted in modified_pose");
  }

  return mover;
}

MoverOP RosettaScriptsParser::parse_protocol_tag(TagCOP protocol_tag)
{
  Pose temp_pose;
  bool modified_pose = false;

  MoverOP mover =  generate_mover_for_protocol(temp_pose, modified_pose, protocol_tag);

  if (modified_pose)
  {
    utility_exit_with_message("parse_protocol_tag resulted in modified_pose");
  }

  return mover;
}

MoverOP RosettaScriptsParser::generate_mover_for_protocol(Pose & pose, bool & modified_pose, TagCOP tag)
{
  protocols::rosetta_scripts::ParsedProtocolOP protocol( new protocols::rosetta_scripts::ParsedProtocol );

  Movers_map movers;
  protocols::filters::Filters_map filters;
  basic::datacache::DataMap data; // abstract objects, such as scorefunctions, to be used by filter and movers

  MoverOP mover;

  typedef std::pair< std::string const, MoverOP > StringMover_pair;
  typedef std::pair< std::string const, protocols::filters::FilterOP > StringFilter_pair;
  //typedef std::pair< std::string const, ScoreFunctionOP > StringScorefxn_pair;
//setting up some defaults
  protocols::filters::FilterOP true_filter = new protocols::filters::TrueFilter;
  protocols::filters::FilterOP false_filter = new protocols::filters::FalseFilter;
  filters.insert( StringFilter_pair( "true_filter", true_filter ) );
  filters.insert( StringFilter_pair( "false_filter", false_filter ) );

  MoverOP null_mover = new protocols::moves::NullMover();
  movers.insert( StringMover_pair( "null", null_mover) );

// default scorefxns
  ScoreFunctionOP commandline_sfxn = core::scoring::getScoreFunction();
  ScoreFunctionOP talaris2013 = ScoreFunctionFactory::create_score_function(TALARIS_2013);
  ScoreFunctionOP score12 = ScoreFunctionFactory::create_score_function( PRE_TALARIS_2013_STANDARD_WTS, SCORE12_PATCH );
  ScoreFunctionOP docking_score = ScoreFunctionFactory::create_score_function( PRE_TALARIS_2013_STANDARD_WTS, DOCK_PATCH );
  ScoreFunctionOP soft_rep = ScoreFunctionFactory::create_score_function( SOFT_REP_DESIGN_WTS );
  ScoreFunctionOP docking_score_low = ScoreFunctionFactory::create_score_function( "interchain_cen" );
  ScoreFunctionOP score4L = ScoreFunctionFactory::create_score_function( "cen_std", "score4L" );

  data.add( "scorefxns", "commandline", commandline_sfxn );
  data.add( "scorefxns", "talaris2013", talaris2013 );
  data.add( "scorefxns", "score12", score12 );
  data.add( "scorefxns", "score_docking", docking_score );
  data.add( "scorefxns", "soft_rep", soft_rep );
  data.add( "scorefxns", "score_docking_low", docking_score_low );
  data.add( "scorefxns", "score4L", score4L );
  //default scorefxns end

  /// Data Loaders
  std::set< std::string > non_data_loader_tags;
  non_data_loader_tags.insert( "MOVERS" );
  non_data_loader_tags.insert( "APPLY_TO_POSE" );
  non_data_loader_tags.insert( "FILTERS" );
  non_data_loader_tags.insert( "PROTOCOLS" );
  non_data_loader_tags.insert( "OUTPUT" );

  /// Load in data into the basic::datacache::DataMap object.  All tags beside those listed
  /// in the non_data_loader_tags set are considered DataLoader tags.
  TagCOPs const all_tags = tag->getTags();
  for ( Size ii = 0; ii < all_tags.size(); ++ii ) {
    using namespace parser;
    TagCOP iitag = all_tags[ ii ];
    if ( non_data_loader_tags.find( iitag->getName() ) != non_data_loader_tags.end() ) continue;
    DataLoaderOP loader = DataLoaderFactory::get_instance()->newDataLoader( iitag->getName() );
    loader->load_data( pose, iitag, data );
  }
  if ( !tag->hasTag("MOVERS") )
    utility_exit_with_message("parser::protocol file must specify MOVERS section");
  if ( !tag->hasTag("PROTOCOLS") )
      utility_exit_with_message("parser::protocol file must specify PROTOCOLS section");

  foreach( TagCOP const curr_tag, all_tags ){
    ///// APPLY TO POSE
    if ( curr_tag->getName() == "APPLY_TO_POSE" ) { // section is not mandatory
      /// apply to pose may affect all of the scorefxn definitions below, so it is called first.
      TagCOPs const apply_tags( curr_tag->getTags() );

      foreach(TagCOP apply_tag_ptr, apply_tags){
        std::string const mover_type( apply_tag_ptr->getName() );
        MoverOP new_mover( MoverFactory::get_instance()->newMover( apply_tag_ptr, data, filters, movers, pose ) );
        runtime_assert( new_mover );
        new_mover->apply( pose );
        TR << "Defined and applied mover of type " << mover_type << std::endl;
        bool const name_exists( movers.find( mover_type ) != movers.end() );
        if ( name_exists ) {
          utility_exit_with_message( "Can't apply_to_pose the same mover twice" + mover_type );
        }

        modified_pose = true;
      } /// done applyt_tag_it

    }//  fi apply_to_pose


    TR.flush();

    ////// Filters
    if ( curr_tag->getName() == "FILTERS" ) {

      foreach(TagCOP tag_ptr, curr_tag->getTags()){
        std::string const type( tag_ptr->getName() );
        if ( ! tag_ptr->hasOption("name") )
          utility_exit_with_message("Can't define unnamed Filter of type " + type );
        std::string const user_defined_name( tag_ptr->getOption<std::string>("name") );
        bool const name_exists( filters.find( user_defined_name ) != filters.end() );
        if ( name_exists ) {
          TR.Error << "ERROR filter of name \"" << user_defined_name << "\" (with type " << type
            << ") already exists. \n" << tag << std::endl;
          utility_exit_with_message("Duplicate definition of Filter with name " + user_defined_name);
        }
        protocols::filters::FilterOP new_ddf( protocols::filters::FilterFactory::get_instance()->newFilter( tag_ptr, data, filters, movers, pose ) );
        runtime_assert( new_ddf );
        filters.insert( std::make_pair( user_defined_name, new_ddf ) );
        TR << "Defined filter named \"" << user_defined_name << "\" of type " << type << std::endl;
      }
    }

    TR.flush();

    ////// MOVERS
    if( curr_tag->getName() == "MOVERS" ){
      foreach(TagCOP tag_ptr, curr_tag->getTags()){
        std::string const type( tag_ptr->getName() );
        if ( ! tag_ptr->hasOption("name") )
          utility_exit_with_message("Can't define unnamed Mover of type " + type );
        std::string const user_defined_name( tag_ptr->getOption<std::string>("name") );
        bool const name_exists( movers.find( user_defined_name ) != movers.end() );
        if ( name_exists ) {
          TR.Error << "ERROR mover of name \"" << user_defined_name << "\" (with type " << type
            << ") already exists.\n" << tag << std::endl;
          utility_exit_with_message("Duplicate definition of Mover with name " + user_defined_name);
        }
        /// APL -- singleton factory replacement.
        MoverOP new_mover( MoverFactory::get_instance()->newMover( tag_ptr, data, filters, movers, pose ) );
        runtime_assert( new_mover );
        movers.insert( std::make_pair( user_defined_name, new_mover ) );
        TR << "Defined mover named \"" << user_defined_name << "\" of type " << type << std::endl;
      }
    }//fi MOVERS

    TR.flush();
  }// foreach curr_tag

  ////// ADD MOVER FILTER PAIRS
  TagCOP const protocols_tag( tag->getTag("PROTOCOLS") );
  protocol->parse_my_tag( protocols_tag, data, filters, movers, pose );
  TR.flush();

  ////// Set Output options
  if ( tag->hasTag("OUTPUT") ) {
    TagCOP const output_tag( tag->getTag("OUTPUT") );
    protocol->final_scorefxn( rosetta_scripts::parse_score_function( output_tag, data, "commandline" ) );
  } else {
    protocol->final_scorefxn( commandline_sfxn );
  }

  tag->die_for_unaccessed_options_recursively();

  return protocol;
}

///@brief Create a variable substituted version of the input stream, given a StringVectorOption formated list of variables
///to substitiute. Each item in script_vars should be in the form of "variable=value", where "value" is the string to substitiute
///into the input stream whereever the string "%%variable%%" is found in the input.
///
///Having the return value be passed through a parameter is to get around some copy-constructor limitations of C++ streams.
void
RosettaScriptsParser::substitute_variables_in_stream( std::istream & in, utility::options::StringVectorOption const& script_vars, std::stringstream & out){
	using namespace std;
	//Parse variable substitutions
	map<string,string> var_map;
	for( Size ii = 1; ii <= script_vars.size(); ii++ ) {
		Size eq_pos(script_vars[ii].find('='));
		if(eq_pos != string::npos) { //ignore ones without a '='
			var_map[ script_vars[ii].substr(0,eq_pos) ] = script_vars[ii].substr(eq_pos+1);
		}
	}
	//Print parsing for confirmation
	TR << "Variable substitution will occur with the following values: ";
	for( map<string,string>::const_iterator map_it( var_map.begin() ), map_end( var_map.end() );
			 map_it != map_end; ++map_it ){
		TR << "'%%" << map_it->first << "%%'='" << map_it->second << "';  ";
	}
	TR << std::endl;
	var_map[""] = "%%"; //for %%%%->%% substitutions
	//Perform actual variable substitution
	TR << "Substituted script:" << "\n";
	string line;
	while( getline( in, line ) ) {
		Size pos, end, last(0);
		string outline; // empty string
		while( ( pos = line.find("%%", last) ) != string::npos) {
			end = line.find("%%", pos+2);
			if(end == string::npos) break; //Unmatched %% on line - keep as is.
			outline += line.substr( last, pos-last );
			last = end+2;
			string var( line.substr( pos+2, end-(pos+2) ) );
			if( var_map.find( var ) != var_map.end() ) {
				outline += var_map[var];
			}
			else {
				outline += "%%" + var + "%%"; // Silently ignore missing values - will probably cause error later.
			}
		}
		outline += line.substr( last ); // Add the rest of the line, starting at last
		TR << outline << "\n";
		out << outline << "\n";
	}
	TR.flush();
}

///@brief Factories avoid requiring compiler dependencies for all possible constructible derived classes,
///by allowing runtime registration of derived class prototypes. However, this requires
///pre-registration of a derived type with the factory prior to asking the factory to return an
///instance of that type. This method registers those additional derived classes that are available for
///construction in the RosettaScriptsParser context.
/// TO-DO: replace this manual factory registration system with a load-time-initialized singleton scheme (see r32404 for example)
void
RosettaScriptsParser::register_factory_prototypes()
{
	// note: TaskOperations are now registered with a singleton factory at load time using apl's creator/registrator scheme

	// also register some constraint types with the ConstraintFactory (global singleton class)
	// this allows derived non-core constraints to be constructed from string definitions in constraints files
	//using namespace core::scoring::constraints;
	//ConstraintFactory & cstf( ConstraintIO::get_cst_factory() );
	//cstf.add_type( new core::scoring::constraints::SequenceProfileConstraint(
	//	Size(), utility::vector1< id::AtomID >(), NULL ) );

	// register calculators
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

}//jd2
}//protocols
