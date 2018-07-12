// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/sergey/hal_demo.cc
/// @brief HAL client for executing XML-constructible Movers
/// @author Sergey Lyskov

#include <utility/json_utilities.hh>


#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)

#include <devel/init.hh>

#include <basic/options/option.hh>

#include <protocols/network/hal.hh>
#include <protocols/network/util.hh>
#include <protocols/network/ui_mover.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/datacache/DataMap.hh>

#include <protocols/moves/MoverFactory.hh>


#include <basic/Tracer.hh>
#include <basic/options/option_macros.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/xsd_util/util.hh>
#include <utility/tag/Tag.hh>

#include <json.hpp>

#include <algorithm>

static basic::Tracer TR("hal");

using std::string;
using namespace protocols::network;

OPT_KEY( StringVector, filter)

struct MoverParameter
{
	string name, type, default_, description;
};

using MoverParameters = std::vector<MoverParameter>;

/// adapted from file:utility/xsd_util/util.cc output_all_tag_options
MoverParameters collect_tag_info(xmlNode* node)
{
	using namespace utility::xsd_util;

	//options << (level > 1 ? "\n" : "" ) << "\"" << tagname << "\" ";
	//for ( platform::Size i=1; i<level; ++i ) options << "sub-";
	//options << "tag:";

	// if ( level > 1 ) {
	//  options << " ";
	//  generate_human_readable_documentation(node, options);
	//  options << "\n";
	// } else {
	//  options << "\n\n";
	// }

	MoverParameters params;

	for ( xmlNode* subnode = node->children; subnode != nullptr; subnode = subnode->next ) { //Loop through sub-nodes.
		if ( subnode->type == XML_ELEMENT_NODE && !strcmp(reinterpret_cast<const char*>(subnode->name),"attribute") ) {
			std::string const name( get_node_option( subnode, "name" ) );
			std::string const type( get_type_name( get_node_option( subnode, "type" ) ) );
			std::string const default_( get_node_option(subnode, "default") );
			//tags << " " << name << "=(" << type << (def.empty() ? "" : ",\"" + def + "\"" ) << ")";

			//options << "\t" << name << " (" << type << (def.empty() ? "" : ",\"" + def + "\"" ) << "):  ";
			std::stringstream description;
			generate_human_readable_documentation( subnode, description);

			//TR << "\ntag:\t" << name << "\ntype:\t" << type << "\ndefault:\t" << default_ << "\ndescription:\t" << description.str() << std::endl;

			params.emplace_back( MoverParameter{name, type, default_, description.str() } );
		}
	}
	return params;
}

/// adapted from file:utility/xsd_util/util.cc generate_human_readable_recursive
MoverParameters get_element_info_from_xml_node(xmlNode* rootnode, platform::Size /* level */, std::string const &complextype, std::string const &tag_name_to_print)
{
	using namespace utility::xsd_util;

	MoverParameters params;

	for ( xmlNode* cplxtype_node = rootnode->children; cplxtype_node!=nullptr; cplxtype_node = cplxtype_node->next ) {
		if ( cplxtype_node->type == XML_ELEMENT_NODE && !strcmp(reinterpret_cast<const char*>(cplxtype_node->name),"complexType") ) {

			if ( complextype.empty() ||  get_node_option(cplxtype_node, "name") == complextype ) {
				// for ( platform::Size i=1; i<level; ++i ) { tags << "\t"; }

				// tags << "<";
				std::string const tagname( tag_name_to_print.empty() ? get_tag_name(cplxtype_node) : tag_name_to_print );
				// tags << tagname;
				params = collect_tag_info(cplxtype_node);
				// tags << ">\n";

				// //If this is level 1, generate the documentation information.
				// if ( level==1 ) {
				//  for ( xmlNode* annotation_node = cplxtype_node->children; annotation_node!=nullptr; annotation_node = annotation_node->next ) {
				//   if ( !strcmp(reinterpret_cast<const char*>(annotation_node->name),"annotation") ) {
				//    generate_human_readable_documentation( annotation_node, description );
				//   }
				//  }
				// }

				// //Get sub-tags.
				// for ( xmlNode* element_node = cplxtype_node->children; element_node != nullptr; element_node = element_node->next ) {
				//  if ( element_node->type == XML_ELEMENT_NODE ) {
				//   if ( !strcmp(reinterpret_cast<const char*>(element_node->name),"element") ) {
				//    std::string type_name( get_node_option( element_node, "type" ) );
				//    generate_human_readable_recursive( type_name.empty() ? element_node : rootnode, description, tags, options, level+1, type_name, type_name.empty() ? get_node_option( element_node, "name" ) : "" );
				//   } else if ( !strcmp(reinterpret_cast<const char*>(element_node->name),"choice") ) {
				//    for ( xmlNode* element_node2 = element_node->children; element_node2 != nullptr; element_node2 = element_node2->next ) {
				//     if ( element_node2->type == XML_ELEMENT_NODE && !strcmp(reinterpret_cast<const char*>(element_node2->name),"element") ) {
				//      std::string type_name( get_node_option( element_node2, "type" ) );
				//      generate_human_readable_recursive( type_name.empty() ? element_node2 : rootnode, description, tags, options, level+1, type_name, type_name.empty() ? get_node_option( element_node2, "name" ) : "" );
				//     }
				//    }
				//   }
				//  }
				// }

				// for ( platform::Size i=1; i<level; ++i ) { tags << "\t"; }
				// tags << "</" << tagname << ">\n";
			}
		}
	}
	return params;
}


/// adapted from file:utility/xsd_util/util.cc generate_human_readable_summary
MoverParameters get_mover_info(std::string const &mover_name, std::string const &xsd)
{
	//bool const has_name_and_type( !component_name.empty() && !component_type.empty() );
	//std::string const name_and_type( has_name_and_type ? component_type + "_" + component_name + "_type" : "" );
	std::string const name_and_type( "mover_" + mover_name + "_type" );

	// std::stringstream description(""), tags(""), options("");
	// description << "DESCRIPTION:\n\n";
	// tags << "USAGE:\n\n";
	// options << "OPTIONS:\n\n";

	xmlDoc* doc( xmlReadMemory( xsd.c_str(), xsd.length()+1, nullptr, nullptr, 0 ) );

	MoverParameters params( get_element_info_from_xml_node( xmlDocGetRootElement(doc), 1, name_and_type, "" ) );

	xmlFreeDoc(doc); doc=nullptr;
	xmlCleanupParser();

	return params;
}


// bool is_gui_constructable(MoverParameters const & params)
// {
//  //auto known_types = {"bool", "int", "real", "string"};
//  for(auto const & p : params) {
//   //auto it = std::find( known_types.begin(), known_types.end(), p.type);
//   auto it = type_mapping.find(p.type);
//   if( it == type_mapping.end() ) {
//    TR << "Unknown parameter type: " << TR.Red << p.type << std::endl;
//    TR << "           description: " << p.description << std::endl;
//    return false;
//   }
//  }
//  return true;
// }


// provide type mapping RosettaScripts -> UI specification
std::map<string, json> const type_mapping = {
{"bool",                  json( { {_f_type_, _t_boolean_}, } ) },  // {"default", true}, {"optional", true},

{"int",                   json( { {_f_type_, _t_integer_}, } ) },
{"port_range",            json( { {_f_type_, _t_integer_}, {"min", 1024}, {"max", 65535}, } ) },
{"max_packet_size_range", json( { {_f_type_, _t_integer_}, {"min", 1},    {"max", 65535}, } ) },

{"real",                  json( { {_f_type_, _t_float_}, } ) },

{"string",                json( { {_f_type_, _t_string_}, } ) },
// {"int_cslist", "string"},

//{"minimizer_type",        json( { {_f_type_, _t_string_}, } ) },

};


/// return JSON representation for parameter specification, return `null` JSON object if type could not be adequately represented
json convert_to_json(MoverParameter const &p)
{
	auto it = type_mapping.find(p.type);
	if ( it == type_mapping.end() ) {
		// TR << "Unknown parameter type: " << TR.Red << p.type << std::endl;
		// TR << "           description: " << p.description << std::endl;
		return nullptr;
	} else return it->second;
}



///  generate specification for function with given argument list, return `null` JSON object if params contain incompatible types
json generate_function_specification(MoverParameters const &params)
{
	std::map<string, std::function< json (string const&) > > const converters = {
		{ _t_boolean_, [](std::string const &v) { return json( utility::is_true_string(v) ); } },
		{_t_integer_,  [](std::string const &v) { return json( std::stoi(v) ); } },
		{_t_float_,    [](std::string const &v) { return json( std::stod(v) ); } },
		{_t_string_,   [](std::string const &v) { return json( v ); } },
		};

	json f = json::object();

	for ( auto const & param : params ) {
		json arg = convert_to_json(param);
		if ( arg.is_null() ) {
			if ( not param.default_.empty() ) continue;  // if parameter could not be represented in UI but have deault value we can still build function specification
			f = nullptr;
			break;
		}

		auto const banned_arguments = {_f_name_};
		auto it = std::find( banned_arguments.begin(), banned_arguments.end(), param.name);

		if ( it == banned_arguments.end() ) {
			f[param.name] = arg;

			if ( not param.default_.empty() ) f[param.name][_f_default_] = converters.at(arg[_f_type_])(param.default_);

			f[param.name][_f_description_] = param.description;
		}
	}

	f[_f_pose_] = { {_f_type_, _t_pose_}, };

	return f;
}



/// Generate HAL specification
string specification()
{
	std::map<string, MoverParameters> movers;

	auto const banned_movers = {"MetropolisHastings", "MultipleOutputWrapper", "MultiplePoseMover", "ScriptCM"};

	protocols::moves::MoverFactory* mover_factory( protocols::moves::MoverFactory::get_instance() );
	auto const & mover_creators = mover_factory->mover_creator_map();
	for ( auto const & mc : mover_creators ) {
		string const & name = mc.first;

		//TR << name << std::endl;

		auto it = std::find( banned_movers.begin(), banned_movers.end(), name);
		if ( it == banned_movers.end() ) {

			if ( basic::options::option[basic::options::OptionKeys::filter]().size() ) {
				bool skip = true;
				for ( auto & s : basic::options::option[basic::options::OptionKeys::filter]() ) {
					if ( name.find(s) != std::string::npos ) { skip = false; break; }
				}
				if ( skip ) continue;
			}

			//if( name == "Small" or name == "MinMover"  or  name == "Docking"  or name == "PyMOLMover" ) {
			//if( name == "TopologyBrokerMover" ) {
			//TR << "Mover \"" << name << '"' << std::endl;;
			utility::tag::XMLSchemaDefinition xsd;
			mover_factory->provide_xml_schema( name, xsd );

			MoverParameters params = get_mover_info(name, xsd.full_definition() );

			// if( is_gui_constructable(params) ) {
			//  movers.emplace( std::make_pair(name, std::move(params) ) );
			//  TR << TR.Green << "Adding mover: " << name << TR.Reset << std::endl;
			// }
			movers.emplace( std::make_pair(name, std::move(params) ) );

			//TR << xsd.human_readable_summary( name, "mover" );
			//TR << "full definition: " << xsd.full_definition() << std::endl;;
			//TR << "_________" << std::endl;
		}
		// if( name == "MinMover" ) {
		//  utility::tag::XMLSchemaDefinition xsd;
		//  mover_factory->provide_xml_schema( name, xsd );
		//  TR << xsd.human_readable_summary( name, "mover" );
		//  TR << "full definition: " << xsd.full_definition() << std::endl;;
		//  TR << "_________" << std::endl;
		// }
	}

	json specification, functions;

	for ( auto const & mover : movers ) {
		json f = generate_function_specification(mover.second);
		if ( f.is_object() ) {
			//TR << TR.Green << "Adding mover: " << mover.first << TR.Reset << std::endl;
			functions[mover.first] = f;
		}
	}

	specification[_f_functions_] = functions;

	string r;
	nlohmann::json::basic_json::to_msgpack(specification, r);
	return r;

	// {
	//  nlohmann::json j;
	//  // j["functions"] = {
	//  //  {{ "name", "my_mover" },
	//  //  { "arguments",
	//  //  {
	//  //  {"param1", "int"},
	//  //  {"param2", "float"},
	//  //  },
	//  //  }},
	//  //  {{ "name", "test_mover2" },
	//  //  { "arguments",
	//  //  {
	//  //  {"param1", "string"},
	//  //  {"param2", "float"},
	//  //  },
	//  //  }},
	//  //  };
	//  json foo = json::object();
	//  //foo["arguments"] = json::object();
	//  //foo["arguments"]["argument_1"] = { {"type", "integer"}, {"optional", true}, {"default", 1}, };
	//  foo["residue_n"] = { {"type", "integer"}, {"optional", true}, {"default", 1}, {"min", -10}, {"max", 10}, };
	//  foo.emplace("delta angle", json::object( { {"type", "float"}, {"optional", true}, {"default", 1}, {"min", -1}, } ) );
	//  foo.emplace("magnitude",   json::object( { {"type", "float"}, {"optional", true}, {"default", 1}, {"min", -10}, {"max", -2}, } ) );
	//  json bar = json::object();
	//  j["functions"] = { {"foo_mover", foo}, {"bar_mover", bar} };
	//  string r;
	//  nlohmann::json::basic_json::to_msgpack(j, r);
	//  return r;
	// }
}


protocols::moves::MoverOP construct_mover_from_command(json const &j)
{
	//TR << "construct_mover_from_command:" << j.dump(2).c_str();

	string name = j.value(_f_name_, "");

	auto tag = std::make_shared<utility::tag::Tag>();
	tag->setName(name);

	for ( auto arg = j[_f_arguments_].begin(); arg != j[_f_arguments_].end(); ++arg ) {
		if ( arg.key() != _f_pose_ ) {
			auto v = arg.value();
			auto t = v.type();
			if ( t == json::value_t::boolean ) tag->setOption(arg.key(), static_cast<int>(v) );
			if ( t == json::value_t::number_integer  or  t == json::value_t::number_unsigned ) tag->setOption(arg.key(), static_cast<int>(v) );
			if ( t == json::value_t::number_float ) tag->setOption(arg.key(), static_cast<double>(v) );
		}
	}

	TR << "Tag:" << TR.Blue << *tag << TR.Reset << std::endl;

	basic::datacache::DataMap data_map;

	auto mover = protocols::moves::MoverFactory::get_instance()->newMover(tag, data_map, protocols::filters::Filters_map(), protocols::moves::Movers_map(), core::pose::Pose() );
	TR << TR.Blue;  mover->show(TR);  TR << TR.Reset << std::endl;

	return mover;

	// {
	//  string name = "Small";
	//  auto tag = std::make_shared<utility::tag::Tag>();
	//  tag->setName(name);
	//  //tag->setOption("address", "123.456.789.0");
	//  //tag->setOption("port", "1234");
	//  tag->setOption("nmoves", 10);
	//  //TR << "Tag:" << TR.Blue << *tag << TR.Reset << std::endl;
	//  basic::datacache::DataMap data_map;
	//  // protocols::filters::Filters_map filters_map;
	//  // protocols::moves::Movers_map mover_map;
	//  // core::pose::Pose pose;
	//  auto mover = protocols::moves::MoverFactory::get_instance()->newMover(tag, data_map, protocols::filters::Filters_map(), protocols::moves::Movers_map(), core::pose::Pose() );
	//  TR << TR.Blue;  mover->show(TR);  TR << TR.Reset << std::endl;
	// }
}

json hal_executioner(json const &command)
{
	//std::cout << "Rosetta: executing command: " << command << std::endl;
	std::cout << "Rosetta: executing command: " << command[_f_name_] << std::endl;

	nlohmann::json result;

	core::pose::PoseOP pose = protocols::network::bytes_to_pose(command[_f_arguments_].value(_f_pose_, "") );

	if ( pose ) {
		auto mover = construct_mover_from_command(command);
		if ( mover ) {

			//protocols::network::UIMover ui;
			protocols::network::AddUIObserver(*pose);
			mover->apply(*pose);
		}

		auto pose_binary = protocols::network::pose_to_bytes(*pose);

		result[_f_pose_] = pose_binary;
	}

	return result;
}

int main(int argc, char * argv [])
{

	try {
		NEW_OPT( filter, "Specify list of strings to use as `filter` for Mover names", "");

		devel::init(argc, argv);

		{ // creating dummy pose object to trigger database load so later we can create Pose immeditaly
			core::pose::Pose p;
			core::import_pose::pose_from_pdbstring(p, "ATOM     17  N   ILE A   1      16.327  47.509  23.466  1.00  0.00\n");
		}

		if ( basic::options::option[basic::options::OptionKeys::filter]().size() ) {
			TR << "Filtering Movers with: " << TR.Yellow;
			for ( auto & s : basic::options::option[basic::options::OptionKeys::filter]() ) TR << s << " ";
			TR << TR.Reset << std::endl;
		}

		protocols::network::hal(specification, hal_executioner, CommandLineArguments{argc, argv});
		return 0;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

#else // ! ( defined(ZEROMQ) and  defined(_NLOHMANN_JSON_ENABLED_) )

#include <utility/excn/Exceptions.hh>
#include <devel/init.hh>
#include <iostream>

int main(int argc, char * argv [])
{
	try {
		devel::init(argc, argv);

		std::cerr << "HAL app need to be build with extras=serialization! Aborting..." << std::endl;
		return 1;
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

#endif // defined(ZEROMQ) and  defined(_NLOHMANN_JSON_ENABLED_)
