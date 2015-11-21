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

#ifndef INCLUDED_app_javierbq_micro_utils_hh
#define INCLUDED_app_javierbq_micro_utils_hh

#include <core/chemical/AA.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/rms_util.hh>
#include <protocols/moves/DsspMover.hh>

// Utility headers
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include <jsoncpp/json/json.h>
#include <jsoncpp/json/value.h>

#include <curlpp/cURLpp.hpp>
#include <curlpp/Easy.hpp>
#include <curlpp/Options.hpp>
#include <curlpp/Infos.hpp>

#define foreach BOOST_FOREACH


namespace cyrus {
namespace micro {
namespace utils {

static THREAD_LOCAL basic::Tracer TR( "cyrus_micro_utils" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

core::scoring::ScoreFunctionOP parse_score_function(const Json::Value& json);

core::pose::Pose get_pose_from_json(const Json::Value& json)
{
	using namespace core::io::silent;
	std::string pdb_data = json.get("pdb_data","").asString();
	std::string pose_data = json.get("pose","").asString();

	core::pose::Pose pose;
	if(!pose_data.empty() && pdb_data.empty())
	{
		TR << "loading pose from silent structure" << std::endl;
		std::istringstream data;
		data.str(pose_data);
		SilentFileData silent_data;
		utility::vector1<std::string> tags;
		tags.push_back("serialized_pose");
		silent_data.read_stream(data, tags, true);
		silent_data.get_structure("serialized_pose").fill_pose(pose);

	}
	else if(pose_data.empty() && !pdb_data.empty())
	{
		TR << "loading pose from pdb" << std::endl;
		core::import_pose::pose_from_pdbstring(pose, pdb_data);
	}
	else if(pose_data.empty() && pdb_data.empty())
	{
		TR.Error << "You need to specify an input source!" << std::endl;
		throw std::runtime_error("Missing pose input in json");
	}
	else
	{
		TR.Error << "Too much input! pass only pdb_data or pose in json" << std::endl;
		throw std::runtime_error("Ambiguous pose input in json");
	}

	return pose;
}

void
dump_secstruct(
               pose::Pose const & pose,
               std::ostream & out
               ) {
    
    Size n_helix = 0, n_sheet = 0, n_loop = 0;
    for (Size ires=1; ires<pose.n_residue(); ++ires) {
        char secstruct=pose.secstruct()[ires-1];
        Size jres = ires;
        while (jres<pose.n_residue() && pose.secstruct()[jres-1] == secstruct) {
            jres += 1;
        }
        if (secstruct == 'E') {
            n_sheet += 1;
        }
        ires = jres;
    }

    char const chainids []="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890";
    Size i_sheet = 0;
    for (Size ires=1; ires<pose.n_residue(); ++ires) {
        char secstruct=pose.secstruct()[ires-1];
        Size jres = ires;
        while (jres<pose.n_residue() && pose.secstruct()[jres-1] == secstruct) {
            jres += 1;
        }
        
        using namespace ObjexxFCL::format;
        std::string sec_str;
        if (secstruct == 'H') {
            n_helix += 1;
            sec_str = "HELIX " + I(4, n_helix) + " H" + I(2,2,n_helix)
            + " " + pose.residue(ires).name3() +" "+ chainids[pose.residue(ires).chain()-1] + I(5, ires) + " "
            + " " + pose.residue(jres).name3() +" "+ chainids[pose.residue(jres).chain()-1] + I(5, jres) + " "
            + "01                               " + I(5,jres-ires+1); // assuming standard helix for now
        }
        else if (secstruct == 'E') {
            i_sheet += 1;
            sec_str = "SHEET " + I(4, i_sheet) + " B" + I(2,2,i_sheet) + I(2,n_sheet)
            + " " + pose.residue(ires).name3() +" "+ chainids[pose.residue(ires).chain()-1] + I(4, ires) + " "
            + " " + pose.residue(jres).name3() +" "+ chainids[pose.residue(jres).chain()-1] + I(4, jres) + " "
            + " 0";
        }
        else if (secstruct == 'L') {
            n_loop += 1;
            sec_str = "TURN  " + I(4, n_loop) + " T" + I(2,2,n_loop)
            + " " + pose.residue(ires).name3() +" "+ chainids[pose.residue(ires).chain()-1] + I(4, ires) + " "
            + " " + pose.residue(jres).name3() +" "+ chainids[pose.residue(jres).chain()-1] + I(4, jres) + " ";
        }
        out << sec_str << std::endl;
        
        ires = jres;
    }
}

void add_pose_to_result_value(const Json::Value& task, Json::Value& result, core::pose::Pose & pose, core::pose::PoseOP input_pose = 0)
{
	using namespace core::io::silent;
	core::scoring::ScoreFunctionOP scorefxn = parse_score_function(task);
	core::Real score(scorefxn->score(pose));
	std::ostringstream data;
	Json::Value pose_object;
	protocols::moves::DsspMover dssp;
	dssp.apply(pose);
	pose_object["secstruct"] = pose.secstruct();

        dump_secstruct(pose, data);

	if(input_pose)
	{
		core::Real rmsd = core::scoring::calpha_superimpose_pose(pose, *input_pose);
		TR.Debug << "RMSD: " <<rmsd <<  std::endl << std::flush;
		pose_object["rms"] = rmsd;
	}
	bool return_pose  = task.get("pose", Json::nullValue) != Json::nullValue;
	if(return_pose)
	{
		TR.Debug << "Adding silent struct to results" << std::endl;
		SilentFileData silent_data;
		BinarySilentStruct silent_struct(pose, "serialized_pose");
		silent_struct.print_header(data);
		silent_data.write_silent_struct(silent_struct, data);
	}
	else
	{
		TR.Debug << "Adding pdb data to results" << std::endl;
		core::io::pdb::dump_pdb(pose, data);
	}

    TR << data.str();

	pose_object["data"] = data.str();
	pose_object["sequence"] = pose.sequence();
	pose_object["score"] = score;

	Json::Value& pose_list = result["poses"];
	if(pose_list == Json::nullValue)
		pose_list = Json::Value();
	pose_list.append(pose_object);
}

core::kinematics::MoveMapOP parse_movemap(const Json::Value& json)
{
	const Json::Value mm_json = json.get("movemap", Json::nullValue);
	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	if(mm_json != Json::nullValue)
	{
		mm->set_chi(mm_json.get("chi", true).asBool());
		mm->set_bb(mm_json.get("bb", true).asBool());
	}
	else
	{
		mm->set_chi(true);
		mm->set_bb(true);
	}
	return mm;
}

core::pack::task::PackerTaskOP parse_mutations(const core::pose::Pose pose, const Json::Value& json)
{
    core::pack::task::PackerTaskOP packertask =  core::pack::task::TaskFactory::create_packer_task(pose);
    const Json::Value mut_json = json.get("mutations", Json::nullValue);
    TR << "Parsing mutations" << std::endl;
    if(mut_json != Json::nullValue)
    {
       Json::Value::Members keys = mut_json.getMemberNames();
       TR << keys.size() << " mutated positions." << std::endl;

       for(int i=1; i <= pose.n_residue(); i++) {
	  std::string pos = boost::lexical_cast<std::string>(i);
	  if(std::find(keys.begin(), keys.end(), pos) != keys.end()) {
           std::string mutations = mut_json.get(pos, "").asString();
	   TR << i << "\t" << mutations<< std::endl;
           utility::vector1< bool > res_for_i( core::chemical::num_canonical_aas, false );
           foreach(char r, mutations) {
               res_for_i[ core::chemical::aa_from_oneletter_code( r ) ] = true;
           }
               packertask->nonconst_residue_task( i ).restrict_absent_canonical_aas( res_for_i );


	  }else
             packertask->nonconst_residue_task(i).restrict_to_repacking();
	}

    } else {
       TR << "No mutations passed" << std::endl;
    }
    return packertask;
}

core::scoring::ScoreFunctionOP parse_score_function(const Json::Value& json)
{
	const Json::Value scorefxn_json =   json.get("score_function", Json::nullValue);
	std::string weights = scorefxn_json.get("weights", "talaris2013").asString();

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function(weights);

	return scorefxn;
}

Json::Value get_work_queue(std::string const & endpoint)
{

	curlpp::Easy request;
	std::ostringstream os;
	curlpp::options::WriteStream ws(&os);
	request.setOpt(ws);
	curlpp::options::Url url(endpoint);
	request.setOpt(url);
	request.perform();

	Json::Value payload;
	Json::Reader reader;
	bool parsingSuccessful = reader.parse( os.str(), payload, false );
	payload["response"] = boost::lexical_cast<int>(curlpp::infos::ResponseCode::get(request));
	return payload;
}

} // namespace utils
} // namespace micro
} // namespace cyrus
#endif INCLUDED_app_javierbq_micro_utils_hh
