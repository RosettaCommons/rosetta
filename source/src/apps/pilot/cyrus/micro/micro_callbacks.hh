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
#include <protocols/moves/DsspMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>



#include <apps/pilot/cyrus/micro/micro_utils.hh>
#include <jsoncpp/json/json.h>


// CALLBACKS

namespace cyrus {
namespace micro {
namespace callbacks {

static THREAD_LOCAL basic::Tracer TR( "cyrus_micro_callbacks" );

using utils::get_pose_from_json;
using utils::add_pose_to_result_value;

Json::Value pdb_to_pose(const Json::Value&  task) {
    using namespace core::io::silent;
    std::string pdb_data = task.get("pdb_data","").asString();
    std::string pdb_id = task.get("pose_id","").asString();
    core::pose::Pose pose;
    TR << "Loading PDB" << std::endl;
    core::import_pose::pose_from_pdbstring(pose, pdb_data);
    TR << "n_residues: " << pose.n_residue() <<std::endl;

    Json::Value results = Json::Value(Json::objectValue);
    add_pose_to_result_value(task, results, pose);
    return results;
}

Json::Value pose_to_pdb(const Json::Value&  task) {
    using namespace core::io::silent;

    core::pose::Pose pose = get_pose_from_json(task);
    protocols::moves::DsspMover dssp;
    dssp.apply(pose);


    std::ostringstream data_out;
    core::io::pdb::dump_pdb(pose, data_out);

    Json::Value results = Json::Value(Json::objectValue);
    results["pdb_data"] = data_out.str();
    results["seq"] = pose.sequence();
    results["ss"] = pose.secstruct();

    return results;
}

Json::Value fast_relax(const Json::Value&  task) {
    using namespace core::io::silent;
    core::pose::Pose pose = get_pose_from_json(task);
    core::pose::PoseOP ref = pose.clone();

    core::scoring::ScoreFunctionOP score_fxn = cyrus::micro::utils::parse_score_function(task);
    protocols::relax::FastRelax relax(score_fxn);
    core::kinematics::MoveMapOP mm = cyrus::micro::utils::parse_movemap(task);
    relax.set_movemap(mm);
    relax.apply(pose);

    Json::Value results = Json::Value(Json::objectValue);
    add_pose_to_result_value(task, results, pose, ref);
    return results;
}


Json::Value repack(const Json::Value&  task) {
    using namespace core::io::silent;
    core::pose::Pose pose = get_pose_from_json(task);

    core::scoring::ScoreFunctionOP score_fxn = cyrus::micro::utils::parse_score_function(task);
    core::pack::task::PackerTaskOP packertask =  core::pack::task::TaskFactory::create_packer_task(pose);
    packertask->restrict_to_repacking();

    protocols::simple_moves::PackRotamersMover packer(score_fxn, packertask);
    packer.apply(pose);

    Json::Value results = Json::Value(Json::objectValue);
    add_pose_to_result_value(task, results, pose);
    return results;
}


Json::Value design(const Json::Value&  task) {
    using namespace core::io::silent;
    core::pose::Pose pose = get_pose_from_json(task);

    core::scoring::ScoreFunctionOP score_fxn = cyrus::micro::utils::parse_score_function(task);
    core::pack::task::PackerTaskOP packertask =  utils::parse_mutations(pose, task);

    protocols::simple_moves::PackRotamersMover packer(score_fxn, packertask);
    packer.apply(pose);

    Json::Value results = Json::Value(Json::objectValue);
    add_pose_to_result_value(task, results, pose);
    return results;
}

Json::Value minimize(const Json::Value&  task) {
    using namespace core::io::silent;
    core::pose::Pose pose = get_pose_from_json(task);
    core::pose::PoseOP ref = pose.clone();
    core::scoring::ScoreFunctionOP scorefxn = cyrus::micro::utils::parse_score_function(task);
    core::kinematics::MoveMapOP mm = cyrus::micro::utils::parse_movemap(task);

    protocols::simple_moves::MinMover min( mm, scorefxn, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
    min.apply(pose);

    Json::Value results = Json::Value(Json::objectValue);
    add_pose_to_result_value(task, results, pose, ref);
    return results;
}

} // namespace callbacks
} // namespace micro
} // namespace cyrus
