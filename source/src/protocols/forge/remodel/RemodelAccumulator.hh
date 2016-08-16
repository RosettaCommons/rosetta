// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/forge/remodel/RemodelAccumulator.hh
/// @brief
/// @author Possu Huang ( possu@uw.edu )


#ifndef INCLUDED_protocols_forge_remodel_RemodelAccumulator_hh
#define INCLUDED_protocols_forge_remodel_RemodelAccumulator_hh

#include <protocols/cluster/cluster.hh>
//#include <protocols/forge/remodel/RemodelData.hh>
#include <protocols/forge/remodel/RemodelWorkingSet.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace remodel {

class RemodelAccumulator: public protocols::moves::Mover {

private: // typedefs

	typedef protocols::moves::Mover Super;

public: // typedefs

	typedef core::Real Real;
	typedef core::Size Size;

	typedef core::kinematics::MoveMap MoveMap;
	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionFactory ScoreFunctionFactory;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::pack::task::PackerTask PackerTask;
	typedef protocols::moves::MoverOP MoverOP;
	typedef protocols::forge::remodel::RemodelData RemodelData;
	typedef utility::pointer::shared_ptr< protocols::cluster::ClusterPhilStyle > ClusterPhilStyleOP;


public: //constructor/destructor

	RemodelAccumulator();

	//RemodelAccumulator(RemodelData remodeldata);
	RemodelAccumulator(RemodelWorkingSet & working_model);

	virtual
	~RemodelAccumulator();

public: // virtual constructors

	virtual
	MoverOP clone() const;

	virtual
	MoverOP fresh_instance() const;

public: // options

public:

	void keep_top_pose(core::Size num_to_keep);
	void cluster_pose();
	void cluster_loop();

	using Super::apply;

	virtual void apply( Pose & pose );
	virtual void apply( Pose & pose, core::Real & score );
	virtual Pose pop();
	virtual void clear();
	virtual core::Size size();
	virtual std::string get_name() const;


private: // data

	bool cluster_switch_;  // check cluster radius for on/off state.  0 radius is no clustering

	//RemodelData remodel_data_; // design mode determined in here
	RemodelWorkingSet working_model_;
	//PackerTask packer_task_;
	ScoreFunctionOP sfxn_;
	//MoveMap move_map_;
	ClusterPhilStyleOP cluster_;
	std::multimap<core::Real,core::pose::PoseOP> pose_store_;

public: // accessors
	core::Size recover_checkpoint();

	void write_checkpoint(core::Size progress_point);

	bool cluster_switch();

	//ScoreFunction const & scorefunction() const;

	//MoveMap const & movemap() const;

	//PackerTask const & packertask() const;

	void run_cluster();
	void shrink_cluster(core::Size num_top);
	std::vector<core::pose::PoseOP>  clustered_best_poses();
	std::vector<core::pose::PoseOP>  contents_in_pose_store();
	std::vector<core::pose::PoseOP>  clustered_top_poses(core::Size count);
};

} // namespace remodel
} // namespace forge
} // namespace protocols

#endif /* INCLUDED_protocols_forge_remodel_RemodelAccumulator_HH */
