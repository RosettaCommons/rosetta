// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody/LHSnugFitLegacy.hh
/// @brief Build a homology model of an antibody
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#ifndef INCLUDED_protocols_antibody_LHSnugFitLegacy_hh
#define INCLUDED_protocols_antibody_LHSnugFitLegacy_hh






#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>

#include <protocols/antibody_legacy/LHSnugFitLegacy.fwd.hh>



using namespace core;
namespace protocols {
namespace antibody {

class LHSnugFitLegacy: public moves::Mover {


public:

	/// @brief default constructor
	LHSnugFitLegacy();

	/// @brief constructor with arguments
	LHSnugFitLegacy(loops::LoopsOP loops_in );
	LHSnugFitLegacy(antibody::AntibodyInfoOP antibody_in );
	LHSnugFitLegacy(antibody::AntibodyInfoOP antibody_in, bool camelid );

	virtual protocols::moves::MoverOP clone() const;

	/// @brief default destructor
	~LHSnugFitLegacy();

	void set_default();


	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	void set_task_factory(core::pack::task::TaskFactoryCOP tf) {
		tf_ = new core::pack::task::TaskFactory(*tf);
	}

private:

	AntibodyInfoOP ab_info_;

	bool user_defined_;
	bool benchmark_;
	bool is_camelid_;
	loops::LoopsOP all_loops_;
	std::string min_type_;
	Real rot_mag_;
	Real trans_mag_;
	Real temperature_;

	void init(loops::LoopsOP loops_in, bool camelid);

	void setup_objects();

	void snugfit_mcm_protocol( core::pose::Pose & pose_in, loops::Loops loops_in );

	//packer task
	core::pack::task::TaskFactoryOP tf_;

};







} // namespace antibody
} // namespace protocols

#endif








