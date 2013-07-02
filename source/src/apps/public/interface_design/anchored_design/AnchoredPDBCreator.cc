// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/public/interface_design/anchored_design/AnchoredPDBCreator.cc
/// @brief Anchored Design pre-protocol to make its inputs
/// @author Steven Lewis

// Unit Headers
#include <protocols/analysis/LoopAnalyzerMover.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>

#include <protocols/moves/Mover.hh>

#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/grafting/AnchoredGraftMover.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("apps.public.interface_design.anchored_design.AnchoredPDBCreator");

namespace basic{ namespace options{ namespace OptionKeys{
namespace AnchoredPDBCreator{
basic::options::FileOptionKey const scaffold_loop("AnchoredPDBCreator::scaffold_loop");
basic::options::FileOptionKey const scaffold_pdb("AnchoredPDBCreator::scaffold_pdb");
basic::options::FileOptionKey const anchor_pdb("AnchoredPDBCreator::anchor_pdb");
basic::options::FileOptionKey const target_pdb("AnchoredPDBCreator::target_pdb");
basic::options::IntegerOptionKey const APDBC_cycles("AnchoredPDBCreator::APDBC_cycles");
}//AnchoredPDBCreator
}}}//basic::options::OptionKeys



///@brief APDBC mover
class APDBCMover : public protocols::moves::Mover {

private:

public:
	APDBCMover() : Mover()
	{
		//set up our three poses
		using namespace basic::options::OptionKeys::AnchoredPDBCreator;
		core::import_pose::pose_from_pdb( scaffold, basic::options::option[scaffold_pdb].value() );
		core::import_pose::pose_from_pdb( anchor, basic::options::option[anchor_pdb].value() );
		core::import_pose::pose_from_pdb( target, basic::options::option[target_pdb].value() );

		read_in_insert_loop_info();
	}

	virtual ~APDBCMover(){};

	void read_in_insert_loop_info()
	{
		std::string filename(basic::options::option[ basic::options::OptionKeys::AnchoredPDBCreator::scaffold_loop].value());

		//find the file, open it, handle error
		utility::io::izstream scaffold_loop( filename );
		if ( !scaffold_loop ) {
			Error() << "Can't open scaffold insert loop specification file, looked for: " << filename << std::endl;
			Error() << "use -AnchoredPDBCreator::scaffold_loop <filename> to specify" << std::endl;
			utility_exit();
		}

		core::Size PDBloopstart, PDBinsertstart, PDBloopend;
		char chain;

		scaffold_loop >> chain >> PDBloopstart >> PDBinsertstart >> PDBloopend;

		//if the stream fails, bad formatting
		if( scaffold_loop.fail() ){
			Error() << "Can't parse insert file.  Using the PDB numbering, format is:\n    chain insert_loop_start insert_point insert_loop_end" <<std::endl;
			utility_exit();
		}

		if (chain == '_') chain = ' ';

		core::pose::PDBPoseMap const & pose_map(scaffold.pdb_info()->pdb2pose());
		insert_loop_start = pose_map.find(chain, PDBloopstart);
		insert_loop_end = pose_map.find(chain, PDBloopend);
		insert_point = pose_map.find(chain, PDBinsertstart);

		if( !( insert_loop_start < insert_point)
			&& (insert_point < insert_loop_end) ){
			utility_exit_with_message("insert/loop definitions must obey insert_loop_start < insert_point < insert_loop_end");
		}

	}//read_in_insert_loop_info

	virtual
	void
	apply( core::pose::Pose & pose ){

		//domain insertion
		//using protocols::toolbox::pose_manipulation::insert_pose_into_pose;
		core::Size const cycles(basic::options::option[ basic::options::OptionKeys::AnchoredPDBCreator::APDBC_cycles ].value());


		//core::pose::Pose combined(insert_pose_into_pose(
		//	scaffold,
		//	anchor,
		//	insert_loop_start,
		//	insert_point,
		//	loop_end,
		//	cycles));

		protocols::grafting::AnchoredGraftMover grafter = protocols::grafting::AnchoredGraftMover(insert_point, insert_point+1);
		Size const Nter_flexibility = insert_point-insert_loop_start+1;
		Size const Cter_flexibility = insert_loop_end-insert_point;

		//(jadolfbr) Set AnchoredGraftMover options.
		grafter.set_scaffold_flexibility(Nter_flexibility, Cter_flexibility);
        grafter.set_cycles(cycles);
        core::pose::Pose piece(anchor);
		grafter.set_piece(piece, 0, 0);
		grafter.set_mintype("dfpmin_armijo");//mintype from pose_into_pose
		Pose combined(scaffold);//Copy the scaffold into combined before passing to Mover.
		grafter.apply(combined);


		//check on the loop quality
		protocols::loops::Loops loops;
		loops.add_loop(insert_loop_start, grafter.get_Cter_loop_end(), insert_point);
		protocols::analysis::LoopAnalyzerMover LAM(loops);
		LAM.apply(combined);

		using namespace protocols::jd2;
		JobCOP job_me( JobDistributor::get_instance()->current_job() );
		//JobDistributor::get_instance()->job_outputter()->other_pose( job_me, combined, "closed");

		//superimpose combined onto anchor
		utility::vector1< core::Size > positions;
		core::Size const anchorlength = anchor.total_residue();
		for(core::Size i(1); i<=anchorlength; ++i) positions.push_back(i);
		using protocols::toolbox::pose_manipulation::superimpose_pose_on_subset_CA;
		superimpose_pose_on_subset_CA(combined, anchor, positions, insert_point);

		//JobDistributor::get_instance()->job_outputter()->other_pose( job_me, combined, "moved");

		core::pose::Pose total(target);
		//add combined to target
		total.append_residue_by_jump(combined.residue(1), total.total_residue(), "CA", "CA", true);
		for (core::Size i = 2; i <= combined.total_residue(); ++i){
			//TR << i << std::endl;
			total.append_residue_by_bond(combined.residue(i));
		}

		pose = total;
		return;
	}

	virtual std::string get_name() const { return "APDBCMover"; }


private:
	core::Size insert_loop_start;
	core::Size insert_point;
	core::Size insert_loop_end;

	core::pose::Pose scaffold;
	core::pose::Pose anchor;
	core::pose::Pose target;
};

typedef utility::pointer::owning_ptr< APDBCMover > APDBCMoverOP;

int main( int argc, char* argv[] )
{
	try {
	using basic::options::option;
	using namespace basic::options::OptionKeys::AnchoredPDBCreator;
	option.add( scaffold_loop, "scaffold anchor loop location file").def("scaffold_loop");
	option.add( scaffold_pdb, "scaffold pdb location file").def("scaffold.pdb");
	option.add( anchor_pdb, "anchor pdb location file").def("anchor.pdb");
	option.add( target_pdb, "target pdb location file").def("target.pdb");
	option.add( APDBC_cycles, "loop closure cycles").def(500);

	devel::init(argc, argv);

	protocols::jd2::JobDistributor::get_instance()->go(new APDBCMover);

	TR << "************************d**o**n**e**************************************" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

 	return 0;
}
