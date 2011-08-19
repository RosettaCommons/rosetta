#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/init.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/util.hh>//option.hh>
#include <protocols/loops/LoopRelaxMover.fwd.hh>
#include <protocols/loops/LoopRelaxMover.hh>
#include <core/kinematics/FoldTree.hh>
#include <numeric/random/random.hh>
#include <protocols/loops/Loops.hh>
#include <core/chemical/VariantType.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/conformation/util.hh>
#include <core/pose/disulfide_util.hh>
#include <protocols/loops/LoopRelaxMover.fwd.hh>
#include <devel/init.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/conformation/Conformation.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
//#include <protocols/jobdist/standard_mains.hh>
//#include <protocols/moves/Mover.hh>

#include <ObjexxFCL/format.hh>

// C++ headers
//#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>


using basic::T;
using basic::Warning;
using basic::Error;
namespace ss_loop_builder {
	basic::options::IntegerOptionKey start_cys("ss_loop_builder:start");
	basic::options::IntegerOptionKey   end_cys("ss_loop_builder:end");
	basic::options::IntegerOptionKey   ncycles("ss_loop_builder:ncycles");
}
std::string int2string(int n){
	      std::stringstream Num;
	      std::string str;
	      Num << n;
	      str = Num.str();
	      return str;
}
void disulfide(core::pose::Pose& pose) {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	protocols::loops::LoopRelaxMover loop_relax_mover;
	loop_relax_mover.remodel("perturb_kic");
	loop_relax_mover.refine("refine_kic");
	loop_relax_mover.copy_sidechains(false);
	Size pepsize = pose.n_residue();

}

int main( int argc, char** argv ){

	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	option.add( ss_loop_builder::start_cys,
								"lower seq. cys to form the loop in pose numbering").def(-1);
		option.add( ss_loop_builder::end_cys,"higher seq. cys to form the loop in pose numbering").def(-1);
		option.add( ss_loop_builder::ncycles,"Num. cycles of relax/loopmodel and S-S samples").def(1);
	core::init(argc,argv);
	basic::Tracer TR("protocols.moves.CycPep");
	pose::Pose workpose;
	core::import_pose::pose_from_pdb(workpose, basic::options::start_file());
	option[ out::file::fullatom ].value(true);
	core::scoring::ScoreFunctionOP scorefxn(
				core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::STANDARD_WTS, core::scoring::SCORE12_PATCH ));
	Size start_cys = option[ ss_loop_builder::start_cys ]();
	Size end_cys   = option[ ss_loop_builder::end_cys ]();
	int const nstruct = option[ OptionKeys::out::nstruct ];
	Size cycles = option[ ss_loop_builder::ncycles ]();

	protocols::loops::Loops loops;
	loops.add_loop(start_cys+1,end_cys-1,(start_cys+end_cys)/2);

	core::conformation::form_disulfide(workpose.conformation(), start_cys, end_cys);
	protocols::loops::LoopRelaxMover loop_relax_mover;
	//loop_relax_mover.scorefxns(scorefxn,scorefxn);
	loop_relax_mover.remodel("perturb_kic");
	loop_relax_mover.refine("refine_kic");
	loop_relax_mover.loops(loops);
	loop_relax_mover.copy_sidechains(false);

	for (Size i = 0; i<cycles; ++i) {
		loop_relax_mover.apply(workpose);
	}
	workpose.dump_pdb("disulfide_loop.pdb");
}

//int main(int argc, char** argv) {
//	using namespace core;
//	using namespace basic::options;
//	using namespace basic::options::OptionKeys;
//
//	core::init(argc,argv);
//	core::scoring::ScoreFunctionOP scorefxn_ = core::scoring::getScoreFunction();
//	protocols::loops::Loops myloop;
//	myloop.add_loop(2,6,4);
//	protocols::loops::LoopRelaxMover loop_relax_mover;
//	loop_relax_mover.remodel("perturb_kic");
//	loop_relax_mover.refine("refine_kic");
//	loop_relax_mover.loops(myloop);
//	loop_relax_mover.copy_sidechains(false);
//	basic::Tracer TR("protocols.moves.CycPep");
//	pose::Pose mypose;
//	core::import_pose::pose_from_pdb(mypose, options::start_file());
//	option[ out::file::fullatom ].value(true);
//	disulfide(mypose);
//	return 0;
//	Size pepsize = mypose.n_residue();
//	pose::Pose temp = mypose;
////	temp.dump_pdb("temp_0.pdb");
//	core::pose::remove_lower_terminus_type_from_pose_residue(temp,1);
//	core::pose::remove_upper_terminus_type_from_pose_residue(temp,pepsize);
////	temp.dump_pdb("temp_1.pdb");
////
////	for (Size i=1; i<=pepsize-1; ++i){
////		temp.append_polymer_residue_after_seqpos(mypose.residue(i),i,false);
////		temp.dump_pdb("temp_"+int2string(i)+".pdb");
////	}
////
//	TR << "INITIAL Score: " << scorefxn_->score(temp)<<std::endl;
//	return 0;
//	temp.dump_pdb("start.pdb");
//	for (Size i=1; i<=pepsize; ++i){
////		core::pose::remove_lower_terminus_type_from_pose_residue(temp,pepsize);
////		core::pose::remove_upper_terminus_type_from_pose_residue(temp,1);
//		pack::task::PackerTaskOP task;
//		task = pack::task::TaskFactory::create_packer_task( temp );
//		pack::task::operation::NoRepackDisulfides noRepackDisulf;
//		noRepackDisulf.apply(temp, *task);
//		loop_relax_mover.apply(temp);
//		//temp.dump_pdb("test"+int2string(i)+".pdb");
//		//mypose.clear();
//		//mypose.append_residue_by_bond(temp.residue(pepsize));
////		temp.dump_pdb("test_temp_before.pdb");
//
////		if (i==1)
////			temp.dump_pdb("test_temp_after.pdb");
//		pose::Pose tmp2(temp, pepsize, pepsize);
//		core::pose::remove_lower_terminus_type_from_pose_residue(tmp2,1);
//		core::pose::remove_upper_terminus_type_from_pose_residue(tmp2,1);
//		for( Size j=1; j<=pepsize-1; ++j){
//			tmp2.append_polymer_residue_after_seqpos(temp.residue(j),j,false);
//		}
//		//TR << tmp2.sequence() << std::endl;
//
//		core::pose::remove_lower_terminus_type_from_pose_residue(tmp2,pepsize);
//		core::pose::remove_upper_terminus_type_from_pose_residue(tmp2,1);
//		tmp2.dump_pdb("tst2"+int2string(i)+".pdb");
////		for (Size i=1; i<=pepsize; ++i){
////			TR << "res: " << i << " lower: " << tmp2.residue(i).has_variant_type(core::chemical::LOWER_TERMINUS) <<std::endl;
////			TR << "res: " << i << " upper: " << tmp2.residue(i).has_variant_type(core::chemical::UPPER_TERMINUS) <<std::endl;
////		}
//		temp = tmp2;
//		break;
//	}
//
////	for (Size i=1; i<=pepsize; ++i){
////		//loop_relax_mover.apply(mypose);
////
////		mypose.replace_residue(1,temp.residue(pepsize),true);
////		std::cerr<<"TTTTTTT\n";
////		for( Size j=1; j<=pepsize-1; ++j){
////			mypose.replace_residue(j+1,temp.residue(j),true);
////			mypose.conformation().update_polymeric_connection(j);
////			mypose.conformation().update_polymeric_connection(j+1);
////		}
////		temp = mypose;
////	}
//
//	/**
//	 * This is the actual definition of a new edge + cut required!!!
//	 */
////	fd.delete_edge(e);
////	fd.add_edge(14,1,-1);
////	fd.new_jump(5,7,6);
////	mypose.fold_tree(fd);
////	for (Size i =1; i<2; ++i){
////		std::cout<<i<<" "<<mypose.total_residue()-i<<" "<<(mypose.total_residue()-i)/2<<std::endl;
////		myloop.add_loop(3,13,8);
////
////		loop_relax_mover.loops(myloop);
////		loop_relax_mover.apply(mypose);
////
////	}
//	io::pdb::dump_pdb(mypose,"test.pdb");
//
//}
