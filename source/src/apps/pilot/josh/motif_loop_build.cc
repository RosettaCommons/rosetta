#include <iostream>
#include <devel/init.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>

#include <devel/dna/MotifLoopBuild.hh>


int main( int argc, char * argv[] )
{
	try {
		devel::init( argc, argv);
		devel::dna::MotifLoopBuild this_run =  devel::dna::MotifLoopBuild();

		core::pose::Pose pose;
		if ( ! basic::options::option[ basic::options::OptionKeys::in::file::s].user() ) {
			utility_exit_with_message("No input pdb provided. Use \"-s\" name.pdb flag.");
		}
		std::string pdb_name( basic::options::start_file() );
		core::import_pose::pose_from_file(pose, pdb_name, core::import_pose::PDB_file);

		this_run.apply(pose);
		//pose.dump_pdb("final_model.pdb");
		std::cout << "******Sucessfull Run!******\n\n";
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
