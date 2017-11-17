//////////////////////////////
// EXE
//////////////////////////////
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/mc.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <core/io/silent/SilentFileData.hh>
#include <iostream>

#include <protocols/moves/mc_convergence_checks/HPool.hh>

int main(int argc, char *argv[])
{

	try {

	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::moves::mc_convergence_checks;

	devel::init(argc, argv);
	HPool_RMSD_OP hp = new HPool_RMSD(option[mc::known_structures]);

	//hp->debug();

    io::silent::SilentFileData sfd;
    sfd.read_file( option[ in::file::silent ]()[1] );

    std::string tag;
    core::Real best_rmsd;

    for ( io::silent::SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it )
    {
        //hp->evaluate(**it, option[ cluster::K_radius ]()[1], tag, best_rmsd);
        hp->evaluate(**it, option[ cluster::K_threshold ](), tag, best_rmsd);
        std::cout << "decoy:" << it->decoy_tag() << " " << tag << " " << best_rmsd << std::endl;
    }

    std::cout << "number of centers: " << hp->size() << std::endl;

	return 0;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

