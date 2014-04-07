#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/not_universal_main.hh>
//#include <protocols/moves/Mover.hh>
#include <protocols/abinitio/AbrelaxMover.hh>

#include <core/options/option.hh>
#include <core/options/keys/OptionKeys.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/score.OptionKeys.gen.hh>
#include <core/options/option_macros.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/saxs/PDDFEnergy.hh>

#include <core/pose/Pose.hh>


#include <devel/init.hh>
#include <core/types.hh>
#include <core/util/Tracer.hh>

#include <utility/excn/Exceptions.hh>

using namespace core;

static core::util::Tracer trFilteredAbrelax("FilteredAbrelax");


int main( int argc, char * argv [] ) {
    try {
    using namespace protocols;
    using namespace protocols::jobdist;
    using namespace protocols::moves;
    using namespace core::options;
    using namespace core::options::OptionKeys;

//    register_options();
    devel::init(argc, argv);

    try {
        protocols::abinitio::AbrelaxMover runMe;
	not_universal_main( runMe );
    } catch ( utility::excn::EXCN_Base &excn ) {
    	std::cerr<< excn.msg() <<std::endl;
    }
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }

    return 0;
}
