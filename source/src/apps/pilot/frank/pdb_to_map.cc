/// @file
/// @brief


#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>


int
main( int argc, char * argv [] ) {
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		devel::init(argc, argv);

		utility::vector1<core::pose::PoseOP> models = core::import_pose::poseOPs_from_pdbs(option[OptionKeys::in::file::s]());
		core::Real mapreso = option[OptionKeys::edensity::mapreso]();
		core::Real gridspacing = option[OptionKeys::edensity::grid_spacing]();
		core::scoring::electron_density::ElectronDensity edens( models, mapreso, gridspacing );
		std::string outfile = option[OptionKeys::edensity::mapfile]();
		edens.writeMRC(outfile);
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
