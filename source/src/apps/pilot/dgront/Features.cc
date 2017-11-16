#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/moves/Mover.hh>

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
#include <basic/Tracer.hh>

using namespace core;

static basic::Tracer trRescorePDDF("RescorePDDF");


OPT_1GRP_KEY( Boolean, saxs, show_pddf )

void register_options() {
	using namespace core::options;
	using namespace core::options::OptionKeys;

	OPT(in::file::native);
	OPT(in::file::s);
	OPT(in::file::residue_type_set);
	OPT(out::nooutput);
	OPT(in::file::silent);
}


class PrintFeatures : public protocols::moves::Mover {
public:

	PrintFeatures() {
	}

	virtual void apply( core::pose::Pose & pose ) {

		for ( Size i=1; i<=pose.size(); i++ ) {
			char tmp = torsion2big_bin(pose.phi(i),
				pose.psi(i),pose.omega(i));
			std::cout<<tmp;
		}
		std::cout<<'\n';
	}

private:
	Size n_;

	char torsion2big_bin(core::Real const phi,  core::Real const psi,  core::Real const omega) {

		if ( std::abs( omega ) < 90 ) {
			return 'O'; // cis-omega
		} else if ( phi >= 0.0 ) {
			if ( -100 < psi && psi <= 100 ) {
				return 'G'; // alpha-L
			} else {
				return 'E'; // E
			}
		} else {
			if ( -125 < psi && psi <= 50 ) {
				return 'A'; // helical
			} else {
				return 'B'; // extended
			}
		}
		return 0;
	}

};


int main( int argc, char * argv [] ) {
	try {
		using namespace protocols;
		using namespace protocols::jobdist;
		using namespace protocols::moves;
		using namespace core::options;
		using namespace core::options::OptionKeys;

		register_options();
		devel::init(argc, argv);

		PrintFeatures f;
		not_universal_main( f );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
