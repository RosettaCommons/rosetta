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
#include <core/scoring/saxs/SAXSEnergy.hh>
#include <core/scoring/saxs/SAXSEnergyCreatorFA.hh>

#include <core/pose/Pose.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>

using namespace core;

static basic::Tracer trRotateSAXS("RotateSAXS");


void register_options() {
	using namespace core::options;
	using namespace core::options::OptionKeys;

	OPT(in::file::native);
	OPT(in::file::s);
	OPT(in::file::residue_type_set);
	OPT(out::nooutput);
	OPT(score::saxs::q_min);
	OPT(score::saxs::q_max);
	OPT(score::saxs::q_step);
	OPT(score::saxs::custom_ff);
	OPT(score::saxs::ref_spectrum);
}


class RotateSAXS : public protocols::moves::Mover {
public:

	RotateSAXS() {

		step_ = 3.6;
		ires_ = 50;
		using namespace core::options;
		using namespace core::options::OptionKeys;

		core::chemical::ResidueTypeSetCAP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
			core::options::option[ in::file::residue_type_set ]());

		scorefxn_ = core::scoring::get_score_function();

		if ( core::options::option[ in::file::residue_type_set ]() == "fa_standard" ) {
			is_fa_ = true;
			scorefxn_->set_weight( core::scoring::saxs_fa_score, 1.0 );
		} else {
			is_fa_ = false;
			scorefxn_->set_weight( core::scoring::saxs_cen_score, 1.0 );
		}
	}

	virtual ~RotateSAXS() {}

	virtual void apply( core::pose::Pose & pose ) {

		Real phi_is = pose.phi(ires_);
		Real psi_is = pose.psi(ires_);
		for ( Size i=1; i<=100; i++ ) {
			for ( Size j=1; j<=100; j++ ) {
				pose.set_phi(ires_,phi_is + i*step_);
				pose.set_psi(ires_,psi_is + j*step_);
				core::Real score = (*scorefxn_)(pose);
				core::scoring::EnergyMap emap = pose.energies().total_energies();

				std::cerr<< phi_is + i*step_ << " " << psi_is + j*step_ << " "<< emap[ core::scoring::saxs_fa_score ] << std::endl;
			}
		}
	}
	virtual std::string get_name() const { return "RotateSAXS"; }

private:
	Real step_;
	Size ires_;
	core::scoring::ScoreFunctionOP scorefxn_;
	bool is_fa_;
};

int main( int argc, char * argv [] ) {
	try {
		using namespace protocols;
		using namespace protocols::jobdist;
		using namespace protocols::moves;

		register_options();
		devel::init(argc, argv);

		// configure score function


		RotateSAXS debay;

		not_universal_main( debay );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
