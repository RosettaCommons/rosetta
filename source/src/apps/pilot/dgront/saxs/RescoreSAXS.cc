#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/moves/Mover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/Pose.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>


using namespace core;

static THREAD_LOCAL basic::Tracer trRescoreSAXS( "RescoreSAXS" );


void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

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


class RescoreSAXS : public protocols::moves::Mover {
public:

	RescoreSAXS() {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		core::chemical::ResidueTypeSetCAP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
			basic::options::option[ in::file::residue_type_set ]());

		scorefxn_ = core::scoring::get_score_function();

		if ( basic::options::option[ in::file::residue_type_set ]() == "fa_standard" ) {
			is_fa_ = true;
			scorefxn_->set_weight( core::scoring::saxs_fa_score, 1.0 );
		} else {
			is_fa_ = false;
			scorefxn_->set_weight( core::scoring::saxs_cen_score, 1.0 );
		}
	}

	virtual ~RescoreSAXS() {}

	virtual void apply( core::pose::Pose & pose ) {

		core::Real score = (*scorefxn_)(pose);
		core::scoring::EnergyMap emap = pose.energies().total_energies();
		if ( is_fa_ ) {
			trRescoreSAXS.Warning << "SAXS  score: " << emap[ core::scoring::saxs_fa_score ]<< std::endl;
		} else {
			trRescoreSAXS.Warning << "SAXS  score: " << emap[ core::scoring::saxs_cen_score ]<<std::endl;
		}
		trRescoreSAXS.Warning << "Total score: " << score<<std::endl;
	}
	virtual std::string get_name() const { return "RescoreSAXS"; }

private:
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


		RescoreSAXS debay;

		not_universal_main( debay );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
