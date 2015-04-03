#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/moves/Mover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>


#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/scoring/methods/saxs/PDDFEnergy.hh>

#include <core/pose/Pose.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>

using namespace core;

static thread_local basic::Tracer trRescorePDDF( "RescorePDDF" );


OPT_1GRP_KEY( Boolean, saxs, show_pddf )

void register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  OPT(in::file::native);
  OPT(in::file::s);
  OPT(in::file::residue_type_set);
  OPT(out::nooutput);
  OPT(score::saxs::ref_pddf);
  OPT(score::saxs::custom_ff);

  NEW_OPT( saxs::show_pddf  , "prints PDDF on stdout", "" );
}


class RescorePDDF : public protocols::moves::Mover {
public:

	RescorePDDF() {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		scorefxn_ = core::scoring::get_score_function();
		scorefxn_->set_weight( core::scoring::pddf_score, 1.0 );
	}

	virtual ~RescorePDDF() {}

	virtual void apply( core::pose::Pose & pose ) {
		core::Real score = (*scorefxn_)(pose);
		trRescorePDDF << "PDDF score: "<< score<<std::endl;
		//core::scoring::EnergyMap emap = pose.energies().total_energies();  // unused ~Labonte
	}

	virtual std::string get_name() const { return "RescorePDDF"; }

private:
    core::scoring::ScoreFunctionOP scorefxn_;
};


class PrintPDDF : public protocols::moves::Mover {
public:

	PrintPDDF() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		pddf_score_ = new protocols::scoring::methods::saxs::PDDFEnergy();
	}

	virtual ~PrintPDDF() {}

	virtual void apply( core::pose::Pose & pose ) {

		pddf_score_->compute_pddf_without_ff(pose);
		utility::vector1<Real> & d = pddf_score_->get_dist_bins();
		utility::vector1<Real> & pddf = pddf_score_->get_pddf();
		for(Size i=1;i<=d.size();i++) {
			std::cout << d[i]<<" "<<pddf[i]<<std::endl;
		}
	}
 	virtual std::string get_name() const { return "PrintPDDF"; }


private:
    protocols::scoring::methods::saxs::PDDFEnergy* pddf_score_;
};

int main( int argc, char * argv [] ) {
    try {
    using namespace protocols;
    using namespace protocols::jobdist;
    using namespace protocols::moves;
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    register_options();
    devel::init(argc, argv);

    // configure score function


    if(option[saxs::show_pddf].user()) {
      PrintPDDF printPddf;
      not_universal_main( printPddf );
    } else {
      RescorePDDF rescore;
      not_universal_main( rescore );
    }
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
    return 0;
}
