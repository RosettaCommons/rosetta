// AUTO-REMOVED #include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/moves/Mover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/EnergyMap.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/saxs/SAXSEnergyCEN.hh>
#include <core/scoring/saxs/SAXSEnergyFA.hh>
#include <core/scoring/saxs/SAXSEnergyCreator.hh>
#include <core/scoring/saxs/SAXSEnergyCreatorFA.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

#include <protocols/jobdist/Jobs.hh>
#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>


using namespace core;

static basic::Tracer trRescoreSAXS("ComputeSAXSSpectrum");


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
}


class ComputeSAXSSpectrum : public protocols::moves::Mover {
public:

    ComputeSAXSSpectrum() {

    using namespace basic::options;
    using namespace basic::options::OptionKeys;

	if(basic::options::option[ in::file::residue_type_set ]() == "fa_standard") {
	    saxs_ = new core::scoring::saxs::SAXSEnergyFA();
	} else {
	    saxs_ = new core::scoring::saxs::SAXSEnergyCEN();
	}
	mdl_cnt_ = 1;
    }

    virtual ~ComputeSAXSSpectrum() {}

	virtual void apply( core::pose::Pose & pose ) {

	    saxs_->total_energy(pose);
	    utility::vector1<Real> q = saxs_->get_q();
	    utility::vector1<Real> I = saxs_->get_pose_intensities();

	    std::cout << mdl_cnt_<<" "<<0.0<<" "<<std::setw(15)<<std::setprecision(8)<<
		    saxs_->compute_zero_intensity()<<std::endl;
	    for(Size i=1;i<=q.size();i++) {
		std::cout << mdl_cnt_<<" "<<q[i]<<" "<<std::setw(15)<<std::setprecision(8)<<I[i]<<std::endl;
	    }
	    mdl_cnt_++;
	}

	virtual std::string get_name() const { return "ComputeSAXSSpectrum"; }

private:
    Size mdl_cnt_;
    core::scoring::saxs::SAXSEnergy* saxs_;
};

int main( int argc, char * argv [] ) {
    try {
    using namespace protocols;
    using namespace protocols::jobdist;
    using namespace protocols::moves;

    register_options();
    devel::init(argc, argv);

    ComputeSAXSSpectrum debay;

    not_universal_main( debay );

    } catch ( utility::excn::EXCN_Base const & e ) {
                             std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                }
       return 0;
}
