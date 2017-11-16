#include <protocols/jobdist/standard_mains.hh>
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
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/scoring/methods/saxs/PDDFEnergy.hh>

#include <core/pose/Pose.hh>
#include <core/id/NamedAtomID.hh>
#include <numeric/model_quality/rms.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

#include <utility/pointer/owning_ptr.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <fstream>
#include <iostream>

#include <utility/excn/Exceptions.hh>

using namespace core;
using namespace ObjexxFCL;

static basic::Tracer trAllCrmsd( "AllCrmsd" );


OPT_1GRP_KEY( Boolean, saxs, show_pddf )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	OPT(in::file::native);
	OPT(in::file::s);
	OPT(in::file::residue_type_set);
	OPT(out::nooutput);
}

class XYZ : public FArray2D_double, public utility::pointer::ReferenceCount {
public:
	XYZ(core::Size j,core::Size i,core::Real d) : FArray2D_double(j,i,d){}
};

typedef utility::pointer::owning_ptr<XYZ> XYZOP;

class CollectCoordinates : public protocols::moves::Mover {
public:

	CollectCoordinates(utility::vector1<XYZOP> & xyz) : xyz_(xyz) {}

	virtual ~CollectCoordinates() {}

	virtual void apply( core::pose::Pose & pose ) {

		Size len = pose.size();
		XYZOP c = new XYZ(3,len,0.0);
		copy_coordinates(pose,c);
		xyz_.push_back(c);
	}

private:
	utility::vector1<XYZOP> & xyz_;
	void copy_coordinates(pose::Pose& in_pose,XYZOP dst) {

		Size len = in_pose.size();
		for ( core::Size i = 1; i <= len; i++ ) {
			id::NamedAtomID idCA("CA", i);
			PointPosition const& xyz = in_pose.xyz(idCA);
			for ( core::Size d = 1; d <= 3; ++d ) {
				(*dst)(d, i) = xyz[d - 1];
			}
		}
	}
};


class AllCrmsd {
public:

	AllCrmsd(utility::vector1<XYZOP> & xyz) : xyz_(xyz) {}

	Size calculate() {

		std::ofstream crmsd_out;
		crmsd_out.open ("all_crmsd");
		Size cnt = 0;
		for ( Size i=2; i<=xyz_.size(); i++ ) {
			for ( Size j=1; j<i; j++ ) {
				Real val = numeric::model_quality::rms_wrapper(xyz_[i]->size2(),*xyz_[i],*xyz_[j]);
				crmsd_out<<i<<" "<<j<<" "<<val<<std::endl;
				cnt++;
			}
		}
		crmsd_out.close();
		return cnt;
	}

private:
	utility::vector1<XYZOP> & xyz_;
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

		utility::vector1<XYZOP> xyz;

		time_t time_start = time(NULL);
		CollectCoordinates cc(xyz);
		not_universal_main( cc );
		time_t time_end = time(NULL);
		trAllCrmsd << xyz.size() << " poses cached in "<<(time_end - time_start)<<" seconds"<<std::endl;

		AllCrmsd all(xyz);
		time_start = time(NULL);
		Size cnt = all.calculate();
		time_end = time(NULL);
		trAllCrmsd << "Computed " << cnt <<" crmsd values in "<<(time_end - time_start)<<" seconds"<<std::endl;
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
