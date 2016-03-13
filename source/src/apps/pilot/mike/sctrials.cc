// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


#include <protocols/jd2/JobDistributor.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/Energies.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DiagnosticData.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/jd2/ScoreMap.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>
#include <string>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <ObjexxFCL/string.functions.hh>

using namespace ObjexxFCL;

using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;
using namespace protocols;
using namespace protocols::moves;
using namespace core::scoring;

using utility::vector1;

static THREAD_LOCAL basic::Tracer TR( "sctrials" );

class ScTrials : public moves::Mover {

public:
	ScTrials();

	~ScTrials();

	virtual MoverOP clone() const;
	virtual MoverOP fresh_instance() const;

	virtual void apply( Pose & pose );
	virtual std::string get_name() const;
	virtual void test_move( Pose & pose )
	{
		apply(pose);
	}

private:
	ScoreFunctionOP score_function_;
	std::map< std::string, core::Real > score_map_;
	bool verbose_;
	std::string scorefile_;
};

ScTrials::ScTrials() :
	Mover( "rama_test_mover" ),
	score_function_( get_score_function() )
{}

ScTrials::~ScTrials() {}

MoverOP ScTrials::clone() const {
	return new ScTrials( *this );
}
MoverOP ScTrials::fresh_instance() const {
	return new ScTrials;
}

void
ScTrials::apply( Pose & pose ) {
	using namespace pose;
	using datacache::CacheableDataType::SCORE_MAP;

	core::pose::Pose start_pose = pose;
	core::Real start_score = (*score_function_)(pose);
	std::cout << "StartScore: " << start_score << std::endl;
	const core::Real limit=4.0;
	core::Size total_rotamers = 0;
	core::Real total_free_energy = 0;
	core::Real total_enthalpy    = 0;
	core::Real total_minus_T_times_entropy = 0;
//	RotamerLibraryScratchSpaceOP = new core::pack::dunbrack::RotamerLibraryScratchSpace();

	TR << protocols::jd2::JobDistributor::get_instance()->current_output_name();
	for( core::Size ir=1; ir <= pose.total_residue(); ++ ir ){
		EnergyMap emap = pose.energies().onebody_energies(ir);

		core::pack::dunbrack::RotamerLibrary const & rotamer_library_ = * core::pack::dunbrack::RotamerLibrary::get_instance();
		core::pack::rotamers::SingleResidueRotamerLibraryCAP residue_rotamer_library(
    	rotamer_library_.get_rsd_library(pose.residue( ir ).type())
		);
		core::pack::dunbrack::SingleResidueDunbrackLibraryCAP dun_rotlib(
		          dynamic_cast< core::pack::dunbrack::SingleResidueDunbrackLibrary const * >
							          ( residue_rotamer_library.get() ));

		if ( dun_rotlib == 0 ) continue;

		utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > sample_data = dun_rotlib->get_all_rotamer_samples( pose.phi( ir ), pose.psi(ir ) );


		TR << "NRotamers: " << sample_data.size() << std::endl;

		total_rotamers += sample_data.size();


		std::vector< core::Real > energies;

		for( core::Size j = 1; j <= sample_data.size(); j++ ){
			pose = start_pose;

			core::Size  nchis = sample_data[j].nchi();
			core::pack::dunbrack::Real4 chis = sample_data[j].chi_mean();

			// make a random chi angle change
			for( core::Size ichi=1; ichi <= pose.residue( ir ).nchi(); ichi ++ ){
				pose.set_chi( ichi, ir, chis[ichi] );
			}

			core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.001, true, false );
			core::kinematics::MoveMap final_mm;
			final_mm.set_chi( true );
			final_mm.set_bb( false );
			core::optimization::AtomTreeMinimizer().run( pose, final_mm, *score_function_, options );
			core::Real final_score = (*score_function_)(pose);
			energies.push_back( final_score );
		}

		core::Real lowest_energy = *( energies.begin() );
		for( std::vector< core::Real >::const_iterator it = energies.begin();
					it != energies.end();
					++ it ){

				if( *it < lowest_energy ) lowest_energy = *it;
		}

		core::Real exp_sum = 0;
		core::Real avene = 0;
		core::Real kT = 1.0;
	 	for( std::vector< core::Real >::const_iterator it = energies.begin();
				 it != energies.end();
				 ++ it ){

			core::Real ediff = *it - lowest_energy;
			TR << ir << "  Emin: " << lowest_energy << "  " <<  " Ediff: " << ediff << "  " << exp( -ediff / kT ) << std::endl;
			exp_sum += exp( -ediff / kT );
			avene += ediff * exp( -ediff / kT );
		}
		avene /= sample_data.size();
		core::Real free_energy = -kT * log( exp_sum);
		core::Real enthalpy = avene;
		core::Real minus_T_times_entropy = free_energy - avene;

		TR << "ExpSum: "
		   << " ir : " << ir << " "
			 << " sz : " << sample_data.size() << " "
			 << " exp: " << exp_sum << " "
			 << " fe : " << free_energy << "  "
			 << " ent: " << enthalpy    << "  "
			 << " -TS: " << minus_T_times_entropy << "  "
			 << "   " << std::endl;

		total_free_energy            += free_energy;
  	total_enthalpy               += enthalpy;
  	total_minus_T_times_entropy  += minus_T_times_entropy;
	}

	pose = start_pose;
	(*score_function_)(pose);

	core::pose::setPoseExtraScore( pose, "free_energy", total_free_energy );
	core::pose::setPoseExtraScore( pose, "enthalpy",    total_enthalpy );
	core::pose::setPoseExtraScore( pose, "entropy",     total_minus_T_times_entropy );
	core::pose::setPoseExtraScore( pose, "total_rot",   total_rotamers );
	TR << "Total Rotamers: " << total_rotamers <<  std::endl;
}

std::string
ScTrials::get_name() const {
	return "sctrials";
}


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
    try {
    	using namespace protocols::jobdist;
    	using namespace protocols::moves;
    	using namespace scoring;

    	devel::init(argc, argv);

    	MoverOP protocol = new ScTrials();
    	protocols::jd2::JobDistributor::get_instance()->go( protocol );
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}
