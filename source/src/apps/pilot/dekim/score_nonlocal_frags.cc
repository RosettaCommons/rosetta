// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file David E Kim
/// @brief


// libRosetta headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/scoring/ResidualDipolarCoupling.hh>
#include <core/conformation/Conformation.hh>

// Project headers
#include <core/scoring/rms_util.hh>
#include <numeric/model_quality/rms.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/util.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>

#include <sstream>
#include <fstream>

//Auto Headers
#include <protocols/simple_filters/DdgFilter.hh>


static THREAD_LOCAL basic::Tracer TR( "main" );

using namespace protocols::moves;
using namespace core::scoring;
// using namespace basic::options;
// using namespace basic::options::OptionKeys;


class MyScoreMover : public Mover {
public:
	MyScoreMover();

	virtual void apply( core::pose::Pose& pose );
	std::string get_name() const { return "NonLocalFragsScoreMover"; }

	virtual MoverOP clone() const {
		return MoverOP( new MyScoreMover( *this ) );
	}

	virtual MoverOP fresh_instance() const {
		return MoverOP( new MyScoreMover );
	}

private:
	core::scoring::ScoreFunctionOP sfxn_;
	bool has_native_;
	core::pose::Pose native_;
	bool has_rdc_;
	core::scoring::ResidualDipolarCouplingOP rdc_file_;
	core::scoring::ResidualDipolarCoupling::RDC_lines rdc_raw_data_;

};

MyScoreMover::MyScoreMover():
	has_native_(false),
	has_rdc_(false)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;

	// get sfxn_ and add constraints if defined
	sfxn_ = scoring::get_score_function();

	if ( option[ in::file::native ].user() ) {
		// read native structure
		import_pose::pose_from_file( native_, option[ basic::options::OptionKeys::in::file::native ] , core::import_pose::PDB_file);
		has_native_ = true;
	}

	if ( option[ OptionKeys::in::file::rdc ].user() ) {
		rdc_file_ = core::scoring::ResidualDipolarCouplingOP( new core::scoring::ResidualDipolarCoupling );
		rdc_raw_data_ = rdc_file_->get_RDC_data();
		has_rdc_ = true;
	}
}

void MyScoreMover::apply( core::pose::Pose& pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;

	pose::Pose orig_pose = pose;

	std::string input_tag = protocols::jd2::get_current_job()->input_tag();
	// tag format example:  12asA_1_17_28_132
	utility::vector1< std::string >  tag = utility::string_split( input_tag, '_');
	if ( tag.size() != 5 ) {
		utility_exit_with_message( "Pose " + input_tag + " format error!" );
	}
	core::Size frag1_start = atoi(tag[2].c_str());
	core::Size frag2_start = atoi(tag[3].c_str());
	// fragment size should be pose total residue / 2
	if ( pose.size() % 2 != 0 ) {
		utility_exit_with_message( "Pose " + input_tag + " total_residue() % 2 != 0" );
	}
	Size frag_size = pose.size() / 2;

	TR.Debug << "scoring: " << input_tag << std::endl;

	if ( !pose.is_fullatom() ) {
		core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
	}

	pose.energies().clear();
	pose.data().clear();

	//   Detect disulfides
	pose.conformation().detect_disulfides();
	pose.conformation().detect_bonds();

	utility::vector1< bool > needToRepack( pose.size(), true );
	core::pack::task::PackerTaskOP taskstd = core::pack::task::TaskFactory::create_packer_task( pose );
	taskstd->restrict_to_repacking();
	taskstd->or_include_current(true);
	taskstd->restrict_to_residues( needToRepack );
	protocols::simple_moves::PackRotamersMover pack1( sfxn_, taskstd );
	pack1.apply( pose );

	// quick SC minimization
	core::optimization::AtomTreeMinimizer mzr;
	core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 1e-5, true, false );
	core::kinematics::MoveMapOP mm_min( new core::kinematics::MoveMap() );
	mm_min->set_bb( false );
	mm_min->set_chi( true );
	mzr.run( pose, *mm_min, *sfxn_, options );

	(*sfxn_)( pose );

	// RELAX non-local fragment pair
	// setup relax protocol for sub pose (frag pair)
	protocols::relax::RelaxProtocolBaseOP sub_pose_relax_protocol = protocols::relax::generate_relax_from_cmd();
	kinematics::MoveMapOP mm = sub_pose_relax_protocol->get_movemap();
	mm->set_jump(true); // set jumps movable
	pose::Pose relax_sub_pose = pose;
	sub_pose_relax_protocol->apply( relax_sub_pose );

	// check relaxed pose
	Real relaxed_rmsd = scoring::CA_rmsd( relax_sub_pose, pose );
	Real relaxed_score = (*sfxn_)( relax_sub_pose );
	pose::Pose relaxed_pose = relax_sub_pose;

	core::pose::setPoseExtraScore( pose, "nlf_rlx_rms", relaxed_rmsd );
	core::pose::setPoseExtraScore( pose, "nlf_rlx", relaxed_score );

	// DDG of relaxed fragment pair
	protocols::simple_filters::DdgFilter ddg = protocols::simple_filters::DdgFilter( 1000, sfxn_, 1, 5);
	Real relaxed_ddg_score = ddg.compute( relax_sub_pose );
	core::pose::setPoseExtraScore( pose, "nlf_ddg", relaxed_ddg_score );

	if ( has_rdc_ ) {
		Size weight = 1;
		Size max = 6*frag_size; // 3 rdcs per residue (maximum) should be given a weight of 1
		core::scoring::ResidualDipolarCoupling::RDC_lines rdc_data_given_fragment_lines;
		for ( core::scoring::ResidualDipolarCoupling::RDC_lines::const_iterator it = rdc_raw_data_.begin();
				it != rdc_raw_data_.end(); ++it ) {
			if ( it->res1() >= frag1_start && it->res1() <= frag1_start + frag_size - 1 && it->res2() >= frag1_start &&
					it->res2() <= frag1_start + frag_size - 1 ) {
				// first fragment of pair
				core::scoring::RDC selected ( it->res1() - frag1_start + 1, it->atom1(), it->res2() - frag1_start + 1,
					it->atom2(), it->Jdipolar());
				rdc_data_given_fragment_lines.push_back( selected );
			} else if ( it->res1() >= frag2_start && it->res1() <= frag2_start + frag_size - 1 &&
					it->res2() >= frag2_start && it->res2() <= frag2_start + frag_size - 1 ) {
				// second fragment of pair
				core::scoring::RDC selected ( it->res1() - frag2_start + 1 + frag_size, it->atom1(),
					it->res2() - frag2_start + 1 + frag_size, it->atom2(), it->Jdipolar());
				rdc_data_given_fragment_lines.push_back( selected );
			}
		}
		Real rdc_score = 0.0;
		if ( rdc_data_given_fragment_lines.size() < frag_size ) {
			rdc_score = 30.0; // if there's only 0.5 rdcs per residue, give it a big penalty
		} else {
			weight = ((-29/max)*rdc_data_given_fragment_lines.size()) + 30;  // linear weight
			if ( weight > 30 ) weight = 30;
			pose::Pose rdc_pose = orig_pose;
			core::scoring::ResidualDipolarCouplingOP rdc_data_given_fragment( new core::scoring::ResidualDipolarCoupling ( rdc_data_given_fragment_lines ) );
			//rdc_data_given_fragment->show(std::cout);
			//std::cout << std::endl;
			core::scoring::store_RDC_in_pose( rdc_data_given_fragment, rdc_pose );
			rdc_score = weight*rdc_data_given_fragment->compute_dipscore( rdc_pose );
		}
		core::pose::setPoseExtraScore( pose, "rdc", rdc_score );
		core::pose::setPoseExtraScore( pose, "rdcs", rdc_data_given_fragment_lines.size() );
	}

	// rmsd to native
	if ( has_native_ ) {
		// get pose CA coords
		std::vector< core::Vector > orig_pose_coords;
		std::vector< core::Vector > relaxed_pose_coords;
		std::vector< core::Vector > native_pose_coords;
		for ( Size i=1; i<=orig_pose.size(); i++ ) {
			orig_pose_coords.push_back( orig_pose.residue(i).xyz("CA") );
		}
		for ( Size i=1; i<=relaxed_pose.size(); i++ ) {
			relaxed_pose_coords.push_back( relaxed_pose.residue(i).xyz("CA") );
		}
		for ( Size i=0; i<frag_size; i++ ) {
			// get pose residue CA
			Size respos = frag1_start+i;
			native_pose_coords.push_back( native_.residue(respos).xyz("CA") );
		}
		for ( Size i=0; i<frag_size; i++ ) {
			// get pose residue CA
			Size respos = frag2_start+i;
			native_pose_coords.push_back( native_.residue(respos).xyz("CA") );
		}
		int const natoms = orig_pose_coords.size();
		ObjexxFCL::FArray2D< core::Real > p1a( 3, natoms );  // orig pose
		ObjexxFCL::FArray2D< core::Real > p2a( 3, natoms );  // relaxed pose
		ObjexxFCL::FArray2D< core::Real > p3a( 3, natoms );  // native pose
		for ( int i = 0; i < natoms; ++i ) {
			for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
				p1a(k+1,i+1) = orig_pose_coords[i][k];
				p2a(k+1,i+1) = relaxed_pose_coords[i][k];
				p3a(k+1,i+1) = native_pose_coords[i][k];
			}
		}
		// calculate rms of native to original pose
		core::Real rms_orig_native = numeric::model_quality::rms_wrapper( natoms, p1a, p3a );
		core::Real rms_relax_native = numeric::model_quality::rms_wrapper( natoms, p2a, p3a );
		core::pose::setPoseExtraScore( pose, "rms", rms_orig_native );
		core::pose::setPoseExtraScore( pose, "rms_rlx", rms_relax_native );
	}
}

int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols;
		using namespace protocols::jd2;

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core;

		jd2::register_options();

		// initialize core
		devel::init(argc, argv);

		//MyScoreMover* scoremover = new MyScoreMover;
		MoverOP scoremover( new MyScoreMover );

		using namespace protocols::jd2;

		// Make sure the default JobOutputter is SilentJobOutputter to ensure that when score_jd2
		// is called with default arguments is prints a proper scorefile and not the hacky thing that
		// the  JobOutputter scorefile() function produces (which for example skips Evaluators!!)

		// Set up a job outputter that writes a scorefile and no PDBs and no Silent Files.
		SilentFileJobOutputterOP jobout( new SilentFileJobOutputter );
		jobout->set_write_no_structures();
		jobout->set_write_separate_scorefile(true);

		// If the user chooses something else, then so be it, but by default score(_jd2) should only create a score
		// file and nothing else.
		protocols::jd2::JobDistributor::get_instance()->set_job_outputter( JobDistributorFactory::create_job_outputter( jobout ));

		JobDistributor::get_instance()->go( scoremover );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
		std::cout << "Exception: " << std::endl;
		excn.show( std::cout ); //so its also seen in a >LOG file
		return -1;
	}
	return 0;
}
