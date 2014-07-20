// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief protocols for folding into density
/// @detailed
/// @author Frank DiMaio

// keep first
#include <core/scoring/cryst/PhenixInterface.hh>

#include <protocols/cryst/cryst_movers.hh>
#include <protocols/cryst/cryst_movers_creator.hh>

#include <protocols/comparative_modeling/coord_util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>

#include <core/util/cryst_util.hh>
#include <core/scoring/cryst/util.hh>
#include <core/scoring/cryst/XtalMLEnergy.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDB_Info.hh>
#include <core/pose/Remarks.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>

/////////////
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/types.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/AtomTreeMultifunc.hh>

#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/symmetry/SymMinimizerMap.hh>
#include <core/optimization/symmetry/SymAtomTreeMultifunc.hh>

#include <core/optimization/CartesianMinimizerMap.hh>
#include <core/optimization/CartesianMultifunc.hh>
#include <core/optimization/CartesianMinimizer.hh>

///////////

#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>
#include <numeric/model_quality/rms.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>


using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace cryst {

static basic::Tracer TR("protocols.cryst.cryst_movers");

using namespace protocols;
using namespace core;
  using namespace kinematics;
  using namespace scoring;


//////////////////
//////////////////
/// creators


std::string
SetCrystWeightMoverCreator::keyname() const {
	return SetCrystWeightMoverCreator::mover_name();
}

protocols::moves::MoverOP
SetCrystWeightMoverCreator::create_mover() const {
	return new SetCrystWeightMover;
}

std::string
SetCrystWeightMoverCreator::mover_name() {
	return "SetCrystWeight";
}

///

std::string
RecomputeDensityMapMoverCreator::keyname() const {
	return RecomputeDensityMapMoverCreator::mover_name();
}

protocols::moves::MoverOP
RecomputeDensityMapMoverCreator::create_mover() const {
	return new RecomputeDensityMapMover;
}

std::string
RecomputeDensityMapMoverCreator::mover_name() {
	return "RecomputeDensityMap";
}

///

std::string
LoadDensityMapMoverCreator::keyname() const {
	return LoadDensityMapMoverCreator::mover_name();
}

protocols::moves::MoverOP
LoadDensityMapMoverCreator::create_mover() const {
	return new LoadDensityMapMover;
}

std::string
LoadDensityMapMoverCreator::mover_name() {
	return "LoadDensityMap";
}

///

std::string
FitBfactorsMoverCreator::keyname() const {
	return FitBfactorsMoverCreator::mover_name();
}

protocols::moves::MoverOP
FitBfactorsMoverCreator::create_mover() const {
	return new FitBfactorsMover;
}

std::string
FitBfactorsMoverCreator::mover_name() {
	return "FitBfactors";
}

///

std::string
UpdateSolventMoverCreator::keyname() const {
	return UpdateSolventMoverCreator::mover_name();
}

protocols::moves::MoverOP
UpdateSolventMoverCreator::create_mover() const {
	return new UpdateSolventMover;
}

std::string
UpdateSolventMoverCreator::mover_name() {
	return "UpdateSolvent";
}

///

std::string
TagPoseWithRefinementStatsMoverCreator::keyname() const {
	return TagPoseWithRefinementStatsMoverCreator::mover_name();
}

protocols::moves::MoverOP
TagPoseWithRefinementStatsMoverCreator::create_mover() const {
	return new TagPoseWithRefinementStatsMover;
}

std::string
TagPoseWithRefinementStatsMoverCreator::mover_name() {
	return "TagPoseWithRefinementStats";
}

///

std::string
SetRefinementOptionsMoverCreator::keyname() const {
	return SetRefinementOptionsMoverCreator::mover_name();
}

protocols::moves::MoverOP
SetRefinementOptionsMoverCreator::create_mover() const {
	return new SetRefinementOptionsMover;
}

std::string
SetRefinementOptionsMoverCreator::mover_name() {
	return "SetRefinementOptions";
}

//////////////////
//////////////////
/// movers

void SetCrystWeightMover::apply( core::pose::Pose & pose ) {
	if (weight_ < 0) {
		using namespace core::optimization;
		using namespace core::optimization::symmetry;

		// report!
		core::scoring::ScoreFunctionOP rosetta_scorefxn = new core::scoring::ScoreFunction();
		core::scoring::ScoreFunctionOP xtal_scorefxn = new core::scoring::ScoreFunction();

		rosetta_scorefxn->set_weight( core::scoring::xtal_ml, 0.0 );
		xtal_scorefxn->set_weight( core::scoring::xtal_ml, 1.0 );

		core::kinematics::MoveMap move_map;
		move_map.set_bb  ( true );
		move_map.set_chi ( true );
		move_map.set_jump( true );
		if (core::pose::symmetry::is_symmetric(pose)) {
			core::conformation::symmetry::SymmetricConformation const & symm_conf (
					dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation() ) );
			core::conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

			// symmetrize scorefunct & movemap
			rosetta_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *rosetta_scorefxn );
			xtal_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *xtal_scorefxn );
			core::pose::symmetry::make_symmetric_movemap( pose, move_map );
		}

		// compute gradients using both scorefunctions
		CartesianMinimizerMap min_map;
		min_map.setup( pose, move_map );
		Multivec vars( min_map.ndofs() ), dExtal_dvars;
		min_map.copy_dofs_from_pose( pose, vars );

		(*xtal_scorefxn)(pose);  // score pose first
		xtal_scorefxn->setup_for_minimizing( pose, min_map );
		CartesianMultifunc f_xtal( pose, min_map, *xtal_scorefxn, false, false );
		f_xtal.dfunc( vars, dExtal_dvars );

		utility::vector1< Multivec > dEros_dvars(1);
		for (int ii=0; ii<1; ++ii) {
			rosetta_scorefxn->set_weight( core::scoring::fa_atr       , (ii<=1 || ii==2)? score_function_ref_->get_weight(core::scoring::fa_atr) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::fa_rep       , (ii<=1 || ii==3)? score_function_ref_->get_weight(core::scoring::fa_rep) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::fa_sol       , (ii<=1 || ii==4)? score_function_ref_->get_weight(core::scoring::fa_sol) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::fa_intra_rep , (ii<=1 || ii==5)? score_function_ref_->get_weight(core::scoring::fa_intra_rep) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::fa_pair      , (ii<=1 || ii==6)? score_function_ref_->get_weight(core::scoring::fa_pair) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::fa_dun       , (ii<=1 || ii==7)? score_function_ref_->get_weight(core::scoring::fa_dun) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::hbond_lr_bb  , (ii<=1 || ii==8)? score_function_ref_->get_weight(core::scoring::hbond_lr_bb) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::hbond_sr_bb  , (ii<=1 || ii==9)? score_function_ref_->get_weight(core::scoring::hbond_sr_bb) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::hbond_bb_sc  , (ii<=1 || ii==10)? score_function_ref_->get_weight(core::scoring::hbond_bb_sc) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::hbond_sc     , (ii<=1 || ii==11)? score_function_ref_->get_weight(core::scoring::hbond_sc) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::p_aa_pp      , (ii<=1 || ii==12)? score_function_ref_->get_weight(core::scoring::p_aa_pp) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::dslf_ss_dst  , (ii<=1 || ii==13)? score_function_ref_->get_weight(core::scoring::dslf_ss_dst) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::dslf_cs_ang  , (ii<=1 || ii==14)? score_function_ref_->get_weight(core::scoring::dslf_cs_ang) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::dslf_ss_dih  , (ii<=1 || ii==15)? score_function_ref_->get_weight(core::scoring::dslf_ss_dih) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::dslf_ca_dih  , (ii<=1 || ii==16)? score_function_ref_->get_weight(core::scoring::dslf_ca_dih) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::pro_close    , (ii<=1 || ii==17)? score_function_ref_->get_weight(core::scoring::pro_close) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::rama         , (ii<=1 || ii==18)? score_function_ref_->get_weight(core::scoring::rama) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::omega        , (ii<=1 || ii==19)? score_function_ref_->get_weight(core::scoring::omega) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::cart_bonded_angle , (ii==0 || ii==20)?
				score_function_ref_->get_weight(core::scoring::cart_bonded_angle)+score_function_ref_->get_weight(core::scoring::cart_bonded) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::cart_bonded_length , (ii==0 || ii==21)?
				score_function_ref_->get_weight(core::scoring::cart_bonded_length)+score_function_ref_->get_weight(core::scoring::cart_bonded) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::cart_bonded_torsion , (ii==0 || ii==22)?
				score_function_ref_->get_weight(core::scoring::cart_bonded_torsion)+score_function_ref_->get_weight(core::scoring::cart_bonded) : 0.0 );

			(*rosetta_scorefxn)(pose);  // score pose first
			rosetta_scorefxn->setup_for_minimizing( pose, min_map );
			CartesianMultifunc f_ros( pose, min_map, *rosetta_scorefxn, false, false );
			f_ros.dfunc( vars, dEros_dvars[ii+1] );
		}

		// report
		// GRAD  name id atom_name total total_no_cartbond fa_atr fa_rep fa_sol fa_intra_rep fa_pair fa_dun hbond_lr_bb hbond_sr_bb hbond_bb_sc hbond_sc p_aa_pp dslf_ss_dst dslf_cs_ang dslf_ss_dih dslf_ca_dih pro_close rama omega cart_bonded_angle cart_bonded_length cart_bonded_torsion
		for (core::Size counter=0; counter<dEros_dvars[1].size()/3; ++counter) {
			id::AtomID id = min_map.get_atom( counter+1 );

			core::conformation::Residue const & rsd_i = pose.residue( id.rsd() );
			core::Size pdb_res = pose.pdb_info()->number( id.rsd() );
			char pdb_chain = pose.pdb_info()->chain( id.rsd() );

			if (id.atomno() <= rsd_i.natoms()) {
				std::cout << "gradient " << rsd_i.name3() << " " << pdb_res << " " << pdb_chain<< " " << rsd_i.atom_name( id.atomno() )
					<< ": " << dEros_dvars[1][3*counter+1] << "," << dEros_dvars[1][3*counter+2] << "," << dEros_dvars[1][3*counter+3]
					<< std::endl;
			}
		}
	} else if (autoset_wt_) {
		// if autoset_wt_ guess at cryst wt
		//    use score_function_ref_ to do the scaling
		core::Real auto_weight = 0;
		if (cartesian_) {
			if (mm_) {
				auto_weight = core::util::getMLweight_cart( *score_function_ref_, pose, *mm_  );
			} else {
				auto_weight = core::util::getMLweight_cart( *score_function_ref_, pose );
			}
		} else {
			if (mm_) {
				auto_weight = core::util::getMLweight( *score_function_ref_, pose, *mm_  );
			} else {
				auto_weight = core::util::getMLweight( *score_function_ref_, pose );
			}
		}
		score_function_->set_weight( core::scoring::xtal_ml, auto_weight * weight_scale_ );
		TR << "set xtal_weight to " << score_function_->get_weight( core::scoring::xtal_ml ) << std::endl;
	} else {
		TR << "override automatic weight with wt=" << weight_ * weight_scale_ << std::endl;
		score_function_->set_weight( core::scoring::xtal_ml, weight_ * weight_scale_ );
	}
}

void SetCrystWeightMover::parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap &data,
			filters::Filters_map const & /*filters*/,
			moves::Movers_map const & /*movers*/,
			core::pose::Pose const & pose ) {

	score_function_ = protocols::rosetta_scripts::parse_score_function( tag, data );

	if ( tag->hasOption("scorefxn_ref") ) {
		std::string const scorefxn_key( tag->getOption<std::string>("scorefxn_ref") );
		if ( ! data.has( "scorefxns", scorefxn_key ) )
			utility_exit_with_message("ScoreFunction " + scorefxn_key + " not found in basic::datacache::DataMap.");
		score_function_ref_ = data.get< ScoreFunction* >( "scorefxns", scorefxn_key );
	} else {
		score_function_ref_ = score_function_->clone();
	}

	if ( tag->hasOption("MoveMap") ) {
		protocols::rosetta_scripts::parse_movemap( tag, pose, mm_ );
	}

	// also can specify defaults
	if ( tag->hasOption("jump") ) {
		if ( ! mm_ ) mm_ = new MoveMap;
		utility::vector1<std::string> jumps = utility::string_split( tag->getOption<std::string>( "jump" ), ',' );
		// string 'ALL' makes all jumps movable
		if (jumps.size() == 1 && (jumps[1] == "ALL" || jumps[1] == "All" || jumps[1] == "all") ) {
				mm_->set_jump( true );
		} else {
			for( std::vector<std::string>::const_iterator it = jumps.begin(); it != jumps.end(); ++it ) {
				Size const value = std::atoi( it->c_str() ); // convert to C string, then convert to integer, then set a Size (phew!)
				mm_->set_jump( value, true );
			}
		}
	}
	if ( tag->hasOption("chi") ) {
		bool const value( tag->getOption<bool>("chi") );
		if ( ! mm_ ) mm_ = new MoveMap;
		mm_->set_chi(value);
	}
	if ( tag->hasOption("bb") ) {
		bool const value( tag->getOption<bool>("bb") );
		if ( ! mm_ ) mm_ = new MoveMap;
		mm_->set_bb(value);
	}
	if ( tag->hasOption("bondangle") ) {
		bool const value( tag->getOption<bool>("bondangle") );
		if ( ! mm_ ) mm_ = new MoveMap;
		mm_->set( core::id::THETA, value );
	}
	if ( tag->hasOption("bondlength") ) {
		bool const value( tag->getOption<bool>("bondlength") );
		if ( ! mm_ ) mm_ = new MoveMap;
		mm_->set( core::id::D, value );
	}

	if ( tag->hasOption("weight") ) {
		weight_ = tag->getOption<core::Real>("weight");
		autoset_wt_ = false;
	}
	weight_scale_ = tag->getOption<core::Real>("weight_scale", 1.0);
	cartesian_ = tag->getOption<bool>("cartesian", false);
}

////

void RecomputeDensityMapMover::apply( core::pose::Pose & pose )
{
	using namespace core::scoring::electron_density;
	std::string newMap = core::scoring::cryst::getPhenixInterface().calculateDensityMap( pose, !keep_sidechains_ );
	/*ElectronDensity& densityMap =*/ getDensityMap(newMap, true);
}

void RecomputeDensityMapMover::parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & /*data*/,
			filters::Filters_map const & /*filters*/,
			moves::Movers_map const & /*movers*/,
			core::pose::Pose const & /*pose*/ )
{
	keep_sidechains_ = tag->getOption<bool>("sidechains", true);
}

////
void LoadDensityMapMover::apply( core::pose::Pose & /*pose*/ ) {
	using namespace core::scoring::electron_density;

	ElectronDensity& densityMap = getDensityMap(mapfile_, true);

	//fpd make these selectable?
	densityMap.setWindow(window_);
	densityMap.setSCscaling(sc_scale_);
}

void LoadDensityMapMover::parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & /*data*/,
			filters::Filters_map const & /*filters*/,
			moves::Movers_map const & /*movers*/,
			core::pose::Pose const & /*pose*/ )
{
	mapfile_ = tag->getOption<std::string>("mapfile");
	sc_scale_ = tag->getOption<core::Real>("sc_scale", 1.0);
	window_ = tag->getOption<core::Size>("window", 3);
}

////

void FitBfactorsMover::apply( core::pose::Pose & pose ) {
	if (adp_strategy_ == "randomize") {
		randomize_bs( pose );
	} else {
		core::scoring::cryst::getPhenixInterface().set_adp_strategy( adp_strategy_ );
		core::scoring::cryst::getPhenixInterface().fitBfactors( pose );
	}

	// reinitialize fmodel object by rescoring pose
	core::scoring::cryst::getPhenixInterface().getScore( pose );
}

void FitBfactorsMover::parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & /*data*/,
			filters::Filters_map const & /*filters*/,
			moves::Movers_map const & /*movers*/,
			core::pose::Pose const & /*pose*/ )
{
	std::string adp_strat = tag->getOption<std::string>("adp_strategy", "individual");
	runtime_assert( adp_strat=="individual" || adp_strat=="group" || adp_strat=="randomize" );

	if ( adp_strat=="group") {
		Size group_adp_strat = tag->getOption<Size>("group_adp_mode", 1);
		runtime_assert( group_adp_strat==1 || group_adp_strat==2 );
		if (group_adp_strat==1)
			adp_strategy_ = "group1";
		else
			adp_strategy_ = "group2";
	} else if ( adp_strat=="randomize") {
		adp_strategy_ = "randomize";
		b_min_ = tag->getOption<core::Real>("b_min", 5.0);
		b_max_ = tag->getOption<core::Real>("b_max", 5.0);
	} else {
		adp_strategy_ = "individual";
	}
}

void FitBfactorsMover::randomize_bs( core::pose::Pose & pose ) {
	for (core::Size resid=1; resid<=pose.total_residue(); ++resid) {
		if (pose.residue(resid).aa() == core::chemical::aa_vrt) continue;
		core::conformation::Residue const &rsd_i = pose.residue(resid);
		for (core::Size atmid=1; atmid<=rsd_i.natoms(); ++atmid) {
			core::Real B = numeric::random::uniform()*(b_max_-b_min_) + b_min_;
			pose.pdb_info()->temperature( resid, atmid, B);
		}
	}
}

////

void UpdateSolventMover::apply( core::pose::Pose & pose ) {
	// make sure fmodel is initialized by scoring the pose
	core::scoring::cryst::getPhenixInterface().getScore( pose );

	//fcalc
	if( update_fcalc_ ) {
		core::scoring::cryst::getPhenixInterface().updateFcalc();
	}

	//mask
	if (optimize_mask_ && optimize_params_) {
		core::scoring::cryst::getPhenixInterface().optimizeSolvParamsAndMask();
	} else if (optimize_mask_) {
		core::scoring::cryst::getPhenixInterface().optimizeSolventMask();
	} else if( optimize_params_ ) {
		core::scoring::cryst::getPhenixInterface().optimizeSolvParams();
	} else if( update_mask_ ) {
		core::scoring::cryst::getPhenixInterface().updateSolventMask();
	}
}

void UpdateSolventMover::parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & /*data*/,
			filters::Filters_map const & /*filters*/,
			moves::Movers_map const & /*movers*/,
			core::pose::Pose const & /*pose*/ )
{
	update_mask_ = tag->getOption<bool>("update_mask", 1);
	update_fcalc_ = tag->getOption<bool>("update_fcalc", 1);
	optimize_mask_ = tag->getOption<bool>("optimize_mask", 0);
	optimize_params_ = tag->getOption<bool>("optimize_params", 0);
}

////

void TagPoseWithRefinementStatsMover::apply( core::pose::Pose & pose ) {
	// make sure fmodel is initialized by scoring the pose
	/*core::Real xraytgt =*/ core::scoring::cryst::getPhenixInterface().getScore( pose );

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	//fpd if the pose is symmetric use a symmetric scorefunction
	if (core::pose::symmetry::is_symmetric(pose))
		scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *scorefxn );
	core::Real score = (*scorefxn)(pose);
	core::pose::RemarkInfo remark;

	std::ostringstream oss;

	oss << "[ " << tag_ << " ] R=" << core::scoring::cryst::getPhenixInterface().getInfoLine();

	//oss.str("");
	//oss << "[ " << tag_ << " ] xray_tgt = " << xraytgt << "   score = " << score;
	oss << " sc=" << score;

	core::pose::Pose pose_asu;
	if (basic::options::option[ basic::options::OptionKeys::in::file::native ].user() || dump_pose_) {
		core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);
	}

	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		// rms to native
		if (!get_native_pose()) protocols::jd2::set_native_in_mover(*this);
		core::pose::PoseCOP native = get_native_pose();

		core::sequence::SequenceOP model_seq ( new core::sequence::Sequence( pose_asu.sequence(),  "model",  1 ) );
		core::sequence::SequenceOP native_seq( new core::sequence::Sequence( native->sequence(), "native", 1 ) );
		core::sequence::SequenceAlignment aln = align_naive(model_seq,native_seq);

		int n_atoms;
		ObjexxFCL::FArray2D< core::Real > p1a, p2a;
		protocols::comparative_modeling::gather_coords( pose_asu, *native, aln, n_atoms, p1a, p2a );
		core::Real rms = numeric::model_quality::rms_wrapper( n_atoms, p1a, p2a );
		oss << " rms=" << rms;
	}

	TR << oss.str() << std::endl;
	remark.num = 1;	remark.value = oss.str();
	pose.pdb_info()->remarks().push_back( remark );

	if (dump_pose_) {
		std::string out_name = protocols::jd2::JobDistributor::get_instance()->current_output_name() + "_" + tag_ + ".pdb";
		pose_asu.dump_pdb( out_name );
	}
}


void TagPoseWithRefinementStatsMover::parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & /*data*/,
			filters::Filters_map const & /*filters*/,
			moves::Movers_map const & /*movers*/,
			core::pose::Pose const & /*pose*/ )
{
	tag_ = tag->getOption<std::string>("tag", "");
	dump_pose_ = tag->getOption<bool>("dump", 0);
}

////
void SetRefinementOptionsMover::apply( core::pose::Pose & /*pose*/ ) {
	TR << "Setting options:" << std::endl;

	if (res_high_ != 0 || res_low_ != 0) {
		TR << "  res_high:" << res_high_ << std::endl;
		TR << "  res_low:" << res_low_ << std::endl;
		core::scoring::cryst::getPhenixInterface().setResLimits(res_high_, res_low_);
	}
	if (twin_law_.length() != 0) {
		TR << "  twin_law:" << twin_law_ << std::endl;
		core::scoring::cryst::getPhenixInterface().setTwinLaw(twin_law_);
	}
	if (algo_.length() != 0) {
		TR << "  algorithm:" << algo_ << std::endl;
		core::scoring::cryst::getPhenixInterface().setAlgorithm(algo_);
	}
	if (target_.length() != 0) {
		TR << "  target:" << std::endl;
		core::scoring::cryst::getPhenixInterface().set_target_function(target_);
	}

	if (setmap_type_) {
		TR << "  map_type:" << map_type_ << std::endl;
		core::scoring::cryst::getPhenixInterface().set_map_type ( map_type_ );
	}

	if (cif_files_.size() > 0) {
		TR << "Passing cif files to phenix refine:" << std::endl;
		for (core::Size i=1; i<=cif_files_.size(); ++i) {
			TR << "  " << cif_files_[i] << std::endl;
		}
		core::scoring::cryst::getPhenixInterface().set_cif_files ( cif_files_ );
	}
}


void SetRefinementOptionsMover::parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & /*data*/,
			filters::Filters_map const & /*filters*/,
			moves::Movers_map const & /*movers*/,
			core::pose::Pose const & /*pose*/ )
{
	res_high_ = tag->getOption<core::Real>("res_high", 0.0);
	res_low_ = tag->getOption<core::Real>("res_low", 0.0);
	twin_law_ = tag->getOption<std::string>("twin_law", "");
	algo_ = tag->getOption<std::string>("algorithm", "");
	target_ = tag->getOption<std::string>("target", "");
	setmap_type_ = tag->hasOption("map_type");
	map_type_ = tag->getOption<std::string>("map_type", "");

	std::string allcifs = tag->getOption<std::string>("cifs", "");
	if (allcifs != "")
		cif_files_ = utility::string_split( allcifs, ',' );
}

}
}
