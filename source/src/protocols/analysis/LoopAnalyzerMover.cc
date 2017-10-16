// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/analysis/LoopAnalyzerMover.cc
/// @brief LoopAnalyzerMover implementation - class for in-depth loop structure quality analysis
/// @author Steven Lewis

// Unit Headers
#include <protocols/analysis/LoopAnalyzerMover.hh>
#include <protocols/analysis/LoopAnalyzerMoverCreator.hh>

// Package Headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Project Headers
#include <core/kinematics/FoldTree.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/symmetry/util.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/VariantType.hh>

#include <protocols/jd2/util.hh>

#include <protocols/loop_modeling/utilities/rosetta_scripts.hh>

// Utility Headers
#include <ObjexxFCL/FArray1D.hh> //necessary for fold tree tricks
#include <ObjexxFCL/FArray2D.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

// C++ Headers
#include <sstream>
#include <iomanip>
#include <limits>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.analysis.LoopAnalyzerMover" );

namespace protocols {
namespace analysis {

///stupid helper function needed because ternary operator does not allow variable return types
std::ostream & which_ostream( std::ostream & ost, std::ostream & oss, bool const tracer){
	if ( tracer ) return ost;
	return oss;
}

//This should be std::numeric_limits<core::Real>::lowest(), but that's cxx11
core::Real const REAL_FAKE_MIN = -1000000;


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////LoopAnalyzerMover////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
LoopAnalyzerMover::LoopAnalyzerMover( protocols::loops::Loops const & loops, bool const tracer ) :
	Mover(),
	loops_(protocols::loops::LoopsCOP( new protocols::loops::Loops(loops) )),
	tracer_(tracer),
	sf_(/* NULL */),
	chbreak_sf_(/* NULL */),
	total_score_(0),
	max_rama_(REAL_FAKE_MIN),
	max_chainbreak_(REAL_FAKE_MIN),
	max_omega_(REAL_FAKE_MIN),
	max_pbond_(REAL_FAKE_MIN)
{
	protocols::moves::Mover::type( "LoopAnalyzer" );
	set_sf();
}

LoopAnalyzerMover::~LoopAnalyzerMover() = default;

//Isn't there some fancy way to have ctors call each other to deduplicate this code?
LoopAnalyzerMover::LoopAnalyzerMover() :
	Mover(),
	tracer_(true), // WAG, who knows what user wants
	sf_(/* NULL */),
	chbreak_sf_(/* NULL */),
	total_score_(0),
	max_rama_(REAL_FAKE_MIN),
	max_chainbreak_(REAL_FAKE_MIN),
	max_omega_(REAL_FAKE_MIN),
	max_pbond_(REAL_FAKE_MIN)
{
	protocols::moves::Mover::type( "LoopAnalyzer" );
	set_sf();
}

LoopAnalyzerMover::LoopAnalyzerMover( LoopAnalyzerMover const & rhs ) :
	//utility::pointer::ReferenceCount(),
	Mover()
{
	*this = rhs;
}

// lhs.scorefxn_minimizer_        = rhs.scorefxn_minimizer_       ;


LoopAnalyzerMover & LoopAnalyzerMover::operator=( LoopAnalyzerMover const & rhs ) {
	if ( this == &rhs ) return *this;

	loops_ = protocols::loops::LoopsCOP( new protocols::loops::Loops(*(rhs.loops_) ) );
	tracer_ = rhs.tracer_;
	positions_ = rhs.positions_; //this is useless data
	sf_ = rhs.sf_->clone();
	chbreak_sf_ = rhs.chbreak_sf_->clone();
	scores_ = rhs.scores_; //useless
	total_score_ = rhs.total_score_;
	max_rama_ = rhs.max_rama_;
	max_chainbreak_ = rhs.max_chainbreak_;
	max_omega_ = rhs.max_omega_;
	max_pbond_ = rhs.max_pbond_;

	return *this;
}


void LoopAnalyzerMover::set_sf(){
	using namespace core::scoring;

	sf_ = ScoreFunctionOP( new ScoreFunction );
	sf_->set_weight( rama, 1.0 );
	sf_->set_weight( omega, 1.0 );
	sf_->set_weight( fa_dun, 1.0 );
	//sf_->set_weight( p_aa_pp, 1.0 );
	//sf_->set_weight( chainbreak, 1.0 );
	sf_->set_weight( peptide_bond, 1.0 );

	chbreak_sf_ = ScoreFunctionOP( new ScoreFunction );
	chbreak_sf_->set_weight( chainbreak, 20.0 );

	return;
}

/// @brief reset stored data
void LoopAnalyzerMover::reset() {
	scores_.clear();
	max_rama_ = REAL_FAKE_MIN;
	max_chainbreak_ = REAL_FAKE_MIN;
	max_omega_ = REAL_FAKE_MIN;
	max_pbond_ = REAL_FAKE_MIN;
}

/// @details LoopAnalyzerMover
void LoopAnalyzerMover::apply( core::pose::Pose & input_pose )
{
	TR << "running LoopAnalyzerMover" << std::endl;

	reset();

	//prep pose
	core::pose::Pose pose;
	//if symmetric, desymmetrize; otherwise just make a local copy
	if ( !basic::options::option[ basic::options::OptionKeys::symmetry::symmetry_definition ].user() ) {
		pose = input_pose; //protecting input pose from our chainbreak changes (and its energies object)
	} else {
		core::pose::symmetry::extract_asymmetric_unit(input_pose, pose, false); //don't need VIRTs
	}

	loops_->verify_against(pose);

	find_positions(pose);
	calculate_all_chainbreaks(pose);
	(*sf_)(pose);

	//make output
	std::ostringstream results_oss;
	std::ostream & results = which_ostream(TR, results_oss, tracer_); //easy swap between tracer/job output
	core::Real total_rama(0), total_omega(0), total_peptide_bond(0), total_chbreak(0);

	results << "LoopAnalyzerMover: unweighted bonded terms and angles (in degrees)" << std::endl
		<< "position phi_angle psi_angle omega_angle peptide_bond_C-N_distance rama_score omega_score dunbrack_score peptide_bond_score chainbreak_score" << std::endl
		<< " pos phi_ang psi_ang omega_ang pbnd_dst    rama  omega_sc dbrack pbnd_sc   cbreak" << std::endl;

	for ( core::Size i(1); i <= positions_.size(); ++i ) {
		core::Size const res(positions_[i]);
		//peptide_bond distance
		core::Real const pbnd_dist(input_pose.residue(res).atom("C").xyz().distance(input_pose.residue(res+1).atom("N").xyz()));
		//phi/psi/omega
		//  results.precision(3);
		results << std::setw(4) << res
			<< std::setw(8) << std::setprecision(4) << pose.phi(res)
			<< std::setw(8) << std::setprecision(4) << pose.psi(res)
			<< std::setw(10) << std::setprecision(4) << pose.omega(res)
			<< std::setw(9) << std::setprecision(4) << pbnd_dist;

		using namespace core::scoring;
		EnergyMap const & emap(pose.energies().residue_total_energies(res));
		results << std::setw(8) << std::setprecision(3) << emap[rama]
			//<< std::setw(8) << std::setprecision(4) << emap[p_aa_pp]
			<< std::setw(10) << std::setprecision(3) << emap[omega]
			<< std::setw(7) << std::setprecision(3) << emap[fa_dun]
			<< std::setw(8) << std::setprecision(3) << emap[peptide_bond]
			<< std::setw(9) << std::setprecision(3) << scores_[i]
			<< std::setprecision(6) << std::endl;

		total_rama += emap[rama];
		max_rama_=std::max(max_rama_, emap[rama]);

		total_omega += emap[omega];
		max_omega_=std::max(max_omega_, emap[omega]);

		total_peptide_bond += emap[peptide_bond];
		max_pbond_=std::max(max_pbond_, emap[peptide_bond]);

		total_chbreak += scores_[i];
		max_chainbreak_=std::max(max_chainbreak_,scores_[i]);
	}//for all loop positions

	results << "total_rama " << total_rama << std::endl;
	results << "total_omega " << total_omega << std::endl;
	results << "total_peptide_bond " << total_peptide_bond << std::endl;
	results << "total_chainbreak " << total_chbreak << std::endl;
	total_score_ = total_rama + total_omega + total_peptide_bond + total_chbreak;
	results << "total rama+omega+peptide bond+chainbreak " << total_score_ << std::endl;


	if ( !tracer_ ) {
		protocols::jd2::add_string_to_current_job(results_oss.str());
		//store the loop_total where jd2 silent file can get it
		protocols::jd2::add_string_real_pair_to_current_job("LAM_total", total_score_);
	}

	return;
}//LoopAnalyzerMover::apply


void
LoopAnalyzerMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	set_use_tracer(tag->getOption< bool >( "use_tracer", false ) );
	set_loops(protocols::loop_modeling::utilities::parse_loops_from_tag(tag));
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
LoopAnalyzerMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new LoopAnalyzerMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
LoopAnalyzerMover::clone() const
{
	return protocols::moves::MoverOP( new LoopAnalyzerMover( *this ) );
}

// XRW TEMP std::string
// XRW TEMP LoopAnalyzerMover::get_name() const {
// XRW TEMP  return "LoopAnalyzerMover";
// XRW TEMP }

//////////////////////////getters, setters/////////////////////
/// @brief set loops object, because public setters/getters are a rule
void LoopAnalyzerMover::set_loops( protocols::loops::LoopsCOP loops ) { loops_ = loops; }

/// @brief get loops object, because public setters/getters are a rule
protocols::loops::LoopsCOP const & LoopAnalyzerMover::get_loops( void ) const { return loops_; }

core::Real
LoopAnalyzerMover::get_total_score() const {
	return total_score_;
}

core::Real
LoopAnalyzerMover::get_max_rama() const {
	return max_rama_;
}

core::Real
LoopAnalyzerMover::get_max_omega() const {
	return max_omega_;
}
core::Real
LoopAnalyzerMover::get_max_pbond() const {
	return max_pbond_;
}

core::Real
LoopAnalyzerMover::get_max_chainbreak() const {
	return max_chainbreak_;
}

utility::vector1<core::Real>
LoopAnalyzerMover::get_chainbreak_scores() {
	return scores_;
}

void LoopAnalyzerMover::find_positions( core::pose::Pose const & pose ){
	positions_.clear();
	for ( auto const & it : *loops_ ) {
		core::Size start(it.start());
		core::Size end(it.stop());
		//in here we are about to grow the analysis area by moving start and end each "out" one residue, unless we hit cases where we cannot/should not (because there is not residue there), or we're at a terminus (in which case the next residue is not the same chain), or we're adjacent to the terminus (in which case the chainbreak on the terminus cannot be calculated).
		if ( !pose.residue(start).is_terminus() //if this is not a terminus,
				&& start != 1 //AND there are residues before it,
				/*&& !pose.residue(start-1).is_terminus()*/ ) { //if the previous residue is not a terminus (actually ok)
			--start;   //grow the analysis area
		}
		if ( pose.residue(end).is_terminus() //if this is a terminus,
				|| end == pose.size() ) { //or if there are no residues after it,
			--end; //move end IN one, so that the i+1 activity in add/remove variant type in calculate_all_chainbreaks won't fail
		} else if ( !pose.residue(end+1).is_terminus()  //if the next residue is not a terminus (note this end+1 is a safe check because if end=pose.size we triggered the preceding if instead, and this is never checked)
				&& pose.size() != (end+1) ) { //and the next residue is not the end of the pose
			++end;  //grow the analysis area
		} //implicitly - if "end" is not a terminus but end+1 is - do nothing.
		TR << "LoopAnalyzerMover will consider these positions (Rosetta numbering) - remember that it includes an extra residue on both sides of each loop, conditions permitting: ";
		for ( core::Size i(start); i <= end; ++i ) {
			TR << i << " ";
			positions_.push_back(i);
		}
		TR << std::endl;
	}//for all loops
	scores_.clear();
	scores_.resize(positions_.size(), 0);
	return;
}//find_positions

void LoopAnalyzerMover::calculate_all_chainbreaks( core::pose::Pose & pose )
{
	using namespace core::chemical;
	using core::kinematics::Edge;

	//remove all chainbreak variants within loop
	core::Size const numpos(positions_.size());
	for ( core::Size i(1); i<=numpos; ++i ) {
		core::Size const pos = positions_[i];
		core::pose::remove_variant_type_from_pose_residue(pose, CUTPOINT_LOWER, pos);
		core::pose::remove_variant_type_from_pose_residue(pose, CUTPOINT_UPPER, pos);
	}

	//create a chainbreak at the desired position, score it, store the score, and revert changes
	for ( core::Size i(1); i<=numpos; ++i ) {
		core::Size const pos = positions_[i];

		//create chainbreak variant
		core::pose::add_variant_type_to_pose_residue(pose, CUTPOINT_LOWER, pos);
		runtime_assert_string_msg((pose.size() >= pos+1), "LoopAnalyzerMover: cannot analyze the last position of a pose, as residue i+1 does not exist");
		core::pose::add_variant_type_to_pose_residue(pose, CUTPOINT_UPPER, pos+1);

		//create cut in fold tree
		core::kinematics::FoldTree ft(pose.size());
		ft.clear(); //remove initial edge
		ft.add_edge(Edge(1, pos, Edge::PEPTIDE));
		ft.add_edge(Edge(pos+1, pose.size(), Edge::PEPTIDE));
		ft.add_edge(Edge(pos, pos+1, 1)); //a jump and therefore a cut
		ft.delete_self_edges(); //needed if loop ends adjacent to a pose's end
		ft.reorder(1);
		runtime_assert_string_msg(ft.check_fold_tree(), "FoldTree in LoopAnalyzerMover invalid: aborting");
		pose.fold_tree(ft);

		//score the sucker
		scores_[i] = (*chbreak_sf_)(pose);

		//undo that chainbreak variant
		core::pose::remove_variant_type_from_pose_residue(pose, CUTPOINT_LOWER, pos);
		core::pose::remove_variant_type_from_pose_residue(pose, CUTPOINT_UPPER, pos+1);
	}

	//revert to innocuous fold tree
	pose.fold_tree(core::kinematics::FoldTree(pose.size()));

	return;
}//calculate_all_chainbreaks


////////////////////creator/////////////////////////////

// XRW TEMP std::string
// XRW TEMP LoopAnalyzerMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return LoopAnalyzerMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP LoopAnalyzerMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new LoopAnalyzerMover() );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP LoopAnalyzerMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "LoopAnalyzerMover";
// XRW TEMP }

std::string LoopAnalyzerMover::get_name() const {
	return mover_name();
}

std::string LoopAnalyzerMover::mover_name() {
	return "LoopAnalyzerMover";
}

std::string LoopAnalyzerMover_subelement_ct_name( std::string const & name ) {
	return "LoopAnalyzerMover_Loop_subelement_" + name + "Type";
}

void LoopAnalyzerMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default( "use_tracer", xsct_rosetta_bool, "report output to a Tracer", "false" );

	//empty SubelementList for MoveMap util function
	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( & LoopAnalyzerMover_subelement_ct_name );
	loop_modeling::utilities::append_subelement_and_attributes_for_parse_loops_from_tag(xsd, subelements, attlist);

	moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd,
		mover_name(),
		"This mover reports scores and statistics useful for judging the quality of loops",
		attlist,
		subelements);

	//protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "This mover reports scores and statistics useful for judging the quality of loops", attlist );
}

std::string LoopAnalyzerMoverCreator::keyname() const {
	return LoopAnalyzerMover::mover_name();
}

protocols::moves::MoverOP
LoopAnalyzerMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoopAnalyzerMover );
}

void LoopAnalyzerMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoopAnalyzerMover::provide_xml_schema( xsd );
}



}//analysis
}//protocols
