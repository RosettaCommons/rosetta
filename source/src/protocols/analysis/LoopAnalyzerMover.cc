// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/analysis/LoopAnalyzerMover.cc
/// @brief LoopAnalyzerMover implementation - class for in-depth loop structure quality analysis
/// @author Steven Lewis

// Unit Headers
#include <protocols/analysis/LoopAnalyzerMover.hh>

// Package Headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Project Headers
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/pose/Pose.hh>


#include <core/chemical/VariantType.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh> //CENTROID, FA_STANDARD

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Utility Headers
#include <ObjexxFCL/FArray1D.hh> //necessary for fold tree tricks
#include <ObjexxFCL/FArray2D.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <numeric/xyz.io.hh>

// C++ Headers
#include <sstream>
#include <iomanip>

#include <core/chemical/AtomType.hh>
#include <core/pose/util.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.analysis.LoopAnalyzerMover" );

namespace protocols{
namespace analysis{

///stupid helper function needed because ternary operator does not allow variable return types
std::ostream & which_ostream( std::ostream & ost, std::ostream & oss, bool const tracer){
	if(tracer) return ost;
	return oss;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////LoopAnalyzerMover////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
LoopAnalyzerMover::LoopAnalyzerMover( protocols::loops::Loops const & loops, bool const tracer ) :
	Mover(),
	loops_(new protocols::loops::Loops(loops)),
	tracer_(tracer),
	sf_(NULL),
	chbreak_sf_(NULL),
	total_score_(0),
	max_rama_(-100000),
	max_chainbreak_(-100000),
	max_omega_(-100000),
	max_pbond_(-100000)
{
	protocols::moves::Mover::type( "LoopAnalyzer" );
	set_sf();
}

LoopAnalyzerMover::~LoopAnalyzerMover() {}

///@brief do not use a default constructor with this class - function exists as part of the remove #include drive
LoopAnalyzerMover::LoopAnalyzerMover() : tracer_(false) { utility_exit_with_message("do not use default constructor of class LoopAnalyzerMover"); }

LoopAnalyzerMover::LoopAnalyzerMover( LoopAnalyzerMover const & rhs ) :
	//utility::pointer::ReferenceCount(),
	Mover(),
	loops_(new protocols::loops::Loops(*(rhs.loops_))),
	tracer_(rhs.tracer_),
	positions_(rhs.positions_), //this is useless data
	sf_(rhs.sf_->clone()),
	chbreak_sf_(rhs.chbreak_sf_->clone()),
	scores_(rhs.scores_) //useless
{}

void LoopAnalyzerMover::set_sf(){
	using namespace core::scoring;

	sf_ = new ScoreFunction;
	sf_->set_weight( rama, 1.0 );
	sf_->set_weight( omega, 1.0 );
	sf_->set_weight( fa_dun, 1.0 );
	//sf_->set_weight( p_aa_pp, 1.0 );
	//sf_->set_weight( chainbreak, 1.0 );
	sf_->set_weight( peptide_bond, 1.0 );

	chbreak_sf_ = new ScoreFunction;
	chbreak_sf_->set_weight( chainbreak, 20.0 );

	return;
}

///@details LoopAnalyzerMover is mostly a container for other movers for the anchored design protocol.
void LoopAnalyzerMover::apply( core::pose::Pose & input_pose )
{
	TR << "running LoopAnalyzerMover" << std::endl;

	//prep pose
	core::pose::Pose pose(input_pose); //protecting input pose from our chainbreak changes (and its energies object)
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

	for( core::Size i(1); i <= positions_.size(); ++i ){
		core::Size const res(positions_[i]);
		//peptide_bond distance
		core::Real const pbnd_dist(input_pose.residue(res).atom("C").xyz().distance(input_pose.residue(res+1).atom("N").xyz()));
		//phi/psi/omega
		//		results.precision(3);
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


	if(!tracer_){
		protocols::jd2::JobDistributor::get_instance()->current_job()->add_string(results_oss.str());
		//store the loop_total where jd2 silent file can get it
		protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair("LAM_total", total_score_);
	}

	return;
}//LoopAnalyzerMover::apply

std::string
LoopAnalyzerMover::get_name() const {
	return "LoopAnalyzerMover";
}

core::Real
LoopAnalyzerMover::get_total_score() const{
	return total_score_;
}

core::Real
LoopAnalyzerMover::get_max_rama() const{
	return max_rama_;
}

core::Real
LoopAnalyzerMover::get_max_omega() const{
	return max_omega_;
}
core::Real
LoopAnalyzerMover::get_max_pbond() const{
	return max_pbond_;
}

core::Real
LoopAnalyzerMover::get_max_chainbreak() const{
	return max_chainbreak_;
}

utility::vector1<core::Real>
LoopAnalyzerMover::get_chainbreak_scores(){
	return scores_;
}

void LoopAnalyzerMover::find_positions( core::pose::Pose const & pose ){
	positions_.clear();
	for( protocols::loops::Loops::const_iterator it=loops_->begin(), it_end=loops_->end(); it != it_end; ++it ){
		core::Size start(it->start());
		core::Size end(it->stop());
		if(!pose.residue(start).is_terminus() && start != 1) --start;
		if(!pose.residue(end).is_terminus() && end != pose.total_residue()) ++end;
		for( core::Size i(start); i <= end; ++i ) positions_.push_back(i);
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
	for(core::Size i(1); i<=numpos; ++i){
		core::Size const pos = positions_[i];
		core::pose::remove_variant_type_from_pose_residue(pose, CUTPOINT_LOWER, pos);
		core::pose::remove_variant_type_from_pose_residue(pose, CUTPOINT_UPPER, pos);
	}

	//create a chainbreak at the desired position, score it, store the score, and revert changes
	for(core::Size i(1); i<=numpos; ++i){
		core::Size const pos = positions_[i];

		//create chainbreak variant
		core::pose::add_variant_type_to_pose_residue(pose, CUTPOINT_LOWER, pos);
		core::pose::add_variant_type_to_pose_residue(pose, CUTPOINT_UPPER, pos+1);

		//create cut in fold tree
		core::kinematics::FoldTree ft(pose.total_residue());
		ft.clear(); //remove initial edge
		ft.add_edge(Edge(1, pos, Edge::PEPTIDE));
		ft.add_edge(Edge(pos+1, pose.total_residue(), Edge::PEPTIDE));
		ft.add_edge(Edge(pos, pos+1, 1)); //a jump and therefore a cut
		ft.reorder(1);
		pose.fold_tree(ft);

		//score the sucker
		scores_[i] = (*chbreak_sf_)(pose);

		//undo that chainbreak variant
		core::pose::remove_variant_type_from_pose_residue(pose, CUTPOINT_LOWER, pos);
		core::pose::remove_variant_type_from_pose_residue(pose, CUTPOINT_UPPER, pos+1);
	}

	//revert to innocuous fold tree
	pose.fold_tree(core::kinematics::FoldTree(pose.total_residue()));

	return;
}//calculate_all_chainbreaks


}//analysis
}//protocols
