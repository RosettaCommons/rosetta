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
/// @author Liz Kellogg

#include <protocols/moves/Mover.hh>
#include <protocols/ddg/ddGMover.hh>

#include <core/types.hh>


#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <protocols/scoring/Interface.hh> //added ek for interface ddgs

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

//constraint stuff
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/Func.hh>

#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>

#include <basic/options/util.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/ddg.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>


#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <protocols/rigid/RigidBodyMover.hh>

#include <fstream>
#include <iostream>
#include <sstream>

//new includes for rotamer constraints

// C++ headers
#include <cstdlib>
#include <string>
#include <sys/stat.h>

#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "protocols.moves.ddGMover" );

namespace protocols {
namespace ddg {


using namespace core;
using namespace core::scoring;

typedef std::vector<double> ddGs;

ddGMover::ddGMover() :
	Mover("ddGMover"),
	num_decoys_used_in_calculations_(20), //default to number of iterations
	mean_(false),
	min_(true),
	min_cst_(true),
	nbr_cutoff_(8.0),
	restrict_to_nbrhood_(false),
	interface_ddg_(false),
	scorefxn_(/* 0 */),
	min_cst_sfxn_(/* 0 */),
	min_cst_sfxn_no_cst_weight_( /* 0 */ ),
	cst_set_(core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet() )),
	min_cst_set_wt_(core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet() )),
	min_cst_wt_types_(),
	min_cst_set_mut_(core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet() )),
	min_cst_mut_types_(),
	repack_cst_set_wt_(core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet() )),
	repack_wt_types_(),
	repack_cst_set_mut_(core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet() )),
	repack_mut_types_(),
	residues_to_mutate_(0,core::chemical::aa_unk),
	num_iterations_(20),
	dmp_pdb_(false),
	dbg_output_(false),
	wt_components_(1,1,-999.99),
	wt_unbound_components_(1,1,-999.99),//only used for interface ddg
	mutant_components_(1,1,-999.99),
	mutant_unbound_components_(1,1,-999.99),//only used for interface ddg
	natives_(),
	mutants_()
{}

ddGMover::ddGMover(
	core::scoring::ScoreFunctionOP s,
	core::scoring::ScoreFunctionOP m,
	utility::vector1<core::chemical::AA> res_to_mutate
) :
	Mover("ddGMover"),
	num_decoys_used_in_calculations_(20),
	mean_(false),
	min_(true),
	min_cst_(true),
	nbr_cutoff_(8.0),
	restrict_to_nbrhood_(false),
	interface_ddg_(false),
	scorefxn_(s),
	min_cst_sfxn_( m->clone() ),
	min_cst_sfxn_no_cst_weight_( m->clone() ),
	cst_set_(core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet() )),
	min_cst_set_wt_(core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet() )),
	min_cst_wt_types_(),
	min_cst_set_mut_(core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet() )),
	min_cst_mut_types_(),
	repack_cst_set_wt_(core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet() )),
	repack_wt_types_(),
	repack_cst_set_mut_(core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet() )),
	repack_mut_types_(), //these arrays only because
	//theres no easy way to get constraint info after the fact
	residues_to_mutate_(res_to_mutate),
	num_iterations_(20),
	dmp_pdb_(false),
	dbg_output_(false),
	wt_components_(1,1,-999.99),
	wt_unbound_components_(1,1,-999.99),//only used for interface ddg
	mutant_components_(1,1,-999.99),
	mutant_unbound_components_(1,1,-999.99),//only used for interface ddg
	natives_(),
	mutants_()
{}

ddGMover::~ddGMover(){}


utility::vector1<int>
ddGMover::find_nbrs(
	pose::Pose & p,
	utility::vector1<int> & mutation_position,
	double radii
){
	utility::vector1<int> nbrs;
	for ( core::Size i = 1; i <= mutation_position.size(); ++i ) {
		numeric::xyzVector<double> i_pos;
		if ( p.residue(mutation_position[i]).name1() == 'G' ) {
			i_pos = p.residue(mutation_position[i]).xyz(" CA ");
		} else {
			i_pos = p.residue(mutation_position[i]).xyz(" CB ");
		}
		for ( core::Size j = 1; j <= p.total_residue(); j++ ) {
			numeric::xyzVector<double> j_pos;
			if ( p.residue(j).name1() == 'G' ) {
				j_pos = p.residue(j).xyz(" CA ");
			} else {
				j_pos = p.residue(j).xyz(" CB ");
			}
			if ( i_pos.distance(j_pos) <= radii ) {
				nbrs.push_back(j);
			}
		}
	}
	return nbrs;
}

bool
ddGMover::is_complete(ObjexxFCL::FArray2D<double> to_check){
	bool is_complete=true;
	for ( int i = to_check.l1(); i <= to_check.u1(); i++ ) {
		for ( int j = to_check.l2(); j <= to_check.u2(); j++ ) {
			if ( to_check(i,j) == -999.99 ) {
				is_complete=false;
			}
		}
	}
	return is_complete;
}

double
ddGMover::sum(ddGs &scores_to_sum)
{
	double sum=0;
	for ( unsigned int i =0; i<scores_to_sum.size(); i++ ) {
		sum+=scores_to_sum[i];
	}
	return sum;
}

double
ddGMover::average( utility::vector1<double> &scores_to_average)
{
	double sum = 0;
	for ( unsigned int i =1; i<=scores_to_average.size(); i++ ) {
		sum+=scores_to_average[i];
	}
	return (sum/scores_to_average.size());
}

int
ddGMover::store_energies(
	ObjexxFCL::FArray2D< double > & two_d_e_arrays,
	core::scoring::ScoreFunction & s,
	pose::Pose &p,
	int next_index,
	int size_to_expect
)
{
	p.update_residue_neighbors();
	s(p);

	//hack. we don't want to consider atom-pair constraints in computing ddgs!
	using core::scoring::score_type_from_name;
	core::scoring::ScoreType const atom_pair_cst = score_type_from_name("atom_pair_constraint");
	core::Real saved_weight = 0;
	if ( s.get_weight( atom_pair_cst) > 0 ) {
		saved_weight = s.get_weight(atom_pair_cst);
	}
	s.set_weight( atom_pair_cst, (core::Real)0);


	//all this to determine how many non-zero weights there are
	int num_score_components = 0;
	using core::scoring::EnergyMap;
	using core::scoring::ScoreType;
	for ( EnergyMap::const_iterator it = s.weights().begin(); it != s.weights().end(); ++it ) {
		if ( *it != 0.0 ) {
			num_score_components++;
		}
	}
	//
	if ( (two_d_e_arrays.u1() != num_score_components) ||
			(two_d_e_arrays.u2() != size_to_expect) ) {
		two_d_e_arrays.dimension(num_score_components,size_to_expect,-999.99);
	}
	int current_score_component=0;
	int j =1;
	for ( EnergyMap::const_iterator i = (s.weights()).begin(), end = s.weights().end(); i != end; ++i ) {
		//get score component of pose, then store in next slot of two_d_e_arrays
		current_score_component++;
		if ( *i != 0 ) {
			two_d_e_arrays(j++,next_index)=(*i) *
				((p.energies()).total_energies())[ScoreType(current_score_component)];
		}
	}

	s.set_weight(atom_pair_cst,saved_weight);

	return 0;
}

int
ddGMover::average_score_components(
	ObjexxFCL::FArray2D< double > &scores_to_average,
	utility::vector1<double> &averaged_scores
)
{
	if ( num_decoys_used_in_calculations_ > num_iterations_ ) {
		num_decoys_used_in_calculations_ = num_iterations_;
	} //if user specifies more decoys to use than produced, simply set to the number
	//of decoys produced.

	if ( num_decoys_used_in_calculations_ == num_iterations_ ) { //if we do a straight average, things are easy
		averaged_scores = utility::vector1<double>(scores_to_average.u1()-scores_to_average.l1()+1);
		for ( int i = scores_to_average.l1(); i <= scores_to_average.u1(); i++ ) {
			double sum_score_component = 0;
			for ( int j = scores_to_average.l2(); j <= scores_to_average.u2(); j++ ) {
				sum_score_component += scores_to_average(i,j);
			}
			averaged_scores[i]=
				sum_score_component/(scores_to_average.u2()-scores_to_average.l2()+1);
		}
		return 0;
	} else {
		typedef utility::vector1<double> score_components;
		typedef std::pair<double,score_components> enrgs;

		utility::vector1<enrgs> scores_to_sort;

		for ( int i = scores_to_average.l2(); i <= scores_to_average.u2(); i++ ) {
			double total_score = 0; utility::vector1<double> sc;
			enrgs e;
			for ( int j = scores_to_average.l1(); j <= scores_to_average.u1(); j++ ) {
				total_score += scores_to_average(j,i);
				sc.push_back(scores_to_average(j,i));
			}
			e.first = total_score;
			e.second = sc;
			scores_to_sort.push_back(e);

		}
		utility::vector1<enrgs> sorted;

		while ( sorted.size() < num_decoys_used_in_calculations_ ) {
			double min = 1000; int min_index = -1;
			for ( unsigned int i=1; i<= scores_to_sort.size(); i++ ) {
				if ( scores_to_sort[i].first < min ) {
					min = scores_to_sort[i].first;
					min_index = i;
				}
			}
			if ( min_index != -1 ) {
				sorted.push_back(scores_to_sort[min_index]);
				scores_to_sort.erase(scores_to_sort.begin()+(min_index-1));
			}
		}

		for ( unsigned int i =1; i <= sorted[1].second.size(); i++ ) {
			double sum_score_component = 0.0;
			for ( unsigned int j = 1; j <= sorted.size(); j++ ) {
				sum_score_component += sorted[j].second[i];
			}
			averaged_scores.push_back(sum_score_component/(sorted.size()));
		}

		return 0;
	}
}

void
ddGMover::calculate_interface_unbound_energy(
	core::pose::Pose & p,
	core::scoring::ScoreFunctionOP s,
	core::pack::task::PackerTaskOP pt
)
{

	using namespace protocols::moves;
	using namespace core::pack;
	//from sarel. thank you!
	int const rb_jump(1);
	rigid::RigidBodyTransMoverOP separate_partners( new rigid::RigidBodyTransMover(p,rb_jump) );
	separate_partners->step_size(1000.0);
	separate_partners->apply(p);
	core::pack::pack_rotamers(p,(*s),pt);
}

void
ddGMover::setup_repack_constraints(
	pose::Pose & pose,
	ScoreFunctionOP sfxn,
	bool /*all_but_mut_site*/,
	utility::vector1<char> & /*constraints_at_pos*/
)
{
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace basic::options::OptionKeys;
	using namespace basic::options;

	double const CONSTRAINT_WEIGHT = 1.0;

	pose.remove_constraints();
	//if member object constraint set is non-zero, then
	//transfer constraints and be done with it!
	utility::vector1<int> mutations;
	for ( unsigned int i=1; i <= residues_to_mutate_.size(); i++ ) {
		if ( residues_to_mutate_[i] != core::chemical::aa_unk ) {
			mutations.push_back(i);
		}
	} //i is mutation position

	sfxn->set_weight(atom_pair_constraint, CONSTRAINT_WEIGHT);

}

void
ddGMover::setup_constraints(
	pose::Pose & pose,
	ScoreFunctionOP sfxn,
	float const cst_tol,
	bool /*all_but_mut_site*/,
	utility::vector1<char> & /*constraints_at_pos*/
)
{
	//create constraints for all residues
	//type: HARMONIC
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace basic::options::OptionKeys;
	using namespace basic::options;
	static float const CA_cutoff(9.0);
	double const CONSTRAINT_WEIGHT = 1.0;

	pose.remove_constraints();

	//if member object constraint set is non-zero, then
	//transfer constraints and be done with it!
	//if(cst_set_->has_constraints()){ //depreciated. replaced by specific constraintsets for each phase of optimization
	// TR << "constraint set already exists, using member object" << std::endl;
	// pose.constraint_set(cst_set_);
	//}//if we defined constraints outside of the program, use those instead
	if ( basic::options::option[basic::options::OptionKeys::constraints::cst_file].user() ) {
		core::scoring::constraints::ConstraintSetOP cstset( new
			core::scoring::constraints::ConstraintSet() );
		cstset = core::scoring::constraints::ConstraintIO::read_constraints(
			option[basic::options::OptionKeys::constraints::cst_file][1],cstset,
			pose);
		pose.constraint_set(cstset);
	} else { //otherwise, define your own constraints
		int nres = pose.total_residue();
		for ( int i = 1; i <= nres; i++ ) {
			if ( pose.residue(i).is_protein() ) {
				Vector const CA_i( pose.residue(i).xyz(" CA "));
				for ( int j = 1; j < i; j++ ) {
					if ( pose.residue(j).is_protein() ) {
						Vector const CA_j(pose.residue(j).xyz(" CA "));
						Real const CA_dist = (CA_i - CA_j).length();
						if ( CA_dist < CA_cutoff ) {
							core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc(CA_dist, cst_tol) );
							ConstraintCOP cst( ConstraintOP( new AtomPairConstraint(
								core::id::AtomID(pose.residue(i).atom_index(" CA "),i),
								core::id::AtomID(pose.residue(j).atom_index(" CA "),j),
								fx) ) );
							pose.add_constraint(cst);
						}
					}
				}
			}
		}
	}
	sfxn->set_weight(atom_pair_constraint, CONSTRAINT_WEIGHT);
	TR.Debug<< "size of cst set is " << pose.constraint_set()->get_all_constraints().size() << std::endl;
}

void
ddGMover::minimize_with_constraints(
	pose::Pose & pose,
	ScoreFunctionOP s,
	core::pack::task::PackerTaskOP pt
)
{
	using namespace basic::options::OptionKeys;
	using namespace basic::options;

	core::optimization::AtomTreeMinimizer min_struc;
	float minimizer_tol = 0.000001;
	core::optimization::MinimizerOptions options(
		"dfpmin_armijo_nonmonotone", minimizer_tol,
		true /*use_nb_list*/,
		false /*deriv_check_in*/,
		true /*deriv_check_verbose_in*/);

	options.nblist_auto_update( true );
	options.max_iter(5000); //otherwise, they don't seem to converge
	core::kinematics::MoveMap mm;
	mm.set_bb(true);
	mm.set_chi(true);

	if ( basic::options::option[OptionKeys::ddg::sc_min_only]() ) {
		mm.set_bb(false);
		//minimizer_tol = 0.0001;
	}


	if ( basic::options::option[OptionKeys::ddg::ramp_repulsive]() ) {
		//set scorefxn fa_rep to 1/10 of original weight and then minimize
		ScoreFunctionOP one_tenth_orig(s->clone());
		//reduce_fa_rep(0.1,*one_tenth_orig);
		//min_struc.run(p,mm,s,options);

		one_tenth_orig->set_weight( core::scoring::fa_rep, s->get_weight(core::scoring::fa_rep)*0.1 );
		min_struc.run(pose,mm,*one_tenth_orig,options);

		//then set scorefxn fa_rep to 1/3 of original weight and then minimize
		ScoreFunctionOP one_third_orig(s->clone());
		//reduce_fa_rep(0.33,*one_third_orig);
		one_third_orig->set_weight( core::scoring::fa_rep, s->get_weight(core::scoring::fa_rep)*0.3333);

		min_struc.run(pose,mm,*one_third_orig,options);
		//then set scorefxn fa_rep to original weight and then minimize
	}

	(*s)(pose);
	s->show(std::cout,pose);
	pose::Pose before(pose);
	min_struc.run(pose,mm,*s,options);
	//(*s)(pose);
	while ( std::abs((*s)(pose)-(*s)(before)) > 1 ) {
		//std::cout << "running another iteration of minimization. difference is: " << ((*s)(pose)-(*s)(before)) << std::endl;
		before=pose;

		if ( basic::options::option[OptionKeys::ddg::pack_until_converge]() ) {
			core::pack::pack_rotamers(pose,(*(this->score_function())),pt);
			//use packing score function, not minimization score function
		}

		if ( basic::options::option[OptionKeys::ddg::ramp_repulsive]() ) {
			//set scorefxn fa_rep to 1/10 of original weight and then minimize
			ScoreFunctionOP one_tenth_orig(s->clone());
			//reduce_fa_rep(0.1,*one_tenth_orig);
			//min_struc.run(p,mm,s,options);

			one_tenth_orig->set_weight(core::scoring::fa_rep, s->get_weight(core::scoring::fa_rep)*0.1);
			min_struc.run(pose,mm,*one_tenth_orig,options);

			//then set scorefxn fa_rep to 1/3 of original weight and then minimize
			ScoreFunctionOP one_third_orig(s->clone());
			//reduce_fa_rep(0.33,*one_third_orig);
			one_third_orig->set_weight(core::scoring::fa_rep, s->get_weight(core::scoring::fa_rep)*0.3333);
			min_struc.run(pose,mm,*one_third_orig,options);
		}
		//then set scorefxn fa_rep to original weight and then minimize
		min_struc.run(pose,mm,*s,options);

	}


	s->show(std::cout, pose);
}

bool
sort_numerically_ascending(double a, double b){
	return a < b;
}

/// @details APL Note that utility::arg_min does exactly this without the awful 99999999 bug
int
min_index(utility::vector1<double>scores)
{
	double min=9999999;
	int min_index=-1;
	for ( unsigned int i=1; i<=scores.size(); i++ ) {
		if ( scores[i] < min ) {
			min_index=i;
			min = scores[i];
		}
	}
	return min_index;
}

void
ddGMover::neighbor_cutoff(double cutoff){
	nbr_cutoff_ = cutoff;
}

void
ddGMover::restrict_to_nbrs(bool truefalse){
	restrict_to_nbrhood_ = truefalse;
}

void
ddGMover::set_minimization_score_function( core::scoring::ScoreFunctionOP s){
	min_cst_sfxn_ = s->clone();
	min_cst_sfxn_no_cst_weight_ = s->clone();
}

void
ddGMover::score_function( core::scoring::ScoreFunctionOP s){
	scorefxn_=s;
}

void
ddGMover::num_iterations( int num ){
	num_iterations_=num;
}

void
ddGMover::dump_pdbs(bool truefalse){
	dmp_pdb_=truefalse;
}

void
ddGMover::debug_output(bool truefalse){
	dbg_output_=truefalse;
}

void
ddGMover::is_interface_ddg(bool truefalse){
	interface_ddg_=truefalse;
}

void
ddGMover::wt_score_components(ObjexxFCL::FArray2D<double> wsc){
	wt_components_=wsc;
}

void
ddGMover::wt_unbound_score_components(ObjexxFCL::FArray2D<double> wusc){
	wt_unbound_components_=wusc;
}

void
ddGMover::mutant_score_components(ObjexxFCL::FArray2D<double> msc){
	mutant_components_=msc;
}

void
ddGMover::residues_to_mutate(utility::vector1<core::chemical::AA> residues){
	residues_to_mutate_=residues;
}

void
ddGMover::set_min_cst(bool truefalse){
	min_cst_ = truefalse;
}

void
ddGMover::set_mean(bool truefalse){
	mean_ = truefalse;
}

void
ddGMover::set_min(bool truefalse){
	min_  = truefalse;
	if ( min_ ) {
		num_decoys_used_in_calculations_ = 1;
	}
}

void
ddGMover::set_num_decoys_used_in_calculations(core::Real num_lowe_used){
	num_decoys_used_in_calculations_ = num_lowe_used;
}

Real
ddGMover::neighbor_cutoff(){
	return nbr_cutoff_;
}

bool
ddGMover::restrict_to_nbrs(){
	return restrict_to_nbrhood_;
}

core::scoring::ScoreFunctionOP
ddGMover::score_function(){
	return scorefxn_;
}

core::scoring::ScoreFunctionOP
ddGMover::minimization_score_function(){
	return min_cst_sfxn_;
}

int
ddGMover::num_iterations(){return num_iterations_;}

bool
ddGMover::is_interface_ddg(){return interface_ddg_;}

ObjexxFCL::FArray2D<double>
ddGMover::wt_score_components(){
	return wt_components_;
}

ObjexxFCL::FArray2D<double>
ddGMover::wt_unbound_score_components(){
	return wt_unbound_components_;
}

ObjexxFCL::FArray2D<double>
ddGMover::mutant_score_components(){
	return mutant_components_;
}

utility::vector1<core::chemical::AA>
ddGMover::residues_to_mutate(){
	return residues_to_mutate_;
}

utility::vector1<double>
ddGMover::get_wt_min_score_components(){
	//std::cout << "computing minimum energy components" << std::endl;
	utility::vector1<double> wt_min_score_components;
	utility::vector1<double> total_scores;
	//sum all components for each structure
	for ( int i = wt_components_.l2(); i <= wt_components_.u2(); i++ ) {
		double total_score = 0.0;
		for ( int j = wt_components_.l1(); j <= wt_components_.u1(); j++ ) {
			total_score += wt_components_(j,i);
		}
		//std::cout << "total score computed for wt number " << i << " is " << total_score << std::endl;
		total_scores.push_back(total_score);
	}
	//num_score_components, size_to_expect
	//find minimum index
	int min_index = -1;
	double min_score = 10000;
	for ( unsigned int i = 1; i <= total_scores.size(); i++ ) {
		if ( total_scores[i] < min_score ) {
			min_index = i;
			min_score = total_scores[i];
		}
	}
	assert(min_index != -1);
	//std::cout << "min_index is " << min_index << " and corresponds to an energy " << min_score  << std::endl;
	//return score components of min energy structure
	for ( int i = wt_components_.l1(); i <= wt_components_.u1(); i++ ) {
		wt_min_score_components.push_back(wt_components_(i,min_index));
	}
	return wt_min_score_components;
}

utility::vector1<double>
ddGMover::get_wt_averaged_score_components(){
	//std::cout << "computing wt averaged score components" << std::endl;
	utility::vector1<double> wt_averaged_score_components;
	if ( !interface_ddg_ ) {
		average_score_components(wt_components_,wt_averaged_score_components);
	} else {
		assert((wt_components_.u1()-wt_components_.l1())== (wt_unbound_components_.u1()-wt_unbound_components_.l1()) &&
			(wt_components_.u2()-wt_components_.l2()) == (wt_unbound_components_.u2()-wt_unbound_components_.l2()));
		ObjexxFCL::FArray2D<double> dG_bound_unbound((wt_components_.u1()-wt_components_.l1()+1),
			(wt_components_.u2()-wt_components_.l2()+1),-999.99);
		for ( int i =wt_components_.l1(); i<=wt_components_.u1(); i++ ) {
			for ( int j = wt_components_.l2(); j<= wt_components_.u2(); j++ ) {
				dG_bound_unbound(i,j)=wt_components_(i,j)- wt_unbound_components_(i,j);
			}
		}
		average_score_components(dG_bound_unbound,wt_averaged_score_components);
	}
	return wt_averaged_score_components;
}

utility::vector1<double>
ddGMover::get_mutant_min_score_components(){
	//std::cout << "computing minimum energy score components for mutant" << std::endl;
	utility::vector1<double> mutant_min_score_components;
	utility::vector1<double> total_scores;
	//sum all components for each structure
	for ( int i = mutant_components_.l2(); i <= mutant_components_.u2(); i++ ) {
		double total_score = 0.0;
		for ( int j = mutant_components_.l1(); j <= mutant_components_.u1(); j++ ) {
			total_score += mutant_components_(j,i);
		}
		//std::cout << "total score computed for mut number " << i << " is " << total_score << std::endl;
		total_scores.push_back(total_score);
	}
	//num_score_components, size_to_expect
	//find minimum index
	int min_index = -1;
	double min_score = 10000;
	for ( unsigned int i = 1; i <= total_scores.size(); i++ ) {
		if ( total_scores[i] < min_score ) {
			min_index = i;
			min_score = total_scores[i];
		}
	}
	assert(min_index != -1);
	//std::cout << "MUT min_index is " << min_index << " and corresponds to an energy " << min_score  << std::endl;
	//return score components of min energy structure
	for ( int i = mutant_components_.l1(); i <= mutant_components_.u1(); i++ ) {
		mutant_min_score_components.push_back(mutant_components_(i,min_index));
	}
	return mutant_min_score_components;
}

utility::vector1<double>
ddGMover::get_mutant_averaged_score_components(){
	//std::cout << "computing average mutant score components" << std::endl;
	utility::vector1<double> mutant_averaged_score_components;
	if ( !interface_ddg_ ) {
		average_score_components(mutant_components_, mutant_averaged_score_components);
	} else {
		assert(
			(mutant_components_.u1()-mutant_components_.l1()) == (mutant_unbound_components_.u1()-mutant_unbound_components_.l1()) &&
			(mutant_components_.u2()-mutant_components_.l2()) == (mutant_unbound_components_.u2()-mutant_unbound_components_.l2()));
		ObjexxFCL::FArray2D<double> dG_bound_unbound(
			(mutant_components_.u1()-mutant_components_.l1()+1),
			(mutant_components_.u2()-mutant_components_.l2()+1),
			-999.99);
		for ( int i = mutant_components_.l1(); i <= mutant_components_.u1(); ++i ) {
			for ( int j = mutant_components_.l2(); j <= mutant_components_.u2(); ++j ) {
				dG_bound_unbound(i,j) = mutant_components_(i,j) - mutant_unbound_components_(i,j);
			}
		}
		average_score_components(dG_bound_unbound,mutant_averaged_score_components);
	}
	return mutant_averaged_score_components;
}

utility::vector1<double>
ddGMover::get_delta_energy_components()
{
	utility::vector1<double> wt;
	utility::vector1<double> mut;
	if ( mean_ ) {
		wt = this->get_wt_averaged_score_components();
		mut = this->get_mutant_averaged_score_components();
	} else if ( min_ ) {
		wt = this->get_wt_min_score_components();
		mut = this->get_mutant_min_score_components();
	}
	utility::vector1<double> delta_energy;
	assert(wt.size() == mut.size());
	for ( unsigned int i=1; i<=wt.size(); i++ ) {
		delta_energy.push_back(mut[i]-wt[i]);
	}
	return delta_energy;
}

void
ddGMover::get_scorefunction_header(
	core::scoring::ScoreFunctionOP sfxn,
	utility::vector1<std::string> & components
)
{
	//utility::vector1<std::string> components;
	core::scoring::ScoreType const atom_pair_cst = score_type_from_name("atom_pair_constraint"); //atom-pair-csts not considered in ddg score
	if ( sfxn != 0 ) {
		int score_component=0;
		for ( EnergyMap::const_iterator i = (sfxn->weights()).begin(), end = sfxn->weights().end(); i != end; ++i ) {
			score_component++;
			if ( *i != 0 && ScoreType(score_component) != atom_pair_cst ) {
				components.push_back(name_from_score_type(ScoreType(score_component)));
			}
		}
	}
	// return components;
}

std::string
ddGMover::mutation_label( pose::Pose const & pose ) const
{
	std::string mutation_label="";
	if ( residues_to_mutate_.size() != pose.total_residue() ) {
		return mutation_label; //residues_to_mutate_ hasn't been initialized
	} else {
		for ( unsigned int i=1; i <= residues_to_mutate_.size(); i++ ) {
			if ( residues_to_mutate_[i] != core::chemical::aa_unk ) {
				std::ostringstream q;
				q << i;
				mutation_label = mutation_label + (pose.residue(i)).name1() + q.str() + core::chemical::oneletter_code_from_aa(residues_to_mutate_[i]);
				//std::cout << "at position " << i << " the mutation label now becomes " << mutation_label << std::endl; //DEBUG
			}
		}
	}
	return mutation_label;
}

double
ddGMover::get_wt_averaged_totals(){
	utility::vector1<double> wt = this->get_wt_averaged_score_components();
	double sum=0;
	for ( unsigned int i =1; i<=wt.size(); i++ ) {
		sum=sum+wt[i];
	}
	return sum;
}

double
ddGMover::get_wt_min_totals(){
	utility::vector1<double> wt = this->get_wt_min_score_components();
	double sum = 0;
	for ( unsigned int i =1; i <= wt.size(); i++ ) {
		sum = sum + wt[i];
	}
	return sum;
}


double
ddGMover::get_mutant_averaged_totals(){
	utility::vector1<double> mut = this->get_mutant_averaged_score_components();
	double sum=0;
	for ( unsigned int i =1; i<=mut.size(); i++ ) {
		sum=sum+mut[i];
	}
	return sum;
}

double
ddGMover::get_mutant_min_totals(){
	utility::vector1<double> mut = this->get_mutant_min_score_components();
	double sum = 0;
	for ( unsigned int i =1; i <= mut.size(); i++ ) {
		sum = sum + mut[i];
	}
	return sum;
}

double
ddGMover::ddG(){
	if ( mean_ ) {
		return get_mutant_averaged_totals() - get_wt_averaged_totals();
	} else if ( min_ ) {
		return get_mutant_min_totals() - get_wt_min_totals();
	} else {
		std::cerr << "neither mean or min set, reverting to default" << std::endl;
		min_=true; mean_=false;
		return get_mutant_min_totals() - get_wt_min_totals();
	}
}

bool
ddGMover::is_wt_calc_complete(){
	return is_complete(wt_components_);
}

bool
ddGMover::is_mutant_calc_complete(){
	return is_complete(mutant_components_);
}

bool
ddGMover::is_properly_initialized(
	pose::Pose & pose
)
{
	bool is_initialized=true;
	//check scorefxn
	//check residues_to_mutate_
	if ( residues_to_mutate_.size() != pose.total_residue() ) {
		is_initialized=false;
	}
	//check num_iterations_
	//check dmp_pdbs
	//check dbg_output
	return is_initialized;
}

bool
ddGMover::get_min_cst(){
	return min_cst_;
}

bool
ddGMover::use_mean(){
	return mean_;
}

bool
ddGMover::use_min(){
	return min_;
}

core::Real
ddGMover::get_num_decoys_used_in_calculations(){
	return num_decoys_used_in_calculations_;
}

void
ddGMover::relax_wildtype_structure(
	core::pose::Pose & pose,
	protocols::scoring::Interface & protein_interface,
	std::string const & wt_traj,
	bool output_silent
)
{
	using namespace pose;
	using namespace scoring;
	using namespace conformation;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::task;


	//debug output?
	if ( dbg_output_ ) {
		TR.Debug << "weights being used: " <<
			scorefxn_->weights() << "\n";
	}
	min_cst_wt_types_.resize(pose.total_residue(),'X');
	min_cst_mut_types_.resize(pose.total_residue(),'X');
	repack_wt_types_.resize(pose.total_residue(),'X');
	repack_mut_types_.resize(pose.total_residue(),'X');

	//store constraints for later usage
	if ( basic::options::option[OptionKeys::ddg::use_rotamer_constraints_to_native]() ) {
		setup_repack_constraints(pose,scorefxn_,false,repack_wt_types_); //not going to use now, but later
		repack_cst_set_wt_ = core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet(*(pose.constraint_set())) );
		//
		setup_repack_constraints(pose,scorefxn_,true,repack_mut_types_);
		repack_cst_set_mut_ = core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet(*(pose.constraint_set())) );
	}

	if ( basic::options::option[OptionKeys::ddg::min_cst]() ) {
		if ( basic::options::option[OptionKeys::ddg::harmonic_ca_tether].user() ) {
			setup_constraints( pose, min_cst_sfxn_, basic::options::option[OptionKeys::ddg::harmonic_ca_tether](), false, min_cst_wt_types_);
			min_cst_set_wt_ = core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet(*(pose.constraint_set())) );
			//
			setup_constraints( pose, min_cst_sfxn_, basic::options::option[OptionKeys::ddg::harmonic_ca_tether](), true, min_cst_mut_types_);
			min_cst_set_mut_ = core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet(*(pose.constraint_set())) );

		} else {
			setup_constraints(pose,min_cst_sfxn_,0.5,false,min_cst_wt_types_);
			min_cst_set_wt_ = core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet(*(pose.constraint_set())) );
			//
			setup_constraints(pose,min_cst_sfxn_,0.5,true,min_cst_mut_types_);
			min_cst_set_mut_ = core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet(*(pose.constraint_set())) );
		}
	}
	pose.remove_constraints(); //only reason we added constraints was to create a storage set of constraints
	//for later use

	if ( option[OptionKeys::ddg::opt_input_structure].user() && option[OptionKeys::ddg::opt_input_structure]() ) {
		PackerTaskOP norepacking(pack::task::TaskFactory::create_packer_task(pose));
		norepacking->temporarily_fix_everything();
		pose.constraint_set(min_cst_set_wt_);
		minimize_with_constraints(pose,min_cst_sfxn_,norepacking);
		pose.dump_pdb("input_pose.pdb");
	}

	// initialize the scoring function stuff
	(*scorefxn_)(pose);
	/// Now handled automatically.  scorefxn_->accumulate_residue_total_energies(pose);

	if ( interface_ddg_ ) {
		protein_interface.distance(10.0);
		protein_interface.calculate(pose);
	}

	bool dG_native_calculated = this->is_wt_calc_complete();
	if ( !dG_native_calculated &&
			!(basic::options::option[OptionKeys::ddg::mut_only].user() && (basic::options::option[OptionKeys::ddg::mut_only]())) ) {

		if ( dbg_output_ ) {
			TR << "dG_wildtype hasn't been calculated yet! " <<
				"starting to calculate dG for native structure\n";
		}

		//start filehandler for recording the native structure logfile
		std::ofstream record_trajectories;
		if ( wt_traj != "" ) {
			record_trajectories.open(wt_traj.c_str());
		} else {
			TR << "wt_traj is emtpy, not recording trajectories.\n";
		}

		core::io::silent::SilentFileData sfd;
		std::string silentfilename;
		if ( basic::options::option[out::file::silent].user() ) {
			silentfilename = basic::options::option[out::file::silent]();
			sfd.set_filename(silentfilename);
		} else {
			silentfilename = "wt_"+this->mutation_label(pose)+".out";
			sfd.set_filename(silentfilename); //setup in case we dump in form of silent file
		}

		//end add in interface specific stuff

		//PACKERTASK SETUP FOR NATIVE START
		//measure dG of input structure by repacking n times and taking the average score
		//this corresponds to wt bound energies of protein complex
		pack::task::PackerTaskOP repack_native(pack::task::TaskFactory::create_packer_task(pose));
		//add option here to restrict to local region
		//find mutation position
		utility::vector1<int> mutations;
		for ( unsigned int i=1; i <= residues_to_mutate_.size(); i++ ) {
			if ( residues_to_mutate_[i] != core::chemical::aa_unk ) {
				mutations.push_back(i);
			}
		} //i is mutation position
		//find all neighbors within 8 angstrom
		if ( restrict_to_nbrhood_ ) {
			if ( mutations.size() == 0 ) {
				TR << "FATAL ERROR no mutations specified, which are needed to determine neighbors" << std::endl;
				exit(1);
			}

			utility::vector1< bool > neighborhood = neighborhood_of_mutations( pose, mutations );
			// Now: report which residues will be repacked
			for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				if ( neighborhood[ ii ] ) {
					TR << "residue " << ii << " is within range of a mutating residue and will be repacked" << std::endl;
					repack_native->nonconst_residue_task(ii).restrict_to_repacking();
				} else {
					TR << "residue " << ii << " is not within range of any mutating residue and will NOT be repacked" << std::endl;
					repack_native->nonconst_residue_task(ii).prevent_repacking();
				}
			}

		} else {
			repack_native->restrict_to_repacking();
		}

		for ( unsigned int j =1; j <= pose.total_residue(); j++ ) {
			initialize_rotamer_behavior_for_residue_level_task( repack_native->nonconst_residue_task(j) );
			if ( interface_ddg_ ) {
				if ( ! protein_interface.is_interface( j ) ) {
					repack_native->nonconst_residue_task(j).prevent_repacking();
				}
			}
		}
		initialize_task_level_behavior( repack_native );

		//PACKERTASK SETUP FOR NATIVE END
		utility::vector1< std::pair < Real,std::string > > results;
		utility::vector1<core::pose::PoseOP> poses;
		pose::Pose temporary_pose = pose;

		temporary_pose.constraint_set(repack_cst_set_wt_);
		pack::pack_rotamers_loop(temporary_pose,(*scorefxn_),repack_native,num_iterations_,results,poses);

		for ( unsigned int i=1; i<=poses.size(); i++ ) {

			std::ostringstream q;
			q << i;

			pose::Pose resulting_pose=(*poses[i]);
			resulting_pose.remove_constraints();
			if ( min_cst_ ) {
				resulting_pose.constraint_set(min_cst_set_wt_);

				minimize_with_constraints(resulting_pose,min_cst_sfxn_,repack_native);
				( *min_cst_sfxn_no_cst_weight_) ( resulting_pose ); // now, report the score in the absence of constraints
				store_energies(wt_components_, (*min_cst_sfxn_no_cst_weight_), resulting_pose,i,num_iterations_);


			} else {
				resulting_pose.constraint_set(repack_cst_set_wt_);
			}

			TR << i <<
				" score before mutation: residue ";
			if ( !min_cst_ ) {
				TR << (*scorefxn_)(resulting_pose) << " " << resulting_pose.energies().total_energies().weighted_string_of( scorefxn_->weights() ) << std::endl;
			} else {
				TR << (*min_cst_sfxn_no_cst_weight_)(resulting_pose) << " " << resulting_pose.energies().total_energies().weighted_string_of( min_cst_sfxn_no_cst_weight_->weights()) << std::endl;
			}

			if ( !interface_ddg_ && !min_cst_ ) {
				store_energies(wt_components_, (*scorefxn_), resulting_pose,i,num_iterations_);
				//end store components
			} else if ( interface_ddg_ ) {
				assert(resulting_pose.chain(resulting_pose.total_residue())-pose.chain(1) !=0);
				//dGinterface = Gbound-Gunbound
				//double debug_bound_score = (*scorefxn_)(temporary_pose);
				store_energies(wt_components_, (*scorefxn_), resulting_pose,i,num_iterations_);
				pose::Pose temporary_unbound_pose = resulting_pose;
				calculate_interface_unbound_energy(temporary_unbound_pose,scorefxn_,repack_native);
				//store in wt_components_
				store_energies(wt_unbound_components_, (*scorefxn_),temporary_unbound_pose,i,num_iterations_);
				//double debug_unbound_score = (*scorefxn_)(temporary_unbound_pose);
			}
			//store repacked pose
			natives_.push_back(resulting_pose);
			if ( dmp_pdb_ && !restrict_to_nbrhood_ && !output_silent ) {
				std::string dump_repacked_wt = "repacked_wt_round_" + q.str() + ".pdb";
				resulting_pose.dump_pdb(dump_repacked_wt);
			} else if ( (restrict_to_nbrhood_ && dmp_pdb_) || output_silent ) {
				//if you optimize the local neighborhood, automatically dumps out in silent file form
				//otherwise there'll be too many files!!!
				core::io::silent::BinarySilentStruct ss(resulting_pose,"repacked_wt_round_" + q.str());
				sfd.write_silent_struct(ss,silentfilename,false);
			}
		}
	}


	//std::string mutation_label= this->mutation_label(pose);

}

void
ddGMover::setup_packer_task_for_mutations(
	core::pose::Pose const & pose,
	protocols::scoring::Interface const & protein_interface,
	utility::vector1<int> const & mutations,
	core::pack::task::PackerTaskOP packer_task
)
{
	using namespace pose;
	using namespace scoring;
	using namespace conformation;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::task;


	//PACKERTASK SETUP FOR MUTANT
	utility::vector1< bool > neighborhood;
	if ( restrict_to_nbrhood_ ) neighborhood = neighborhood_of_mutations( pose, mutations );

	for ( unsigned int i = 1; i<=pose.total_residue(); i++ ) {
		if ( residues_to_mutate_[i] != core::chemical::aa_unk ) {
			//contains_mutation=true;
			//add option here to restrict to local neighborhood
			utility::vector1<bool> restrict_to_aa(20,false);
			restrict_to_aa[(int)residues_to_mutate_[i]]=true;
			packer_task->nonconst_residue_task(i).restrict_absent_canonical_aas(restrict_to_aa);
			initialize_rotamer_behavior_for_residue_level_task( packer_task->nonconst_residue_task(i) );
		} else {
			if ( restrict_to_nbrhood_ ) {
				if ( neighborhood[i] ) {
					TR << "MUT residue " << i << " is within the cutoff distance to a mutation and will be repacked " << std::endl;
					initialize_rotamer_behavior_for_residue_level_task( packer_task->nonconst_residue_task(i) );
					packer_task->nonconst_residue_task(i).restrict_to_repacking();
				} else {
					TR << "MUT residue " << i << " is not within the cutoff distance to any mutation and will NOT be repacked " << std::endl;
					packer_task->nonconst_residue_task(i).prevent_repacking();
				}
			} else { //not a mutation position, not restricted by distance to mutation site, so allow all to repack
				packer_task->nonconst_residue_task(i).restrict_to_repacking();
				initialize_rotamer_behavior_for_residue_level_task( packer_task->nonconst_residue_task(i) );
			}
		}
	}


	if ( interface_ddg_ ) { // further restrict non-interface aa's if we're in interface_ddg mode
		for ( core::Size i = 1; i<= pose.total_residue(); ++i ) {
			if ( ! protein_interface.is_interface( i ) ) {
				packer_task->nonconst_residue_task(i).prevent_repacking();
			}
		}
	}
	//PACKERTASK SETUP FOR MUTANT END
}

/// @details Experimental: add extra options controlling how the packer behaves
void
ddGMover::initialize_task_level_behavior(
	core::pack::task::PackerTaskOP packer_task
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	packer_task->initialize_extra_rotamer_flags_from_command_line();
	if ( option[ packing::multi_cool_annealer ].user()  ) {
		packer_task->or_multi_cool_annealer( true );
		packer_task->increase_multi_cool_annealer_history_size( option[ packing::multi_cool_annealer ] );
	}
}


/// @details determine if we need to do any work, or if we have checkpoints saving us from some work
bool
ddGMover::is_any_pdb_empty(
	core::pose::Pose const & pose,
	std::string const & mutant_traj,
	bool const output_silent
) const
{
	using namespace pose;
	using namespace scoring;
	using namespace conformation;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::task;

	bool any_pdb_empty = false;
	struct stat filesize;
	if ( !dmp_pdb_ && stat(mutant_traj.c_str(), &filesize) == 0 ) {
		std::cout << "mutation " << this->mutation_label(pose) << " already has a trajectory file. skipping!....." << std::endl;
		any_pdb_empty = false;
	} else if ( !restrict_to_nbrhood_ && !output_silent ) {
		for ( int i =1; i <= num_iterations_; i++ ) {
			std::ostringstream q;
			q << i;
			std::string file_to_check = ("mut_" + this->mutation_label(pose) + "_"
				+ "round_" + q.str() + ".pdb");
			if ( stat(file_to_check.c_str(), &filesize) == 0 && !basic::options::option[OptionKeys::ddg::suppress_checkpointing]() ) {
				//std::cout << "file " << file_to_check << " exists. checking to see if all files non-empty" << std::endl;
				if ( filesize.st_size == 0 ) {
					any_pdb_empty = true;
					std::cout << "file " << file_to_check << " is empty! re-doing calculations " << std::endl;
				}
			} else {
				any_pdb_empty = true;
			}
		}
	} else {
		//silentfile specific behavior when local repacks activated
		std::string file_to_check = "mut_"+this->mutation_label(pose)+".out";;
		if ( !basic::options::option[OptionKeys::ddg::suppress_checkpointing]() ) {

			std::cout << "checking for decoys in silent file: " << file_to_check << std::endl;
			if ( stat(file_to_check.c_str(), &filesize) == 0 ) {
				std::cout << "silent file " << file_to_check << " exists. checking how many tags exist" << std::endl;
				core::io::silent::SilentFileData sfd(file_to_check,false,false,"binary");
				sfd.read_file(file_to_check);
				utility::vector1<std::string> tags = sfd.read_tags_fast(file_to_check);
				if ( /**num tags**/ tags.size() < core::Size(num_iterations_) ) {
					std::cout << "number of tags " << tags.size() << " is less than number of iterations " << num_iterations_ << " re-doing calculations " << std::endl;
					any_pdb_empty = false;
				}
			} else {
				std::cout << "mutant silent file " << file_to_check << " doesn't exist. starting calculations" << std::endl;
				any_pdb_empty=true;
			}
		}
	}
	return any_pdb_empty;
}


utility::vector1< bool >
ddGMover::neighborhood_of_mutations(
	core::pose::Pose const & pose,
	utility::vector1< int > const & mutations
) const
{
	utility::vector1< bool > neighbors( pose.total_residue(), false );
	Real rad2 = nbr_cutoff_ * nbr_cutoff_;

	for ( core::Size  ii = 1; ii <= mutations.size(); ++ii ) {
		Size iiresid = mutations[ ii ];
		core::Vector ii_pos;
		if ( pose.residue(iiresid).name1() == 'G' ) {
			ii_pos = pose.residue(iiresid).xyz(" CA ");
		} else {
			ii_pos = pose.residue(iiresid).xyz(" CB ");
		}
		for ( core::Size jj = 1; jj <= pose.total_residue(); ++jj ) {
			if ( neighbors[ jj ] ) continue;
			core::Vector jj_pos;
			/// Gee -- this will crash if we're dealing with something other than a protein
			if ( pose.residue(jj).name1() == 'G' ) {
				jj_pos = pose.residue(jj).xyz(" CA ");
			} else {
				jj_pos = pose.residue(jj).xyz(" CB ");
			}
			//if c-beta is within x angstrom, set movemap(i) true
			if ( ii_pos.distance_squared(jj_pos) <= rad2 ) {
				neighbors[ jj ] = true;
			}
		}//for all jj
	} // for all mutation positions

	return neighbors;
}

void
ddGMover::initialize_rotamer_behavior_for_residue_level_task(
	core::pack::task::ResidueLevelTask & rltask
) const
{
	rltask.or_include_current(true);
	rltask.or_ex1(true);
	rltask.or_ex2(true);
}


void
ddGMover::apply(core::pose::Pose & pose)
{
	using namespace pose;
	using namespace scoring;
	using namespace conformation;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::task;

	/// This mover operates in two phases:
	/// first, it relaxes the wild type structure in the region of the mutation, then
	/// second, it relaxes the mutant structure in the region of the mutation.  It is critical
	/// that both relaxation protocols be as similar as possible, excluding of course,
	/// the introduction of the non-native amino acid in the protocol modeling the mutant.

	std::string wt_traj = "wt_traj";
	std::string mutant_traj = "mutant_traj"+this->mutation_label(pose);

	// set fa_max_dis to 9.0
	// APL NOTE: this does not work if the Etable has already been created, and it probably has.
	// basic::options::option[ score::fa_max_dis ](9.0);

	protocols::scoring::Interface protein_interface(1); // need to define this outside if, to prevent compilation errors
	//check for properly initialized variables

	//initialize output_silent variable. should go as boolean in private class data?
	bool output_silent = false;
	//force output of silent file
	if ( basic::options::option[OptionKeys::ddg::output_silent].user() ) {
		output_silent = basic::options::option[OptionKeys::ddg::output_silent]();
	}


	if ( this->is_properly_initialized(pose) ) {
		relax_wildtype_structure( pose, protein_interface, wt_traj, output_silent );
	} //for now can only store one. will be last listed mutation

	//also initialize packertask using this array
	pack::task::PackerTaskOP packer_task(pack::task::TaskFactory::create_packer_task(pose));

	bool contains_mutation=false;
	utility::vector1<int> mutations;

	//TR << "restricting to local mutant neighborhood" << std::endl;
	//find mutation position
	//bug?//utility::vector1<int> mutations;
	for ( unsigned int j=1; j <= residues_to_mutate_.size(); j++ ) {
		//bool contains_mutation = false;
		if ( residues_to_mutate_[j] != core::chemical::aa_unk ) {
			contains_mutation = true;
			mutations.push_back(j);
		}
	} //j is mutation position

	setup_packer_task_for_mutations( pose, protein_interface, mutations, packer_task );
	initialize_task_level_behavior( packer_task );

	if ( !contains_mutation ) {
		return; // no mutation specified, so no work to be done.
	}

	//apply all mutations in packertask to pose
	std::ofstream record_mutant_trajectories;

	bool any_pdb_empty = is_any_pdb_empty( pose, mutant_traj, output_silent );

	if ( any_pdb_empty || basic::options::option[OptionKeys::ddg::suppress_checkpointing]() ) {
		TR << "[DEBUG] file does not exist. "<< mutant_traj << ". Creating logfile now.\n";

		record_mutant_trajectories.open( mutant_traj.c_str()  );
		if ( dbg_output_ ) {
			record_mutant_trajectories << "beginning logfile" <<std::endl;
		}

		core::io::silent::SilentFileData sfd;
		std::string silentfilename;
		if ( basic::options::option[out::file::silent].user() ) {
			silentfilename = basic::options::option[out::file::silent]();
			sfd.set_filename(silentfilename);
		} else {
			silentfilename = "mut_"+this->mutation_label(pose)+".out";
			sfd.set_filename(silentfilename); //setup in case we dump in form of silent file
		}

		utility::vector1< std::pair< Real,std::string > > results;
		utility::vector1< pose::PoseOP > poses;


		pose::Pose temporary_pose = pose;
		temporary_pose.remove_constraints();
		temporary_pose.constraint_set(repack_cst_set_mut_);

		pack::pack_rotamers_loop(temporary_pose,(*scorefxn_),packer_task,num_iterations_,results,poses);

		for ( unsigned int i = 1; i <= poses.size(); i++ ) {
			pose::Pose resulting_pose = (*poses[i]);
			resulting_pose.remove_constraints();
			if ( min_cst_ ) {
				resulting_pose.constraint_set(min_cst_set_mut_);
				minimize_with_constraints(resulting_pose,min_cst_sfxn_,packer_task);
				( * min_cst_sfxn_no_cst_weight_)( resulting_pose ); // now compute the score in the absence of the constranits.
				store_energies(mutant_components_, (*min_cst_sfxn_no_cst_weight_),resulting_pose, i,num_iterations_);
			} else {
				resulting_pose.constraint_set(repack_cst_set_mut_);
			}
			Real final_score = -1;
			if ( !min_cst_ ) {
				final_score=(*scorefxn_)( resulting_pose ) ;
			} else {
				final_score=(*min_cst_sfxn_no_cst_weight_)(resulting_pose);
			}

			//store scores
			if ( !interface_ddg_ && !min_cst_ ) {
				store_energies(mutant_components_, (*scorefxn_),resulting_pose, i,num_iterations_);
			} else if ( interface_ddg_ ) {
				assert(pose.chain(pose.total_residue())-pose.chain(1) !=0);
				//dGinterface = Gbound-Gunbound
				//double debug_bound = (*scorefxn_)(resulting_pose);
				store_energies(mutant_components_, (*scorefxn_), resulting_pose,i,num_iterations_);
				pose::Pose temporary_unbound_pose = resulting_pose;
				calculate_interface_unbound_energy(temporary_unbound_pose,scorefxn_,packer_task);
				//store in wt_components_
				//    double debug_unbound = (*scorefxn_)(temporary_unbound_pose);
				store_energies(mutant_unbound_components_, (*scorefxn_),temporary_unbound_pose,i,num_iterations_);
				//TR << "debug unbound " << debug_unbound << " debug bound is " << debug_bound << " so dg is " << (debug_bound-debug_unbound) << std::endl;
			}

			TR << i <<
				" score after mutation: residue " << this->mutation_label(pose) << " "
				<< final_score << " ";
			if ( !min_cst_ ) {
				TR << resulting_pose.energies().total_energies().weighted_string_of( scorefxn_->weights() ) << std::endl;
			} else {
				TR << resulting_pose.energies().total_energies().weighted_string_of( min_cst_sfxn_no_cst_weight_->weights()) << std::endl;
			}

			//output the mutated pdb
			std::ostringstream q;
			q << i;
			//record_mutant_trajectories << "round " << q.str()
			TR << "round " << q.str()
				<< " mutate " << this->mutation_label(pose)
				<< " "  << (final_score) << std::endl;

			if ( dmp_pdb_ && !restrict_to_nbrhood_ && !output_silent ) {
				std::string output_pdb = "mut_" +
					this->mutation_label(pose) + "_" + "round_" + q.str() + ".pdb";
				resulting_pose.dump_pdb(output_pdb);
			} else if ( (dmp_pdb_ && restrict_to_nbrhood_) || output_silent ) {
				core::io::silent::BinarySilentStruct ss(resulting_pose,"mut_" + this->mutation_label(pose) + "_round_" + q.str());
				sfd.write_silent_struct(ss,silentfilename,false);
			}
			mutants_.push_back(resulting_pose);
		}

		TR << "mutate " << this->mutation_label(pose) << //DEBUG
			" " << " "  << " wildtype_dG is: "
			<< this->get_wt_min_totals() << " and mutant_dG is: "
			<< this->get_mutant_min_totals() << " ddG is: " << (this->get_mutant_min_totals()-this->get_wt_min_totals())
			<< std::endl;
		record_mutant_trajectories.close();


	}
}


std::string
ddGMover::get_name() const {
	return "ddGMover";
}

} //moves
} //protocols
