// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

 //////////////////////////////////////////////
 ///
 /// @file PseudocontactShiftData.cc
 ///
 /// @brief  Hold the PCS data on which the SVD will be applied
 ///
 /// @details
 ///
 /// @param
 ///
 /// @return
 ///
 /// @remarks
 ///
 /// @references C Schmitz et.al. J Mol Biol. Mar 9, 2012; 416(5): 668–677 ; Yagi H et.al Structure, 2013, 21(6):883-890
 ///
 /// @authorv Christophe Schmitz , Kala Bharath Pilla
 ///
 ////////////////////////////////////////////////


// Unit headers
#include <protocols/scoring/methods/pcsTs4/PseudocontactShiftData.hh>

// Package headers
#include <protocols/scoring/methods/pcsTs4/PseudocontactShiftInput.hh>
#include <protocols/scoring/methods/pcsTs4/PseudocontactShiftTensor.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/PCSTS4.OptionKeys.gen.hh>

// Utility headers

// Numeric headers
#include <numeric/constants.hh>

// Objexx headers
#include <ObjexxFCL/FArray1D.hh>

// C++ headers
#include <iostream>
#include <iomanip>

namespace protocols{
namespace scoring{
namespace methods{
namespace pcsTs4{

using namespace ObjexxFCL;

static thread_local basic::Tracer TR_pcs_d_p_l_Ts4( "protocols.scoring.methods.pcsTs4.PCS_data_per_lanthanides_Ts4" );
static thread_local basic::Tracer TR_pcs_d_Ts4( "protocols.scoring.methods.pcsTs4.PCS_data_Ts4" );

PCS_data_per_lanthanides_Ts4::~PCS_data_per_lanthanides_Ts4(){
}

PCS_data_per_lanthanides_Ts4::PCS_data_per_lanthanides_Ts4():
	filename_(""), weight_(0)
{
	utility_exit_with_message( "You shouldn't call the empty constructor for PCS_data_per_lanthanides_Ts4 class" );
}

PCS_data_per_lanthanides_Ts4::PCS_data_per_lanthanides_Ts4(PCS_data_per_lanthanides_Ts4 const &other):
	filename_(other.filename_), svd_s_(other.svd_s_),  weight_(other.weight_)
{
	n_pcs_ = other.n_pcs_;
	A_index_ = other.A_index_;
	fstyle_A_ = other.fstyle_A_;
	fstyle_b_ = other.fstyle_b_;
	normalization_1_ = other.normalization_1_;
	normalization_2_ = other.normalization_2_;
	normalization_3_ = other.normalization_3_;
	normalization_factor_ = other.normalization_factor_;
}

core::Real
PCS_data_per_lanthanides_Ts4::get_weight() const{
	return weight_;
}

core::Real
PCS_data_per_lanthanides_Ts4::get_normalization_factor() const{
	return normalization_factor_;
}


PCS_data_per_lanthanides_Ts4 &
PCS_data_per_lanthanides_Ts4::operator=( PCS_data_per_lanthanides_Ts4 const & other )
{
	if ( this != &other ) {
		n_pcs_ = other.n_pcs_;
		A_index_ = other.A_index_;
		fstyle_A_ = other.fstyle_A_;
		fstyle_b_ = other.fstyle_b_;
		svd_s_ = other.svd_s_;
		normalization_1_ = other.normalization_1_;
		normalization_2_ = other.normalization_2_;
		normalization_3_ = other.normalization_3_;
		normalization_factor_ = other.normalization_factor_;
	}
	return *this;
}

PCS_data_Ts4::PCS_data_Ts4(){
	utility_exit_with_message( "You shouldn't call the empty constructor for PCS_data_Ts4 class" );
}

PCS_data_Ts4::~PCS_data_Ts4(){
}

PCS_data_Ts4 &
PCS_data_Ts4::operator=( PCS_data_Ts4 const &other )
{
	if ( this != &other ) {
		n_lanthanides_ = other.n_lanthanides_;
		n_pcs_spin_ = other.n_pcs_spin_;
		PCS_data_line_all_spin_ = other.PCS_data_line_all_spin_;
		PCS_data_per_lanthanides_all_ = other.PCS_data_per_lanthanides_all_;
		A_all_ = other.A_all_;
		X_all_ = other.X_all_;
		Z_all_ = other.Z_all_;
		Y_all_ = other.Y_all_;
	}
	return *this;
}

PCS_data_Ts4::PCS_data_Ts4(PCS_data_Ts4 const &other):
	CacheableData()
{
	n_lanthanides_ = other.n_lanthanides_;
	n_pcs_spin_ = other.n_pcs_spin_;
	PCS_data_line_all_spin_ = other.PCS_data_line_all_spin_;
	PCS_data_per_lanthanides_all_ = other.PCS_data_per_lanthanides_all_;
	A_all_ = other.A_all_;
	X_all_ = other.X_all_;
	Y_all_ = other.Y_all_;
	Z_all_ = other.Z_all_;
}

basic::datacache::CacheableDataOP
PCS_data_Ts4::clone() const {
	return basic::datacache::CacheableDataOP( new PCS_data_Ts4( *this ) );
}

void
fill_A_line(utility::vector1<core::Real> & A_line,
						core::Real const xM,
						core::Real const yM,
						core::Real const zM,
						core::Real const x,
						core::Real const y,
						core::Real const z){

	using namespace core;

	core::Real value_1_4_PI_r5;
	core::Real r2;
	core::Real r;
	core::Real r5;
	core::Real x2, y2, z2;
	core::Real v_x, v_y, v_z;

	static const core::Real FACT_USI_PRECALC_FOR_A_3( (10000.0/12.0/ core::Real( numeric::constants::d::pi ) ) * 3.0 );

	v_x = x - xM;
	v_y = y - yM;
	v_z = z - zM;
	x2 = v_x * v_x;
	y2 = v_y * v_y;
	z2 = v_z * v_z;

	r2 = x2 + y2 + z2;
	r = sqrt(r2);
	r5 = r2 * r2 * r;

	value_1_4_PI_r5 = FACT_USI_PRECALC_FOR_A_3 / r5;

	A_line[1] = value_1_4_PI_r5 * (x2 - z2);
	A_line[2] = value_1_4_PI_r5 * (2.0 * v_x * v_y);
	A_line[3] = value_1_4_PI_r5 * (2.0 * v_x * v_z);
	A_line[4] = value_1_4_PI_r5 * (y2 - z2);
	A_line[5] = value_1_4_PI_r5 * (2.0 * v_y * v_z);
}

void
PCS_data_per_lanthanides_Ts4::update_my_A_matrix(utility::vector1< utility::vector1<core::Real> > & A_all){
	core::Size i, j;
	for(i = 1; i <= A_index_.size(); ++i){
		for(j = 1; j <= 5; ++j){
			fstyle_A_(i, j) = A_all[A_index_[i]][j];
		}
	}
	svd_s_.set_matrix_A(fstyle_A_);
}

core::Real
PCS_data_per_lanthanides_Ts4::calculate_tensor_and_cost_with_svd(PCS_tensor_Ts4 &PCS_t){

	core::Real score;

	//#define FAST_SCORING
	//FAST SCORING is actually not faster (nor slower)
	//FAST SCORING gives the PCS score without actually solving Ax = b
	//(but still, the matrix has to be decomposed.
  //So FAST_SCORING doesn't give the tensor parameters...
	//FAST SCORING might become handy when working on unassigned data, that's why I keep it.

	svd_s_.run_decomp_svd();

#ifdef FAST_SCORING
	score = svd_s_.run_score_svd_without_solving();
	PCS_t.reset_tensor(0, 0, 0, 0, 0); //Junk tensor to avoid warning at compilation
#else
	svd_s_.run_solve_svd();
	score = svd_s_.run_score_svd_on_matrix(fstyle_A_);
	// The 2 following calls are not that necessary for just scoring
	//	FArray1D<core::Real> f (svd_s_.get_svd_solution());

	utility::vector1<core::Real> const & f (svd_s_.get_svd_solution());

	//	PCS_t.reset_tensor(f(1), f(2), f(3), f(4), f(5));
	PCS_t.reset_tensor(f[1], f[2], f[3], f[4], f[5]);
#endif

	//	std::cerr << weight_ << "applyed as individual weight in svd" << std::endl;

	score *= weight_;

	return(score/normalization_factor_);
}

utility::vector1<PCS_data_per_lanthanides_Ts4> &
PCS_data_Ts4::get_pcs_data_per_lanthanides_all() {
	return (PCS_data_per_lanthanides_all_);
}

const utility::vector1<PCS_line_data_Ts4> &
PCS_data_Ts4::get_PCS_data_line_all_spin() const{
	return (PCS_data_line_all_spin_);
}

const utility::vector1<PCS_data_per_lanthanides_Ts4>&
PCS_data_Ts4::get_pcs_data_per_lanthanides_all() const {
	return (PCS_data_per_lanthanides_all_);
}

void
PCS_data_Ts4::update_matrix_A(){
	core::Size i;
	for( i = 1; i <= PCS_data_per_lanthanides_all_.size(); ++i){
		PCS_data_per_lanthanides_all_[i].update_my_A_matrix(A_all_);
	}
}

void
PCS_data_Ts4::update_matrix_A_all(core::Real const X,
															core::Real const Y,
															core::Real const Z){
	core::Size i;
	core::Real x, y, z;

	for( i = 1; i <= n_pcs_spin_; ++i){
		x = X_all_[i];
		y = Y_all_[i];
		z = Z_all_[i];
		fill_A_line(A_all_[i], X, Y, Z, x, y, z);
	}
	update_matrix_A();
}


//To be called each time the pose is changed
void
PCS_data_Ts4::update_X_Y_Z_all(core::pose::Pose const & pose){
	core::Size i, res;
	std::string at;
	numeric::xyzVector< core::Real > coo;

	for( i = 1; i <= PCS_data_line_all_spin_.size(); ++i){
		res = PCS_data_line_all_spin_[i].residue_num();
		if(res > pose.total_residue()){
			std::cerr << "Error: Couldn't find residue " << res << std::endl;
			std::cerr << "Numbering residue within Rosetta match the sequence provided as input" << std::endl;
			std::cerr << "Make sure the numbering between the sequence and the PseudocontactShift (npc) input file match" << std::endl;
			utility_exit_with_message("Check your pdb and PseudocontactShift (npc) input file");
		}

		at = PCS_data_line_all_spin_[i].atom_name();
		if( ! pose.residue(res).has(at)){
			std::cerr << "Error: Couldn't find the atom " << at << " in residue " << res << std::endl;
			std::cerr << "Numbering residue within Rosetta match the sequence provided as input" << std::endl;
			std::cerr << "Make sure the numbering between the sequence and the PseudocontactShift (npc) input file match" << std::endl;
			std::cerr << "Use only PCS for the backbone for abinitio." << std::endl;
			std::cerr << "only N, CA, C, O, CB, H and CEN" << std::endl;

			utility_exit_with_message("Check your pdb and PseudocontactShift (npc) input file");
		}
		coo = pose.residue(res).atom(at).xyz();
		X_all_[i] = coo.x();
		Y_all_[i] = coo.y();
		Z_all_[i] = coo.z();
	}
}

core::Size
PCS_data_per_lanthanides_Ts4::get_n_pcs() const{
	return (n_pcs_);
}


const utility::vector1<core::Size> &
PCS_data_per_lanthanides_Ts4::get_A_index() const{
	return(A_index_);
}

const FArray1D< core::Real > &
PCS_data_per_lanthanides_Ts4::get_fstyle_b() const{
	return(fstyle_b_);
}


core::Size
PCS_data_Ts4::where_is_line(PCS_line_data_Ts4 & P_l_d){
	utility::vector1<PCS_line_data_Ts4>::iterator it;
	core::Size index;
	index = 1;
	for ( it = PCS_data_line_all_spin_.begin(); it != PCS_data_line_all_spin_.end(); ++it ) {
		if ((*it).residue_num() == P_l_d.residue_num()){
			if ((*it).atom_name() == P_l_d.atom_name()){
				return(index);
			}
		}
		index = index + 1;
	}
	return (0);
}


PCS_data_per_lanthanides_Ts4::PCS_data_per_lanthanides_Ts4(std::string filename,
																									 core::Real const weight,
																									 utility::vector1< PCS_line_data_Ts4 > & PCS_d_l_a):
		filename_(filename), svd_s_(basic::svd::SVD_Solver(PCS_d_l_a.size(), 5)), weight_(weight)
{
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::Size i;
	utility::vector1<PCS_line_data_Ts4>::iterator it;

	n_pcs_ =  PCS_d_l_a.size();
	A_index_.resize(n_pcs_);

	core::Size M, N;
	M = n_pcs_;
	N = 5;

	svd_s_ = basic::svd::SVD_Solver(M, N);
	utility::vector1<core::Real> vec_temp;
	fstyle_A_.dimension(M, N);
	fstyle_b_.dimension(M);

	normalization_1_ = 0;
	core::Real sum (0);

	for (i = 1; i <= PCS_d_l_a.size(); ++i){
		PCS_line_data_Ts4 PCS_l_d = PCS_d_l_a[i];
		vec_temp.push_back(PCS_l_d.PCS_experimental());
		fstyle_b_(i) = PCS_l_d.PCS_experimental();
		normalization_1_ += PCS_l_d.PCS_experimental() * PCS_l_d.PCS_experimental();
		sum += PCS_l_d.PCS_experimental();
	}
	svd_s_.set_vector_b(vec_temp);
	normalization_1_ = sqrt(normalization_1_);

	normalization_3_ = sqrt(normalization_1_ / PCS_d_l_a.size());

	core::Real average(sum/PCS_d_l_a.size());

	normalization_2_ = 0;
	for (i = 1; i <= PCS_d_l_a.size(); ++i){
		normalization_2_ += (fstyle_b_(i) - average) *  (fstyle_b_(i) - average);
	}
	normalization_2_ = sqrt( normalization_2_ / PCS_d_l_a.size() );

	TR_pcs_d_p_l_Ts4 << "n_pcs: " << n_pcs_ << std::endl;
	TR_pcs_d_p_l_Ts4 << "Normalization 1: " << normalization_1_ << std::endl;
	TR_pcs_d_p_l_Ts4 << "Normalization 2: " << normalization_2_ << std::endl;
	TR_pcs_d_p_l_Ts4 << "Normalization 3: " << normalization_3_ << std::endl;

	if( option[ basic::options::OptionKeys::PCSTS4::normalization_id ].user() ){
		core::Size norma_id (option[ basic::options::OptionKeys::PCSTS4::normalization_id ]());

		switch (norma_id){
		case 1:{
			TR_pcs_d_p_l_Ts4 << "Using Normalization '1': " << normalization_1_ << ". Normalize each data set by SQRT( SUM( PCScalc(i)^2 ) ) " << std::endl;
			normalization_factor_ = normalization_1_;
			break;
		}
		case 2:{
			TR_pcs_d_p_l_Ts4 << "Using Normalization '2': " << normalization_2_ << ". Normalize each data set by Standard Deviation" << std::endl;
			normalization_factor_ = normalization_2_;
			break;
		}
		case 3:{
			TR_pcs_d_p_l_Ts4 << "Using Normalization '3': " << normalization_3_ << ". Normalize each data set by SQRT( SUM( PCScalc(i)^2 ) / N) " << std::endl;
			normalization_factor_ = normalization_3_;
			break;
		}
		default:{
			TR_pcs_d_p_l_Ts4 << "Normalization '"<< norma_id << "' id not recognized " << std::endl;
			utility_exit_with_message("You should use a valid normalization id ('1' or '2' or '3')");
			break;
		}
		}
	}
	else{
		TR_pcs_d_p_l_Ts4 << "Normalization NOT used in calculations " << std::endl;
		normalization_factor_ = 1.0;
	}
}

std::string
PCS_data_per_lanthanides_Ts4::get_filename() const {
	return (filename_);
}

void
PCS_data_per_lanthanides_Ts4::set_A_index(core::Size j,
																			core::Size index){
	A_index_[j] = index;
}

PCS_data_Ts4::PCS_data_Ts4(PCS_data_input_Ts4 & P_d_i){
	std::map< std::string, PCS_file_data_Ts4 >::iterator it;

	std::string filename;
	PCS_line_data_Ts4 * P_l_d_temp;
	core::Size index, i, j;

	n_pcs_spin_ = 0;
	n_lanthanides_ = 0;

	std::map< std::string, PCS_file_data_Ts4 > & P_f_a_d = P_d_i.get_PCS_data_input_reference();

	for ( it = P_f_a_d.begin(); it != P_f_a_d.end(); ++it ) {
		filename = it->first;

		PCS_file_data_Ts4 & P_f_d = (it->second);
		n_lanthanides_ = n_lanthanides_ + 1;
		TR_pcs_d_Ts4 << "Filename " << filename  <<std::endl;
		core::Real weight(P_f_d.get_weight());

		utility::vector1<PCS_line_data_Ts4>::iterator it2;
		utility::vector1<PCS_line_data_Ts4> & PCS_d_l_a = P_f_d.get_PCS_data_line_all_reference();

		PCS_data_per_lanthanides_Ts4 P_d_p_l = PCS_data_per_lanthanides_Ts4(filename, weight, PCS_d_l_a );

		j = 1;
		for ( it2 = PCS_d_l_a.begin(); it2 != PCS_d_l_a.end(); ++it2 ){
			P_l_d_temp = &(*it2);
			index = where_is_line(*P_l_d_temp);
			if(index == 0){
				PCS_data_line_all_spin_.push_back(*it2);
				n_pcs_spin_ = n_pcs_spin_ + 1;
				P_d_p_l.set_A_index(j, n_pcs_spin_);
			}
			else{
				P_d_p_l.set_A_index(j, index);
			}
			j = j + 1;
		}
		PCS_data_per_lanthanides_all_.push_back(P_d_p_l);
	}
	TR_pcs_d_Ts4 << "Total spin independent: " <<	n_pcs_spin_ << std::endl;

	A_all_.resize(n_pcs_spin_);
	for(i = 1; i <= n_pcs_spin_; i++){
		(A_all_[i]).resize(5);
	}
	X_all_.resize(n_pcs_spin_);
	Y_all_.resize(n_pcs_spin_);
	Z_all_.resize(n_pcs_spin_);
}

PCS_data_Ts4::PCS_data_Ts4(PCS_data_input_Ts4 & P_d_i, utility::vector1< bool > const exclude_residues ){
	std::map< std::string, PCS_file_data_Ts4 >::iterator it;

	std::string filename;
	PCS_line_data_Ts4 * P_l_d_temp;
	core::Size index, i, j;

	n_pcs_spin_ = 0;
	n_lanthanides_ = 0;

	std::map< std::string, PCS_file_data_Ts4 > & P_f_a_d = P_d_i.get_PCS_data_input_reference();

	for ( it = P_f_a_d.begin(); it != P_f_a_d.end(); ++it ) {
		filename = it->first;

		PCS_file_data_Ts4 & P_f_d = (it->second);
		n_lanthanides_ = n_lanthanides_ + 1;
		TR_pcs_d_Ts4 << "Filename " << filename  <<std::endl;
		core::Real weight(P_f_d.get_weight());

		utility::vector1< PCS_line_data_Ts4 >::iterator it2;
		utility::vector1< PCS_line_data_Ts4 > & file_PCS_d_l_a = P_f_d.get_PCS_data_line_all_reference();
		utility::vector1< PCS_line_data_Ts4 > PCS_d_l_a;

		for ( it2 = file_PCS_d_l_a.begin(); it2 != file_PCS_d_l_a.end(); ++it2 ){
			bool use_residue = false;

			if ( it2->residue_num() > exclude_residues.size() ) {
				use_residue = true;
			} else {
				if ( exclude_residues[it2->residue_num()] == false ) {
					use_residue = true;
				}
			}

			if (use_residue == true) {
				PCS_d_l_a.push_back(*it2);
			}
		}

		PCS_data_per_lanthanides_Ts4 P_d_p_l = PCS_data_per_lanthanides_Ts4(filename, weight, PCS_d_l_a );

		j = 1;
		for ( it2 = PCS_d_l_a.begin(); it2 != PCS_d_l_a.end(); ++it2 ){
			P_l_d_temp = &(*it2);
			index = where_is_line(*P_l_d_temp);
			if(index == 0){
				PCS_data_line_all_spin_.push_back(*it2);
				n_pcs_spin_ = n_pcs_spin_ + 1;
				P_d_p_l.set_A_index(j, n_pcs_spin_);
			}
			else{
				P_d_p_l.set_A_index(j, index);
			}
			j = j + 1;
		}
		PCS_data_per_lanthanides_all_.push_back(P_d_p_l);
	}
	TR_pcs_d_Ts4 << "Total spin independent: " <<	n_pcs_spin_ << std::endl;

	A_all_.resize(n_pcs_spin_);
	for(i = 1; i <= n_pcs_spin_; i++){
		(A_all_[i]).resize(5);
	}
	X_all_.resize(n_pcs_spin_);
	Y_all_.resize(n_pcs_spin_);
	Z_all_.resize(n_pcs_spin_);
}


core::Size
PCS_data_Ts4::get_n_lanthanides() const{
	return (n_lanthanides_);
}

const utility::vector1<core::Real> &
PCS_data_Ts4::get_X_all() const{
	return(X_all_);
}

const utility::vector1<core::Real> &
PCS_data_Ts4::get_Y_all() const{
	return(Y_all_);
}

const utility::vector1<core::Real> &
PCS_data_Ts4::get_Z_all() const{
	return(Z_all_);
}

std::ostream &
operator<<(std::ostream& out, const PCS_data_Ts4 & P_d){
	core::Size i;
	out << "n lanthanides experiment: " << P_d.get_n_lanthanides() << std::endl;
	out << "n independent spins in total: " << P_d.n_pcs_spin_ << std::endl;
	for (i = 1 ; i <= P_d.n_lanthanides_; ++i){
		out << P_d.PCS_data_per_lanthanides_all_[i] << std::endl;
	}

	return out;
}

std::ostream &
operator<<(std::ostream& out, const PCS_data_per_lanthanides_Ts4 &PCS_d_p_l){
	core::Size i, j;
	out << std::setprecision(4) << std::endl;
	out << "  Filename : " << PCS_d_p_l.get_filename() << std::endl;
	out << "  Number of pcs : " << PCS_d_p_l.get_n_pcs() << std::endl;
	out << "  b vector (pcs values)  : " << std::endl;

	for ( i = 1; i <= PCS_d_p_l.get_n_pcs(); ++i){
		out << "  " << std::setw(4) << i << ":"<< std::setw(10) <<PCS_d_p_l.fstyle_b_(i)  << std::endl;
	}
	out << std::endl;

	out << "  A_index  : " << std::endl;
	for ( i = 1; i <= PCS_d_p_l.get_n_pcs(); ++i){
		out << "  " << std::setw(4) << i << ":"  << std::setw(4) <<PCS_d_p_l.A_index_[i];
		out << std::endl;
	}


	out << "  Matrix A:"<< std::endl;
	for( i = 1; i <= PCS_d_p_l.n_pcs_ ; ++i ){
		out << "  " << std::setw(4) << i << ":";
		for (j = 1; j <= 5; ++j){
			out << " " << std::setw(12) << PCS_d_p_l.fstyle_A_( i, j);
		}
		out << std::endl;
	}
	return out;
}

}//namespcacs pcsTs4
}//namespace methods
}//namespace scoring
}//namespace protocols
