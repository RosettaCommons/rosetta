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
 /// @file protocols/scoring/methods/pcs2/PcsDataCenter.cc
 ///
 /// @brief
 ///
 /// @details
 ///
 /// @param
 ///
 /// @return
 ///
 /// @remarks
 ///
 /// @references
 ///
 /// @authorv Christophe Schmitz
 ///
 ////////////////////////////////////////////////


// Unit headers
#include <protocols/scoring/methods/pcs2/PcsDataCenter.hh>

// Package headers
#include <protocols/scoring/methods/pcs2/PcsInputFile.hh>
#include <protocols/scoring/methods/pcs2/PcsInputCenter.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>

// Numeric headers
//#include <numeric/constants.hh>

// Objexx headers

// C++ headers
#include <iostream>

#include <utility/vector1.hh>


namespace protocols{
namespace scoring{
namespace methods{
namespace pcs2{

static thread_local basic::Tracer TR_PcsDataCenter( "protocols.scoring.methods.pcs.PcsDataCenter" );

PcsDataCenter::PcsDataCenter(){
	utility_exit_with_message( "You shouldn't call the empty constructor for PcsDataCenter class" );
}

PcsDataCenter::~PcsDataCenter(){
}

PcsDataCenter &
PcsDataCenter::operator=( PcsDataCenter const &other )
{
	//	TR_PcsDataCenter << " = called" << std::endl;
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

PcsDataCenter::PcsDataCenter(PcsDataCenter const &other):
	ReferenceCount()
{
	//	TR_PcsDataCenter << " () called" << std::endl;
	n_lanthanides_ = other.n_lanthanides_;
	n_pcs_spin_ = other.n_pcs_spin_;
	PCS_data_line_all_spin_ = other.PCS_data_line_all_spin_;
	PCS_data_per_lanthanides_all_ = other.PCS_data_per_lanthanides_all_;
	A_all_ = other.A_all_;
	X_all_ = other.X_all_;
	Y_all_ = other.Y_all_;
	Z_all_ = other.Z_all_;
}

core::Real
fill_A_line_fast(utility::vector1<core::Real> & A_line,
								 core::Real const xM,
								 core::Real const yM,
								 core::Real const zM,
								 core::Real const x,
								 core::Real const y,
								 core::Real const z){


	core::Real r2;
	core::Real r;
	core::Real r5;
	core::Real x2, y2, z2;
	core::Real v_x, v_y, v_z;

	v_x = x - xM;
	v_y = y - yM;
	v_z = z - zM;
	x2 = v_x * v_x;
	y2 = v_y * v_y;
	z2 = v_z * v_z;

	r2 = x2 + y2 + z2;
	r = sqrt(r2);
	r5 = r2 * r2 * r;

	A_line[1] = (x2 - z2);
	A_line[2] = (2.0 * v_x * v_y);
	A_line[3] = (2.0 * v_x * v_z);
	A_line[4] = (y2 - z2);
	A_line[5] = (2.0 * v_y * v_z);

	return(r5);
}

void
fill_A_line_slow(utility::vector1<core::Real> & A_line,
								 core::Real const xM,
								 core::Real const yM,
								 core::Real const zM,
								 core::Real const x,
								 core::Real const y,
								 core::Real const z){


	//	static const core::Real FACT_USI_PRECALC_FOR_A_3( (10000.0/12.0/ core::Real( numeric::constants::d::pi ) ) * 3.0 );

	core::Real v_x ( x - xM);
	core::Real v_y ( y - yM);
	core::Real v_z ( z - zM);
	core::Real x2  ( v_x * v_x);
	core::Real y2  ( v_y * v_y);
	core::Real z2  ( v_z * v_z);

	core::Real r2 (x2 + y2 + z2);
	core::Real r (sqrt(r2));
	core::Real r5 (r2 * r2 * r);

	core::Real value_1_4_PI_r5 (OTHER_FACT_USI_PRECALC_FOR_A_3 / r5);

	v_x *= 2.0;

	A_line[1] = value_1_4_PI_r5 * (x2 - z2);
	A_line[2] = value_1_4_PI_r5 * (v_x * v_y);
	A_line[3] = value_1_4_PI_r5 * (v_x * v_z);
	A_line[4] = value_1_4_PI_r5 * (y2 - z2);
	A_line[5] = value_1_4_PI_r5 * (2.0 * v_y * v_z);
}


utility::vector1<PcsDataLanthanide> &
PcsDataCenter::get_pcs_data_per_lanthanides_all() {
	return (PCS_data_per_lanthanides_all_);
}

const utility::vector1<PcsInputLine> &
PcsDataCenter::get_PCS_data_line_all_spin() const{
	return (PCS_data_line_all_spin_);
}

const utility::vector1<PcsDataLanthanide>&
PcsDataCenter::get_pcs_data_per_lanthanides_all() const {
	return (PCS_data_per_lanthanides_all_);
}

void
PcsDataCenter::update_matrix_A(){
	core::Size i;
	for( i = 1; i <= PCS_data_per_lanthanides_all_.size(); ++i){
		PCS_data_per_lanthanides_all_[i].update_my_A_matrix_for_svd(A_all_);
	}
}

void
PcsDataCenter::update_matrix_A_cstyle(){
	core::Size i;
	for( i = 1; i <= PCS_data_per_lanthanides_all_.size(); ++i){
		PCS_data_per_lanthanides_all_[i].update_my_A_matrix_for_cstyle(A_all_);
	}
}


void
PcsDataCenter::update_matrix_A_all(core::Real const X,
																	 core::Real const Y,
																	 core::Real const Z){
	core::Size i;

	for( i = 1; i <= n_pcs_spin_; ++i){
		fill_A_line_slow(A_all_[i], X, Y, Z, X_all_[i], Y_all_[i], Z_all_[i]);
	}
}

void
PcsDataCenter::update_matrix_A_all_for_svd(core::Real const X,
																					 core::Real const Y,
																					 core::Real const Z){
	core::Size i;
	core::Real x, y, z;

	for( i = 1; i <= n_pcs_spin_; ++i){
		x = X_all_[i];
		y = Y_all_[i];
		z = Z_all_[i];
		fill_A_line_slow(A_all_[i], X, Y, Z, x, y, z);
	}
	update_matrix_A();
}

void
PcsDataCenter::update_matrix_A_all_for_cstyle(core::Real const X,
																							core::Real const Y,
																							core::Real const Z){
	core::Size i;
	core::Real x, y, z;

	for( i = 1; i <= n_pcs_spin_; ++i){
		x = X_all_[i];
		y = Y_all_[i];
		z = Z_all_[i];
		fill_A_line_slow(A_all_[i], X, Y, Z, x, y, z);
	}
	update_matrix_A_cstyle();
}


//To be called each time the pose is changed
void
PcsDataCenter::update_X_Y_Z_all(core::pose::Pose const & pose){
	core::Size i, res;
	std::string at;
	numeric::xyzVector< core::Real > coo;

	for( i = 1; i <= PCS_data_line_all_spin_.size(); ++i){
		res = PCS_data_line_all_spin_[i].get_residue_num();
		if(res > pose.total_residue()){
			std::cerr << "Error: Couldn't find residue " << res << std::endl;
			std::cerr << "Numbering residue within Rosetta match the sequence provided as input" << std::endl;
			std::cerr << "Make sure the numbering between the sequence and the PseudocontactShift (npc) input file match" << std::endl;
			utility_exit_with_message("Check your pdb and PseudocontactShift (npc) input file");
		}

		at = PCS_data_line_all_spin_[i].get_atom_name();
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
PcsDataCenter::where_is_line(PcsInputLine & pcs_i_l){
	utility::vector1<PcsInputLine>::iterator it;
	core::Size index;
	index = 1;
	for ( it = PCS_data_line_all_spin_.begin(); it != PCS_data_line_all_spin_.end(); ++it ) {
		if ((*it).get_residue_num() == pcs_i_l.get_residue_num()){
			if ((*it).get_atom_name() == pcs_i_l.get_atom_name()){
				return(index);
			}
		}
		index = index + 1;
	}
	return (0);
}

	PcsDataCenter::PcsDataCenter(PcsInputCenter & pcs_i_c, core::Size start, core::Size end, core::Real individual_scale){

	//	TR_PcsDataCenter << " constructor called" << std::endl;

	std::map< std::string, PcsInputFile >::iterator it;

	std::string filename;
	PcsInputLine * pcs_i_l_temp;
	core::Size index, i, j;

	n_pcs_spin_ = 0;
	n_lanthanides_ = 0;

	std::map< std::string, PcsInputFile > & pcs_i_f_all = pcs_i_c.get_PcsInputFile_all();

	for ( it = pcs_i_f_all.begin(); it != pcs_i_f_all.end(); ++it ) {
		filename = it->first;

		PcsInputFile & pcs_i_f = (it->second);
		n_lanthanides_ = n_lanthanides_ + 1;
		//		TR_PcsDataCenter << "Filename " << filename  <<std::endl;
		core::Real weight(pcs_i_f.get_weight());

		utility::vector1<PcsInputLine>::iterator it2;
		utility::vector1<PcsInputLine> & pcs_i_l_all = pcs_i_f.get_PcsInputLine_all();

		PcsDataLanthanide pcs_d_l = PcsDataLanthanide(filename, weight, pcs_i_l_all, start, end, individual_scale);

		j = 1;
		for ( it2 = pcs_i_l_all.begin(); it2 != pcs_i_l_all.end(); ++it2 ){

			pcs_i_l_temp = &(*it2);
			if(do_I_skip(*pcs_i_l_temp, start, end)){
				continue;
			}

			index = where_is_line(*pcs_i_l_temp);
			if(index == 0){
				PCS_data_line_all_spin_.push_back(*it2);
				n_pcs_spin_ = n_pcs_spin_ + 1;
				pcs_d_l.set_A_index(j, n_pcs_spin_);
			}
			else{
				pcs_d_l.set_A_index(j, index);
			}
			j = j + 1;
		}
		PCS_data_per_lanthanides_all_.push_back(pcs_d_l);
	}
	//	TR_PcsDataCenter << "Total spin independent: " <<	n_pcs_spin_ << std::endl;

	A_all_.resize(n_pcs_spin_);
	for(i = 1; i <= n_pcs_spin_; i++){
		(A_all_[i]).resize(5);
	}
	X_all_.resize(n_pcs_spin_);
	Y_all_.resize(n_pcs_spin_);
	Z_all_.resize(n_pcs_spin_);
}


core::Size
PcsDataCenter::get_n_lanthanides() const{
	return (n_lanthanides_);
}


const	utility::vector1< utility::vector1<core::Real> > &
PcsDataCenter::get_A_all() const{
	return(A_all_);
}

const utility::vector1<core::Real> &
PcsDataCenter::get_X_all() const{
	return(X_all_);
}

const utility::vector1<core::Real> &
PcsDataCenter::get_Y_all() const{
	return(Y_all_);
}

const utility::vector1<core::Real> &
PcsDataCenter::get_Z_all() const{
	return(Z_all_);
}

std::ostream &
operator<<(std::ostream& out, const PcsDataCenter & me){
	core::Size i;
	out << "n lanthanides experiment: " << me.get_n_lanthanides() << std::endl;
	out << "n independent spins in total: " << me.n_pcs_spin_ << std::endl;
	for (i = 1 ; i <= me.n_lanthanides_; ++i){
		out << me.PCS_data_per_lanthanides_all_[i] << std::endl;
	}

	return out;
}

}//namespcacs PCS
}//namespace methods
}//namespace scoring
}//namespace protocols
