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
/// @file PseudocontactShiftData.hh
///
/// @brief  Hold the PCS data on which the SVD will be applyed
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


#ifndef INCLUDED_protocols_scoring_methods_pcsTs3_PseudocontactShiftData_hh
#define INCLUDED_protocols_scoring_methods_pcsTs3_PseudocontactShiftData_hh

#include <protocols/scoring/methods/pcsTs3/PseudocontactShiftData.fwd.hh>

// Package headers
#include <protocols/scoring/methods/pcsTs3/PseudocontactShiftInput.fwd.hh>
#include <protocols/scoring/methods/pcsTs3/PseudocontactShiftTensor.fwd.hh>

// Project headers
#include <basic/datacache/CacheableData.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Numeric headers
#include <basic/svd/SVD_Solver.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

// C++ headers
#include <iostream>

#ifdef PYROSETTA
#include <protocols/scoring/methods/pcsTs3/PseudocontactShiftInput.hh>
#include <protocols/scoring/methods/pcsTs3/PseudocontactShiftTensor.hh>
#endif


namespace protocols {
namespace scoring {
namespace methods {
namespace pcsTs3 {

class PCS_data_per_lanthanides_Ts3{
private:
	std::string const filename_;
	core::Size n_pcs_;

	utility::vector1<core::Size> A_index_; //index on the giant matrix A to build all the smalls matrix A_

	ObjexxFCL::FArray2D< core::Real > fstyle_A_; //We are going to SVD Ax = b
	ObjexxFCL::FArray1D< core::Real > fstyle_b_; //I should make this one const

	basic::svd::SVD_Solver svd_s_;
	core::Real const weight_;
	core::Real normalization_1_; // SQRT(SUMi( PCS(calc,i)^2 ) )
	core::Real normalization_2_; // Standard deviation
	core::Real normalization_3_; // SQRT(SUMi( PCS(calc,i)^2 )/N )
	core::Real normalization_factor_;
public:
	//PCS_data_per_lanthanides_Ts3(std::string, PCS_file_data_Ts3 & P_f_d, core::Real const weight);

	PCS_data_per_lanthanides_Ts3(std::string, core::Real const weight, utility::vector1< PCS_line_data_Ts3 > & PCS_d_l_a);

private:
	PCS_data_per_lanthanides_Ts3();

public:
	~PCS_data_per_lanthanides_Ts3();

	PCS_data_per_lanthanides_Ts3(PCS_data_per_lanthanides_Ts3 const &other);

	PCS_data_per_lanthanides_Ts3 &
	operator=( PCS_data_per_lanthanides_Ts3 const & other );

	void
	set_A_index(core::Size j, core::Size n_pcs_spin_);

	void
	update_my_A_matrix(utility::vector1< utility::vector1<core::Real> > & A_all);

	core::Real
	get_weight() const;

	/*
	core::Real
	get_normalization_1() const;

	core::Real
	get_normalization_2() const;

	core::Real
	get_normalization_3() const;
	*/

	core::Real
	get_normalization_factor() const;


	//For Debugging purpose
	/*
	void
	print_index_A() const;
	*/


	std::string
	get_filename() const;

	core::Size
	get_n_pcs() const;

	utility::vector1<core::Size> const &
	get_A_index() const;

	ObjexxFCL::FArray1D< core::Real > const &
	get_fstyle_b() const;

	friend
	std::ostream &
	operator << ( std::ostream& out, const PCS_data_per_lanthanides_Ts3 &PCS_d_p_l );

	core::Real
	calculate_tensor_and_cost_with_svd(PCS_tensor_Ts3 &PCS_t);

	//core::Real calculate_tensor_and_cost_with_svd_precalc(PCS_tensor_Ts3 &PCS_t);
};


class PCS_data_Ts3 : public basic::datacache::CacheableData {
private:
	core::Size n_lanthanides_;
	core::Size n_pcs_spin_;
	utility::vector1<PCS_line_data_Ts3> PCS_data_line_all_spin_;
	utility::vector1<PCS_data_per_lanthanides_Ts3> PCS_data_per_lanthanides_all_;
	utility::vector1< utility::vector1<core::Real> > A_all_;
	utility::vector1<core::Real> X_all_;
	utility::vector1<core::Real> Y_all_;
	utility::vector1<core::Real> Z_all_;

public:
	PCS_data_Ts3();

	~PCS_data_Ts3();

	PCS_data_Ts3(PCS_data_input_Ts3 & P_d_i);

	PCS_data_Ts3(PCS_data_input_Ts3 & P_d_i, utility::vector1< bool > const exclude_residues );

	PCS_data_Ts3(PCS_data_Ts3 const &other);

	PCS_data_Ts3 &
	operator=( PCS_data_Ts3 const & src );

	virtual basic::datacache::CacheableDataOP
	clone() const;

	//void update_matrix_fstyle_A();

	core::Size
	get_n_lanthanides() const;

	utility::vector1<core::Real> const &
	get_X_all() const;

	utility::vector1<core::Real> const &
	get_Y_all() const;

	utility::vector1<core::Real> const &
	get_Z_all() const;

	core::Size
	where_is_line(PCS_line_data_Ts3 & P_l_d);

	void
	update_X_Y_Z_all(core::pose::Pose const & pose); //To be called each time the pose is changed

	void
	update_matrix_A_all(core::Real const X,
		core::Real const Y,
		core::Real const Z);

	//void print_matrix_A_all() const;

	//void svd_matrix_A_all();

	utility::vector1<PCS_data_per_lanthanides_Ts3>&
	get_pcs_data_per_lanthanides_all();

	const utility::vector1<PCS_line_data_Ts3> &
	get_PCS_data_line_all_spin() const;

	const utility::vector1<PCS_data_per_lanthanides_Ts3>&
	get_pcs_data_per_lanthanides_all() const;

	friend std::ostream &
	operator<<(std::ostream& out, const PCS_data_Ts3 & P_d);

private:
	void
	update_matrix_A();
};

void
fill_A_line(utility::vector1<core::Real> & A_line,
	core::Real const xM,
	core::Real const yM,
	core::Real const zM,
	core::Real const x,
	core::Real const y,
	core::Real const z
);

}//namespace pcsTs3
}//namespace methods
}//namespace scoring
}//namespace protocols
#endif
