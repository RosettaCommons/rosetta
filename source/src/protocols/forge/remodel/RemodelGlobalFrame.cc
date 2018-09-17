// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/remodel/RemodelLoopMover.cc
/// @brief  Loop modeling protocol based on routines from Remodel and EpiGraft
///         packages in Rosetta++.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)
/// @author Possu Huang (possu@u.washington.edu)

// unit headers
#include <protocols/forge/remodel/RemodelGlobalFrame.hh>
#include <protocols/forge/remodel/RemodelRotamerLinks.hh>
#include <protocols/forge/methods/util.hh>

// package headers

// project headers
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.hh>
#include <core/util/disulfide_util.hh>
#include <core/conformation/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/io/Remarks.hh>
#include <core/io/StructFileRep.hh>
#include <core/pose/util.hh> // for pdbinfo
#include <core/id/AtomID.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/io/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/symmetry/SetupForSymmetryMover.hh>
#include <protocols/constraints_additional/BindingSiteConstraint.hh>
#include <protocols/constraints_additional/COMCoordinateConstraint.hh>

// numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//external
#include <Eigen/Dense>

// boost headers

// C++ headers
#include <iostream>
#include <cmath>

using namespace basic::options;
using namespace Eigen;

namespace protocols {
namespace forge {
namespace remodel {

// Tracer instance for this file
// Named after the original location of this code
static basic::Tracer TR( "protocols.forge.remodel.RemodelGlobalFrame" );

// RNG


// @brief default constructor
RemodelGlobalFrame::RemodelGlobalFrame()
{
	// has to reinitialize state before apply
	//state_.clear();
}

/// @brief value constructor
RemodelGlobalFrame::RemodelGlobalFrame(RemodelData const & remodel_data,
	RemodelWorkingSet const & working_model,
	ScoreFunctionOP const & sfxn)
{

	remodel_data_ = remodel_data;
	seg_size = remodel_data.blueprint.size();
	working_model_ = working_model;
	score_fxn_ = sfxn->clone();
	left_handed_ = 0;
	/*
	if(option[OptionKeys::remodel::repeat_structure].user()){
	Size repeatCount =option[OptionKeys::remodel::repeat_structure];
	for (Size rep = 0; rep < repeatCount ; rep++){
	for (std::set< core::Size >::iterator it = uup.begin(); it != uup.end(); ++it){
	//DEBUG
	// std::cout << *it + remodel_data.blueprint.size()*rep << std::endl;
	// std::cout << "manger size"  << working_model.manager.union_of_intervals_containing_undefined_positions().size() <<  std::endl;
	// std::cout << *it  << std::endl;
	if ( !(*it+remodel_data.blueprint.size()*rep > (remodel_data.blueprint.size()*repeatCount)) ){ //Extrapolation of positions shouldn't go beyond the length of pose
	und_pos.insert(*it + remodel_data.blueprint.size()*rep);
	}
	}
	}
	}
	else {
	und_pos = working_model.manager.union_of_intervals_containing_undefined_positions();
	}

	*/

}

//simpel constructor
RemodelGlobalFrame::RemodelGlobalFrame(Size segment_size){
	seg_size = segment_size;
	left_handed_ = 0;
}

/// @brief copy constructor

/// @brief default destructor
RemodelGlobalFrame::~RemodelGlobalFrame()= default;

/// @brief clone this object
protocols::moves::MoverOP
RemodelGlobalFrame::clone() const {
	return protocols::moves::MoverOP( new RemodelGlobalFrame( *this ) );
}

/// @brief create this type of object
protocols::moves::MoverOP
RemodelGlobalFrame::fresh_instance() const {
	return protocols::moves::MoverOP( new RemodelGlobalFrame() );
}

void RemodelGlobalFrame::apply( core::pose::Pose & pose )
{
	using namespace basic::options;

	get_helical_params(pose);

}

void RemodelGlobalFrame::get_helical_params( core::pose::Pose & pose ) {
	using Eigen::MatrixXd;
	using namespace Eigen;
	using namespace std;
	using namespace basic::options;

	//Size numRes = pose.size();  // unused ~Labonte

	// dumping information into PDB header
	core::pose::PDBInfoOP temp_pdbinfo( new core::pose::PDBInfo(pose,true) );

	core::io::RemarkInfo remark;
	//capture stream
	std::stringstream capture_stream;

	if ( option[OptionKeys::remodel::repeat_structure].user() ) {}
	else {
		TR << "only applicable in Repeat mode";
		exit(0);
	}

	// get the size of each repeating sector
	//Size seg_size = remodel_data_.blueprint.size();
	//TR << "debug: seg_size = " << seg_size << std::endl;

	MatrixXf A(3, seg_size);
	MatrixXf B(3, seg_size);

	for ( Size i = 1; i <= seg_size ; i++ ) {
		core::Vector coord = pose.residue(i).xyz("CA");
		A(0,i-1) = coord[0];
		A(1,i-1) = coord[1];
		A(2,i-1) = coord[2];

		//TR << "resA " << i << " has " << coord[0] << " " << coord[1]  << " " << coord[2] << std::endl;
	}
	for ( Size i = seg_size+1; i <= seg_size*2 ; i++ ) {
		core::Vector coord = pose.residue(i).xyz("CA");
		B(0,i-1-seg_size) = coord[0];
		B(1,i-1-seg_size) = coord[1];
		B(2,i-1-seg_size) = coord[2];
		//TR << "resB " << i << " has " << coord[0] << " " << coord[1]  << " " << coord[2] << std::endl;
	}


	Vector3f c_A(A.row(0).mean(), A.row(1).mean(), A.row(2).mean());
	Vector3f c_B(B.row(0).mean(), B.row(1).mean(), B.row(2).mean());
	MatrixXf x_A(A.rows(),A.cols());
	MatrixXf x_B(B.rows(),B.cols());

	for ( int i=0; i<A.cols(); i++ ) {
		x_A.col(i)=A.col(i)-c_A;
		x_B.col(i)=B.col(i)-c_B;
	}

	Matrix3f cov= (x_B * x_A.transpose()) / x_A.cols();
	JacobiSVD<MatrixXf> svd(cov, ComputeFullU | ComputeFullV);

	Matrix3f Rt=svd.matrixU() * svd.matrixV().transpose();
	Matrix3f R;
	R<< 1,0,0, 0,1,0, 0,0,Rt.determinant();

	Matrix3f H= svd.matrixU() * R * svd.matrixV().transpose();
	capture_stream <<"Rotation matrix:\n"<<H<<endl;

	double omega = acos((H.trace()-1)/2);
	Matrix3f I = Matrix3f::Identity();
	Matrix3f N = 0.5*(H+H.transpose()) - cos(omega)*I;

	Vector3f hN;
	Real scalar = 0;
	Real max_scalar = -10000000;
	for ( Size i = 0; i<=2; i++ ) {
		scalar = N.col(i).norm();
		if ( scalar > max_scalar ) {
			max_scalar = scalar;
			hN = N.col(i)/N.col(i).norm();
		}
	}

	double sin_omega = (H(1,0)-H(0,1)) / (2*hN(2));
	TR.Debug <<"sin_omega "<<sin_omega<<endl;
	if ( sin_omega < 0 ) hN = -1 * hN;


	Vector3f t = c_B - H*c_A;
	double L = t.dot(hN) ;
	double rise= std::abs(L);

	capture_stream<<"getParams:N:\n"<<N<<endl;
	capture_stream<<"getParams:helical axis:\n"<<hN<<endl;


	Matrix3f Ncross;
	Ncross << 0,-1*hN(2),hN(1), hN(2),0,-1*hN(0), -1*hN(1),hN(0),0 ;
	Matrix3f R0t= (1-cos(omega))*I - sin(omega)*Ncross;
	Vector3f R0 = R0t.inverse() * (t-L*hN);

	TR.Debug << "R0" << std::endl << R0 << std::endl;
	capture_stream <<"R0 (helical axis position):\n"<<R0<<endl;
	/*
	cout<<"helical axis:\n"<<hN<<endl<<endl;
	cout<<"t:\n"<<t<<endl<<endl;
	cout<<"L:\n"<<L<<endl<<endl;
	cout<<"R0t:\n"<<R0t<<endl<<endl;
	cout<<"R0t_inv:\n"<<R0t.inverse()<<endl<<endl;
	cout<<"helical axis:\n"<<hN<<endl<<endl;
	*/

	Vector3f pA= (c_A-R0)-(hN*(hN.dot(c_A-R0)));
	Vector3f pB= (c_B-R0)-(hN*(hN.dot(c_B-R0)));


	double direction = L * hN.dot(pA.cross(pB));
	left_handed_ = 1;
	if ( direction > 0 ) left_handed_ = -1;
	capture_stream<<"left handed: "<<left_handed_<<endl;

	double radius=pA.norm();

	capture_stream <<"L: " << L << endl << "rise: "<<rise<<endl<<"radius: "<<radius<<endl<<"omega: "<<omega<<endl;

	rise_=rise;
	radius_=radius;
	omega_=omega;

	std::cout << capture_stream.str()<< endl;

	while ( capture_stream.good() ) {
		std::string one_line;
		getline(capture_stream, one_line, '\n');
		remark.value = one_line ;
		temp_pdbinfo->remarks().push_back( remark );
	}

	pose.pdb_info(temp_pdbinfo);

}


// extract parameters from repeating segments and align the first to X axis on XZ plane
void RemodelGlobalFrame::align_segment( core::pose::Pose & pose ) {
	using Eigen::MatrixXd;
	using namespace Eigen;
	using namespace std;
	using namespace numeric;
	using namespace basic::options;
	using namespace core::pose::symmetry;
	using namespace core::scoring::constraints;
	using namespace protocols;

	TR.Debug << "align seg 1" << std::endl;


	bool is_sym = false;
	if ( is_symmetric( pose ) ) {  // if symmetrical, has to extract monomer for the transformation
		is_sym = true;
		ConstraintSetOP pose_cst_set( new ConstraintSet( *pose.constraint_set() ) );
		Pose junk_for_copy;
		extract_asymmetric_unit( pose, junk_for_copy, false);
		pose=junk_for_copy;
		pose.constraint_set(pose_cst_set);
		pose.pdb_info()->obsolete(true);
	}


	//Size numRes = pose.size();  // unused ~Labonte

	if ( option[OptionKeys::remodel::repeat_structure].user() ) {}
	else {
		TR << "only applicable in Repeat mode";
		exit(0);
	}

	//////////////////////////////////
	//  compute helical axis with CA
	//////////////////////////////////

	// get the size of each repeating sector
	//moved this to constructor
	//Size seg_size = remodel_data_.blueprint.size();

	TR.Debug << "align seg 2" << std::endl;
	TR.Debug << "pose length" << pose.size() << std::endl;
	TR.Debug << "seg_size" << seg_size << std::endl;

	MatrixXf A(3,seg_size);
	MatrixXf B(3,seg_size);

	for ( Size i = 1; i <= seg_size ; i++ ) {
		core::Vector coord = pose.residue(i).xyz("CA");
		A(0,i-1) = coord[0];
		A(1,i-1) = coord[1];
		A(2,i-1) = coord[2];

	}
	for ( Size i = seg_size+1; i <= seg_size*2 ; i++ ) {
		core::Vector coord = pose.residue(i).xyz("CA");
		B(0,i-1-seg_size) = coord[0];
		B(1,i-1-seg_size) = coord[1];
		B(2,i-1-seg_size) = coord[2];
	}

	TR.Debug << "align seg 3" << std::endl;

	Vector3f c_A(A.row(0).mean(), A.row(1).mean(), A.row(2).mean());
	Vector3f c_B(B.row(0).mean(), B.row(1).mean(), B.row(2).mean());
	MatrixXf x_A(A.rows(),A.cols());
	MatrixXf x_B(B.rows(),B.cols());

	for ( int i=0; i<A.cols(); i++ ) {
		x_A.col(i)=A.col(i)-c_A;
		x_B.col(i)=B.col(i)-c_B;
	}

	Matrix3f cov= (x_B * x_A.transpose()) / x_A.cols();
	JacobiSVD<MatrixXf> svd(cov, ComputeFullU | ComputeFullV);

	Matrix3f Rt=svd.matrixU() * svd.matrixV().transpose();
	Matrix3f R;
	R<< 1,0,0, 0,1,0, 0,0,Rt.determinant();

	Matrix3f H= svd.matrixU() * R * svd.matrixV().transpose();

	TR.Debug << "H: " << H << std::endl;

	double omega = acos((H.trace()-1)/2);
	Matrix3f I = Matrix3f::Identity();
	Matrix3f N = 0.5*(H+H.transpose()) - cos(omega)*I;

	Vector3f hN;
	Real scalar = 0;
	Real max_scalar = -10000000;
	for ( Size i = 0; i<=2; i++ ) {
		scalar = N.col(i).norm();
		if ( scalar > max_scalar ) {
			max_scalar = scalar;
			hN = N.col(i)/N.col(i).norm();
		}
	}

	TR.Debug << "N: " << N << std::endl;
	TR.Debug << "helical axis: " << hN << std::endl;

	double sin_omega = (H(1,0)-H(0,1)) / (2*hN(2));

	TR.Debug<<"sin_omega "<<sin_omega<<endl;

	if ( sin_omega < 0 ) hN = -1 * hN;

	Vector3f t = c_B - H*c_A;
	double L = t.dot(hN) ;
	//double rise= std::abs(L);  // unused ~Labonte


	Matrix3f Ncross;
	Ncross << 0,-1*hN(2),hN(1), hN(2),0,-1*hN(0), -1*hN(1),hN(0),0 ;
	Matrix3f R0t= (1-cos(omega))*I - sin(omega)*Ncross;
	Vector3f R0 = R0t.inverse() * (t-L*hN);

	Vector3f pA= (c_A-R0)-(hN*(hN.dot(c_A-R0)));
	Vector3f pB= (c_B-R0)-(hN*(hN.dot(c_B-R0)));

	double direction = L * hN.dot(pA.cross(pB));
	left_handed_ = 1;
	if ( direction<0 ) left_handed_ = -1;

	double radius=pA.norm();

	radius_ = radius;

	TR.Debug << "align seg 4" << std::endl;

	// make residue list, needed for rotational transformation
	std::list <core::Size> residue_list;
	for ( Size ires=1; ires<= pose.size(); ++ires ) {
		if ( !pose.residue(ires).is_protein() ) continue;
		residue_list.push_back(ires);
	}


	//////////////////////////////////
	//  align axis to Z
	//////////////////////////////////
	//3x1
	Vector3f cA2axis= R0 + hN*(hN.dot(c_A-R0));
	double u(hN(0)),v(hN(1)),w(hN(2));
	// check! u,v,w should not be parallel to Z-axis (0,0,1)

	//translation
	for ( int i=0; i < A.cols(); i++ ) A.col(i) -= cA2axis;
	for ( int i=0; i < B.cols(); i++ ) B.col(i) -= cA2axis;

	// also need to translate the full coords

	//elements needed for transformation
	xyzMatrix< core::Real > Rid; //identity matrix for translation
	xyzVector< core::Real > preT = xyzVector< core::Real >(cA2axis(0),cA2axis(1),cA2axis(2));
	xyzVector< core::Real > postT = xyzVector< core::Real >(0,0,0);
	identity_matrix( Rid );
	protocols::forge::methods::apply_transformation( pose, residue_list, Rid, preT, postT);

	TR.Debug << "align seg 5" << std::endl;

	Matrix3f Txz(3,3);
	Txz<< u/sqrt(u*u+v*v), v/sqrt(u*u+v*v),0,
		-v/sqrt(u*u+v*v),u/sqrt(u*u+v*v),0,
		0               ,0              ,1;

	Matrix3f Tz(3,3);
	Tz<< w/sqrt(u*u+v*v+w*w)            ,0          ,-sqrt(u*u+v*v)/sqrt(u*u+v*v+w*w),
		0                              ,1          ,0,
		sqrt(u*u+v*v)/sqrt(u*u+v*v+w*w),0          ,w/sqrt(u*u+v*v+w*w);

	Matrix3f rot2Z=Tz*Txz;

	//rotation to z-axis
	A = (rot2Z * A);
	B = (rot2Z * B);


	//figure out the new axis with the CA set of coordinates
	double a(A.row(0).mean()), b(A.row(1).mean());
	//double c(A.row(2).mean());

	Matrix3f Tx;
	Tx<< a/sqrt(a*a+b*b), b/sqrt(a*a+b*b),0,
		-b/sqrt(a*a+b*b), a/sqrt(a*a+b*b),0,
		0               ,0              ,1;

	//rotate the CA coord set to x-axis, needed for generating ideal repeats
	A = (Tx * A);
	B = (Tx * B);

	//combine the rotation matrices for one step action
	//Matrix3f fullTx=rot2Z*Tx;
	Matrix3f fullTx=Tx * rot2Z;

	xyzMatrix< core::Real > xyzFull;
	matrix3f_to_xyzMatrix(fullTx, xyzFull);

	// in this step no translation
	preT =  xyzVector< core::Real >(0,0,0);
	postT = xyzVector< core::Real >(0,0,0);
	protocols::forge::methods::apply_transformation( pose, residue_list, xyzFull, preT, postT);

	for ( Size i = 1; i <= seg_size ; i++ ) {
		core::Vector coord = pose.residue(i).xyz("CA");
		A(0,i-1) = coord[0];
		A(1,i-1) = coord[1];
		A(2,i-1) = coord[2];
	}

	/*  only need to update A, this B thing just for debugging.
	for (Size i = seg_size+1; i <= 2*seg_size ; i++){
	core::Vector coord = pose.residue(i).xyz("CA");
	B(0,i-1-seg_size) = coord[0];
	B(1,i-1-seg_size) = coord[1];
	B(2,i-1-seg_size) = coord[2];
	}
	*/

	if ( is_sym ) { //re-symmetrize
		ConstraintSetOP pose_cst_set( new ConstraintSet( *pose.constraint_set() ) );
		symmetry::SetupForSymmetryMover pre_mover1;
		pre_mover1.apply( pose );
		pose.constraint_set(pose_cst_set);
		pose.pdb_info()->obsolete(true);
	}

	TR.Debug << "align seg 6" << std::endl;

}

void
RemodelGlobalFrame::setup_helical_constraint(Pose & pose){
	using Eigen::MatrixXd;
	using namespace Eigen;
	using core::id::AtomID;
	using namespace numeric;
	using namespace core::scoring::constraints;
	using namespace core::pose;
	using namespace basic::options;
	using namespace core::chemical;


	align_segment(pose);

	// now all coords (A) are along z-axis and the first c_A is in (radius,0,0)
	//////////////////////////////////////////////////
	TR.Debug << "setup RGF cst 1" << std::endl;
	//cache the constraints currently in pose
	ConstraintSetOP input_pose_cst_set;
	//input_pose_cst_set = new ConstraintSet( *pose.constraint_set() );
	//native_cst_set = new ConstraintSet(*input_pose_cst_set);

	//duplicate the native_set and add the new ones to it.
	input_pose_cst_set = ConstraintSetOP( new ConstraintSet(*native_cst_set) );

	// Extract the first segment from pose for generating the ideal decoy

	utility::vector1< core::Size > residue_indices;
	for ( Size i = 1; i <= seg_size; ++i ) {
		residue_indices.push_back(i);
	}
	PoseOP singleton_pose( new Pose );
	core::chemical::ResidueTypeSetCOP type_set = pose.residue_type_set_for_pose();
	core::io::pose_from_pose( *singleton_pose, pose, *type_set, residue_indices);

	// However the constraint type maps the same atoms, so we actually need yet another pose, with two copies.
	// the strategy is to use the singleton to figure out the coordinates, and dump the xyz to the second segment.
	// these's going to be an issue with connectivity....

	// expand residue indices to double
	for ( Size i = seg_size+1; i <= 2*seg_size; ++i ) {
		residue_indices.push_back(i);
	}

	PoseOP double_pose( new Pose );
	core::io::pose_from_pose( *double_pose, pose, *type_set, residue_indices);

	TR.Debug << "setup RGF cst 2" << std::endl;
	// make residue list, needed for rotational transformation
	std::list <core::Size> single_residue_list;
	std::list <core::Size> double_residue_list;
	std::list <core::Size> fullpose_residue_list;

	for ( Size ires=1; ires<= singleton_pose->size(); ++ires ) {
		if ( !singleton_pose->residue(ires).is_protein() ) continue;
		single_residue_list.push_back(ires);
	}
	for ( Size ires=1; ires<= double_pose->size(); ++ires ) {
		if ( !double_pose->residue(ires).is_protein() ) continue;
		double_residue_list.push_back(ires);
	}
	for ( Size ires=1; ires<= pose.size(); ++ires ) {
		if ( !pose.residue(ires).is_protein() ) continue;
		fullpose_residue_list.push_back(ires);
	}

	// make atomID list, use CA for now.   if wanted to adjust the number of atoms involved, do it here.
	utility::vector1< AtomID > atms;

	// we want to setup constraints for two copies, so make lists spanning the first two segments.
	for ( Size i = 1; i <= 2*seg_size ; i++ ) {
		AtomID id = AtomID(double_pose->residue(i).atom_index("CA"), i);
		atms.push_back(id);
	}

	// generating new coordinates following pre-defined helical parameters
	// pre-defined helical parameters, inputs from outside
	double t_rise =option[OptionKeys::remodel::helical_rise];
	double t_radius =option[OptionKeys::remodel::helical_radius];
	double t_omega =option[OptionKeys::remodel::helical_omega];

	//if(left_handed_ == -1) t_omega *= -1;

	Matrix3f Rz;
	Rz<<cos(t_omega), -sin(t_omega),0,
		sin(t_omega),  cos(t_omega),0,
		0           ,  0           ,1;

	//elements needed for transformation
	xyzMatrix< core::Real > Rid; //identity matrix for translation
	xyzVector< core::Real > preT;
	xyzVector< core::Real > postT;
	identity_matrix( Rid );
	//translate
	preT = xyzVector< core::Real >(0,0,0);
	postT = xyzVector< core::Real >(t_radius-radius_,0,0);

	TR.Debug << "setup RGF cst 3" << std::endl;
	protocols::forge::methods::apply_transformation( *singleton_pose, single_residue_list, Rid, preT, postT );
	//also have to move the doulbe_pose to the new frame.
	protocols::forge::methods::apply_transformation( *double_pose, double_residue_list, Rid, preT, postT );
	//as well as the full frame
	//commented out now, so constraints are applied to structure still
	//along Z; otherwise when workign with symmetry, it causes serious
	//issues with translation
	//protocols::forge::methods::apply_transformation( pose, fullpose_residue_list, Rid, preT, postT );

	// package Rz
	xyzMatrix< core::Real > xyzFull;
	matrix3f_to_xyzMatrix(Rz, xyzFull);
	postT = xyzVector< core::Real >(0,0,t_rise);

	protocols::forge::methods::apply_transformation( *singleton_pose, single_residue_list, xyzFull, preT, postT );

	//update second segment in the double pose


	TR.Debug << "setup RGF cst 4" << std::endl;
	for ( Size i = 1; i <= seg_size; i++ ) {
		Size natoms = double_pose->residue(i+seg_size).natoms(); // singleton may have more terminal atoms
		for ( Size j = 1; j <= natoms ; j++ ) {
			double_pose->set_xyz( AtomID( j, seg_size+i /* offset res number */ ), singleton_pose->residue(i).xyz(j));
		}
	}

	TR.Debug << "setup RGF cst 5" << std::endl;
	input_pose_cst_set->add_constraint(ConstraintCOP( ConstraintOP( new protocols::constraints_additional::BindingSiteConstraint( atms, *double_pose) ) ));

	TR.Debug << "setup RGF cst 6" << std::endl;
	pose.constraint_set(input_pose_cst_set);
	TR.Debug << "setup RGF cst 7" << std::endl;

	//get_helical_params(pose);

}

void
RemodelGlobalFrame::setup_CM_helical_constraint(Pose & pose){
	using Eigen::MatrixXd;
	using namespace Eigen;
	using core::id::AtomID;
	using namespace core;
	using namespace numeric;
	using namespace core::scoring::constraints;
	using namespace core::pose;
	using namespace basic::options;
	using namespace core::chemical;

	//align_segment(pose);


	// now all coords (A) are along z-axis and the first c_A is in (radius,0,0)
	//////////////////////////////////////////////////
	TR.Debug << "setup RGF cst 1" << std::endl;
	//cache the constraints currently in pose
	ConstraintSetOP input_pose_cst_set;
	//input_pose_cst_set = new ConstraintSet( *pose.constraint_set() );
	//native_cst_set = new ConstraintSet(*input_pose_cst_set);

	//duplicate the native_set and add the new ones to it.
	input_pose_cst_set = ConstraintSetOP( new ConstraintSet(*native_cst_set) );

	Size repeat_number = option[OptionKeys::remodel::repeat_structure];

	MatrixXf COMs( 3, repeat_number );

	for ( Size i = 0; i < repeat_number ; i++ ) {
		xyzVector< Real > com1 = compute_center_of_mass( pose, i*seg_size+1, (i+1)*seg_size);
		COMs.col(i) << com1[0], com1[1], com1[2];
	}

	/*
	//3 x N, column is x, y and z
	COMs << com1[0], com2[0], com3[0], com4[0], com5[0], com6[0],
	com1[1], com2[1], com3[1], com4[1], com5[1], com6[1],
	com1[2], com2[2], com3[2], com4[2], com5[2], com6[2];
	*/
	double t_rise =option[OptionKeys::remodel::helical_rise];
	double t_radius =option[OptionKeys::remodel::helical_radius];
	double t_omega =option[OptionKeys::remodel::helical_omega];
	double sd = option[OptionKeys::remodel::COM_sd];
	double tolerance = option[OptionKeys::remodel::COM_tolerance];

	MatrixXf COM_target = ideal_COMs(t_rise , t_radius , t_omega, repeat_number);

	Matrix3f H = rot_mat(COMs,COM_target);
	//    cout <<"CM Rotation matrix:\n"<<H<<endl<<endl;

	Vector3f c_COMs(COMs.row(0).mean(), COMs.row(1).mean(), COMs.row(2).mean());
	Vector3f c_COM_target(COM_target.row(0).mean(), COM_target.row(1).mean(), COM_target.row(2).mean());
	xyzVector< Real >  preT (c_COMs(0), c_COMs(1), c_COMs(2));   // preT
	xyzVector< Real >  postT (c_COM_target(0), c_COM_target(1), c_COM_target(2));   // preT

	using namespace core::pose::symmetry;
	using namespace protocols;

	// build residue list for transformation
	std::list <core::Size> fullpose_residue_list;
	for ( Size ires=1; ires<= seg_size*repeat_number; ++ires ) {
		if ( !pose.residue(ires).is_protein() ) continue;
		fullpose_residue_list.push_back(ires);
	}

	// package H
	xyzMatrix< core::Real > xyzFull;
	matrix3f_to_xyzMatrix(H, xyzFull);

	protocols::forge::methods::apply_transformation( pose, fullpose_residue_list, xyzFull , preT, postT );

	/* debug
	// re-assess the RMSD after transformation
	com1 = compute_center_of_mass( pose, 1, seg_size);
	com2 = compute_center_of_mass( pose, seg_size+1, 2*seg_size);
	com3 = compute_center_of_mass( pose, 2*seg_size+1, 3*seg_size);
	com4 = compute_center_of_mass( pose, 3*seg_size+1, 4*seg_size);

	MatrixXf newCOMs(3,4);
	//3 x N, column is x, y and z
	newCOMs << com1[0], com2[0], com3[0], com4[0],
	com1[1], com2[1], com3[1], com4[1],
	com1[2], com2[2], com3[2], com4[2];

	//get_RMSD( newCOMs, COM_target);
	*/

	for ( Size i = 0; i < repeat_number ; i++ ) {

		//get coordinates for the moving target
		utility::vector1< id::AtomID > ID_1s; //pull on first segment
		for ( Size ires = 1; ires <= seg_size; ++ires ) {
			conformation::Residue const &rsd = pose.residue(ires + seg_size*i); // offset
			for ( Size atmno=1 ; atmno <= rsd.natoms(); ++atmno ) {
				//std::cout << ires << " " << atmno << std::endl;
				ID_1s.push_back( id::AtomID(atmno, ires+ seg_size*i) ); // offset
			}
		}

		// repackage
		Vector idealCM_1( COM_target(0,i), COM_target(1,i), COM_target(2,i));

		//constraint
		input_pose_cst_set->add_constraint( ConstraintCOP( ConstraintOP( new protocols::constraints_additional::COMCoordinateConstraint( ID_1s, idealCM_1, sd, tolerance ) ) ));
	}

	pose.constraint_set(input_pose_cst_set);


}

/*  CODE BASE FOR 4 COM superposition
void
RemodelGlobalFrame::setup_CM_helical_constraint(Pose & pose){
using Eigen::MatrixXd;
using namespace Eigen;
using core::id::AtomID;
using namespace core;
using namespace numeric;
using namespace core::scoring::constraints;
using namespace core::pose;
using namespace basic::options;
using namespace core::chemical;

//align_segment(pose);


// now all coords (A) are along z-axis and the first c_A is in (radius,0,0)
//////////////////////////////////////////////////
TR.Debug << "setup RGF cst 1" << std::endl;
//cache the constraints currently in pose
ConstraintSetOP input_pose_cst_set;
//input_pose_cst_set = new ConstraintSet( *pose.constraint_set() );
//native_cst_set = new ConstraintSet(*input_pose_cst_set);

//duplicate the native_set and add the new ones to it.
input_pose_cst_set = new ConstraintSet(*native_cst_set);

xyzVector< Real > com1 = compute_center_of_mass( pose, 1, seg_size);
xyzVector< Real > com2 = compute_center_of_mass( pose, seg_size+1, 2*seg_size);
xyzVector< Real > com3 = compute_center_of_mass( pose, 2*seg_size+1, 3*seg_size);
xyzVector< Real > com4 = compute_center_of_mass( pose, 3*seg_size+1, 4*seg_size);

MatrixXf COMs(3,4);
//3 x N, column is x, y and z
COMs << com1[0], com2[0], com3[0], com4[0],
com1[1], com2[1], com3[1], com4[1],
com1[2], com2[2], com3[2], com4[2];

double t_rise =option[OptionKeys::remodel::helical_rise];
double t_radius =option[OptionKeys::remodel::helical_radius];
double t_omega =option[OptionKeys::remodel::helical_omega];

MatrixXf COM_target = ideal_COMs(t_rise , t_radius , t_omega, 4);
//get_RMSD(COMs,COM_target);

Matrix3f H = rot_mat(COMs,COM_target);
//    cout <<"CM Rotation matrix:\n"<<H<<endl<<endl;

//repackage COM_target

Vector idealCM_1( COM_target(0,0), COM_target(1,0), COM_target(2,0));
Vector idealCM_2( COM_target(0,1), COM_target(1,1), COM_target(2,1));
Vector idealCM_3( COM_target(0,2), COM_target(1,2), COM_target(2,2));
Vector idealCM_4( COM_target(0,3), COM_target(1,3), COM_target(2,3));


Vector3f c_COMs(COMs.row(0).mean(), COMs.row(1).mean(), COMs.row(2).mean());
Vector3f c_COM_target(COM_target.row(0).mean(), COM_target.row(1).mean(), COM_target.row(2).mean());
xyzVector< Real >  preT (c_COMs(0), c_COMs(1), c_COMs(2));   // preT
xyzVector< Real >  postT (c_COM_target(0), c_COM_target(1), c_COM_target(2));   // preT


// build residue list for transformation
std::list <core::Size> fullpose_residue_list;
for ( Size ires=1; ires<= pose.size(); ++ires ) {
if ( !pose.residue(ires).is_protein() ) continue;
fullpose_residue_list.push_back(ires);
}

// package H
xyzMatrix< core::Real > xyzFull;
matrix3f_to_xyzMatrix(H, xyzFull);

protocols::forge::methods::apply_transformation( pose, fullpose_residue_list, xyzFull , preT, postT );

// re-assess the RMSD after transformation
com1 = compute_center_of_mass( pose, 1, seg_size);
com2 = compute_center_of_mass( pose, seg_size+1, 2*seg_size);
com3 = compute_center_of_mass( pose, 2*seg_size+1, 3*seg_size);
com4 = compute_center_of_mass( pose, 3*seg_size+1, 4*seg_size);

MatrixXf newCOMs(3,4);
//3 x N, column is x, y and z
newCOMs << com1[0], com2[0], com3[0], com4[0],
com1[1], com2[1], com3[1], com4[1],
com1[2], com2[2], com3[2], com4[2];

//get_RMSD( newCOMs, COM_target);

//get coordinates for the moving target
utility::vector1< id::AtomID > ID_1s; //pull on first segment
utility::vector1< id::AtomID > ID_2s; //pull on second segment
utility::vector1< id::AtomID > ID_3s; //pull on second segment
utility::vector1< id::AtomID > ID_4s; //pull on second segment
utility::vector1< id::AtomID > ID_5s; //pull on second segment
utility::vector1< id::AtomID > ID_6s; //pull on second segment
for( Size ires = 1; ires <= seg_size; ++ires){
conformation::Residue const &rsd = pose.residue(ires);
for(Size atmno=1 ; atmno <= rsd.natoms(); ++atmno ){
//std::cout << ires << " " << atmno << std::endl;
ID_1s.push_back( id::AtomID(atmno, ires) );
}
}
for( Size ires = 1; ires <= seg_size; ++ires){
conformation::Residue const &rsd = pose.residue(ires+seg_size); // offset here
for(Size atmno=1 ; atmno <= rsd.natoms(); ++atmno ){
//std::cout << ires << " " << atmno << std::endl;
ID_2s.push_back( id::AtomID(atmno, ires+seg_size) ); //offset here
}
}
for( Size ires = 1; ires <= seg_size; ++ires){
conformation::Residue const &rsd = pose.residue(ires+2*seg_size); // offset here
for(Size atmno=1 ; atmno <= rsd.natoms(); ++atmno ){
//std::cout << ires << " " << atmno << std::endl;
ID_3s.push_back( id::AtomID(atmno, ires+2*seg_size) ); //offset here
}
}
for( Size ires = 1; ires <= seg_size; ++ires){
conformation::Residue const &rsd = pose.residue(ires+3*seg_size); // offset here
for(Size atmno=1 ; atmno <= rsd.natoms(); ++atmno ){
//std::cout << ires << " " << atmno << std::endl;
ID_4s.push_back( id::AtomID(atmno, ires+3*seg_size) ); //offset here
}
}

//constraint
input_pose_cst_set->add_constraint( new protocols::constraints_additional::COMCoordinateConstraint( ID_1s, idealCM_1 ));
input_pose_cst_set->add_constraint( new protocols::constraints_additional::COMCoordinateConstraint( ID_2s, idealCM_2 ));
input_pose_cst_set->add_constraint( new protocols::constraints_additional::COMCoordinateConstraint( ID_3s, idealCM_3 ));
input_pose_cst_set->add_constraint( new protocols::constraints_additional::COMCoordinateConstraint( ID_4s, idealCM_4 ));

pose.constraint_set(input_pose_cst_set);


}
*/

void
RemodelGlobalFrame::restore_original_cst(Pose & pose){
	TR.Debug << "pre_restore" << std::endl;
	pose.constraint_set(native_cst_set);
	TR.Debug << "post_restore" << std::endl;
}

void
RemodelGlobalFrame::set_native_cst_set( ConstraintSet const & cst_set){
	native_cst_set = ConstraintSetOP( new ConstraintSet(cst_set) );
}

void
RemodelGlobalFrame::set_native_cst_set( Pose const & pose ){
	native_cst_set = ConstraintSetOP( new ConstraintSet( *pose.constraint_set() ) );
}

void
RemodelGlobalFrame::set_segment_size( Size segment_size ){
	seg_size = segment_size;
}

void
RemodelGlobalFrame::matrix3f_to_xyzMatrix(Eigen::Matrix3f const & Re, numeric::xyzMatrix< core::Real> & R){
	R.xx(Re(0,0));R.xy(Re(0,1));R.xz(Re(0,2));
	R.yx(Re(1,0));R.yy(Re(1,1));R.yz(Re(1,2));
	R.zx(Re(2,0));R.zy(Re(2,1));R.zz(Re(2,2));
}

void
RemodelGlobalFrame::identity_matrix( numeric::xyzMatrix< core::Real> & R ){
	R.xx(1);R.xy(0);R.xz(0);
	R.yx(0);R.yy(1);R.yz(0);
	R.zx(0);R.zy(0);R.zz(1);
}


std::string
RemodelGlobalFrame::get_name() const {
	return "RemodelGlobalFrame";
}


void RemodelGlobalFrame::scorefunction( ScoreFunctionOP const & sfxn) {
	score_fxn_ = sfxn->clone();
}


numeric::xyzVector< core::Real >
compute_center_of_mass( core::pose::Pose const &  pose, core::Size range_start, core::Size range_stop){

	using namespace core;

	Real sum_x = 0;
	Real sum_y = 0;
	Real sum_z = 0;

	Size span = range_stop-range_start+1;

	for ( Size i = range_start; i <= range_stop ; i++ ) {
		core::Vector coord = pose.residue(i).xyz("CA");
		sum_x = sum_x + coord[0];
		sum_y = sum_y + coord[1];
		sum_z = sum_z + coord[2];
	}

	//cen of mass
	numeric::xyzVector< core::Real > com(sum_x/span, sum_y/span, sum_z/span);
	return com;
}

Matrix3f rot_mat(MatrixXf &A,MatrixXf &B)
{
	// A,B is 3 x N

	Vector3f c_A(A.row(0).mean(), A.row(1).mean(), A.row(2).mean());
	Vector3f c_B(B.row(0).mean(), B.row(1).mean(), B.row(2).mean());
	MatrixXf x_A(A.rows(),A.cols());
	MatrixXf x_B(B.rows(),B.cols());

	for ( int i=0; i < A.cols(); i++ ) {
		x_A.col(i) = A.col(i) - c_A;
		x_B.col(i) = B.col(i) - c_B;
	}

	Matrix3f cov= (x_B * x_A.transpose()) / x_A.cols();
	JacobiSVD<MatrixXf> svd(cov, ComputeFullU | ComputeFullV);

	Matrix3f Rt=svd.matrixU() * svd.matrixV().transpose();
	Matrix3f R;
	R<< 1,0,0, 0,1,0, 0,0,Rt.determinant();

	Matrix3f H = svd.matrixU() * R * svd.matrixV().transpose();
	return H;
}

MatrixXf ideal_COMs(double rise, double radius, double omega, int unitn){

	Matrix3f Rz;  // rot_mat along Z
	Rz<< cos(omega),-sin(omega),0, sin(omega),cos(omega),0, 0,0,1;
	Vector3f RISE(0,0,rise);
	Vector3f COM1(radius,0,0);

	MatrixXf COMs(3,unitn);   // 3 x N

	COMs.col(0) = COM1;
	Vector3f COMi = COM1;
	for ( int i=1; i<unitn; i++ ) { // n-1
		COMi = Rz * COMi + RISE;
		COMs.col(i) = COMi;
	}
	return COMs;
}

double get_RMSD(MatrixXf &A,MatrixXf &B){

	double rmsd=0;
	for ( int i=0; i < A.cols(); i++ ) {
		double d2 = (A.col(i)-B.col(i)).squaredNorm();
		rmsd += d2;
	}
	rmsd /= A.cols();
	rmsd = sqrt (rmsd);
	//  cout<<"RGF rmsd:"<<rmsd<<endl;
	return rmsd;
}


} // remodel
} // forge
} // protocol
