// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/RemodelLoopMover.cc
/// @brief  Loop modeling protocol based on routines from Remodel and EpiGraft
///         packages in Rosetta++.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)
/// @author Possu Huang (possu@u.washington.edu)

// unit headers
#include <protocols/forge/remodel/RemodelGlobalFrame.hh>
// AUTO-REMOVED #include <protocols/forge/remodel/RemodelMover.hh>
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
#include <core/pose/Remarks.hh>
#include <core/pose/util.hh> // for pdbinfo
#include <core/id/AtomID.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/io/pdb/file_data.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/constraints_additional/BindingSiteConstraint.hh>

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
#include <math.h>

using namespace basic::options;

namespace protocols {
namespace forge{
namespace remodel{

// Tracer instance for this file
// Named after the original location of this code
static basic::Tracer TR( "protocols.forge.remodel.RemodelGlobalFrame" );

// RNG
//static numeric::random::RandomGenerator RG( 342342 ); // magic number, don't change


// @brief default constructor
RemodelGlobalFrame::RemodelGlobalFrame()
	: op_user_remodel_repeat_structure_(option[OptionKeys::remodel::repeat_structure].user()),
		op_remodel_repeat_structure_(option[OptionKeys::remodel::repeat_structure]),
		op_remodel_helical_rise_(option[OptionKeys::remodel::helical_rise]),
		op_remodel_helical_radius_(option[OptionKeys::remodel::helical_radius]),
		op_remodel_helical_omega_(option[OptionKeys::remodel::helical_omega])
{
// has to reinitialize state before apply
	//state_.clear();
}

/// @brief value constructor
RemodelGlobalFrame::RemodelGlobalFrame(RemodelData const & remodel_data,
																			 RemodelWorkingSet const & working_model,
																			 ScoreFunctionOP const & sfxn)
	: op_user_remodel_repeat_structure_(option[OptionKeys::remodel::repeat_structure].user()),
		op_remodel_repeat_structure_(option[OptionKeys::remodel::repeat_structure]),
		op_remodel_helical_rise_(option[OptionKeys::remodel::helical_rise]),
		op_remodel_helical_radius_(option[OptionKeys::remodel::helical_radius]),
		op_remodel_helical_omega_(option[OptionKeys::remodel::helical_omega])
{

  remodel_data_ = remodel_data;
	seg_size = remodel_data.blueprint.size();
	working_model_ = working_model;
	score_fxn_ = sfxn->clone();
	left_handed_ = 0;
/*
	if (op_user_remodel_repeat_structure_){
		Size repeatCount = op_remodel_repeat_structure_;
		for (Size rep = 0; rep < repeatCount ; rep++){
			for (std::set< core::Size >::iterator it = uup.begin(); it != uup.end(); ++it){
			//DEBUG
			//	std::cout << *it + remodel_data.blueprint.size()*rep << std::endl;
			//	std::cout << "manger size"  << working_model.manager.union_of_intervals_containing_undefined_positions().size() <<  std::endl;
			//	std::cout << *it  << std::endl;
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
RemodelGlobalFrame::~RemodelGlobalFrame(){}

/// @brief clone this object
protocols::moves::MoverOP
RemodelGlobalFrame::RemodelGlobalFrame::clone() const {
	return new RemodelGlobalFrame( *this );
}

/// @brief create this type of object
protocols::moves::MoverOP
RemodelGlobalFrame::RemodelGlobalFrame::fresh_instance() const {
	return new RemodelGlobalFrame();
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

	//Size numRes = pose.total_residue();  // unused ~Labonte

	// dumping information into PDB header
  core::pose::PDBInfoOP temp_pdbinfo( new core::pose::PDBInfo(pose,true));

  core::pose::RemarkInfo remark;
	//capture stream
	std::stringstream capture_stream;

	if (op_user_remodel_repeat_structure_){
	}
	else {
		TR << "only applicable in Repeat mode";
		exit(0);
	}

	// get the size of each repeating sector
	//Size seg_size = remodel_data_.blueprint.size();
	//TR << "debug: seg_size = " << seg_size << std::endl;

	MatrixXf A(3, seg_size);
	MatrixXf B(3, seg_size);

	for (Size i = 1; i <= seg_size ; i++){
		core::Vector coord = pose.residue(i).xyz("CA");
		A(0,i-1) = coord[0];
		A(1,i-1) = coord[1];
		A(2,i-1) = coord[2];

		//TR << "resA " << i << " has " << coord[0] << " " << coord[1]  << " " << coord[2] << std::endl;
	}
	for (Size i = seg_size+1; i <= seg_size*2 ; i++){
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

	for(int i=0;i<A.cols();i++){
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
	for (Size i = 0; i<=2; i++){
		scalar = N.col(i).norm();
		if (scalar > max_scalar){
			max_scalar = scalar;
			hN = N.col(i)/N.col(i).norm();
		}
	}

	double sin_omega = (H(1,0)-H(0,1)) / (2*hN(2));
  TR.Debug <<"sin_omega "<<sin_omega<<endl;
  if( sin_omega < 0) hN = -1 * hN;


	Vector3f t = c_B - H*c_A;
	double L = t.dot(hN) ;
	double rise=abs(L);

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
  if(direction > 0) left_handed_ = -1;
  capture_stream<<"left handed: "<<left_handed_<<endl;

	double radius=pA.norm();

	capture_stream <<"L: " << L << endl << "rise: "<<rise<<endl<<"radius: "<<radius<<endl<<"omega: "<<omega<<endl;

	std::cout << capture_stream.str()<< endl;

	while (capture_stream.good()){
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

TR.Debug << "align seg 1" << std::endl;

	//Size numRes = pose.total_residue();  // unused ~Labonte

	if (op_user_remodel_repeat_structure_){
	}
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
TR.Debug << "pose length" << pose.total_residue() << std::endl;
TR.Debug << "seg_size" << seg_size << std::endl;

	MatrixXf A(3,seg_size);
	MatrixXf B(3,seg_size);

	for (Size i = 1; i <= seg_size ; i++){
		core::Vector coord = pose.residue(i).xyz("CA");
		A(0,i-1) = coord[0];
		A(1,i-1) = coord[1];
		A(2,i-1) = coord[2];

	}
	for (Size i = seg_size+1; i <= seg_size*2 ; i++){
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

	for(int i=0;i<A.cols();i++){
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
	for (Size i = 0; i<=2; i++){
		scalar = N.col(i).norm();
		if (scalar > max_scalar){
			max_scalar = scalar;
			hN = N.col(i)/N.col(i).norm();
		}
	}

TR.Debug << "N: " << N << std::endl;
TR.Debug << "helical axis: " << hN << std::endl;

	double sin_omega = (H(1,0)-H(0,1)) / (2*hN(2));

TR.Debug<<"sin_omega "<<sin_omega<<endl;

  if( sin_omega < 0) hN = -1 * hN;

	Vector3f t = c_B - H*c_A;
	double L = t.dot(hN) ;
	//double rise=abs(L);  // unused ~Labonte


	Matrix3f Ncross;
	Ncross << 0,-1*hN(2),hN(1), hN(2),0,-1*hN(0), -1*hN(1),hN(0),0 ;
	Matrix3f R0t= (1-cos(omega))*I - sin(omega)*Ncross;
	Vector3f R0 = R0t.inverse() * (t-L*hN);

  Vector3f pA= (c_A-R0)-(hN*(hN.dot(c_A-R0)));
  Vector3f pB= (c_B-R0)-(hN*(hN.dot(c_B-R0)));

  double direction = L * hN.dot(pA.cross(pB));
  left_handed_ = 1;
  if(direction<0) left_handed_ = -1;

	double radius=pA.norm();

	radius_ = radius;

TR.Debug << "align seg 4" << std::endl;

	// make residue list, needed for rotational transformation
  std::list <core::Size> residue_list;
  for ( Size ires=1; ires<= pose.total_residue(); ++ires ) {
    if ( !pose.residue(ires).is_protein() ) continue;
    residue_list.push_back(ires);
  }


	//////////////////////////////////
	//  align axis to Z
	//////////////////////////////////
	//3x1
   Vector3f cA2axis= R0 + hN*(hN.dot(c_A-R0));
   double u,v,w;
   u=hN(0),v=hN(1),w=hN(2);
   // check! u,v,w should not be parallel to Z-axis (0,0,1)

	 //translation
	 for(int i=0; i < A.cols(); i++) A.col(i) -= cA2axis;
	 for(int i=0; i < B.cols(); i++) B.col(i) -= cA2axis;

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
   double a,b;
	 //double c;
   a=A.row(0).mean(),b=A.row(1).mean(); //,c=A.row(2).mean();

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

	for (Size i = 1; i <= seg_size ; i++){
		core::Vector coord = pose.residue(i).xyz("CA");
		A(0,i-1) = coord[0];
		A(1,i-1) = coord[1];
		A(2,i-1) = coord[2];
	}

/*	 only need to update A, this B thing just for debugging.
	for (Size i = seg_size+1; i <= 2*seg_size ; i++){
		core::Vector coord = pose.residue(i).xyz("CA");
		B(0,i-1-seg_size) = coord[0];
		B(1,i-1-seg_size) = coord[1];
		B(2,i-1-seg_size) = coord[2];
	}
*/

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
				input_pose_cst_set = new ConstraintSet(*native_cst_set);

				// Extract the first segment from pose for generating the ideal decoy

				utility::vector1< core::Size > residue_indices;
				for(Size i = 1; i <= seg_size; ++i){
    			residue_indices.push_back(i);
  			}
				PoseOP singleton_pose = new Pose;
				core::chemical::ResidueTypeSet const & typeSet = (pose.residue(1).residue_type_set());
				core::io::pdb::pose_from_pose( *singleton_pose, pose, typeSet, residue_indices);

				// However the constraint type maps the same atoms, so we actually need yet another pose, with two copies.
				// the strategy is to use the singleton to figure out the coordinates, and dump the xyz to the second segment.
				// these's going to be an issue with connectivity....

				// expand residue indices to double
				for(Size i = seg_size+1; i <= 2*seg_size; ++i){
    			residue_indices.push_back(i);
  			}

				PoseOP double_pose = new Pose;
				core::io::pdb::pose_from_pose( *double_pose, pose, typeSet, residue_indices);

TR.Debug << "setup RGF cst 2" << std::endl;
				// make residue list, needed for rotational transformation
				std::list <core::Size> single_residue_list;
				std::list <core::Size> double_residue_list;
				std::list <core::Size> fullpose_residue_list;

				for ( Size ires=1; ires<= singleton_pose->total_residue(); ++ires ) {
						if ( !singleton_pose->residue(ires).is_protein() ) continue;
						single_residue_list.push_back(ires);
				}
				for ( Size ires=1; ires<= double_pose->total_residue(); ++ires ) {
						if ( !double_pose->residue(ires).is_protein() ) continue;
						double_residue_list.push_back(ires);
				}
				for ( Size ires=1; ires<= pose.total_residue(); ++ires ) {
						if ( !pose.residue(ires).is_protein() ) continue;
						fullpose_residue_list.push_back(ires);
				}

				// make atomID list, use CA for now.   if wanted to adjust the number of atoms involved, do it here.
				utility::vector1< AtomID > atms;

				// we want to setup constraints for two copies, so make lists spanning the first two segments.
				for (Size i = 1; i <= 2*seg_size ; i++){
					AtomID id = AtomID(double_pose->residue(i).atom_index("CA"), i);
					atms.push_back(id);
				}

        // generating new coordinates following pre-defined helical parameters
        // pre-defined helical parameters, inputs from outside
        double t_rise = op_remodel_helical_rise_;
        double t_radius = op_remodel_helical_radius_;
        double t_omega = op_remodel_helical_omega_;

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
				protocols::forge::methods::apply_transformation( pose, fullpose_residue_list, Rid, preT, postT );

				// package Rz
				xyzMatrix< core::Real > xyzFull;
				matrix3f_to_xyzMatrix(Rz, xyzFull);
				postT = xyzVector< core::Real >(0,0,t_rise);

				protocols::forge::methods::apply_transformation( *singleton_pose, single_residue_list, xyzFull, preT, postT );

				//update second segment in the double pose


TR.Debug << "setup RGF cst 4" << std::endl;
				for (Size i = 1; i <= seg_size; i++){
					Size natoms = double_pose->residue(i+seg_size).natoms(); // singleton may have more terminal atoms
					for (Size j = 1; j <= natoms ; j++){
						double_pose->set_xyz( AtomID( j, seg_size+i /* offset res number */ ), singleton_pose->residue(i).xyz(j));
					}
				}

TR.Debug << "setup RGF cst 5" << std::endl;
			  input_pose_cst_set->add_constraint(new protocols::constraints_additional::BindingSiteConstraint( atms, *double_pose));

TR.Debug << "setup RGF cst 6" << std::endl;
				pose.constraint_set(input_pose_cst_set);
TR.Debug << "setup RGF cst 7" << std::endl;

				//get_helical_params(pose);

}


void
RemodelGlobalFrame::restore_original_cst(Pose & pose){
TR.Debug << "pre_restore" << std::endl;
	pose.constraint_set(native_cst_set);
TR.Debug << "post_restore" << std::endl;
}

void
RemodelGlobalFrame::set_native_cst_set( ConstraintSet const & cst_set){
	native_cst_set = new ConstraintSet(cst_set);
}

void
RemodelGlobalFrame::set_native_cst_set( Pose const & pose ){
	native_cst_set = new ConstraintSet( *pose.constraint_set() );
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




} // remodel
} // forge
} // protocol
