// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/RepeatParameterFilter.cc
/// @brief can filter a pose for the rise, run and omega.  Many of the functions have been stolen from RepeatGlobalFrame but made const to work with filter obejcts. 
/// @author TJ Brunette

//Unit Headers
#include <protocols/simple_filters/RepeatParameterFilter.hh>
#include <protocols/simple_filters/RepeatParameterFilterCreator.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//external
#include <Eigen/Dense>

//Project Headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

namespace protocols{
namespace simple_filters {
        
using namespace core;
using namespace core::scoring;
using core::Real;
        
static basic::Tracer TR( "protocols.simple_filters.RepeatParameterFilter" );
        
protocols::filters::FilterOP RepeatParameterFilterCreator::create_filter() const { return protocols::filters::FilterOP( new RepeatParameterFilter ); }
        
std::string RepeatParameterFilterCreator::keyname() const { return "RepeatParameter"; }
        
RepeatParameterFilter::RepeatParameterFilter() : Filter( "RepeatParameter" ) {}
        
RepeatParameterFilter::~RepeatParameterFilter() {}
    
bool RepeatParameterFilter::apply( core::pose::Pose const & pose ) const {
    if(!filter_)
        return(true);
    else{
        std::string handedness;
        Real rise;
        Real radius;
        Real omega;
        calculate_helical_parameters(pose, handedness, rise, radius, omega);
        if(param_type_=="handedness"){ //never filtering on handedness 
                return(true);
        }
        if(param_type_ == "rise"){
            if((rise >= min_) && (rise <= max_))
                return(true);
            else
                return(false);
        }
        if(param_type_ == "radius"){
            if((radius >= min_) && (radius <= max_))
                return(true);
            else
                return(false);
        }
        if(param_type_ == "omega"){
            if((omega >= min_) && (omega <= max_))
                return(true);
            else
                return(false);
        }
    }
    throw utility::excn::EXCN_RosettaScriptsOption("Execution should not have gotten here!");
    return(true);
}

    
core::Real RepeatParameterFilter::report_sm(const core::pose::Pose & pose ) const{
    std::string handedness;
    Real rise;
    Real radius;
    Real omega;
    calculate_helical_parameters(pose, handedness, rise, radius, omega);
    if(param_type_== "handedness"){
        if(handedness=="R"){
	            return(1.0); //right handed
				}
        else{
            return(-1.0); //left handed
				}
		}
    if(param_type_=="rise")
        return(rise);
    if(param_type_=="radius")
        return(radius);
    if(param_type_=="omega")
        return(omega);
    throw utility::excn::EXCN_RosettaScriptsOption("Execution should not have gotten here. ");
    return(true);
}
    
void RepeatParameterFilter::report( std::ostream & out,const core::pose::Pose & pose ) const{
    std::string handedness;
    Real rise;
    Real radius;
    Real omega;
    calculate_helical_parameters(pose, handedness, rise, radius, omega);
    out << "rise:" << rise << " radius:" << radius << " omega:" << omega << " handedness:" << handedness << std::endl;
    
    
    }
    
void RepeatParameterFilter::calculate_helical_parameters( core::pose::Pose const & pose, std::string & handedness, Real & rise_out, Real & radius_out, Real & omega_out) const {
    core::pose::PoseOP tmpPose = pose.clone();
		calculate_helical_parameters_helper(*tmpPose, handedness, rise_out, radius_out,omega_out);
}

void RepeatParameterFilter::apply_transformation(core::pose::Pose & mod_pose,std::list <core::Size> const & residue_list, numeric::xyzMatrix< core::Real > const & R, numeric::xyzVector< core::Real > const & preT, numeric::xyzVector< core::Real > const & postT) const {
    using namespace ObjexxFCL;
    // translate xx2 by COM and fill in the new ref_pose coordinates
		utility::vector1< core::id::AtomID > ids;
    utility::vector1< numeric::xyzVector<core::Real> > positions;
        
    for (std::list<core::Size>::const_iterator it = residue_list.begin(); it != residue_list.end(); ++it) {
        core::Size ires = *it;
        for ( core::Size iatom=1; iatom<= mod_pose.residue_type(ires).natoms(); ++iatom ) { // use residue_type to prevent internal coord update
            ids.push_back(core::id::AtomID(iatom,ires));
            positions.push_back(postT + (R*( mod_pose.xyz(core::id::AtomID(iatom,ires)) - preT )));
        }
    }
    mod_pose.batch_set_xyz(ids,positions);
}

void RepeatParameterFilter::matrix3f_to_xyzMatrix(Eigen::Matrix3f const & Re, numeric::xyzMatrix< core::Real> & R) const{
    R.xx(Re(0,0));R.xy(Re(0,1));R.xz(Re(0,2));
    R.yx(Re(1,0));R.yy(Re(1,1));R.yz(Re(1,2));
    R.zx(Re(2,0));R.zy(Re(2,1));R.zz(Re(2,2));
}
    
void RepeatParameterFilter::identity_matrix( numeric::xyzMatrix< core::Real> & R ) const{
        R.xx(1);R.xy(0);R.xz(0);
        R.yx(0);R.yy(1);R.yz(0);
        R.zx(0);R.zy(0);R.zz(1);
}
    
void RepeatParameterFilter::calculate_helical_parameters_helper( core::pose::Pose const & pose, std::string & handedness, Real & rise_out, Real & radius_out, Real & omega_out) const {
    using Eigen::MatrixXd;
    using namespace Eigen;
    using namespace std;
    //copy from RemodelGlobalFrame because A. lack of const correctness, #remodel is in protocols.b.5 and I wanted this code to live in simple filters which is in protocols.3
		Size seg_size = pose.total_residue()/numb_repeats_;
		Size startResOffset = seg_size*(startAtRepeat_-1);
    MatrixXf A(3, seg_size);
    MatrixXf B(3, seg_size);
		for (Size i = 1; i <= seg_size; i++){
				Size offset =i+startResOffset;
        core::Vector coord = pose.residue(offset).xyz("CA");
        A(0,i-1) = coord[0];
        A(1,i-1) = coord[1];
        A(2,i-1) = coord[2];

    }
    for (Size i = seg_size+1; i <= seg_size*2 ; i++){
				Size offset =i+startResOffset;
        core::Vector coord = pose.residue(offset).xyz("CA");
        B(0,i-1-seg_size) = coord[0];
        B(1,i-1-seg_size) = coord[1];
        B(2,i-1-seg_size) = coord[2];
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

		double acos_fix_tmp = (H.trace()-1)/2;
		if(acos_fix_tmp <= -1.0)
			acos_fix_tmp = -1.0;
		if(acos_fix_tmp >= 1.0)
			acos_fix_tmp = 1.0;
		double omega = acos(acos_fix_tmp);
    Matrix3f I = Matrix3f::Identity();
    Matrix3f N = 0.5*(H+H.transpose()) - cos(omega)*I;
        
    
		Vector3f hN(0,0,0);//initializiation for the testing server  
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
        
    Matrix3f Ncross;
    Ncross << 0,-1*hN(2),hN(1), hN(2),0,-1*hN(0), -1*hN(1),hN(0),0 ;
    Matrix3f R0t= (1-cos(omega))*I - sin(omega)*Ncross;
    Vector3f R0 = R0t.inverse() * (t-L*hN);
    TR.Debug << "R0" << std::endl << R0 << std::endl;
    Vector3f pA= (c_A-R0)-(hN*(hN.dot(c_A-R0)));
    Vector3f pB= (c_B-R0)-(hN*(hN.dot(c_B-R0)));
        
        
    double direction = L * hN.dot(pA.cross(pB));
    if(direction > 0)
        handedness = "R";
    else
        handedness = "L";
    double radius=pA.norm();
    rise_out = rise;
    radius_out = radius;
    omega_out = omega;
}
    
    
    
    
void RepeatParameterFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & ){
    if (!tag->hasOption("numb_repeats"))
        throw utility::excn::EXCN_RosettaScriptsOption("RepeatPrarameter filter requires the number of repeats be entered with numb_repeats tag");
    numb_repeats_ = tag->getOption<Size>("numb_repeats");
		startAtRepeat_ = tag->getOption<Size>("start_at_repeat",1);
		if(startAtRepeat_+1>numb_repeats_)
				throw utility::excn::EXCN_RosettaScriptsOption("start_at_repeat must be atleast 1 smaller then the number of repeats ");
		param_type_ = tag->getOption<std::string>("param_type");
    if(!(param_type_=="radius" || param_type_=="rise" || param_type_=="omega" || param_type_=="handedness"))
        throw utility::excn::EXCN_RosettaScriptsOption("RepeatPrarameter filter requires one of 4 param_type = radius,rise,omega or handedness");
    if(param_type_=="radius" || param_type_=="rise" || param_type_=="omega"){
        min_=-9999;
        max_=9999;
        if(tag->hasOption("min")){
            min_ = tag->getOption<Real>("min");
            filter_=true;
        }
        if(tag->hasOption("max")){
            max_ = tag->getOption<Real>("max");
            filter_=true;
        }
    }
		std::string output_param_type = param_type_+"_r";  //r stands for repeat. There are other flags
    set_user_defined_name( output_param_type);
}
}//simple_filters
}//protocols
