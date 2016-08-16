// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/brunette/evalRepeats
///
/// @brief  analyzes helix oritentation in repeat proteins

/// @usage:

/// @author TJ Brunette


// Utility Headers
#include <basic/Tracer.hh>
#include <numeric/conversions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

#include <numeric/random/random.hh>

// Core Headers
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.hh>

#include <core/id/NamedAtomID.hh>
#include <core/id/SequenceMapping.hh>

#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <core/types.hh>

#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/jumping/util.hh>

#include <protocols/hybridization/TMalign.hh>
#include <protocols/hybridization/util.hh>

#include <devel/init.hh>

#include <utility/vector1.hh>
//external
#include <Eigen/Dense>

//basic & utility
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <utility/io/ozstream.hh>
#include <iostream>
#include <utility/io/ozstream.hh>

#include <ObjexxFCL/format.hh>

using namespace ObjexxFCL::format;
using utility::vector1;
using core::Size;
using core::Real;

static THREAD_LOCAL basic::Tracer TR( "extractNativeRepeats" );

void get_change_in_distance(const core::pose::Pose& pose,vector1<Real> & running_average_distance, vector1<Real> & change_in_distance, vector1<Real> & magnitude_of_change){
    Size local_smoothing_residues = 4;
    Size region_smoothing_residues = 20; //I have no clue how to set this.
    //distance seems to work better if it's proportional to the protein size
    numeric::xyzVector< core::Real >  first_res_pt = pose.xyz(core::id::NamedAtomID("CA", 1));
    numeric::xyzVector< core::Real >  last_res_pt = pose.xyz(core::id::NamedAtomID("CA", pose.total_residue()));
    Real max_dist_from_center_of_mass_along_axis = 80;
    if(pose.total_residue()<150)
        max_dist_from_center_of_mass_along_axis = 30; //seems to get better results having the point set closer in. Still miss several cases.
    numeric::xyzVector< core::Real > center_of_mass =  get_center_of_mass(pose);
    numeric::xyzVector< core::Real > random_pt;
    vector1<Real> running_total_distance;
    vector1<Real> local_running_average_smoothed_distance;
    vector1<Real> region_running_average_smoothed_distance;
    for(Size ii = 1; ii <= 25; ++ii){
        random_pt.x(center_of_mass.x()+numeric::random::random_range(1,max_dist_from_center_of_mass_along_axis));
        random_pt.y(center_of_mass.y()+numeric::random::random_range(1,max_dist_from_center_of_mass_along_axis));
        random_pt.z(center_of_mass.z()+numeric::random::random_range(1,max_dist_from_center_of_mass_along_axis));
        for(Size jj=1; jj<=pose.total_residue(); ++jj){
            numeric::xyzVector< core::Real >  pose_pt = pose.xyz(core::id::NamedAtomID("CA", jj));
            Real tmp_dist = random_pt.distance(pose_pt);
            if(ii == 1)
                running_total_distance.push_back(tmp_dist);
            else
                running_total_distance[jj] = running_total_distance[jj]+tmp_dist;
        }
    }
    for(Size jj=1; jj<=pose.total_residue(); ++jj){
        running_average_distance.push_back(running_total_distance[jj]/25);
    }
    for(Size ii=1; ii<=pose.total_residue(); ++ii){
        Real local_tmpSum = 0;
        Real local_count = 0;
        Real region_tmpSum = 0;
        Real region_count = 0;
        for(int jj=(int)ii-(int)local_smoothing_residues/2; jj<=(int)ii+(int)local_smoothing_residues/2; jj++){
            if(jj>=1 && jj<=(int)pose.total_residue()){
                local_tmpSum += running_average_distance[jj];
                local_count += 1;
            }
        }
        local_running_average_smoothed_distance.push_back(local_tmpSum/local_count);
        for(int jj=(int)ii-(int)region_smoothing_residues/2; jj<=(int)ii+(int)region_smoothing_residues/2; jj++){
            if(jj>=1 && jj<=(int)pose.total_residue()){
                region_tmpSum += running_average_distance[jj];
                region_count += 1;
            }
        }
        region_running_average_smoothed_distance.push_back(region_tmpSum/region_count);
    }

    change_in_distance.push_back(0);
    for(Size ii=2; ii<=pose.total_residue(); ++ii){
        change_in_distance.push_back(local_running_average_smoothed_distance[ii]-local_running_average_smoothed_distance[ii-1]);
    }
    for(Size ii=1; ii<=pose.total_residue(); ++ii){
        magnitude_of_change.push_back(std::abs(local_running_average_smoothed_distance[ii]-region_running_average_smoothed_distance[ii]));
    }
}

vector1<Size> get_inflection_points(vector1<Real> change_in_distance,vector1<Real> change_in_magnitude){
   //gets the max inflection points within a range
    vector1<Size> inflection_points;
    vector1<Size> filtered_inflection_points;
    int NEAREST_INFLECTION_POINT = 7; //eliminate points which are too close
    //step 1 find inflection points
    inflection_points.push_back(1);//first and last residues are automatically
    for(Size ii=2; ii<=change_in_distance.size()-1; ++ii){
        if((change_in_distance[ii-1]>=0 && change_in_distance[ii]<=0)||(change_in_distance[ii-1]<=0 && change_in_distance[ii]>=0))
            inflection_points.push_back(ii);
    }
    inflection_points.push_back(change_in_distance.size());//last residue
    filtered_inflection_points.push_back(1);//first residue
    for(Size ii=2; ii<=inflection_points.size()-1; ++ii){
        bool max_pt = true;
        for(Size kk=2; kk<=inflection_points.size()-1; ++kk){
            if(ii!=kk){
                int res1 = (int)inflection_points[ii];
                int res2 = (int)inflection_points[kk];
                if(std::abs(res2-res1) <= NEAREST_INFLECTION_POINT){//only worry about points close together
                    if(change_in_magnitude[res2] >= change_in_magnitude[res1])
                        max_pt = false;
                }
            }
        }
        if(max_pt == true){
            filtered_inflection_points.push_back(inflection_points[ii]);
        }
    }
    filtered_inflection_points.push_back(change_in_distance.size());//last residue
    return(filtered_inflection_points);
}

Size get_nearest_loop_to_helix(core::pose::Pose& pose, Size starting_pt){
    Size max_dist = numeric::max(starting_pt-1,pose.total_residue()-starting_pt);
    for(int ii=0; ii<(int)max_dist; ++ii){
        if((pose.secstruct(starting_pt-ii-1) == 'L') && (pose.secstruct(starting_pt-ii) == 'H'))
            return(starting_pt-ii);
        if((pose.secstruct(starting_pt+ii) == 'L') && (pose.secstruct(starting_pt+ii+1) == 'H'))
            return(starting_pt+ii+1);
    }
    return(pose.total_residue());
}


vector1<Size> get_loop_to_helix_points(core::pose::Pose& pose,vector1<Size> inflection_points){
    vector1<Size> loop_to_helix_points;
    for(Size ii=1; ii<=inflection_points.size(); ++ii){
        Size tmp_loop_to_helix = get_nearest_loop_to_helix(pose,inflection_points[ii]);
        if(std::find(loop_to_helix_points.begin(), loop_to_helix_points.end(), tmp_loop_to_helix)==loop_to_helix_points.end()){
            loop_to_helix_points.push_back(tmp_loop_to_helix);
        }
    }
    return(loop_to_helix_points);
}


Size find_next_repeat_inflection_point(const core::pose::Pose& pose, vector1<Size> inflection_points, Size starting_inflection_pt,Size MIN_REPEAT_SIZE){
    Real min_distance_above_threshold = 999;
    Size position = 999;
    for(Size ii=starting_inflection_pt+1; ii<=inflection_points.size(); ++ii){
        Size res1 = inflection_points[starting_inflection_pt];
        Size res2 = inflection_points[ii];
        if(res2-res1>MIN_REPEAT_SIZE){
            numeric::xyzVector< core::Real >  res1_xyz = pose.xyz(core::id::NamedAtomID("CA", res1));
            numeric::xyzVector< core::Real >  res2_xyz = pose.xyz(core::id::NamedAtomID("CA", res2));
            Real tmp_distance = res1_xyz.distance(res2_xyz);
            if(tmp_distance<min_distance_above_threshold){
                min_distance_above_threshold=tmp_distance;
                position = ii;
            }
        }
    }
    return(position);
}


vector1<Size> gather_repeat(const core::pose::Pose& pose, vector1<Size> inflection_points,Size max_res_variance,Size starting_inflection_pt,Size numb_repeats){
    Size MIN_REPEAT_SIZE = 15;
    //step1 get potential repeat
    bool failed=false;
    Size last_inflection_pt = starting_inflection_pt;
    vector1 <Size> repeat_positions;
    repeat_positions.push_back(inflection_points[starting_inflection_pt]);
    for(Size ii=1; ii<=numb_repeats && !failed; ++ii){
        Size next_pt = find_next_repeat_inflection_point(pose,inflection_points,last_inflection_pt,MIN_REPEAT_SIZE);
        if(next_pt == 999)
            failed = true;
        else{
            last_inflection_pt = next_pt;
            repeat_positions.push_back(inflection_points[next_pt]);
        }
    }
    if(repeat_positions.size() != numb_repeats+1){
        vector1 <Size> repeat_positions_to_return;
        repeat_positions_to_return.push_back(99999);
        return(repeat_positions_to_return);
    }
    //step3 check that the repeats unit length is < max_res_variance
    //step3a get avg length between repeats
    Size total_dist = 0;
    for(Size ii = 2; ii <= repeat_positions.size(); ++ii)
        total_dist += repeat_positions[ii]-repeat_positions[ii-1];
    Real avg_dist = (Real)total_dist/(Real)(repeat_positions.size()-1);
    //step3b see if any diverge from the repeat by more then max_variance
    bool within_res_variance = true;
    for(Size ii = 2; ii <= repeat_positions.size(); ++ii){
        Real tmp_dist = repeat_positions[ii]-repeat_positions[ii-1];
        if((tmp_dist > avg_dist+max_res_variance) || (tmp_dist < avg_dist-max_res_variance))
            within_res_variance = false;
    }
    if(!within_res_variance){
        vector1 <Size> repeat_positions_to_return;
        repeat_positions_to_return.push_back(99999);
        return(repeat_positions_to_return);
    }
    return(repeat_positions);
}


void matrix3f_to_xyzMatrix(Eigen::Matrix3f const & Re, numeric::xyzMatrix< core::Real> & R){
    R.xx(Re(0,0));R.xy(Re(0,1));R.xz(Re(0,2));
    R.yx(Re(1,0));R.yy(Re(1,1));R.yz(Re(1,2));
    R.zx(Re(2,0));R.zy(Re(2,1));R.zz(Re(2,2));
}

void identity_matrix( numeric::xyzMatrix< core::Real> & R ){
        R.xx(1);R.xy(0);R.xz(0);
        R.yx(0);R.yy(1);R.yz(0);
        R.zx(0);R.zy(0);R.zz(1);
}

void calculate_helical_parameters( core::pose::Pose const & pose, Size startRep1, Size endRep1, Size startRep2, Size endRep2, std::string & handedness, Real & rise_out, Real & radius_out, Real & omega_out){
    using Eigen::MatrixXd;
    using namespace Eigen;
    using namespace std;
    //copy from RemodelGlobalFrame because A. lack of const correctness, #remodel is in protocols.b.5 and I wanted this code to live in simple filters which is in protocols.3
	Size seg1_size = endRep1-startRep1+1;
    Size seg2_size = endRep2-startRep2+1;
    Size max_size;
    if(seg1_size > seg2_size)
        max_size = seg1_size;
    else
        max_size = seg2_size;
    MatrixXf A(3, max_size);
    MatrixXf B(3, max_size);
    for (Size i = 1; i <= max_size ; i++){
    if(i+startRep1-1<=endRep1){
            core::Vector coord = pose.residue(i+startRep1-1).xyz("CA");
            A(0,i-1) = coord[0];
            A(1,i-1) = coord[1];
            A(2,i-1) = coord[2];
        }
        else{
            A(0,i-1) = 0;
            A(1,i-1) = 0;
            A(2,i-1) = 0;
        }
    }
    for (Size i = 1; i <= max_size ; i++){
         if(i+startRep2-1<=endRep2){
            core::Vector coord = pose.residue(i+startRep2-1).xyz("CA");
            B(0,i-1) = coord[0];
            B(1,i-1) = coord[1];
            B(2,i-1) = coord[2];
         }
         else{
            B(0,i-1) = 0;
            B(1,i-1) = 0;
            B(2,i-1) = 0;
        }
    }
    Vector3f c_A(A.row(0).mean(), A.row(1).mean(), A.row(2).mean());
    Vector3f c_B(B.row(0).mean(), B.row(1).mean(), B.row(2).mean());
    MatrixXf x_A(A.rows(),A.cols());
    MatrixXf x_B(B.rows(),B.cols());
    for(int i=0;i<A.cols();i++){
        x_A.col(i)=A.col(i)-c_A;
    }
     for(int i=0;i<B.cols();i++){
        x_B.col(i)=B.col(i)-c_B;
    }
    Matrix3f cov= (x_B * x_A.transpose()) / x_A.cols();
    JacobiSVD<MatrixXf> svd(cov, ComputeFullU | ComputeFullV);
    //JacobiSVD<MatrixXf> svd(cov, ComputeThinU | ComputeThinV);
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

Real tm_superimpose(core::pose::Pose & pose, core::pose::Pose & ref_pose){
    using namespace core::id;
    using namespace protocols::hybridization;
    std::list <Size> pose_residue_list; // all residue numbers in ipose, used for transformation after alignment
    for (Size ires = 1; ires <= pose.total_residue(); ++ires) {
        pose_residue_list.push_back(ires);
    }
    std::list <Size> ref_pose_residue_list; // all residue numbers in ipose, used for transformation after alignment
    for (Size ires = 1; ires <= ref_pose.total_residue(); ++ires) {
        ref_pose_residue_list.push_back(ires);
    }
    TMalign tm_align;
    tm_align.apply(pose,ref_pose,pose_residue_list,ref_pose_residue_list);
    core::id::AtomID_Map< core::id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose, core::id::BOGUS_ATOM_ID );
    core::Size n_mapped_residues=0;
	tm_align.alignment2AtomMap(pose, ref_pose, pose_residue_list,ref_pose_residue_list, n_mapped_residues, atom_map);
    core::Size normalize_length = pose.total_residue() < ref_pose.total_residue() ? pose.total_residue() : ref_pose.total_residue();
    core::Real TMscore = tm_align.TMscore(normalize_length);
    utility::vector1< core::Real > aln_cutoffs;
	aln_cutoffs.push_back(2);
	aln_cutoffs.push_back(1.5);
	aln_cutoffs.push_back(1.0);
	aln_cutoffs.push_back(0.5);
	core::Real min_coverage = 0.2;
    partial_align(pose,ref_pose, atom_map, pose_residue_list, true, aln_cutoffs, min_coverage);
    return(TMscore);
}

void calculate_helical_parameters_helper(core::pose::Pose const & pose, Size startRep1, Size endRep1, Size startRep2, Size endRep2, std::string & handedness, Real & rise_out, Real & radius_out, Real & omega_out, Real & tm_out ){
    using namespace protocols;
    using namespace core::id;
    //Create new pose with repeat 1 superimposed onto repeat 2.
    //step1: cutout repeat 1
    vector1<Size> slice_res1;
    vector1<Size> slice_res2;
    core::pose::Pose repeat1;
    core::pose::Pose repeat2;
    for(Size ii=startRep1; ii<=endRep1; ++ii)
        slice_res1.push_back(ii);
    pdbslice(repeat1,pose,slice_res1);
    for(Size ii=startRep2; ii<=endRep2; ++ii)
        slice_res2.push_back(ii);
    pdbslice(repeat2,pose,slice_res2);
    core::pose::Pose repeat1_clone = *repeat1.clone();
    Size max_res = numeric::min(repeat1.total_residue(), repeat2.total_residue());
   /* core::id::AtomID_Map< core::id::AtomID > atom_map;
    core::pose::initialize_atomid_map( atom_map, repeat1, BOGUS_ATOM_ID );
   	for (Size ii=1; ii<=max_res; ++ii){
        core::id::AtomID const id1(repeat1_clone.residue(ii).atom_index("CA"),ii );
        core::id::AtomID const id2(repeat2.residue(ii).atom_index("CA"),ii);
        atom_map[id1]=id2;
    }
    core::scoring::superimpose_pose(repeat1_clone,repeat2,atom_map);
    */
    tm_out = tm_superimpose(repeat1_clone,repeat2);
    append_pose_to_pose(repeat1,repeat1_clone,false);
    Size res_per_repeat = repeat1.total_residue()/2;
    calculate_helical_parameters(repeat1,1,res_per_repeat,res_per_repeat+1,repeat1.total_residue(),handedness,rise_out,radius_out, omega_out);
}

Real calculate_sheet_pct(core::pose::Pose& pose,vector1<Size> repeat_positions){
    Size total_sheet = 0;
    Size total_residue = 0;
    for(Size ii=repeat_positions[0]; ii < repeat_positions[repeat_positions.size()-1]; ++ii){
        if(pose.secstruct(ii) == 'E')
            total_sheet++;
        total_residue++;
    }
    return((Real)total_sheet/(Real)total_residue);
}

Real calculate_rmsd_variance(core::pose::Pose & pose,Size start1_res, Size end1_res, Size start2_res, Size end2_res){
    using namespace core::sequence;
	using namespace core::scoring;
	using namespace core::id;
	vector1<Size> mod_pose_positions;
	vector1<Size> ref_pose_positions;
	for (Size ii=start1_res; ii<= end1_res; ++ii){
		mod_pose_positions.push_back(ii);
	}
    for (Size ii=start2_res; ii<= end2_res; ++ii){
		ref_pose_positions.push_back(ii);
	}
	core::kinematics::FoldTree f_mod(mod_pose_positions.size());
	core::kinematics::FoldTree f_ref(ref_pose_positions.size());
    core::pose::Pose mod_pose;
    core::pose::Pose ref_pose;
	core::pose::create_subpose(pose,mod_pose_positions,f_mod,mod_pose);
	core::pose::create_subpose(pose,ref_pose_positions,f_ref,ref_pose);
	return(core::scoring::CA_rmsd(mod_pose,ref_pose));
}

Real calculate_tm_variance(core::pose::Pose & pose,Size start1_res, Size end1_res, Size start2_res, Size end2_res){
    using namespace core::sequence;
	using namespace core::scoring;
	using namespace core::id;
	vector1<Size> mod_pose_positions;
	vector1<Size> ref_pose_positions;
	for (Size ii=start1_res; ii<= end1_res; ++ii){
		mod_pose_positions.push_back(ii);
	}
    for (Size ii=start2_res; ii<= end2_res; ++ii){
		ref_pose_positions.push_back(ii);
	}
	core::kinematics::FoldTree f_mod(mod_pose_positions.size());
	core::kinematics::FoldTree f_ref(ref_pose_positions.size());
    core::pose::Pose mod_pose;
    core::pose::Pose ref_pose;
	core::pose::create_subpose(pose,mod_pose_positions,f_mod,mod_pose);
	core::pose::create_subpose(pose,ref_pose_positions,f_ref,ref_pose);
    Real tm_score = tm_superimpose(mod_pose,ref_pose);
	return(tm_score);
}


Real calculate_max_rmsd_variance(core::pose::Pose & pose, vector1<Size> repeat_positions){
    Real max_variance = 0;
    for(Size ii =1; ii<=repeat_positions.size()-2; ++ii){
        Real tmp_variance = calculate_rmsd_variance(pose,repeat_positions[ii],repeat_positions[ii+1]-1,repeat_positions[ii+1],repeat_positions[ii+2]-1);
        if(tmp_variance > max_variance){
            max_variance = tmp_variance;
        }
    }
    return(max_variance);
}


Real calculate_min_tm_variance(core::pose::Pose & pose, vector1<Size> repeat_positions){
    Real min_variance = 999;
    for(Size ii =1; ii<=repeat_positions.size()-2; ++ii){
        Real tmp_variance = calculate_tm_variance(pose,repeat_positions[ii],repeat_positions[ii+1]-1,repeat_positions[ii+1],repeat_positions[ii+2]-1);
        if(tmp_variance < min_variance){
            min_variance = tmp_variance;
        }
    }
    return(min_variance);
}


vector1 <vector1 <Size> >  gather_repeats(core::pose::Pose& pose, vector1<Size> inflection_points){
    Size max_res_variance = 3;
    Size numb_repeats = 3;
    Size current_starting_inflection_pt = 1;
    vector1 <vector1 <Size> > all_repeats;
    while(current_starting_inflection_pt <= (inflection_points.size()-2*numb_repeats)){
        vector1 <Size> repeat_positions = gather_repeat(pose,inflection_points,max_res_variance,current_starting_inflection_pt,numb_repeats);
        Real sheet_pct = calculate_sheet_pct(pose,repeat_positions);
        if(repeat_positions[1] == 99999||sheet_pct > 0){
            current_starting_inflection_pt += 1;
        }
        else{
            Real sheet_pct = calculate_sheet_pct(pose,repeat_positions);
            Real tm_variance = calculate_min_tm_variance(pose,repeat_positions);
            Real rmsd_variance = calculate_max_rmsd_variance(pose,repeat_positions);
            std::cout << "rmsd" << rmsd_variance << std::endl;
            std::cout << "tm" << tm_variance << std::endl;
            if((sheet_pct == 0) && (tm_variance>0.6)){
                std::cout << "tm_taken" << tm_variance << std::endl;
                all_repeats.push_back(repeat_positions);
            }
            current_starting_inflection_pt += 1; //may want to jump to put more complicated logic in now to skip intermediate repeats
          }
    }
    return(all_repeats);
}

void gather_info_all_repeats(core::pose::Pose& pose,vector1 <vector1 <Size> > all_repeats,vector1 <vector1 <Size> > & selected_repeats, vector1<std::string> & handednesses_v, vector1<Real> & radius_v, vector1<Real> & rise_v, vector1<Real> & omega_v,vector1<Real> & rmsd_v,vector1<Real> & tm_v){
    //just average rise,run,omega rather then looking for canonical repeat.
    std::string handedness;
    Real rise;
    Real radius;
    Real omega;
    Real tm;
    vector1<Size> repeatStarts;
    for(Size ii=1; ii<=all_repeats.size(); ++ii){
        Real tmp_rmsd = calculate_max_rmsd_variance(pose,all_repeats[ii]);
        for(Size kk=1; kk<=all_repeats[ii].size()-2; ++kk){
            if(std::find(repeatStarts.begin(), repeatStarts.end(), all_repeats[ii][kk])==repeatStarts.end()){
                repeatStarts.push_back(all_repeats[ii][kk]);
                selected_repeats.push_back(all_repeats[ii]);
                calculate_helical_parameters_helper(pose,all_repeats[ii][kk], all_repeats[ii][kk+1]-1, all_repeats[ii][kk+1], all_repeats[ii][kk+2]-1, handedness, rise, radius, omega,tm);
                handednesses_v.push_back(handedness);
                radius_v.push_back(radius);
                std::cout << "radius" << radius << ",start:" << all_repeats[ii][kk] << std::endl;
                rise_v.push_back(rise);
                omega_v.push_back(omega);
                tm_v.push_back(tm);
                rmsd_v.push_back(tmp_rmsd);
            }
        }
    }
}

void gather_info_canonical_repeats(core::pose::Pose& pose,vector1 <vector1 <Size> > all_repeats,vector1 <vector1 <Size> > & selected_repeats, vector1<std::string> & handednesses_v, vector1<Real> & radius_v, vector1<Real> & rise_v, vector1<Real> & omega_v,vector1<Real> & rmsd_v,vector1<Real> & tm_v){
    //just average rise,run,omega rather then looking for canonical repeat.
    std::string handedness;
    Real rise;
    Real radius;
    Real omega;
    Real tm;
    Real rise_total = 0;
    Real radius_total = 0;
    Real omega_total = 0;
    Real tm_total;
    Real rmsd_total;
    Size angle_ct = 0;
    using std::map;
    Size canonical_repeat = 0;
    Real max_tm_variance = 0;
    for(Size ii=1; ii<=all_repeats.size(); ++ii){
        Real tmp_tm_variance = calculate_min_tm_variance(pose, all_repeats[ii]);
        if(tmp_tm_variance > max_tm_variance){
            max_tm_variance = tmp_tm_variance;
            canonical_repeat = ii;
        }
    }
    selected_repeats.push_back(all_repeats[canonical_repeat]);
    for(Size kk=1; kk<=all_repeats[canonical_repeat].size()-2; ++kk){
        calculate_helical_parameters_helper(pose,all_repeats[canonical_repeat][kk], all_repeats[canonical_repeat][kk+1]-1, all_repeats[canonical_repeat][kk+1], all_repeats[canonical_repeat][kk+2]-1, handedness, rise, radius, omega,tm);
        rise_total += rise;
        radius_total += radius;
        omega_total += omega;
        tm_total += tm;
        rmsd_total += calculate_max_rmsd_variance(pose,all_repeats[canonical_repeat]);
        angle_ct += 1;
    }
    /* From when I was doing averaging.
    selected_repeats.push_back(all_repeats[canonical_repeat]);
    for(Size ii=1; ii<=all_repeats.size(); ++ii){
        for(Size kk=1; kk<=all_repeats[ii].size()-2; ++kk){
            calculate_helical_parameters_helper(pose,all_repeats[ii][kk], all_repeats[ii][kk+1]-1, all_repeats[ii][kk+1], all_repeats[ii][kk+2]-1, handedness, rise, radius, omega,tm);
            rise_total += rise;
            radius_total += radius;
            omega_total += omega;
            std::cout << "TM:" << tm << std::endl;
            tm_total += tm;
            rmsd_total += calculate_max_rmsd_variance(pose,all_repeats[ii]);
            angle_ct += 1;
        }
    }*/

    handednesses_v.push_back(handedness);
    radius_v.push_back(radius_total/(Real)angle_ct);
    rise_v.push_back(rise_total/(Real)angle_ct);
    omega_v.push_back(omega_total/(Real)angle_ct);
    tm_v.push_back(tm_total/(Real)angle_ct);
    rmsd_v.push_back(rmsd_total/(Real)angle_ct);
}

void output_all_repeats(vector1 <vector1 <Size> > all_repeats, std::string id,  vector1<std::string> & handedness, vector1<Real> & radii, vector1<Real> & rises, vector1<Real> & omegas, vector1<Real> rmsd, vector1<Real> tm, std::ostream& output){
    for(Size ii=1; ii<= all_repeats.size(); ++ii){
        output << id << " ";
        for(Size kk=1; kk<=all_repeats[ii].size(); ++kk){
            output << all_repeats[ii][kk] << ",";
        }
        output << rmsd[ii] <<",";
        output << tm[ii] << ",";
        output << handedness[ii] <<",";
        output << radii[ii] << ",";
        output << rises[ii] << ",";
        output << omegas[ii];
        output << std::endl;
    }
}


int main( int argc, char * argv [] ) {
    try {
    using namespace core::chemical;
    using namespace core::import_pose::pose_stream;
    using core::import_pose::pose_from_file;
    using namespace core::scoring;
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    devel::init(argc, argv);
    //create vector of input poses.
    MetaPoseInputStream input = streams_from_cmd_line();
    core::pose::Pose current_pose;
   utility::io::ozstream all_out("all_repeats_out.txt");
   utility::io::ozstream canonical_out("canonical_repeats.txt");
   utility::io::ozstream change("dist_change.txt");
   utility::io::ozstream inflect("inflection_pts.txt");
    while(input.has_another_pose()){
        input.fill_pose(current_pose,*rsd_set_from_cmd_line().lock());
        std::string tag = core::pose::tag_from_pose(current_pose);
        protocols::jumping::assign_ss_dssp(current_pose);
        std::cout << "**********************working on " << tag << std::endl;
        vector1<Real> change_in_distance;
        vector1<Real> magnitude_of_change;
        vector1<Real> running_average_distance;
        get_change_in_distance(current_pose,running_average_distance,change_in_distance,magnitude_of_change);
        for(Size ii=1; ii<=running_average_distance.size(); ++ii)
           change <<running_average_distance[ii] << std::endl;
        vector1<Size> inflection_points = get_inflection_points(change_in_distance,magnitude_of_change);
        vector1<Size> loop_to_helix_points = get_loop_to_helix_points(current_pose,inflection_points);
        for(Size ii=1; ii<=inflection_points.size(); ++ii)
            inflect << inflection_points[ii] << std::endl;
        inflect << "---------------------------" << std::endl;
        std::cout << loop_to_helix_points.size() << std::endl;
        for(Size ii=1; ii<=loop_to_helix_points.size(); ++ii)
            inflect << loop_to_helix_points[ii] << std::endl;
        vector1<vector1<Size> > all_repeats = gather_repeats(current_pose,loop_to_helix_points);
        vector1<std::string> handedness;
        vector1<Real> radius;
        vector1<Real> rise;
        vector1<Real> omega;
        vector1<Real> rmsd;
        vector1<Real> tm;
        vector1<vector1<Size> > output_repeats;
        if(all_repeats.size() > 0){
            gather_info_canonical_repeats(current_pose,all_repeats,output_repeats,handedness,radius,rise,omega,rmsd,tm);
            output_all_repeats(output_repeats,tag,handedness,radius,rise,omega,rmsd,tm,canonical_out);
            handedness.clear();
            radius.clear();
            rise.clear();
            omega.clear();
            rmsd.clear();
            tm.clear();
            output_repeats.clear();
            gather_info_all_repeats(current_pose,all_repeats,output_repeats,handedness,radius,rise,omega,rmsd,tm);
            output_all_repeats(output_repeats,tag,handedness,radius,rise,omega,rmsd,tm,all_out);
        }
    }
    }catch ( utility::excn::EXCN_Base const & e ) {
        std::cout << "caught exception " << e.msg() << std::endl;
        return -1;
    }
return 0;
}

