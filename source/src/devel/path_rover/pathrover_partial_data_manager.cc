
////////////////////////////////////////////////////////////////////////////////
//     pathways_partial_data_manager.cc: module for handling partial data info
//
//     see pathways_partial_data_manager.h
/////////////////////////////////////////////////////////////////////////////////


// Rosetta Headers
#include "pack.h"
#include "pathways.h"
#include "pathways_planners.h"
#include "pose.h"
#include "pose_io.h"
#include "pose_rotamer_trials.h"
#include "prof.h"
#include "random_numbers.h"
#include "rotamer_trials.h"
#include "score.h"
#include "score_ns.h"

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// Utility Headers
//#include "src/utility/basic_sys_util.h"
//#include "src/utility/io/ozstream.h"

// C++ Headers
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <string.h>
#include <string>
#include <sstream>
#include <math.h>
#include "pathways_partial_data_manager.h"
#include "pathways_DOFs_manager.h"


namespace pathways{

	bool by_id(PD_match p,PD_match q) { return (p.id < q.id);}
  //bool by_id(PD_rmsd_match p,PD_rmsd_match q) { return (p.id < q.id);}

 double max2(double x, double y){
  return ((x > y)? x : y);
}
	double pythag(double a, double b)
{
  double absa, absb;

  absa = fabs(a);
  absb = fabs(b);
  if(absa > absb)
    return (absa * sqrt(1.0 + ((absb/absa)*(absb/absa))));
  else
    return (absb == 0.0 ? 0.0 : absb * sqrt(1.0+((absa/absb)*(absa/absb))));
}

	double torsion4(Vec3d v1,
                    Vec3d v2,
                    Vec3d v3,
                    Vec3d v4) {

    Vec3d n1,n2;
    double l1,l2,ct;
    double sign,angle;

    Vec3d v12= v1-v2;
    Vec3d v32= v3-v2;
    Vec3d v34= v3-v4;

	// normal to plane 1
    n1=v12.cross(v32);

    // normal to plane 2
    n2=v34.cross(v32);

	l1=std::sqrt(n1.dot(n1));
	l2=std::sqrt(n2.dot(n2));

    // cos of the angle between the normals to the 2 planes
    ct=n1.dot(n2)/(l1*l2);

    if (ct >  1.0)
        ct = 1.0;
    if (ct < (-1.0))
        ct = -1.0;

	angle = std::acos(ct);

    sign = v32[0] * (n1[2] * n2[1] - n1[1] * n2[2])
         + v32[1] * (n1[0] * n2[2] - n1[2] * n2[0])
         + v32[2] * (n1[1] * n2[0] - n1[0] * n2[1]);

    if (sign < 0.0)
        angle = -angle;

    angle = (angle>0.0) ? 3.14159265359-angle : -(3.14159265359+angle);
	// angle is degree [0..360]
	angle = angle*180/3.14159265359;
    if (angle < 0) angle +=360;


    return angle;

} // end of torsion


bool svd3d(double a[3][3], double w[3], double v[3][3])
{
  int flag, i, its, j, jj, k, l, nm=0;
  double anorm, c, f, g, h, s, scale,  x, y, z;
  double rv1[3];

  g = scale = anorm = 0.0;

  // Householder reduction to bidiagonal form

  for (i = 0; i < 3; ++i) {
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;
    if(i < 3) {
      for(k=i; k<3; ++k) scale += fabs(a[k][i]);
      if(scale) {
	for(k=i; k<3; ++k) {
	  a[k][i] /= scale;
	  s += a[k][i] * a[k][i];
	}
	f = a[i][i];
	g = -sign(sqrt(s),f);
	h = f * g - s;
	a[i][i] = f-g;
	for ( j= l; j < 3; ++j) {
	  for(s= 0.0, k= i; k<3; ++k) s+= a[k][i] * a[k][j];
	  f = s/h;
	  for(k=i; k<3;++k) a[k][j] += f * a[k][i];
	}
	for(k=i;k<3;++k) a[k][i] *= scale;
      }
    }
    w[i] = scale * g;
    g = s = scale = 0.0;
    if(i<2) {
      for (k = l; k < 3; ++k) scale += fabs(a[i][k]);
      if(scale) {
	for(k= l; k<3; ++k) {
	  a[i][k] /= scale;
	  s += a[i][k] * a[i][k];
	}
	f = a[i][l];
	g = -sign(sqrt(s),f);
	h= f*g-s;
	a[i][l] = f-g;
	for(k = l; k < 3; ++k) rv1[k] = a[i][k]/h;
	for ( j = l; j < 3; ++j) {
	  for (s = 0.0, k= l; k < 3; ++k) s += a[j][k] * a[i][k];
	  for (k = l; k < 3; ++k) a[j][k] += s * rv1[k];
	}
	for(k = l; k < 3; ++k) a[i][k] *= scale;
      }
    }
    anorm = max2(anorm, (fabs(w[i]) + fabs(rv1[i])));
  }

  // Accumulation of right-hand transformations

  for (i = 2; i >= 0; --i) {
    if (i < 2) {
      if(g) {
	for (j = l; j < 3; ++j)
	  v[j][i] = (a[i][j]/a[i][l])/g;
	for(j = l; j < 3; ++j) {
	  for (s = 0.0, k = l; k < 3; ++k) s += a[i][k] * v[k][j];
	  for(k = l; k < 3; ++k)  v[k][j] += s * v[k][i];
	}
      }
      for(j = l; j < 3; ++j) v[i][j] = v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g= rv1[i];
    l = i;
  }

  // Accumulation of left-hand transformations

  for(i = 2; i >= 0; --i) {
    l = i + 1;
    g = w[i];
    for(j = l; j < 3; ++j) a[i][j] = 0.0;
    if(g) {
      g = 1.0/g;
      for ( j = l; j < 3; ++j) {
	for ( s= 0.0, k = l; k < 3; ++k) s+= a[k][i] * a[k][j];
	f = (s/a[i][i]) * g;
	for( k = i; k < 3; ++k) a[k][j] += f * a[k][i];
      }
      for(j = i; j < 3; ++j) a[j][i] *= g;
    }
    else for (j = i; j < 3; ++j) a[j][i] = 0.0;
    ++a[i][i];
  }

  // Diagonalization of bidiagonal form

  for(k = 2; k >= 0; --k)   {
    for(its= 1; its <= 30; ++its) {
      flag = 1;
      for(l=k; l >= 0; --l) {
	nm = l-1;
	if((double)(fabs(rv1[l]) + anorm) == (double)anorm) {
	  flag = 0;
	  break;
	}
	if((double)(fabs(w[nm]) + anorm) == (double)anorm) break;
      }
      if(flag) {
	c = 0.0;
	s = 1.0;
	for(i = l; i <= k; ++i) {
	  f = s * rv1[i];
	  rv1[i] = c * rv1[i];
	  if ((double)(fabs(f) + anorm) == (double)anorm) break;
	  g = w[i];
	  h = pythag(f,g);
	  w[i] = h;
	  h = 1.0/h;
	  c = g * h;
	  s = -f * h;
	  for( j = 0; j < 3; ++j) {
	    y = a[j][nm];
	    z = a[j][i];
	    a[j][nm] = y * c + z * s;
	    a[j][i] = z * c - y * s;
	  }
	}
      }
      z = w[k];
      if( l == k) {
	if ( z < 0.0) {
	  w[k] = -z;
	  for ( j = 0; j < 3; ++j) v[j][k] = (-v[j][k]);
	}
	break;
      }
      if(its == 30) return true;
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y-z) * (y + z) + (g - h) * (g + h))/(2.0*h*y);
      g = pythag(f,1.0);
      f= ((x-z) * (x+z) + h * ((y/(f + sign(g,f))) - h))/x;
      c=s = 1.0;
      for(j = l; j <= nm; ++j) {
	i = j + 1;
	g = rv1[i];
	y = w[i];
	h = s * g;
	g = c * g;
	z = pythag(f,h);
	rv1[j] = z;
	c = f/z;
	s = h/z;
	f = x * c + g * s;
	g = g * c - x * s;
	h = y * s;
	y *= c;
	for ( jj= 0; jj < 3; ++jj) {
	  x = v[jj][j];
	  z = v[jj][i];
	  v[jj][j] = x * c + z * s;
	  v[jj][i] = z * c - x * s;
	}
	z = pythag(f,h);
	w[j] = z;
	if(z) {
	  z = 1.0/z;
	  c = f * z;
	  s = h * z;
	}
	f = (c * g) + (s * y);
	x = (c * y) - (s * g);
	for(jj = 0; jj < 3; ++jj) {
	  y = a[jj][j];
	  z = a[jj][i];
	  a[jj][j] = y * c + z * s;
	  a[jj][i] = z * c - y * s;
	}
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }
  return false;
}


void build_correlation_matrix(std::vector<Vec3d>& points, double* mat) {
  double x_mean,y_mean,z_mean,var,var_x,var_y,var_z;
  var=x_mean=y_mean=z_mean=var_x=var_y=var_z=0.0;
  unsigned int i;
  for(i = 0 ; i < points.size() ; i++) {
    x_mean+=(points[i])[0];
    y_mean+=(points[i])[1];
    z_mean+=(points[i])[2];
  }

  x_mean/=((double)points.size());
  y_mean/=((double)points.size());
  z_mean/=((double)points.size());

  for(i = 0 ; i < points.size() ; i++) {
    var_x+=(x_mean-(points[i])[0])*(x_mean-(points[i])[0]);
    var_y+=(y_mean-(points[i])[1])*(y_mean-(points[i])[1]);
    var_z+=(z_mean-(points[i])[2])*(z_mean-(points[i])[2]);
  }

  var_x/=((double)points.size());
  var_y/=((double)points.size());
  var_z/=((double)points.size());

  mat[0]=var_x;mat[3]=var_y;mat[5]=var_z;

  var=0.0;
  for(i = 0 ; i < points.size() ; i++) {
      var+=((points[i])[0]-x_mean)*((points[i])[1]-y_mean);
  }

  var=var/((double)points.size());
  mat[1]=var;

  var=0.0;
  for(i = 0 ; i < points.size() ; i++) {
    var+=((points[i])[0]-x_mean)*((points[i])[2]-z_mean);
  }

  var=var/((double)points.size());
  mat[2]=var;

  var=0.0;
  for(i = 0 ; i < points.size() ; i++) {
    var+=((points[i])[1]-y_mean)*((points[i])[2]-z_mean);
  }

  var=var/points.size();
  mat[4]=var;
}


void compute_eigen_vector(double *a, double* eigenvector) {
  double fullMatrix[3][3];
  double eigenValuesMatrix[3];
  double eigenVectorMatrix[3][3];

  fullMatrix[0][0] = (double) a[0];
  fullMatrix[1][0] = (double) a[1];
  fullMatrix[2][0] = (double) a[2];

  fullMatrix[0][1] = (double) a[1];
  fullMatrix[1][1] = (double) a[3];
  fullMatrix[2][1] = (double) a[4];

  fullMatrix[0][2] = (double) a[2];
  fullMatrix[1][2] = (double) a[4];
  fullMatrix[2][2] = (double) a[5];

  svd3d(fullMatrix, eigenValuesMatrix, eigenVectorMatrix);

  int maxEigenValueIndex = 0;
  for (int i = 1 ; i < 3 ; i++) {
    if (eigenValuesMatrix[i] > eigenValuesMatrix[maxEigenValueIndex]) {
      maxEigenValueIndex = i;
    }
  }

  eigenvector[0] = eigenVectorMatrix[0][maxEigenValueIndex];
  eigenvector[1] = eigenVectorMatrix[1][maxEigenValueIndex];
  eigenvector[2] = eigenVectorMatrix[2][maxEigenValueIndex];
}


double Partial_Data::compute_match(pose_ns::Pose const& p)
{
  //std::cout<<"computing distance with respect to partial data"<<std::endl;
	double measure = 0, tmp_measure = 0;
	unsigned int i,j;
    double MAX_HH = 100,//_owner->get_src_pose().get_0D_score(pose_ns::HB_SRBB)*1000,
    MAX_SHEET = 10;

    if (false){//_helix_flag){
	  pose_ns::Pose trg;
      pose_ns::Pose src;
      src = _owner->get_src_pose();
      trg = _owner->get_trg_partial_data_pose();

      FArray1D_bool allow_repack1(src.total_residue(), true);
	  FArray1D_bool allow_repack2(trg.total_residue(), true);
	  trg.repack( allow_repack2, true/*include_current*/ );
      src.repack( allow_repack1, true/*include_current*/ );

      /*if (_owner->is_full_atom()){

        trg.set_fullatom_flag( false );

	  }
      if (_owner->is_full_atom()){
      src.set_fullatom_flag( false );
      }*/

      std::cout<<"==================================================\n";
	  std::cout<<"SHEET compute target score = ";
      //trg.score( score4); src.score( score4);
	  //std::cout<<trg.get_0D_score(pose_ns::SHEET)<<std::endl;
      //std::cout<<"SHEET compute source score = ";
      //std::cout<<src.get_0D_score(pose_ns::SHEET)<<std::endl;

      trg.score( score12); src.score( score12);
      trg.show_scores(std::cout); std::cout<<std::endl;std::cout<<"SOURCS = "<<std::endl;
      src.show_scores(std::cout); std::cout<<std::endl;

	}

	if (_angle_flag){ // degrees
		assert( _LMS_indices_A.size() == _LMS_indices_B.size() );
		for(unsigned int i = 0; i < _LMS_indices_A.size(); i++){
			//assert( _LMS_indices_A[i].size() == _LMS_indices_B[i].size() );
			std::vector<Vec3d> pa1 = extract_points_from_pose_to_Vec3d(p,_LMS_indices_A[i]);
			std::vector<Vec3d> pa2 = extract_points_from_pose_to_Vec3d(p,_LMS_indices_B[i]);
			Line l1 = Line().line_fitting(pa1);
			Line l2 = Line().line_fitting(pa2);
			Vec3d p1,p2,p3,p4;
			p1 = l1.project_point(pa1[pa1.size()-1]);
			l1.compute_segment_dist(l2,p2,p3);
			p4 = l2.project_point(pa2[pa2.size()-1]);
			tmp_measure = torsion4(p1,p2,p3,p4); // degrees (0..360)
			double line_measure = std::abs(tmp_measure - _line_angles[i]);
			measure += line_measure * weight_angle;
			std::cout<<"Line Angle : " << line_measure << " weight " << weight_angle << std::endl;
		}
	}
	if (_centroid_flag){
		assert(_centroid_indices_A.size() == _centroid_indices_B.size() );
		for(unsigned int i = 0; i < _centroid_indices_A.size(); i++){
			//assert( _centroid_indices_A[i].size() == _centroid_indices_B[i].size() );
			std::vector<Vec3d> pa1 = extract_points_from_pose_to_Vec3d(p,_centroid_indices_A[i]);
			std::vector<Vec3d> pa2 = extract_points_from_pose_to_Vec3d(p,_centroid_indices_B[i]);
			Vec3d cen1 = compute_center_mass(pa1);
			Vec3d cen2 = compute_center_mass(pa2);
			Vec3d cen = cen1 - cen2;
			tmp_measure = cen.norm();
			double centroid_measure = std::abs(tmp_measure -_centroid_distances[i]);
			measure += centroid_measure * weight_centroid;
			std::cout << "Centroid: " << centroid_measure << " weight " << weight_centroid << std::endl;
		}
	}

	if (_dist_line_flag){ // ok
		assert(_dist_line_indices_A.size() == _dist_line_indices_B.size() );
		for(unsigned int i = 0; i < _dist_line_indices_A.size(); i++){
			//assert( _dist_line_indices_A[i].size() == _dist_line_indices_B[i].size() );
			std::vector<Vec3d> pa1 = extract_points_from_pose_to_Vec3d(p,_dist_line_indices_A[i]);
			std::vector<Vec3d> pa2 = extract_points_from_pose_to_Vec3d(p,_dist_line_indices_B[i]);
			Line l1 = Line().line_fitting(pa1);
			Line l2 = Line().line_fitting(pa2);
			tmp_measure = l1.distance_line(l2);
			double dist_line_measure = std::abs(tmp_measure -_line_distances[i]);
			measure += dist_line_measure * weight_dist_line;
			std::cout << "Line Dist: " << dist_line_measure << " weight " << weight_dist_line << std::endl;
		}
	}
	if (_match){ // unique case where we compare 2 diff poses
		assert( _source_indices.size() == _partial_data_indices.size() );
		for(unsigned int i = 0; i < _source_indices.size(); i++){
		  //std::cout<<"match "<<i<<std::endl;
			assert( _source_indices[i].size() == _partial_data_indices[i].size() );
		    FArray2D_double source = extract_points_from_pose(p,_source_indices[i]);
		    FArray2D_double target = extract_points_from_pose(_owner->get_trg_partial_data_pose(),_partial_data_indices[i]);
		    //measure += svd_CA_match(source,_target[i],_source_indices[i].size())*weight_match;
		    double match_measure = svd_CA_match(source,target,_source_indices[i].size());
		    measure += match_measure * weight_match;
    		std::cout << "Match: " << match_measure << " weight " << weight_match << std::endl;
                   //std::vector<Vec3d> pa1 = extract_points_from_pose_to_Vec3d(p,_source_indices[i]);
                   //std::vector<Vec3d> pa2 =
			extract_points_from_pose_to_Vec3d(_owner->get_trg_partial_data_pose(),_partial_data_indices[i]);
                   //measure = std::abs((pa1[0] - pa1[1]).norm() -  (pa2[0] - pa2[1]).norm());

        }

	}

    if (_match_rmsd){
		assert( _source_rmsd_indices.size() == _partial_data_rmsd_indices.size() );
		double match_rmsd_measure = rmsd_CA_match(p);
		measure += match_rmsd_measure * weight_match;
		std::cout << "Match RMSD: " << match_rmsd_measure << " weight " << weight_match << std::endl;
	}

	if (_beta_flag){

	  //std::cout<<p.get_0D_score(pose_ns::SHEET)<<std::endl;
          //pose_ns::Pose trg;
          //trg = _owner->get_trg_partial_data_pose();
          //trg.score( score4);
          //std::cout<<trg.get_0D_score(pose_ns::SHEET)<<std::endl;
          //exit(0);
        if (_owner->is_full_atom()){
            pose_ns::Pose centroid_pose;
	        centroid_pose = p;
	        centroid_pose.set_fullatom_flag( false );
            centroid_pose.score( score4);
            //measure += (MAX_SHEET - centroid_pose.get_0D_score(pose_ns::SHEET));
            measure += centroid_pose.get_0D_score(pose_ns::SSPAIR);//SHEET);
	  }

      else{
            //p.score(score4);
            //measure += (MAX_SHEET - p.get_0D_score(pose_ns::SHEET));
            measure += p.get_0D_score(pose_ns::SSPAIR);//SHEET);
	  }
	  //std::cout<<"beta sheet score = "<<p.get_0D_score(pose_ns::SHEET)<<std::endl;


	  /*for(unsigned int i = 0; i < _beta_indices.size(); i++){
	FArray2D_double beta = extract_points_from_pose(p,_beta_indices[i]);
	// tmp_measure = check score for beta???
	measure += tmp_measure*weight_beta;
	}*/
	}
	if (_helix_flag){

	  //double new_score = p.get_0D_score(pose_ns::HB_SRBB)+ p.get_0D_score(pose_ns::HB_LRBB); //HB_LRBB ,  HB_SC
      double alpha_score = p.get_0D_score(pose_ns::HB_SC)+ p.get_0D_score(pose_ns::HB_LRBB); //HB_LRBB ,  HB_SC
      weight_alpha = 10;
      measure += alpha_score*weight_alpha;
      std::cout << "Match Helix: " << alpha_score << " weight " << weight_alpha << std::endl;


	/*
		for(unsigned int i = 0; i < _alpha_indices.size(); i++){
           double tmp_measure_psi =0; // to ask Barak if in rosetta (0..360) -40 i.e. 320
		   double tmp_measure_phi =0; // -60 i.e. 360
           for(unsigned int j = 0; j < _alpha_indices[i].size(); j++){
			   tmp_measure_psi += std::pow(p.psi(_alpha_indices[i][j]) - CANONICAL_HELIX_PSI,2);
			   tmp_measure_phi += std::pow(p.phi(_alpha_indices[i][j]) - CANONICAL_HELIX_PHI,2);
		   }
		   measure += (std::sqrt(tmp_measure_psi+tmp_measure_phi))*weight_alpha;
		}
	*/
	}
	return measure;
}


Vec3d Partial_Data::compute_center_mass(std::vector<Vec3d> const & points){
	Vec3d res(0,0,0);
	for(unsigned int i = 0; i < points.size(); i++){
	  res = res + points[i];
	}
    res = res / ((double)points.size());
	return res;
}

  std::vector<Vec3d> Partial_Data::extract_points_from_pose_to_Vec3d
	  (pose_ns::Pose const & p,std::vector<int>& index_list)
  {

	std::vector<Vec3d> pa;
	bool local_debug = false;
    FArray3D_float const & Epos ( p.Eposition() );
    //TODO: do we want full?
    //FArray3D_float const & full_coord ( p.full_coord() );
   if(local_debug){
    std::cout<<"Extract:\n";
   }
	for (unsigned int i = 0; i < index_list.size(); i++ )
	{   assert( index_list[i] < p.total_residue() );
	    static int const atom_index_protein ( 2 ); // CA
		   pa.push_back(Vec3d(Epos(1,atom_index_protein,index_list[i]),
			              Epos(2,atom_index_protein,index_list[i]),
				      Epos(3,atom_index_protein,index_list[i])));
                   if(local_debug){
		      std::cout<<"Vec "<<index_list[i]<<" = "<<pa[pa.size()-1]<<std::endl;
		   }
	}
	return pa;

  }


  FArray2D_double Partial_Data::extract_points_from_pose(pose_ns::Pose const & p,std::vector<int>& index_list)
  {
	  FArray2D_double pa(3,index_list.size()); // TODO: is this the correct way to initialize? check compiler warning
	  //using namespace param_aa;
	  FArray3D_float const & Epos ( p.Eposition() );

	  for (unsigned int i = 0; i < index_list.size(); i++ )
	  {
		  assert( index_list[i] < p.total_residue() );
	      for (int k = 1; k <= 3; k++) {
		     static int const atom_index_protein ( 2 ); // CA
		     pa(k,i+1) = Epos( k, atom_index_protein, index_list[i]);
	      }
	  }
	  return pa;
  }


  double  Partial_Data::svd_CA_match( FArray2D_double& pa1, FArray2D_double& pa2,int natoms)
  {
    //using namespace param_aa;
    if (natoms == 0) return 0.0;
    // calc rms
    FArray1D_double const ww( natoms, 1.0 );
    FArray2D_double uu( 3, 3 );
    double ctx;
    findUU( pa1, pa2, ww, natoms, uu, ctx );
    float fast_rms;
    calc_rms_fast( fast_rms, pa1, pa2, ww, natoms, ctx );
    return fast_rms;
  }


  double  Partial_Data::rmsd_CA_match(const pose_ns::Pose& pose)
  {

      double sum_square_dev=0;
      int const CA = 2;//atom code for C-Alpha
      const FArray3D_float & Epos1( pose.Eposition() ); // set of backbone coordinates
      static const FArray3D_float & Epos2( _owner->get_trg_partial_data_pose().Eposition() );

      for(unsigned int i = 0; i < _partial_data_rmsd_indices.size(); i++){
        int res_id_s = _source_rmsd_indices[i];
        int res_id_t = _partial_data_rmsd_indices[i];
	for(int coord=1 ; coord <= 3 ; coord++) // X,Y,Z
	  sum_square_dev +=
	    pow(Epos1(coord,CA,res_id_s) - Epos2(coord,CA,res_id_t), 2);
      }
     return std::sqrt(sum_square_dev / _partial_data_rmsd_indices.size());

  }


  void Partial_Data::read_partial_data()
  {
    // set the flags to true - of the partial data that should be concidered
    // TODO: maybe we want to use some different types of partial data simultaniously
    // should be loaded from a file
    // case MATCH target:[A,15-20;B,14;16] source:[20-30;2;4]

    pose_ns::Pose const& src = _owner->get_src_pose(); // we are working on the source no matter what
    pdbres_to_poseres_pose(src,_map_pdbres_to_pose_source);

    // if read match {
    _match = true;
    //if (_owner->get_trg_partial_data_pose()){
    //pose_ns::Pose const& trg = _owner->get_trg_partial_data_pose();
    //pdbres_to_poseres_pose(trg,_map_pdbres_to_pose_target);
	  //}
    // read the residues + chains taht match in target and src:
    // currently - just testing
    test();
    //}

    // if read angle {
    //    _angle_flag = true;
    //    update _LMS_indices_A _LMS_indices_B _angles
    //}

    // if read centroid {
    //    _centroid_flag = true;
    //    update _centroid_indices_A _centroid_indices_B _distances
    //}

    // if read alpha {
    //    _helix_flag = true;
    //    update _alpha_indices
    //}

    // if read beta {
    //    _beta_flag = true;
    //    update _beta_indices
    //}

    // else{ std::cout<<"unknown partial_data_info"<<std::cout; exit(0);}


  }


  Partial_Data::Partial_Data(pathways::Pathways * owner) {
    _dist_line_flag = false; _angle_flag = false; _centroid_flag = false; _match = false; _beta_flag = false; _helix_flag = false; _match_rmsd = false;
    weight_match = 1; weight_angle = 0.1; weight_centroid = 1; weight_alpha = 1; weight_beta = 1,weight_dist_line = 1;
    _owner = owner;
	pose_ns::Pose const& src = _owner->get_src_pose(); // we are working on the source no matter what
    pdbres_to_poseres_pose(src,_map_pdbres_to_pose_source);

    //read_partial_data(); - TODO:remove
  }


  // initilize the mapping of the source template and target template (if the last is given)
  void Partial_Data::pdbres_to_poseres_pose(pose_ns::Pose const & _template_pose,t_map_pdbres_to_pose& map_pdbres_to_pose)
  {
    using namespace pdb;
    Pdb_info const & pdb_info = _template_pose.pdb_info();
    Pdbres_id pdbres;
    for(int i=1; i < _template_pose.total_residue(); i++)
      {
	pdbres.chain = pdb_info.res_chain(i);
	pdbres.res_id = pdb_info.pdb_res_num(i);
	pdbres.insert_letter = ' ';//pdb_info.pdb_insert_let(i);
	assert(map_pdbres_to_pose.find(pdbres) == map_pdbres_to_pose.end());
	map_pdbres_to_pose[pdbres] = i;
      }
  }
  // Looks up index in pose of a residue, according to PDB indexing
  int Partial_Data::pdbres_to_poseres(Pdbres_id pdbres_id,t_map_pdbres_to_pose& map_pdbres_to_pose)
  {
    t_map_pdbres_to_pose::const_iterator find_iter
      = map_pdbres_to_pose.find(pdbres_id);
    assert(find_iter != map_pdbres_to_pose.end());
    if(find_iter == map_pdbres_to_pose.end())
      return -1; // not found
    return find_iter->second;
  }


// --------------------------------------------------------------------------------

    //struct PD_line_angle{ char chain1;char chain2;int from1,from2,to1,to2;double angle;
  void Partial_Data::set_line_angle_recs(std::vector<PD_line_angle> line_angle_recs){

	  _angle_flag = true;

	  unsigned int i,j;
	  for(i =0; i < line_angle_recs.size(); i++){
		     std::vector<int> res_a;
			 std::vector<int> res_b;

	         char chain_1 = line_angle_recs[i].chain1;
             char chain_2 = line_angle_recs[i].chain2;
             unsigned int start_res_1 = line_angle_recs[i].from1;
             unsigned int end_res_1 = line_angle_recs[i].to1;
             unsigned int start_res_2 = line_angle_recs[i].from2;
			 unsigned int end_res_2 = line_angle_recs[i].to2;
			 Pdbres_id pdbres;
             for(j=start_res_1; j <= end_res_1;j++){
                      pdbres.chain = chain_1;
                      pdbres.res_id = j;
                      pdbres.insert_letter = ' ';
                      res_a.push_back(pdbres_to_poseres
						 (pdbres, _map_pdbres_to_pose_source));
			 }
			 for(j=start_res_2; j <= end_res_2;j++){
                      pdbres.chain = chain_2;
                      pdbres.res_id = j;
                      pdbres.insert_letter = ' ';
                      res_b.push_back(pdbres_to_poseres
						 (pdbres, _map_pdbres_to_pose_source));
			 }
             _LMS_indices_A.push_back(res_a);
			 _LMS_indices_B.push_back(res_b);
			 _line_angles.push_back(line_angle_recs[i].angle);
    }

  }

void Partial_Data::set_line_distance_recs(std::vector<PD_line_distance> line_distance_recs)
{
	_dist_line_flag = true;

 unsigned int i,j;
	  for(i =0; i < line_distance_recs.size(); i++){
		     std::vector<int> res_a;
			 std::vector<int> res_b;

	         char chain_1 = line_distance_recs[i].chain1;
             char chain_2 = line_distance_recs[i].chain2;
             unsigned int start_res_1 = line_distance_recs[i].from1;
             unsigned int end_res_1 = line_distance_recs[i].to1;
             unsigned int start_res_2 = line_distance_recs[i].from2;
			 unsigned int end_res_2 = line_distance_recs[i].to2;
			 Pdbres_id pdbres;
             for(j=start_res_1; j <= end_res_1;j++){
                      pdbres.chain = chain_1;
                      pdbres.res_id = j;
                      pdbres.insert_letter = ' ';
                      res_a.push_back(pdbres_to_poseres
						 (pdbres, _map_pdbres_to_pose_source));
			 }
			 for(j=start_res_2; j <= end_res_2;j++){
                      pdbres.chain = chain_2;
                      pdbres.res_id = j;
                      pdbres.insert_letter = ' ';
                      res_b.push_back(pdbres_to_poseres
						 (pdbres, _map_pdbres_to_pose_source));
			 }
             _dist_line_indices_A.push_back(res_a);
			 _dist_line_indices_B.push_back(res_b);
			 _line_distances.push_back(line_distance_recs[i].distance);
    }
}


void Partial_Data::set_centroid_distance_recs(std::vector<PD_centroid_distance> centroid_distance_recs)
{
	 _centroid_flag = true;
 unsigned int i,j;
	  for(i =0; i < centroid_distance_recs.size(); i++){
		     std::vector<int> res_a;
			 std::vector<int> res_b;

	         char chain_1 = centroid_distance_recs[i].chain1;
             char chain_2 = centroid_distance_recs[i].chain2;
             unsigned int start_res_1 = centroid_distance_recs[i].from1;
             unsigned int end_res_1 = centroid_distance_recs[i].to1;
             unsigned int start_res_2 = centroid_distance_recs[i].from2;
			 unsigned int end_res_2 = centroid_distance_recs[i].to2;
			 Pdbres_id pdbres;
             for(j=start_res_1; j <= end_res_1;j++){
                      pdbres.chain = chain_1;
                      pdbres.res_id = j;
                      pdbres.insert_letter = ' ';
                      res_a.push_back(pdbres_to_poseres
						 (pdbres, _map_pdbres_to_pose_source));
			 }
			 for(j=start_res_2; j <= end_res_2;j++){
                      pdbres.chain = chain_2;
                      pdbres.res_id = j;
                      pdbres.insert_letter = ' ';
                      res_b.push_back(pdbres_to_poseres
						 (pdbres, _map_pdbres_to_pose_source));
			 }
             _centroid_indices_A.push_back(res_a);
			 _centroid_indices_B.push_back(res_b);
			 _centroid_distances.push_back(centroid_distance_recs[i].distance);
    }
}

//char chain;	int from,to;
void Partial_Data::set_form_alpha_recs(std::vector<PD_form_alpha> form_alpha_recs){
	  _helix_flag = true;

unsigned int i,j;
	  for(i =0; i < form_alpha_recs.size(); i++){
		     std::vector<int> res;

	         char chain = form_alpha_recs[i].chain;
             unsigned int start_res = form_alpha_recs[i].from;
             unsigned int end_res = form_alpha_recs[i].to;
			 Pdbres_id pdbres;
             for(j=start_res; j <= end_res;j++){
                      pdbres.chain = chain;
                      pdbres.res_id = j;
                      pdbres.insert_letter = ' ';
                      res.push_back(pdbres_to_poseres
						 (pdbres, _map_pdbres_to_pose_source));
			 }
			 _alpha_indices.push_back(res);
	  }
}

//char chain_s;	char chain_t;	int from_s,from_t,to_s,to_t;
void Partial_Data::set_match_recs(std::vector<PD_match> match_recs){
  unsigned int i,j;

  _match = true;

  //if (_owner->get_trg_partial_data_pose()){
  //pose_ns::Pose const& trg = _owner->get_trg_partial_data_pose();
  //if (!_match_rmsd){
  //  pdbres_to_poseres_pose(trg,_map_pdbres_to_pose_target);
  //}

  std::sort(match_recs.begin(), match_recs.end(), by_id);

  int ref_id =match_recs[0].id; int curr_id = ref_id;
  std::vector<int> source_index_entity,partial_data_index_entity;

  //std::cout<<"Ids = " <<match_recs[0].id<<" "<<" Ids ="<<match_recs[1].id<<" "<<std::endl;
  //std::cout<<"done"<<std::endl;

	// push all residues matches to _source_indices and _partial_data_indices
	// matches are grouped by match-id (even if with more than one line)
	for(i =0; i < match_recs.size(); i++){
		curr_id =  match_recs[i].id;
		if(curr_id != ref_id){ //starting new match
			_source_indices.push_back(source_index_entity);
		    _partial_data_indices.push_back(partial_data_index_entity);
			source_index_entity.clear();
			partial_data_index_entity.clear();
			ref_id = curr_id;
			std::cout<<"curr id = "<< curr_id<<std::endl;
		}
		char chain_s = match_recs[i].chain_s;
        char chain_t = match_recs[i].chain_t;
        unsigned int start_res_s = match_recs[i].from_s;
        unsigned int end_res_s = match_recs[i].to_s;
        unsigned int start_res_t = match_recs[i].from_t;
		unsigned int end_res_t = match_recs[i].to_t;
		Pdbres_id pdbres;
        for(j=start_res_s; j <= end_res_s;j++){
        	pdbres.chain = chain_s;
            pdbres.res_id = j;
            pdbres.insert_letter = ' ';
            int pose_res = pdbres_to_poseres
				(pdbres, _map_pdbres_to_pose_source);
            source_index_entity.push_back(pose_res);
            std::cout << "SOURCE: " << chain_s << j << " => " << pose_res << std::endl;
		}
		for(j=start_res_t; j <= end_res_t;j++){
        	pdbres.chain = chain_t;
            pdbres.res_id = j;
            pdbres.insert_letter = ' ';
            int pose_res = pdbres_to_poseres
				(pdbres, _map_pdbres_to_pose_target);
		    std::cout << "PARTIAL: " << chain_t << j << " => " << pose_res << std::endl;
            partial_data_index_entity.push_back(pose_res);
		}
    }
    // flush last match
    if(source_index_entity.size()>0){
		_source_indices.push_back(source_index_entity);
		_partial_data_indices.push_back(partial_data_index_entity);
		source_index_entity.clear();
		partial_data_index_entity.clear();
		std::cout<<"curr id = "<< curr_id<<std::endl;
    }

//  //------------------- tmp -------------------
//          pose_ns::Pose const& trg = _owner->get_trg_partial_data_pose();
//          double tmp_measure = -1;
//         	       trg.dump_pdb("./output/check_TRG_partial.pdb"); // DEBUG
//  std::cout<<"#matching groups = "<<_partial_data_indices[0].size()<<" "<<_partial_data_indices[1].size()<<std::endl;
//  std::cout << "got here 0" << std::endl;
//  std::vector<Vec3d> pa1 = extract_points_from_pose_to_Vec3d(trg,_partial_data_indices[0]);
//  std::vector<Vec3d> pa2 = extract_points_from_pose_to_Vec3d(trg,_partial_data_indices[1]);
//  std::cout << "got here 1" << std::endl;
//  Line l1 = Line().line_fitting(pa1);
//  Line l2 = Line().line_fitting(pa2);
//  Vec3d p1,p2,p3,p4;
//  std::cout << "got here 2" << std::endl;
//  p1 = l1.project_point(pa1[pa1.size()-1]);
//  l1.compute_segment_dist(l2,p2,p3);
//  std::cout << "got here 3" << std::endl;
//  p4 = l2.project_point(pa2[pa2.size()-1]);
//  tmp_measure = torsion4(p1,p2,p3,p4); // degrees (0..360)
//  std::cout<<"Angle between helides = "<<tmp_measure<<std::endl;
//
//  Vec3d cen1 = compute_center_mass(pa1);
//  Vec3d cen2 = compute_center_mass(pa2);
//  Vec3d cen = cen1 - cen2;
//  std::cout<<"Centers = "<<cen1<<" "<<cen2<<std::endl;
//  tmp_measure = cen.norm();
//  std::cout<<"Distance betweeb helix centroids = "<<tmp_measure<<std::endl;
//
//  tmp_measure = l1.distance_line(l2);
//  std::cout<<"Distance between helices = "<<tmp_measure<<std::endl;
//
//  exit(0);

//Angle between helides = 295.077
//Distance betweeb helix centroids = 17.8758
//Distance between helices = 10.081


  //RESULTS:
  //Angle between helides = 289.965
  //Distance betweeb helix centroids = 18.2784
  //Distance between helices = 8.8846

  //-------------------------------------------


  //trg.dump_pdb("./output/PARTIAL_TARGET.pdb");

	// TODO: what are these lines good for? (Barak)
  _target.resize(_partial_data_indices.size());
  for(i =0; i < _partial_data_indices.size(); i++){
    _target[i] = extract_points_from_pose(_owner->get_trg_partial_data_pose(),_partial_data_indices[i]);
  }

}


void Partial_Data::set_partial_data(){
  pose_ns::Pose const& trg = _owner->get_trg_partial_data_pose();
  pdbres_to_poseres_pose(trg,_map_pdbres_to_pose_target);
}


void Partial_Data::set_match_rmsd_recs(std::vector<PD_match_rmsd> match_rmsd_recs){
  unsigned int i,j;

  _match_rmsd = true;

  //pose_ns::Pose const& trg = _owner->get_trg_partial_data_pose();
  //if (!_match){
  //  std::cout<<"MATCH = false %%%%%%%%%%%%%%%%%%%%%%%"<<std::endl;
  //  pdbres_to_poseres_pose(trg,_map_pdbres_to_pose_target);
  //}

  for(i =0; i < match_rmsd_recs.size(); i++){

	     char chain_s = match_rmsd_recs[i].chain_s;
             char chain_t = match_rmsd_recs[i].chain_t;
             unsigned int start_res_s = match_rmsd_recs[i].from_s;
             unsigned int end_res_s = match_rmsd_recs[i].to_s;
             unsigned int start_res_t = match_rmsd_recs[i].from_t;
	     unsigned int end_res_t = match_rmsd_recs[i].to_t;
	     Pdbres_id pdbres;
             assert(std::abs(start_res_s - end_res_s) == std::abs(start_res_t - end_res_t));

	     for(j=start_res_s; j <= end_res_s;j++){
	       pdbres.chain = chain_s;
	       pdbres.res_id = j;
	       pdbres.insert_letter = ' ';
	       _source_rmsd_indices.push_back(pdbres_to_poseres
		 (pdbres, _map_pdbres_to_pose_source));

	     }
	     for(j=start_res_t; j <= end_res_t;j++){
	       pdbres.chain = chain_t;
	       pdbres.res_id = j;
	       pdbres.insert_letter = ' ';
	       _partial_data_rmsd_indices.push_back(pdbres_to_poseres
		 (pdbres, _map_pdbres_to_pose_target));
	     }
  }
   std::cout<<"End set_match_rmsd_recs"<<std::endl;
   //trg.dump_pdb("./output/PARTIAL_TARGET.pdb");

}


// read teh residues from the file and add the corresponding indices to
  void Partial_Data::test(){
    // suppose we have the following line (cesT)
    //MATCH target:[A,3-12;A,95-111] source:[A,5-14;A,110-126]

    char chain_s = 'A';
    char chain_t = 'A';
    int start_res_source1 = 5;
    int end_res_source1 = 14;
    int start_res_target1 = 3;
    int end_res_target1 = 12;

    int start_res_source2 = 110;
    int end_res_source2 = 126;
    int start_res_target2 = 95;
    int end_res_target2 = 111;
    /*
    int i;
    Pdbres_id pdbres;
    for(i=start_res_source1; i <= end_res_source1;i++){
      pdbres.chain = chain_s;
      pdbres.res_id = start_res_source1;
      pdbres.insert_letter = ' ';
      _source_indices.push_back(pdbres_to_poseres(pdbres, _map_pdbres_to_pose_source));
    }

    for(i=start_res_source2; i <= end_res_source2;i++){
      pdbres.chain = chain_s;
      pdbres.res_id = start_res_source2;
      pdbres.insert_letter = ' ';
      _source_indices.push_back(pdbres_to_poseres(pdbres, _map_pdbres_to_pose_source));
    }

    for(i=start_res_target1; i <= end_res_target1;i++){
      pdbres.chain = chain_t;
      pdbres.res_id = start_res_target1;
      pdbres.insert_letter = ' ';
      _partial_data_indices.push_back(pdbres_to_poseres(pdbres, _map_pdbres_to_pose_target));
    }

    for(i=start_res_target2; i <= end_res_target2;i++){
      pdbres.chain = chain_t;
      pdbres.res_id = start_res_target2;
      pdbres.insert_letter = ' ';
      _partial_data_indices.push_back(pdbres_to_poseres(pdbres, _map_pdbres_to_pose_target));
    }
    */


  }


// class Line
// ----------------------------------------------------------------

  	  Line::Line(const Vec3d& point, const Vec3d& direction){
	    _point = point;
        _direction = direction;
	  }
	  Line::Line(double pointCoordinates[3], double directionCoordinates[3]){
	     _point = Vec3d(pointCoordinates);
         _direction = Vec3d(directionCoordinates);
	  }

	  Vec3d Line::project_point(Vec3d& p){
		  double t;

		  double directionNorm2 = _direction.norm2();

		  if(directionNorm2 != 0.0 ) {
			  t = ( _direction[0] * (p[0] - _point[0]) +
				    _direction[1] * (p[1] - _point[1]) +
				    _direction[2] * (p[2] - _point[2]) ) / directionNorm2;

			  return Vec3d(t * _direction[0] + _point[0],
				  t * _direction[1] + _point[1],
				  t * _direction[2] + _point[2]);
		  }
		  else {
		    std::cout << "ERROR: line direction = zero" << std::endl;
			  exit(1);
		  }

	  }


  void Line::compute_segment_dist(const Line& line,Vec3d& p1,Vec3d& p2)
  {
	  Vec3d u = _direction;
      Vec3d v = line.get_direction();
      Vec3d w0 = _point - line.get_point();

  double a = u * u;        // always >= 0
  double b = u * v;
  double c = v * v;        // always >= 0
  double d = u * w0;
  double e = v * w0;

  double D = a*c - b*b;    // always >= 0

  // Compute the line parameters of the two closest points
  double sc, tc;

  if (D < 0.001) {
    // the lines are almost parallel
    sc = 0.0;
    tc = (b>c ? d/b : e/c);   // use the largest denominator
  } else {
    sc = (b*e - c*d) / D;
    tc = (a*e - b*d) / D;
  }

  // compute the two closest points
  p1 =  _point + (sc * u);
  p2 =  line.get_point() + (tc * v);


  }


   // distance between two lines
  // |[(p0-q0),u,v]| / |uxv|
  // |[(p0-q0),u,v]| - is the determinat of a 3x3 matrix whose columns
  // are the vectors: (p0-q0), u and v.
  //  x - is the cross product
  double Line::distance_line(const Line& line)
  {
	  Vec3d subPoint = _point-line.get_point();

	  Vec3d vrow[3];
      vrow[0] = subPoint;
      vrow[1] = line.get_direction();
      vrow[2] = _direction;

	  double determinant =
        vrow[1-1][1-1] * (vrow[2-1][2-1]*vrow[3-1][3-1] - vrow[2-1][3-1]*vrow[3-1][2-1])
      - vrow[1-1][2-1] * (vrow[2-1][1-1]*vrow[3-1][3-1] - vrow[2-1][3-1]*vrow[3-1][1-1])
      + vrow[1-1][3-1] * (vrow[2-1][1-1]*vrow[3-1][2-1] - vrow[2-1][2-1]*vrow[3-1][1-1]);


	  Vec3d verticalVector = _direction & line.get_direction();

	  double verticalVectorNorm = verticalVector.norm();
	  if (verticalVectorNorm != 0) {
		  return fabs(determinant) / verticalVectorNorm;
	  } else {
		  // The lines are parallel.
		  // So, the distance between the two lines is equal
		  // to the distance between any point on one of the
		  // lines to the other.
		  return distance_point(line.get_point());
	  }
  }


  double Line::distance_point(const Vec3d& point) {
      Vec3d subPoint = point - _point;
      Vec3d verticalVector = subPoint & _direction;
	  return std::abs(verticalVector.norm() / _direction.norm());
  }

  // min-square-distance line
  // best fit to the given set of points
  Line Line::line_fitting(std::vector<Vec3d>& points){
  double eigenvector[3];
  double mat[6];

  build_correlation_matrix(points, mat);

  compute_eigen_vector(mat, eigenvector);

  double x_mean,y_mean,z_mean;
  x_mean = y_mean = z_mean = 0.0;

  for(unsigned int i = 0 ; i < points.size() ; i++) {
    x_mean+=(points[i])[0];
    y_mean+=(points[i])[1];
    z_mean+=(points[i])[2];
  }

  x_mean /= (double)points.size();
  y_mean /= (double)points.size();
  z_mean /= (double)points.size();

  double pointCoordinates[3];
  pointCoordinates[0] = x_mean - eigenvector[0];
  pointCoordinates[1] = y_mean - eigenvector[1];
  pointCoordinates[2] = z_mean - eigenvector[2];

  double norm=std::sqrt(eigenvector[0]*eigenvector[0]+
		  eigenvector[1]*eigenvector[1]+
		  eigenvector[2]*eigenvector[2]);

  eigenvector[0] /= norm;
  eigenvector[1] /= norm;
  eigenvector[2] /= norm;

  return Line(pointCoordinates, eigenvector);


  }


} // namespace pathways


