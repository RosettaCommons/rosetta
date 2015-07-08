////////////////////////////////////////////////////////////////////////////////
//     pathways_partial_data_manager.h:
//     This class is for managing partial data about target conformation
//
//     Author: Barak Raveh & Angela Enosh
//     Created: 27/12/2007
//
// methods list:
//     TODO: fill this in
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_devel_path_rover_pathrover_partial_data_manager_hh
#define INCLUDED_devel_path_rover_pathrover_partial_data_manager_hh

#include <ObjexxFCL/ObjexxFCL.hh>
#include "pathways.h"
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include "pdb.h"
#include "pose.h"
#include "random_numbers.h"
#include "score_data.h"
#include "symmetry_info.h"
#include "orient_rms.h"
#include "pathways.h"
#include "pathways_parameters.h"


//#define M_PI 3.14159265359
#define CANONICAL_HELIX_PSI 320
#define CANONICAL_HELIX_PHI 300

// TODO: divide this to inheritence hierarchy? where different partial data are differnet objects derived
//       from the same abstract class? BUT - perhaps we should allow compund classes like this as well?
//       i.e. how do we allow the programmers to add partial info classes of their own, and still maintain
//       a joint weighted score?
//       Perhaps just join into a "score" function, like / in Rosetta Constraints framework?
namespace pathways {

  class Pathways;
  // partial_data_file
  // MATCH target A 15-20 B 14 16 source A 20-30 2 4
  // meaning owner has a target and [A,15-20...] of the target match to [20-30;2;4] of teh source
  // LINES_ANGLE A 15 20 A 20 30 = 45 (degrees) usually will follow the below line
  // CENTROIDS_DISTANCE A 15 20  B 20 30 = 5 (angstrom)
  // LINES_DISTANCE A 15 20  B 20 30 = 5 (angstrom)
  // FORM_ALPHA B 15 20

  class Vec3d{
  public:
	  std::vector< double > _point_3d;
	  Vec3d(){
		  _point_3d.resize(3);
		  _point_3d[0]=0;
		  _point_3d[1]=0;
		  _point_3d[2]=0;

	  };
	  Vec3d(double point[3]){
		  _point_3d.resize(3);
		  _point_3d[0]=point[0];
		  _point_3d[1]=point[1];
		  _point_3d[2]=point[2];

	  };

	  Vec3d(double x,double y,double z){
		  _point_3d.resize(3);
		  _point_3d[0]=x;
		  _point_3d[1]=y;
		  _point_3d[2]=z;

	  };

	  double dot(const Vec3d &v)
	  {
		  return (_point_3d[0]*v._point_3d[0] +
			      _point_3d[1]*v._point_3d[1] +
				  _point_3d[2]*v._point_3d[2]);

	  }

	  Vec3d cross(const Vec3d &v)
	  {
		  Vec3d r(_point_3d[1]*v._point_3d[2] - _point_3d[2]*v._point_3d[1],
			      _point_3d[2]*v._point_3d[0] - _point_3d[0]*v._point_3d[2],
			      _point_3d[0]*v._point_3d[1] - _point_3d[1]*v._point_3d[0]);

		  return r;

	  }


	  double norm() { return std::sqrt(_point_3d[0]*_point_3d[0]+
		                          _point_3d[1]*_point_3d[1]+
								  _point_3d[2]*_point_3d[2]);
	  }

	  double norm2() { return (_point_3d[0]*_point_3d[0]+
		                          _point_3d[1]*_point_3d[1]+
				   _point_3d[2]*_point_3d[2]);}

	  friend Vec3d operator&(const Vec3d& p1, const Vec3d& p2) {
         return Vec3d(p1._point_3d[1]*p2._point_3d[2]-p1._point_3d[2]*p2._point_3d[1],
		   p1._point_3d[2]*p2._point_3d[0]-p1._point_3d[0]*p2._point_3d[2],
		   p1._point_3d[0]*p2._point_3d[1]-p1._point_3d[1]*p2._point_3d[0]);
      }

	  friend Vec3d operator-(const Vec3d& p1, const Vec3d& p2) {
         return Vec3d(p1._point_3d[0]-p2._point_3d[0],
			          p1._point_3d[1]-p2._point_3d[1],
					  p1._point_3d[2]-p2._point_3d[2]);
      }
	  friend Vec3d operator+(const Vec3d& p1, const Vec3d& p2) {
         return Vec3d(p1._point_3d[0]+p2._point_3d[0],
			          p1._point_3d[1]+p2._point_3d[1],
					  p1._point_3d[2]+p2._point_3d[2]);
      }

      friend double operator*(const Vec3d& p1, const Vec3d& p2) {
	return (p1._point_3d[0]*p2._point_3d[0] +
		p1._point_3d[1]*p2._point_3d[1] +
		p1._point_3d[2]*p2._point_3d[2]);
      }

     friend Vec3d operator*(const Vec3d& p1, double s) {
       return Vec3d(p1._point_3d[0]*s,
		    p1._point_3d[1]*s,
		    p1._point_3d[2]*s);
      }

    friend Vec3d operator*(double s, const Vec3d& p1) {
       return Vec3d(p1._point_3d[0]*s,
		    p1._point_3d[1]*s,
		    p1._point_3d[2]*s);
      }
    /*  friend Vec3d& operator+=(const Vec3d& p) {
          _point_3d[0]+=p._point_3d[0];
		  _point_3d[1]+=p._point_3d[1];
		  _point_3d[2]+=p._point_3d[2];
          return *this;
	  }*/

	  friend Vec3d operator/(const Vec3d& p,const double n) {
         return Vec3d(p._point_3d[0]/n,p._point_3d[1]/n,p._point_3d[2]/n);
	 }


          double& operator[](unsigned int coord) {return _point_3d[coord];}
	  Vec3d(const std::vector<double>& point_3d):
	  _point_3d( point_3d )
	  {}

    friend std::ostream& operator<<(std::ostream& s,Vec3d& v){
      return s << "(" << v._point_3d[0] << " , " << v._point_3d[1] << " , " <<v._point_3d[2]<<") ";
    }
  };


  class Line{

  Vec3d _point;
  Vec3d _direction;
  public:
	  Line(){}
	  Line(const Vec3d& point, const Vec3d& direction);
	  Line(double pointCoordinates[3], double directionCoordinates[3]);

  inline const Vec3d& get_direction() const {
    return _direction;
  }

  inline const Vec3d& get_point() const {
    return _point;
  }

  //compute the 2 closest points between the 2 lines
  void compute_segment_dist(const Line& line,Vec3d& p1,Vec3d& p2);

  Vec3d project_point(Vec3d& p);
  // distance between two lines
  // |[(p0-q0),u,v]| / |uxv|
  // |[(p0-q0),u,v]| - is the determinat of a 3x3 matrix whose columns
  // are the vectors: (p0-q0), u and v.
  //  x - is the cross product
  double distance_line(const Line& line);

    double distance_point(const Vec3d& point);

  // min-square-distance line
  // best fit to the given set of points
  static Line line_fitting(std::vector<Vec3d>& points);
    //void build_correlation_matrix( std::vector<Vec3d>& points, double mat[6]);

    //void compute_eigen_vector(double a[6], double eigenvector[3]);
};


  class Partial_Data{

    //MATCH
    std::vector<std::vector<int> > _source_indices;
    std::vector<std::vector<int> > _partial_data_indices;


     //MATCH RMSD
    std::vector<int>  _source_rmsd_indices;
    std::vector<int>  _partial_data_rmsd_indices;

    // all the below options applied on the source template - the target is ignored
    //ANGLE
    std::vector<std::vector<int> > _LMS_indices_A; // the indices in the source to define the least mean sqaure line
    std::vector<std::vector<int> > _LMS_indices_B; // the indices in the source to define the least mean sqaure line
    std::vector<double> _line_angles; // in degrees

    //CENTROID
    std::vector<std::vector<int> > _centroid_indices_A; // the indices in the source to define the centroid
    std::vector<std::vector<int> > _centroid_indices_B; // the indices in the source to define the centroid
    std::vector<double> _centroid_distances;

	// distance between LMS lines
	std::vector<std::vector<int> > _dist_line_indices_A;
	std::vector<std::vector<int> > _dist_line_indices_B;
	std::vector<double> _line_distances;


    //Alpha,Beta
    std::vector<std::vector<int> > _alpha_indices;
    std::vector<std::vector<int> > _beta_indices;

    // general purpose
    pathways::Pathways* _owner;
    bool _dist_line_flag,_angle_flag,_centroid_flag,_match,_match_rmsd,_beta_flag,_helix_flag;
    typedef std::map<Pdbres_id, int> t_map_pdbres_to_pose; // map from pdb indexing to pose indexing
    t_map_pdbres_to_pose _map_pdbres_to_pose_source; // map from pdb indexing to pose indexing
    t_map_pdbres_to_pose _map_pdbres_to_pose_target; // map from pdb indexing to pose indexing
    double weight_dist_line,weight_match,weight_angle,weight_centroid,weight_alpha,weight_beta;
    std::vector<FArray2D_double> _target;

  public:

    void set_alpha_flag(){_helix_flag = true; }
    void set_beta_flag(){_beta_flag = true; }
    double compute_match(pose_ns::Pose const& p);
    FArray2D_double extract_points_from_pose(pose_ns::Pose const & p,std::vector<int>& index_list);
	std::vector<Vec3d> extract_points_from_pose_to_Vec3d
	  (pose_ns::Pose const & p,std::vector<int>& index_list);
	Vec3d compute_center_mass(std::vector<Vec3d> const & points);

    double  svd_CA_match( FArray2D_double& pa1, FArray2D_double& pa2,int natoms);
    double  rmsd_CA_match(const pose_ns::Pose& pose);
    void read_partial_data();
    // read the residues from the file and add the corresponding indices to
    void test();
    Partial_Data() {}

    Partial_Data(pathways::Pathways* owner);


    // initilize the mapping of the source template and target template (if the last is given)
    void pdbres_to_poseres_pose(pose_ns::Pose const & _template_pose,t_map_pdbres_to_pose& map_pdbres_to_pose);
    // Looks up index in pose of a residue, according to PDB indexing
    int pdbres_to_poseres(Pdbres_id pdbres_id,t_map_pdbres_to_pose& map_pdbres_to_pose);

	void set_line_angle_recs(std::vector<PD_line_angle> line_angle_recs);
	void set_line_distance_recs(std::vector<PD_line_distance> line_distance_recs);
	void set_centroid_distance_recs(std::vector<PD_centroid_distance> centroid_distance_recs);
    void set_form_alpha_recs(std::vector<PD_form_alpha> form_alpha_recs);
    void set_match_recs(std::vector<PD_match> match_recs);
    void set_match_rmsd_recs(std::vector<PD_match_rmsd> match_rmsd_recs);
    void set_partial_data();
  };


} // namespace

#endif
