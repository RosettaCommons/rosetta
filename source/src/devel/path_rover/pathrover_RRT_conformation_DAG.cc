////////////////////////////////////////////////////////////////////////////////
//     pathways_RRT_conformation_tree.h:
//     This class is for storing a tree of conformation pathways, for use by the
//     RRT (Rapidly Exploring Random Trees) Algorithm.
//     It contains data structures for a tree of conformations and its nodes
//
//     Originally Created in R++: 3/12/2007 (Barak Raveh & Angela Enosh)
//     Transformed to mini: 10/11/2009 (Barak Raveh)
//
// methods list:
//     TODO: fill this in
//
/////////////////////////////////////////////////////////////////////////////////

#include "pose.h"
#include "pathways_RRT_conformation_DAG.h"
#include "orient_rms.h"

namespace pathways
{

RRT_node::~RRT_node()
{
  if(_custom_info) delete _custom_info;
}


  void
  RRT_node::set_custom_info
  (Node_custom_info* ci)
  {
    if(_custom_info)
      delete _custom_info;
    _custom_info = ci;
  }


  double
  CA_RMSD_functor::operator ()(RRT_node const* node1, RRT_node const* node2)
  {
    // try to work with caching, if not - reproduces poses
    using namespace pose_ns;
    using namespace std;

    //Node_custom_info* p = NULL;
    //Node_rmsd_info* q = dynamic_cast< Node_rmsd_info* >(p);
    double sum_square_dev = 0.0;
    int n_res;
    Node_rmsd_info const* rmsd_info1 =
      dynamic_cast< Node_rmsd_info const* >(node1->get_custom_info());
    Node_rmsd_info const* rmsd_info2 =
      dynamic_cast< Node_rmsd_info const* >(node2->get_custom_info());
    if(rmsd_info1 && rmsd_info2)
      // caching is on
      {
	assert(rmsd_info1->_points.size() == rmsd_info2->_points.size());
	int n = rmsd_info1->_points.size();
	n_res = int(n / 3); // 3D
	for(int i=0 ; i < n; i++)
	  sum_square_dev +=
	    pow(rmsd_info1->_points[i] - rmsd_info2->_points[i], 2);
      }
    else
      // no caching - reproduce nodes poses
    {
      int const CA = 2;//atom code for C-Alpha // TODO: make nicer
      pose_ns::Pose* pose1 = node1->produce_pose(); // TODO: make sure this is the quickest
                                                    //(i.e. "=" is invoked only once, inliding internally
                                                    // by get_pose()
      pose_ns::Pose* pose2 = node2->produce_pose();
      n_res = pose1->size();
      assert(n_res == pose2->size());
      const FArray3D_float & Epos1( pose1->Eposition() ); // set of backbone coordinates
      const FArray3D_float & Epos2( pose2->Eposition() );
      // compute sum-square of all coordinates of all residues
      for(int res_id = 1 ; res_id <= n_res ; res_id++ ){
	for(int coord=1 ; coord <= 3 ; coord++) // X,Y,Z
	  sum_square_dev +=
	    pow(Epos1(coord,CA,res_id) - Epos2(coord,CA,res_id), 2);
      }
      delete pose1;
      delete pose2;
    }
    // return RMSD ( = per residue square deviation)
    return std::sqrt(sum_square_dev / n_res);
}


// L2-norm between dofs vectro
double Dofs_vector_L2_norm_functor::operator()(RRT_node const* node1, RRT_node const* node2)
{
  using namespace std;
  std::vector<double> const& v1 = node1->get_dofs_vector();
  std::vector<double> const& v2 = node2->get_dofs_vector();
  assert(v1.size() = v2.size());
  double sum_sqr_dev = 0.0;
  for(unsigned int i = 0; i < v1.size();i++)
    {
      // calculate "periodic" difference (minimal distance in 360 degrees range)
      double diff_i = v1[i] - v2[i];
      if (diff_i < - 180.0)
	diff_i += 2*180.0;
      else if (diff_i > 180.0)
	diff_i = diff_i - 2*180.0;
      sum_sqr_dev += pow(diff_i, 2);
    }
  return std::sqrt(sum_sqr_dev / v1.size()); // root of mean

}


// initialize object to contain all CA coordinates of input_pose
// into _points array, as caching for RMSD calculations
Node_rmsd_info::Node_rmsd_info(pose_ns::Pose const& input_pose)
{
  int const CA = 2;//atom code for C-Alpha // TODO: use libRosetta keys
  int const n_res = input_pose.size();
  FArray3D_float const& Epos( input_pose.Eposition() ); // set of backbone coordinates
  for(int res_id=1; res_id <= n_res; res_id++)
    {
      for(int coord = 1; coord <=3; coord++)
	_points.push_back( Epos(coord, CA, res_id) );
    }
}

} // namespace pathways
