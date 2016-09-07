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

#ifndef INCLUDED_devel_path_rover_pathrover_RRT_conformation_DAG_hh
#define INCLUDED_devel_path_rover_pathrover_RRT_conformation_DAG_hh


#include <core/pose/Pose.fwd.hh>
#include "pathways_DOFs_manager.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

namespace pathways {


// *************** Forward Declarations: ***************
class Distance_functor; // TODO: use functor ? template? see http://www.newty.de/fpt/functor.html
class CA_RMSD_functor;
class Node_custom_info;
class Node_rmsd_info;
class RRT_conformation_DAG;
class RRT_node;

// ********************** Utility Stuff: ***********************/

// custom info that can be attached to a tree node // TODO: should we template this?
struct Node_custom_info{
 virtual ~Node_custom_info() = default; // fake - just want to make it virtual // TODO is this necessary?
};

// info needed to quickly calculate RMSD between nodes
struct Node_rmsd_info :
  public Node_custom_info
{
  Node_rmsd_info(core::pose::Pose const& pose);

  // a vector of 3D coordinates to be compared for RMSD
  // (points are concatenated to 1D array of size 3 x n)
  // TODO: make it nicer...
  std::vector< float > _points;
};

// abstract functor class for computing distances between nodes
class Distance_functor
{
 public:
  virtual ~Distance_functor() = default;
  virtual double operator ()(RRT_node const* node1, RRT_node const* node2) = 0;
};

// functor for calculating RMSD between all CAs of poses of each node
class CA_RMSD_functor : public Distance_functor
{
 public:
  double operator ()(RRT_node const* node1, RRT_node const* node2) override;
};

// L-2 norm (euclidean) over DOFs vector, e.g. internal coordinatesy
class Dofs_vector_L2_norm_functor : public Distance_functor
{
 public:
  double operator ()(RRT_node const* node1, RRT_node const* node2) override;
};

// ******************** class RRT_node: ****************
class RRT_node{
  friend class RRT_conformation_DAG;
public:


 bool _is_added;
  // Constructor for an unconnected (yet) node
  // Node saves only DOFs vector for the pose using
  // p_dofs_manager->get_dofs_values_vector()
  //
  // params:
  //  &pose - pose of this node (used to extract DOFs etc.)
  //  &template_pose - ref to a template that is identical
  //                    to pose, except at DOFs
  //  &dofs_manager - ref to a DOFs manager that will be
  //                   used to extract / apply DOFs vector on "pose"
  //  is_score_valid - if true, assumes input pose has a pre-computed
  //                   & valid score (default FALSE)
  //  custom_info - any extra info attached to node
  //                (default NULL)
  //
  // NOTE: node is responsible for deleting custom_info
  // NOTE2: &template_pose and &dofs_manager are expected to exist
  //        throught the lifetime of the node
  RRT_node(
	   core::pose::Pose const& pose,
	   core::pose::Pose const& template_pose,
	   DOFs_manager& dofs_manager,
	   bool is_score_valid = false,
	   Node_custom_info* custom_info = nullptr)
    : _template_pose(template_pose),
      _dofs_manager(dofs_manager)
    {
      _dofs_vector = _dofs_manager.get_dofs_values_vector(pose);
      _level = -1; // replace -1 with UNDEFINED define
      //_selected_points = selected_points; // TODO: handle this!!!
      _orig_src_root_id =  -1; // TODO: replace -1 with UNDEFINED define
      _custom_info = custom_info;
      _owner_DAG = nullptr;
      _is_score_valid = is_score_valid;
      if(_is_score_valid)
	_score = pose.energies().total_energy();
      _is_added = false;
    }

  // Constructor for an unconnected (yet) node
  // Input is dofs_vactor + template
  // Real pose is the combination of the two
  //
  // params
  //  dofs_vector - the values of the free DOFs
  //                (see dofs_manager.get_dofs_values_vector())
  //  &template_pose - node pose is extracted by applying
  //                   dofs_vector to template
  //  &dofs_manager - ref to a DOFs manager that will be
  //                   used to extract / apply DOFs vector on "pose"
  //  custom info - any extra info attached to node
  //                (default NULL)
  //
  // NOTE: node is responsible for deleting custom_info
  // NOTE2: score is not initialized when using this cntr

  RRT_node(
	   std::vector<double> dofs_vector,
	   core::pose::Pose const& template_pose,
	   DOFs_manager& dofs_manager,
	   Node_custom_info* custom_info = nullptr)
    : _template_pose(template_pose),
    _dofs_manager(dofs_manager)
    {
      _dofs_vector = dofs_vector;
      _level = -1; // replace -1 with UNDEFINED define
      //_selected_points = selected_points; // TODO: handle this!!!
      _orig_src_root_id =  -1; // TODO: replace -1 with UNDEFINED define
      _custom_info = custom_info;
      _owner_DAG = nullptr;
      _is_score_valid = false;
      _is_added = false;
    }


  ~RRT_node();


  // returns a new pose that is modified at DOFs according to
  // template_pose using dofs_manager
  // (template pose is assumed identical except at the DOFs)
  //
  // NOTE: calling code is responsible for deleting the new pose
  core::pose::Pose* produce_pose() const
    {
      utility::vector1< conformation::ResidueOP > residues;
      for ( Size i=1; i<= pose.total_residue(); ++i ) {
	residues.push_back( pose.residue(i).clone() );
      }
      pose.
//       core::pose::Pose* p_return_pose = new core::pose::Pose();
//       *p_return_pose = _template_pose;
//       _dofs_manager.apply_dofs_values_vector(*p_return_pose, _dofs_vector);

//       return p_return_pose;
    }

  void apply_node_dofs_vector(core::pose::Pose& pose)
  {
	  _dofs_manager.apply_dofs_values_vector(pose, _dofs_vector);
  }

  int get_orig_src_root_id() const
    {
      return _orig_src_root_id;
    }

  RRT_conformation_DAG* get_owner_tree()
    {
      return  _owner_DAG;
    }

  // NOTE: user is no more responsible for deleting ci - ownership is by the node
  void set_custom_info(Node_custom_info* ci);

  Node_custom_info* get_custom_info()
    { return _custom_info; }

  Node_custom_info const* get_custom_info() const
    {	return _custom_info; }


  double compute_dist
    ( RRT_node * other_node,
      Distance_functor* dist_functor)
    {
      return (*dist_functor)(this, other_node); // distance between this node and other node
    }

  double compute_dist_CA_rmsd(RRT_node * other_node)
    {
      CA_RMSD_functor rmsd_functor;
      return compute_dist(other_node, &rmsd_functor);
    }

  double compute_dist_dofs_vector(RRT_node * other_node)
    {
      Dofs_vector_L2_norm_functor dofs_dist_functor;
      return compute_dist(other_node, &dofs_dist_functor);
    }

  int getLevel()
    { return _level; }

  void set_score(double score)
    {
      _score = score;
      _is_score_valid = true;
    }

  double calc_score_pose(core::pose::Score_weight_map& weight_map)
    {
      core::pose::Pose* p_scoring_pose = produce_pose();
      this->set_score( p_scoring_pose->score(weight_map) );
      delete p_scoring_pose;

      return _score;
    }

  // NOTE: assumes score was pre-calculated by pre-invoking
  // score(), otherwise might lead to unexpected results
  float get_score() const
    {
      assert(_is_score_valid);
      return _score; //_pose.get_0D_score(core::pose::SCORE);
      // TODO: add return value for invalid / non-updated score
    }

  std::vector<double> const& get_dofs_vector() const
    { return _dofs_vector; }

  std::vector< RRT_node* >& get_children()
    { return _children; }

  std::vector< RRT_node* >& get_parents()
    { return _parents; }


 private:
  	// add a child to the children list, and define this node as the parent of "child"
	// (NOTE: private = this can be used only internally or by friend classes, such as the DAG)
	void add_child(RRT_node* child){
		assert(child);
		child->_parents.push_back(this);
		// if this is the first parent of child, update _orig_src_root_id
		if(child->_parents.size() == 1){
			child->_orig_src_root_id = _orig_src_root_id;
			child-> _owner_DAG = _owner_DAG; // TODO: do we need this???
			child->_level = _level + 1;
			//std::cout << "add_child(): CHILD LEVEL = " << child->_level << std::endl;
		}
		_children.push_back(child);
	}
 private:
	core::pose::Pose const& _template_pose;
	DOFs_manager& _dofs_manager;
	std::vector<double> _dofs_vector; // assignments to free DOFs
					  // (assuming template for
					  // other DOFs)
	bool _is_score_valid;
	double _score;
	int _level; // distance from root (if more then one, refers to _orig_src_root_id)
	std::vector< RRT_node* > _parents;
	std::vector< RRT_node* > _children; // TODO: all neighbours?
	int _orig_src_root_id; // the root id of the first parent of this node (i.e. the original tree for which it belongs)
  //	double _score, _score_path; // saving energy calculations TODO: adaptations for rosetta? caching info? do we need this - it might be already inside Pose
  //std::vector<Vector3> _selected_points;  // Some selcted point to compute distances
											// between two nodes in the conformation graph // TODO: vector3...

	RRT_conformation_DAG* _owner_DAG; // TODO: do we need this???
	Node_custom_info* _custom_info; // custom node info like RMSD, etc.
}; // RRT Node
  /** *********************************************  **/


// class for trees of conformations (in RRT algorithm, etc.)
class RRT_conformation_DAG{
 private:
  typedef std::vector< RRT_node * > t_vec_nodes;
  typedef std::vector< std::vector< RRT_node * > > t_vec2D_nodes;
  typedef std::map<int, core::pose::Pose> t_map_templates;
 private:
  t_vec_nodes _roots;
  t_vec2D_nodes _nodes; // list of all nodes in the DAG, for each root:
                        // 1st dimension = which root
                        // 2nd dimension = all the nodes for this root (might overlap after connecting trees)
  int _active_root; //which tree is actively grown now? 0 => _roots[0]
  t_vec2D_nodes _close_target_nodes; // the closest nodes in both trees

public:

  RRT_conformation_DAG()
    = default;
  void print(){
    std::cout<<"******************* "<< _nodes[0].size()<<"   "<<_nodes[1].size()<<std::endl;
   }

   int get_active_id()
   	{ return _active_root; }

   void set_active(RRT_node* root)
     { _active_root = get_root_id(root); }

   void set_active(int id)
     {
       assert(id < (int)_roots.size());
       _active_root = id;
     }

   // resets active tree, s.t. invoking set_next_active()
   // will activate the first tree
   // NOTE: after initialize active - no valid tree is active yet
   void initialize_active()
     {
       _active_root = -1;
     }

   // advances the active root
   // returns false if passed the last root, true if still iterating
   bool set_next_active()
     {
       _active_root++;
       assert(_active_root <= (int)_roots.size());
       if (_active_root == (int)_roots.size())
	 {
	   _active_root =0;
	   return false;
	 }
       _active_root = _active_root % _roots.size();
       return true; // old cycle

     }

   RRT_node* get_active_root() const{
     assert(_active_root < (int)_roots.size());
     return _roots[_active_root];
   }

   std::vector<RRT_node*>&  get_active_nodes(){
     assert(_active_root < (int)_roots.size());
     return _nodes[_active_root];
   }

   core::pose::Pose const& get_root_template(int root_id) const
     {
      debug_assert ((root_id < (int)_roots.size()) && (root_id >= 0));
       return _roots[root_id]->_template_pose;
     }

   // get ref to template pose of active root
   // (templates are used for memory saving in nodes storage)
   core::pose::Pose const& get_active_template() const
     {
       RRT_node const* root_node = get_active_root();
       return root_node->_template_pose; // can do it on private - we are friends...
     }

   // get size of tree rooted at root_id...
   int get_tree_size(int root_id)
     {
       assert(root_id < _roots.size());
       return _nodes[root_id].size();
     }

   int get_num_roots() const
     { return _roots.size(); }

   // returns the identifier of a root in the roots list, or -1 if not found
   // see also: get_root_by_id()
   int get_root_id(RRT_node* root) const
     {
       // NOTE: very naive implementation...
       int id;
       for(id = 0; id < (int)_roots.size(); id++){
	 if(_roots[id] == root) break;
       }
       if(id == (int)_roots.size())
	 return -1;
       return id;
     }

   std::vector<RRT_node*>& get_node_list_of_root(RRT_node* root){
     int id = get_root_id(root);
     return _nodes[id];
   }

   std::vector< RRT_node * >& get_roots()
     { return _roots; }

   // returns the set of all nodes whose src_root_id == root_id
   std::vector< RRT_node * >& get_root_nodes(int root_id)
     {
       assert(root_id < (int)_roots.size());
       return _nodes[root_id];
     }

   // see also get_root_id
   RRT_node* get_root_by_id(int id){
     assert(id < (int)_roots.size());
     return _roots[id];
   }


   // adds n as a child of parent
   // if parent is NULL, adds n as root
   // NOTE: if n didn't have a parent yet, also sets n's
   //       root_id to parent's root_id
   //
   // IMPLEMENTATION NOTE: pose is stored only for root, and is
   // used there as a template for all other nodes, together with
   // application of DOF vectors
   //
   // returns:
   // root_id for added node ( = n->get_orig_src_root_id()
   // after addition)
   int add_node(RRT_node* n, RRT_node* parent)
     {
       int root_id;
       if(!parent)
	 {
	   // add as root:
	   root_id = _roots.size();
	   _roots.push_back(n);
	   n->_level = 0;
	   n->_orig_src_root_id = root_id; // self
	   n->_owner_DAG = this;
	   std::vector< RRT_node * > tree_nodes; // list for tracking n's children
	   _nodes.push_back(tree_nodes);

	 }
       else
	 {
	   parent->add_child(n);
	   root_id = n->get_orig_src_root_id();
	 }
       // track n in its root children list
       assert(root_id < (int)_roots.size());
       if(!n->_is_added){
           _nodes[root_id].push_back(n);
       }
       n->_is_added = true;

       return root_id;
     }

};


} // namespace pathways

#endif
