////////////////////////////////////////////////////////////////////////////////
//     pathways_DOFs_manager.hh:
//     This class is for wrapping a pose object, and for managing allowed DOFs,
//     access to chains, symmetry relations etc. All algorithmic manipulation of
//     DOFs is channeled to this object, in order to add flexibility.
//
//     Originally Created in R++: 3/12/2007 (Barak Raveh & Angela Enosh)
//     Transformed to mini: 10/11/2009 (Barak Raveh)
//
// methods list:
//     TODO: fill this in
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_devel_path_rover_pathrover_DOFs_manager_hh
#define INCLUDED_devel_path_rover_pathrover_DOFs_manager_hh

#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/random/random.fwd.hh>
//#include "pathways_RRT_conformation_tree.h"
//#include "pdb.h"
//#include "score_data.h"
//#include "symmetry_info.h"
// #include <ObjexxFCL/ObjexxFCL.hh>

#include <iostream>
#include <map>
#include <set>

namespace pathways {


/********************* Pdbres_id struct *******************/
class Pdbres_id {
public:
  char chain; // ' ' = doesn't matter
  int res_id;
  char insert_letter; // ' ' = no insertion
public:
  Pdbres_id() { insert_letter = ' '; }// default
  Pdbres_id(char chain, int res_id, char insert_letter = ' ');
  bool operator<(Pdbres_id const& other) const;
  bool operator!=(Pdbres_id const& other) const;
};

std::ostream& operator<<(std::ostream& os, Pdbres_id const& res);

/********************* DOF_info struct *******************/
// this contains extended information about a DOF - amplitude of
// allowed motion, etc.
class DOF_info {
public:
  int poseres_id; // index of residue in Pose object
  double std_dev; // standard deviation for gaussian sampling
  double max_uni_dev; // maximal deviation for uniform sampling
public:
  DOF_info() = default;
  DOF_info(int poseres_id, double std_dev, double max_uni_dev);
  bool operator<(DOF_info const& other) const;

  // sample ~Uniform[mean-max_dev .. mean+max_dev]
  inline double uniform_sample(double mean) const
  {
    using namespace numeric::random;
    double rand_pert = (uniform() - 0.5) * max_uni_dev; // [-max_dev..+max_dev]
    return mean + rand_pert;
  }

  // sample ~Norm[mean, std_dev^2]
  inline double gaussian_sample(double mean) const
  {
    return mean + numeric::random::gaussian() * std_dev;
  }
};


/********************* DOFs_manager struct *******************/

//     This class is for wrapping a pose object, and for managing allowed DOFs,
//     access to chains, symmetry relations etc. All algorithmic manipulation of
//     DOFs is channeled to this object, in order to add flexibility.
class DOFs_manager {
public:
  // Typedefs
  // ========
  typedef std::map<Pdbres_id, int> t_map_pdbres_to_pose; // map from pdb indexing to pose indexing

private:
  // Class Variables:
  // ================
  core::pose::PoseCOP _template_pose; // reference pose serving as template for all samples (i.e. same residues and folding tree)
  t_map_pdbres_to_pose _map_pdbres_to_pose; // map from pdb indexing to pose indexing

  // DOFs list: each is a list of residue ids (according to pose indexing) for which the
  // relevant DOFs are considerd free to move
  // (Note: lists are for efficient access to small # of DOFs, unlike a boolean array like in Pose)
  std::set<DOF_info> _free_phi_list;
  std::set<DOF_info> _free_psi_list;
  core::kinematics::MoveMapOP _movemap;
  int _n_free_dofs; // total number of free torsions

  bool _is_initialized;

public:
  // Public Methods:
  // ===============

  // this ctr requires initialization before manager can be used
  // by calling DOFs_manager::initialize()
  DOFs_manager() : _is_initialized(false) {}

  // Params:
  // * template_pose - this is a template pose structure
  // * in_movemap - read flexible DOFs from movemap
  //
  //  The DOFs_manager object should only be applied to Pose objects
  //  That are "similar" to the template - same # of residues, chains, etc.
  // Note: constructor loads "template_pose.allow_dofs_move" arrays to initialize
  //       DOFs vectors.
  //       - "this->set_dofs_dof()" methods can be used to add more DOFs.
  //
  DOFs_manager(core::pose::Pose& template_pose, core::kinematics::MoveMapCOP in_movemap = nullptr);

  // does all constructor stuff, see ctr with same params for documentation
  // in_movemap is used to set flexible DOFs (if not NULL)
  void initialize(core::pose::Pose& template_pose, core::kinematics::MoveMapCOP in_movemap);

  // Following functions set phi/psi/both("residue") DOF as moveable
  // (including all symmetry units)
  //
  // Params:
  // * chain / pdbres/ pdbres_id: alternative ways to denote a specific PDB residue,
  //   default insertion letter is ' '
  // * std_dev - std variance for gaussian sampling of DOF
  // * max_uni_dev - maximal deviation for uniform sampling of DOF
  void set_phi_dof(char chain, int pdbres, double std_dev, double max_uni_dev = 180.0F);

  void set_phi_dof(Pdbres_id pdbres_id, double std_dev, double max_uni_dev = 180.0F);

  // (see documentation of set_phi_dof)
  void set_psi_dof(char chain, int pdbres, double std_dev, double max_uni_dev = 180.0F);

  void set_psi_dof(Pdbres_id pdbres_id, double std_dev, double max_uni_dev = 180.0F);

  // (see documentation of set_phi_dof)
  void set_residue_dofs(char chain, int pdbres, double std_dev, double max_uni_dev = 180.0F); // both phi and psi, = 2 DOFs

  void set_residue_dofs(Pdbres_id pdbres_id, double std_dev, double max_uni_dev = 180.0F); // both phi and psi, = 2 DOFs

  // apply uniform sampling on all DOFs
  // with [-+ max_dev] deviation from current value in "pose" for each DOF
  // TODO: for rotations etc., max_dev is a matrix?
  void apply_uniform_sample_all(core::pose::Pose& pose); // apply uniform sampling on all DOFs for pose

  // apply uniform sampling on one DOF only
  // with [-+ max_dev] deviation from current value in "pose"
  // TODO: for rotations etc., max_dev is a matrix?
  void apply_uniform_sample_random_DOF(core::pose::Pose& pose);

  // apply gaussian sampling on all DOFs, respective to current values of "pose"
  // using variancs of DOF, and the current "pose" value as mean
  // (i.e. Norm(current value, var)
  void apply_gaussian_sample_all(core::pose::Pose& pose); //

  // apply gaussian sampling on one DOF only
  // using variancs of DOF, and the current "pose" value as mean
  // (i.e. Norm(current value, var)
  void apply_gaussian_sample_random_DOF(core::pose::Pose& pose);

  // get a vector of all DOFs values, in the same format
  // as in get_dofs_values_vector()
  //
  // NOTE: apply careful usage!
  // Vector is according to current DOFs defined for DOFs_manager -
  // if using apply_dofs_vector() later on, make sure the state of
  // the DOFs manager didn't change between calls!
  // TODO: what to do with matrices and arrays? perhaps serialize them in a uniform fashion?
  std::vector<double> const get_dofs_values_vector(core::pose::Pose const& pose) const;

  // changes the free DOFs according to the input vector
  // (same format and size as the one returned by get_dofs_vector)
  //
  // NOTE: see get_dofs_vector()
  // TODO: what to do with matrices and arrays? perhaps serialize them in a uniform fashion?
  void apply_dofs_values_vector(core::pose::Pose& pose, std::vector<double> const& values_vector);

  // Looks up index in pose of a residue, according to PDB indexing
  // i.e. chain, res_id and insertion letter
  //
  // Params:
  //  <chain>: ' ' = NULL chain
  //  <pdb_res>: residue id as in PDB file
  //  <insert_letter>: PDB insertion letter ; ' ' = no insertion (default)
  //
  // Returns:
  //  Index of specified residue inside template pose object
  //  (-1 if didn't find).
  int pdbres_to_poseres(char chain, int pdb_res, char insert_letter = ' ');

  // Looks up index in pose of a residue, according to PDB indexing
  // i.e. chain, res_id and insertion letter
  //
  // Params:
  //  <pdbres>: identifies PDB residue [input]
  //    (' ' = NULL chain / insertion letter)
  //
  // Returns:
  //  Index of specified residue inside template pose object
  //  (-1 if didn't find).
  int pdbres_to_poseres(Pdbres_id);


private:
  // go over pose and build a full map from pdb indexing to pose,
  // i.e. inverse mapping to pose.pdb_info()
  void build_pdbres_to_pose_mapping();

  // updates this object DOFs according to this->_template_pose->allow_dofs_move arrays
  void import_flexible_dofs_from_movemap_torsionids
   (core::kinematics::MoveMapCOP new_movemap,
    double default_var = 180.0,
    double default_max_dev = 180.0);

  void update_dofs_number()
  {	this->_n_free_dofs = _free_phi_list.size() + _free_psi_list.size(); }
};


} // namespace DOFs_manager

#endif
