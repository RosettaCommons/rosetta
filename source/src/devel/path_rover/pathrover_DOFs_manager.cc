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


#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <numeric/angle.functions.hh>
// # #include "symmetry_info.h"
#include "pathrover_DOFs_manager.hh"

#include <set>

namespace pathways
{

/**********************************************************/
/********************* Pdbres_id struct *******************/
/**********************************************************/

// ctr
Pdbres_id::Pdbres_id
(char chain, int res_id, char insert_letter)
{
  this->chain = chain;
  this->res_id = res_id;
  this->insert_letter = insert_letter;
}

// straightforward less operator
bool Pdbres_id::operator<(Pdbres_id const& other) const
{
  if(chain != other.chain)
    return chain < other.chain;
  if(res_id != other.res_id)
    return res_id < other.res_id;
  return insert_letter < other.insert_letter;
}

// straightforward less operator
bool Pdbres_id::operator!=(Pdbres_id const& other) const
{
  if(chain != other.chain)
    return true;
  if(res_id != other.res_id)
    return true;
  return insert_letter != other.insert_letter;
}

std::ostream& operator<<(std::ostream& os, Pdbres_id const& res){
  os << res.chain << res.res_id << res.insert_letter;
  return os;
}

/**********************************************************/
/********************* DOF_info struct *******************/
/**********************************************************/


DOF_info::DOF_info
	(int poseres_id, double std_dev, double max_uni_dev)
{
  this->poseres_id = poseres_id;
  this->std_dev = std_dev;
  this->max_uni_dev = max_uni_dev;
}

bool DOF_info::operator<(DOF_info const& other) const
{
  return (poseres_id  <  other.poseres_id);
}


/*******************************************************/
/********************* DOFs_manager ********************/
/*******************************************************/

// ctr
DOFs_manager::DOFs_manager(core::pose::Pose& template_pose, core::kinematics::MoveMapCOP in_movemap)
{
  this->initialize(template_pose, in_movemap);
}

// initialization of manager so it can be used based on a certain pose
void
DOFs_manager::initialize
(core::pose::Pose& template_pose, core::kinematics::MoveMapCOP in_movemap)
{
  // clear arrays:
  _is_initialized = false; // until succesful finish
  if(_movemap)
    _movemap->clear();
  else
    _movemap=in_movemap->clone();
  _free_phi_list.clear();
  _free_psi_list.clear();
  _map_pdbres_to_pose.clear();
  this->_n_free_dofs = 0;
  this->_template_pose = template_pose;
  if(in_movemap)
    import_flexible_dofs_from_movemap_torsionids(in_movemap);
  // build PDB ID indexing (and perhaps DOFs lists) based on template pose
  build_pdbres_to_pose_mapping();
  _is_initialized = true;
}

// see .hh for more info
// Following functions set phi/psi/both("residue") DOF as moveable
// (including all symmetry units)
void
DOFs_manager::set_phi_dof
(char chain, int pdbres, double std_dev, double max_uni_dev)
{
  assert(_is_initialized);
  Pdbres_id pdbres_id (chain, pdbres);
  set_phi_dof(pdbres_id, std_dev, max_uni_dev);
}

// see .hh for doc
// NOTE: no need to touch symmetry - this is handled by pose!
void DOFs_manager::set_phi_dof
	(Pdbres_id pdbres_id, double std_dev, double max_uni_dev)
{
  using namespace core::kinematics;
  using namespace core::id;
  assert(_is_initialized);
  int poseres_id = pdbres_to_poseres(pdbres_id);
  std::cout << "Set_phi_dof: [Poseres id = " << poseres_id << "], for residue "
	    << pdbres_id.chain << pdbres_id.res_id << std::endl;
  assert(poseres_id != -1);
  DOF_info dof_info(poseres_id, std_dev, max_uni_dev);
  if(this->_free_phi_list.find(dof_info) != this->_free_phi_list.end())
    return;
  this->_free_phi_list.insert(dof_info);
  _movemap->set( TorsionID( poseres_id, BB, 1 ), true); // TODO: symmetry?
  this->update_dofs_number();
}

// see .h for doc
void DOFs_manager::set_psi_dof
	(char chain, int pdbres, double std_dev, double max_uni_dev)
{
  assert(_is_initialized);
  Pdbres_id pdbres_id (chain, pdbres);
  set_psi_dof(pdbres_id, std_dev, max_uni_dev);
}

// see .h for doc
// NOTE: no need to touch symmetry - this is handled by pose!
void DOFs_manager::set_psi_dof
	(Pdbres_id pdbres_id, double std_dev, double max_uni_dev)
{
  using namespace core::id;
  assert(_is_initialized);
  int poseres_id = pdbres_to_poseres(pdbres_id);
  assert(poseres_id != -1);
  DOF_info dof_info(poseres_id, std_dev, max_uni_dev);
  if(this->_free_psi_list.find(dof_info) != this->_free_psi_list.end())
    return;
  this->_free_psi_list.insert(dof_info);
  _movemap->set( TorsionID( poseres_id, BB, 2 ), true); // TODO: symmetry?
  this->update_dofs_number();
}

// see .h for doc
void DOFs_manager::set_residue_dofs
	(char chain, int pdbres, double std_dev, double max_uni_dev) // both phi and psi, = 2 DOFs
{
	assert(_is_initialized);
	Pdbres_id pdbres_id (chain, pdbres);
	set_residue_dofs(pdbres_id, std_dev, max_uni_dev);
}

// see .h for doc
void DOFs_manager::set_residue_dofs
	(Pdbres_id pdbres_id, double std_dev, double max_uni_dev) // both phi and psi, = 2 DOFs
{
	assert(_is_initialized);
	set_phi_dof(pdbres_id, std_dev, max_uni_dev);
	set_psi_dof(pdbres_id, std_dev, max_uni_dev);
}

// apply uniform sampling on all DOFs
// with [-+ max_dev] deviation from current value for each
// TODO: for rotations etc., max_dev is a matrix?
void DOFs_manager::apply_uniform_sample_all(core::pose::Pose& pose) // apply uniform sampling on all DOFs for pose
{
  using namespace std;
  bool local_debug = false;//true;
  assert(_is_initialized);
  // iteratoe over phi dofs:
  std::set<DOF_info>::const_iterator iter;
  for(iter = _free_phi_list.begin();	iter != _free_phi_list.end(); iter++)
    {
    double new_val = numeric::nonnegative_principal_angle_degrees(
	iter->uniform_sample(pose.phi(iter->poseres_id)) );
      if(local_debug)
	cout << "applying phi sampling to phi(" << iter->poseres_id << ")"
	     << "from " << pose.phi(iter->poseres_id)
	     << "to " << new_val << std::endl;
      pose.set_phi(iter->poseres_id, new_val);
    }
  // iterate over psi dofs:
  for(iter = _free_psi_list.begin(); iter != _free_psi_list.end(); iter++)
    {
      if(local_debug)
	cout << "applying psi sampling to psi(" << iter->poseres_id << ")" << endl;
    double new_val = numeric::nonnegative_principal_angle_degrees(
	iter->uniform_sample(pose.psi(iter->poseres_id)) );
      pose.set_psi(iter->poseres_id, new_val);
    }
}

// apply uniform sampling on one DOF only
// with [-+ max_dev] deviation from current value
// TODO: for rotations etc., max_dev is a matrix?
// TODO: very naive implementation, imporvement requires switching from set data struct to vector or map
void DOFs_manager::apply_uniform_sample_random_DOF(core::pose::Pose& pose)
{
  assert(_is_initialized);
  std::set<DOF_info>::const_iterator iter;
  // choose DOF randomly
  unsigned int dof_index = int(this->_n_free_dofs * numeric::random::uniform()) + 1;
  // see if phi or psi:
  if(dof_index <= _free_phi_list.size()){ // phi dof
    // very naive finding of Kth-element // TODO: better be changed - inefficient
    iter = _free_phi_list.begin();
    for( unsigned int i=1; i < dof_index; i++)
      iter++;
    double new_val = numeric::nonnegative_principal_angle_degrees(
      iter->uniform_sample(pose.phi(iter->poseres_id)) );
    pose.set_phi(iter->poseres_id, new_val);
  }else { // psi dof
    dof_index -= _free_phi_list.size(); // reindex in psi list
    // very naive finding of Kth-element // TODO: better be changed - inefficient
    iter = _free_psi_list.begin();
    for(unsigned int i=1; i < dof_index; i++)
      iter++;
    double new_val = numeric::nonnegative_principal_angle_degrees(
      iter->uniform_sample(pose.psi(iter->poseres_id)) );
    pose.set_psi(iter->poseres_id, new_val);
  }
}

// apply gaussian sampling on all DOFs,
// using variancs of DOF, and the current value as mean
// (i.e. Norm(current value, var)
void DOFs_manager::apply_gaussian_sample_all(core::pose::Pose& pose)
{
  assert(_is_initialized);
  // iteratoe over phi dofs:
  std::set<DOF_info>::const_iterator iter;
  for(iter = _free_phi_list.begin(); iter != _free_phi_list.end(); iter++) {
    double new_val = numeric::nonnegative_principal_angle_degrees(
      iter->gaussian_sample(pose.phi(iter->poseres_id)) );
    pose.set_phi(iter->poseres_id, new_val);
  }
  // iterate over psi dofs:
  for(iter = _free_psi_list.begin(); iter != _free_psi_list.end(); iter++) {
    double new_val = numeric::nonnegative_principal_angle_degrees(
      iter->gaussian_sample(pose.psi(iter->poseres_id)) );
    pose.set_psi(iter->poseres_id, new_val);
  }
}

// apply gaussian sampling on one DOF only,
// using variancs of DOF, and the current value as mean
// (i.e. Norm(current value, var)
void DOFs_manager::apply_gaussian_sample_random_DOF(core::pose::Pose& pose)
{
  assert(_is_initialized);
  std::set<DOF_info>::const_iterator iter;
  // choose DOF randomly
  unsigned int dof_index = int(this->_n_free_dofs * numeric::random::uniform()) + 1;
  // see if phi or psi:
  if(dof_index <= _free_phi_list.size()){ // phi dof
    // very naive finding of Kth-element // TODO: better be changed - inefficient
    iter = _free_phi_list.begin();
    for(unsigned int i=1; i < dof_index; i++)
      iter++;
    double new_val = numeric::nonnegative_principal_angle_degrees(
      iter->uniform_sample(pose.phi(iter->poseres_id)) );
    pose.set_phi(iter->poseres_id, new_val);
  }else { // psi dof
    dof_index -= _free_phi_list.size(); // reindex in psi list
    // very naive finding of Kth-element // TODO: better be changed - inefficient
    iter = _free_psi_list.begin();
    for(unsigned int i=1; i < dof_index; i++)
      iter++;
    double new_val = numeric::nonnegative_principal_angle_degrees(
      iter->uniform_sample(pose.psi(iter->poseres_id)) );
    pose.set_psi(iter->poseres_id, new_val);
  }
}

// get a vector of all DOFs
//
// NOTE: apply careful usage!
// Vector is according to current DOFs defined for DOFs_manager -
// if using apply_dofs_values_vector() later on, make sure the state of
// the DOFs manager didn't change between calls!
std::vector<double> const
DOFs_manager::get_dofs_values_vector(core::pose::Pose const& pose) const
{
  assert(_is_initialized);
  std::vector<double> return_vector(this->_n_free_dofs);
  int k = 0;
  std::set<DOF_info>::const_iterator iter;
  // iterate over phi dofs:
  for(iter = _free_phi_list.begin(); iter != _free_phi_list.end(); iter++) {
    return_vector[k++] = pose.phi(iter->poseres_id);
  }
  for(iter = _free_psi_list.begin(); iter != _free_psi_list.end(); iter++) {
    return_vector[k++] = pose.psi(iter->poseres_id);
  }
  return return_vector; // TODO: efficiency? copies the vectror... can use an output byref parameter instead
}

// changes the free DOFs according to the input vector
// (same format and size as the one returned by get_dofs_vector)
//
// IMPORTANT NOTE: see get_dofs_vector()
void
DOFs_manager::apply_dofs_values_vector(core::pose::Pose& pose, std::vector<double> const& values_vector)
{
  assert(_is_initialized);
  int k = 0;
  std::set<DOF_info>::const_iterator iter;
  // iterate over phi dofs:
  for(iter = _free_phi_list.begin(); iter != _free_phi_list.end(); iter++) {
    pose.set_phi(iter->poseres_id, values_vector[k++]);
  }
  for(iter = _free_psi_list.begin(); iter != _free_psi_list.end(); iter++) {
    pose.set_psi(iter->poseres_id, values_vector[k++]);
  }
} // TODO: not implemented!!! TODO: Is this still true?

// Looks up index in pose of a residue, according to PDB indexing
// i.e. chain, res_id and insertion letter
int DOFs_manager::pdbres_to_poseres(char chain, int pdb_res, char insert_letter)
{
  assert(_is_initialized);
  // look up in map
  Pdbres_id pdbres_id(chain, pdb_res, insert_letter);
  return pdbres_to_poseres(pdbres_id);
}

// Looks up index in pose of a residue, according to PDB indexing
//
// Params:
//  * pdbres_id - pdb residue identifier
//
// Returns:
//  Index of specified residue inside this->_pose object
//  i.e. reference pose as specified in constructor.
//  Return -1 if didn't find.
//
// Implementation details: use map
int DOFs_manager::pdbres_to_poseres(Pdbres_id pdbres_id)
{
  assert(_is_initialized);
  _template_pose->pdb_info()->pdb2pose
    (pdbres_id.chain, pdbres_id.res_id, pdbres_id.insert_letter);
  t_map_pdbres_to_pose::const_iterator find_iter
    = this->_map_pdbres_to_pose.find(pdbres_id);
  if(find_iter == this->_map_pdbres_to_pose.end())
    return -1; // not found
  return find_iter->second;
}


/********************************************************************************/
/************************ DOFs_manager: Private methods *************************/
/********************************************************************************/

// Builds a mapping from PDB residue indexing to pose indexing.
// This is the inverse of pose.pdb_info() mapping.
// (does not require pre-initialization)
void DOFs_manager::build_pdbres_to_pose_mapping()
{
  core::pose::PDBInfoCOP pdb_info = _template_pose->pdb_info();
  Pdbres_id pdbres;
  for(core::Size i=1; i < _template_pose->total_residue(); i++)
    {
      pdbres.chain = pdb_info->chain(i);
      pdbres.res_id = pdb_info->number(i);
      pdbres.insert_letter = pdb_info->icode(i);
      assert(_map_pdbres_to_pose.find(pdbres) == _map_pdbres_to_pose.end()); // TODO: can we assume no too pdb_ids are identical?
      _map_pdbres_to_pose[pdbres] = i;
    }
}

// updates this object DOFs according to movemap of backbone torsion ids
// (does not require pre-initialization)
//
// IMPLEMENTATION NOTE:
// Symmetry is handled implicitly, since if it was previously incorporated into pose, then DOFs are symmetric as well
void DOFs_manager::import_flexible_dofs_from_movemap_torsionids
(core::kinematics::MoveMapCOP new_movemap,
	 double default_var,
	 double default_max_dev)
{
  using namespace core::kinematics;
  using namespace core::id;

  _movemap = new_movemap->clone();
  // prepate a template with correct var & max_dev
  DOF_info template_dof(0, default_var, default_max_dev);
  // TorsionID
  for ( MoveMap::TorsionID_Map::const_iterator
	  i = new_movemap->torsion_id_begin() ;
	i != new_movemap->torsion_id_end();
	++i )
    {
      if(!i->second) continue;
      template_dof.poseres_id = i->first.rsd();
      if(i->first.torsion() == 1 && i->first.type() == BB)
	this->_free_phi_list.insert(template_dof);
      if(i->first.torsion() == 2 && i->first.type() == BB)
	this->_free_psi_list.insert(template_dof);
    }
}

} // namespace pathways
