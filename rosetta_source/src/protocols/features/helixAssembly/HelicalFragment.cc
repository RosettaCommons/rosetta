// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /rosetta/rosetta_source/src/protocols/features/helixAssembly/HelicalFragment.cc
///
/// @brief

/// @author Tim jacobs

//Unit Headers
#include <protocols/features/helixAssembly/HelicalFragment.hh>

//Core
#include <core/conformation/Residue.hh>

//Utility
#include <utility/string_util.hh>

//Numeric
#include <numeric/xyzVector.hh>
#include <numeric/HomogeneousTransform.hh>

namespace protocols {
namespace features {
namespace helixAssembly {

HelicalFragment::HelicalFragment(){}

HelicalFragment::HelicalFragment(core::Size start, core::Size end):
start_(start),
end_(end),
pdb_source_(""),
direction_(true)
//residue_list_(end+1)
//residue_map_()
{}

HelicalFragment::HelicalFragment(core::Size start, core::Size end, bool direction):
start_(start),
end_(end),
pdb_source_(""),
direction_(direction)
//residue_list_(end+1)
//residue_map_()
{}

HelicalFragment::~HelicalFragment(){}

core::Size HelicalFragment::get_start() const
{
    return start_;
}

core::Size HelicalFragment::get_end() const
{
    return end_;
}

core::Size HelicalFragment::get_size() const
{
    if(end_ == start_){return 0;}
    return end_-start_+1;
}

std::string HelicalFragment::get_pdb_source() const
{
    return pdb_source_;
}

bool HelicalFragment::get_direction() const
{
    return direction_;
}

void HelicalFragment::set_pdb_source(std::string pdb_source_)
{
    this->pdb_source_ = pdb_source_;
}

void HelicalFragment::set_direction(bool direction_)
{
    this->direction_ = direction_;
}

//@Brief
//void HelicalFragment::insertResiduesFromPose(const core::pose::Pose & pose, core::Size start, core::Size end){
//
//  //ensure that we are adding as many residues from the pose as we have elements in this HelicalFragment
//  assert(get_size()==(end-start+1));
//
//  core::Size foo = start_;
//  for(core::Size i=start; i <= end; ++i){
//      core::conformation::Residue ros_res = pose.residue(i);
//
//      std::vector<NativeAtom> atoms;
//      for(core::Size j=1; j<=ros_res.atoms().size(); ++j){
//          numeric::xyzVector<double> ros_atom_xyz = ros_res.atom(j).xyz();
//          NativeAtom atom(ros_atom_xyz.x(),ros_atom_xyz.y(),ros_atom_xyz.z());
//          atoms.push_back(atom);
//      }
//      NativeResidue native_res(ros_res.name(), atoms);
//      residue_map_[foo].push_back(native_res);
//      ++foo;
//  }
//}

//std::string HelicalFragment::print() const{
//
//  std::string output = "Fragment(" + utility::to_string(start_) + "," + utility::to_string(end_) + ")\n";
//
//  std::map<core::Size, std::vector<NativeResidue> >::const_iterator it;
//  for(it=residue_map_.begin(); it != residue_map_.end(); ++it){
//      output += "Residue " + utility::to_string((*it).first) + "\n";
//
//      for(core::Size i=0; i<(*it).second.size(); ++i){
//          output += (*it).second[i].print();
//      }
//  }
//  return output;
//}

//void HelicalFragment::insertResiduesFromPose(const core::pose::Pose & matching_pose,
//		core::Size start, core::Size end, const core::pose::Pose & bundle_pose){
//
//  //ensure that we are adding as many residues from the pose as we have elements in this HelicalFragment
//  assert(get_size()==(end-start+1));
//
//  core::Size counter = start_;
//  for(core::Size i=start; i <= end; ++i){
//      core::conformation::Residue matching_res = matching_pose.residue(i);
//      core::conformation::Residue bundle_res = bundle_pose.residue(counter);
//
//      numeric::xyzVector< core::Real > halfpoint_input = 0.5 * (matching_res.atom("N").xyz() + matching_res.atom("C").xyz() );
//      numeric::HomogeneousTransform< core::Real > input_frame( matching_res.atom("N").xyz(), halfpoint_input, matching_res.atom("CA").xyz() );
//
//      numeric::xyzVector< core::Real > halfpoint_output = 0.5 * (bundle_res.atom("N").xyz() + bundle_res.atom("C").xyz() );
//      numeric::HomogeneousTransform< core::Real > output_frame( bundle_res.atom("N").xyz(), halfpoint_output, bundle_res.atom("CA").xyz() );
//
//      std::cout << "Working on residue: " << i << "(" << matching_res.type().name() << ")/" << counter << "(" << bundle_res.type().name() << ")" << std::endl;
//
//      std::vector<NativeAtom> atoms;
////      std::vector<NativeAtom> atoms2;
//      for(core::Size j=1; j<=matching_res.atoms().size(); ++j){
//    	  //If this atom is a backbone atom we don't want to transform the coordinates at all
//    	  if(matching_res.atom_is_backbone(j)){
//    		  if(j <= bundle_res.atoms().size()){//This should only happen with res type derivatives (like c-term etc) I think....
//				  numeric::xyzVector< core::Real > bundle_atom_xyz(bundle_res.atom(j).xyz());
//				  NativeAtom atom(j, bundle_atom_xyz.x(),bundle_atom_xyz.y(),bundle_atom_xyz.z());
//				  atoms.push_back(atom);
//    		  }
//    	  }
//    	  //If this atom is a sidechain atom we need to change the atom xyz to match be in the coordinate frame of the
//    	  //new residue
//    	  else{
//    		  //create a Homogeneous transform object in the global coordinate frame with the point set to the current atom
//    		  numeric::xyzVector< core::Real > matching_atom_xyz_transformed((output_frame * input_frame.inverse()) *
//    				  matching_res.atom(j).xyz());
//    		  NativeAtom atom(j, matching_atom_xyz_transformed.x(),
//    				  matching_atom_xyz_transformed.y(),matching_atom_xyz_transformed.z());
//    		  atoms.push_back(atom);
//    	  }
////    	  numeric::xyzVector< core::Real > ros_atom_old_xyz(ros_res_old.atom(j).xyz());
////    	  NativeAtom atom2(j, ros_atom_old_xyz.x(),ros_atom_old_xyz.y(),ros_atom_old_xyz.z());
////    	  atoms2.push_back(atom2);
//      }
//      NativeResidue native_res(matching_res.type().name(), atoms);
////      NativeResidue native_res2(ros_res_old.type().name(), atoms2);
//      residue_list_[counter].push_back(native_res);
////      residue_list_[counter].push_back(native_res2);
//      ++counter;
//  }
//  std::cout << "Done with round" << std::endl;
//}
//
//std::string HelicalFragment::print() const{
//
//  std::string output = "Fragment(" + utility::to_string(start_) + "," + utility::to_string(end_) + ")\n";
//
//  for(core::Size j=start_; j<residue_list_.size(); ++j){
//      output += "RESNUM " + utility::to_string(j) + "\n";
//
//      for(core::Size i=0; i<residue_list_[j].size(); ++i){
//          output += residue_list_[j][i].print();
//      }
//  }
//  return output;
//}
    
} //namespace helixAssembly
} //namespace features
} //namespace protocols

