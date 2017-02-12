// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/FreePeptide.hh
/// @brief FreePeptide represents a free peptide
/// @author xingjiepan (xingjiepan@gmail.com)


#ifndef INCLUDED_protocols_backbone_moves_local_backbone_mover_FreePeptide_hh
#define INCLUDED_protocols_backbone_moves_local_backbone_mover_FreePeptide_hh

// Project headers
#include <protocols/backbone_moves/local_backbone_mover/FreePeptide.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// Core headers
#include <core/id/AtomID.hh>
#include <core/pose/Pose.fwd.hh>

// Protocol headers
#include <protocols/backbone_moves/local_backbone_mover/types.hh>


namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {

///@brief Residue is a helper class for FreePeptide
class Residue : public utility::pointer::ReferenceCount {

public:

	Residue(Size seqpos, 
		Real n_ca_bond, Real ca_c_bond, Real c_n_bond,
		Real n_ca_c_angle, Real ca_c_n_angle, Real c_n_ca_angle,
		Real phi, Real psi, Real omega) :
		seqpos_(seqpos), 
		n_ca_bond_(n_ca_bond), ca_c_bond_(ca_c_bond), c_n_bond_(c_n_bond),
		n_ca_c_angle_(n_ca_c_angle), ca_c_n_angle_(ca_c_n_angle), c_n_ca_angle_(c_n_ca_angle),
		phi_(phi), psi_(psi), omega_(omega)
	{}

	// Constructor from pose
	Residue(Size seqpos, core::pose::Pose const& pose);

	// Apply internal coordinates of this residue to the pose
	void apply_to_pose(core::pose::Pose &pose);

	void seqpos( Size seqpos ){
		seqpos_ = seqpos;
	}

	Size seqpos(){
		return seqpos_;
	}

	void n_ca_bond( Real n_ca_bond ){
		n_ca_bond_ = n_ca_bond;
	}

	Real n_ca_bond(){
		return n_ca_bond_;
	}

	void ca_c_bond( Real ca_c_bond ){
		ca_c_bond_ = ca_c_bond;
	}

	Real ca_c_bond(){
		return ca_c_bond_;
	}

	void c_n_bond( Real c_n_bond ){
		c_n_bond_ = c_n_bond;
	}

	Real c_n_bond(){
		return c_n_bond_;
	}

	void n_ca_c_angle( Real n_ca_c_angle ){
		n_ca_c_angle_ = n_ca_c_angle;
	}

	Real n_ca_c_angle(){
		return n_ca_c_angle_;
	}

	void ca_c_n_angle( Real ca_c_n_angle ){
		ca_c_n_angle_ = ca_c_n_angle;
	}

	Real ca_c_n_angle(){
		return ca_c_n_angle_;
	}

	void c_n_ca_angle( Real c_n_ca_angle ){
		c_n_ca_angle_ = c_n_ca_angle;
	}

	Real c_n_ca_angle(){
		return c_n_ca_angle_;
	}

	void phi( Real phi ){
		phi_ = phi;
	}

	Real phi(){
		return phi_;
	}

	void psi( Real psi ){
		psi_ = psi;
	}

	Real psi(){
		return psi_;
	}

	void omega( Real omega ){
		omega_ = omega;
	}

	Real omega(){
		return omega_;
	}

	void n_xyz( xyzVector <Real> n_xyz){
		n_xyz_ = n_xyz;
	}

	xyzVector <Real> n_xyz(){
		return n_xyz_;
	}

	void ca_xyz( xyzVector <Real> ca_xyz){
		ca_xyz_ = ca_xyz;
	}

	xyzVector <Real> ca_xyz(){
		return ca_xyz_;
	}

	void c_xyz( xyzVector <Real> c_xyz){
		c_xyz_ = c_xyz;
	}

	xyzVector <Real> c_xyz(){
		return c_xyz_;
	}

private:

	Size seqpos_;
	
	Real n_ca_bond_;
	Real ca_c_bond_;
	Real c_n_bond_;

	Real n_ca_c_angle_;
	Real ca_c_n_angle_;
	Real c_n_ca_angle_;

	Real phi_;
	Real psi_;
	Real omega_;

	xyzVector <Real> n_xyz_;
	xyzVector <Real> ca_xyz_;
	xyzVector <Real> c_xyz_;
};


///@brief FreePeptide represents a free peptide
class FreePeptide : public utility::pointer::ReferenceCount {

	friend class ::LocalBackboneMoverTests;

public: // Functions from code template generator

	FreePeptide(core::pose::Pose const& pose, Size pivot1, Size pivot2);
	FreePeptide(FreePeptide const &) = default;

	virtual ~FreePeptide();

	FreePeptideOP clone() const;

public:

	/// @brief Total number of atoms
	/// @detail Only consider atoms in the free peptide,
	/// The N of the first residue and C of the last residue
	/// are excluded
	Size number_of_atoms(){
		return 3 * (pivot2_ - pivot1_ - 1) - 2;
	}

	/// @brief Total number of free peptide residues.
	Size number_of_residues(){
		return pivot2_ - pivot1_ - 1;
	}

	/// @brief Getter of the first pivot
	Size pivot1(){
		return pivot1_;
	}

	/// @brief Getter of the second pivot
	Size pivot2(){
		return pivot2_;
	}

	/// @brief Getter of the n_ca_bond
	Real n_ca_bond(Size seqpos){
		return residues_[res_id(seqpos)].n_ca_bond();
	}

	/// @brief Setter of the n_ca_bond
  void n_ca_bond(Size seqpos, Real value){
		residues_[res_id(seqpos)].n_ca_bond(value);
		xyz_updated_ = false;
	}

	/// @brief Getter of the ca_c_bond
	Real ca_c_bond(Size seqpos){
		return residues_[res_id(seqpos)].ca_c_bond();
	}

	/// @brief Setter of the ca_c_bond
  void ca_c_bond(Size seqpos, Real value){
		residues_[res_id(seqpos)].ca_c_bond(value);
		xyz_updated_ = false;
	}

	/// @brief Getter of the c_n_bond
	Real c_n_bond(Size seqpos){
		return residues_[res_id(seqpos)].c_n_bond();
	}

	/// @brief Setter of the c_n_bond
  void c_n_bond(Size seqpos, Real value){
		residues_[res_id(seqpos)].c_n_bond(value);
		xyz_updated_ = false;
	}

	/// @brief Getter of the n_ca_c_angle
	Real n_ca_c_angle(Size seqpos){
		return residues_[res_id(seqpos)].n_ca_c_angle();
	}

	/// @brief Setter of the n_ca_c_angle
  void n_ca_c_angle(Size seqpos, Real value){
		residues_[res_id(seqpos)].n_ca_c_angle(value);
		xyz_updated_ = false;
	}

	/// @brief Getter of the ca_c_n_angle
	Real ca_c_n_angle(Size seqpos){
		return residues_[res_id(seqpos)].ca_c_n_angle();
	}

	/// @brief Setter of the ca_c_n_angle
  void ca_c_n_angle(Size seqpos, Real value){
		residues_[res_id(seqpos)].ca_c_n_angle(value);
		xyz_updated_ = false;
	}

	/// @brief Getter of the c_n_ca_angle
	Real c_n_ca_angle(Size seqpos){
		return residues_[res_id(seqpos)].c_n_ca_angle();
	}

	/// @brief Setter of the c_n_ca_angle
  void c_n_ca_angle(Size seqpos, Real value){
		residues_[res_id(seqpos)].c_n_ca_angle(value);
		xyz_updated_ = false;
	}

	/// @brief Getter of the phi
	Real phi(Size seqpos){
		return residues_[res_id(seqpos)].phi();
	}

	/// @brief Setter of the phi
  void phi(Size seqpos, Real value){
		residues_[res_id(seqpos)].phi(value);
		xyz_updated_ = false;
	}

	/// @brief Getter of the psi
	Real psi(Size seqpos){
		return residues_[res_id(seqpos)].psi();
	}

	/// @brief Setter of the psi
  void psi(Size seqpos, Real value){
		residues_[res_id(seqpos)].psi(value);
		xyz_updated_ = false;
	}

	/// @brief Getter of the omega
	Real omega(Size seqpos){
		return residues_[res_id(seqpos)].omega();
	}

	/// @brief Setter of the omega
  void omega(Size seqpos, Real value){
		residues_[res_id(seqpos)].omega(value);
		xyz_updated_ = false;
	}

	/// @brief Getter of the xyz coordinates of N
	xyzVector <Real> n_xyz(Size seqpos){
		if(!xyz_updated_){ update_xyz_coords(); }
		return residues_[res_id(seqpos)].n_xyz();
	}

	/// @brief Getter of the xyz coordinates of CA
	xyzVector <Real> ca_xyz(Size seqpos){
		if(!xyz_updated_){ update_xyz_coords(); }
		return residues_[res_id(seqpos)].ca_xyz();
	}

	/// @brief Getter of the xyz coordinates of C
	xyzVector <Real> c_xyz(Size seqpos){
		if(!xyz_updated_){ update_xyz_coords(); }
		return residues_[res_id(seqpos)].c_xyz();
	}

	/// @brief Apply the coordinates of the free peptide to a pose.
	/// @detail The internal coordinates from pivot1 - 1 to pivot2 + 1
	/// will be apply to the pose. Torsions in the two gaps should be
	/// further reset to close the segment.
	void apply_to_pose(core::pose::Pose &pose);

	/// @brief Update the xyz coordinates of the residues
	void update_xyz_coords();

	/// @brief Translate the free peptide by the vector v
	void translate(xyzVector <Real> v);

	/// @brief Rotate the free peptide by the rotation matrix M
	/// @detail The center of rotation is the CA atom of the 
	/// first residue of the free peptide
	void rotate(xyzMatrix <Real> M);

	/// @brief Align the stubs to the old stubs.
	void align();

private: // Methods

	/// @brief Translate the seqpos to the id within the free peptide
	Size res_id(Size seqpos){
		assert(seqpos > pivot1_ - 3 && seqpos < pivot2_ + 3);
		return seqpos - pivot1_ + 3;
	}

	/// @brief Translate the id of a residue within the free peptide to its seqpos
	Size seqpos(Size res_id){
		assert(res_id < pivot2_ - pivot1_ + 6);
		return res_id + pivot1_ - 3;
	}

private: // Member data
	
	Size pivot1_;	
	Size pivot2_;

	// Residues from pivot1 - 2 to pivot2 + 2
	// Those residues beyond pivots exist because they provide
	// access points for modifying bond lengths and angles within gaps.
	// Also stubs coordinates could be easily read by the GapCloser.
	vector1<Residue> residues_; // Residues within the free peptide

	// Coordinates of the CA, C of the first residue and N, CA
	// of the last residue. Used by the align() function.
	vector1 < xyzVector<Real> > old_stubs_;
	xyzVector<Real> old_stub_centor_of_mass_;

	// Reference coordinates are the coordinates of C atom of the first pivot and N, CA atoms of the first  
	// residue of the free peptide. Translation and rotation of the free peptide are
	// realized by changing these coordinates.
	vector1 <core::id::AtomID> reference_atom_ids_;
	vector1 < xyzVector<Real> > reference_coordinates_;

	// Indicator telling if the xyz coords of the residues are updated
	bool xyz_updated_ = false;

};


} //protocols
} //backbone_moves
} //local_backbone_mover



#endif //INCLUDED_protocols_backbone_moves_local_backbone_mover_FreePeptide_hh





