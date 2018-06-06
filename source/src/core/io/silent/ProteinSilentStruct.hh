// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/silent/ProteinSilentStruct.hh
///
/// @brief Representation of rosetta++ protein silent-file structures.
/// @author James Thompson, Mike Tyka

#ifndef INCLUDED_core_io_silent_ProteinSilentStruct_hh
#define INCLUDED_core_io_silent_ProteinSilentStruct_hh

// mini headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>

#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace silent {


template <class T>
class ProteinSilentStruct_Template : public SilentStruct {

	typedef SilentStruct Parent;

public:
	ProteinSilentStruct_Template(
		SilentFileOptions const & opts,
		core::pose::Pose const & pose,
		std::string tag = "empty_tag",
		bool fa = false
	) :
		Parent( opts ),
		fullatom_( fa )
	{
		bJumps_use_IntraResStub_ = false;
		symminfo_ = core::conformation::symmetry::SymmetryInfoOP( new core::conformation::symmetry::SymmetryInfo() );
		symminfo_->set_use_symmetry(false);
		fill_struct( pose, tag );
	} // ProteinSilentStruct

	/// @brief Constructors.
	ProteinSilentStruct_Template(
		SilentFileOptions const & opts,
		Size const nres_in
	) :
		Parent( opts )
	{
		nres( nres_in );
		fullatom_  = false;
		decoy_tag( "empty_tag" );
		resize( nres_in );
		//fold_tree_ = core::kinematics::FoldTree();
		bJumps_use_IntraResStub_ = false;
		symminfo_ = new core::conformation::symmetry::SymmetryInfo();
		symminfo_->set_use_symmetry(false);
	}

	ProteinSilentStruct_Template( SilentFileOptions const & opts ) :
		Parent( opts )
	{
		nres( 0 );
		decoy_tag( "empty_tag" );
		fullatom ( false );
		bJumps_use_IntraResStub_ = false;
		symminfo_ = core::conformation::symmetry::SymmetryInfoOP( new core::conformation::symmetry::SymmetryInfo() );
		symminfo_->set_use_symmetry(false);
	}

	// ProteinSilentStruct_Template(
	//  core::pose::Pose const & pose,
	//  std::string tag = "empty_tag",
	//  bool fa = false
	// );

	/// @brief Returns a new ProteinSilentStruct with a copy of the information
	/// in this ProteinSilentStruct.
	virtual SilentStructOP clone() const {
		return SilentStructOP( new ProteinSilentStruct_Template<T>( *this ) );
	}

	// destructor
	~ProteinSilentStruct_Template() {}

	/// @brief Test if this ProteinSilentStruct is equal to the given
	/// ProteinSilentStruct_Template<T> in terms of conformation. Doesn't check energies.
	ProteinSilentStruct_Template<T> & operator= (
		ProteinSilentStruct_Template<T> const & src
	);

	/// @brief Tells this ProteinSilentStruct object to initialize itself from
	//the given set of lines.
	virtual bool init_from_lines(
		utility::vector1< std::string > const & lines,
		SilentFileData & container
	);

	/// @brief Configure the conformation of the given Pose with the
	/// conformational data within this ProteinSilentStruct.
	/// Calls pose.clear() and rebuilds Pose from scratch using FA_STANDARD
	/// residue types.
	virtual void fill_pose(
		core::pose::Pose & pose,
		bool const metapatches = true
	) const;

	/// @brief Configure the conformation of the given Pose with the
	/// conformational data within this ProteinSilentStruct. Calls pose.clear()
	/// and rebuilds Pose from scratch using the user-specified residue types.
	virtual void fill_pose(
		core::pose::Pose & pose,
		core::chemical::ResidueTypeSet const & residue_set,
		bool const metapatches = true
	) const;

	/// @brief opposite of fill_pose
	virtual void fill_struct(
		core::pose::Pose const & pose,
		std::string tag = "empty_tag"
	);

	/// @brief Prints the conformation information within this
	// ProteinSilentStruct to the given std::ostream.
	virtual void print_conformation( std::ostream & output ) const;

	/// @brief data getters/setters
	bool fullatom() const {
		return fullatom_;
	}

	virtual void fullatom( bool fullatom ) {
		fullatom_ = fullatom;
	}

	Real phi( Size const seqpos ) const {
		return phi_[seqpos];
	}

	Real psi( Size const seqpos ) const {
		return psi_[seqpos];
	}

	Real omega( Size const seqpos ) const {
		return omega_[seqpos];
	}

	char secstruct( Size const seqpos ) const {
		return secstruct_[seqpos];
	}

	// not safe, therefore deprecated.

	Real chi( Size const seqpos, Size const chi_num ) const;

	/// @brief returns the number of chis at this position.
	Size n_chi( Size const seqpos ) const;

	numeric::xyzVector<T> const & coords( Size const seqpos ) const {
		return coords_[seqpos];
	}

	utility::vector1< numeric::xyzVector<T> > const & coords() const {
		return coords_;
	}

	void phi( Size const seqpos, Real const phi ) {
		phi_[seqpos] = phi;
	}

	void psi( Size const seqpos, Real const psi ) {
		psi_[seqpos] = psi;
	}

	void omega( Size const seqpos, Real const omega ) {
		omega_[seqpos] = omega;
	}

	void secstruct( Size const seqpos, char const ss ) {
		secstruct_[seqpos] = ss;
	}

	utility::vector1< Size > const & chain_endings() const {
		return chain_endings_;
	}

	/// @brief set the list of chain endings
	/// @remarks All positions in the list must be strictly less than the
	/// number of residues in the data.  If this condition is not met the
	/// routine will fail-fast, so remember to resize() properly prior to
	/// calling this function.
	void chain_endings( utility::vector1< Size > const & endings );

	/// @brief add a chain ending to the list
	/// @remarks All positions in the list must be strictly less than the
	///  number of residues in the data.  If this condition is not met
	///  the routine will fail-fast, so remember to resize() properly prior
	///  to calling this function.
	void add_chain_ending( Size const seqpos );

	void chi( Size const seqpos, utility::vector1< T > const & chis );

	void chi( Size const seqpos, Size const chi_num, Real const chi );

	Size max_chi() const;

	void coords( Size const seqpos,  numeric::xyzVector<T> const & coords ) {
		coords_[seqpos] = coords;
	}

	void fold_tree( kinematics::FoldTree const & f ) {
		fold_tree_ = f;
	}

	kinematics::FoldTree const& fold_tree( ) const {
		return fold_tree_;
	}

	//lin Symmetry
	// @lin - move these to the .cc file so you can only include SymmetryInfo.fwd.hh!
	bool is_symmetric() const { return symminfo_->get_use_symmetry(); }

	void symmetry_info( core::conformation::symmetry::SymmetryInfo & s ) {
		symminfo_ = core::conformation::symmetry::SymmetryInfoOP( new core::conformation::symmetry::SymmetryInfo( s ) );
	}

	core::conformation::symmetry::SymmetryInfoCOP symmetry_info( ) const {
		return symminfo_;
	}

	void add_jump( kinematics::Jump const & jump ) {
		jumps_.push_back( jump.rt() );
	}

	void add_rt( kinematics::RT const & rt ) {
		jumps_.push_back( rt );
	}

	/// @brief returns the number of jumps held in this container.
	Size njumps() const {
		return jumps_.size();
	}

	// it's really odd that this function is called jump, but returns an RT.
	kinematics::RT const & jump( Size const jump_num ) const {
		return jumps_[ jump_num ];
	}

	// @brief returns the positions of the CA atoms in this
	// ProteinSilentStruct. Useful for RMS calculations.
	virtual ObjexxFCL::FArray2D< Real > get_CA_xyz() const;

	// model quality-related methods.
	virtual Real CA_rmsd( ProteinSilentStruct_Template<T> other_pss );

	/// @brief calculates the RMSD between the C-alpha atoms of a Pose built from
	/// the torsions in this ProteinSilentStruct and the C-alpha atoms from this
	/// ProteinSilentStruct.
	virtual Real get_debug_rmsd();

	static bool is_single_precision();

protected:
	/// @brief parse the chain endings string from an input stream
	void parse_chain_endings( std::istream & stream );

	/// @brief return the chain endings string
	std::string chain_endings_str() const;

private: // private member functions
	/// @brief Re-dimension the storage capacity of this ProteinSilentStruct to
	/// the given number of residues.
	void resize( Size const nres_in );

	void resize_chi();
public:
	virtual core::Size mem_footprint() const;
protected:
	const static Size max_chi_ = 4; // maximum number of chis for the classic rosetta++ silent-file format
	bool fullatom_;

	typename utility::vector1< char > secstruct_;
	typename utility::vector1<  T > phi_;
	typename utility::vector1<  T > psi_;
	typename utility::vector1<  T > omega_;
	typename utility::vector1<  numeric::xyzVector < T > > coords_;
	utility::vector1< kinematics::RT > jumps_;
	bool bJumps_use_IntraResStub_;
	kinematics::FoldTree fold_tree_;
	core::conformation::symmetry::SymmetryInfoOP symminfo_;
	utility::vector1< Size > chain_endings_;

private:
	typename utility::vector1< utility::vector1< T > > chi_;

}; // class ProteinSilentStruct_Template


} // namespace silent
} // namespace io
} // namespace core

#endif

// I will be removing this #include in not too long; if you need to use
// this class, you'll also have to #include the .tmpl.hh file.
// #include <core/io/silent/ProteinSilentStruct.tmpl.hh>
