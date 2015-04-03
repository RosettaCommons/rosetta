// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file

#ifndef INCLUDED_protocols_toolbox_DecoySetEvaluation_hh
#define INCLUDED_protocols_toolbox_DecoySetEvaluation_hh

#include <protocols/toolbox/DecoySetEvaluation.fwd.hh>

#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2P.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <core/io/silent/SilentFileData.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {


	/// @brief const_iterator class for SilentFileData container.
template< typename const_iterator, typename SomeOP >
class OP_const_iterator {
	typedef SomeOP value_type;
	typedef value_type* pointer;
	typedef value_type& reference;
	typedef std::ptrdiff_t       difference_type;
	typedef std::bidirectional_iterator_tag iterator_category;

public:
		/// @brief empty constructor
	OP_const_iterator() {}

	/// @brief Constructor, given an iterator into the Structure_Map.
	OP_const_iterator( const_iterator s_iter ) {
		it_ = s_iter;
	}

	~OP_const_iterator() {}

	OP_const_iterator& operator=( const OP_const_iterator& src ) {
		it_ = src.it_;
		return (*this);
	}

	bool operator==( const OP_const_iterator& other ) {
		return ( it_ == other.it_ );
	}

	bool operator!=( const OP_const_iterator& other ) {
		return ( it_ != other.it_ );
	}

	OP_const_iterator& operator++() {
		it_++;
		return (*this);
	}

	OP_const_iterator& operator--() {
		it_--;
		return (*this);
	}

	value_type operator->() const {
		return *it_;
	}

	value_type operator*() const {
		return *it_;
	}

private:
	const_iterator it_; // keeps track of my place in a Structure_Map
}; // class iterator


//template< typename Container >
//dereference_iterator< Container > dereference_iterator( typename Container::const_iterator );

class DecoySetEvaluation : public utility::pointer::ReferenceCount {
public:
  DecoySetEvaluation();
	virtual ~DecoySetEvaluation();

  static void register_options();


  void reserve( core::Size n_decoys_ );

  void push_back( core::pose::Pose& pose );
	//void push_back_CA_xyz( ObjexxFCL::FArray2D< core::Real > const&, core::Size nres );
	void push_back_CA_xyz( ObjexxFCL::FArray2_double const&, core::Size nres );

	void pop_back_CA_xyz( );

	//get all CA from all structures in silent-file
	void push_back_CA_xyz_from_silent_file( core::io::silent::SilentFileData const& sfd, bool store_energies );

	template< typename SilentStructIterator >
	void push_back_CA_xyz_from_silent_file(
         core::Size n_decoys,
				 SilentStructIterator begin,
				 SilentStructIterator end,
				 bool store_energies
	);

  void rmsf( utility::vector1< core::Real >& results );

	core::Size n_decoys_max() const {
		return n_decoys_max_;
	}
  core::Size n_decoys() const {
    return n_decoys_;
  }
	core::Size n_atoms() const {
		return n_atoms_;
	}

	void set_n_atom( core::Size natoms );

	ObjexxFCL::FArray3D_double const& coords() const {
		return coords_;
	}

	void compute_distance_matrix( ObjexxFCL::FArray2D_double& ) const;

	//view on coords of structure i
	ObjexxFCL::FArray2P_double coords( core::Size i ) const {
		return ObjexxFCL::FArray2P_double( coords_( 1, 1, i ), 3, n_atoms_ );
	}

	void rmsf( ObjexxFCL::FArray1_double& result );

	core::Real rmsd(
									ObjexxFCL::FArray1_double const& weights,
									ObjexxFCL::FArray2_double& xx_ref,
									ObjexxFCL::FArray2_double& xx
	) const;
	//returns the index of the structure used as reference  (icenter)
	core::Size wRMSD( core::Real sigma2, core::Real tolerance, ObjexxFCL::FArray1_double& weights );

	void set_weights( ObjexxFCL::FArray1_double const& weights );

	void superimpose( core::Size icenter = 1 ); //take icenter decoy as reference

 	core::pose::Pose const& ref_pose() {
 		return ref_pose_;
 	}

	void clear() {
		n_decoys_max_ = 0;
		n_decoys_ = 0;
		n_atoms_ = 0;
	}

	void center_structure( core::Size i );
	void center_structure( core::Size i, ObjexxFCL::FArray1_double const& weights );
	void center_all( ObjexxFCL::FArray1_double const& weights );

	void create_dist_constraints_median( core::scoring::constraints::ConstraintSet& cst_set ) const;
	void create_dist_constraints( core::scoring::constraints::ConstraintSet& cst_set ) const;
	void create_xyz_constraints_median(
		 core::scoring::constraints::ConstraintSet& cst_set,
		 core::pose::Pose const& ref_pose,
		 core::Size root
	) const;
	//	void dump_coords() const;

	core::Size find_closest_to_average( ObjexxFCL::FArray2_double& average_structure ) const;
	void compute_average_structure( ObjexxFCL::FArray2_double& average_structure ) const;

	utility::vector1< core::Real > const& all_energies() const { return all_energies_; };
	void set_all_energies( utility::vector1< core::Real > const& all_energies ) {
		store_energies_ = true;
		all_energies_ = all_energies;
	};

private:
	void superimpose( ObjexxFCL::FArray1_double const& weights, core::Size icenter = 1 );
	core::Real rmsf( core::Size pos );

	void prepare_push_back( core::Size nres );

private:
	core::pose::Pose ref_pose_;
	ObjexxFCL::FArray1D_double COM;
  core::Size n_decoys_;
  core::Size n_atoms_;
  ObjexxFCL::FArray3D_double coords_;
  ObjexxFCL::FArray2D_double ref_structure_;
  core::Size n_decoys_max_;
  ObjexxFCL::FArray1D_double weights_;

	utility::vector1< core::Real > all_energies_;
	bool store_energies_;

	static bool options_registered_;
};

}
}

#endif
