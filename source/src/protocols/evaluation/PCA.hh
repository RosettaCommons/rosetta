// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_evaluation_PCA_hh
#define INCLUDED_protocols_evaluation_PCA_hh


// Unit Headers
#include <protocols/evaluation/PCA.fwd.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/id/AtomID.fwd.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


// C++ headers

namespace protocols {
namespace evaluation {

class PCA : public utility::pointer::ReferenceCount {
  typedef ObjexxFCL::FArray2D< core::Real > CoordVector;
  typedef utility::vector1< core::id::AtomID > IndexVector;
// this code is stolen from GROMACS. Hence some unfamiliar definitions
  typedef core::Real matrix[3][3];
  typedef core::Real rvec[3];
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~PCA();
  typedef utility::vector1< core::Real > ProjectionVector;

  void read_eigvec_file( std::string fn, core::pose::Pose const& pose, int nvec = -1 );
  void eval( core::pose::Pose const& pose, ProjectionVector& proj );
  void show( std::ostream& );
private:
  void read_structure ( std::istream& data, core::pose::Pose const& pose, IndexVector& ind, CoordVector& x, std::string endtag );
  void fill_coordinates( core::pose::Pose const& pose, IndexVector const& ind, CoordVector & x );
  void reset_x( Size n, CoordVector& x, rvec trans );
  void jacobi(double a[6][6],double d[],double v[6][6],int *nrot);
  void calc_fit_R(int natoms,rvec *xp,rvec const* x,matrix R);
  void rotate_vec(int natoms,rvec *x,matrix R);
  void add_vec( int natoms, rvec *x, rvec transvec );
  void oprod(const rvec a,const rvec b,rvec c);
  CoordVector xref_;
  CoordVector xav_;
  IndexVector ifit_;
  IndexVector ipca_;
  ObjexxFCL::FArray3D< core::Real > eigvec_;
  Size nfit_, npca_, nvec_;
};

}
}

#endif
