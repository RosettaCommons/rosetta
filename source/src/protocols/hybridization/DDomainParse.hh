// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief A reimplementation of DDomain domain parsing algorithm, Protein Sci.16,947--955 (2007)
/// @author Frank DiMaio

#ifndef INCLUDED_protocols_hybridization_DDomainParse_hh
#define INCLUDED_protocols_hybridization_DDomainParse_hh

#include <protocols/loops/Loops.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Edge.hh>

#include <core/types.hh>

#include <numeric/xyzVector.hh>

#include <list>

namespace protocols {
//namespace comparative_modeling {
namespace hybridization {


class DDomainParse {
public:
	DDomainParse( ) {
		pcut_ = 0.81;
		hcut_ = 0.18;
		length_ = 38;
	}

	DDomainParse( core::Real pcut_in, core::Real hcut_in, core::Size length_in) {
		pcut_ = pcut_in;
		hcut_ = hcut_in;
		length_ = length_in;
	}

	// domain splitting
	utility::vector1< loops::Loops > split( core::pose::Pose const &templ );

private:
	// functions
	void pulldomain( int isize0, int ist, int ilast, utility::vector1< core::Real > &resect, int &idip);
	void findpos( int mdom, int ndom, int id, utility::vector1<int> &ipdom, int &ist, int &ilast);
	void findpos2( int ist, int ilast, int imid, int i, int &ista, int &ilasta);
	void small_big( int const& mdom, int const& ipdd, utility::vector1<int> &ipdom);
	void distance( int ist, int ilast, utility::vector1< numeric::xyzVector<core::Real> > const &xs, utility::vector1< utility::vector1< core::Real > > &dij);
	void ddomain_pot( core::Real pw, int ist, int ilast,
		utility::vector1< utility::vector1< core::Real > > const &dij, utility::vector1< core::Real > &resect);

	// parameters
	core::Real pcut_,hcut_;
	core::Size length_;

	core::Size mseq_;
	core::Size nseq_;
};


} // hybridize
//} // comparative_modeling
} // protocols

#endif

