// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/motifs/Motif.cc
/// @brief Implementation of interaction motifs

#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>

#include <core/chemical/ResidueType.hh>

#include <core/kinematics/Jump.hh>

#include <basic/Tracer.hh>


#include <protocols/motifs/SingleMotif.hh>

#include <sstream>

#include <utility/vector1.hh>


namespace protocols {
namespace motifs {

static basic::Tracer mt( "protocols.motifs.Motif", basic::t_info );

SingleMotif::SingleMotif(
	std::string const & resname1,
	std::string const & res1_atom1,
	std::string const & res1_atom2,
	std::string const & res1_atom3,
	std::string const & resname2,
	std::string const & res2_atom1,
	std::string const & res2_atom2,
	std::string const & res2_atom3,
	core::kinematics::Jump const & orientation
) : Motif(
	resname1,
	res1_atom1,
	res1_atom2,
	res1_atom3,
	resname2,
	res2_atom1,
	res2_atom2,
	res2_atom3,
	orientation
	)
{}

SingleMotif::SingleMotif(
	core::pose::Pose const & pose,
	Size const residue_position_1,
	char const chain1,
	std::string const & res1_atom1_name,
	std::string const & res1_atom2_name,
	std::string const & res1_atom3_name,
	Size const residue_position_2,
	char const chain2,
	std::string const & res2_atom1_name,
	std::string const & res2_atom2_name,
	std::string const & res2_atom3_name
) : Motif(
	pose,
	residue_position_1,
	chain1,
	res1_atom1_name,
	res1_atom2_name,
	res1_atom3_name,
	residue_position_2,
	chain2,
	res2_atom1_name,
	res2_atom2_name,
	res2_atom3_name
	)
{}

SingleMotif::SingleMotif(
	core::pose::Pose const & pose,
	Size const residue_position_1,
	std::string const & res1_atom1_name,
	std::string const & res1_atom2_name,
	std::string const & res1_atom3_name,
	Size const residue_position_2,
	std::string const & res2_atom1_name,
	std::string const & res2_atom2_name,
	std::string const & res2_atom3_name
) : Motif(
	pose,
	residue_position_1,
	res1_atom1_name,
	res1_atom2_name,
	res1_atom3_name,
	residue_position_2,
	res2_atom1_name,
	res2_atom2_name,
	res2_atom3_name
	)
{}

SingleMotif::SingleMotif(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2
) : Motif(
	res1,
	res2
	)
{}

SingleMotif::SingleMotif( SingleMotif const & ) = default;

// SingleMotif constructor for ligands
SingleMotif::SingleMotif(
	std::string const & resname1,
	std::string const & res1_atom1,
	std::string const & res1_atom2,
	std::string const & res1_atom3,
	std::string const & res2_atom1,
	std::string const & res2_atom2,
	std::string const & res2_atom3,
	core::kinematics::Jump const & orientation
) : Motif(
	resname1,
	res1_atom1,
	res1_atom2,
	res1_atom3,
	res2_atom1,
	res2_atom2,
	res2_atom3,
	orientation
	)
{}

std::ostream & operator <<(
	std::ostream & os,
	const SingleMotif & motif
)
{
	os << motif.print();
	return os;
}

}
}
