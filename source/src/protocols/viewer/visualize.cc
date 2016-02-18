// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/visualize.cc
/// @brief  3-D visualizations of FoldTree and AtomTree in kinemage format
/// @author Ian W. Davis

// Unit headers
#include <protocols/viewer/visualize.hh>

// Package headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>

// Rosetta headers
#include <core/types.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// ObjexxFCL headers

// C++ Headers
#include <iostream>
#include <fstream>

#include <core/chemical/AtomType.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace viewer {


//////////////////////////////////////////////////////////////////////////////
// helper functions private to this file


void print_node(
	std::ostream & out,
	int residue_num,
	int atom_num,
	core::conformation::Conformation const & conf,
	std::string extras = "" //< P for points, color, width/radius, etc.
)
{
	// atom_num is often 0 in fold tree, means no specific atom.
	// might as well use the first one:
	if ( atom_num == 0 ) atom_num = 1;
	core::conformation::Residue const & res = conf.residue(residue_num);
	core::chemical::ResidueType const & res_type = conf.residue_type(residue_num);
	core::conformation::Atom const & atom = res.atom(atom_num);
	core::chemical::AtomType const & atom_type = res.atom_type(atom_num);
	// This info appears when you click on the point
	out << "{" << res_type.name3() << " " << res.seqpos()
		<< " " << res.atom_name(atom_num) << " (" << atom_type.name() << ")"
		<< "}";
	// Color, width, etc. followed by coordinates
	out << extras;
	out << " " << atom.xyz().x() << " " << atom.xyz().y() << " " << atom.xyz().z() << "\n";
}
void print_node(
	std::ostream & out,
	int residue_num,
	std::string atom_name,
	core::conformation::Conformation const & conf,
	std::string extras = "" //< P for points, color, width/radius, etc.
)
{
	// atom_num is often 0 in fold tree, means no specific atom.
	// might as well use the first one:
	core::conformation::Residue const & res = conf.residue(residue_num);

	int atom_num;
	if ( atom_name == "" ) {
		atom_num = 1;
	} else {
		atom_num = res.atom_index( atom_name );
	}
	print_node( out, residue_num, atom_num, conf, extras );
}


void print_interres_bond(
	std::ostream & out,
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::conformation::Conformation const & conf
)
{
	print_node(out, rsd1.seqpos(), rsd1.connect_atom(rsd2), conf, "P");
	print_node(out, rsd2.seqpos(), rsd2.connect_atom(rsd1), conf);
}


void dump_residue_kinemage(
	std::ostream & out,
	core::conformation::Residue const & rsd,
	core::conformation::Conformation const & conf
)
{
	// intra-residue connections
	// do residues in different (~random) colors to help distinguish them
	int const num_colors = 6;
	std::string colors[num_colors] = {"pinktint", "peachtint", "yellowtint", "greentint", "bluetint", "lilactint"};
	std::string color = colors[ rsd.seqpos() % num_colors ];
	out << "@vectorlist {} color= " << color << " width= 1 master= {intra-res}\n";
	for ( core::Size atom_i = 1; atom_i <= rsd.natoms(); ++atom_i ) {
		core::conformation::Residue::AtomIndices const & nbrs = rsd.nbrs(atom_i);
		for ( core::conformation::Residue::AtomIndices::const_iterator j = nbrs.begin(), end_j = nbrs.end(); j != end_j; ++j ) {
			core::Size atom_j = *j;
			if ( atom_j <= atom_i ) continue; // so we draw each bond just once, not twice
			print_node(out, rsd.seqpos(), atom_i, conf, "P");
			print_node(out, rsd.seqpos(), atom_j, conf);
		}
	}
	// inter-residue connections
	// there *has* to be a better way of getting next/prev residue...
	out << "@vectorlist {} color= gray width= 1 master= {inter-res}\n";
	core::chemical::ResidueType const & res_type = rsd.type();
	if ( rsd.seqpos() > 1 && rsd.is_bonded( conf.residue(rsd.seqpos()-1) ) ) {
		print_interres_bond(out, rsd, conf.residue(rsd.seqpos()-1), conf);
	}
	if ( (core::Size)rsd.seqpos() < conf.size() && rsd.is_bonded( conf.residue(rsd.seqpos()+1) ) ) {
		print_interres_bond(out, rsd, conf.residue(rsd.seqpos()+1), conf);
	}
	for ( core::Size i = 1; i <= res_type.n_possible_residue_connections(); ++i ) {
		print_interres_bond(out, rsd, conf.residue( rsd.residue_connection_partner(i) ), conf);
	}
}


/// @brief A better way to visualize the structure is use Prekin
/// or KiNG's File | Import | Molecules command.
void dump_structure_kinemage(
	std::ostream & out,
	core::conformation::Conformation const & conf
)
{
	out << "@subgroup {by residue} dominant\n";
	for ( core::Size i = 1; i <= conf.size(); ++i ) {
		dump_residue_kinemage(out, conf.residue(i), conf);
	}
}


void dump_foldtree_kinemage(
	std::ostream & out,
	core::kinematics::FoldTree const & fold_tree,
	core::conformation::Conformation const & conf
)
{
	out << "@arrowlist {true} color= gold width=3 radius= 0.6 off\n";
	core::kinematics::FoldTree::const_iterator i = fold_tree.begin(), i_end = fold_tree.end();
	for ( ; i != i_end; ++i ) {
		//std::cout << i->start() << "," << i->start_atom() << " --> " << i->stop() << "," << i->stop_atom() << std::endl;
		print_node(out, i->start(), i->start_atom(), conf, "P");
		if ( i->is_jump() ) print_node(out, i->stop(), i->stop_atom(), conf, "width6");
		else                print_node(out, i->stop(), i->stop_atom(), conf);
	}

	out << "@arrowlist {res-by-res} color= lime radius= 0.6\n";
	i = fold_tree.begin(), i_end = fold_tree.end();
	for ( ; i != i_end; ++i ) {
		if ( i->is_jump() ) {
			print_node(out, i->start(), i->start_atom(), conf, "P");
			print_node(out, i->stop(), i->stop_atom(), conf, "width4");
		} else {
			print_node(out, i->start(), 0, conf, "P");
			// start res num may be greater or less than stop res num:
			int dir = (i->start() < i->stop() ? 1 : -1);
			for ( int j = i->start()+dir, j_end = i->stop()+dir; j != j_end; j+=dir ) {
				print_node(out, j, 0, conf);
			}
		}
	}
}


void visit_atomtree_node(
	std::ostream & out,
	core::kinematics::tree::Atom const & katom,
	core::conformation::Conformation const & conf
)
{
	// Easier to just do point-line all the time than to try and see if
	// previous line was drawn to our parent (it rarely will be).

	if ( katom.parent().get() != NULL ) {
		core::id::AtomID const & p_atom_id = katom.parent()->atom_id();
		int p_residue_num = p_atom_id.rsd();
		int p_atom_num = p_atom_id.atomno();
		print_node(out, p_residue_num, p_atom_num, conf, "P");
	}

	core::id::AtomID const & atom_id = katom.atom_id();
	int residue_num = atom_id.rsd();
	int atom_num = atom_id.atomno();
	if ( katom.is_jump() ) print_node(out, residue_num, atom_num, conf, "width4");
	else                 print_node(out, residue_num, atom_num, conf);

	// Recursively visit child atoms
	core::kinematics::tree::Atom::Atoms_ConstIterator i = katom.atoms_begin(), i_end = katom.atoms_end();
	for ( ; i != i_end; ++i ) visit_atomtree_node(out, **i, conf);
}


void dump_atomtree_kinemage(
	std::ostream & out,
	core::kinematics::AtomTree const & atom_tree,
	core::conformation::Conformation const & conf
)
{
	out << "@arrowlist {true} color= orange\n";
	core::kinematics::tree::Atom const & root = *( atom_tree.root() );
	visit_atomtree_node(out, root, conf);
}


//////////////////////////////////////////////////////////////////////////////
// public functions

void
dump_pose_kinemage(
	std::string const & filename,
	core::pose::Pose const & pose
)
{
	std::ofstream out (filename.c_str());
	if ( !out.good() ) {
		basic::Error() << "Can't open kinemage file " << filename << std::endl;
		return;
	}
	out << "@text\n";
	out << "View this file with KiNG or Mage from http://kinemage.biochem.duke.edu\n";
	out << "@kinemage 1\n";
	out << "@onewidth\n";
	out << "@group {structure}\n";
	dump_structure_kinemage(out, pose.conformation());
	out << "@group {fold tree} animate\n";
	dump_foldtree_kinemage(out, pose.fold_tree(), pose.conformation());
	out << "@group {atom tree} animate\n";
	dump_atomtree_kinemage(out, pose.atom_tree(), pose.conformation());
	out.close();
}


} // namespace viewer
} // namespace protocols

