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

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/types.hh>
//XRW_B_T1
//#include <core/coarse/Translator.fwd.hh>
//XRW_E_T1
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/GeneralEvent.fwd.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.fwd.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomPointer.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/AtomWithDOFChange.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/DomainMap.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/types.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
//#include <core/optimization/MinimizerMap.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <basic/MetricValue.fwd.hh>
// AUTO-REMOVED #include <basic/OStream.fwd.hh>
#include <utility/stream_util.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector0.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/io/all.fwd.hh>
#include <utility/io/icstream.fwd.hh>
#include <utility/io/irstream.fwd.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/ocstream.fwd.hh>
#include <utility/io/orstream.fwd.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/PausableSignalHub.fwd.hh>
#include <utility/signals/PausableSignalHub.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/BodyPosition.fwd.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/Quaternion.fwd.hh>
#include <numeric/all.fwd.hh>
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzTriple.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray.all.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>


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
	if (atom_num == 0) atom_num = 1;
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
	if (atom_name == "") {
		atom_num = 1;
	}	else {
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
	for(core::Size atom_i = 1; atom_i <= rsd.natoms(); ++atom_i) {
		core::conformation::Residue::AtomIndices const & nbrs = rsd.nbrs(atom_i);
		for(core::conformation::Residue::AtomIndices::const_iterator j = nbrs.begin(), end_j = nbrs.end(); j != end_j; ++j) {
			core::Size atom_j = *j;
			if(atom_j <= atom_i) continue; // so we draw each bond just once, not twice
			print_node(out, rsd.seqpos(), atom_i, conf, "P");
			print_node(out, rsd.seqpos(), atom_j, conf);
		}
	}
	// inter-residue connections
	// there *has* to be a better way of getting next/prev residue...
	out << "@vectorlist {} color= gray width= 1 master= {inter-res}\n";
	core::chemical::ResidueType const & res_type = rsd.type();
	if (rsd.seqpos() > 1 && rsd.is_bonded( conf.residue(rsd.seqpos()-1) )) {
		print_interres_bond(out, rsd, conf.residue(rsd.seqpos()-1), conf);
	}
	if ((core::Size)rsd.seqpos() < conf.size() && rsd.is_bonded( conf.residue(rsd.seqpos()+1) )) {
		print_interres_bond(out, rsd, conf.residue(rsd.seqpos()+1), conf);
	}
	for(core::Size i = 1; i <= res_type.n_residue_connections(); ++i) {
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
	for(core::Size i = 1; i <= conf.size(); ++i) {
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
	for( ; i != i_end; ++i) {
		//std::cout << i->start() << "," << i->start_atom() << " --> " << i->stop() << "," << i->stop_atom() << std::endl;
			print_node(out, i->start(), i->start_atom(), conf, "P");
		if ( i->is_jump() ) print_node(out, i->stop(), i->stop_atom(), conf, "width6");
		else                print_node(out, i->stop(), i->stop_atom(), conf);
	}

	out << "@arrowlist {res-by-res} color= lime radius= 0.6\n";
	i = fold_tree.begin(), i_end = fold_tree.end();
	for( ; i != i_end; ++i) {
		if ( i->is_jump() ) {
			print_node(out, i->start(), i->start_atom(), conf, "P");
			print_node(out, i->stop(), i->stop_atom(), conf, "width4");
		} else {
			print_node(out, i->start(), 0, conf, "P");
			// start res num may be greater or less than stop res num:
			int dir = (i->start() < i->stop() ? 1 : -1);
			for(int j = i->start()+dir, j_end = i->stop()+dir; j != j_end; j+=dir ) {
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

	if (katom.parent() != NULL) {
		core::id::AtomID const & p_atom_id = katom.parent()->atom_id();
		int p_residue_num = p_atom_id.rsd();
		int p_atom_num = p_atom_id.atomno();
		print_node(out, p_residue_num, p_atom_num, conf, "P");
	}

	core::id::AtomID const & atom_id = katom.atom_id();
	int residue_num = atom_id.rsd();
	int atom_num = atom_id.atomno();
	if (katom.is_jump()) print_node(out, residue_num, atom_num, conf, "width4");
	else                 print_node(out, residue_num, atom_num, conf);

	// Recursively visit child atoms
	core::kinematics::tree::Atom::Atoms_ConstIterator i = katom.atoms_begin(), i_end = katom.atoms_end();
	for( ; i != i_end; ++i) visit_atomtree_node(out, **i, conf);
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
	std::string const filename,
	core::pose::Pose const & pose
)
{
	std::ofstream out (filename.c_str());
	if (!out.good()) {
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

