// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/GrowLigand.cc
///
/// @brief
/// @Gordon Lemmon

#include <protocols/ligand_docking/GrowLigand.hh>
#include <protocols/ligand_docking/GrowLigandCreator.hh>

#include <protocols/ligand_docking/LigandDesign.hh> // For helper functions.  Refactor this
#include <core/pose/util.hh>


#include <protocols/moves/Mover.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

// option key includes

// AUTO-REMOVED #include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <utility/tag/Tag.hh>

#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
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
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/types.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverCreator.hh>
#include <protocols/moves/MoverStatus.hh>
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
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/options/StringOption.hh>
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
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <utility/tag/Tag.fwd.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <numeric/random/random.fwd.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <numeric/random/uniform.fwd.hh>
#include <numeric/random/uniform.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end
#define foreach BOOST_FOREACH

namespace protocols {
namespace ligand_docking {

static basic::Tracer grow_ligand_tracer("protocols.ligand_docking.GrowLigand", basic::t_debug);

std::string
GrowLigandCreator::keyname() const
{
	return GrowLigandCreator::mover_name();
}

protocols::moves::MoverOP
GrowLigandCreator::create_mover() const {
	return new GrowLigand;
}

std::string
GrowLigandCreator::mover_name()
{
	return "GrowLigand";
}

GrowLigand::GrowLigand():
		Mover("GrowLigand"),
		chain_("")
{
	set_fragments();
}

GrowLigand::GrowLigand(std::string chain):
		Mover("GrowLigand"),
		chain_(chain)
{
	set_fragments();
}

GrowLigand::GrowLigand(GrowLigand const & that):
	    //utility::pointer::ReferenceCount(),
		protocols::moves::Mover( that ),
		chain_(that.chain_),
		fragments_(that.fragments_)
{}

GrowLigand::~GrowLigand() {}

void
GrowLigand::set_fragments(){
	core::chemical::ResidueSelector rs;
	rs.set_property("FRAGMENT");
	core::chemical::ChemicalManager *cm= core::chemical::ChemicalManager::get_instance();
	core::chemical::ResidueTypeSetCAP rsd_set= cm->residue_type_set( core::chemical::FA_STANDARD );
	core::chemical::AtomTypeSetCAP atom_set= cm->atom_type_set( core::chemical::FA_STANDARD );
	core::Size atom_set_size= atom_set->n_atomtypes();
	fragments_.assign(atom_set_size, utility::vector1< std::pair<core::conformation::ResidueCOP, core::Size> >() );
	core::chemical::ResidueTypeCAPs fragment_types= rs.select( *rsd_set );
	grow_ligand_tracer<< fragment_types.size()<< " fragment_types"<< std::endl;

	foreach(core::chemical::ResidueTypeCAP fragment_type, fragment_types){
		core::conformation::Residue* temp= new core::conformation::Residue(*fragment_type, true);
		for(core::Size i=1; i<= temp->n_residue_connections(); ++i){
			core::chemical::ResidueConnection const & res_conn= temp->residue_connection(i);
			int atom_index_number= res_conn.atomno();
			int atom_type_index= temp->atom_type_index(atom_index_number);
			std::pair<core::conformation::ResidueCOP, core::Size> pair(temp, i);
			fragments_[atom_type_index].push_back(pair);
		}

		grow_ligand_tracer<< "frag_name: "<< temp->name()<< std::endl;
	}
}

protocols::moves::MoverOP GrowLigand::clone() const {
	return new GrowLigand( *this );
}

protocols::moves::MoverOP GrowLigand::fresh_instance() const {
	return new GrowLigand;
}

std::string GrowLigand::get_name() const{
	return "GrowLigand";
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
GrowLigand::parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & /*datamap*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "GrowLigand" ) {
		grow_ligand_tracer << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	if ( tag->hasOption("chain") ) {
		chain_ = tag->getOption<std::string>("chain");
	}else{
		utility_exit_with_message("HeavyAtom filter needs a 'chain' option");
	}
}

void
GrowLigand::apply( core::pose::Pose & pose )
{
	assert(!fragments_.empty());
	assert(chain_.size() == 1);

	utility::vector1<core::Size> unconnected_residues;
	{
		core::Size const chain_id= core::pose::get_chain_id_from_chain(chain_, pose);;
		core::Size const start = pose.conformation().chain_begin(chain_id);
		core::Size const end = pose.conformation().chain_end(chain_id);
		unconnected_residues=find_unconnected_residues(pose, start, end);
	}

	numeric::random::random_permutation(unconnected_residues, numeric::random::RG);// shuffle vector

	foreach(core::Size grow_from, unconnected_residues){
		core::Size grow_from_connection= random_connection(&pose.residue(grow_from));
		core::Size atom_type_index= pose.residue(grow_from).residue_connection(grow_from_connection).atom_type_index();
		utility::vector1< std::pair<core::conformation::ResidueCOP, core::Size > >
				residue_connection_pairs= fragments_[atom_type_index];
		if( residue_connection_pairs.size() > 0 ){
			std::pair<core::conformation::ResidueCOP, core::Size > const growth= numeric::random::random_element(residue_connection_pairs);
			bool build_ideal_geometry= true;
			pose.append_residue_by_bond(*(growth.first), build_ideal_geometry, growth.second, grow_from, grow_from_connection);
			return;
		}
	}

}

void GrowLigand::fragments_to_string() const{
	for(core::Size i=1; i <= fragments_.size(); ++i){
		utility::vector1< std::pair<core::conformation::ResidueCOP, core::Size> >::const_iterator  begin= fragments_[i].begin();
		for(; begin != fragments_[i].end(); ++begin){
			core::conformation::Residue const & res= *(begin->first);
			core::Size connect_id= begin->second;
			std::string name= res.name();
			grow_ligand_tracer<< "atom_type, res_name, connection"<< i << " "<< name << " "<< connect_id<< std::endl;
		}
	}
}

void GrowLigand::add_scores_to_job(
	core::pose::Pose & /*pose*/
)
{
}

} // namespace ligand_docking
} // namespace protocols
