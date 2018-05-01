// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraph.hh
/// @brief The BuriedUnsatPenaltyGraph class and its related Node and Edge classes.
/// @author Vikram K. Mulligan (vmullig@uw.edu).


#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraph.hh>
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphOptions.hh>

#include <basic/Tracer.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/types.hh>
#include <core/select/util/burial_utilities.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/PackerTask_.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/MirrorSymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Boost Headers
#include <utility/graph/unordered_object_pool.hpp>
#include <utility/pointer/memory.hh>
#include <boost/pool/pool.hpp>

static basic::Tracer TR( "core.pack.guidance_scoreterms.buried_unsat_penalty.graph.BuriedUnsatPenaltyGraph" );

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace buried_unsat_penalty {
namespace graph {

/*******************************************************************************
class BuriedUnsatPenaltyGraphHbond
*******************************************************************************/

/// @brief Default constructor.
BuriedUnsatPenaltyGraphHbond::BuriedUnsatPenaltyGraphHbond() :
	first_node_is_the_acceptor_(false),
	donor_group_index_(0),
	acceptor_group_index_(0),
	energy_(0),
	lower_numbered_node_symmetry_copy_index_(1),
	higher_numbered_node_symmetry_copy_index_(1)
{}

/// @brief Initialization constructor.
BuriedUnsatPenaltyGraphHbond::BuriedUnsatPenaltyGraphHbond(
	bool const first_node_is_the_acceptor,
	core::Size const donor_group_index,
	core::Size const acceptor_group_index,
	core::Real const energy,
	core::Size const lower_numbered_node_symmetry_copy_index,
	core::Size const higher_numbered_node_symmetry_copy_index
) :
	first_node_is_the_acceptor_(first_node_is_the_acceptor),
	donor_group_index_( donor_group_index ),
	acceptor_group_index_( acceptor_group_index ),
	energy_(energy),
	lower_numbered_node_symmetry_copy_index_(lower_numbered_node_symmetry_copy_index),
	higher_numbered_node_symmetry_copy_index_(higher_numbered_node_symmetry_copy_index)
{
	runtime_assert_string_msg( donor_group_index_ > 0, "Error in BuriedUnsatPenaltyGraphHbond constructor: the donor group index must be positive." );
	runtime_assert_string_msg( acceptor_group_index_ > 0, "Error in BuriedUnsatPenaltyGraphHbond constructor: the acceptor group index must be positive." );
	runtime_assert_string_msg( lower_numbered_node_symmetry_copy_index_ > 0, "Error in BuriedUnsatPenaltyGraphHbond constructor: the lower node symmetry copy index must be positive." );
	runtime_assert_string_msg( higher_numbered_node_symmetry_copy_index_ > 0, "Error in BuriedUnsatPenaltyGraphHbond constructor: the higher node symmetry copy index must be positive." );
}

/// @brief Copy constructor.
BuriedUnsatPenaltyGraphHbond::BuriedUnsatPenaltyGraphHbond( BuriedUnsatPenaltyGraphHbond const &src ) :
	first_node_is_the_acceptor_(src.first_node_is_the_acceptor_),
	donor_group_index_(src.donor_group_index_),
	acceptor_group_index_(src.acceptor_group_index_),
	energy_(src.energy_),
	lower_numbered_node_symmetry_copy_index_(src.lower_numbered_node_symmetry_copy_index_),
	higher_numbered_node_symmetry_copy_index_(src.higher_numbered_node_symmetry_copy_index_)
{}

/// @brief Destructor.
BuriedUnsatPenaltyGraphHbond::~BuriedUnsatPenaltyGraphHbond() = default;

/// @brief Assignment operator.
BuriedUnsatPenaltyGraphHbond &
BuriedUnsatPenaltyGraphHbond::operator= ( BuriedUnsatPenaltyGraphHbond const & src ) {
	first_node_is_the_acceptor_ = src.first_node_is_the_acceptor_;
	donor_group_index_ = src.donor_group_index_;
	acceptor_group_index_ = src.acceptor_group_index_;
	energy_ = src.energy_;
	lower_numbered_node_symmetry_copy_index_ = src.lower_numbered_node_symmetry_copy_index_;
	higher_numbered_node_symmetry_copy_index_ = src.higher_numbered_node_symmetry_copy_index_;
	return *this;
}

/// @brief Get the symmetry copy index of the node that's the acceptor.
/// @details Returns 1 in the asymmetric case.
core::Size
BuriedUnsatPenaltyGraphHbond::acceptor_symmetry_copy_index() const {
	if ( first_node_is_the_acceptor_ ) return lower_numbered_node_symmetry_copy_index_;
	return higher_numbered_node_symmetry_copy_index_;
}

/// @brief Get the symmetry copy index of the node that's the donor.
/// @details Returns 1 in the asymmetric case.
core::Size
BuriedUnsatPenaltyGraphHbond::donor_symmetry_copy_index() const {
	if ( first_node_is_the_acceptor_ ) return higher_numbered_node_symmetry_copy_index_;
	return lower_numbered_node_symmetry_copy_index_;
}

/*******************************************************************************
class BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup
*******************************************************************************/

/// @brief Options constructor.
/// @details Note that protons are initialized to an empty list, and must be added later.
BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup::BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup(
	bool const is_acceptor,
	bool const is_counted,
	core::Size const heavyatom_index,
	core::conformation::Residue const &residue
) :
	is_acceptor_(is_acceptor),
	is_counted_(is_counted),
	heavyatom_index_(heavyatom_index),
	proton_indices_(),
	max_accepted_hbond_count_( determine_max_accepted_hbond_count( is_acceptor, heavyatom_index, residue ) )
{}

/// @brief Destructor.
BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup::~BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup() {}

/// @brief Assignment operator.
BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup &
BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup::operator=(
	BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup const &src
) {
	is_acceptor_ = src.is_acceptor_;
	is_counted_ = src.is_counted_;
	heavyatom_index_ = src.heavyatom_index_;
	proton_indices_ = src.proton_indices_;
	max_accepted_hbond_count_ = src.max_accepted_hbond_count_;
	return (*this);
}

/// @brief Add a proton to the list of protons in this group.
void
BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup::add_proton_index( core::Size const index_in ) {
	runtime_assert( !proton_indices_.has_value(index_in) );
	proton_indices_.push_back(index_in);
}

/// @brief Determine the maximum number of hydrogen bonds that a group can accept, based on the identity of the heavyatom.
/// @details Oxygens can accept 2; nitrogens can accept 1.  For now, the rule is that simple.  At some point, we might look up the information from a database
/// lookup table based on the details of the oxygen type or nitrogen type or whatnot.
/// @note Returns 0 if "is_acceptor" is false.  This is a static function.
core::Size
BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup::determine_max_accepted_hbond_count(
	bool const is_acceptor,
	core::Size const heavyatom_index,
	core::conformation::Residue const & residue
) {
	if ( !is_acceptor ) return 0;
	return ( residue.type().atom( heavyatom_index ).element() == core::chemical::element::O ? 2 : 1 );
}

/*******************************************************************************
class BuriedUnsatPenaltyNodeData
*******************************************************************************/

/// @brief Default constructor.
BuriedUnsatPenaltyNodeData::BuriedUnsatPenaltyNodeData():
	utility::pointer::ReferenceCount(),
	residue_position_(0),
	rotamer_index_(0),
	donor_acceptor_groups_(),
	group_heavyatom_index_to_group_index_(),
	donor_acceptor_groups_intrares_hbonds_donated_(),
	donor_acceptor_groups_intrares_hbonds_accepted_(),
	donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated_(),
	donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted_(),
	options_(nullptr),
	hbond_options_(nullptr)
{}

/// @brief Destructor.
BuriedUnsatPenaltyNodeData::~BuriedUnsatPenaltyNodeData() {}

/// @brief Create a copy by owning pointer and return the owning pointer to the copy.
BuriedUnsatPenaltyNodeDataOP
BuriedUnsatPenaltyNodeData::clone() const {
	return utility::pointer::make_shared< BuriedUnsatPenaltyNodeData >( *this );
}

/// @brief Set the residue position.
void
BuriedUnsatPenaltyNodeData::residue_position(
	core::Size const &setting
) {
	residue_position_ = setting;
}


/// @brief Set the residue position.
void
BuriedUnsatPenaltyNodeData::rotamer_index(
	core::Size const &setting
) {
	rotamer_index_ = setting;
}

/// @brief Set the options pointer.
/// @details This copies the input pointer; it doesn't clone the object.
void
BuriedUnsatPenaltyNodeData::options(
	BuriedUnsatPenaltyGraphOptionsCOP options_in
) {
	debug_assert( options_in != nullptr );
	options_ = options_in;
}

/// @brief Set the hbonds options pointer.
/// @details This copies the input pointer; it doesn't clone the object.
void
BuriedUnsatPenaltyNodeData::hbond_options(
	core::scoring::hbonds::HBondOptionsCOP options_in
) {
	debug_assert( options_in != nullptr );
	hbond_options_ = options_in;
}


/*******************************************************************************
class BuriedUnsatPenaltyNode
*******************************************************************************/

/// @brief Node constructor.
BuriedUnsatPenaltyNode::BuriedUnsatPenaltyNode( utility::graph::Graph * owner, platform::Size const node_id ) :
	utility::graph::Node( owner, node_id ),
	stored_data_( utility::pointer::make_shared< BuriedUnsatPenaltyNodeData >() ),
	donated_hbond_count_(),
	accepted_hbond_count_()
{}

/// @brief Destructor.
BuriedUnsatPenaltyNode::~BuriedUnsatPenaltyNode() = default;

//////////////////////
// These functions are public, but are only intended to be called by BuriedUnsatPenaltyGraph class:
//////////////////////

/// @brief Initialize this node, setting its residue_position, rotamer_index,
void
BuriedUnsatPenaltyNode::initialize_node(
	core::Size const residue_position,
	core::Size const rotamer_index,
	core::conformation::ResidueCOP residue,
	core::pose::Pose const &pose,
	BuriedUnsatPenaltyGraphOptionsCOP options,
	core::scoring::hbonds::HBondOptionsCOP hbond_options,
	core::scoring::hbonds::HBondDatabase const &hbond_database,
	bool const is_symmetric
) {
	core::conformation::symmetry::SymmetricConformationCOP const symmconf( is_symmetric ? utility::pointer::static_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( pose.conformation_ptr() ) : nullptr  );
	core::conformation::symmetry::SymmetryInfoCOP const symminfo( is_symmetric ? symmconf->Symmetry_Info() : nullptr );
	core::Size const nres( is_symmetric ? symminfo->num_independent_residues() : pose.total_residue() );

	runtime_assert( residue_position > 0 && residue_position <= nres );
	runtime_assert( rotamer_index > 0 );
	runtime_assert( options != nullptr );

	BuriedUnsatPenaltyNodeDataOP new_stored_data( stored_data_->clone() );
	stored_data_ = new_stored_data; //Const from nonconst.  We can continue to modify the object pointed to using new_stored_data; when it goes out of scope, stored_data_ is the COP to that object.

	new_stored_data->residue_position( residue_position );
	new_stored_data->rotamer_index( rotamer_index );
	new_stored_data->options( options );
	new_stored_data->hbond_options( hbond_options );

	utility::vector1< BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup > &donor_acceptor_groups( new_stored_data->donor_acceptor_groups() );
	donor_acceptor_groups.clear();

	//Loop through atoms in this residue and detect all the donor and acceptor groups:
	for ( core::Size ia(1), iamax(residue->natoms()); ia<=iamax; ++ia ) {
		core::chemical::AtomType const & atomtype( residue->atom_type(ia) );
		if ( atomtype.is_hydrogen() ) continue; // Hydrogens are handled in second pass.
		if ( atomtype.is_acceptor() ) {
			donor_acceptor_groups.push_back( BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup( true, is_buried( ia, *residue, pose ), ia, *residue  ) );
			new_stored_data->group_heavyatom_index_to_group_index()[ia] = donor_acceptor_groups.size();
		}
	}

	//Loop through atoms in this residue again and add polar hydrogens to groups:
	for ( core::Size ia(1), iamax(residue->natoms()); ia<=iamax; ++ia ) {
		core::chemical::AtomType const & atomtype( residue->atom_type(ia) );
		if ( atomtype.is_polar_hydrogen() ) {
			core::Size const parent_atom( residue->icoor(ia).stub_atom1().atomno() );
			debug_assert( parent_atom != 0 );
			core::Size group_index;
			if ( new_stored_data->group_heavyatom_index_to_group_index().find(parent_atom) == new_stored_data->group_heavyatom_index_to_group_index().end() ) {
				donor_acceptor_groups.push_back( BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup( false, is_buried(parent_atom, *residue, pose), parent_atom, *residue ) );
				new_stored_data->group_heavyatom_index_to_group_index()[parent_atom] = donor_acceptor_groups.size();
				group_index = donor_acceptor_groups.size();
			} else {
				group_index = new_stored_data->group_heavyatom_index_to_group_index().at(parent_atom);
			}
			donor_acceptor_groups[ group_index ].add_proton_index( ia );
		}
	}

	//Detect intra-residue hydrogen bonds:
	detect_intra_residue_hydrogen_bonds( *residue, hbond_database, new_stored_data );
	//In the symmetric case, detect hydrogen bonds between this residue and its symmetry mates:
	if ( is_symmetric ) {
		detect_intra_residue_hydrogen_bonds_symmetric( residue_position, residue, pose, symminfo, hbond_database, *hbond_options, new_stored_data );
	}

	donated_hbond_count_.resize( donor_acceptor_groups.size() );
	accepted_hbond_count_.resize( donor_acceptor_groups.size() );
	for ( core::Size i(1), imax(donor_acceptor_groups.size()); i<=imax; ++i ) {
		donated_hbond_count_[i] = stored_data_->donor_acceptor_groups_intrares_hbonds_donated()[i] + stored_data_->donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated()[i];
		accepted_hbond_count_[i] = stored_data_->donor_acceptor_groups_intrares_hbonds_accepted()[i]+ stored_data_->donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted()[i];
	}

	TR.Debug << "Initialized " << (is_symmetric ? "symmetric " : " " ) << "node " << get_node_index() << " (position " << stored_data_->residue_position() << ", rotamer " << stored_data_->rotamer_index() << ") with " << stored_data_->donor_acceptor_groups().size() << " donor/acceptor groups." << std::endl;
}

/// @brief Given another node, set this node to copy that one.
/// @details Note: this does NOT copy the rotamer index, but instead sets it to 1.  The donated_hbond_count_
/// and accepted_hbond_count_ vars are initialized to the intra-residue donated and intra-residue accepted
/// counts, respectively (i.e. inter-residue donated and inter-residue accepted counts have *not* yet been
/// added in).
void
BuriedUnsatPenaltyNode::copy_from(
	utility::graph::Node const * other_node
) {
	BuriedUnsatPenaltyNode const &other( *(static_cast< BuriedUnsatPenaltyNode const * >(other_node)) );
	stored_data_ = other.stored_data_; //Const from const
	clear_hbond_counts();
}

/// @brief For every donor/acceptor group, reset the counts of hydrogen bonds to/from that group.
/// @details The donated_hbond_count_ and accepted_hbond_count_ vars are initialized to the
/// intra-residue donated and intra-residue accepted counts, respectively (i.e. inter-residue
/// donated and inter-residue accepted counts have *not* yet been added in).
void
BuriedUnsatPenaltyNode::clear_hbond_counts() {
	core::Size const ngroups( stored_data_->donor_acceptor_groups().size() );
	donated_hbond_count_.resize( ngroups );
	accepted_hbond_count_.resize( ngroups );
	for ( core::Size i(1); i<=ngroups; ++i ) {
		//Note that, in the following two lines, donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated_[i] and donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted_[i] will each be 0 in the asymmetric case.
		donated_hbond_count_[i] = stored_data_->donor_acceptor_groups_intrares_hbonds_donated()[i] + stored_data_->donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated()[i];
		accepted_hbond_count_[i] = stored_data_->donor_acceptor_groups_intrares_hbonds_accepted()[i] + stored_data_->donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted()[i];
	}
}

/// @brief Add one to the number of accepted hydrogen bonds for the Nth donor/acceptor group in this node.
void
BuriedUnsatPenaltyNode::increment_accepted_hbond_count_for_group(
	core::Size const group_index
) {
	debug_assert(group_index > 0 && group_index <= accepted_hbond_count_.size());
	++(accepted_hbond_count_[group_index]);
}

/// @brief Add one to the number of donated hydrogen bonds for the Nth donor/acceptor group in this node.
void
BuriedUnsatPenaltyNode::increment_donated_hbond_count_for_group(
	core::Size const group_index
) {
	debug_assert(group_index > 0 && group_index <= donated_hbond_count_.size());
	++(donated_hbond_count_[group_index]);
}

/// @brief After hydrogen bonds from edges have been counted, report the number of donor/acceptor groups in each of several categories.
/// @details The integer instances passed in are incremented appropriately by this function.  (So if there is an incoming count of 5 unsaturated acceptors, and
/// the current node has 3 more, the unsat_acceptor_count will end up being 8).
/// @param[out] unsat_acceptor_count The number of acceptor (and not donor) groups that are unsatisfied.
/// @param[out] unsat_donor_count The number of donor (and not acceptor) groups that are unsatisfied.
/// @param[out] unsat_acceptor_and_donor_count The number of groups that are both donors and acceptors (e.g. hydroxyls) that are unsatisfied (i.e. lack either a donated hbond or an accepted hbond).
/// @param[out] oversat_acceptor_count The number of acceptor (and not donor) groups that are oversatisfied.
/// @param[out] oversat_donor_count The number of donor (and not acceptor) groups that are oversatisfied.  (This generally doesn't happen).
/// @param[out] oversat_acceptor_and_donor_count The number of groups that are both donors and acceptors (e.g. hydroxyls) that are unsatisfied (i.e. have either too many donated hbonds or too many accepted hbonds).
void
BuriedUnsatPenaltyNode::increment_counts(
	core::Size &unsat_acceptor_count,
	core::Size &unsat_donor_count,
	core::Size &unsat_acceptor_and_donor_count,
	core::Size &oversat_acceptor_count,
	core::Size &oversat_donor_count,
	core::Size &oversat_acceptor_and_donor_count
) const {
	utility::vector1< BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup > const & donor_acceptor_groups( stored_data_->donor_acceptor_groups() );

	for ( core::Size i(1), imax(donor_acceptor_groups.size()); i<=imax; ++i ) { //Loop through all donor/acceptor groups
		BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup const & curgroup( donor_acceptor_groups[i] );
		core::Size const donated_hbond_count( donated_hbond_count_[i] );
		core::Size const accepted_hbond_count( accepted_hbond_count_[i] );
		if ( !curgroup.is_counted() ) continue;
		if ( curgroup.is_acceptor() ) {
			if ( curgroup.is_donor() ) { //Acceptor AND donor
				if ( donated_hbond_count == 0 && accepted_hbond_count == 0 ) ++unsat_acceptor_and_donor_count;
				if ( donated_hbond_count > curgroup.max_donated_hbond_count() || accepted_hbond_count > curgroup.max_accepted_hbond_count() ) ++oversat_acceptor_and_donor_count;
			} else { //Acceptor only, but NOT donor
				if ( accepted_hbond_count == 0 ) ++unsat_acceptor_count;
				if ( accepted_hbond_count > curgroup.max_accepted_hbond_count() ) ++oversat_acceptor_count;
			}
		} else {
			if ( curgroup.is_donor() ) { //Donor only, but not acceptor
				if ( donated_hbond_count == 0 ) ++unsat_donor_count;
				if ( donated_hbond_count > curgroup.max_donated_hbond_count() ) ++oversat_donor_count;
			} else {
				utility_exit_with_message( "Error in BuriedUnsatPenaltyNode::increment_counts(): A hydrogen bond donor/acceptor group was encountered that was neither donor nor acceptor.  This should be impossilbe." );
			}
		}
	} //Loop through all donor/acceptor groups
}

/// @brief After hydrogen bonds from edges have been counted, report the number of donor/acceptor groups in each of several categories.
/// @details The integer instances passed in are DECREMENTED appropriately by this function.  (So if there is an incoming count of 5 unsaturated acceptors, and
/// the current node accounts for 3 unsaturated acceptors, the unsat_acceptor_count will drop to 2).
/// @param[out] unsat_acceptor_count The number of acceptor (and not donor) groups that are unsatisfied.
/// @param[out] unsat_donor_count The number of donor (and not acceptor) groups that are unsatisfied.
/// @param[out] unsat_acceptor_and_donor_count The number of groups that are both donors and acceptors (e.g. hydroxyls) that are unsatisfied (i.e. lack either a donated hbond or an accepted hbond).
/// @param[out] oversat_acceptor_count The number of acceptor (and not donor) groups that are oversatisfied.
/// @param[out] oversat_donor_count The number of donor (and not acceptor) groups that are oversatisfied.  (This generally doesn't happen).
/// @param[out] oversat_acceptor_and_donor_count The number of groups that are both donors and acceptors (e.g. hydroxyls) that are unsatisfied (i.e. have either too many donated hbonds or too many accepted hbonds).
void
BuriedUnsatPenaltyNode::decrement_counts(
	core::Size &unsat_acceptor_count,
	core::Size &unsat_donor_count,
	core::Size &unsat_acceptor_and_donor_count,
	core::Size &oversat_acceptor_count,
	core::Size &oversat_donor_count,
	core::Size &oversat_acceptor_and_donor_count
) const {
	utility::vector1< BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup > const & donor_acceptor_groups( stored_data_->donor_acceptor_groups() );

	for ( core::Size i(1), imax(donor_acceptor_groups.size()); i<=imax; ++i ) { //Loop through all donor/acceptor groups
		BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup const & curgroup( donor_acceptor_groups[i] );
		core::Size const donated_hbond_count( donated_hbond_count_[i] );
		core::Size const accepted_hbond_count( accepted_hbond_count_[i] );
		if ( !curgroup.is_counted() ) continue;
		if ( curgroup.is_acceptor() ) {
			if ( curgroup.is_donor() ) { //Acceptor AND donor
				if ( donated_hbond_count == 0 && accepted_hbond_count == 0 ) {
					runtime_assert( unsat_acceptor_and_donor_count > 0 );
					--unsat_acceptor_and_donor_count;
				}
				if ( donated_hbond_count > curgroup.max_donated_hbond_count() || accepted_hbond_count > curgroup.max_accepted_hbond_count() ) {
					runtime_assert( oversat_acceptor_and_donor_count > 0 );
					--oversat_acceptor_and_donor_count;
				}
			} else { //Acceptor only, but NOT donor
				if ( accepted_hbond_count == 0 ) {
					runtime_assert( unsat_acceptor_count > 0 );
					--unsat_acceptor_count;
				}
				if ( accepted_hbond_count > curgroup.max_accepted_hbond_count() ) {
					runtime_assert( oversat_acceptor_count > 0 );
					--oversat_acceptor_count;
				}
			}
		} else {
			if ( curgroup.is_donor() ) { //Donor only, but not acceptor
				if ( donated_hbond_count == 0 ) {
					runtime_assert( unsat_donor_count > 0 );
					--unsat_donor_count;
				}
				if ( donated_hbond_count > curgroup.max_donated_hbond_count() ) {
					runtime_assert( oversat_donor_count > 0 );
					--oversat_donor_count;
				}
			} else {
				utility_exit_with_message( "Error in BuriedUnsatPenaltyNode::increment_counts(): A hydrogen bond donor/acceptor group was encountered that was neither donor nor acceptor.  This should be impossilbe." );
			}
		}
	} //Loop through all donor/acceptor groups
}

//////////////////////
// Private functions
//////////////////////

/// @brief Determine whether an atom is buried by the method of sidechain neighbors.
bool
BuriedUnsatPenaltyNode::is_buried(
	core::Size const atom_index,
	core::conformation::Residue const & residue,
	core::pose::Pose const &pose
) const {
	BuriedUnsatPenaltyGraphOptions const &options( *(stored_data_->options()) );
	return core::select::util::determine_whether_point_is_buried(
		residue.xyz( atom_index ),
		pose,
		options.angle_exponent(),
		options.angle_shift_factor(),
		options.dist_exponent(),
		options.dist_midpoint(),
		options.burial_threshold()
	);
}

/// @brief Given a residue, detect intramolecular hydrogen bonds.
/// @details Populates the donor_acceptor_groups_intrares_hbonds_donated_ and donor_acceptor_groups_intrares_hbonds_accepted_ lists.
/// @note This version is for the ASYMMETRIC and SYMMETRIC cases.  It only detects the intraresidue hbonds.  In the symmetric case,
/// the detect_intra_residue_hydrogen_bonds_symmetric() function should be called immediately after this to detect the interresidue hbonds
/// between this residue and its symmetry mates.
void
BuriedUnsatPenaltyNode::detect_intra_residue_hydrogen_bonds(
	core::conformation::Residue const &residue,
	core::scoring::hbonds::HBondDatabase const &hbond_database,
	BuriedUnsatPenaltyNodeDataOP new_stored_data
) {

	//Initialize:
	utility::vector1< BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup > const & donor_acceptor_groups( new_stored_data->donor_acceptor_groups() );
	utility::vector1< core::Size > & donor_acceptor_groups_intrares_hbonds_donated( new_stored_data->donor_acceptor_groups_intrares_hbonds_donated() );
	utility::vector1< core::Size > & donor_acceptor_groups_intrares_hbonds_accepted( new_stored_data->donor_acceptor_groups_intrares_hbonds_accepted() );
	std::map< core::Size, core::Size > & group_heavyatom_index_to_group_index( new_stored_data->group_heavyatom_index_to_group_index() );

	donor_acceptor_groups_intrares_hbonds_donated.clear();
	donor_acceptor_groups_intrares_hbonds_donated.resize( donor_acceptor_groups.size(), 0 );
	donor_acceptor_groups_intrares_hbonds_accepted.clear();
	donor_acceptor_groups_intrares_hbonds_accepted.resize( donor_acceptor_groups.size(), 0 );

	utility::vector1< core::Size > &donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted(new_stored_data->donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted() );
	utility::vector1< core::Size > &donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated(new_stored_data->donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated() );
	donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted.clear();
	donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated.clear();
	donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted.resize( donor_acceptor_groups.size(), 0 );
	donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated.resize( donor_acceptor_groups.size(), 0 );

	core::scoring::hbonds::HBondSet hbondset( *( new_stored_data->hbond_options() ) );
	core::scoring::hbonds::identify_intra_res_hbonds( hbond_database, residue, 0, false, hbondset, false, false, false, false );

	if ( hbondset.nhbonds() > 0 ) {
		core::Real const hbond_energy_threshold( new_stored_data->options()->hbond_energy_threshold() );
		for ( core::Size i(1), imax(hbondset.nhbonds()); i<=imax; ++i ) {
			core::scoring::hbonds::HBond const &curhbond( hbondset.hbond(i) );
			if ( curhbond.energy() < hbond_energy_threshold ) {
				core::Size const don_hatm_parent( residue.icoor( curhbond.don_hatm() ).stub_atom1().atomno() );
				debug_assert( group_heavyatom_index_to_group_index.find(don_hatm_parent) != group_heavyatom_index_to_group_index.end() );
				debug_assert( group_heavyatom_index_to_group_index.find(curhbond.acc_atm()) != group_heavyatom_index_to_group_index.end() );
				core::Size const donor_group( group_heavyatom_index_to_group_index.at(don_hatm_parent) );
				core::Size const acceptor_group( group_heavyatom_index_to_group_index.at(curhbond.acc_atm()) );
				++donor_acceptor_groups_intrares_hbonds_donated[donor_group];
				++donor_acceptor_groups_intrares_hbonds_accepted[acceptor_group];
			}
		}
	}
}

/// @brief Given a residue, detect intramolecular hydrogen bonds.
/// @details Populates the donor_acceptor_groups_intrares_hbonds_donated_ and donor_acceptor_groups_intrares_hbonds_accepted_ lists.
/// @note This version is for the SYMMETRIC case.  It detects interresidue hbonds between this residue and its symmetry mates.  This should be
/// called after detect_intra_residue_hydrogen_bonds().
void
BuriedUnsatPenaltyNode::detect_intra_residue_hydrogen_bonds_symmetric(
	core::Size const res_index,
	core::conformation::ResidueCOP rotamer,
	core::pose::Pose const &pose,
	core::conformation::symmetry::SymmetryInfoCOP const symminfo,
	core::scoring::hbonds::HBondDatabase const &hbond_database,
	core::scoring::hbonds::HBondOptions const &hbond_options,
	BuriedUnsatPenaltyNodeDataOP new_stored_data
) {
	runtime_assert( symminfo->bb_is_independent(res_index) );
	bool const rotamer_matches_pose( rotamer == pose.conformation().residue_cop(res_index) );

	utility::vector1< core::Size > &donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted(new_stored_data->donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted() );
	utility::vector1< core::Size > &donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated(new_stored_data->donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated() );
	donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted.clear();
	donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated.clear();
	core::Size const ngroups( new_stored_data->donor_acceptor_groups().size() );
	donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted.resize( ngroups, 0 );
	donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated.resize( ngroups, 0 );

	// Get the symmetric conformation.  Note: no type-checking except in debug mode.
	debug_assert( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( pose.conformation_ptr() ) != nullptr );
	core::conformation::symmetry::SymmetricConformationCOP symmconf( utility::pointer::static_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( pose.conformation_ptr() ) );

	// Get mirror symmetry info, if available.  Will be nullptr if not a mirror symmetric pose.
	core::conformation::symmetry::MirrorSymmetricConformationCOP mirrorsymmconf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::MirrorSymmetricConformation const >( pose.conformation_ptr() ) /*Will be nullptr if this is symmetric but not mirror symmetric.*/ );

	utility::vector1< core::Size > symm_copies( symminfo->bb_clones(res_index) );
	for ( core::Size i(1), imax(symm_copies.size()); i<=imax; ++i ) {
		core::Size const symm_copy( symm_copies[i] );
		debug_assert ( symm_copy != res_index );
		core::conformation::ResidueCOP rotamer_copy(nullptr);
		if ( rotamer_matches_pose ) {
			rotamer_copy = pose.conformation().residue_cop( symm_copy );
		} else {
			core::conformation::ResidueOP rotamer_copy_nonconst = ( mirrorsymmconf != nullptr && mirrorsymmconf->res_is_mirrored( symm_copy ) ? rotamer->clone_flipping_chirality( *( symmconf->residue_type_set_for_conf( rotamer->type().mode() ) ) ) : rotamer->clone() );
			for ( core::Size ia(1), iamax(rotamer_copy_nonconst->natoms()); ia<=iamax; ++ia ) {
				rotamer_copy_nonconst->set_xyz(ia, symmconf->apply_transformation_norecompute( rotamer->xyz(ia), res_index, symm_copy ) );
			}
			rotamer_copy_nonconst->update_actcoord();
			rotamer_copy_nonconst->seqpos(symm_copy); //for correct node indexing later
			rotamer_copy = rotamer_copy_nonconst;
		}
		core::scoring::hbonds::HBondSet hbondset( hbond_options, *rotamer, *rotamer_copy, hbond_database );
		if ( hbondset.nhbonds() > 0 ) {
			for ( core::Size j(1), jmax( hbondset.nhbonds() ); j<=jmax; ++j ) {
				core::scoring::hbonds::HBond const & curhbond( hbondset.hbond(j) );
				runtime_assert( curhbond.acc_res() == res_index || curhbond.acc_res() == symm_copy );
				runtime_assert( curhbond.don_res() == res_index || curhbond.don_res() == symm_copy );
				if ( curhbond.acc_res() == curhbond.don_res() ) continue; //We cover intraresidue hbonds elsewhere.
				//Note: since the residue identity is the same between acceptor and donor, we can be more casual about the atom index lookups.  The symmetry
				//means that it doesn't matter which residue index is the acceptor and which is the donor.
				bool const acc_group_is_original( curhbond.acc_res() == res_index );
				core::Size const acc_group( new_stored_data->group_heavyatom_index_to_group_index().at( curhbond.acc_atm() ) );
				core::Size const don_group( new_stored_data->group_heavyatom_index_to_group_index().at( rotamer->icoor( curhbond.don_hatm() ).stub_atom1().atomno() ) );
				if ( acc_group_is_original ) {
					++(new_stored_data->donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_accepted()[acc_group]);
				} else {
					++(new_stored_data->donor_acceptor_groups_equivalent_res_symmetry_copy_hbonds_donated()[don_group]);
				}
			}
		}
	}
}

/*******************************************************************************
class BuriedUnsatPenaltyEdgeData
*******************************************************************************/


/// @brief Constructor.
BuriedUnsatPenaltyEdgeData::BuriedUnsatPenaltyEdgeData() :
	utility::pointer::ReferenceCount(),
	hbonds_list_()
{}

/// @brief Initialize a newly-created edge data object.
void
BuriedUnsatPenaltyEdgeData::add_hbond(
	bool const lower_numbered_node_is_acceptor,
	core::Size const acceptor_group,
	core::Size const donor_group,
	core::Real const hbond_energy,
	core::Size const lower_numbered_node_symmetry_copy_index,
	core::Size const higher_numbered_node_symmetry_copy_index
) {
	hbonds_list_.emplace_back( BuriedUnsatPenaltyGraphHbond( lower_numbered_node_is_acceptor, donor_group, acceptor_group, hbond_energy, lower_numbered_node_symmetry_copy_index, higher_numbered_node_symmetry_copy_index ) );
	runtime_assert( hbonds_list_[hbonds_list_.size()].acceptor_group() > 0 );
	runtime_assert( hbonds_list_[hbonds_list_.size()].donor_group() > 0 );
}

/*******************************************************************************
class BuriedUnsatPenaltyEdge
*******************************************************************************/

/// @brief Edge constructor.
/// @details Note that edge_data_ starts out as nullptr; a BuriedUnsatPenaltyEdgeData object must be created and passed in before the edge can be used!
BuriedUnsatPenaltyEdge::BuriedUnsatPenaltyEdge( utility::graph::Graph * owner, platform::Size const first_node_ind, platform::Size const second_node_ind ) :
	utility::graph::Edge( owner, first_node_ind, second_node_ind ),
	edge_data_(nullptr)
{}

/// @brief Edge copy-like constructor.
BuriedUnsatPenaltyEdge::BuriedUnsatPenaltyEdge(
	utility::graph::Graph * owner,
	BuriedUnsatPenaltyEdge const &src
) :
	utility::graph::Edge( owner, src.get_first_node_ind(), src.get_second_node_ind() ),
	edge_data_(src.edge_data_)
{}

/// @brief Destructor.
BuriedUnsatPenaltyEdge::~BuriedUnsatPenaltyEdge() = default;

/// @brief Initialize this edge from another.
void
BuriedUnsatPenaltyEdge::copy_from(
	utility::graph::Edge const * src
) {
	debug_assert( dynamic_cast< BuriedUnsatPenaltyEdge const * >( src ) != nullptr );
	BuriedUnsatPenaltyEdge const &other( *( static_cast<BuriedUnsatPenaltyEdge const *>(src) ) );
	edge_data_ = other.edge_data_;
}

/// @brief Set the data object for this edge.
/// @details Data object pointer is copied; object is NOT cloned.
void
BuriedUnsatPenaltyEdge::set_edge_data(
	BuriedUnsatPenaltyEdgeDataCOP edge_data_in
) {
	debug_assert( edge_data_in != nullptr );
	edge_data_ = edge_data_in;
}


/*******************************************************************************
class BuriedUnsatPenaltyGraph
*******************************************************************************/

/// @brief Options constructor.
/// @details Note: this stores the owning pointer to the options; it doesn't clone them.  The hbond options *are* cloned and modified, though.
BuriedUnsatPenaltyGraph::BuriedUnsatPenaltyGraph( BuriedUnsatPenaltyGraphOptionsCOP options, core::scoring::hbonds::HBondOptionsCOP hbond_options ) :
	parent( ),
	always_rotamer_one_(false),
	residuepos_rotamerindex_to_nodeindex_(),
	residue_memory_address_to_nodeindex_(),
	nodeindex_to_residue_memory_address_(),
	options_(options),
	hbond_options_(nullptr),
	bunsat_edge_pool_( new boost::unordered_object_pool< BuriedUnsatPenaltyEdge >( 256 ) )
{
	core::scoring::hbonds::HBondOptionsOP hbond_options_copy( new core::scoring::hbonds::HBondOptions(*hbond_options) );
	configure_hbond_options( *hbond_options_copy );
	hbond_options_ = core::scoring::hbonds::HBondOptionsCOP( hbond_options_copy );
	if ( TR.Debug.visible() ) {
		options->show(TR.Debug);
		TR.Debug.flush();
	}
}

/// @brief Nodecount constructor with options.
/// @details Note: this stores the owning pointer to the options; it doesn't clone them.  The hbond options *are* cloned and modified, though.
BuriedUnsatPenaltyGraph::BuriedUnsatPenaltyGraph(
	platform::Size const num_nodes,
	BuriedUnsatPenaltyGraphOptionsCOP options,
	core::scoring::hbonds::HBondOptionsCOP hbond_options
) :
	parent( ),
	always_rotamer_one_(false),
	residuepos_rotamerindex_to_nodeindex_(),
	residue_memory_address_to_nodeindex_(),
	nodeindex_to_residue_memory_address_(),
	options_(options),
	hbond_options_(nullptr),
	bunsat_edge_pool_( new boost::unordered_object_pool< BuriedUnsatPenaltyEdge >( 256 ) )
{
	set_num_nodes(num_nodes);
	core::scoring::hbonds::HBondOptionsOP hbond_options_copy( new core::scoring::hbonds::HBondOptions(*hbond_options) );
	configure_hbond_options( *hbond_options_copy );
	hbond_options_ = core::scoring::hbonds::HBondOptionsCOP( hbond_options_copy );
	if ( TR.Debug.visible() ) {
		options->show(TR.Debug);
		TR.Debug.flush();
	}
}

/// @brief Copy constructor.
BuriedUnsatPenaltyGraph::BuriedUnsatPenaltyGraph(
	BuriedUnsatPenaltyGraph const &src
) :
	parent( ),
	always_rotamer_one_(src.always_rotamer_one_),
	residuepos_rotamerindex_to_nodeindex_(src.residuepos_rotamerindex_to_nodeindex_),
	residue_memory_address_to_nodeindex_(src.residue_memory_address_to_nodeindex_),
	nodeindex_to_residue_memory_address_(src.nodeindex_to_residue_memory_address_),
	options_(src.options_),
	hbond_options_(src.hbond_options_),
	bunsat_edge_pool_( new boost::unordered_object_pool< BuriedUnsatPenaltyEdge >( 256 ) )
{
	parent::operator = ( src );
}


/// @brief Destructor.
BuriedUnsatPenaltyGraph::~BuriedUnsatPenaltyGraph() {
	parent::delete_everything();
	delete bunsat_edge_pool_;
	bunsat_edge_pool_ = nullptr;
}

//////////////////////
// Public, static helper functions:
//////////////////////

/// @brief Set up options for hbond detection.
void
BuriedUnsatPenaltyGraph::configure_hbond_options(
	core::scoring::hbonds::HBondOptions &hbondoptions
) {
	hbondoptions.use_hb_env_dep(false);
	hbondoptions.use_hb_env_dep_DNA(false);
	hbondoptions.smooth_hb_env_dep(false);
	hbondoptions.bb_donor_acceptor_check(false);
	hbondoptions.exclude_DNA_DNA(false);
	hbondoptions.exclude_self_hbonds(false);
	hbondoptions.exclude_intra_res_protein(false);
	hbondoptions.exclude_intra_res_RNA(false);
}

//////////////////////
// Public member functions:
//////////////////////

/// @brief Needed override from base class.
void
BuriedUnsatPenaltyGraph::delete_edge(
	utility::graph::Edge * edge
) {
	debug_assert( dynamic_cast< BuriedUnsatPenaltyEdge * >( edge ) != nullptr );
	bunsat_edge_pool_->destroy( static_cast< BuriedUnsatPenaltyEdge * >( edge ) );
}

/// @brief Needed override from base class.
void
BuriedUnsatPenaltyGraph::delete_node(
	utility::graph::Node * node
) {
	debug_assert( dynamic_cast< BuriedUnsatPenaltyNode * >( node ) != nullptr );
	delete static_cast< BuriedUnsatPenaltyNode * >( node );
}

/// @brief Set whether this is a graph that just stores one rotamer per position.
void
BuriedUnsatPenaltyGraph::set_always_rotamer_one(
	bool const setting
) {
	always_rotamer_one_ = setting;
}

/// @brief Provide Pymol commands to colour the pose grey, non-buried donor and acceptor groups cyan, and buried acceptor
/// and donor groups orange.  Useful for debugging degree of burial.
/// @details To use, pass in a pose.  If this graph contains residues corresponding to those in the pose, commands for colouring
/// them will be written out.
void
BuriedUnsatPenaltyGraph::provide_pymol_commands_to_show_groups(
	std::ostream &out,
	core::pose::Pose const &pose
) const {
	out << "Pymol commands to colour the pose:\n";
	out << "color grey\n";
	core::conformation::Conformation const &conf( pose.conformation() );
	for ( core::Size ir(1), irmax(pose.total_residue()); ir<=irmax; ++ir ) {
		core::conformation::ResidueCOP curres( conf.residue_cop(ir) );
		if ( has_node_corresponding_to_residue( curres ) ) {
			BuriedUnsatPenaltyNode const &curnode( *( static_cast< BuriedUnsatPenaltyNode const * >( get_node( get_node_index( curres ) ) ) ) );
			for ( core::Size igroup(1), igroupmax( curnode.num_donor_acceptor_groups() ); igroup<=igroupmax; ++igroup ) {
				BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup const &curgroup( curnode.donor_acceptor_group(igroup) );
				std::string const cur_colour( curgroup.is_counted() ? "orange" : "cyan" );
				out << "color " << cur_colour << ", resi " << ir << " AND name " << utility::strip( curres->atom_name( curgroup.heavyatom_index() ) );
				for ( core::Size ih(1), ihmax( curgroup.proton_indices().size() ); ih<=ihmax; ++ih ) {
					out << "+" << utility::strip( curres->atom_name( curgroup.proton_indices()[ih] ) );
				}
				out << "\n";
			}
		}
	}
}

/// @brief Initialize a BuriedUnsatPenaltyGraph from a pose, for scoring.
void
BuriedUnsatPenaltyGraph::initialize_graph_for_scoring(
	core::pose::Pose const &pose
) {
	//I THINK that this will work if I just pass an empty RotamerSets object to the other initializer.
	core::pack::rotamer_set::RotamerSets dummy_rotamer_set;
	initialize_graph_for_packing( pose, dummy_rotamer_set, true );
}

/// @brief Initialize a BuriedUnsatPenaltyGraph from a pose and a residue set, for packing.
void
BuriedUnsatPenaltyGraph::initialize_graph_for_packing(
	core::pose::Pose const &pose,
	core::pack::rotamer_set::RotamerSets const &rotamersets,
	bool const only_scoring /*=false*/
) {

	bool const is_symmetric( core::pose::symmetry::is_symmetric(pose) );

	// Get symmetry info, if available
	core::conformation::symmetry::SymmetricConformationCOP const symmconf( is_symmetric ? utility::pointer::static_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( pose.conformation_ptr() ) : nullptr  );
	core::conformation::symmetry::SymmetryInfoCOP const symminfo( is_symmetric ? symmconf->Symmetry_Info() : nullptr );

	// Get mirror symmetry info, if available.
	core::conformation::symmetry::MirrorSymmetricConformationCOP mirrorsymmconf( nullptr );
	if ( symminfo != nullptr ) {
		debug_assert( symmconf != nullptr ); //Needs to be true.
		mirrorsymmconf = utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::MirrorSymmetricConformation const >( symmconf ); //Will be nullptr if this is symmetric but not mirror symmetric.
	}

	core::Size const nres( is_symmetric ? symminfo->num_independent_residues() : pose.total_residue() );

	core::scoring::hbonds::HBondDatabaseCOP hbdata( core::scoring::hbonds::HBondDatabase::get_database( hbond_options_->params_database_tag() ) );

	utility::vector1< bool > input_rotamer_is_in_rotamerset( nres, false ); //Is the current rotamer in the rotamerset?

	// First figure out the number of nodes that we'll need in the graph.
	{ //Scope for count
		core::Size counter(0);
		for ( core::Size i(1); i<=nres; ++i ) {
			if ( !only_scoring && rotamersets.has_rotamer_set_for_residue(i) ) {
				core::pack::rotamer_set::RotamerSet const &rotset( *(rotamersets.rotamer_set_for_residue(i)) );
				counter += rotset.num_rotamers();
				core::conformation::ResidueCOP curposeres( pose.conformation().residue_cop(i) );
				for ( core::Size j(1), jmax(rotset.num_rotamers()); j<=jmax; ++j ) {
					core::conformation::ResidueCOP currotamer( rotset.rotamer(j) );
					if ( currotamer == curposeres ) {
						input_rotamer_is_in_rotamerset[i] = true;
						break;
					}
				}
				if ( !input_rotamer_is_in_rotamerset[i] ) ++counter;
			} else {
				++counter;
			}
		}
		set_num_nodes( counter );
	} //End scope for count

	// Second, since nodes contain essential information, set up each node.
	{ //Scope for node setup
		core::Size counter(0);
		for ( core::Size i(1); i<=nres; ++i ) {
			bool const has_rotset( !only_scoring && rotamersets.has_rotamer_set_for_residue(i) );
			for ( core::Size j(1), jmax( has_rotset ? rotamersets.rotamer_set_for_residue(i)->num_rotamers() + static_cast<core::Size>( !input_rotamer_is_in_rotamerset[i] ) : 1 ); j<=jmax; ++j ) {
				++counter;
				//The following stores the position, rotamer index, and all hydrogen bond donors and acceptors in the node:
				//It also stores intra-residue hydrogen bonds.
				initialize_node(
					counter,
					i,
					j,
					(has_rotset ?
					( (!input_rotamer_is_in_rotamerset[i] && j == jmax) ? pose.conformation().residue_cop(i) : rotamersets.rotamer_set_for_residue(i)->rotamer(j))
					:
					pose.conformation().residue_cop(i))
					, pose,
					*hbdata,
					is_symmetric
				);
			}
		}
		debug_assert( counter == num_nodes() );
	}



	// Third, consider inter-residue hydrogen bonds.  We did intraresidue hydrogen bonds separately, above:
	{ //Scope for inter-residue hydrogen bonds.
		core::Real const hbond_energy_threshold( options_->hbond_energy_threshold() );
		for ( core::Size i(1); i<=nres; ++i ) { //Loop through all residues in pose or asymmetric unit, for first of pair to compare
			bool const i_has_rotamers( !only_scoring && rotamersets.has_rotamer_set_for_residue(i) );
			for ( core::Size j(1); j<i; ++j ) { //Loop through all residues in pose or asymmetric unit, for second of pair to compare
				bool const j_has_rotamers( !only_scoring && rotamersets.has_rotamer_set_for_residue(j) );
				for ( core::Size ii(1), iimax( i_has_rotamers ? rotamersets.rotamer_set_for_residue(i)->num_rotamers() + static_cast<core::Size>( !input_rotamer_is_in_rotamerset[i] ) : 1 ); ii<=iimax; ++ii ) { //Loop through all rotamers for position i
					core::conformation::Residue const & res1( i_has_rotamers && ( input_rotamer_is_in_rotamerset[i] || ii != iimax ) ? *(rotamersets.rotamer_set_for_residue(i)->rotamer(ii)) : pose.residue(i) );
					for ( core::Size jj(1), jjmax( j_has_rotamers ? rotamersets.rotamer_set_for_residue(j)->num_rotamers() + static_cast<core::Size>( !input_rotamer_is_in_rotamerset[j] ) : 1 ); jj<=jjmax; ++jj ) { //Loop through all rotamers for position j
						for ( core::Size jsymm(1), jsymmmax( is_symmetric ? symminfo->num_bb_clones() + 1 : 1 ); jsymm<=jsymmmax; ++jsymm ) { //Loop through all symmetry copies of the second residue.
							core::Size const jsymmindex( is_symmetric ? ( jsymm < jsymmmax ? symminfo->bb_clones(j)[jsymm] : j ) : j );
							core::conformation::ResidueCOP res2;
							if ( !is_symmetric ) {
								if ( j_has_rotamers  && ( input_rotamer_is_in_rotamerset[j] || jj != jjmax ) ) {
									res2 = rotamersets.rotamer_set_for_residue(j)->rotamer(jj);
								} else {
									res2 = pose.conformation().residue_cop(j);
								}
							} else { // if is_symmetric
								if ( j_has_rotamers && ( input_rotamer_is_in_rotamerset[j] || jj != jjmax ) ) {
									core::conformation::ResidueCOP residue_to_clone( rotamersets.rotamer_set_for_residue(j)->rotamer(jj) );
									core::conformation::ResidueOP newres( mirrorsymmconf != nullptr && mirrorsymmconf->res_is_mirrored( jsymmindex ) ? residue_to_clone->clone_flipping_chirality( *( symmconf->residue_type_set_for_conf( residue_to_clone->type().mode() ) ) ) : residue_to_clone->clone() );
									if ( jsymmindex != j ) {
										for ( core::Size ia(1), iamax(newres->natoms()); ia<=iamax; ++ia ) {
											newres->set_xyz(ia, symmconf->apply_transformation_norecompute( residue_to_clone->xyz(ia), j, jsymmindex ) );
										}
										newres->update_actcoord();
									}
									newres->seqpos(j); //for correct node indexing later
									res2 = newres; //nonconst to const
								} else {
									core::conformation::ResidueOP newres( pose.conformation().residue_cop(jsymmindex)->clone() );
									newres->seqpos(j); //for correct node indexing later
									res2 = newres; //nonconst to const
								}
							}

							if ( res1.nbr_atom_xyz().distance_squared(res2->nbr_atom_xyz()) <= std::pow(2.0*res1.nbr_radius()+2.0*res2->nbr_radius(), 2) ) {
								core::scoring::hbonds::HBondSet hbondset( *hbond_options_, res1, *res2, *hbdata );
								if ( hbondset.nhbonds() > 0 ) {
									BuriedUnsatPenaltyEdge * newly_added_edge(
										get_edge_exists( i, ii, j, jj )
										?
										find_edge(i, ii, j, jj)
										:
										add_edge(i, ii, j, jj)
									);
									BuriedUnsatPenaltyEdgeDataOP newly_added_edge_data( utility::pointer::make_shared< BuriedUnsatPenaltyEdgeData >() );
									newly_added_edge->set_edge_data(newly_added_edge_data); // Can continue to access the data using the OP above; COP is now stored.
									for ( core::Size ibond(1), ibondmax(hbondset.nhbonds()); ibond<=ibondmax; ++ibond ) {
										core::scoring::hbonds::HBond const &curhbond( hbondset.hbond(ibond) );
										if ( curhbond.energy() <= hbond_energy_threshold ) {
											runtime_assert( curhbond.acc_res() == i || curhbond.acc_res() == j );
											runtime_assert( curhbond.don_res() == i || curhbond.don_res() == j );
											if ( curhbond.acc_res() == curhbond.don_res() ) continue; //Intraresidue hbonds are handled elsewhere.
											BuriedUnsatPenaltyNode * node1( static_cast< BuriedUnsatPenaltyNode*>( get_node( get_node_index(i, ii) ) ) );
											BuriedUnsatPenaltyNode * node2( static_cast< BuriedUnsatPenaltyNode*>( get_node( get_node_index(j, jj) ) ) );
											bool const lower_numbered_node_is_acceptor( curhbond.acc_res() == std::min( i, j ) );
											bool const node1_is_acceptor( curhbond.acc_res() == i );
											debug_assert( lower_numbered_node_is_acceptor || curhbond.acc_res() == std::max( i, j ) );
											core::Size const acceptor_group_index( node1_is_acceptor ? node1->get_donor_acceptor_group_from_heavyatom_index(curhbond.acc_atm()) : node2->get_donor_acceptor_group_from_heavyatom_index(curhbond.acc_atm()) );
											core::Size const don_hatm_parent( node1_is_acceptor ? res2->icoor( curhbond.don_hatm() ).stub_atom1().atomno() : res1.icoor( curhbond.don_hatm() ).stub_atom1().atomno() );
											core::Size const donor_group_index( node1_is_acceptor ? node2->get_donor_acceptor_group_from_heavyatom_index(don_hatm_parent) : node1->get_donor_acceptor_group_from_heavyatom_index(don_hatm_parent) );
											newly_added_edge_data->add_hbond( lower_numbered_node_is_acceptor, acceptor_group_index, donor_group_index, curhbond.energy(), is_symmetric ? (jsymm < jsymmmax ? jsymm + 1 : 1) : 1, 1 );
											if ( TR.Debug.visible() ) {
												TR.Debug << "Added hbond between position " << i << ", rotamer " << ii << " and position " << j << ", rotamer " << jj;
												if ( is_symmetric ) {
													TR.Debug << " symmcopy " << (jsymm < jsymmmax ? jsymm + 1 : 1);
												}
												TR.Debug << "." << std::endl;
											}
										}
									} //Loop through all hydrogen bonds between pair
								}
							}

						} //Loop through all symmetry copies of the second residue.
					} //Loop through all rotamers in 2nd res
				} //Loop through all rotamers in 1st res
			} //Inner loop through all residues in pose
		} //Outer loop through all residues in pose
	} //Scope for inter-residue hydrogen bonds
}


/// @brief Given this BuriedUnsatPenaltyGraph with some number of nodes, iterate through each node and update the
/// internally-stored counts for unsats and oversats based on the edges connected to that node.
/// @details Calls compute_unsats_for_node();
void
BuriedUnsatPenaltyGraph::compute_unsats_all_nodes() {
	for ( core::Size i(1), imax( num_nodes() ); i<=imax; ++i ) {
		compute_unsats_for_node( i );
	}
}

/// @brief Given two lists (one of changed nodes, one of their partners), update the internally-stored counts for unsats and oversats for
/// those nodes only.
/// @details calls compute_unsats_for_node().
void
BuriedUnsatPenaltyGraph::compute_unsats_changed_nodes(
	utility::vector1< core::Size > const & changed_node_indices,
	utility::vector1< core::Size > const & changed_node_partners
) {
	for ( core::Size i(1), imax(changed_node_indices.size()); i<=imax; ++i ) {
		compute_unsats_for_node( changed_node_indices[i] );
	}
	for ( core::Size j(1), jmax( changed_node_partners.size()); j<=jmax; ++j ) {
		compute_unsats_for_node( changed_node_partners[j] );
	}
}

/// @brief Given this BuriedUnsatPenaltyGraph with some number of nodes and the index of a node, update the
/// internally-stored counts for unsats and oversats based on the edges connected to that node.
void
BuriedUnsatPenaltyGraph::compute_unsats_for_node(
	core::Size const node_index
) {
	BuriedUnsatPenaltyNode & curnode( *( static_cast< BuriedUnsatPenaltyNode * >( get_node( node_index ) ) ) );
	curnode.clear_hbond_counts();

	for ( utility::graph::EdgeListConstIterator it( curnode.const_edge_list_begin() ); it != curnode.const_edge_list_end(); ++it ) {
		BuriedUnsatPenaltyEdge const & curedge( *( static_cast< BuriedUnsatPenaltyEdge const * >( *it ) ) );
		bool const this_node_is_first_in_edge( curedge.get_first_node_ind() == node_index );
		debug_assert( this_node_is_first_in_edge || curedge.get_second_node_ind() == node_index); //Should be true.
		for ( core::Size i(1), imax(curedge.n_hbonds()); i<=imax; ++i ) {
			BuriedUnsatPenaltyGraphHbond const & curhbond( curedge.hbond(i) );
			bool const this_node_is_acceptor( this_node_is_first_in_edge ? curhbond.first_node_is_the_acceptor() : !curhbond.first_node_is_the_acceptor() );
			if ( this_node_is_acceptor ) {
				core::Size const acceptor_group_index( curhbond.acceptor_group() );
				BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup const & acceptor_group( curnode.donor_acceptor_group(acceptor_group_index) );
				if ( curhbond.acceptor_symmetry_copy_index() == 1 || curhbond.donor_symmetry_copy_index() == 1 ) {
					if ( acceptor_group.is_counted() ) curnode.increment_accepted_hbond_count_for_group(acceptor_group_index);
				}
			} else {
				core::Size const donor_group_index( curhbond.donor_group() );
				BuriedUnsatPenaltyGraphHbondDonorAcceptorGroup const & donor_group( curnode.donor_acceptor_group(donor_group_index) );
				if ( curhbond.donor_symmetry_copy_index() == 1 || curhbond.acceptor_symmetry_copy_index() == 1 ) {
					if ( donor_group.is_counted() ) curnode.increment_donated_hbond_count_for_group(donor_group_index);
				}
			}
		} //Loop through hbonds for this edge.
	} //Loop through edges
}

/// @brief Given an index of a node in this graph, an owning pointer to another graph, and a node index in the other graph, copy the node
/// from the other graph to the node in this graph, flush the edges that were connected to the node in this graph, and copy those edges from
/// the other graph that can be connected to nodes in this graph.
/// @details Note that the logic for determining whether an edge from the other graph can be related to this graph is based on ResidueCOPs.
void
BuriedUnsatPenaltyGraph::copy_node_and_connected_edges(
	core::Size const node_index_in_this_graph,
	BuriedUnsatPenaltyGraph const & other_graph,
	core::Size const node_index_in_other_graph
) {
	debug_assert( dynamic_cast< BuriedUnsatPenaltyNode * >( this->get_node( node_index_in_this_graph ) ) != nullptr );
	BuriedUnsatPenaltyNode & this_node( *( static_cast< BuriedUnsatPenaltyNode * >( this->get_node( node_index_in_this_graph ) ) ) );
	this_node.drop_all_edges(); //Delete all edges in the current node.

	// Erase entries in current maps:
	if ( nodeindex_to_residue_memory_address_.count( node_index_in_this_graph ) ) {
		core::conformation::ResidueCOP curres( nodeindex_to_residue_memory_address_.at(node_index_in_this_graph) );
		residue_memory_address_to_nodeindex_.erase( curres );
		nodeindex_to_residue_memory_address_.erase( node_index_in_this_graph );
	}
	std::pair< core::Size, core::Size > residuepos( this_node.residue_position(), always_rotamer_one_ ? 1 : this_node.rotamer_index() );
	if ( residuepos_rotamerindex_to_nodeindex_.count( residuepos ) ) {
		residuepos_rotamerindex_to_nodeindex_.erase( residuepos );
	}

	utility::graph::Node const * other_node_raw( other_graph.get_node( node_index_in_other_graph ) );
	debug_assert( dynamic_cast< BuriedUnsatPenaltyNode const * >( other_node_raw ) != nullptr );
	BuriedUnsatPenaltyNode const & other_node( *( static_cast< BuriedUnsatPenaltyNode const * >( other_node_raw ) ) );
	this_node.copy_from( other_node_raw );

	//Update entries in current maps:
	core::conformation::ResidueCOP new_res( other_graph.nodeindex_to_residue_memory_address(node_index_in_other_graph) );
	nodeindex_to_residue_memory_address_[ node_index_in_this_graph ] = new_res;
	residue_memory_address_to_nodeindex_[ new_res ] = node_index_in_this_graph;
	residuepos_rotamerindex_to_nodeindex_[ std::pair< core::Size, core::Size >( this_node.residue_position(), always_rotamer_one_ ? 1 : this_node.rotamer_index() ) ] = node_index_in_this_graph;

	for ( utility::graph::EdgeListConstIterator it(other_node.const_edge_list_begin()); it!=other_node.const_edge_list_end(); ++it ) {
		debug_assert( dynamic_cast< BuriedUnsatPenaltyEdge const * >( *it ) != nullptr );
		BuriedUnsatPenaltyEdge const &curedge( *( static_cast< BuriedUnsatPenaltyEdge const * >( *it ) ) );
		core::Size const second_node_index_in_other_graph( curedge.get_other_ind( node_index_in_other_graph ) ); //The index of the node to which this edge connects the node in the other graph.
		core::conformation::ResidueCOP res_pointer( other_graph.nodeindex_to_residue_memory_address_.at(second_node_index_in_other_graph) ); //Get the memory address of the residue corresponding to that node.
		if ( this->residue_memory_address_to_nodeindex_.find(res_pointer) != this->residue_memory_address_to_nodeindex_.end() ) { //If we also have a node corresponding to that residue...
			core::Size const second_node_index_in_this_graph( this->residue_memory_address_to_nodeindex_.at(res_pointer) );
			this->add_edge( node_index_in_this_graph, second_node_index_in_this_graph, curedge); //...then we add an edge in this graph.
		}
	}

}

//////////////////////
// Private functions
//////////////////////

/// @brief Add an edge, representing an interresidue hydrogen bonded interaction (consisting of one or more interresidue hydrogen bonds), to the graph.
/// @details Returns a pointer to the newly-created edge.
/// @note In the graph base class, the first node index is always numbered lower than the second node index for an edge.  The swap happens automatically.
BuriedUnsatPenaltyEdge *
BuriedUnsatPenaltyGraph::add_edge(
	core::Size const seqpos1,
	core::Size const rotamer_index1,
	core::Size const seqpos2,
	core::Size const rotamer_index2
) {
	runtime_assert( seqpos1 != seqpos2 );
	core::Size const node_index1(get_node_index( seqpos1, rotamer_index1 ) );
	core::Size const node_index2(get_node_index( seqpos2, rotamer_index2 ) );
	return static_cast< BuriedUnsatPenaltyEdge * >( parent::add_edge( node_index1, node_index2 ) );
}

/// @brief Add an edge, representing an interresidue hydrogen bonded interaction (consisting of one or more interresidue hydrogen bonds), to the graph.
/// @details Returns a pointer to the newly-created edge.  This version copies an edge from another graph.
/// @note In the graph base class, the first node index is always numbered lower than the second node index for an edge.  The swap happens automatically.
BuriedUnsatPenaltyEdge *
BuriedUnsatPenaltyGraph::add_edge(
	core::Size const node_index1,
	core::Size const node_index2,
	BuriedUnsatPenaltyEdge const & other_edge
) {
	runtime_assert( node_index1 != node_index2 );
#ifndef NDEBUG
	BuriedUnsatPenaltyEdge * new_edge( static_cast< BuriedUnsatPenaltyEdge * >( utility::graph::Graph::add_edge( node_index1, node_index2 ) ) );
#else
	BuriedUnsatPenaltyEdge * new_edge( dynamic_cast< BuriedUnsatPenaltyEdge * >( utility::graph::Graph::add_edge( node_index1, node_index2 ) ) );
	debug_assert( new_edge != nullptr );
#endif
	new_edge->copy_from(&other_edge);
	return new_edge;
}

/// @brief Determine whether an edge exists in the graph, by seqpos and rotamer index.
bool
BuriedUnsatPenaltyGraph::get_edge_exists(
	core::Size const seqpos1,
	core::Size const rotamer_index1,
	core::Size const seqpos2,
	core::Size const rotamer_index2
) const {
	core::Size const node1( get_node_index(seqpos1, rotamer_index1) );
	core::Size const node2( get_node_index(seqpos2, rotamer_index2) );
	return parent::get_edge_exists( node1, node2 );
}

/// @brief Retrieve an edge by seqpos and rotamer index.
/// @details Nonconst version.
BuriedUnsatPenaltyEdge *
BuriedUnsatPenaltyGraph::find_edge(
	core::Size const seqpos1,
	core::Size const rotamer_index1,
	core::Size const seqpos2,
	core::Size const rotamer_index2
) {
	core::Size const node1( get_node_index(seqpos1, rotamer_index1) );
	core::Size const node2( get_node_index(seqpos2, rotamer_index2) );
	debug_assert( dynamic_cast< BuriedUnsatPenaltyEdge * >( parent::find_edge(node1, node2) ) != nullptr  );
	return static_cast< BuriedUnsatPenaltyEdge * >( parent::find_edge(node1, node2) );
}

/// @brief Retrieve an edge by seqpos and rotamer index.
/// @details Const version.
BuriedUnsatPenaltyEdge const *
BuriedUnsatPenaltyGraph::find_edge(
	core::Size const seqpos1,
	core::Size const rotamer_index1,
	core::Size const seqpos2,
	core::Size const rotamer_index2
) const {
	core::Size const node1( get_node_index(seqpos1, rotamer_index1) );
	core::Size const node2( get_node_index(seqpos2, rotamer_index2) );
	debug_assert( dynamic_cast< BuriedUnsatPenaltyEdge const * >( parent::find_edge(node1, node2) ) != nullptr  );
	return static_cast< BuriedUnsatPenaltyEdge const * >( parent::find_edge(node1, node2) );
}

/// @brief factory method for node creation, defined by derived graph
/// classes, called by the base class
utility::graph::Node*
BuriedUnsatPenaltyGraph::create_new_node(
	platform::Size node_index
) {
	return new BuriedUnsatPenaltyNode( this, node_index );
}

/// @brief factory method for edge creation, defined by derived graph
/// classes, called by the base class
utility::graph::Edge*
BuriedUnsatPenaltyGraph::create_new_edge(
	platform::Size index1,
	platform::Size index2
) {
	return bunsat_edge_pool_->construct( this, index1, index2 );
}

/// @brief This is also needed for edge creation, when copying graphs.
utility::graph::Edge *
BuriedUnsatPenaltyGraph::create_new_edge(
	utility::graph::Edge const * example_edge
) {
	return bunsat_edge_pool_->construct(
		this,
		static_cast< BuriedUnsatPenaltyEdge const & > (*example_edge)
	);
}

/// @brief Initialize a node to represent a given residue position and rotamer index, and store in it all of the relevant hydrogen bond donors and acceptors.
/// @param[in] node_index The index of the node to initialize.
/// @param[in] residue_position The index of the residue that this node represents in the pose.
/// @param[in] rotamer_index The index of the rotamer for this residue that this node represents.
/// @param[in] residue The residue object itself, for extracting hbond donor/acceptor information.
/// @param[in] pose The pose, for context.
/// @param[in] hb_data The hydrogen bonding database object.
/// @param[in] is_symmetric Is this a symmetric pose?  (We figure this out once and only once, to avoid repeated dynamic_casts.)
void
BuriedUnsatPenaltyGraph::initialize_node(
	core::Size const node_index,
	core::Size const residue_position,
	core::Size const rotamer_index,
	core::conformation::ResidueCOP residue,
	core::pose::Pose const &pose,
	core::scoring::hbonds::HBondDatabase const &hb_data,
	bool const is_symmetric
) {
	core::conformation::symmetry::SymmetricConformationCOP const symmconf( is_symmetric ? utility::pointer::static_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( pose.conformation_ptr() ) : nullptr  );
	core::Size const nres( is_symmetric ? symmconf->Symmetry_Info()->num_independent_residues() : pose.total_residue() );

	runtime_assert( residue_position > 0 && residue_position <= nres );

	//Store the (position, rotamer index) pair as a map key pointing to the current node:
	std::pair< core::Size, core::Size > const position_rotamer( residue_position, always_rotamer_one_ ? 1 : rotamer_index );
	runtime_assert( residuepos_rotamerindex_to_nodeindex_.find( position_rotamer ) == residuepos_rotamerindex_to_nodeindex_.end() );
	residuepos_rotamerindex_to_nodeindex_[position_rotamer] = node_index;
	runtime_assert( residue_memory_address_to_nodeindex_.find( residue ) == residue_memory_address_to_nodeindex_.end() );
	residue_memory_address_to_nodeindex_[ residue ] = node_index;
	runtime_assert( nodeindex_to_residue_memory_address_.find(node_index) == nodeindex_to_residue_memory_address_.end() );
	nodeindex_to_residue_memory_address_[ node_index ] = residue;

	//Get and initialize the current node:
	BuriedUnsatPenaltyNode* cur_node(static_cast< BuriedUnsatPenaltyNode * >( this->get_node(node_index) ));
	cur_node->initialize_node( residue_position, rotamer_index, residue, pose, options_, hbond_options_, hb_data, is_symmetric );
}

} //graph
} //buried_unsat_penalty
} //guidance_scoreterms
} //pack
} //core
