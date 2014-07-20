// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/Pose.cc
/// @brief  Pose class
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov

// Unit headers
#include <core/pose/Pose.hh>

// Package headers
#include <core/pose/util.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/pose/PDB_Info.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/pose/metrics/PoseMetricContainer.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/RingConformerSet.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/id/NamedAtomID.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/io/pdb/file_data.hh>

// Basic headers
#include <basic/init.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/prof.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/vector0.hh>

// C++ headers
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/unordered_map.hpp>


namespace core {
namespace pose {

/// @details default init function
void Pose::init(void)
{
	#ifdef PYROSETTA
		// Sanity check: check if core::init was called already and abort otherwise with helpful message...
		if( !basic::was_init_called() ) utility_exit_with_message("Attempt to initialize Pose object before core::init was called detected… Have you forgot to call core::init?");
	#endif

	conformation_ = new Conformation();

	// have the Pose observe it's Conformation for XYZ changes
	// we discard the Link because we own the Conformation
	conformation_->attach_xyz_obs( &Pose::on_conf_xyz_change, this );

	energies_ = new scoring::Energies();
	energies_->set_owner( this );

	data_cache_ = new BasicDataCache( datacache::CacheableDataType::num_cacheable_data_types );

	observer_cache_ = new ObserverCache( datacache::CacheableObserverType::num_cacheable_data_types, *this );

	metrics_ = new metrics::PoseMetricContainer;
	metrics_->attach_to( *this );
}


/// @details default constructor
Pose::Pose() :
	pdb_info_( NULL ),
	constraint_set_( 0 )
{
	init();
}

/// @details destructor -- > kill data on the heap
Pose::~Pose()
{
	//std::cout << "Pose dstor" << std::endl;
	notify_destruction_obs( DestructionEvent( this ) );
	clear();
}

/// @brief copy constructor
Pose::Pose( Pose const & src ) :
	ReferenceCount ( src )
{
	*this = src;
}

/// @brief partial copy constructor
Pose::Pose( Pose const & src, Size begin, Size const end):
	pdb_info_( NULL ),
	constraint_set_( 0 )
{
	init();
	utility::vector1< core::Size > residue_indices;
	for(; begin <= end; ++begin){
		residue_indices.push_back(begin);
	}
	core::io::pdb::pose_from_pose(*this, src, residue_indices);
}

/// @brief copy assignment
Pose &
Pose::operator=( Pose const & src )
{
	PROF_START ( basic::POSE_COPY );

	ConformationOP old_conf = conformation_; // track for observer transfer

	if ( conformation_ && conformation_->same_type_as_me( *src.conformation_, true ) ) {
		(*conformation_) = (*src.conformation_);
	} else {
		conformation_ = src.conformation_->clone();
		conformation_->attach_xyz_obs( &Pose::on_conf_xyz_change, this );
	}

	if ( energies_ && energies_->same_type_as_me( *src.energies_ ) ) {
		(*energies_) = (*src.energies_);
	} else {
		energies_ = src.energies_->clone();
		energies_->set_owner( this );
	}

	data_cache_  = new BasicDataCache( datacache::CacheableDataType::num_cacheable_data_types );
	*data_cache_ = *(src.data_cache_);

	observer_cache_ = new ObserverCache( datacache::CacheableObserverType::num_cacheable_data_types, *this );
	*observer_cache_ = *src.observer_cache_;

	this->pdb_info( src.pdb_info_ );

	metrics_ = new metrics::PoseMetricContainer( *src.metrics_ );
	metrics_->attach_to( *this );

	if( constraint_set_ ) constraint_set_->detach_from_conformation();

	if ( src.constraint_set_ ) {
		constraint_set_ = src.constraint_set_->clone();
		constraint_set_->attach_to_conformation( conformation_.get() );
	} else {
		constraint_set_ = 0;
	}

	num_stacks_ = src.num_stacks_;

	// transfer remaining observers that honor the Conformation::TRANSFER
	// event after everything else is done
	if ( old_conf ) {
		conformation_->receive_observers_from( *old_conf );
		old_conf.reset_to_null(); // force clear
	}

	PROF_STOP ( basic::POSE_COPY );

	return *this;
}

/// @brief clone the conformation
PoseOP
Pose::clone() const
{
	return new Pose( *this );
}

kinematics::FoldTree const &
Pose::fold_tree() const
{
	return conformation_->fold_tree();
}

void
Pose::fold_tree( kinematics::FoldTree const & fold_tree_in )
{
	conformation_->fold_tree( fold_tree_in );
}

/// @brief get the atom_tree
kinematics::AtomTree const &
Pose::atom_tree() const
{
	return conformation_->atom_tree();
}

int
Pose::chain( Size const seqpos ) const
{
	PyAssert( (seqpos<=total_residue()), "Pose::chain( Size const seqpos ): variable seqpos is out of range!" );
	return residue( seqpos ).chain();
}

/// @details  Note that we do not clone the input conformation -- we just take it directly. This could be unsafe (?)
///  but it's more efficient. Maybe we want to switch to cloning... Of course we already hand out nonconst refs to
///  our conformation, which is a little unsafe anyhow.
/// @warning Classes observing the Pose's old conformation will not automatically
///  be re-attached/listening to the new Conformation.  Please pay special attention
///  if you have a CacheableObserver in the ObserverCache that listens to
///  a Pose's Conformation.  The prior PDB_Info, ConstraintSet and Energies will be
///  cleared as well.
void
Pose::set_new_conformation( conformation::ConformationCOP new_conformation )
{

	/// drop stuff
	pdb_info_.reset_to_null();
	constraint_set_.reset_to_null();
	observer_cache_->detach();

	// initialize new
	energies_ = new scoring::Energies();
	energies_->set_owner( this );

	data_cache_ = new BasicDataCache( datacache::CacheableDataType::num_cacheable_data_types );

	metrics_ = new metrics::PoseMetricContainer;
	metrics_->attach_to( *this );

	/// clone and reassign the pointer
	conformation_ = new_conformation->clone();
	conformation_->attach_xyz_obs( &Pose::on_conf_xyz_change, this );

	/* OPERATIONS AFTER THIS POINT NEED TO HAPPEN AFTER THE CONFORMATION
     HAS BEEN SWAPPED */

}

/// @details This function allow us to attach an Energies object from a derived class. What are the proper
///	checks to do for this? This should probably be a protected function, since we do not want this function
/// to be used regularly
void
Pose::set_new_energies_object( scoring::EnergiesOP energies )
{
  energies_ = energies->clone();
	// Set the owner of the new energies object
	energies_->set_owner( this );
}

/// @details splits the current pose into several poses containing only a single chain each.
utility::vector1<PoseOP>
Pose::split_by_chain() const
{
	utility::vector1<PoseOP> singlechain_poses;

	for (core::Size i = 1; i <= conformation_->num_chains(); i++)
	{
		singlechain_poses.push_back( split_by_chain(i) );
	}

	return singlechain_poses;
}

/// @details Returns a pose containing only the given chain.
PoseOP
Pose::split_by_chain(Size const chain_id) const
{

	core::pose::PoseOP chain_pose = new Pose(*this);
	Size chain_begin, chain_end, delete_begin, delete_end;

	chain_begin = chain_pose->conformation().chain_begin( chain_id );
	chain_end = chain_pose->conformation().chain_end( chain_id );

	// if there is only one chain in the pose do nothing and return a copy of the pose
	if ((conformation_->num_chains() == 1) && (chain_id == 1))
	{

	}
	// if this is the first chain, delete chain_end to the end of the pose
	else if (chain_begin == 1) {
		delete_begin = chain_end + 1;
		delete_end = chain_pose->total_residue();
		chain_pose->conformation().delete_residue_range_slow( delete_begin, delete_end );
	}
	// if this is the last chain, delete the start of the pose to chain_begin
	else if ( chain_end == chain_pose->total_residue() ) {
		delete_begin = 1;
		delete_end = chain_begin - 1;
		chain_pose->conformation().delete_residue_range_slow( delete_begin, delete_end );
	}
	// otherwise, make two deletes around the chain of interest
	else {
		delete_begin = 1;
		delete_end = chain_begin - 1;
		chain_pose->conformation().delete_residue_range_slow( delete_begin, delete_end );
		// just deleted residues --> renumbering pose, so need to reset deletion mask
		delete_begin = chain_end - chain_begin + 2;
		delete_end = chain_pose->total_residue();
		chain_pose->conformation().delete_residue_range_slow( delete_begin, delete_end );
	}

	// disulfides
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	if ( option[ in::detect_disulf ].user() ?
			option[ in::detect_disulf ]() : // detect_disulf true
			chain_pose->is_fullatom() ) // detect_disulf default but fa pose
	{
		chain_pose->conformation().detect_disulfides();
	}

	// restore broken pdb_info to new pose ~ Labonte
	chain_pose->pdb_info()->obsolete(false);

	return chain_pose;
}


// TODO: Move to util.cc.
/// @details  This method updates the pose chain IDs to match the chain IDs
/// found in pdb_info().  In some applications, it is more intuitive to change
/// pdb chain ID letters than it is to change pose chain IDs.  This method
/// adds chain endings between pdb chains and re-derives the pose chain IDs.
/// Currently, disconnected segments with the same pdb chain ID character are
/// treated as separate pose chains, e.g., it is possible for pose chains 1,
/// 3, and 5 to all be chain X.  In the future, I will add a flag to force a
/// one-to-one correspondence between the two chain designations, e.g., if
/// residues 6 through 10 are chain B and residues 1 through 5 AND residues 11
/// through 15 are chain A, then the pose will be reordered to place all res-
/// idues with the same pdb chain ID into a single pose chain. (This is how it
/// works when a pose is loaded from a pdb file.)  I personally have needed
/// use of both functionalities.  I have chosen to create this as a separate
/// method, rather than a part of a change_pdb_chain_ID_of_range(), to avoid
/// multiple calls to Conformation.rederive_chain_ids().  Thus, this method
/// should be called once after all modifications to pdb_info() have been
/// made. ~ Labonte
///
/// See also:
///  PDB_Info.chain()
///  PDB_Info.set_resinfo()
///  Pose.split_by_chain()
///  Conformation.rederive_chain_IDs()
///  Conformation.rederive_chain_endings()
void
Pose::update_pose_chains_from_pdb_chains()
{
	// Declare a vector for storing new (between-residue) chain endings.
	utility::vector1<Size> new_endings;

	char last_pdb_chain = pdb_info_->chain(1);
	char current_pdb_chain;

	for (Size i = 1; i <= conformation_->size(); ++i) {
		current_pdb_chain = pdb_info_->chain(i);
		if (current_pdb_chain != last_pdb_chain) {
			new_endings.push_back(i - 1);
			last_pdb_chain = current_pdb_chain;
		}

	// (chain_endings() includes a call to Conformer.rederive_chain_IDs.)
	conformation_->chain_endings(new_endings);
	}
}


void
Pose::metric( std::string const & calculator_name, std::string const & key, basic::MetricValueBase & val ) const
{ metrics_->get(calculator_name, key, val, *this); return; }

std::string
Pose::print_metric( std::string const & calculator_name, std::string const & key ) const
{ return metrics_->print(calculator_name, key, *this); }

void
Pose::append_residue_by_jump(
	conformation::Residue const & new_rsd,
	Size const jump_anchor_residue,
	std::string const& jump_anchor_atom, // = "",
	std::string const& jump_root_atom, // = "",
	bool const start_new_chain //= false
)
{
	//PyAssert( (jump_anchor_residue>0) && (jump_anchor_residue<total_residue()), "Pose::append_residue_by_jump( ...Size const jump_anchor_residue... ): variable jump_anchor_residue is out of range!" );    // check later: may be fixed in conformation
	energies_->clear(); // TEMPORARY
	conformation_->append_residue_by_jump( new_rsd, jump_anchor_residue, jump_anchor_atom, jump_root_atom, start_new_chain );
}

void
Pose::append_residue_by_bond(
	conformation::Residue const & new_rsd,
	bool const build_ideal_geometry, // = false,
	int const connection, // = 0,
	Size const anchor_residue, // = 0,
	int const anchor_connection, // = 0,
	bool const start_new_chain, // = false
	bool const lookup_bond_length // = false
)
{
	//PyAssert( (anchor_residue>0) && (anchor_residue<=total_residue()), "Pose::append_residue_by_bond( ...Size const anchor_residue... ): variable anchor_residue is out of range!" );    // check later: may be fixed in conformation
	energies_->clear(); // TEMPORARY
	conformation_->append_residue_by_bond( new_rsd, build_ideal_geometry, connection, anchor_residue, anchor_connection, start_new_chain, lookup_bond_length);
}

void
Pose::insert_residue_by_jump(
	Residue const & new_rsd_in,
	Size const seqpos, // desired seqpos of new_rsd
	Size anchor_pos, // in the current sequence numbering, ie before insertion of seqpos
	std::string const& anchor_atomno, // = "",
	std::string const& root_atomno // = ""
)
{
	//PyAssert( (seqpos>0), "Pose::insert_residue_by_jump( ...Size const seqpos... ): variable seqpos is out of range!" );    // check later:
	PyAssert( (anchor_pos<=total_residue()), "Pose::insert_residue_by_jump( ...Size anchor_pos... ): variable anchor_pos is out of range!" );    // check later:
	energies_->clear(); // TEMPORARY
	conformation_->insert_residue_by_jump( new_rsd_in, seqpos, anchor_pos, anchor_atomno, root_atomno );
}

void
Pose::insert_residue_by_bond(
                             Residue const & new_rsd_in,
                             Size const seqpos, // desired seqpos of new_rsd
                             Size anchor_pos, // in the current sequence numbering, ie before insertion of seqpos
                             bool const build_ideal_geometry, // = false,
                             std::string const& anchor_atom, // could be ""
                             std::string const& root_atom, // ditto
                             bool new_chain,
                             bool const lookup_bond_length // default false
)
{
    //PyAssert( (seqpos>0), "Pose::insert_residue_by_jump( ...Size const seqpos... ): variable seqpos is out of range!" );    // check later:
    PyAssert( (anchor_pos<=total_residue()), "Pose::insert_residue_by_jump( ...Size anchor_pos... ): variable anchor_pos is out of range!" );    // check later:
    energies_->clear(); // TEMPORARY
    conformation_->insert_residue_by_bond( new_rsd_in, seqpos, anchor_pos, build_ideal_geometry, anchor_atom, root_atom, new_chain, lookup_bond_length );
}

void
Pose::replace_residue(
	Size const seqpos,
	Residue const & new_rsd_in,
	bool const orient_backbone
)
{
	PyAssert( (seqpos<=total_residue()), "Pose::replace_residue( ...Size const seqpos... ): variable seqpos is out of range!" );    // check later: may become unecessary
	conformation_->replace_residue( seqpos, new_rsd_in, orient_backbone );
}

void
Pose::replace_residue(
	int const seqpos,
	Residue const & new_rsd_in,
	utility::vector1< std::pair< std::string, std::string > > const & atom_pairs
)
{
	PyAssert( (seqpos<=static_cast<int>(total_residue())), "Pose::replace_residue( ...Size const seqpos... ): "
			"variable seqpos is out of range!" );    // check later: may become unnecessary
	conformation_->replace_residue( seqpos, new_rsd_in, atom_pairs );
}

void
Pose::append_polymer_residue_after_seqpos(
	Residue const & new_rsd,
	Size const seqpos,
	bool const build_ideal_geometry
)
{
	PyAssert( (seqpos<=total_residue()), "Pose::append_polymer_residue_after_seqpos( ...Size const seqpos... ): variable seqpos is out of range!" );    // check later: may become unecessary
	energies_->clear(); // TEMPORARY
	conformation_->append_polymer_residue_after_seqpos( new_rsd, seqpos, build_ideal_geometry );
}

void
Pose::prepend_polymer_residue_before_seqpos(
	Residue const & new_rsd,
	Size const seqpos,
	bool const build_ideal_geometry
)
{
	PyAssert( (seqpos<=total_residue()), "Pose::prepend_polymer_residue_before_seqpos( ...Size const seqpos... ): variable seqpos is out of range!" );    // check later:
	energies_->clear(); // TEMPORARY
	conformation_->prepend_polymer_residue_before_seqpos( new_rsd, seqpos, build_ideal_geometry );
}

void
Pose::append_pose_by_jump(
	Pose const & src,
	Size const jump_anchor_residue,
	std::string const & jump_anchor_atom,
	std::string const & jump_root_atom)
{
	core::Size old_n_residue = n_residue();

	conformation().insert_conformation_by_jump(
			src.conformation(),
			n_residue() + 1,
			num_jump() + 1,
			jump_anchor_residue,
			0, // Set default anchor jump number
			jump_anchor_atom,
			jump_root_atom);

	if(pdb_info().get() != NULL && src.pdb_info().get() != NULL) {
		pdb_info()->copy(
				*src.pdb_info(),
				1,
				src.pdb_info()->nres(),
				old_n_residue + 1);
	}
}

void
Pose::delete_polymer_residue( Size const seqpos )
{
	PyAssert( (seqpos<=total_residue()), "Pose::delete_polymer_residue( Size const seqpos ): variable seqpos is out of range!" );
	energies_->clear(); // TEMPORARY
	conformation_->delete_polymer_residue( seqpos );
}


void
Pose::copy_segment(
	Size const size,
	Pose const & src,
	Size const begin,
	Size const src_begin
)
{
	conformation_->copy_segment( size, src.conformation(), begin, src_begin );
	// now copy any other data
}


Size
Pose::total_residue() const
{
	return conformation_->size();
}

Size
Pose::n_residue() const
{
	return conformation_->size();
}

bool
Pose::empty() const
{
	return conformation_->empty();
}

Size
Pose::num_jump() const
{
	return conformation_->fold_tree().num_jump();
}

chemical::AA
Pose::aa( Size const seqpos ) const
{
	PyAssert( (seqpos<=total_residue()), "Pose::aa( Size const seqpos ): variable seqpos is out of range!" );
	return conformation_->aa( seqpos );
}

char
Pose::secstruct( Size const seqpos ) const
{
	PyAssert( (seqpos<=total_residue()), "Pose::secstruct( Size const seqpos ): variable seqpos is out of range!" );
	return conformation_->secstruct(seqpos);
}

std::string
Pose::secstruct() const {
	std::string ss="";
	for ( Size i = 1; i <= total_residue(); i++ )	{
		ss += secstruct( i );
	}
	return ss;
}

void
Pose::set_secstruct( Size const seqpos, char const setting )
{
	PyAssert( (seqpos<=total_residue()), "Pose::set_secstruct( Size const seqpos , char const setting ): variable seqpos is out of range!" );
	// check variable "setting" to ensure it is logical?
	conformation_->set_secstruct( seqpos, setting );
}


// Sequence accessors
std::string
Pose::sequence() const
{
	std::string seq;
	for ( Size i=1; i<= conformation_->size(); ++i ) {
		seq += residue(i).name1();
	}
	return seq;
}

std::string
Pose::annotated_sequence( bool show_all_variants ) const
{
	using namespace core::chemical;

	std::string seq;
	for ( Size i=1; i<= conformation_->size(); ++i ) {
		char c = residue(i).name1();
		std::string rn = residue( i ).name();
		bool not_cyd = (rn.find( "CYD" ) == std::string::npos); // awatkins: captures all cyd variations
		seq += c;
		if ( ( !oneletter_code_specifies_aa(c) || name_from_aa( aa_from_oneletter_code(c) ) != residue(i).name() )
		  && ( show_all_variants || not_cyd ) ) {
			seq = seq + '[' + residue(i).name() + ']';
		}
	}
	return seq;
}

std::string
Pose::chain_sequence(core::Size const chain_in) const
{
	using namespace std;

	assert(chain_in <= conformation_->num_chains());
	PyAssert((chain_in <= conformation_->num_chains()),
			"Pose::chain_sequence(core::Size const chain_in): variable chain_in is out of range!");

	stringstream seq(stringstream::out);

	Size begin = conformation_->chain_begin(chain_in);
	Size end = conformation_->chain_end(chain_in);

	if (!residue(begin).is_carbohydrate()) {
		for (Size i = begin; i <= end; ++i) {
			seq << residue(i).name1();
		}
	} else /*is carbohydrate*/ {
		// Carbohydrate sequences are listed in the opposite direction as they are numbered.
		for (Size i = end; i >= begin; --i) {
			seq << residue(i).carbohydrate_info()->short_name();
			if (i != begin) {
				seq << "(";
				seq << residue(i).carbohydrate_info()->anomeric_carbon();
			}
		}
	}
	return seq.str();
}


// Residue at position seqpos
Pose::Residue const &
Pose::residue(
	Size const seqpos
) const
{
	PyAssert( (seqpos<=total_residue()), "Pose::residue( Size const seqpos ): variable seqpos is out of range!" );
	return conformation_->residue( seqpos );
}

chemical::ResidueType const &
Pose::residue_type(
	Size const seqpos
) const
{
	PyAssert( (seqpos<=total_residue()), "Pose::residue_type( Size const seqpos ): variable seqpos is out of range!" );
	return conformation_->residue_type( seqpos );
}


// backbone torsions
// peptides and saccharides

/// @details  For proteins and peptoids, phi is defined as C(n-1)-N(n)-CA(n)-C(n).\n
/// For aldopyranoses, phi is defined as O5(n)-C1(n)-OX(n-1)-CX(n-1),
/// where X is the position of the glycosidic linkage.\n
/// For aldofuranoses, phi is defined as O4(n)-C1(n)-OX(n-1)-CX(n-1).\n
/// For 2-ketopyranoses, phi is defined as O6(n)-C2(n)-OX(n-1)-CX(n-1).\n
/// For 2-ketofuranoses, phi is defined as O5(n)-C2(n)-OX(n-1)-CX(n-1).\n
/// Et cetera...
Real
Pose::phi( Size const seqpos ) const
{
	using namespace id;

	assert( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() );
	PyAssert( (seqpos<=total_residue()), "Pose::phi( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() ),
		"Pose::phi( Size const seqpos ): residue seqpos is not part of a protein, peptoid, or carbohydrate!" );

	if ( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() ) {
		return residue(seqpos).mainchain_torsion(phi_torsion);
	} else /*is carbohydrate*/ {
		return carbohydrates::get_glycosidic_torsion(phi_torsion, *this, seqpos);
	}
}

/// @details  For proteins, phi is defined as C(n-1)-N(n)-CA(n)-C(n).\n
/// For aldopyranoses, phi is defined as O5(n)-C1(n)-OX(n-1)-CX(n-1),
/// where X is the position of the glycosidic linkage.\n
/// For aldofuranoses, phi is defined as O4(n)-C1(n)-OX(n-1)-CX(n-1).\n
/// For 2-ketopyranoses, phi is defined as O6(n)-C2(n)-OX(n-1)-CX(n-1).\n
/// For 2-ketofuranoses, phi is defined as O5(n)-C2(n)-OX(n-1)-CX(n-1).\n
/// Et cetera...
void
Pose::set_phi( Size const seqpos, Real const setting )
{
	using namespace id;

	assert( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() );
	PyAssert( (seqpos<=total_residue()), "Pose::set_phi( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate()),
		"Pose::set_phi( Size const seqpos , Real const setting ): residue seqpos is not part of a protein, peptoid, or carbohydrate!" );

	if ( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() ) {
		conformation_->set_torsion( TorsionID( seqpos, BB, phi_torsion ), setting );
	} else /*is carbohydrate*/ {
		carbohydrates::set_glycosidic_torsion(phi_torsion, *this, seqpos, setting);
	}
}

/// @details  For proteins, psi is defined as N(n)-CA(n)-C(n)-N(n+1).\n
/// For saccharides, psi is defined as: C(anomeric)(n)-OX(n-1)-CX(n-1)-CX-1(n-1),\n
/// where X is the position of the glycosidic linkage.
Real
Pose::psi( Size const seqpos ) const
{
	using namespace id;

	assert( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() );
	PyAssert( (seqpos<=total_residue()), "Pose::psi( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate()),
		"Pose::psi( Size const seqpos ): residue seqpos is not part of a protein, peptoid, or carbohydrate!" );

	if ( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() ) {
		return residue(seqpos).mainchain_torsion(psi_torsion);
	} else /*is carbohydrate*/ {
		return carbohydrates::get_glycosidic_torsion(psi_torsion, *this, seqpos);
	}
}

/// @details  For proteins, psi is defined as N(n)-CA(n)-C(n)-N(n+1).\n
/// For saccharides, psi is defined as: C(anomeric)(n)-OX(n-1)-CX(n-1)-CX-1(n-1),\n
/// where X is the position of the glycosidic linkage.
void
Pose::set_psi( Size const seqpos, Real const setting )
{

	using namespace id;

	assert( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() );
	PyAssert( (seqpos<=total_residue()), "Pose::set_psi( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate()),
		"Pose::set_psi( Size const seqpos , Real const setting ): residue seqpos is not part of a protein, peptoid, or carbohydrate!" );

	if ( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() ) {
		conformation_->set_torsion( TorsionID( seqpos, BB, psi_torsion ), setting);
	} else /*is carbohydrate*/ {
		carbohydrates::set_glycosidic_torsion(psi_torsion, *this, seqpos, setting);
	}
}

/// @details  For proteins, omega is defined as CA(n)-C(n)-N(n+1)-CA(n+1).
/// For carbohydrates glycosylated at an exocyclic position,
/// omega of residue n is defined as OX(n-1)-CX(n-1)-CX-1(n-1)-CX-2(n-1),
/// where X is the position of the glycosidic linkage.  (Note that every atom
/// defining this torsion comes from the previous residue!)
Real Pose::omega( Size const seqpos ) const
{
	using namespace id;

	assert( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() );
	PyAssert( (seqpos<=total_residue()), "Pose::omega( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() ),
		"Pose::omega( Size const seqpos ): residue seqpos is not part of a protein,peptoid, or carbohydrate!" );

	if ( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() ) {
		return residue(seqpos).mainchain_torsion(omega_torsion);
	} else /*is carbohydrate*/ {
		return carbohydrates::get_glycosidic_torsion(omega_torsion, *this, seqpos);
	}
}

/// @details  For proteins, omega is defined as CA(n)-C(n)-N(n+1)-CA(n+1).
/// For carbohydrates glycosylated at an exocyclic position,
/// omega of residue n is defined as OX(n-1)-CX(n-1)-CX-1(n-1)-OX-1(n-1),
/// where X is the position of the glycosidic linkage.  (Note that every atom
/// defining this torsion comes from the previous residue!)
void
Pose::set_omega( Size const seqpos, Real const setting )
{
	using namespace id;

	assert( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() );
	PyAssert( (seqpos<=total_residue()), "Pose::set_omega( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() ),
		"Pose::set_omega( Size const seqpos , Real const setting ): residue seqpos is not part of a protein, peptoid, or carbohydrate!" );

	if ( residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() ) {
		conformation_->set_torsion( TorsionID( seqpos, BB, omega_torsion ),  setting);
	} else /*is carbohydrate*/ {
		carbohydrates::set_glycosidic_torsion(omega_torsion, *this, seqpos, setting);
	}
}

// nucleic acids

Real
Pose::alpha( Size const seqpos ) const{
	assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=total_residue()), "Pose::alpha( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::alpha( Size const seqpos ): residue seqpos is not part of a Nucleic Acid!" );
	return torsion( id::TorsionID( seqpos, id::BB, 1 ) );
}

void
Pose::set_alpha( Size const seqpos, Real const setting )
{
	assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=total_residue()), "Pose::set_alpha( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::set_alpha( Size const seqpos , Real const setting ): residue seqpos is not part of a Nucleic Acid!" );
	conformation_->set_torsion( TorsionID( seqpos, id::BB, 1 ), setting );
}

Real
Pose::beta( Size const seqpos ) const{
	assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=total_residue()), "Pose::beta( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::beta( Size const seqpos ): residue seqpos is not part of a Nucleic Acid!" );
	return torsion( id::TorsionID( seqpos, id::BB, 2 ) );
}

void
Pose::set_beta( Size const seqpos, Real const setting )
{
	assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=total_residue()), "Pose::set_beta( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::set_beta( Size const seqpos , Real const setting ): residue seqpos is not part of a Nucleic Acid!" );
	conformation_->set_torsion( TorsionID( seqpos, id::BB, 2 ), setting );
}

Real
Pose::gamma( Size const seqpos ) const{
	assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=total_residue()), "Pose::gamma( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::gamma( Size const seqpos ): residue seqpos is not part of a Nucleic Acid!" );
	return torsion( id::TorsionID( seqpos, id::BB, 3 ) );
}

void
Pose::set_gamma( Size const seqpos, Real const setting )
{
	assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=total_residue()), "Pose::set_gamma( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::set_gamma( Size const seqpos , Real const setting ): residue seqpos is not part of a Nucleic Acid!" );
	conformation_->set_torsion( TorsionID( seqpos, id::BB, 3 ), setting );
}

Real
Pose::delta( Size const seqpos ) const{
	assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=total_residue()), "Pose::delta( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::delta( Size const seqpos ): residue seqpos is not part of a Nucleic Acid!" );
	return torsion( id::TorsionID( seqpos, id::BB, 4 ) );
}

void
Pose::set_delta( Size const seqpos, Real const setting )
{
	assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=total_residue()), "Pose::set_delta( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::set_delta( Size const seqpos , Real const setting ): residue seqpos is not part of a Nucleic Acid!" );
	conformation_->set_torsion( TorsionID( seqpos, id::BB, 4 ), setting );
}

Real
Pose::epsilon( Size const seqpos ) const{
	assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=total_residue()), "Pose::epsilon( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::epsilon( Size const seqpos ): residue seqpos is not part of a Nucleic Acid!" );
	return torsion( id::TorsionID( seqpos, id::BB, 5 ) );
}

void
Pose::set_epsilon( Size const seqpos, Real const setting )
{
	assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=total_residue()), "Pose::set_epsilon( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::set_epsilon( Size const seqpos , Real const setting ): residue seqpos is not part of a Nucleic Acid!" );
	conformation_->set_torsion( TorsionID( seqpos, id::BB, 5 ), setting );
}

Real
Pose::zeta( Size const seqpos ) const{
	assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=total_residue()), "Pose::zeta( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::zeta( Size const seqpos ): residue seqpos is not part of a Nucleic Acid!" );
	return torsion( id::TorsionID( seqpos, id::BB, 6 ) );
}

void
Pose::set_zeta( Size const seqpos, Real const setting )
{
	assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=total_residue()), "Pose::set_zeta( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::set_zeta( Size const seqpos , Real const setting ): residue seqpos is not part of a Nucleic Acid!" );
	conformation_->set_torsion( TorsionID( seqpos, id::BB, 6 ), setting );
}

// sidechain torsions
// peptides and saccharides

Real
Pose::chi(
	int const chino,
	Size const seqpos
) const
{
	PyAssert( (seqpos<=total_residue()), "Pose::chi( int const chino , Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_protein() || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() ),
		"Pose::chi( int const chino , Size const seqpos ): residue seqpos is not part of a protein, peptoid, or carbohydrate!" );
	PyAssert( (chino>0) && (chino<=static_cast<int>(residue(seqpos).nchi())),
		"Pose::chi( int const chino , Size const seqpos ): variable chino innappropriate for this residue!" );

	return residue( seqpos ).chi( chino );
}

void
Pose::set_chi(
	int const chino,
	Size const seqpos,
	Real const setting
)
{
	PyAssert( (seqpos<=total_residue()), "Pose::set_chi( int const chino , Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_protein()  || residue_type(seqpos).is_peptoid() || residue_type(seqpos).is_carbohydrate() || residue_type(seqpos).is_ligand() ),
			"Pose::set_chi( int const chino , Size const seqpos , Real const setting ): residue seqpos is not part of a protein, peptoid, ligand or carbohydrate!" );
	PyAssert( (chino>0) && (chino<=static_cast<int>(residue(seqpos).nchi())),
		"Pose::set_chi( int const chino , Size const seqpos ): variable chino innappropriate for this residue!" );

	conformation_->set_torsion( TorsionID(seqpos, id::CHI, chino), setting);
}

// nucleic acids

Real
Pose::chi( Size const seqpos ) const{
	assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=total_residue()), "Pose::chi( Size const seqpos ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::chi( Size const seqpos ): residue seqpos is not part of a Nucleic Acid!" );
	return torsion( id::TorsionID( seqpos, id::CHI, 1 ) );
}

void
Pose::set_chi( Size const seqpos, Real const setting )
{
	assert( residue_type( seqpos ).is_NA() );
	PyAssert( (seqpos<=total_residue()), "Pose::set_chi( Size const seqpos, Real const setting ): variable seqpos is out of range!" );
	PyAssert( (residue_type(seqpos).is_NA()), "Pose::set_chi( Size const seqpos , Real const setting ): residue seqpos is not part of a Nucleic Acid!" );
	conformation_->set_torsion( TorsionID( seqpos, id::CHI, 1 ), setting );
}

// Set the given residue's ring conformation, if appropriate.
void
Pose::set_ring_conformation(uint const seqpos, core::chemical::RingConformer const & conformer)
{
	using namespace std;
	using namespace id;
	using namespace numeric;

	Residue const & res = residue(seqpos);

	assert(res.type().is_cyclic());
	PyAssert((seqpos<=total_residue()),
			"Pose::set_ring_conformation(uint const seqpos, core::chemical::RingConformer const &"
			"conformer): variable seqpos is out of range!");
	PyAssert((res.type().is_cyclic()),
			"Pose::set_ring_conformation(uint const seqpos, core::chemical::RingConformer const &"
			"conformer): residue seqpos is not a cyclic residue!" );

	// First, set the nus, which DEFINE the ideal ring conformer.
	Size n_nus = res.type().nu_atoms().size();
	for (uint i = 1; i <= n_nus; ++i) {
		set_torsion(TorsionID(seqpos, NU, i), conformer.nu_angles[i]);
	}

	// Then, set the taus, which result from ring strain.
	Size n_taus = n_nus + 1;  // There will always be one more tau stored in the conformer than nus.
	for (uint i = 1; i < n_taus; ++i) {
		// The reference atoms for the bond angle can be extracted from those used for the corresponding nu angle.
		// For example, nu2 is defined as C1-C2-C3-C4, and tau 1 is defined as C1-C2-C3.
		AtomID ref1(res.type().nu_atoms(i)[1], seqpos);
		AtomID ref2(res.type().nu_atoms(i)[2], seqpos);
		AtomID ref3(res.type().nu_atoms(i)[3], seqpos);
		conformation_->set_bond_angle(ref1, ref2, ref3, conversions::radians(conformer.tau_angles[i]));
	}

	// Since one fewer nus are stored than taus, we need the LAST 3 reference atoms from the last nu, instead of the
	// 1st 3 atoms as we used above.
	AtomID ref1(res.type().nu_atoms(n_nus)[2], seqpos);
	AtomID ref2(res.type().nu_atoms(n_nus)[3], seqpos);
	AtomID ref3(res.type().nu_atoms(n_nus)[4], seqpos);
	conformation_->set_bond_angle(ref1, ref2, ref3, conversions::radians(conformer.tau_angles[n_taus]));

	// Finally, fix the virtual alignment in the special case of saccharide residues.
	if (res.is_carbohydrate()) {
		carbohydrates::align_virtual_atoms_in_carbohydrate_residue(*this, seqpos);  // TODO: Where does this belong?
	}
}

/////////////////////////////////////////////////////////////////////////////
// generic torsion-angle access

/// @brief  get the torsion angle identified by id
Real
Pose::torsion( TorsionID const & id ) const
{
	return conformation_->torsion( id );
}

/// @brief  set the torsion angle identified by id
void
Pose::set_torsion( TorsionID const & id, Real const setting )
{
	conformation_->set_torsion( id, setting );
}


/////////////////////////////////////////////////////////////////////////////
// jumps

///
void
Pose::set_jump(
	int const jump_number,
	const kinematics::Jump & new_jump
)
{
	conformation_->set_jump( jump_number, new_jump );
}

///
kinematics::Jump const &
Pose::jump( int const jump_number ) const
{
	return conformation_->jump( jump_number );
}

///
void
Pose::set_jump(
	AtomID const & id,
	const kinematics::Jump & new_jump
)
{
	conformation_->set_jump( id, new_jump );
}

///
kinematics::Jump const &
Pose::jump( AtomID const & id ) const
{
	return conformation_->jump( id );
}

///////////////////////////////////////////////////////////////////////////
// access atomtree dof's

/// get the value of the atomtree degree of freedom (DOF)
Real
Pose::dof( DOF_ID const & id ) const
{
	return conformation_->dof( id );
}


/// set the value of the atomtree degree of freedom (DOF)
void
Pose::set_dof( DOF_ID const & id, Real const setting )
{
	conformation_->set_dof( id, setting );
}

bool
Pose::has_dof(
	DOF_ID const & did
) const
{
	if( Size(did.rsd()) > total_residue() || did.rsd() < 1 ) return false;
	if( Size(did.atomno()) > residue(did.rsd()).natoms() || did.atomno() < 1 ) return false;
	if( id::PHI == did.type() || id::THETA == did.type() || id::D == did.type() )
		if( conformation_->atom_tree().atom(did.atom_id()).is_jump() ) return false;
	// TODO SHEFFLER MAKE THIS RIGHT!!!!!
	return true;
}


/// get the location of an atom
PointPosition const &
Pose::xyz( AtomID const & id ) const
{
	return conformation_->xyz( id );
}

/// get the location of an atom
PointPosition const &
Pose::xyz( NamedAtomID const & id ) const
{
	return conformation_->residue(id.rsd()).xyz(id.atom());
}

/// set the location of an atom
void
Pose::set_xyz( AtomID const & id, PointPosition const & point )
{
	conformation_->set_xyz( id, point );
}

/// set the location of an atom
void
Pose::set_xyz(
	NamedAtomID const & id,
	PointPosition const & point
)
{
	conformation_->set_xyz( named_atom_id_to_atom_id(id, *this), point );
}

/// set the locations of a vector of atoms
void
Pose::batch_set_xyz( utility::vector1< AtomID > const & ids, utility::vector1< PointPosition > const & points )
{
	conformation_->batch_set_xyz( ids, points );
}


/// get the locations of a vector of atoms
void
Pose::batch_get_xyz( utility::vector1< AtomID > const & ids, utility::vector1< PointPosition > & points ) const
{
	conformation_->batch_get_xyz( ids, points );
}

kinematics::Stub
Pose::stub_from_id(
	id::NamedStubID const& id
)
{
	return conformation_->stub_from_id( named_stub_id_to_stub_id( id, *this ) );
}

void
Pose::center()
{

	PointPosition cog(0,0,0);

	Size count=0;
	for( Size ir = 1; ir <= total_residue(); ir++){
		for( Size at = 1; at <= residue( ir ).natoms(); at++){
			cog += xyz( AtomID( at, ir ) );
			count++;
		}
	}

	//std::cout << cog.x() << std::endl;
	cog /= (Real) count;

	//std::cout << cog.x() << std::endl;

	for( Size ir = 1; ir <= total_residue(); ir++){
		for( Size at = 1; at <= residue( ir ).natoms(); at++){
			set_xyz(  AtomID( at, ir ),  xyz( AtomID( at, ir )) - cog  );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
/// @details transfers domain map information into the Energies object, and then
/// resets the domain map information from the Conformation object
void
Pose::update_residue_neighbors()
{
	runtime_assert( total_residue() != 0 ); // Avoid crash on empty Pose (e.g. if PDB occupancy is all zero)
	residue( total_residue() ); // completely unnecessary temporary hack to force refold if nec.

	if ( conformation_->structure_moved() ) {
		energies_->structure_has_moved( total_residue() );
		conformation_->reset_structure_moved();
	} else if ( energies_->residue_neighbors_updated() ) {
		return;
	}

	// figure out what's changed since the residue neighbor update
	kinematics::DomainMap domain_map( total_residue() );
	conformation_->update_domain_map( domain_map );

	energies_->update_residue_neighbors( domain_map, *this );
	if ( energies_->discard_conformation_domain_map() ) conformation_->reset_move_data();
}

///////////////////////////////////////////////////////////////////////////////
///@details called by the ScoreFunction at the start of scoring.  If the score
/// function has changed since the last round of scoring, then cached energies
/// may have become invalidated -- the Energies object makes that decision.
void
Pose::scoring_begin(
	scoring::ScoreFunction const & sfxn
)
{
	// notify of any structure changes
	if ( conformation_->structure_moved() ) {
		conformation_->reset_structure_moved();
		energies_->structure_has_moved( total_residue() );
	}

	energies_->scoring_begin( sfxn, *this );
	// figure out what's changed since the last score evaluation
	/// and update the Energies object
	update_residue_neighbors();
}

void
Pose::scoring_end( scoring::ScoreFunction const & scorefxn )
{
	// reset the data that describes what has moved
	//conformation_->reset_move_data();
	energies_->scoring_end( scorefxn );
	notify_energy_obs( EnergyEvent( this ) );
}

/// @brief called by PairEPotential
void
Pose::update_actcoords()
{
	conformation_->update_actcoords();
}

void
Pose::update_actcoord( Size resid )
{
	conformation_->update_actcoord( resid );
}

void
Pose::update_orbital_coords( Size resid )
{
	conformation_->update_orbital_coords( resid );
}


/// @brief Applies a transform of the form Rx + v, where R is a rotation
/// matrix, V is a vector, and x is the original position in xyz space.
void
Pose::apply_transform_Rx_plus_v(
	numeric::xyzMatrix< Real > const & R,
	Vector const & v
) {
	for ( Size i = 1; i <= total_residue(); ++i ) {
		for ( Size j = 1; j <= residue_type(i).natoms(); ++j ) {
			AtomID id( j, i );
			set_xyz( id, R * xyz(id) + v );
			//apply_transform_Rx_plus_v( R, v );
		}
	}
}

void
Pose::clear()
{
	conformation_->clear();
	energies_->clear();
	constraint_set_ = 0;
	metrics_->clear();
	data_cache_->clear();
	observer_cache_->clear();
	pdb_info( NULL ); // will check for existence, remove observers, etc
}

///////////////////////////////////////////////////////////////////////////////
/// @details save pose data to file with supplied file_name
///
bool
Pose::dump_pdb(std::string const &file_name, std::string const & tag) const
{
	return core::io::pdb::FileData::dump_pdb(*this, file_name, tag);
}

/// @brief  Dump a pdbfile with some score info at the end.
void
Pose::dump_scored_pdb(
	std::string const &file_name,
	scoring::ScoreFunction const & scorefxn,
	std::string const & tag
) {
	Real const total_score( scorefxn( *this ) );
	/// Now handled automatically.  scorefxn.accumulate_residue_total_energies( *this );
	std::ofstream out( file_name.c_str() );
	core::io::pdb::FileData::dump_pdb( *this, out, tag );
	// verbose output
	out << "END\n";
	std::string secstruct;
	for ( Size i=1; i<= total_residue(); ++i ) secstruct += conformation_->secstruct(i);
	out << "SS: " << secstruct << '\n';
	out << "SCORE_INFO:\n";
	out << "TOTAL_SCORE: " << total_score << '\n';
	scoring::EnergyMap const & wts( scorefxn.weights() );
	out << "WTS: " << wts.show_nonzero() << '\n';
	out << "TOTAL_WTD: " << this->energies().total_energies().weighted_string_of( wts ) << '\n';
	for ( Size i=1; i<= total_residue(); ++i ) {
		out << "RSD_WTD: " << i << ' ' << this->energies().residue_total_energies( i ).weighted_string_of( wts ) << '\n';
	}
	scorefxn.show( out );
	out.close();
}

void Pose::dump_pdb(std::ostream & out, std::string const & tag) const
{
	return core::io::pdb::FileData::dump_pdb(*this, out, tag);
}

/// @brief for writing a specified subset of residues in pdb format
void
Pose::dump_pdb(
	std::ostream & out,
	utility::vector1< Size > const & residue_indices,
	std::string const & tag
) const
{
	return core::io::pdb::FileData::dump_pdb( *this, out, residue_indices, tag );
}


Pose::ConstraintSetCOP
Pose::constraint_set() const
{
	if ( constraint_set_ == 0 ) {
		return new scoring::constraints::ConstraintSet; // create an empty constraint set
	}
	return constraint_set_;
}

scoring::constraints::ConstraintCOP
Pose::add_constraint( scoring::constraints::ConstraintCOP cst )
{
	energies_->clear();
	if ( constraint_set_ == 0 ) {
		constraint_set_ = new scoring::constraints::ConstraintSet; // create an empty constraint set the first time it's asked for
		constraint_set_->attach_to_conformation( conformation_.get() );
	}
	scoring::constraints::ConstraintCOP new_cst( cst->clone() );
	constraint_set_->add_constraint( new_cst );
	return( new_cst );
}

scoring::constraints::ConstraintCOPs
Pose::add_constraints( scoring::constraints::ConstraintCOPs csts )
{
	energies_->clear();
	if ( constraint_set_ == 0 ) {
		constraint_set_ = new scoring::constraints::ConstraintSet; // create an empty constraint set the first time it's asked for
		constraint_set_->attach_to_conformation( conformation_.get() );
	}
	using namespace scoring::constraints;
	ConstraintCOPs new_csts;
	for( ConstraintCOPs::const_iterator cst_it = csts.begin(); cst_it != csts.end(); ++cst_it )
		new_csts.push_back( (*cst_it)->clone() );
	constraint_set_->add_constraints( new_csts );
	return( new_csts );
}

bool
Pose::remove_constraint(
	scoring::constraints::ConstraintCOP cst,
	bool object_comparison )
{
	if ( constraint_set_ == 0 ) return false;
	energies_->clear();
	return constraint_set_->remove_constraint(cst, object_comparison);
}

bool
Pose::remove_constraints(
	scoring::constraints::ConstraintCOPs csts,
	bool object_comparison )
{
	if ( constraint_set_ == 0 ) return false;
	energies_->clear();
	return constraint_set_->remove_constraints(csts, object_comparison);
}

bool
Pose::remove_constraints(){
	if ( constraint_set_ == 0 ) return false;
	energies_->clear();
	constraint_set_->clear();
	return true;
}

void
Pose::constraint_set( ConstraintSetOP constraint_set )
{
	energies_->clear();

	if( constraint_set_ != 0 ) constraint_set_->detach_from_conformation();

	if( constraint_set != 0 ){
		constraint_set_ = constraint_set->clone();
		constraint_set_->attach_to_conformation( conformation_.get() );
	}
	else constraint_set_ = constraint_set;
}

void Pose::transfer_constraint_set( const pose::Pose &pose ){
	energies_->clear();

	if( constraint_set_ != 0 ) constraint_set_->detach_from_conformation();

	if( pose.constraint_set_ != 0 ){
		constraint_set_ = pose.constraint_set_->clone();
		constraint_set_->attach_to_conformation( conformation_.get() );
	}
	else constraint_set_ = pose.constraint_set_;
}

void
Pose::clear_stacking_map()
{
	num_stacks_.clear();
}

void
Pose::record_stacking_interaction( core::Size const & resid ) const
{
	num_stacks_[ resid ]++;
}

const boost::unordered_map< core::Size, core::Size >&
Pose::get_stacking_map() const
{
	return num_stacks_;
}



//////////////////////////////// PDB_Info methods /////////////////////////////////////


/// @brief get pdb info (const)
/// @return NULL if no PDB_Info instance exists, the pdb info instance otherwise
PDB_InfoCOP
Pose::pdb_info() const
{
	assert( pdb_info_ ? pdb_info_->nres() == total_residue() : true );
	return pdb_info_;
}

/// @brief get pdb info
/// @return NULL if no PDB_Info instance exists, the pdb info instance otherwise
PDB_InfoOP
Pose::pdb_info()
{
	assert( pdb_info_ ? pdb_info_->nres() == total_residue() : true );
	return pdb_info_;
}

/// @brief copy new pdb info into this Pose
/// @param[in] new_info  the new pdb info to copy, pass NULL if you want to zero
///  the existence of pdb info inside this Pose
/// @return the prior pdb info instance
PDB_InfoOP
Pose::pdb_info( PDB_InfoOP new_info )
{
	if ( pdb_info_ ) {
		pdb_info_->detach_from();
	}

	PDB_InfoOP prior_pdb_info = pdb_info_;

	if ( new_info ) {
		pdb_info_ = new PDB_Info( *new_info ); // make a copy
		pdb_info_->attach_to( *conformation_ );
	} else {
		pdb_info_.reset_to_null();
	}

	PyAssert( pdb_info_ ? pdb_info_->nres() == total_residue() : true, "Invalid PDB_Info, pdb_info_->nres() != total_residue()" );

	assert( pdb_info_ ? pdb_info_->nres() == total_residue() : true );

	return prior_pdb_info;
}

bool Pose::is_fullatom() const {
	return conformation_->is_fullatom();
}

bool Pose::is_centroid() const {
	return conformation_->is_centroid();
}

/// @brief notify DestructionEvent observers
/// @remarks called only upon destruction of the Pose
void
Pose::notify_destruction_obs( DestructionEvent const & e ) {
	destruction_obs_hub_( e );
}


/// @brief notify GeneralEvent observers
/// @remarks should only be called when there are no other suitable event types
///  since specific event notifications will automatically fire a GeneralEvent signal
void
Pose::notify_general_obs( GeneralEvent const & e ) {
	general_obs_hub_( e );
}


/// @brief notify EnergyEvent observers
/// @param e the event
/// @param fire_general fire a GeneralEvent afterwards? default true
void
Pose::notify_energy_obs( EnergyEvent const & e, bool const fire_general ) {
	energy_obs_hub_( e );
	if ( fire_general ) {
		notify_general_obs( e );
	}
}

/// @brief notify ConformationEvent observers
/// @param e the event
/// @param fire_general fire a GeneralEvent afterwards? default true
void
Pose::notify_conformation_obs( ConformationEvent const & e, bool const fire_general ) {
	conformation_obs_hub_( e );
	if ( fire_general ) {
		notify_general_obs( e );
	}
}

/// @brief upon receiving a conformation::signals::XYZEvent
void
Pose::on_conf_xyz_change( core::conformation::signals::XYZEvent const & ) {
	notify_conformation_obs( ConformationEvent( this ) );
}


std::ostream & operator << ( std::ostream & os, Pose const & pose)
{
	PDB_InfoCOP p = pose.pdb_info();
	if( p ) {
		os << "PDB file name: "<< p->name() << std::endl;
	}
	os << "Total residues:" << pose.total_residue() << std::endl;
	os << "Sequence: " << pose.sequence() << std::endl;
	os << "Fold tree:" << std::endl << pose.fold_tree();
	return os;
}

#ifdef USELUA
PoseSP clone( PoseSP p ) {
	PoseSP pp( new Pose( *p ) );
	return pp;
}

void lregister_Pose( lua_State * lstate ) {
	luabind::module(lstate, "core")
	[
		luabind::namespace_("pose")
		[
			luabind::class_<Pose>("Pose")
				.def("total_residue", &Pose::total_residue)
				.def("sequence", &Pose::sequence)
				.def("clone", (PoseSP (*) (PoseSP)) &clone )
		]
	];
}
#endif

} // pose
} // core
