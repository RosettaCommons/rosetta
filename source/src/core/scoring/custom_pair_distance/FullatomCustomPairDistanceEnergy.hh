// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/FullatomCustomPairDistanceEnergy.hh
/// @brief  Custom atom pair distance energy
/// @author David E Kim


#ifndef INCLUDED_core_scoring_custom_pair_distance_FullatomCustomPairDistanceEnergy_hh
#define INCLUDED_core_scoring_custom_pair_distance_FullatomCustomPairDistanceEnergy_hh

// Unit Headers
#include <core/scoring/custom_pair_distance/FullatomCustomPairDistanceEnergy.fwd.hh>

// Package headers

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/kinematics/MinimizerMapBase.fwd.hh>

#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/RestypeDestructionEvent.fwd.hh>

// Utility headers
#include <numeric/interpolation/Histogram.hh>

// C++ headers
#include <map>
#include <list>

#include <utility/OrderedTuple.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace custom_pair_distance {

struct atoms_and_func_struct
{
	Size resA_atom_index_;
	Size resB_atom_index_;
	DistanceFuncCOP func_;
};

struct resatom_and_func_struct
{
	Size res_index_;
	Size atom_index_;
	DistanceFuncCOP func_;
};


class AtomPairFuncList : public utility::pointer::ReferenceCount
{
public:
	AtomPairFuncList();
	virtual ~AtomPairFuncList();

	std::list< atoms_and_func_struct > const & ats_n_func_list() const {
		return interacting_atompair_list_;
	}

	void
	add_interaction( atoms_and_func_struct const & interacting_pair );

private:
	std::list< atoms_and_func_struct > interacting_atompair_list_;
};

typedef utility::pointer::shared_ptr< AtomPairFuncList > AtomPairFuncListOP;
typedef utility::pointer::shared_ptr< AtomPairFuncList const > AtomPairFuncListCOP;


class PairFuncMap {
public:

	PairFuncMap() = default;

	// If you copy/assign this object, you need to detach/reattach
	// the ResidueType destruction observers for the appropriate object
	PairFuncMap( PairFuncMap const & ) = delete;
	PairFuncMap & operator=( PairFuncMap const & ) = delete;

	~PairFuncMap();

	core::Size
	size() const { return pair_map_.size(); }

	// Ideally this should pass back a COP, but the setup function needs it to be modifiable
	AtomPairFuncListOP
	get_atom_pair_func_list( core::chemical::ResidueType const & rsdtype1, core::chemical::ResidueType const & rsdtype2 ) const;

	AtomPairFuncListCOP
	get_atom_pair_func_list( core::conformation::Residue const & rsd1, core::conformation::Residue const & rsd2 ) const;

	void
	set_atom_pair_func( core::chemical::ResidueType const & rsdtype1,
		core::chemical::ResidueType const & rsdtype2,
		AtomPairFuncListOP func_list );

	void restype_destruction_observer( core::chemical::RestypeDestructionEvent const & event );

private:

	// The raw pointers here are registered in the destruction observer for their respective ResidueType
	typedef std::pair< chemical::ResidueType const *, chemical::ResidueType const * > ResTypePair;

	std::map< ResTypePair, AtomPairFuncListOP > pair_map_;

};

typedef utility::pointer::shared_ptr< PairFuncMap > PairFuncMapOP;
typedef utility::pointer::shared_ptr< PairFuncMap const > PairFuncMapCOP;

class FullatomCustomPairDistanceEnergy : public methods::ContextIndependentTwoBodyEnergy {
public:
	typedef methods::ContextIndependentTwoBodyEnergy parent;

public:
	typedef utility::fixedsizearray1< Size, 2 > ResAtomIndex;
	typedef utility::OrderedTuple< ResAtomIndex > ResAtomIndexTuple;
	typedef std::map< ResAtomIndexTuple, std::list< resatom_and_func_struct > > ResAtomIndexFuncMap;



	FullatomCustomPairDistanceEnergy();

	FullatomCustomPairDistanceEnergy( FullatomCustomPairDistanceEnergy const & src );

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;


	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	// necessary pure virtual functions
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const {}

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const {}

	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	virtual
	Distance
	atomic_interaction_cutoff() const;

	/*virtual
	void
	setup_for_minimizing(
	pose::Pose & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map
	) const;*/

	/*virtual
	void
	eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const & sfxn,
	EnergyMap const & emap,
	Vector & F1,
	Vector & F2
	) const;*/

	/// @brief Returns true if there are atoms on residue 1 & 2 that have interactions defined between them,
	/// and false otherwise.  This avoids unnecessary calls to residue_pair_energy_ext during minimization
	/// between pairs of atoms that have no interactons
	virtual
	bool
	defines_score_for_residue_pair(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		bool res_moving_wrt_eachother
	) const;

	virtual
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const;

	/// @brief Returns true.  This class takes advantage of the opportunity to store sequence-specific data
	/// in a ResiduePairMinimiazationObject listing the atom-pair interactions for a particular residue pair
	/// and to extract that data from cache for easy use inside residue_pair_energy_ext.
	virtual
	bool
	use_extended_residue_pair_energy_interface() const;

	/// @brief Extract the cached std::list< atoms_and_func_struct > object from the ResPairMinimizationData
	/// object that was stored for this particular residue pair in setup_for_minimizing_for_residue_pair, and
	/// use this cached list to evaluate the atom-pair interactions.
	virtual
	void
	residue_pair_energy_ext(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	/// @brief Find the list of atom-pair interactions from the pair_and_func_map_ for the given
	/// residue-type pair and cache that list in the input ResPairMinimizationData object.
	virtual
	void
	setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & minmap,
		ResSingleMinimizationData const & res1_data_cache,
		ResSingleMinimizationData const & res2_data_cache,
		ResPairMinimizationData & data_cache
	) const;


	/// @brief Evaluate all f1/f2 derivative vector pairs for all interacting atoms
	/// on the input residue pair.
	virtual
	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData const & min_data,
		pose::Pose const & pose, // provides context
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

private:

	bool
	interaction_defined(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2 ) const;

	void
	set_pair_and_func_map();

	// DATA

private:

	PairFuncMapOP pair_and_func_map_;
	Real max_dis_;
	virtual
	core::Size version() const;
};


// stolen from Spencer's centroid disulfide stuff
class DistanceFunc : public func::Func
{
public:
	DistanceFunc( std::string const & name );
	virtual ~DistanceFunc();
	func::FuncOP clone() const;
	virtual bool operator == ( func::Func const & rhs ) const;
	virtual bool same_type_as_me( func::Func const & other ) const;
	virtual Real func( Real const ) const;
	virtual Real dfunc( Real const ) const;
	virtual Real max_dis() const;
	virtual Real min_dis() const;
private:
	numeric::interpolation::HistogramCOP<Real,Real>::Type scores_hist_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DistanceFunc();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


/*class CacheableAtomPairFuncMap : public basic::datacache::CacheableData
{
public:
typedef utility::fixedsizearray1< Size, 2 > ResAtomIndex;
typedef utility::OrderedTuple< ResAtomIndex > ResAtomIndexTuple;
typedef std::map< ResAtomIndexTuple, std::list< resatom_and_func_struct > > ResAtomIndexFuncMap;

CacheableAtomPairFuncMap(){};
virtual ~CacheableAtomPairFuncMap(){};
virtual basic::datacache::CacheableDataOP clone() const { return new CacheableAtomPairFuncMap( *this ); }
ResAtomIndexFuncMap & map() { return map_; }
ResAtomIndexFuncMap const & map() const { return map_; }
private:
ResAtomIndexFuncMap map_;
virtual
core::Size version() const;
};*/


}
}
}


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_custom_pair_distance_FullatomCustomPairDistanceEnergy )
#endif // SERIALIZATION


#endif
