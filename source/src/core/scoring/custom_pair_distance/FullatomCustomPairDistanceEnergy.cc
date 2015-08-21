// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/custom_pair_distance/FullatomCustomPairDistanceEnergy.cc
/// @brief
/// @author David E Kim


// Unit headers
#include <core/scoring/custom_pair_distance/FullatomCustomPairDistanceEnergy.hh>
#include <core/scoring/custom_pair_distance/FullatomCustomPairDistanceEnergyCreator.hh>


// Package headers
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/MinimizationData.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <basic/database/open.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

//Auto Headers
#include <core/scoring/EnergyMap.hh>
#include <utility/OrderedTuple.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/vector1.hh>
#include <basic/datacache/CacheableData.hh>

// Numeric headers
//#include <numeric/deriv/distance_deriv.hh>

// option key includes

namespace core {
namespace scoring {
namespace custom_pair_distance {


/// @details This must return a fresh instance of the FullatomCustomPairDistanceEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
FullatomCustomPairDistanceEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new FullatomCustomPairDistanceEnergy );
}

ScoreTypes
FullatomCustomPairDistanceEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_cust_pair_dist );
	return sts;
}

////////////////////////////////////////////////////////////////////////////////////////////


class RespairInteractions : public basic::datacache::CacheableData
{
public:
	typedef utility::fixedsizearray1< Size, 2 > ResAtomIndex;
	typedef utility::OrderedTuple< ResAtomIndex > ResAtomIndexTuple;
	typedef std::map< ResAtomIndexTuple, std::list< resatom_and_func_struct > > ResAtomIndexFuncMap;

	RespairInteractions(){}
	virtual ~RespairInteractions(){}
	virtual basic::datacache::CacheableDataOP clone() const { return basic::datacache::CacheableDataOP( new RespairInteractions( *this ) ); }

	void apfc_list( AtomPairFuncListCOP apfclist ) { apfc_list_ = apfclist; }
	AtomPairFuncList const & apfc_list() const { return *apfc_list_; }
private:
	AtomPairFuncListCOP apfc_list_;
};

typedef utility::pointer::shared_ptr< RespairInteractions > RespairInteractionsOP;
typedef utility::pointer::shared_ptr< RespairInteractions const > RespairInteractionsCOP;

////////////////////////////////////////////////////////////////////////////////////////////

AtomPairFuncList::AtomPairFuncList() {}
AtomPairFuncList::~AtomPairFuncList() {}


void
AtomPairFuncList::add_interaction( atoms_and_func_struct const & interacting_pair )
{
	interacting_atompair_list_.push_back( interacting_pair );
}

////////////////////////////////////////////////////////////////////////////////////////////


static thread_local basic::Tracer tr( "core.scoring.custom_pair_distance.FullatomCustomPairDistanceEnergy" );

// constructor
FullatomCustomPairDistanceEnergy::FullatomCustomPairDistanceEnergy() :
	parent( methods::EnergyMethodCreatorOP( new FullatomCustomPairDistanceEnergyCreator ) )
{
	set_pair_and_func_map();
}

FullatomCustomPairDistanceEnergy::FullatomCustomPairDistanceEnergy( FullatomCustomPairDistanceEnergy const & src ):
	parent( src ),
	pair_and_func_map_( src.pair_and_func_map_ ),
	max_dis_( src.max_dis_ )
{
}

/// clone
methods::EnergyMethodOP
FullatomCustomPairDistanceEnergy::clone() const
{
	return methods::EnergyMethodOP( new FullatomCustomPairDistanceEnergy( *this ) );
}


void
FullatomCustomPairDistanceEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & /*pose*/,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	//ResTypePair respair;
	//respair[1] = rsd1.type();
	//respair[2] = rsd2.type();

	PairFuncMap::const_iterator respairiter = find( rsd1, rsd2 );
	if ( respairiter == pair_and_func_map_.end() ) return;
	debug_assert( (*respairiter).second ); // there should be no null-pointing elements in the pair_and_func_map_
	Energy score = 0.0;
	for ( std::list<atoms_and_func_struct>::const_iterator
			atom_func_iter     = (*respairiter).second->ats_n_func_list().begin(),
			atom_func_iter_end = (*respairiter).second->ats_n_func_list().end();
			atom_func_iter != atom_func_iter_end; ++atom_func_iter ) {
		Vector const& atom_a_xyz( rsd1.xyz((*atom_func_iter).resA_atom_index_));
		Vector const& atom_b_xyz( rsd2.xyz((*atom_func_iter).resB_atom_index_));
		score += (*atom_func_iter).func_->func(atom_a_xyz.distance_squared(atom_b_xyz));

		//if ((*atom_func_iter).func_->func(atom_a_xyz.distance_squared(atom_b_xyz)) > 0)
		//tr.Debug << rsd1.name() << " " << rsd1.seqpos() << " " << rsd1.atom_name((*atom_func_iter).resA_atom_index_) << " " <<
		// rsd2.name() << " " << " " << rsd2.seqpos() << " " << rsd2.atom_name((*atom_func_iter).resB_atom_index_) <<
		// " score: " << (*atom_func_iter).func_->func(atom_a_xyz.distance_squared(atom_b_xyz)) << " dist_sq: " << atom_a_xyz.distance_squared(atom_b_xyz) << std::endl;
	}
	emap[ fa_cust_pair_dist ] += score;
}

/// @brief FullatomCustomPairDistanceEnergy distance cutoff
Distance
FullatomCustomPairDistanceEnergy::atomic_interaction_cutoff() const
{
	return max_dis_;
}

bool
FullatomCustomPairDistanceEnergy::defines_score_for_residue_pair(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	bool res_moving_wrt_eachother
) const
{
	return  res_moving_wrt_eachother && interaction_defined( res1, res2 );
}

bool
FullatomCustomPairDistanceEnergy::minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

bool
FullatomCustomPairDistanceEnergy::use_extended_residue_pair_energy_interface() const
{
	return true;
}

void
FullatomCustomPairDistanceEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData const & pair_data,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	debug_assert( rsd1.seqpos() < rsd2.seqpos() );
	debug_assert( utility::pointer::dynamic_pointer_cast< RespairInteractions const > (pair_data.get_data( fa_custom_pair_dist_data ) ));

	RespairInteractions const & respair_intxns( static_cast< RespairInteractions const & > (pair_data.get_data_ref( fa_custom_pair_dist_data ) ));
	Energy score( 0.0 );
	for ( std::list< atoms_and_func_struct >::const_iterator
			atom_func_iter     = respair_intxns.apfc_list().ats_n_func_list().begin(),
			atom_func_iter_end = respair_intxns.apfc_list().ats_n_func_list().end();
			atom_func_iter != atom_func_iter_end; ++atom_func_iter ) {
		Vector const& atom_a_xyz( rsd1.xyz((*atom_func_iter).resA_atom_index_));
		Vector const& atom_b_xyz( rsd2.xyz((*atom_func_iter).resB_atom_index_));
		score += (*atom_func_iter).func_->func(atom_a_xyz.distance_squared(atom_b_xyz));
	}
	emap[ fa_cust_pair_dist ] += score;
}

void
FullatomCustomPairDistanceEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData & pair_data
) const
{
	debug_assert( rsd1.seqpos() < rsd2.seqpos() );

	// update the existing respair_intxns object if one is already present in the pair_data object
	RespairInteractionsOP respair_intxns = utility::pointer::static_pointer_cast< core::scoring::custom_pair_distance::RespairInteractions > ( pair_data.get_data( fa_custom_pair_dist_data ) );
	if ( ! respair_intxns ) respair_intxns = RespairInteractionsOP( new RespairInteractions );

	AtomPairFuncListCOP apfclist = find( rsd1, rsd2 )->second;
	debug_assert( apfclist );
	respair_intxns->apfc_list( apfclist );
	pair_data.set_data( fa_custom_pair_dist_data, respair_intxns );
}


void
FullatomCustomPairDistanceEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const & pair_data,
	pose::Pose const &,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	debug_assert( rsd1.seqpos() < rsd2.seqpos() );
	debug_assert( rsd1.seqpos() < rsd2.seqpos() );
	debug_assert( utility::pointer::dynamic_pointer_cast< RespairInteractions const > (pair_data.get_data( fa_custom_pair_dist_data ) ));

	RespairInteractions const & respair_intxns( static_cast< RespairInteractions const & > (pair_data.get_data_ref( fa_custom_pair_dist_data ) ));
	for ( std::list< atoms_and_func_struct >::const_iterator
			atom_func_iter     = respair_intxns.apfc_list().ats_n_func_list().begin(),
			atom_func_iter_end = respair_intxns.apfc_list().ats_n_func_list().end();
			atom_func_iter != atom_func_iter_end; ++atom_func_iter ) {
		Vector const & atom_a_xyz = rsd1.xyz( atom_func_iter->resA_atom_index_ );
		Vector const & atom_b_xyz = rsd2.xyz( atom_func_iter->resB_atom_index_ );

		Real dist_sq = atom_a_xyz.distance_squared(atom_b_xyz);
		if ( dist_sq < atom_func_iter->func_->min_dis() || dist_sq > atom_func_iter->func_->max_dis() ) continue;

		Vector f1( atom_a_xyz.cross( atom_b_xyz ));
		Vector f2( atom_a_xyz - atom_b_xyz );
		Real const dist( f2.length() );
		if ( dist == 0.0 ) continue;

		Real deriv = (*atom_func_iter).func_->dfunc(dist_sq);
		// WRONG f1 *= ( deriv / dist ) * weights[ fa_cust_pair_dist ];
		// WRONG f2 *= ( deriv / dist ) * weights[ fa_cust_pair_dist ];
		f1 *= deriv * 2 * weights[ fa_cust_pair_dist ];
		f2 *= deriv * 2 * weights[ fa_cust_pair_dist ];

		/*  THIS CODE WORKS!
		Vector f1,f2; Real dis;
		numeric::deriv::distance_f1_f2_deriv( atom_a_xyz, atom_b_xyz, dis, f1, f2 );
		Real const dE_dd2 = (*atom_func_iter).func_->dfunc(dis*dis);
		f1 *= dE_dd2 * 2 * dis * weights[ fa_cust_pair_dist ];
		f2 *= dE_dd2 * 2 * dis * weights[ fa_cust_pair_dist ];*/

		// equal and opposite
		r1_atom_derivs[ atom_func_iter->resA_atom_index_ ].f1() += f1;
		r1_atom_derivs[ atom_func_iter->resA_atom_index_ ].f2() += f2;
		r2_atom_derivs[ atom_func_iter->resB_atom_index_ ].f1() -= f1;
		r2_atom_derivs[ atom_func_iter->resB_atom_index_ ].f2() -= f2;
	}

}


/*void
FullatomCustomPairDistanceEnergy::setup_for_minimizing(
pose::Pose & pose,
ScoreFunction const &,
kinematics::MinimizerMapBase const &
) const
{
using namespace basic::datacache;

if( !pose.data().has( pose::datacache::CacheableDataType::CUSTOM_PAIR_DIST_SCORE_INFO ) ) {
pose.data().set( pose::datacache::CacheableDataType::CUSTOM_PAIR_DIST_SCORE_INFO, new CacheableAtomPairFuncMap );
}
CacheableDataOP dat( pose.data().get_ptr( pose::datacache::CacheableDataType::CUSTOM_PAIR_DIST_SCORE_INFO ) );
CacheableAtomPairFuncMap *cachemap = (CacheableAtomPairFuncMap*)dat();
ResAtomIndexFuncMap & resatomfuncmap(cachemap->map());

// Identify which residue and atom pairs to evaluate in the pose by looping
// through the pair_and_func_map keys (residue pairs to evaluate).
// Save the data in the pose data cache for evaluating derivative.
for (PairFuncMap::const_iterator respairiter = pair_and_func_map_.begin(),
respairiter_end = pair_and_func_map_.end();
respairiter != respairiter_end; ++respairiter) {
ResTypePair respair = (*respairiter).first.data();
for ( Size i=1; i<=pose.total_residue(); ++i ) {
if ( & (pose.residue_type(i)) == respair[1] ) {
for ( Size j=1; j<=pose.total_residue(); ++j ) {
if ( & (pose.residue_type(j)) == respair[2] ) {
for ( std::list<atoms_and_func_struct>::const_iterator iter_a = (*respairiter).second.begin(),
iter_end_a = (*respairiter).second.end(); iter_a != iter_end_a; ++iter_a) {
if ( i == j && (*iter_a).resA_atom_index_ == (*iter_a).resB_atom_index_ ) continue; // skip same atom
// add pair
resatom_and_func_struct resatom_b;
ResAtomIndex rai_a;
rai_a[1] = i;
rai_a[2] = (*iter_a).resA_atom_index_;
resatom_b.res_index_ = j;
resatom_b.atom_index_ = (*iter_a).resB_atom_index_;
resatom_b.func_ = (*iter_a).func_;
resatomfuncmap[rai_a].push_back( resatom_b );
}
}
}
}
}
}
}


void
FullatomCustomPairDistanceEnergy::eval_atom_derivative(
id::AtomID const & id,
pose::Pose const & pose,
kinematics::DomainMap const &, //domain_map,
ScoreFunction const &, // sfxn,
EnergyMap const & emap,
Vector & F1,
Vector & F2
) const
{
using namespace basic::datacache;
ResAtomIndex resatom;
resatom[1] = id.rsd();
resatom[2] = id.atomno();
CacheableDataCOP dat( pose.data().get_const_ptr( pose::datacache::CacheableDataType::CUSTOM_PAIR_DIST_SCORE_INFO ) );
CacheableAtomPairFuncMap const *cachemap = (CacheableAtomPairFuncMap const *)dat();
ResAtomIndexFuncMap const & resatomfuncmap(cachemap->map());
ResAtomIndexFuncMap::const_iterator resatomiter = resatomfuncmap.find( resatom );
if ( resatomiter == resatomfuncmap.end() ) return;

for ( std::list<resatom_and_func_struct>::const_iterator iter_a = (*resatomiter).second.begin(),
iter_end_a = (*resatomiter).second.end(); iter_a != iter_end_a; ++iter_a) {
// determine distance squared
Vector const& atom_a_xyz( pose.residue(id.rsd()).atom(id.atomno()).xyz() );
Vector const& atom_b_xyz( pose.residue((*iter_a).res_index_).atom((*iter_a).atom_index_).xyz() );
Real dist_sq = atom_a_xyz.distance_squared(atom_b_xyz);
if (dist_sq < (*iter_a).func_->min_dis() || dist_sq > (*iter_a).func_->max_dis()) continue;
Vector const f1( atom_a_xyz.cross( atom_b_xyz ));
Vector const f2( atom_a_xyz - atom_b_xyz );
Real const dist( f2.length() );
debug_assert(dist != 0);
Real deriv = (*iter_a).func_->dfunc(dist_sq);
F1 += ( deriv / dist ) * emap[ fa_cust_pair_dist ] * f1;
F2 += ( deriv / dist ) * emap[ fa_cust_pair_dist ] * f2;

}
}*/

FullatomCustomPairDistanceEnergy::PairFuncMap::const_iterator
FullatomCustomPairDistanceEnergy::find(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const
{
	ResTypePair respair;
	respair[1] = & rsd1.type();
	respair[2] = & rsd2.type();
	return pair_and_func_map_.find( respair );
}


bool
FullatomCustomPairDistanceEnergy::interaction_defined(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const
{
	return find( rsd1, rsd2 ) != pair_and_func_map_.end();
}

void
FullatomCustomPairDistanceEnergy::set_pair_and_func_map()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace chemical;

	std::string pairfuncfile = (option[score::fa_custom_pair_distance_file].user()) ?
		option[score::fa_custom_pair_distance_file]() :
		basic::database::full_name( "scoring/score_functions/custom_pair_distance/fa_custom_pair_distance" );

	tr.Debug << "Reading fa_custom_pair_distance_file: " << pairfuncfile << std::endl;

	utility::io::izstream in( pairfuncfile );
	if ( !in.good() ) {
		utility_exit_with_message( "Unable to open fa_custom_pair_distance file: " + pairfuncfile );
	}

	std::string line, residue_type_set, score_function_name;
	utility::vector1< std::string > resA, resB, atomA, atomB;
	max_dis_ = 0.0;

	while ( getline( in, line ) ) {
		if ( line.size() < 1 || line[0] == '#' ) continue;
		std::istringstream l( line );
		std::string tag;
		l >> tag;
		if ( tag == "RESIDUE_TYPE_SET" ) {
			resA.clear(); resB.clear(); atomA.clear(); atomB.clear();
			l >> residue_type_set;
			tr.Debug << "Residue type set: " <<  residue_type_set << std::endl;
			continue;
		}
		if ( tag == "PAIR" ) {
			std::string buf;
			utility::vector1< std::string > pair_cols;
			while ( l >> buf ) {
				pair_cols.push_back( buf );
			}
			if ( pair_cols.size() == 4 ) {
				resA.push_back( pair_cols[1] );
				atomA.push_back( pair_cols[2] );
				resB.push_back( pair_cols[3] );
				atomB.push_back( pair_cols[4] );
				tr.Debug << "Pair: " << pair_cols[1] << " " <<  pair_cols[2] << " " << pair_cols[3] << " " << pair_cols[4] << std::endl;
				continue;
			}
		}
		if ( tag == "SCORE_FUNCTION" ) {
			l >> score_function_name;

			atoms_and_func_struct pair_func;
			pair_func.func_ = DistanceFuncCOP( DistanceFuncOP( new DistanceFunc( score_function_name ) ) );
			if ( pair_func.func_->max_dis() > max_dis_ ) {
				max_dis_ = pair_func.func_->max_dis();
			}
			tr.Debug << "SCORE_FUNCTION: " << score_function_name << " (min: " <<
				pair_func.func_->min_dis() << " max: " << pair_func.func_->max_dis() << ")" << std::endl;

			ResidueTypeSetCOP restype_set =
				ChemicalManager::get_instance()->residue_type_set( residue_type_set );

			// get all possible residue types for each residue pair
			for ( Size i = 1; i <= resA.size(); ++i ) {
				ResidueTypeCOPs const & possible_res_types_a = restype_set->name3_map_DO_NOT_USE( resA[i] );
				ResidueTypeCOPs const & possible_res_types_b = restype_set->name3_map_DO_NOT_USE( resB[i] );
				for ( Size j = 1; j <= possible_res_types_a.size(); ++j ) {
					ResidueTypeCOP const & rsd_type_a = possible_res_types_a[ j ];
					Size atom_index_a = rsd_type_a->atom_index( atomA[i] );
					for ( Size k = 1; k <= possible_res_types_b.size(); ++k ) {
						ResidueTypeCOP const & rsd_type_b = possible_res_types_b[ k ];
						Size atom_index_b = rsd_type_b->atom_index( atomB[i] );
						ResTypePair respair;

						respair[1] = rsd_type_a.get();
						respair[2] = rsd_type_b.get();
						pair_func.resA_atom_index_ =  atom_index_a;
						pair_func.resB_atom_index_ =  atom_index_b;
						AtomPairFuncListOP atpairlist = pair_and_func_map_[ respair ];
						if ( ! atpairlist ) {
							//std::cout << "Adding new residue-pair interaction for " << respair[1]->name() << " " << respair[2]->name() << std::endl;
							atpairlist = AtomPairFuncListOP( new AtomPairFuncList );
							pair_and_func_map_[ respair ] = atpairlist;
						} else {
							//std::cout << "Updating residue-pair interaction for " << respair[1]->name() << " " << respair[2]->name() << std::endl;
						}
						atpairlist->add_interaction( pair_func );

						// save the mirrored pair if the types are not equal so if AB or
						// BA gets passed to residue_pair_energy they'll be evaluated.
						// I tried to prevent having to do this by sorting the pair but
						// couldn't get it to work
						if ( rsd_type_a != rsd_type_b ) {
							respair[1] = rsd_type_b.get();
							respair[2] = rsd_type_a.get();
							pair_func.resA_atom_index_ =  atom_index_b;
							pair_func.resB_atom_index_ =  atom_index_a;
							atpairlist =  pair_and_func_map_[ respair ];
							if ( ! atpairlist ) {
								atpairlist = AtomPairFuncListOP( new AtomPairFuncList );
								pair_and_func_map_[ respair ] = atpairlist;
							}
							atpairlist->add_interaction( pair_func );
						}
					}
				}
			}
		}
	}
	tr << "Added " << pair_and_func_map_.size() << " AtomPairFuncList lists" << std::endl;
}

core::Size
FullatomCustomPairDistanceEnergy::version() const
{
	return 1; // Initial versioning
}
/// NOTE: distance functions should be input in square-distance bins.

DistanceFunc::DistanceFunc( std::string const name ) {
	utility::io::izstream scores_stream;
	basic::database::open( scores_stream, "scoring/score_functions/custom_pair_distance/" + name);
	scores_hist_ = numeric::interpolation::HistogramCOP<Real,Real>::Type( numeric::interpolation::HistogramOP<Real, Real>::Type( new numeric::interpolation::Histogram<Real,Real>( scores_stream() ) ) );
	scores_stream.close();
}

DistanceFunc::~DistanceFunc() {}
Real DistanceFunc::func( Real const dist_sq ) const {
	Real e(0.0);
	if ( dist_sq < scores_hist_->minimum() ||
			dist_sq > scores_hist_->maximum() ) return e;
	scores_hist_->interpolate(dist_sq,e);
	return e;
}

Real DistanceFunc::dfunc( Real const dist_sq ) const {
	Real df(0.0), e(0.0);
	if ( dist_sq < scores_hist_->minimum() ||
			dist_sq > scores_hist_->maximum() ) return df;
	scores_hist_->interpolate(dist_sq,e,df);
	return df;
}

Real DistanceFunc::max_dis() const {
	return scores_hist_->maximum();
}

Real DistanceFunc::min_dis() const {
	return scores_hist_->minimum();
}

}
}
}
