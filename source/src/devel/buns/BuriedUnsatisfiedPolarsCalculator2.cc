// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/buns/BuriedUnsatPolarsFinder2.cc
/// @brief
/// @details
/// @author Kevin Houlihan (khouli@unc.edu)
/// @author Bryan Der

// Unit headers

#include <devel/buns/BuriedUnsatisfiedPolarsCalculator2.hh>

//Project Headers
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
//#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>
#include <protocols/vardist_solaccess/VarSolDRotamerDots.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
//#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <core/conformation/Residue.hh>
#include <utility>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <numeric/conversions.hh>
#include <basic/options/keys/bunsat_calc2.OptionKeys.gen.hh>


// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <basic/MetricValue.hh>

#include <utility/assert.hh>

#include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>

static basic::Tracer TR( "devel.buns.BuriedUnsatisfiedPolarsCalculator2" );


#ifdef    SERIALIZATION
// Project serialization headers
#include <core/id/AtomID_Map.srlz.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace devel {
namespace buns {

using namespace core;
using namespace core::chemical;
using namespace core::conformation;
using namespace core::pose;
using namespace basic::options;
using namespace core::pose::metrics;
using namespace protocols::toolbox::pose_metric_calculators;
using namespace utility;
using Vector = numeric::xyzVector<core::Real>;

BuriedUnsatisfiedPolarsCalculator2::BuriedUnsatisfiedPolarsCalculator2(
	std::string const & weak_bunsat_calc
) :
	all_bur_unsat_polars_( 0 ),
	special_region_bur_unsat_polars_(0),
	name_of_weak_bunsat_calc_( weak_bunsat_calc ),
	layered_sasa_(option[OptionKeys::bunsat_calc2::layered_sasa]),
	generous_hbonds_(option[OptionKeys::bunsat_calc2::generous_hbonds]),
	sasa_burial_cutoff_(option[OptionKeys::bunsat_calc2::sasa_burial_cutoff]),
	AHD_cutoff_(option[OptionKeys::bunsat_calc2::AHD_cutoff]),
	dist_cutoff_(option[OptionKeys::bunsat_calc2::dist_cutoff]),
	hxl_dist_cutoff_(option[OptionKeys::bunsat_calc2::hxl_dist_cutoff]),
	sulph_dist_cutoff_(option[OptionKeys::bunsat_calc2::sulph_dist_cutoff]),
	metal_dist_cutoff_(option[OptionKeys::bunsat_calc2::metal_dist_cutoff])
{
	atom_bur_unsat_.clear();
	residue_bur_unsat_polars_.clear();
	special_region_.clear();
	assert_calculators();
}


BuriedUnsatisfiedPolarsCalculator2::BuriedUnsatisfiedPolarsCalculator2(
	std::string const & weak_bunsat_calc,
	std::set< core::Size > const & special_region
) :
	all_bur_unsat_polars_(0),
	special_region_bur_unsat_polars_(0),
	name_of_weak_bunsat_calc_( weak_bunsat_calc ),
	special_region_( special_region ),
	layered_sasa_(option[OptionKeys::bunsat_calc2::layered_sasa]),
	generous_hbonds_(option[OptionKeys::bunsat_calc2::generous_hbonds]),
	sasa_burial_cutoff_(option[OptionKeys::bunsat_calc2::sasa_burial_cutoff]),
	AHD_cutoff_(option[OptionKeys::bunsat_calc2::AHD_cutoff]),
	dist_cutoff_(option[OptionKeys::bunsat_calc2::dist_cutoff]),
	hxl_dist_cutoff_(option[OptionKeys::bunsat_calc2::hxl_dist_cutoff]),
	sulph_dist_cutoff_(option[OptionKeys::bunsat_calc2::sulph_dist_cutoff]),
	metal_dist_cutoff_(option[OptionKeys::bunsat_calc2::metal_dist_cutoff])
{
	atom_bur_unsat_.clear();
	residue_bur_unsat_polars_.clear();
	assert_calculators();
}


void
BuriedUnsatisfiedPolarsCalculator2::assert_calculators() {
	if ( !CalculatorFactory::Instance().check_calculator_exists( name_of_weak_bunsat_calc_ ) ) {
		// sboyken NEED DIFFERENT NAME THAN "default" OTHERWISE DOES NOT PLAY NICE WHEN CALL MULTIPLE UNSAT FILTERS WITHIN SAME XML!!!
		//  if ( name_of_weak_bunsat_calc_ != "default" ) {
		//   TR <<
		//    "Attention: couldn't find the specified buried unsat calculator ( " <<
		//    name_of_weak_bunsat_calc_ << " ), instantiating default one." << std::endl;
		//  }
		name_of_weak_bunsat_calc_ = "bunsat_calc2_default_weak_bunsat_calc";

		std::string sasa_calc_name("sasa_calc_name");
		if ( !CalculatorFactory::Instance().check_calculator_exists( sasa_calc_name ) ) {
			if ( layered_sasa_ ) {
				using namespace protocols::vardist_solaccess;
				TR << "Registering VarSolDist SASA Calculator" << std::endl;
				CalculatorFactory::Instance().register_calculator(sasa_calc_name,
					VarSolDistSasaCalculatorOP( new VarSolDistSasaCalculator() ));
			} else {
				TR << "Registering SASA Calculator" << std::endl;
				CalculatorFactory::Instance().register_calculator(sasa_calc_name,
					PoseMetricCalculatorOP( new pose::metrics::simple_calculators::SasaCalculatorLegacy() ));
			}
		}
		std::string num_hbonds_calc_name("num_hbonds_calc_name");
		if ( !CalculatorFactory::Instance().check_calculator_exists( num_hbonds_calc_name ) ) {
			CalculatorFactory::Instance().register_calculator(num_hbonds_calc_name,
				PoseMetricCalculatorOP( new NumberHBondsCalculator() ));
		}
		if ( !CalculatorFactory::Instance().check_calculator_exists( name_of_weak_bunsat_calc_ ) ) {
			using namespace protocols::toolbox::pose_metric_calculators;
			TR << "Registering new basic buried unsat calculator." << std::endl;
			// need to set weak_bunsat_calc to legacy to maintain existing BUNS2 behavior
			//   17/09/05 sboyken NEED TO MAKE SURE THAT THIS SETS DEFAULTS CORRECTLY AND MAINTAINS CURRENT BUNS2 BEHAVIOR
			CalculatorFactory::Instance().register_calculator( name_of_weak_bunsat_calc_, PoseMetricCalculatorOP(
				new BuriedUnsatisfiedPolarsCalculator(sasa_calc_name, num_hbonds_calc_name, sasa_burial_cutoff_, false /* generous */, true /* legacy counting */, layered_sasa_ /* vsasa */ )));
			// defining vsasa won't matter because calculators asserted above and weak calc will find them based on the name
		}
	}
}


void
BuriedUnsatisfiedPolarsCalculator2::lookup(
	std::string const & key,
	basic::MetricValueBase * valptr
) const
{

	if ( key == "all_bur_unsat_polars" ) {
		basic::check_cast( valptr, &all_bur_unsat_polars_,
			"all_bur_unsat_polars expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( all_bur_unsat_polars_ );

	} else if ( key == "special_region_bur_unsat_polars" ) {
		basic::check_cast( valptr, &special_region_bur_unsat_polars_,
			"special_region_bur_unsat_polars expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( special_region_bur_unsat_polars_ );

	} else if ( key == "atom_bur_unsat" ) {
		basic::check_cast( valptr, &atom_bur_unsat_,
			"atom_bur_unsat expects to return a id::AtomID_Map< bool >" );
		(static_cast<basic::MetricValue<id::AtomID_Map< bool > > *>(valptr))->set( atom_bur_unsat_ );

	} else if ( key == "residue_bur_unsat_polars" ) {
		basic::check_cast( valptr, &residue_bur_unsat_polars_,
			"residue_bur_unsat_polars expects to return a utility::vector1< Size >" );
		(static_cast<basic::MetricValue<utility::vector1< Size > > *>(valptr))->set( residue_bur_unsat_polars_ );

	} else {
		basic::Error() << "BuriedUnsatisfiedPolarsCalculator2 cannot compute the requested metric " << key << std::endl;
		utility_exit();
	}

} //lookup


std::string
BuriedUnsatisfiedPolarsCalculator2::print( std::string const & key ) const
{
	if ( key == "all_bur_unsat_polars" ) {
		return utility::to_string( all_bur_unsat_polars_ );
	} else if ( key == "special_region_bur_unsat_polars" ) {
		return utility::to_string( special_region_bur_unsat_polars_ );
	} else if ( key == "residue_bur_unsat_polars" ) {
		return utility::to_string( residue_bur_unsat_polars_ );
	}

	basic::Error() << "BuriedUnsatisfiedPolarsCalculator2 cannot compute metric " << key << std::endl;
	utility_exit();
	return "";

} //print


void
BuriedUnsatisfiedPolarsCalculator2::recompute( Pose const & pose )
{
	TR << "Recomputing buried unsats" << std::endl;
	this->show();

	all_bur_unsat_polars_ = 0;
	special_region_bur_unsat_polars_ = 0;

	if ( pose.size() != residue_bur_unsat_polars_.size() ) {
		residue_bur_unsat_polars_.resize( pose.size() );
		atom_bur_unsat_.resize( pose.size() );
	}

	if ( generous_hbonds_ ) {
		generous_hbond();
	}

	TR << "Running basic buried unsat calc" << std::endl;
	basic::MetricValue< id::AtomID_Map< bool > > bunsat_atomid_map;
	pose.metric(name_of_weak_bunsat_calc_, "atom_bur_unsat", bunsat_atomid_map);

	//id::AtomID_Map< bool > bunsat_thorough_atomid_map(bunsat_atomid_map.value());
	//bunsats_thorough_check(pose, bunsat_thorough_atomid_map);
	atom_bur_unsat_ = bunsat_atomid_map.value();

	TR << "Validating buried unsats" << std::endl;

	bunsats_thorough_check(pose, atom_bur_unsat_);


	//    pose.update_residue_neighbors();
	//    scoring::ScoreFunctionOP scorefxn = scoring::get_score_function();
	//    scorefxn->score(pose);
	//    pose.energies();
	//
	//    const Real max_nbr_radius = pose::pose_max_nbr_radius(pose);
	//    //constexpr const Real max_covalent_hbond_length = 1.5;
	//    const Real max_covalent_hbond_length = 1.5;
	//    const Real neighbor_cutoff = 2 * max_nbr_radius + max_covalent_hbond_length + AHdist_threshold;
	//
	//    conformation::PointGraphOP pg = new conformation::PointGraph;
	//    conformation::residue_point_graph_from_conformation(pose.conformation(), *pg);
	//    conformation::find_neighbors<conformation::PointGraphVertexData,
	//        conformation::PointGraphEdgeData>(pg, neighbor_cutoff);
	//
	//    id::AtomID_Map<Real> atom_sasa = vsasa_calc.calculate(pose);
	//
	//    Size nres = pose.size();
	//
	//    scoring::hbonds::HBondSet hbond_set;
	//
	//    const scoring::hbonds::HBondDatabaseCOP hb_database = scoring::hbonds::HBondDatabase::get_database();
	//    const scoring::TenANeighborGraph& tenA_neighbor_graph(pose.energies().tenA_neighbor_graph());
	//
	//    for (Size lowerResNum = 1; lowerResNum <= nres; ++lowerResNum ) {
	//        const conformation::Residue& lowerRes(pose.residue(lowerResNum));
	//        const Size nbl = tenA_neighbor_graph.get_node(lowerResNum)->
	//            num_neighbors_counting_self_static();
	//        for (conformation::PointGraph::UpperEdgeListConstIter
	//                ue  = pg->get_vertex(lowerResNum).const_upper_edge_list_begin(),
	//                ue_end = pg->get_vertex(lowerResNum).const_upper_edge_list_end();
	//                ue != ue_end; ++ue ) {
	//            Size upperResNum = ue->upper_vertex();
	//            const conformation::Residue& upperRes(pose.residue(upperResNum));
	//            const Size nbu = tenA_neighbor_graph.get_node(upperResNum)->
	//                num_neighbors_counting_self_static();
	//            scoring::hbonds::identify_hbonds_1way_AHdist(*hb_database,
	//                lowerRes, upperRes, nbl, nbu, AHdist_threshold, hbond_set);
	//            scoring::hbonds::identify_hbonds_1way_AHdist(*hb_database,
	//                upperRes, lowerRes, nbu, nbl, AHdist_threshold, hbond_set);
	//        }
	//    }
	//
	//    Size buns = 0;
	//    const pose::PDBInfo& pdb_info = *(pose.pdb_info());
	//    for (Size resNum = 1; resNum <= nres; ++resNum) {
	//        const conformation::Residue& res = pose.residue(resNum);
	//        const chemical::ResidueType& res_type = res.type();
	//        for (Size atomNum = 1, natom = pose.residue(resNum).natoms(); atomNum <= natom; ++atomNum){
	//            if (!res.heavyatom_is_an_acceptor(atomNum) && !res.atom_is_polar_hydrogen(atomNum)) continue;
	//            if (pdb_info.temperature(resNum, atomNum) > 30) continue;
	//            if (pdb_info.occupancy(resNum, atomNum) < 1) continue;
	//
	//            id::AtomID at(atomNum, resNum);
	//            const Real vsasa = atom_sasa[at];
	//
	//            const core::chemical::AtomType& atom_type = res_type.atom_type(atomNum);
	//
	//            if (vsasa > burial_cutoff) continue;
	//            utility::vector1<scoring::hbonds::HBondCOP> hbonds = hbond_set.atom_hbonds(at, false /*include only allowed*/);
	//            bool hbonded = false;
	//            for (utility::vector1<scoring::hbonds::HBondCOP>::iterator h = hbonds.begin(), end = hbonds.end();
	//                 h != end; h++) {
	//                if (hb_eval.evaluate(pose, *h)){
	//                    hbonded = true;
	//                    break;
	//                }
	//            }
	//            if (!hbonded) {
	//                ++buns;
	//    // found a bunsat, do stuff here
	//            }
	//        }
	//    }


} //recompute

void
BuriedUnsatisfiedPolarsCalculator2::generous_hbond()
const
{
	TR << "Setting generous_hbonds." << std::endl;
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	scoring::methods::EnergyMethodOptionsOP emopts( new scoring::methods::EnergyMethodOptions(
		scorefxn->energy_method_options() ) );
	emopts->hbond_options().decompose_bb_hb_into_pair_energies(true);
	emopts->hbond_options().use_hb_env_dep(false);
	emopts->hbond_options().bb_donor_acceptor_check(false);
	scorefxn->set_energy_method_options(*emopts);
}

void
BuriedUnsatisfiedPolarsCalculator2::bunsats_thorough_check(
	pose::Pose const & pose,
	id::AtomID_Map< bool > & bunsat_thorough_atomid_map
)
{
	for ( Size res = 1; res <= bunsat_thorough_atomid_map.size(); ++res ) {

		residue_bur_unsat_polars_[res] = 0;

		for ( Size atm = 1; atm <= bunsat_thorough_atomid_map.n_atom(res); ++atm ) {

			//std::string res_debug = pose.pdb_info()->pose2pdb(res);
			//std::string atm_debug = pose.residue(res).atom_name(static_cast< int >(atm));

			// if unsat by weak bunsat calc, do extra checks
			if ( bunsat_thorough_atomid_map(res, atm) ) {
				TR << "Considering candidate buried unsat, residue " << res << ", atom " << atm << "." << std::endl;
				id::AtomID bunsat_candidate_atom_id(atm, res);
				bunsat_thorough_atomid_map.set(bunsat_candidate_atom_id,
					single_bunsat_thorough_check(pose, bunsat_candidate_atom_id));
				//if (!bunsat_thorough_atomid_map(res, atm)) TR << "Rejected a buried unsat, residue " << res_debug << ", atom " << atm_debug << "." << std::endl;
			}
			// if still bunsat, add to counts
			if ( bunsat_thorough_atomid_map(res, atm) ) {
				TR << "Validated a buried unsat, rosetta resi " << res << ", atom " << atm << "." << std::endl;
				residue_bur_unsat_polars_[res]++;
				if ( special_region_.find( res ) != special_region_.end() ) {
					special_region_bur_unsat_polars_++;
				}
			}
		}
	}
}

/*
return true if bunsat_candidate_atom_id in pose passes additional checks,
otherwise return false
*/
bool
BuriedUnsatisfiedPolarsCalculator2::single_bunsat_thorough_check(
	pose::Pose const & pose,
	id::AtomID const & bunsat_candidate_atom_id
)
{
	Size const bunsat_resi = bunsat_candidate_atom_id.rsd();

	TR << "Calculating neighbors for resi " << bunsat_resi << std::endl;

	NeighborsByDistanceCalculator nbr_calc(bunsat_resi);
	basic::MetricValue< std::set< Size > > neighbors;
	nbr_calc.get("neighbors", neighbors, pose);
	//pose.metric(nbr_calc_name, "neighbors", neighbors);

	Residue const & bunsat_rsd(pose.residue(bunsat_resi));
	Size const bunsat_atom_num = bunsat_candidate_atom_id.atomno();

	Vector bunsat_xyz = bunsat_rsd.atom(bunsat_atom_num).xyz();

	bool bunsat_is_donor(bunsat_rsd.atom_type(
		static_cast< int >(bunsat_atom_num)).is_donor());

	bool bunsat_is_acceptor(bunsat_rsd.atom_type(
		static_cast< int >(bunsat_atom_num)).is_acceptor());

	Size num_hbonds = 0;

	if ( bunsat_is_donor ) {
		for ( core::Size test_resi : neighbors.value() ) {
			//std::string test_res_debug = pose.pdb_info()->pose2pdb(test_resi);
			//TR << "checking if donates to " << test_resi << std::endl;

			bunsat_donor_nbr_residue_check(pose, bunsat_candidate_atom_id, bunsat_rsd, bunsat_xyz, test_resi, num_hbonds);
		}
	}
	if ( bunsat_is_acceptor ) {
		for ( core::Size test_resi : neighbors.value() ) {
			bunsat_acc_nbr_residue_check(pose, bunsat_candidate_atom_id, bunsat_rsd, bunsat_xyz, test_resi, num_hbonds);
		}
	}

	Size satisfac_cut = satisfaction_cutoff( bunsat_rsd.type().atom_type( bunsat_atom_num ).name() );
	Size bonded_heavyatoms = bunsat_rsd.n_bonded_neighbor_all_res( bunsat_atom_num ) - bunsat_rsd.type().number_bonded_hydrogens( bunsat_atom_num );
	//all_bur_unsat_hbonds_ += (satisfac_cut - ( bonded_heavyatoms + num_hbonds));
	Size nbonds = bonded_heavyatoms + num_hbonds;
	TR << "Residue " << bunsat_resi << ", atom " << bunsat_atom_num << " has " << nbonds << " bonded heavy atoms." << std::endl;
	if ( nbonds < satisfac_cut ) {
		TR << "nbonds less than satisfac_cut of" << satisfac_cut << std::endl;
		all_bur_unsat_polars_++;
		return true;
	} else {
		return false;
	}
}

/*
return false if a plausible hbond between bunsat_candidate_atom_id as a donor
and test_resi as an acceptor is found, otherwise return true
Why don't we have it return the number of hbonds instead of just true.
Then we can add up the sum of hbonds detected and run the same satisfac_cut check as in
the original BuriedUnsatisfiedHbondCalculator. Here we could be say the atom is fine if
it only makes 1 hbond when it may need to make 2.
*/
void
BuriedUnsatisfiedPolarsCalculator2::bunsat_donor_nbr_residue_check(
	pose::Pose const & pose,
	id::AtomID const & bunsat_candidate_atom_id,
	Residue const & bunsat_rsd,
	Vector const & bunsat_xyz,
	Size const test_resi,
	Size & num_hbonds
)
{
	Size const bunsat_resi = bunsat_candidate_atom_id.rsd();
	Size const bunsat_atom_num = bunsat_candidate_atom_id.atomno();
	std::string const bunsat_atom_name =
		bunsat_rsd.atom_name(static_cast< int >(bunsat_atom_num));

	Residue const & test_rsd(pose.residue(test_resi));

	AtomIndices const & test_accpt_pos = test_rsd.accpt_pos();
	for ( auto test_atom_it = test_accpt_pos.begin();
			test_atom_it < test_accpt_pos.end(); ++test_atom_it ) {
		Size const test_atom_num = *test_atom_it;

		// exclude neighboring backbone-backbone
		std::string const test_atom_name(
			test_rsd.atom_name(static_cast< int >(test_atom_num)));
		if ( adjacent_bbbb_check(bunsat_resi, bunsat_atom_name, test_resi, test_atom_name) ) {
			continue;
		}

		// exclude self sidechain-sidechain
		if ( self_scsc(bunsat_rsd, bunsat_resi, bunsat_atom_num, test_rsd, test_resi, test_atom_num) ) {
			continue;
		}

		Vector test_xyz = test_rsd.atom(test_atom_num).xyz();
		// check for sulfur-bonds
		if ( sulphur_bond_check(test_rsd, test_atom_num, bunsat_xyz, test_xyz) ) {
			TR << "HBond donating to a sulphur atom detected for resi " << test_resi << std::endl;
			num_hbonds++;
		}

		if ( test_rsd.atom_type(static_cast< int >(test_atom_num)).is_acceptor() ) {
			//if(bunsat_xyz.distance(test_xyz) < 3.3) {
			//if (check_AHD_angle(pose, bunsat_candidate_atom_id, id::AtomID( test_atom_num, test_resi)))
			//TR << "Checking if candidate sites donates to acceptor" << std::endl;
			if ( don_geom_check(pose, bunsat_rsd, bunsat_atom_num, bunsat_xyz, test_xyz) ) {
				TR << "HBond donating to a standard acceptor detected for resi " << test_resi << std::endl;
				num_hbonds++;
			}
		}

	}
}

/*
return false if a plausible hbond between bunsat_candidate_atom_id as an acceptor
and test_resi as an donor is found, otherwise return true
*/
void
BuriedUnsatisfiedPolarsCalculator2::bunsat_acc_nbr_residue_check(
	pose::Pose const & pose,
	id::AtomID const & bunsat_candidate_atom_id,
	Residue const & bunsat_rsd,
	Vector const & bunsat_xyz,
	Size const & test_resi,
	Size & num_hbonds
)
{
	Size const bunsat_resi = bunsat_candidate_atom_id.rsd();
	Size const bunsat_atom_num = bunsat_candidate_atom_id.atomno();
	std::string const bunsat_atom_name =
		bunsat_rsd.atom_name(static_cast< int >(bunsat_atom_num));

	Residue const & test_rsd(pose.residue(test_resi));

	AtomIndices const & test_hpos_polar = test_rsd.Hpos_polar();
	for ( auto test_atom_it = test_hpos_polar.begin();
			test_atom_it < test_hpos_polar.end(); ++test_atom_it ) {
		Size const test_atom_num = *test_atom_it;
		// exclude neighboring backbone-backbone
		std::string const test_atom_name(
			test_rsd.atom_name(static_cast< int >(test_atom_num)));
		if ( adjacent_bbbb_check(bunsat_resi, bunsat_atom_name, test_resi, test_atom_name) ) {
			continue;
		}

		// exclude self sidechain-sidechain
		if ( self_scsc(bunsat_rsd, bunsat_resi, bunsat_atom_num, test_rsd, test_resi, test_atom_num) ) {
			continue;
		}

		Vector test_xyz = test_rsd.atom(test_atom_num).xyz();

		// check if coordinating metal
		if ( metal_check(test_rsd, bunsat_xyz, test_xyz) ) {
			//TR << "HBond accepting from a metal ion detected" << std::endl;
			num_hbonds++;
		}
		if ( test_rsd.atom_type(static_cast< int >(test_atom_num)).is_polar_hydrogen() ) {
			//if (check_AHD_angle(pose, bunsat_candidate_atom_id, id::AtomID( test_atom_num, test_resi)))
			if ( acc_geom_check(pose, bunsat_xyz, test_rsd, test_atom_num, test_xyz) ) {
				TR << "HBond accepting from a standard donor detected for resi " << test_resi << std::endl;
				num_hbonds++;
			}
		}
	}
}

/*
return true if bunsat_xyz within h-bond distance to test_xyz and test_rsd is a metal ion,
otherwise return false
*/
bool
BuriedUnsatisfiedPolarsCalculator2::metal_check(
	Residue const & test_rsd,
	Vector const & bunsat_xyz,
	Vector const & test_xyz
) const
{
	std::string test_rsd_name = test_rsd.name();
	if ( test_rsd_name == "CA" || test_rsd_name == "MG" || test_rsd_name == "ZN" ||
			test_rsd_name == "FE" || test_rsd_name == "MN" || test_rsd_name == "NA" ) {
		return (bunsat_xyz.distance(test_xyz) < metal_dist_cutoff_); //metal_dist_cutoff
	} else {
		return false;
	}
}

/*
return true if bunsat_atom_name on bunsat_resi and test_atom_name on test_resi
are adjacent backbone atoms, otherwise false
*/
bool
BuriedUnsatisfiedPolarsCalculator2::adjacent_bbbb_check(
	Size const & bunsat_resi,
	std::string const & bunsat_atom_name,
	Size const & test_resi,
	std::string const & test_atom_name
) const
{
	bool adjacent = std::abs(static_cast< int >(bunsat_resi - test_resi)) <= 1;
	bool bunsat_bb = (bunsat_atom_name == "H" || bunsat_atom_name == "O");
	bool test_bb = (test_atom_name == "H" && test_atom_name == "O");
	return (adjacent && bunsat_bb && test_bb);
}

/*
return true if bunsat_atom_num on bunsat_resi and test_atom_num on test_resi
are both SC atoms in the same residue, otherwise false
*/
bool
BuriedUnsatisfiedPolarsCalculator2::self_scsc(
	Residue const & bunsat_rsd,
	Size const & bunsat_resi,
	Size const & bunsat_atom_num,
	Residue const & test_rsd,
	Size const & test_resi,
	Size const & test_atom_num
) const
{
	bool same_res = bunsat_resi == test_resi;
	bool bunsat_sc = bunsat_atom_num >= bunsat_rsd.first_sidechain_atom();
	bool test_sc = test_atom_num >= test_rsd.first_sidechain_atom();
	return (same_res && bunsat_sc && test_sc);
}

/*
return true if bunsat_xyz is within hbond distance of test_xyz and test_xyz
is a sulphur atom
*/
bool
BuriedUnsatisfiedPolarsCalculator2::sulphur_bond_check(
	Residue const & test_rsd,
	Size const & test_atom_num,
	Vector const & bunsat_xyz,
	Vector const & test_xyz
) const
{
	if ( test_rsd.atom_type(static_cast< int >(test_atom_num)).element() == "S" ) {
		return(bunsat_xyz.distance(test_xyz) < sulph_dist_cutoff_); //suphur_dist_cutoff
	} else {
		return false;
	}
}

/*
return true if bunsat_atom_num in bunsat_rsd as an hbond donor to an acceptor
at test_xyz has acceptable hbond geometry, otherwise return false
*/
bool
BuriedUnsatisfiedPolarsCalculator2::don_geom_check(
	core::pose::Pose const &,
	Residue const & bunsat_rsd,
	Size const & bunsat_atom_num,
	Vector const & bunsat_xyz,
	Vector const & test_xyz
) const
{
	for ( Size bunsat_H_atom_num = bunsat_rsd.attached_H_begin(static_cast< int >(bunsat_atom_num));
			bunsat_H_atom_num <= bunsat_rsd.attached_H_end(static_cast< int >(bunsat_atom_num));
			++bunsat_H_atom_num ) {
		//TR << "checking site attached H for donor geom" << std::endl;
		Vector bunsat_H_xyz = bunsat_rsd.atom(bunsat_H_atom_num).xyz();
		Real AHD_angle = numeric::conversions::degrees(
			angle_of(bunsat_xyz, bunsat_H_xyz, test_xyz));
		Real dist = bunsat_H_xyz.distance(test_xyz);
		//TR << "bunsat_xyz = " << bunsat_xyz[0] << "," << bunsat_xyz[1] << "," << bunsat_xyz[2] << "; bunsat_H_xyz = " << bunsat_H_xyz[0] << "," << bunsat_H_xyz[1] << "," << bunsat_H_xyz[2] << "; test_xyz = " << test_xyz[0] << "," << test_xyz[1] << "," << test_xyz[2] << std::endl;
		//TR << "AH dist = " << dist << ", AHD angle = " << AHD_angle << std::endl;
		//if(dist >= option[OptionKeys::bunsat_calc2::dist_cutoff]) TR << "failed dist check" << std::endl;
		//if (AHD_angle <= option[OptionKeys::bunsat_calc2::AHD_cutoff]) TR << "failed AHD check" << std::endl;
		if ( dist < dist_cutoff_ && AHD_angle > AHD_cutoff_ ) { //dist_cutoff and AHD_cutoff
			return true;
		} else if ( (bunsat_rsd.name() == "SER" || bunsat_rsd.name3() == "THR" || bunsat_rsd.name3() == "TYR")
				&& bunsat_atom_num >= bunsat_rsd.first_sidechain_atom() ) {
			/*
			HXL hydrogen placement is ambiguous, for serine or threonine do
			a second check based only on heavy atom distance
			*/
			Real dist = bunsat_xyz.distance(test_xyz);
			return dist < hxl_dist_cutoff_; // hydroxyl_dist_cutoff
		}
	}
	return false;
}


/// @brief return true if bunsat_atom_num in bunsat_rsd as an hbond acceptor to a donor
/// at test_xyz has acceptable hbond geometry, otherwise return false
bool
BuriedUnsatisfiedPolarsCalculator2::acc_geom_check(
	core::pose::Pose const &,
	Vector const & bunsat_xyz,
	Residue const & test_rsd,
	Size const & test_atom_num,
	Vector const & test_xyz
) const
{
	Vector test_base_xyz = test_rsd.atom(
		(test_rsd.atom_base(static_cast< int >(test_atom_num)))).xyz();
	Real AHD_angle = numeric::conversions::degrees(
		angle_of(bunsat_xyz, test_xyz, test_base_xyz));
	Real dist = bunsat_xyz.distance(test_xyz);
	//dist_cutoff and AHD_cutoff
	if ( dist < dist_cutoff_
			&& AHD_angle > AHD_cutoff_ ) {
		return true;
	} else {
		/*
		HXL hydrogen placement is ambiguous, for serine or threonine do
		a second check based only on heavy atom distance
		*/
		if ( (test_rsd.name() == "SER" || test_rsd.name3() == "THR" || test_rsd.name3() == "TYR")
				&& test_atom_num >= test_rsd.first_sidechain_atom() ) {
			Real dist = bunsat_xyz.distance(test_base_xyz);
			return dist < hxl_dist_cutoff_; // hydroxyl_dist_cutoff
		} else {
			return false;
		}
	}
}

core::Size
BuriedUnsatisfiedPolarsCalculator2::satisfaction_cutoff( std::string atom_type )
{

	//according to jk, buried hydroxyls are often seen making only one hydrogen bond. also, ether oxygens often are bad h-bond acceptors
	if ( atom_type == "OH" ) return 2;

	//backbone oxygens also only have one h-bbond in most secondary structure elements
	else if ( atom_type == "OCbb" ) return 2;

	else if ( atom_type ==  "S" ) return 2;

	//everything else we expect to have 3 bonded/h-bonded neighbours to count as satisfied
	else return 3;


}

void
BuriedUnsatisfiedPolarsCalculator2::show() {
	std::stringstream sstream;
	sstream << "Buns2calc with name_of_weak_bunsat_calc_ = "
		<< name_of_weak_bunsat_calc_ << ", layered_sasa_ = "
		<< layered_sasa_ << ", generous_hbonds_ = "
		<< generous_hbonds_ << ", sasa_burial_cutoff_ = "
		<< sasa_burial_cutoff_ << ", AHD_cutoff_ = "
		<< AHD_cutoff_ << ", dist_cutoff_ = "
		<< dist_cutoff_ << ".";
	TR << sstream.str() << std::endl;
}

} // namespace buns
} // namespace devel

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
devel::buns::BuriedUnsatisfiedPolarsCalculator2::BuriedUnsatisfiedPolarsCalculator2() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
devel::buns::BuriedUnsatisfiedPolarsCalculator2::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::EnergyDependentCalculator >( this ) );
	arc( CEREAL_NVP( all_bur_unsat_polars_ ) ); // Size
	arc( CEREAL_NVP( special_region_bur_unsat_polars_ ) ); // Size
	arc( CEREAL_NVP( atom_bur_unsat_ ) ); // core::id::AtomID_Map<_Bool>
	arc( CEREAL_NVP( residue_bur_unsat_polars_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( name_of_weak_bunsat_calc_ ) ); // std::string
	arc( CEREAL_NVP( special_region_ ) ); // std::set<Size>
	arc( CEREAL_NVP( layered_sasa_ ) ); // _Bool
	arc( CEREAL_NVP( generous_hbonds_ ) ); // _Bool
	arc( CEREAL_NVP( sasa_burial_cutoff_ ) ); // core::Real
	arc( CEREAL_NVP( AHD_cutoff_ ) ); // core::Real
	arc( CEREAL_NVP( dist_cutoff_ ) ); // core::Real
	arc( CEREAL_NVP( hxl_dist_cutoff_ ) ); // core::Real
	arc( CEREAL_NVP( sulph_dist_cutoff_ ) ); // core::Real
	arc( CEREAL_NVP( metal_dist_cutoff_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
devel::buns::BuriedUnsatisfiedPolarsCalculator2::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::EnergyDependentCalculator >( this ) );
	arc( all_bur_unsat_polars_ ); // Size
	arc( special_region_bur_unsat_polars_ ); // Size
	arc( atom_bur_unsat_ ); // core::id::AtomID_Map<_Bool>
	arc( residue_bur_unsat_polars_ ); // utility::vector1<Size>
	arc( name_of_weak_bunsat_calc_ ); // std::string
	arc( special_region_ ); // std::set<Size>
	arc( layered_sasa_ ); // _Bool
	arc( generous_hbonds_ ); // _Bool
	arc( sasa_burial_cutoff_ ); // core::Real
	arc( AHD_cutoff_ ); // core::Real
	arc( dist_cutoff_ ); // core::Real
	arc( hxl_dist_cutoff_ ); // core::Real
	arc( sulph_dist_cutoff_ ); // core::Real
	arc( metal_dist_cutoff_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( devel::buns::BuriedUnsatisfiedPolarsCalculator2 );
CEREAL_REGISTER_TYPE( devel::buns::BuriedUnsatisfiedPolarsCalculator2 )

CEREAL_REGISTER_DYNAMIC_INIT( devel_buns_BuriedUnsatisfiedPolarsCalculator2 )
#endif // SERIALIZATION
