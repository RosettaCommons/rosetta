// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/kevin/BuriedUnsatPolarsFinder2.cc
/// @brief
/// @details
/// @author Kevin Houlihan
/// @author Bryan Der

#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>

#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>
//#include <core/scoring/hbonds/hbonds.hh>
//#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/MetricValue.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>
#include <devel/vardist_solaccess/LoadVarSolDistSasaCalculatorMover.hh>
#include <numeric/xyzVector.hh>

#include <core/scoring/sasa.hh>
#include <numeric/conversions.hh> //degrees-radians

//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static THREAD_LOCAL basic::Tracer TR( "apps.pilot.bder.BuriedUnsatPolarsFinder" );


using namespace basic::options;
using namespace core;
using namespace core::conformation;
using namespace chemical;
using namespace utility;
using namespace core::pose::metrics;
using namespace protocols::toolbox::pose_metric_calculators;

typedef numeric::xyzVector< core::Real > Vector;

basic::options::BooleanOptionKey const use_varsoldist_sasa_calc("use_varsoldist_sasa_calc");
basic::options::BooleanOptionKey const use_regular_sasa_calc("use_regular_sasa_calc");
//basic::options::BooleanOptionKey const add_probe_radius("add_probe_radius");

basic::options::BooleanOptionKey const generous_hbond_settings("generous_hbond_settings");
basic::options::RealOptionKey const sasa_cutoff("bunsat_sasa_cutoff");
basic::options::RealOptionKey const AHD_cutoff("bunsat_AHD_cutoff");
basic::options::RealOptionKey const dist_cutoff("bunsat_dist_cutoff");
basic::options::RealOptionKey const hydroxyl_dist_cutoff("bunsat_hydroxyl_dist_cutoff");
basic::options::RealOptionKey const sulphur_dist_cutoff("bunsat_sulphur_dist_cutoff");
basic::options::RealOptionKey const metal_dist_cutoff("bunsat_metal_dist_cutoff");

/// @brief
class BuriedUnsatPolarsFinder : public protocols::moves::Mover {
public:
	BuriedUnsatPolarsFinder() {}
	
	virtual ~BuriedUnsatPolarsFinder(){};
	
	virtual
	void
	apply(pose::Pose & pose)
	{
		
		assert(option[use_varsoldist_sasa_calc] ||
			   option[use_regular_sasa_calc]);
		
		utility::file::FileName filename(pose.pdb_info()->name());
		std::string pdbname_base = filename.base();
		//std::string pymol_bunsat_selection  = "bunsat_"  + pdbname_base;
		//std::string pymol_bunsat_selection2 = "bunsat2_" + pdbname_base;
		
		using namespace core::scoring;
		
		if(basic::options::option[generous_hbond_settings])
			generous_hbond();
		
		core::scoring::ScoreFunctionOP scorefxn = get_score_function();
		scorefxn->score(pose);
		
		TR << "Registering num_Hbond Calculator" << std::endl;
		std::string num_hbonds_calc_name("num_hbonds_calc_name");
		CalculatorFactory::Instance().remove_calculator(num_hbonds_calc_name);
		CalculatorFactory::Instance().register_calculator(num_hbonds_calc_name,
			new NumberHBondsCalculator());
		
		std::string sasa_calc_name("sasa_calc_name");
		std::string bunsat_calc_name("bunsat_calc_name");
		Real sasa_buried = option[sasa_cutoff]; //sasa_cutoff
		
		CalculatorFactory::Instance().remove_calculator(sasa_calc_name);
		CalculatorFactory::Instance().remove_calculator(bunsat_calc_name);
		
		if(basic::options::option[use_varsoldist_sasa_calc]) {
			TR << "Registering VarSolDist SASA Calculator" << std::endl;
			CalculatorFactory::Instance().register_calculator(sasa_calc_name,
				new devel::vardist_solaccess::VarSolDistSasaCalculator());
		}
		else { // default to vanilla SasaCalculator
			TR << "Registering SASA Calculator" << std::endl;
			CalculatorFactory::Instance().register_calculator(sasa_calc_name,
				new pose::metrics::simple_calculators::SasaCalculator());
		}
		
		CalculatorFactory::Instance().register_calculator(bunsat_calc_name,
			new BuriedUnsatisfiedPolarsCalculator(
				sasa_calc_name, num_hbonds_calc_name, sasa_buried));
		
		basic::MetricValue< id::AtomID_Map< Real > > sasa_atomid_map;
		basic::MetricValue< id::AtomID_Map< bool > > bunsat_atomid_map;
		
		pose.metric(sasa_calc_name, "atom_sasa", sasa_atomid_map);
		pose.metric(bunsat_calc_name, "atom_bur_unsat", bunsat_atomid_map);
		
		//convert atomid-map to vector of atom ids
		vector1< id::AtomID > bunsat_atom_ids;
		
		for (Size res = 1; res <= (bunsat_atomid_map.value()).n_residue(); ++res) {
			for (Size atm = 1; atm <= (bunsat_atomid_map.value()).n_atom(res); ++atm) {
				if ((bunsat_atomid_map.value())(res, atm)) {
					core::id::AtomID atm_id(atm, res);
					bunsat_atom_ids.push_back(atm_id);
				}
			}
		}
		
		id::AtomID_Map< bool > bunsat_thorough_atomid_map(bunsat_atomid_map.value());
		bunsats_thorough_check(pose, bunsat_thorough_atomid_map);
		
		
		
		// TODO Remove all this, just to verify identical outputs
		/*
		 
		std::string pymol_bunsat_selection  = "bunsat_"  + pdbname_base;
		std::string pymol_bunsat_selection2 = "bunsat2_" + pdbname_base;
		
		//Print for Talaris2013 hbond detection only
		TR << "PYMOL_SELECTION: select " << pymol_bunsat_selection << ",";
		
		for(Size i(1); i <= bunsat_atom_ids.size(); ++i) {
			Size rosetta_resnum = bunsat_atom_ids[i].rsd();
			std::string atom_name = pose.residue(rosetta_resnum).atom_name(bunsat_atom_ids[i].atomno());
			
			Size pdb_resnum = pose.pdb_info()->number(rosetta_resnum);
			char chain = pose.pdb_info()->chain(rosetta_resnum);
			
			if(i != bunsat_atom_ids.size()) {
				TR << " " << pdbname_base << " and chain " << chain << " and resi " << pdb_resnum << " and name " << atom_name << " +";
			}
			else {
				TR << " " << pdbname_base << " and chain " << chain << " and resi " << pdb_resnum << " and name " << atom_name ;
			}
		}
		TR << std::endl;
		TR << "PYMOL_SELECTION: alter " << pymol_bunsat_selection << ", vdw=0.5" << std::endl;
		
		
		
		//Print for Talaris2013 hbond detection + distance-AHD hbond detection
		
		TR << "PYMOL_SELECTION_2: select " << pymol_bunsat_selection2 << ",";
		
		for(Size i(1); i <= bunsat_thorough_atom_ids.size(); ++i) {
			Size rosetta_resnum = bunsat_thorough_atom_ids[i].rsd();
			std::string atom_name = pose.residue(rosetta_resnum).atom_name(bunsat_thorough_atom_ids[i].atomno());
			
			Size pdb_resnum = pose.pdb_info()->number(rosetta_resnum);
			char chain = pose.pdb_info()->chain(rosetta_resnum);
			
			if(i != bunsat_thorough_atom_ids.size()) {
				TR << " " << pdbname_base << " and chain " << chain << " and resi " << pdb_resnum << " and name " << atom_name << " +";
			}
			else {
				TR << " " << pdbname_base << " and chain " << chain << " and resi " << pdb_resnum << " and name " << atom_name;
			}
		}
		TR << std::endl;
		TR << "PYMOL_SELECTION_2: alter " << pymol_bunsat_selection2 << ", vdw=0.5" << std::endl;
		
		
		
	  	//This will print a zero for buried_unsat_atom_ids_with_distAHD_check if the dist-AHD option was turned off
	  	TR << pose.pdb_info()->name() << " has " << bunsat_atom_ids.size() << " buns  and  " << bunsat_thorough_atom_ids.size() << " buns_with_distAHD_check" << std::endl;
		 */
		// TODO Remove above
		
		
		
		return;
	}
	
	void
	generous_hbond()
	const
	{
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		scoring::methods::EnergyMethodOptionsOP emopts(
			new scoring::methods::EnergyMethodOptions(
				scorefxn->energy_method_options()));
		emopts->hbond_options().decompose_bb_hb_into_pair_energies(true);
		emopts->hbond_options().use_hb_env_dep(false);
		emopts->hbond_options().bb_donor_acceptor_check(false);
		scorefxn->set_energy_method_options(*emopts);
	}
	
	void
	bunsats_thorough_check(
		pose::Pose const & pose,
		id::AtomID_Map< bool > & bunsat_thorough_atomid_map)
	const
	{
	  	//for(Size i(1); i <= bunsat_atomid_map.size(); ++i) {
			
			
		for (Size res = 1; res <= bunsat_thorough_atomid_map.n_residue(); ++res) {
			for (Size atm = 1; atm <= bunsat_thorough_atomid_map.n_atom(res); ++atm) {
				if (bunsat_thorough_atomid_map(res, atm))
				{
					id::AtomID bunsat_candidate_atom_id(atm, res);
					bunsat_thorough_atomid_map.set(bunsat_candidate_atom_id,
						single_bunsat_thorough_check(pose, bunsat_candidate_atom_id));
				}
			}
		}
	}

	bool
	single_bunsat_thorough_check(
		pose::Pose const & pose,
		id::AtomID const & bunsat_candidate_atom_id)
	const
	{
		Size const bunsat_resi = bunsat_candidate_atom_id.rsd();
		std::string const nbr_calc_name("nbr_calc");
		// remove unnecessary residues and guarantee necessary ones?
		core::Real dist_cutoff(8); //related to dist_cutoff
		CalculatorFactory::Instance().remove_calculator(nbr_calc_name);
		CalculatorFactory::Instance().register_calculator(nbr_calc_name,
			new NeighborsByDistanceCalculator(bunsat_resi, dist_cutoff));
		basic::MetricValue< std::set< Size > > neighbors;
		pose.metric(nbr_calc_name, "neighbors", neighbors);
		
		Residue const & bunsat_rsd(pose.residue(bunsat_resi));
		Size const bunsat_atom_num = bunsat_candidate_atom_id.atomno();
		
		Vector bunsat_xyz = bunsat_rsd.atom(bunsat_atom_num).xyz();
		
		bool bunsat_is_donor(bunsat_rsd.atom_type(
			static_cast< int >(bunsat_atom_num)).is_donor());
		
		if (bunsat_is_donor)
		{
			for(std::set<Size>::const_iterator test_res_it = neighbors.value().begin();
				test_res_it != neighbors.value().end(); ++test_res_it)
			{
				Size const test_resi = *test_res_it;
				if (bunsat_donor_nbr_residue_check(
					pose, bunsat_candidate_atom_id, bunsat_rsd, bunsat_xyz, test_resi))
						return true;
			}
		}
		else
		{
			for(std::set<Size>::const_iterator test_res_it = neighbors.value().begin();
				test_res_it != neighbors.value().end(); ++test_res_it)
			{
				Size const test_resi = *test_res_it;
				if (bunsat_acc_nbr_residue_check(
					pose, bunsat_candidate_atom_id, bunsat_rsd, bunsat_xyz, test_resi))
				{
					return true;
				}
			}
		}
		return false;
	}	
	
	bool
	bunsat_donor_nbr_residue_check(
		pose::Pose const & pose,
		id::AtomID const & bunsat_candidate_atom_id,
		Residue const & bunsat_rsd,
		Vector const & bunsat_xyz,
		Size const test_resi)
	const
	{
		Size const bunsat_resi = bunsat_candidate_atom_id.rsd();
		Size const bunsat_atom_num = bunsat_candidate_atom_id.atomno();
		std::string const bunsat_atom_name =
			bunsat_rsd.atom_name(static_cast< int >(bunsat_atom_num));
		
		Residue const & test_rsd(pose.residue(test_resi));
		
		
		AtomIndices const & test_accpt_pos = test_rsd.accpt_pos();
		for(vector1< Size >::const_iterator test_atom_it = test_accpt_pos.begin();
			test_atom_it < test_accpt_pos.end(); ++test_atom_it)
		{
			Size const test_atom_num = *test_atom_it;
			
			// exclude neighboring backbone-backbone
			std::string const test_atom_name(
				test_rsd.atom_name(static_cast< int >(test_atom_num)));
			if (adjacent_bbbb_check(bunsat_resi, bunsat_atom_name, test_resi, test_atom_name))
			{
				continue;
			}

			// exclude self sidechain-sidechain
			if (self_scsc(bunsat_rsd, bunsat_resi, bunsat_atom_num, test_rsd, test_resi, test_atom_num))
			{
				continue;
			}
			
			Vector test_xyz = test_rsd.atom(test_atom_num).xyz();
			// check for sulfur-bonds
			if (sulphur_bond_check(test_rsd, test_atom_num, bunsat_xyz, test_xyz))
			{
				return true;
			}
			
			if (test_rsd.atom_type(static_cast< int >(test_atom_num)).is_acceptor()){
				//if(bunsat_xyz.distance(test_xyz) < 3.3) {
				//if (check_AHD_angle(pose, bunsat_candidate_atom_id, id::AtomID( test_atom_num, test_resi)))
				if (don_geom_check(pose, bunsat_rsd, bunsat_atom_num, bunsat_xyz, test_xyz))
				{
					return true;
				}
			}
			
		}
							   
		return false;
	}
	
	bool
	bunsat_acc_nbr_residue_check(
		pose::Pose const & pose,
		id::AtomID const & bunsat_candidate_atom_id,
		Residue const & bunsat_rsd,
		Vector const & bunsat_xyz,
		Size const & test_resi)
	const
	{
		Size const bunsat_resi = bunsat_candidate_atom_id.rsd();
		Size const bunsat_atom_num = bunsat_candidate_atom_id.atomno();
		std::string const bunsat_atom_name =
		bunsat_rsd.atom_name(static_cast< int >(bunsat_atom_num));
		
		Residue const & test_rsd(pose.residue(test_resi));
		
		AtomIndices const & test_hpos_polar = test_rsd.Hpos_polar();
		for(vector1< Size>::const_iterator test_atom_it = test_hpos_polar.begin();
			test_atom_it < test_hpos_polar.end(); ++test_atom_it)
		{
			Size const test_atom_num = *test_atom_it;
			
			// exclude neighboring backbone-backbone
			std::string const test_atom_name(
				test_rsd.atom_name(static_cast< int >(test_atom_num)));
			if (adjacent_bbbb_check(bunsat_resi, bunsat_atom_name, test_resi, test_atom_name))
			{
				continue;
			}
				
			// exclude self sidechain-sidechain
			if (self_scsc(bunsat_rsd, bunsat_resi, bunsat_atom_num, test_rsd, test_resi, test_atom_num))
			{
				continue;
			}
				
			Vector test_xyz = test_rsd.atom(test_atom_num).xyz();

			// check if coordinating metal
			if (metal_check(test_rsd, bunsat_xyz, test_xyz))
			{
				return true;
			}
				
			if (test_rsd.atom_type(static_cast< int >(test_atom_num)).is_polar_hydrogen()) {
				
				//if (check_AHD_angle(pose, bunsat_candidate_atom_id, id::AtomID( test_atom_num, test_resi)))
				if (acc_geom_check(pose, bunsat_xyz, test_rsd, test_atom_num, test_xyz))
				{
					return true;
				}
			}
		}
		
		return false;
	}
	
	bool
	metal_check(
		Residue const & test_rsd,
		Vector const & bunsat_xyz,
		Vector const & test_xyz)
	const
	{
		std::string test_rsd_name = test_rsd.name();
		if(test_rsd_name == "CA" || test_rsd_name == "MG" || test_rsd_name == "ZN" ||
		   test_rsd_name == "FE" || test_rsd_name == "MN" || test_rsd_name == "NA")
		{
			return (bunsat_xyz.distance(test_xyz) < option[metal_dist_cutoff]); //metal_dist_cutoff
		}
		else
		{
			return false;
		}
	}
	
	// TODO fill this in
	bool
	AHD_and_dist_check(){return false;}
	
	bool
	adjacent_bbbb_check(
		Size const & bunsat_resi,
		std::string const & bunsat_atom_name,
		Size const & test_resi,
		std::string const & test_atom_name)
	const
	{
		bool adjacent = std::abs(static_cast< int >(bunsat_resi - test_resi)) <= 1;
		bool bunsat_bb = (bunsat_atom_name == "H" || bunsat_atom_name == "O");
		bool test_bb = (test_atom_name == "H" && test_atom_name == "O");
		return (adjacent && bunsat_bb && test_bb);
	}
	
	bool
	self_scsc(
		Residue const & bunsat_rsd,
		Size const & bunsat_resi,
		Size const & bunsat_atom_num,
		Residue const & test_rsd,
		Size const & test_resi,
		Size const & test_atom_num)
	const
	{
		bool same_res = bunsat_resi == test_resi;
		bool bunsat_sc = bunsat_atom_num >= bunsat_rsd.first_sidechain_atom();
		bool test_sc = test_atom_num >= test_rsd.first_sidechain_atom();
		return (same_res && bunsat_sc && test_sc);
	}
	
	bool sulphur_bond_check(
		Residue const & test_rsd,
		Size const & test_atom_num,
		Vector const & bunsat_xyz,
		Vector const & test_xyz)
	const
	{
		if(test_rsd.atom_type(static_cast< int >(test_atom_num)).element() == "S")
		{
			return(bunsat_xyz.distance(test_xyz) < option[sulphur_dist_cutoff]); //suphur_dist_cutoff
		}
		else
		{
			return false;
		}
	}
	
	bool
	don_geom_check(
		core::pose::Pose const & pose,
		Residue const & bunsat_rsd,
		Size const & bunsat_atom_num,
		Vector const & bunsat_xyz,
		Vector const & test_xyz)
	const
	{
		for (Size bunsat_H_atom_num = bunsat_rsd.attached_H_begin(static_cast< int >(bunsat_atom_num));
			 bunsat_H_atom_num <= bunsat_rsd.attached_H_end(static_cast< int >(bunsat_atom_num));
			 ++bunsat_H_atom_num)
		{
			Vector bunsat_H_xyz = bunsat_rsd.atom(
				(bunsat_rsd.atom_base(static_cast< int >(bunsat_H_atom_num)))).xyz();
			Real AHD_angle = numeric::conversions::degrees(
				angle_of(bunsat_xyz, bunsat_H_xyz, test_xyz));
			Real dist = bunsat_H_xyz.distance(test_xyz);
			if (dist < option[dist_cutoff] && AHD_angle > option[AHD_cutoff]) //dist_cutoff and AHD_cutoff
			{
				return true;
			}
		}
		return false;
	}
	
	bool
	acc_geom_check(
		core::pose::Pose const & pose,
		Vector const & bunsat_xyz,
		Residue const & test_rsd,
		Size const & test_atom_num,
		Vector const & test_xyz)
	const
	{
		Vector test_base_xyz = test_rsd.atom(
			(test_rsd.atom_base(static_cast< int >(test_atom_num)))).xyz();
		Real AHD_angle = numeric::conversions::degrees(
			angle_of(bunsat_xyz, test_xyz, test_base_xyz));
		Real dist = bunsat_xyz.distance(test_xyz);
		if (dist < option[dist_cutoff] && AHD_angle > option[AHD_cutoff]) //dist_cutoff and AHD_cutoff
		{
			return true;
		}
		else
		{
			/*
			 HXL hydrogen placement is ambiguous, for serine or threonine do
			 a second check based only on heavy atom distance
			 */
			if ((test_rsd.name() == "SER" || test_rsd.name3() == "THR")
				&& test_atom_num >= test_rsd.first_sidechain_atom())
			{
				Real dist = bunsat_xyz.distance(test_base_xyz);
				return dist < option[hydroxyl_dist_cutoff]; // hydroxyl_dist_cutoff
			}
			else
			{
				return false;
			}
		}
	}
	
	virtual
	std::string
	get_name() const { return "BuriedUnsatPolarsFinder"; }
	
	
	
	/*
	// TODO remove
	
	virtual
	Size  //return 1 if there is an hbond detected, return 0 if not
	check_AHD_angle(
		core::pose::Pose const & pose,
		id::AtomID bunsat_atid,
		id::AtomID hbond_atid)
	const
	{
		Vector bunsat_xyz = pose.residue(bunsat_atid.rsd()).atom(bunsat_atid.atomno()).xyz();
		Vector hbond_xyz = pose.residue(hbond_atid.rsd()).atom(hbond_atid.atomno()).xyz();
		Real dist = bunsat_xyz.distance(hbond_xyz);
		
		
		//Determine hydrogen atom by closest distance to the donor and acceptor
		//A better solution would be to only check the polar hydrogens bonded to the donor heavy atom
		Vector hydrogen_xyz(0, 0, 0);
		Real h_dist(999.9);
		Size closest_h_atomno(99);
		Real closest_combined_h_dist (999.9);
		
		if( pose.residue(bunsat_atid.rsd()).atom_type(bunsat_atid.atomno()).is_donor() ) {
			for( Size i(1); i<=pose.residue(bunsat_atid.rsd()).natoms(); ++i ) {
				if( pose.residue(bunsat_atid.rsd()).atom_type(i).is_polar_hydrogen() ) {
					Vector h_xyz = pose.residue(bunsat_atid.rsd()).atom(i).xyz();
					Real combined_h_dist = h_xyz.distance(hbond_xyz) + h_xyz.distance(bunsat_xyz);
					if(combined_h_dist < closest_combined_h_dist) {
						closest_combined_h_dist = combined_h_dist;
						closest_h_atomno = i;
						hydrogen_xyz = h_xyz;
					}
				}
			}
			h_dist = hydrogen_xyz.distance(bunsat_xyz);
		}
		else if( pose.residue(bunsat_atid.rsd()).atom_type(bunsat_atid.atomno()).is_acceptor() ) {
			for( Size i(1); i<=pose.residue(hbond_atid.rsd()).natoms(); ++i ) {
				if( pose.residue(hbond_atid.rsd()).atom_type(i).is_polar_hydrogen() ) {
					Vector h_xyz = pose.residue(hbond_atid.rsd()).atom(i).xyz();
					Real combined_h_dist = h_xyz.distance(bunsat_xyz) + h_xyz.distance(hbond_xyz);
					if(combined_h_dist < closest_combined_h_dist) {
						closest_combined_h_dist = combined_h_dist;
						closest_h_atomno = i;
						hydrogen_xyz = h_xyz;
					}
				}
			}
			h_dist = hydrogen_xyz.distance(hbond_xyz);
		}
		
		//Measure the AHD angle (Acceptor-Hydrogen-Donor)
		Real AHD_angle = numeric::conversions::degrees(angle_of(bunsat_xyz, hydrogen_xyz, hbond_xyz));
		
		if( AHD_angle > 120.0 && h_dist < 2.4 ) {
			
			TR << "DIST and AHD_angle: " << pose.pdb_info()->name() << " " << bunsat_atid << " " << hbond_atid << "  " << dist << " " << AHD_angle << std::endl;
			return 1;
		}
		else {
			return 0;
		}
	}
	
	 */
	
	
	
private:
	
	// basic::MetricValue< id::AtomID_Map< Real > > atom_sasa_;
	// basic::MetricValue< id::AtomID_Map< Size > > atom_hbonds_;
	
};

typedef utility::pointer::owning_ptr<BuriedUnsatPolarsFinder> BuriedUnsatPolarsFinderOP;

int main(int argc, char* argv[])
{
	try {
		using basic::options::option;
		option.add(use_varsoldist_sasa_calc, "var sol dist sasa calculator").def(true);
		option.add(use_regular_sasa_calc, "regular sasa calculator").def(false);
		//option.add(add_probe_radius, "include_probe_radius_in_atom_radii").def(false);
		
		option.add(generous_hbond_settings, "decompose, allow sc/bb in secstruct, no env_dep").def(true);

		option.add(sasa_cutoff, "SASA cutoff to be considered buried").def(0.01);
		option.add(AHD_cutoff, "AHD angle above which hbonds are considered "
			"for geometry-only based H-bond checking").def(120);
		option.add(dist_cutoff, "acceptor-hydrogen distance below which hbonds "
				   "are considered for geometry-only based H-bond checking").def(3.0);
		option.add(hydroxyl_dist_cutoff, "acceptor to hydroxyl oxygen ditance cutoff"
				   "for detecting h-bonds by heavyatom positions; meant to "
				   "account for ambiguity in hydrogen placement").def(3.5);
		option.add(sulphur_dist_cutoff, "sulphur to hydrogen distance cutoff below "
				   "which an h-bond is detected").def(3.3);
		option.add(metal_dist_cutoff, "metal atom to hydrogen distance cutoff below "
				   "which an h-bond is detected").def(2.7);
		
		
		devel::init(argc, argv);
		protocols::jd2::JobDistributor::get_instance()->go(new BuriedUnsatPolarsFinder);
		
		TR << "Recommended option:   -sasa_calculator_probe_radius 1.2" << std::endl;
		
		TR << "************************d**o**n**e**************************************" << std::endl;
		
	} catch (utility::excn::EXCN_Base const & e) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	
	return 0;
}
