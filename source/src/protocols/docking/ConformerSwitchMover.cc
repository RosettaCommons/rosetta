// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ConformerSwitchMover.cc
/// @brief code for the conformer switch mover in ensemble docking
/// @author Sid Chaudhury
/// @author Modified by Daisuke Kuroda


#include <protocols/docking/ConformerSwitchMover.hh>
#include <protocols/docking/ConformerSwitchMoverCreator.hh>

// Rosetta Headers
#include <basic/datacache/CacheableStringFloatMap.hh>
#include <core/id/types.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/types.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/docking/DockFilters.hh>
#include <protocols/docking/DockingEnsemble.hh>
#include <protocols/scoring/Interface.hh>


// Random number generator
#include <numeric/random/random.hh>
//
#include <string>

#include <basic/Tracer.hh>


#include <core/chemical/ChemicalManager.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/datacache/BasicDataCache.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
using basic::T;
using basic::Error;
using basic::Warning;


static THREAD_LOCAL basic::Tracer TR( "protocols.docking.ConformerSwitchMover" );

namespace protocols {
namespace docking {


ConformerSwitchMover::ConformerSwitchMover() :
	moves::Mover(),
	random_conformer_(false),
	temperature_(0.8)
{
	moves::Mover::type( "ConformerSwitchMover" );
}


std::string
ConformerSwitchMoverCreator::keyname() const
{
	return ConformerSwitchMoverCreator::mover_name();
}

protocols::moves::MoverOP
ConformerSwitchMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ConformerSwitchMover );
}

std::string
ConformerSwitchMoverCreator::mover_name()
{
	return "ConformerSwitchMover";
}

//constructor with arguments
ConformerSwitchMover::ConformerSwitchMover(
	protocols::docking::DockingEnsembleOP ensemble,
	bool random_conformer
) :
	moves::Mover(),
	random_conformer_( random_conformer ),
	temperature_(0.8)
{
	ensemble_ = ensemble;

	lowres_filter_ = protocols::docking::DockingLowResFilterOP( new protocols::docking::DockingLowResFilter() );
}

void ConformerSwitchMover::set_temperature( core::Real temp_in )
{
	temperature_ = temp_in;
}

void ConformerSwitchMover::apply( core::pose::Pose & pose )
{
	using namespace core::pose::datacache;
	core::Size conf_num = 1;

	(*(ensemble_->scorefxn_low()))(pose);

	// the conformer used can either be random, or generated from the probability table of energies of the conformers
	// to use a random conformer, random_conformer must be set to true
	if ( !random_conformer_ && lowres_filter_->apply( pose ) ) {  //only calculate partition function if filters are passed
		GenerateProbTable( pose );
		core::Real rand_num( numeric::random::rg().uniform() );
		for ( Size i = 1; i <= ensemble_->size(); i++ ) {
			if ( (rand_num >= prob_table_[i]) ) conf_num++;
		}
	} else {
		conf_num = numeric::random::rg().random_range( 1, ensemble_->size() );
	}

	TR << "Switching partner with conformer: " << conf_num << std::endl;

	switch_conformer( pose, conf_num );
	ensemble_->set_current_confnum( conf_num );

	// make sure that the pose has ARBITRARY_FLOAT_DATA in the DataCache
	if ( !pose.data().has( ( CacheableDataType::ARBITRARY_FLOAT_DATA ) ) ) {
		using namespace basic::datacache;
		pose.data().set(
			CacheableDataType::ARBITRARY_FLOAT_DATA,
			DataCache_CacheableData::DataOP( new basic::datacache::CacheableStringFloatMap() )
		);
	}

	basic::datacache::CacheableStringFloatMapOP data
		= utility::pointer::dynamic_pointer_cast< basic::datacache::CacheableStringFloatMap >
		( pose.data().get_ptr(CacheableDataType::ARBITRARY_FLOAT_DATA) );

	data->map()[ ensemble_->partner() ] = ensemble_->lowres_reference_energy(conf_num);
	pose.data().set( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA, data );
}

void ConformerSwitchMover::show(std::ostream & output) const
{
	Mover::show(output);
	output << "Temperature:           " << get_temperature() << std::endl <<
		"Use random conformer?: " << (use_random_conformer() ? "True" : "False") << std::endl;
}

void ConformerSwitchMover::GenerateProbTable( core::pose::Pose & pose )
{

	core::pose::Pose complex_pose = pose;
	utility::vector1< core::Real > e_table;
	core::Real partition_sum(0.0);

	prob_table_.clear();

	for ( core::Size i = 1; i <= ensemble_->size(); i++ ) {
		switch_conformer( complex_pose, i );
		// the score takes the reference energy into consideration
		core::Real complex_score = (*(ensemble_->scorefxn_low()))(complex_pose) - ensemble_->lowres_reference_energy(i);
		complex_pose = pose;
		e_table.push_back( complex_score );
	}

	core::Real min_energy(0.0);
	for ( core::Size i = 1; i <= ensemble_->size(); i++ ) {
		if ( e_table[i] <= min_energy ) min_energy = e_table[i];
	}

	for ( core::Size i = 1; i <= ensemble_->size(); i++ ) {
		e_table[i] = std::exp((-1*(e_table[i] - min_energy))/temperature_);
		partition_sum += e_table[i];
	}

	prob_table_.push_back( e_table[1]/partition_sum );
	for ( core::Size i = 2; i <=ensemble_->size(); i++ ) {
		prob_table_.push_back( prob_table_[i-1] + (e_table[i]/partition_sum) );
	}

}

void ConformerSwitchMover::switch_conformer(
	core::pose::Pose & pose,
	core::Size conf_num
)
{
	core::pose::Pose new_conf = ensemble_->get_conformer_cen(conf_num); // new_conf is centroid

	(*(ensemble_->scorefxn_low()))(pose);
	scoring::Interface interface( ensemble_->jump_id() );

	interface.calculate( pose );

	utility::vector1<Size>conf_interface;
	for ( Size i = ensemble_->start_res(); i <= ensemble_->end_res(); i++ ) {
		if ( interface.is_interface(i) ) conf_interface.push_back( i );
	}

	if ( conf_interface.size() <= 5 ) {
		conf_interface.clear();
		for ( Size i = ensemble_->start_res(); i <= ensemble_->end_res(); i++ ) {
			conf_interface.push_back( i );
		}
	}

	core::id::AtomID_Map< core::id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, new_conf, core::id::BOGUS_ATOM_ID ); // maps every atomid to bogus

	for ( Size i = 1; i <= conf_interface.size(); i++ ) {
		Size new_conf_resnum = conf_interface[i]-ensemble_->start_res()+1;
		Size pose_resnum = conf_interface[i];
		core::id::AtomID const id1( new_conf.residue(new_conf_resnum).atom_index("CA"), new_conf_resnum );
		core::id::AtomID const id2( pose.residue(pose_resnum).atom_index("CA"), pose_resnum );
		atom_map[ id1 ] = id2;
	}

	core::scoring::superimpose_pose( new_conf, pose, atom_map );
	pose.copy_segment( ensemble_->conf_size(), new_conf, ensemble_->start_res(), 1);
}

std::string ConformerSwitchMover::get_name() const {
	return ConformerSwitchMoverCreator::mover_name();
}

std::ostream &operator<< (std::ostream & output, ConformerSwitchMover const & mover)
{
	mover.show(output);
	return output;
}

}  // namespace docking
}  // namespace protocols
