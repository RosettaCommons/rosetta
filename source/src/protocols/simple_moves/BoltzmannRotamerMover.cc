// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file BoltzmannRotamerMover.cc
/// @brief implementation of BoltzmannRotamerMover class and functions
/// @author Noah Ollikainen (nollikai@gmail.com)

// Unit headers
#include <protocols/simple_moves/BoltzmannRotamerMover.hh>
#include <protocols/simple_moves/BoltzmannRotamerMoverCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/packer_neighbors.hh>
#include <utility/graph/Graph.hh>
#include <core/conformation/Residue.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.fwd.hh>

// Parser headers
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// Numeric Headers
#include <numeric/random/random.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.simple_moves.BoltzmannRotamerMover" );

bool compare_values(
	const std::pair<core::Size, core::Real> &lhs,
	const std::pair<core::Size, core::Real> &rhs) {
	return lhs.second > rhs.second;
}

namespace protocols {
namespace simple_moves {

// XRW TEMP std::string
// XRW TEMP BoltzmannRotamerMoverCreator::keyname() const {
// XRW TEMP  return BoltzmannRotamerMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP BoltzmannRotamerMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new BoltzmannRotamerMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP BoltzmannRotamerMover::mover_name() {
// XRW TEMP  return "BoltzmannRotamerMover";
// XRW TEMP }

// default constructor
BoltzmannRotamerMover::BoltzmannRotamerMover() : protocols::moves::Mover()
{
	protocols::moves::Mover::type( "BoltzmannRotamer" );
	resnum_ = 0;
	ligand_resnum_ = 0;
	ligand_weight_ = 1.0;
	temperature_ = 1.0;
	bias_sampling_ = true;
	randomize_resnum_ = false;
	bump_check_ = true;
}

// constructor with arguments
BoltzmannRotamerMover::BoltzmannRotamerMover(
	ScoreFunctionCOP scorefxn_in,
	PackerTaskOP & task_in
) : protocols::moves::Mover(), scorefxn_(std::move( scorefxn_in )), factory_( /* NULL */ ), show_packer_task_( false )
{
	protocols::moves::Mover::type( "BoltzmannRotamer" );
	task_ = task_in->clone();
	resnum_ = 0;
	ligand_resnum_ = 0;
	ligand_weight_ = 1.0;
	temperature_ = 1.0;
	bias_sampling_ = true;
	randomize_resnum_ = false;
	bump_check_ = true;
}

// constructor with arguments
BoltzmannRotamerMover::BoltzmannRotamerMover(
	ScoreFunctionCOP scorefxn_in,
	TaskFactoryCOP factory_in
) : protocols::moves::Mover(), scorefxn_(std::move( scorefxn_in )), task_( /* NULL */ ), factory_(std::move( factory_in )), show_packer_task_( false )
{
	protocols::moves::Mover::type( "BoltzmannRotamer" );
	resnum_ = 0;
	ligand_resnum_ = 0;
	ligand_weight_ = 1.0;
	temperature_ = 1.0;
	bias_sampling_ = true;
	randomize_resnum_ = false;
	bump_check_ = true;
}

// copy constructor
BoltzmannRotamerMover::BoltzmannRotamerMover( BoltzmannRotamerMover const & rval ):
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( rval ),
	scorefxn_( rval.scorefxn_ ),
	task_( rval.task_ ),
	factory_( rval.factory_ ),
	show_packer_task_( rval.show_packer_task_ ),
	resnum_( rval.resnum_ ),
	ligand_resnum_( rval.ligand_resnum_ ),
	ligand_weight_( rval.ligand_weight_ ),
	temperature_( rval.temperature_ ),
	bias_sampling_( rval.bias_sampling_ ),
	randomize_resnum_( rval.randomize_resnum_ ),
	bump_check_( rval.bump_check_ )
{}

// destructor
BoltzmannRotamerMover::~BoltzmannRotamerMover()= default;

// clone this object
BoltzmannRotamerMover::MoverOP
BoltzmannRotamerMover::clone() const
{
	return BoltzmannRotamerMover::MoverOP( new protocols::simple_moves::BoltzmannRotamerMover( *this ) );
}

// create this type of object
BoltzmannRotamerMover::MoverOP
BoltzmannRotamerMover::fresh_instance() const
{
	return BoltzmannRotamerMover::MoverOP( new protocols::simple_moves::BoltzmannRotamerMover() );
}

void
BoltzmannRotamerMover::apply( core::pose::Pose & pose )
{
	core::pack::task::PackerTaskOP unedited_task( task(pose)->clone() );
	core::pack::task::PackerTaskOP ptask( task(pose)->clone() );
	if ( show_packer_task_ ) {
		TR << *ptask;
	}

	if ( resnum_ == 0 ) {
		randomize_resnum_ = true;
	}

	if ( randomize_resnum_ ) {
		utility::vector1< core::Size > move_positions;
		for ( core::Size i = 1; i <= ptask->total_residue(); i++ ) {
			if ( ptask->pack_residue(i) || ptask->design_residue(i) ) {
				move_positions.push_back(i);
			}
		}
		core::Size random_index = numeric::random::rg().random_range(1, move_positions.size());
		resnum_ = move_positions[random_index];
	}

	// generate rotamers for position resnum_
	pose.update_residue_neighbors();
	ptask->set_bump_check( bump_check_ );
	ptask->or_include_current( true );
	ptask->temporarily_fix_everything();
	ptask->temporarily_set_pack_residue( resnum_, true );
	scorefxn_->setup_for_packing( pose, ptask->repacking_residues(), ptask->designing_residues() );
	core::conformation::Residue const & res = pose.residue( resnum_ );
	core::pack::rotamer_set::RotamerSetFactory rsf;
	core::pack::rotamer_set::RotamerSetOP rotset = rsf.create_rotamer_set( res );
	rotset->set_resid( resnum_ );
	utility::graph::GraphOP packer_graph = core::pack::create_packer_graph( pose, *scorefxn_, ptask );
	rotset->build_rotamers( pose, *scorefxn_, *ptask, packer_graph );
	utility::vector1< core::PackerEnergy > one_body_energies( rotset->num_rotamers() );
	utility::vector1< utility::vector1< core::PackerEnergy > > two_body_energies( rotset->num_rotamers() );
	utility::vector1< core::Size > packable_neighbors;
	scorefxn_->prepare_rotamers_for_packing( pose, *rotset );

	if ( rotset->num_rotamers() <= 1 ) {
		return;
	}

	rotset->compute_one_and_two_body_energies(
		pose, *scorefxn_, *unedited_task, packer_graph,
		one_body_energies, two_body_energies, packable_neighbors);

	core::Size ligand_neighbor_index = 0;
	for ( core::Size i = 1; i <= packable_neighbors.size(); i++ ) {
		if ( packable_neighbors[i] == ligand_resnum_ ) {
			ligand_neighbor_index = i;
		}
	}

	utility::vector1< utility::vector1<std::pair<core::Size, core::Real> > > boltzmann_factors;
	utility::vector1<core::Real> rotamer_partition_funtions;
	core::Real final_rot_id = 1;
	utility::vector1< core::Real > ligand_bonuses;
	ligand_bonuses.resize( rotset->num_rotamers() );

	// if resnum_ is in a protein, do the following:
	// 1) calculate boltzmann weighted probability for each rotamer
	// 2) use probabilities to select one rotamer for each amino acid
	// 3) calculate boltzmann weighted probabilities for each amino acid
	// 4) use probabilities to select an amino acid
	// 5) replace resnum_ with the selected rotamer / amino acid
	if ( pose.residue(resnum_).is_protein() ) {

		boltzmann_factors.resize( core::chemical::num_canonical_aas );
		rotamer_partition_funtions.resize( core::chemical::num_canonical_aas, 0.0 );

		// iterator over each rotamer
		for ( core::Size i = 1; i <= one_body_energies.size(); i++ ) {
			core::chemical::AA aa_type = (*rotset->rotamer(i)).type().aa();
			core::PackerEnergy init_score = one_body_energies[rotset->id_for_current_rotamer()];
			core::PackerEnergy move_score = one_body_energies[i];

			ligand_bonuses[i] = 0.0;
			if ( ligand_neighbor_index != 0 ) {
				init_score += two_body_energies[rotset->id_for_current_rotamer()][ligand_neighbor_index] * (ligand_weight_ - 1.0);
				move_score += two_body_energies[i][ligand_neighbor_index] * (ligand_weight_ - 1.0);
				ligand_bonuses[i] = two_body_energies[i][ligand_neighbor_index] * (ligand_weight_ - 1.0);
			}

			core::Real boltzmann_factor = std::exp((init_score - move_score) / temperature_);

			if ( boltzmann_factor > 0 ) {
				rotamer_partition_funtions[aa_type] += boltzmann_factor;
				boltzmann_factors[aa_type].push_back(std::make_pair(i, boltzmann_factor));
			}
		}

		utility::vector1<std::pair<core::Size, core::Real> > selected_rotamers;
		core::Real amino_acid_partition_function = 0.0;

		// select rotamer for each amino acid
		core::Size aa_count = 0;
		for ( core::Size aa_type = 1; aa_type <= boltzmann_factors.size(); aa_type++ ) {
			if ( boltzmann_factors[aa_type].size() > 0 ) {
				aa_count += 1;
				core::Real rot = select_rotamer(boltzmann_factors[aa_type], rotamer_partition_funtions[aa_type]);
				selected_rotamers.push_back(std::make_pair(boltzmann_factors[aa_type][rot].first,boltzmann_factors[aa_type][rot].second));
				amino_acid_partition_function += boltzmann_factors[aa_type][rot].second;
			}
		}

		// select an amino acid and get the id of its selected rotamer
		core::Real final_aa = select_rotamer(selected_rotamers, amino_acid_partition_function);
		final_rot_id = selected_rotamers[final_aa].first;
	} else {
		// if resnum_ is not in a protein, do the following:
		// 1) calculate boltzmann weighted probability for each rotamer
		// 2) use probabilities to select one rotamer for each amino acid
		// 3) replace resnum_ with selected rotamer
		boltzmann_factors.resize( 1 );
		rotamer_partition_funtions.resize( 1, 0.0 );

		// iterator over rotamers
		for ( core::Size i = 1; i <= one_body_energies.size(); i++ ) {
			core::PackerEnergy init_score = one_body_energies[rotset->id_for_current_rotamer()];
			core::PackerEnergy move_score = one_body_energies[i];

			ligand_bonuses[i] = 0.0;
			core::Real boltzmann_factor = std::exp((init_score - move_score) / temperature_);
			if ( boltzmann_factor > 0 ) {
				rotamer_partition_funtions[1] += boltzmann_factor;
				boltzmann_factors[1].push_back(std::make_pair(i, boltzmann_factor));
			}
		}

		// select a rotamer and get its id
		core::Real rot = select_rotamer(boltzmann_factors[1], rotamer_partition_funtions[1]);
		final_rot_id = boltzmann_factors[1][rot].first;
	}

	// replace resnum_ with the selected rotamer
	core::conformation::ResidueOP newresidue(  rotset->rotamer( final_rot_id )->clone() );
	pose.replace_residue (resnum_, *newresidue, false );

}

// XRW TEMP std::string
// XRW TEMP BoltzmannRotamerMover::get_name() const {
// XRW TEMP  return "BoltzmannRotamerMover";
// XRW TEMP }

void
BoltzmannRotamerMover::show(std::ostream & output) const
{
	Mover::show(output);
	if ( scorefxn() != nullptr ) {
		output << "Score function: " << scorefxn()->get_name() << std::endl;
	} else { output << "Score function: none" << std::endl; }
}

// helpers

/// @brief select a rotamer based on Boltzmann weighted probabilities
core::Size
BoltzmannRotamerMover::select_rotamer(
	utility::vector1<std::pair<core::Size, core::Real> > const & boltzmann_factors,
	core::Real const & partition_function)
{
	utility::vector1<std::pair<core::Size, core::Real > > boltzmann_probs;
	for ( core::Size rot = 1; rot <= boltzmann_factors.size(); rot++ ) {
		core::Real boltzmann_probability = boltzmann_factors[rot].second / partition_function;
		if ( bias_sampling_ ) {
			boltzmann_probs.push_back(std::make_pair(rot, boltzmann_probability));
		} else {
			boltzmann_probs.push_back(std::make_pair(rot, core::Real(1.0) / core::Real(boltzmann_factors.size())));
		}
	}
	std::sort(boltzmann_probs.begin(), boltzmann_probs.end(), compare_values);
	core::Real rotnum = 0;
	core::Real random_prob = numeric::random::uniform();
	while ( random_prob > 0 ) {
		rotnum++;
		random_prob -= boltzmann_probs[rotnum].second;
		if ( rotnum == boltzmann_probs.size() ) break;
	}
	return boltzmann_probs[rotnum].first;
}

// setters
void BoltzmannRotamerMover::set_score_function( core::scoring::ScoreFunctionCOP sf ) { scorefxn_ = sf; }
void BoltzmannRotamerMover::set_task_factory( core::pack::task::TaskFactoryCOP tf ) { factory_ = tf; }
void BoltzmannRotamerMover::set_resnum( core::Size resnum ) { resnum_ = resnum; }
void BoltzmannRotamerMover::set_ligand_resnum( core::Size ligand_resnum ) { ligand_resnum_ = ligand_resnum; }
void BoltzmannRotamerMover::set_ligand_weight( core::Real ligand_weight ) { ligand_weight_ = ligand_weight; }
void BoltzmannRotamerMover::set_temperature( core::Real temperature ) { temperature_ = temperature; }
void BoltzmannRotamerMover::set_bias_sampling( bool bias_sampling ) { bias_sampling_ = bias_sampling; }
void BoltzmannRotamerMover::set_randomize_resnum( bool randomize_resnum ) { randomize_resnum_ = randomize_resnum; }
void BoltzmannRotamerMover::set_bump_check( bool bump_check ) { bump_check_ = bump_check; }

// getters
core::Size
BoltzmannRotamerMover::get_resnum() const {
	return resnum_;
}
core::Size
BoltzmannRotamerMover::get_ligand_resnum() const {
	return ligand_resnum_;
}
core::Real
BoltzmannRotamerMover::get_ligand_weight() const {
	return ligand_weight_;
}
core::Real
BoltzmannRotamerMover::get_temperature() const {
	return temperature_;
}
bool
BoltzmannRotamerMover::get_bias_sampling() const {
	return bias_sampling_;
}
bool
BoltzmannRotamerMover::get_randomize_resnum() const {
	return randomize_resnum_;
}
bool
BoltzmannRotamerMover::get_bump_check() const {
	return bump_check_;
}


/// @brief read access for derived classes
BoltzmannRotamerMover::ScoreFunctionCOP
BoltzmannRotamerMover::scorefxn() const
{
	return scorefxn_;
}

/// @brief read access for derived classes
BoltzmannRotamerMover::PackerTaskCOP
BoltzmannRotamerMover::task( core::pose::Pose const & pose ) const
{
	//if we have a factory, generate and return a new task
	if ( factory_ ) return factory_->create_task_and_apply_taskoperations( pose );
	//else runtime_assert( task_is_valid( pose ) );

	//else return the unsafe one
	return task_;
}

/// @brief parse xml
void
BoltzmannRotamerMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & )
{
	using core::scoring::ScoreFunction;
	using core::pack::task::operation::TaskOperation;
	using core::pack::task::TaskFactoryOP;
	using core::pack::task::TaskFactory;

	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	show_packer_task_ = tag->getOption<bool>( "show_packer_task", 0 );
	set_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
}

std::string BoltzmannRotamerMover::get_name() const {
	return mover_name();
}

std::string BoltzmannRotamerMover::mover_name() {
	return "BoltzmannRotamerMover";
}

void BoltzmannRotamerMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	rosetta_scripts::attributes_for_parse_task_operations(attlist);
	rosetta_scripts::attributes_for_parse_score_function(attlist);
	attlist + XMLSchemaAttribute( "show_packer_task", xsct_rosetta_bool, "show the PackerTask to be used at the beginning of apply");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Part of CoupledMoves. Replaces a single rotamer based on the Boltzmann probabilities of its rotamers", attlist );
}

std::string BoltzmannRotamerMoverCreator::keyname() const {
	return BoltzmannRotamerMover::mover_name();
}

protocols::moves::MoverOP
BoltzmannRotamerMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new BoltzmannRotamerMover );
}

void BoltzmannRotamerMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BoltzmannRotamerMover::provide_xml_schema( xsd );
}


std::ostream &operator<< (std::ostream &os, BoltzmannRotamerMover const &mover)
{
	mover.show(os);
	return os;
}

} // moves
} // protocols
