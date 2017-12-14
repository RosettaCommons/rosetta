// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/sidechain_moves/SidechainMover.cc
/// @brief implementation of SidechainMover class and functions
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMoverCreator.hh>

// Procols Headers
#include <basic/datacache/DataMap.hh>

// Core Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/id/TorsionID.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <basic/prof.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

// Utility
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>

// C++ Headers
#include <sstream>
#include <utility/fixedsizearray1.hh>

#include <core/id/TorsionID_Range.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <protocols/rosetta_scripts/util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using namespace core;
using namespace core::pose;

static basic::Tracer TR( "protocols.simple_moves.sidechain_moves.SidechainMover" );

namespace protocols {
namespace simple_moves {
namespace sidechain_moves {

// XRW TEMP std::string
// XRW TEMP SidechainMoverCreator::keyname() const {
// XRW TEMP  return SidechainMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SidechainMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SidechainMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SidechainMover::mover_name() {
// XRW TEMP  return "Sidechain";
// XRW TEMP }

SidechainMover::SidechainMover():
	rotamer_library_( * core::pack::dunbrack::RotamerLibrary::get_instance() ),
	pose_(/* 0 */),
	prob_uniform_(0.1),
	prob_withinrot_(0.0),
	prob_random_pert_to_current_(0.0),
	preserve_detailed_balance_(false),
	accept_according_to_dunbrack_(true),
	sample_rotwells_unif_(false),
	change_chi_without_replacing_residue_(false),
	next_resnum_(0),
	last_proposal_density_ratio_(1),
	task_initialized_(false),
	scratch_( core::pack::dunbrack::RotamerLibraryScratchSpaceOP( new core::pack::dunbrack::RotamerLibraryScratchSpace ) ),
	temperature0_(0.56),
	sampling_temperature_(0.56)
{}

SidechainMover::SidechainMover(
	core::pack::dunbrack::RotamerLibrary const & rotamer_library
):
	rotamer_library_(rotamer_library),
	pose_(/* 0 */),
	prob_uniform_(0.1),
	prob_withinrot_(0.0),
	prob_random_pert_to_current_(0.0),
	preserve_detailed_balance_(false),
	accept_according_to_dunbrack_(true),
	sample_rotwells_unif_(false),
	change_chi_without_replacing_residue_(false),
	next_resnum_(0),
	last_proposal_density_ratio_(1),
	task_initialized_(false),
	scratch_( core::pack::dunbrack::RotamerLibraryScratchSpaceOP( new core::pack::dunbrack::RotamerLibraryScratchSpace ) ),
	temperature0_(0.56),
	sampling_temperature_(0.56)
{}

SidechainMover::SidechainMover(
	SidechainMover const & mover
) :
	//utility::pointer::ReferenceCount(),
	protocols::canonical_sampling::ThermodynamicMover(mover),
	rotamer_library_(mover.rotamer_library_),
	packed_residues_(mover.packed_residues_),
	residue_packed_(mover.residue_packed_),
	prob_uniform_(mover.prob_uniform_),
	prob_withinrot_(mover.prob_withinrot_),
	prob_random_pert_to_current_(mover.prob_random_pert_to_current_),
	preserve_detailed_balance_(mover.preserve_detailed_balance_),
	accept_according_to_dunbrack_(mover.accept_according_to_dunbrack_),
	sample_rotwells_unif_(mover.sample_rotwells_unif_),
	change_chi_without_replacing_residue_(mover.change_chi_without_replacing_residue_),
	next_resnum_(mover.next_resnum_),
	last_chi_angles_(mover.last_chi_angles_),
	last_nchi_(mover.last_nchi_),
	last_mutation_(mover.last_mutation_),
	last_uniform_(mover.last_uniform_),
	last_withinrot_(mover.last_withinrot_),
	last_pertrot_(mover.last_pertrot_),
	last_proposal_density_ratio_(mover.last_proposal_density_ratio_),
	task_initialized_(mover.task_initialized_),
	temperature0_(0.56),
	sampling_temperature_(0.56)
{
	if ( mover.task_factory_ ) task_factory_ = core::pack::task::TaskFactoryCOP( core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory(*mover.task_factory_) ) );
	if ( mover.task_ ) task_ = mover.task_->clone();
	if ( mover.pose_ ) pose_ = core::pose::PoseOP( new core::pose::Pose(*mover.pose_) );
	if ( mover.scratch_ ) scratch_ = core::pack::dunbrack::RotamerLibraryScratchSpaceOP( new core::pack::dunbrack::RotamerLibraryScratchSpace(*mover.scratch_) );
}

SidechainMover::~SidechainMover() = default;

protocols::moves::MoverOP
SidechainMover::clone() const
{
	return protocols::moves::MoverOP( new protocols::simple_moves::sidechain_moves::SidechainMover(*this) );
}

void
SidechainMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
)
{
	core::pack::task::TaskFactoryOP new_task_factory( new core::pack::task::TaskFactory );

	if ( tag->hasOption("task_operations") ) {

		std::string const t_o_val( tag->getOption<std::string>("task_operations") );
		using StringVec = utility::vector1<std::string>;
		StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );
		for ( auto const & t_o_key : t_o_keys ) {
			if ( data.has( "task_operations", t_o_key ) ) {
				new_task_factory->push_back( data.get_ptr< core::pack::task::operation::TaskOperation >( "task_operations", t_o_key ) );
			} else {
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "TaskOperation " + t_o_key + " not found in basic::datacache::DataMap.");
			}
		}

	} else {
		using core::pack::task::operation::TaskOperationCOP;
		new_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );
	}

	task_factory_ = new_task_factory;

	set_preserve_detailed_balance( tag->getOption<bool>( "preserve_detailed_balance", preserve_detailed_balance() ) );
	set_prob_uniform( tag->getOption<core::Real>( "prob_uniform", prob_uniform() ) );
	set_prob_withinrot( tag->getOption<core::Real>( "prob_withinrot", prob_withinrot() ) );
	set_prob_random_pert_current( tag->getOption<core::Real>( "prob_random_pert_current", prob_random_pert_current() ) );
	set_change_chi_without_replacing_residue( tag->getOption<core::Real>( "change_chi_without_replacing_residue", change_chi_without_replacing_residue() ) );
}

/// @details
/// Check to make sure that a packer task exists and matches the numer of residues in
/// the given pose. If that isn't the case, create a new one with the task factory.
/// Exits with an error if no task factory exists.
void
SidechainMover::init_task(
	core::pose::Pose const & pose
)
{
	// check to see if a valid task already exists
	if ( !task_ || task_->total_residue() != pose.size() ) {
		// if not, create one using the task factory
		if ( task_factory_ ) {
			set_task(task_factory_->create_task_and_apply_taskoperations( pose ));
		} else {
			utility_exit_with_message("Cannot create task because no task factory is set");
		}
	}
	/// apl -- is it reasonable to first check the task to see if infact we are designing?
	/// copying a pose is expensive -- try to avoid it if possible.
	if ( !pose_ ) { //pose not initialized
		pose_ = core::pose::PoseOP( new core::pose::Pose( pose ) ); //need this in case we're designing as well. ek 4/28/10
	}


}

void
SidechainMover::initialize_simulation(
	core::pose::Pose & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const&, /*metropolis_hastings_mover*/
	core::Size /*cycle*/
)
{
	idealize_sidechains(pose);
}

void
SidechainMover::perturb_chi(
	numeric::random::RandomGenerator & Rand,
	core::Real max_deviation,
	utility::vector1< core::Real > & old_chi,
	utility::vector1< core::Real > & new_chi
){
	new_chi.resize( old_chi.size() );
	core::Real rand;
	for ( unsigned int chi_i = 1; chi_i <= old_chi.size(); chi_i++ ) {
		rand = Rand.uniform();
		//int sign = Rand.random_range(0,1);
		//if(sign == 0){
		// sign = -1;
		//}
		//new_chi[ chi_i ] = basic::periodic_range( ((rand * max_deviation * sign) + old_chi[ chi_i ]) , 360.0 );
		new_chi[ chi_i ] = basic::periodic_range( ( (2.0*rand-1.0)*max_deviation + old_chi[ chi_i ]) , 360.0 );
	}
}

bool
SidechainMover::dunbrack_accept(
	numeric::random::RandomGenerator & Rand,
	core::conformation::Residue & res,
	utility::vector1<core::Real> const & previous_chi_angles,
	utility::vector1<core::Real> const & new_chi_angles
)
{
	core::Real rand = Rand.uniform();
	//conformation::Residue copy(res); // APL: EXPENSIVE!
	res.set_all_chi( previous_chi_angles );
	core::Real prev_chi_angle_score = rotamer_library_.rotamer_energy( res, *scratch_);
	res.set_all_chi( new_chi_angles );

	core::Real new_chi_angle_score = rotamer_library_.rotamer_energy( res, *scratch_);
	if ( new_chi_angle_score > prev_chi_angle_score ) {
		//Real const boltz_factor = (prev_chi_angle_score - new_chi_angle_score)/1.0;
		Real const boltz_factor = (prev_chi_angle_score - new_chi_angle_score)/sampling_temperature_*temperature0_;
		//Real const probability = std::exp(std::max(Real(-40.0),boltz_factor));
		Real const probability = std::exp(boltz_factor);
		if ( rand >= probability ) {
			return false;
		}
	}
	return true;
}


core::conformation::ResidueOP
SidechainMover::make_move( core::conformation::ResidueOP input_residue )
{

	using numeric::conversions::degrees;
	using numeric::conversions::radians;

	// reset the last move tracking information
	last_mutation_ = false;
	last_uniform_ = false;
	last_withinrot_ = false;
	last_pertrot_ = false;

	core::chemical::ResidueType const & previous_residue_type( input_residue->type() );
	utility::vector1<core::Real> previous_chi_angles( input_residue->chi() );

	core::Size resnum = input_residue->seqpos();
	core::pack::task::ResidueLevelTask const & residue_task(task_->residue_task(resnum));
	core::pack::task::ResidueLevelTask::ResidueTypeCOPList const & residue_types(residue_task.allowed_residue_types());

	// select a residue type
	core::chemical::ResidueTypeCOP residue_type;
	do {
		Size const restypenum(numeric::random::rg().random_range(1, residue_types.size()));
		auto iter(residue_types.begin());
		for ( Size i = 1; i < restypenum; ++i ) ++iter;
		residue_type = *iter;
	}
	while (residue_type->aa() == core::chemical::aa_pro); // temporary hack to exlcude sampling of proline sidechains

	last_nchi_ = residue_type->nchi();

	///code note: stl uses lazy resizes; if a vector doesn't have to allocate more space, it doesn't, it just
	/// updates its end-of-memory pointer.  By always performing the resize, the vector size() will reflect
	/// the actual number of chi.
	///if (last_chi_angles_.size() < last_nchi_) last_chi_angles_.resize(last_nchi_);
	last_chi_angles_.resize( last_nchi_ );

	if ( TR.visible( basic::t_debug ) ) {
		TR.Debug << input_residue->name() << " " << resnum << " -> " << residue_type->name() << std::endl;
	}

	if ( residue_type->aa() != (input_residue->type()).aa() ) last_mutation_ = true;

	core::pack::rotamers::SingleResidueRotamerLibraryCOP residue_rotamer_library(
		core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get(*residue_type)
	);
	core::pack::dunbrack::SingleResidueDunbrackLibraryCOP residue_dunbrack_library(
		utility::pointer::dynamic_pointer_cast< core::pack::dunbrack::SingleResidueDunbrackLibrary const >(residue_rotamer_library)
	);

	/// last_chi_angles_ holds degrees for a short period;
	/// Apply uniform sampling to rotamers with a small probability, or with probability 1 if
	/// the residue's library is not a dunbrack library. (Is this what we would want for pSer?).
	/// At the moment, the ligand rotamer library does not respond to the assign_random_rotamer_with_bias
	/// method, so this is an ok check.

	core::Real move_type_prob( numeric::random::rg().uniform() );
	if ( move_type_prob > prob_uniform_  && residue_dunbrack_library ) { //set up rotamer library stuff
		//p_within_rot or p_jump_rots or p_random_pert_to_current
		core::Real const phi( residue_dunbrack_library->get_phi_from_rsd( *input_residue ) );
		core::Real const psi( residue_dunbrack_library->get_psi_from_rsd( *input_residue ) );

		/* //add in option to sample rotamers above a certain threshold uniformly
		utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > most_probable_rotamers;
		Size itr_i = 1;
		core::Real const probability_threshold = 0.01;
		debug_assert( rotamer_sample_data.size() > 0 );
		if ( sample_rotwells_unif_  ){ //uniformly samples rotamer wells
		while((itr_i <= rotamer_sample_data.size())  &&
		(rotamer_sample_data[ itr_i ].probability() > probability_threshold)){
		most_probable_rotamers.push_back( rotamer_sample_data[ itr_i ] );
		itr_i++;
		}
		debug_assert( most_probable_rotamers.size() > 0 );
		rotamer_sample_data = most_probable_rotamers; //replace all rots with ones that pass a certain threshold (1%)
		}*/

		if ( move_type_prob > (prob_uniform_ + prob_withinrot_ + prob_random_pert_to_current_ ) || last_mutation_ ) {
			utility::fixedsizearray1< Real, 5 > bbs;
			bbs[ 1 ] = phi; bbs[ 2 ] = psi;
			utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > rotamer_sample_data(
				residue_dunbrack_library->get_all_rotamer_samples(/*phi, psi*/bbs)
			);
			make_rotwell_jump( rotamer_sample_data ); // returns last_chi_angles_
		} else if ( move_type_prob > (prob_uniform_ + prob_withinrot_) ) {
			//ek added in new option to perturb current chi and evaluate according to dunbrack
			preturb_rot_and_dunbrack_eval( input_residue ); // ?
		} else if ( move_type_prob > prob_uniform_ ) {
			utility::fixedsizearray1< Real, 5 > bbs;
			bbs[ 1 ] = phi; bbs[ 2 ] = psi;
			utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > rotamer_sample_data(
				residue_dunbrack_library->get_all_rotamer_samples(/*phi, psi*/bbs)
			);
			perturb_rot_within_well( rotamer_sample_data, previous_chi_angles );
		}
	} else {
		if ( TR.visible( basic::t_debug ) ) {
			TR.Debug << "making a uniform jump " << std::endl;
		}
		for ( Size ii = 1; ii <= last_nchi_; ++ii ) {
			last_chi_angles_[ ii ] = numeric::random::rg().uniform() * 360.0 - 180.0;
		}
		last_uniform_ = true;
	}
	//at the end, have the new angles assigned to last_chi_angles_

	core::Real proposal_density_reverse(1);
	core::Real proposal_density_forward(1);

	if ( preserve_detailed_balance_ ) {
		proposal_density_reverse = proposal_density(*input_residue, resnum, *residue_type, last_chi_angles_);
	}


	if ( TR.visible( basic::t_debug ) ) {
		TR.Debug << "residue angles were previously: ";

		for ( Size i = 1; i <= input_residue->nchi(); ++i ) {
			TR.Debug << " " << input_residue->chi( i );;
		}
		TR.Debug << std::endl;
	}

	/// set the chi
	for ( Size i = 1; i <= input_residue->nchi(); ++i ) {
		input_residue->set_chi( i, last_chi_angles_[ i ] );
	}

	if ( TR.visible( basic::t_debug ) ) {
		TR.Debug << "Set residue chi angles to:";
		for ( Size i = 1; i <= input_residue->nchi(); ++i ) {
			TR.Debug << " " << last_chi_angles_[ i ];
			TR.Debug << " " << input_residue->chi( i );
		}
		TR.Debug << std::endl;
	}

	// swap in a new residue if necessary
	core::conformation::ResidueOP previous_residue = input_residue;
	core::conformation::ResidueOP new_residue;
	if ( residue_type.get() != & previous_residue->type() ) { // apl, note: pointer comparison is faster
		if ( TR.visible( basic::t_debug ) ) {
			TR.Debug << "previous residue " << previous_residue->type().name() << " doesn't match new res-name " << residue_type->name() << std::endl;
		}
		new_residue =
			core::conformation::ResidueFactory::create_residue(*residue_type,
			*input_residue,
			pose_->conformation(),
			residue_task.preserve_c_beta());
		if ( TR.visible( basic::t_debug ) ) {
			TR.Debug << "type of new residue created " << new_residue->type().name() << std::endl;
		}
		for ( Size ii = 1; ii <= residue_type->nchi(); ii++ ) {
			if ( TR.visible( basic::t_debug ) ) {
				TR.Debug << " setting chi of new residue " << last_chi_angles_[ ii ] << std::endl;
			}
			new_residue->set_chi( ii, numeric::principal_angle_degrees(last_chi_angles_[ ii ] ));
			if ( TR.visible( basic::t_debug ) ) {
				TR.Debug << " input residue chi is now: " << input_residue->chi( ii ) << std::endl;
			}
		}
	} else {
		new_residue = input_residue;
	}
	//core::conformation::ResidueOP new_residue(input_residue);
	//pose_.replace_residue(resnum, *new_residue, false); //moved to outside apply

	/// last_chi_angles_ holds radians after this code executes;
	for ( Size ii = 1; ii <= residue_type->nchi(); ++ii ) {
		last_chi_angles_[ ii ] = numeric::principal_angle_degrees(last_chi_angles_[ ii ]);
		//pose.set_chi( ii, resnum, last_chi_angles_[ ii ] );
		last_chi_angles_[ ii ] = radians( last_chi_angles_[ ii ] );
	}


	if ( preserve_detailed_balance_ ) {
		proposal_density_forward = proposal_density(*input_residue, resnum, previous_residue_type, previous_chi_angles);
	}
	last_proposal_density_ratio_ = proposal_density_reverse / proposal_density_forward;
	if ( TR.visible( basic::t_debug ) ) {
		TR.Debug << "last proposal density ratio set to: "  << last_proposal_density_ratio_ << std::endl;
	}


	update_type();
	/// comments in this section were distracting me .. they have been moved to the end of the file.
	return new_residue;

}

void
SidechainMover::make_rotwell_jump(
	utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > const & rotamer_sample_data
)
{
	Size rotnum = 0;
	if ( TR.visible( basic::t_debug ) ) {
		TR.Debug << "making a rot-well jump " << std::endl;
	}
	// choose a random rotamer according to its probability
	core::Real random_prob = numeric::random::rg().uniform();

	//yuan, 01/06/10
	//modify the rotprob, if want to sample at different temperature
	utility::vector1< core::Real > rot_prob_list(rotamer_sample_data.size());
	core::Real sum_prob=0.0;
	for ( rotnum=rotamer_sample_data.size(); rotnum>0; rotnum-- ) {
		rot_prob_list[rotnum] = rotamer_sample_data[rotnum].probability();
		//modified by temperature
		rot_prob_list[rotnum]=std::pow(rot_prob_list[rotnum],temperature0_/sampling_temperature_);
		sum_prob+=rot_prob_list[rotnum];
	}

	runtime_assert(sum_prob>0.0);
	for ( rotnum=rotamer_sample_data.size(); rotnum>0; rotnum-- ) {rot_prob_list[rotnum]/=sum_prob;}

	rotnum = 0;
	/// go through rotamers in decreasing order by probability and stop when the
	while ( random_prob > 0 ) {
		if ( sample_rotwells_unif_ ) {
			random_prob -= (1.0/rotamer_sample_data.size()); //pick any rotamer with uniform probability
			rotnum++;
		} else {
			//random_prob -= rotamer_sample_data[++rotnum].probability();
			random_prob -= rot_prob_list[++rotnum];
		}
		// loop condition might end up satisfied even if we've walked through all possible rotamers
		// if the chosen random number was nearly 1
		// (and interpolation introduced a tiny bit of numerical noise).
		if ( rotnum == rotamer_sample_data.size() ) break;
	}
	//rotamer_sample_data[rotnum].assign_random_chi(last_chi_angles_,numeric::random::rg());
	rotamer_sample_data[rotnum].assign_random_chi(last_chi_angles_, numeric::random::rg(), std::sqrt(sampling_temperature_/temperature0_));
	for ( unsigned int ii = 1; ii <= last_chi_angles_.size(); ii++ ) {
		last_chi_angles_[ii] = basic::periodic_range( last_chi_angles_[ii], 360.0 );
	}

}

void
SidechainMover::preturb_rot_and_dunbrack_eval( core::conformation::ResidueOP input_residue )
{
	if ( TR.visible( basic::t_debug ) ) {
		TR.Debug << "making a perterb-rot (ek) " << std::endl;
	}
	core::Real const max_deviation = 10.0; //arbitrary_value!
	utility::vector1< core::Real > new_chi_angles;
	utility::vector1< core::Real > previous_chi_angles = input_residue->chi();
	if ( accept_according_to_dunbrack_ ) {
		//compute probability of acceptance based on dunbrack rotamer distribution
		/// apl -- unused -- core::Real rand = 0;
		/// apl -- unused -- core::Real proposed_move_dunbrack_probability = -1;
		/// apl -- unused -- core::Real before_pert_prob = 0;
		do {
			perturb_chi( numeric::random::rg(), max_deviation, previous_chi_angles, new_chi_angles );
			// new_chi_angles.resize( previous_chi_angles.size() );
			// for(unsigned int chi_i = 1; chi_i <= previous_chi_angles.size(); chi_i++)
			//{
			//  new_chi_angles[ chi_i ] = basic::periodic_range( ( (2.0*RG.uniform()-1.0)*max_deviation + previous_chi_angles[ chi_i ]) , 360.0 );
			// }

		}
		while( !dunbrack_accept( numeric::random::rg(), *input_residue, previous_chi_angles, new_chi_angles ) );
		if ( preserve_detailed_balance_ ) {
			TR.Error << "You cannot specify accept_according_to_dunbrack_ and preserve_detailed_balance_ both as true!" << std::endl;
			TR.Error << "Recommend rerunning with 'random perturb current' frequency set to zero." << std::endl;
			utility_exit();
		}
	}
	last_pertrot_ = true;
	last_chi_angles_ = new_chi_angles;

}

void
SidechainMover::perturb_rot_within_well(
	utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > const & rotamer_sample_data,
	utility::vector1<core::Real> const & previous_chi_angles
)
{
	Size rotnum = 0;
	//ek added in new option to perturb current chi and evaluate according to dunbrack 3/30/10
	if ( TR.visible( basic::t_debug ) ) {
		TR.Debug << "making a within-rot-well jump " << std::endl;
	}
	// find the rotamer that has the highest probability of proposing the previous chi angles
	core::Real max_rot_prob = 0;
	for ( core::Size ii = 1; ii <= rotamer_sample_data.size(); ++ii ) {
		core::Real candidate_rot_prob( rotamer_sample_data[ii].chi_probability(previous_chi_angles) );
		if ( candidate_rot_prob > max_rot_prob ) {
			rotnum = ii;
			max_rot_prob = candidate_rot_prob;
		}
	}
	last_withinrot_ = true;
	rotamer_sample_data[rotnum].assign_random_chi(last_chi_angles_,numeric::random::rg());
	for ( unsigned int ii = 1; ii <= last_chi_angles_.size(); ii++ ) {
		last_chi_angles_[ii] = basic::periodic_range(last_chi_angles_[ii],360.0);
	}
	//residue_rotamer_library->assign_random_rotamer_with_bias(
	// pose.residue( resnum ), *scratch_, numeric::random::rg(),
	// last_chi_angles_, true );

}


bool SidechainMover::task_initialized(){
	return task_initialized_;
}

/// @details
void
SidechainMover::apply(
	Pose & pose
)
{
	using numeric::conversions::degrees;
	using numeric::conversions::radians;

	PROF_START( basic::APPLY_SC_MOVE );

	init_task(pose);

	// select a residue to change
	runtime_assert(packed_residues_.size() > 0);

	Size resnum;
	if ( next_resnum_ ) {
		resnum = next_resnum_;
		next_resnum_ = 0;
	} else {
		resnum = packed_residues_[ numeric::random::rg().random_range(1, packed_residues_.size()) ];
	}
	conformation::ResidueOP newresidue( new conformation::Residue( pose.residue( resnum ) ) );
	if ( TR.visible( basic::t_debug ) ) {
		TR.Debug << "old residue chi is : ";
		for ( Size i = 1; i <= newresidue->nchi(); ++i ) TR.Debug << " " << newresidue->chi(i);
		TR.Debug << std::endl;
	}

	conformation::ResidueOP final = make_move( newresidue );
	if ( TR.visible( basic::t_debug ) ) {
		TR.Debug << "new residue chi is : ";
		for ( Size i = 1; i <= pose.residue(resnum).type().nchi(); ++i ) TR.Debug << " " << newresidue->chi(i);
		TR.Debug << std::endl;
	}

	//perform replace_residue by default, or if there was a mutation
	if ( !change_chi_without_replacing_residue() || last_mutation_ ) {
		pose.replace_residue(resnum, *final, true);
	} else { //optionally, perform in-place chi rotations instead of a residue replace.  This causes proper kinematic dependency if the atomtree has atoms that are children of the moving residue.
		core::Size const nchi(final->nchi());
		for ( core::Size i(1); i<=nchi; ++i ) {
			pose.set_chi(i, resnum, final->chi(i));
		}
	}

	if ( TR.visible( basic::t_debug ) ) {
		TR.Debug << "pose replaced " << resnum << " " << pose.residue(resnum).type().name() << " has phi-psi " << pose.phi(resnum) << " " << pose.psi(resnum) << " and chi: ";

		for ( Size i = 1; i <= pose.residue(resnum).type().nchi(); ++i ) TR.Debug << " " << (pose.residue(resnum).chi(i));
		TR.Debug << std::endl;
	}
	PROF_STOP( basic::APPLY_SC_MOVE );
}

// XRW TEMP std::string
// XRW TEMP SidechainMover::get_name() const {
// XRW TEMP  return "SidechainMover";
// XRW TEMP }


core::Real
SidechainMover::proposal_density(
	core::conformation::Residue const & proposed_residue,
	core::Size const proposed_resnum,
	core::chemical::ResidueType const & initial_residue_type,
	utility::vector1<core::Real> const & initial_chi_angles
) const
{
	utility::vector1<core::Real> const & proposed_chi_angles( proposed_residue.chi() );

	core::Real density(0);

	core::pack::rotamers::SingleResidueRotamerLibraryCOP residue_rotamer_library(
		core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get(proposed_residue.type())
	);
	core::pack::dunbrack::SingleResidueDunbrackLibraryCOP residue_dunbrack_library(
		utility::pointer::dynamic_pointer_cast< core::pack::dunbrack::SingleResidueDunbrackLibrary const >(residue_rotamer_library)
	);

	if ( residue_dunbrack_library ) {

		core::Real const phi( residue_dunbrack_library->get_phi_from_rsd( proposed_residue ) );
		core::Real const psi( residue_dunbrack_library->get_psi_from_rsd( proposed_residue ) );

		utility::fixedsizearray1< Real, 5 > bbs;
		bbs[ 1 ] = phi; bbs[ 2 ] = psi;
		utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > rotamer_sample_data(
			residue_dunbrack_library->get_all_rotamer_samples(/*phi, psi*/bbs)
		);

		utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > most_probable_rotamers;
		core::Real probability_threshold = 0.01;
		if ( sample_rotwells_unif_ ) {
			core::Size itr_i = 1;
			while ( itr_i <= rotamer_sample_data.size() &&
					rotamer_sample_data[ itr_i ].probability() > probability_threshold ) {
				most_probable_rotamers.push_back( rotamer_sample_data[ itr_i ] );
				itr_i++;
			}
			rotamer_sample_data = most_probable_rotamers; //replace all rots with most probable ones
		}

		core::Real rot_density(0);

		core::Real max_proposed_rot_prob(0);
		core::Real max_initial_rot_prob(0);

		bool calc_withinrot( prob_withinrot_ && proposed_residue.type().aa() == initial_residue_type.aa() );

		for ( core::Size ii = 1; ii <= rotamer_sample_data.size(); ++ii ) {

			core::Real const proposed_rot_prob( rotamer_sample_data[ii].chi_probability(proposed_chi_angles) );

			if ( sample_rotwells_unif_ ) {
				rot_density += ( 1.0 / rotamer_sample_data.size() ) * proposed_rot_prob;
			} else {
				//rot_density += rotamer_sample_data[ii].probability() * proposed_rot_prob;
				rot_density += exp(log(rotamer_sample_data[ii].probability()) + log(proposed_rot_prob));
			}

			if ( calc_withinrot ) {

				core::Real const initial_rot_prob( rotamer_sample_data[ii].chi_probability(initial_chi_angles) );

				if ( initial_rot_prob > max_initial_rot_prob ) {
					max_initial_rot_prob = initial_rot_prob;
					max_proposed_rot_prob = proposed_rot_prob;
				}
			}
		}

		rot_density *= 1 - prob_uniform_ - prob_withinrot_;

		density += rot_density;

		if ( calc_withinrot ) {

			core::Real const withinrot_density( prob_withinrot_ * max_proposed_rot_prob );

			density += withinrot_density;
		}
	}

	if ( prob_uniform_ || !residue_dunbrack_library ) {

		core::Real uniform_density(1);

		for ( core::Size ii = 1; ii <= proposed_residue.nchi(); ++ii ) {
			uniform_density *= 1.0/360.0;
		}

		if ( residue_dunbrack_library ) uniform_density *= prob_uniform_;

		density += uniform_density;
	}

	// Probability of selecting the given ResidueType
	// This isn't needed now, but should be changed if the way the ResidueTypes
	// are sampled changes.
	core::pack::task::ResidueLevelTask const & residue_task(task_->residue_task(proposed_resnum));
	density *= 1.0/residue_task.allowed_residue_types().size();
	return density;
}

/// @details
void
SidechainMover::test_move(
	Pose &
)
{

}

/// @details
/// all sidechains that might be changed are replaced with ideal coordinates that have
/// the original chi angles
void
SidechainMover::idealize_sidechains(
	core::pose::Pose & pose
)
{
	utility::vector1<Real> chi;

	init_task(pose);

	for ( Size i = 1; i <= packed_residues_.size(); ++i ) {

		Size const resnum(packed_residues_[i]);

		// get the residue type and residue level task
		chemical::ResidueType const & residue_type(pose.residue_type(resnum));
		pack::task::ResidueLevelTask const & residue_task(task_->residue_task(resnum));

		// disable proline idealization for now
		if ( residue_type.aa() == chemical::aa_pro ) continue;

		// save original chi angles
		if ( residue_type.nchi() > chi.size() ) chi.resize(residue_type.nchi());
		for ( Size chinum = 1; chinum <= residue_type.nchi(); ++chinum ) {
			chi[chinum] = pose.chi(chinum, resnum);
		}

		// create a new residue and replace the old one with it
		conformation::ResidueOP new_residue(
			conformation::ResidueFactory::create_residue(residue_type, pose.residue(resnum), pose.conformation(),
			residue_task.preserve_c_beta())
		);
		pose.replace_residue(resnum, *new_residue, false);

		// put original chi angles back
		for ( Size chinum = 1; chinum <= residue_type.nchi(); ++chinum ) {
			pose.set_chi(chinum, resnum, chi[chinum]);
		}
	}
}

core::pack::dunbrack::RotamerLibrary const &
SidechainMover::rotamer_library() const
{
	return rotamer_library_;
}

core::pack::task::TaskFactoryCOP
SidechainMover::task_factory() const
{
	return task_factory_;
}

void
SidechainMover::set_task_factory(
	core::pack::task::TaskFactoryCOP task_factory
)
{
	task_factory_ = task_factory;
}

core::pack::task::PackerTaskCOP
SidechainMover::task() const
{
	return task_;
}

void
SidechainMover::set_task(
	core::pack::task::PackerTaskCOP task
)
{
	using core::pack::task::ResidueLevelTask;

	task_ = task;

	// update the list of residues being packed
	packed_residues_.clear();
	residue_packed_.resize(task_->total_residue());
	for ( Size i = 1; i <= task_->total_residue(); ++i ) {

		// loop over all possible residue types and check if any have chi angles
		residue_packed_[i] = false;
		core::pack::task::ResidueLevelTask const & residue_task(task->residue_task(i));
		for ( auto iter(residue_task.allowed_residue_types_begin());
				iter != residue_task.allowed_residue_types_end(); ++iter ) {

			// temporarily exclude proline residues from side chain sampling
			if ( (*iter)->nchi() > 0 && (*iter)->aa() != core::chemical::aa_pro ) {
				packed_residues_.push_back(i);
				residue_packed_[i] = true;
				break;
			}
		}
	}
}

core::Real
SidechainMover::prob_uniform() const
{
	return prob_uniform_;
}

void
SidechainMover::set_prob_uniform(
	core::Real prob_uniform
)
{
	prob_uniform_ = prob_uniform;
}

core::Real
SidechainMover::prob_withinrot() const
{
	return prob_withinrot_;
}

void
SidechainMover::set_prob_withinrot(
	core::Real prob_withinrot
)
{
	prob_withinrot_ = prob_withinrot;
}

core::Real
SidechainMover::prob_random_pert_current() const{
	return prob_random_pert_to_current_;
}

void
SidechainMover::set_prob_random_pert_current(
	core::Real prob_pert
)
{
	prob_random_pert_to_current_ = prob_pert;
}

bool
SidechainMover::preserve_detailed_balance() const
{
	return preserve_detailed_balance_;
}

void
SidechainMover::set_preserve_detailed_balance(
	bool preserve_detailed_balance
)
{
	preserve_detailed_balance_ = preserve_detailed_balance;
}

bool
SidechainMover::change_chi_without_replacing_residue() const
{
	return change_chi_without_replacing_residue_;
}

void
SidechainMover::set_change_chi_without_replacing_residue(
	bool const change_chi_without_replacing_residue
)
{
	change_chi_without_replacing_residue_ = change_chi_without_replacing_residue;
}

utility::vector1<core::id::TorsionID_Range>
SidechainMover::torsion_id_ranges(
	core::pose::Pose & //pose
)
{
	return utility::vector1<core::id::TorsionID_Range>();
}

utility::vector1<core::id::DOF_ID_Range>
SidechainMover::dof_id_ranges(
	core::pose::Pose & pose
)
{
	Real static const pi(numeric::NumericTraits<Real>::pi());

	init_task(pose);

	utility::vector1<core::id::DOF_ID_Range> range_vector;

	for ( core::Size i = 1; i <= packed_residues_.size(); ++i ) {

		core::Size const resnum(packed_residues_[i]);

		bool all_names_match = true;

		// we only monitor DOFs of residues that won't change primary type
		// check to see if pose residue type has the same name as all allowed residue types
		core::pack::task::ResidueLevelTask const & residueleveltask(task_->residue_task(resnum));
		for ( auto iter(residueleveltask.allowed_residue_types_begin());
				iter != residueleveltask.allowed_residue_types_end(); ++iter ) {

			if ( (*iter)->name3() != pose.residue_type(resnum).name3() ) {
				all_names_match = false;
				break;
			}
		}

		if ( !all_names_match ) break;

		for ( core::Size j = 1; j <= pose.residue_type(resnum).nchi(); ++j ) {

			core::id::TorsionID chi_torsion(resnum, core::id::CHI, j);
			core::id::DOF_ID chi_dof(pose.conformation().dof_id_from_torsion_id(chi_torsion));
			if ( chi_dof.valid() ) range_vector.push_back(core::id::DOF_ID_Range(chi_dof, -pi, pi));
		}
	}

	return range_vector;
}

utility::vector1<core::Size> const &
SidechainMover::packed_residues() const
{
	return packed_residues_;
}

utility::vector1<bool> const &
SidechainMover::residue_packed() const
{
	return residue_packed_;
}

core::Size
SidechainMover::next_resnum() const
{
	return next_resnum_;
}

void
SidechainMover::next_resnum(
	core::Size resnum
)
{
	runtime_assert(resnum == 0 || residue_packed_[resnum]);
	next_resnum_ = resnum;
}

core::Size
SidechainMover::last_nchi() const
{
	return last_nchi_;
}

bool
SidechainMover::last_mutation() const
{
	return last_mutation_;
}

bool
SidechainMover::last_uniform() const
{
	return last_uniform_;
}

bool
SidechainMover::last_withinrot() const
{
	return last_withinrot_;
}

core::Real
SidechainMover::last_proposal_density_ratio()
{
	return last_proposal_density_ratio_;
}

/// @details
/// All move types are prefixed with "sc". Sections are divided by underscores.
/// The next section indicates whether a mutation was made ("mut") or not ("chi").
/// The last section indicates wehter chi sampling was uniform ("unif"), used
/// Dunbrack rotamer statistics ("rot"), or whether no chi angles existed in the
/// placed residue ("none").
void
SidechainMover::update_type()
{
	std::stringstream mt;

	mt << "sc_" << (last_mutation_ ? "mut" : "chi") << "_" << (last_nchi_ ? (last_uniform_ ? "unif" : (last_withinrot_ ? "withinrot" : (last_pertrot_ ? "pertrot" : "rot"))) : "none");

	std::string const new_type(mt.str());
	type(new_type);
}

std::string SidechainMover::get_name() const {
	return mover_name();
}

std::string SidechainMover::mover_name() {
	return "Sidechain";
}

void SidechainMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default("ntrials", xsct_non_negative_integer, "Number of trials.", "10000" )
		+ XMLSchemaAttribute::attribute_w_default("preserve_detailed_balance", xsct_rosetta_bool, "Should the simulation preserve detailed balance?", "1" )
		+ XMLSchemaAttribute::attribute_w_default("prob_uniform", xsct_real, "The probability of uniform chi sampling.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("prob_withinrot", xsct_real, "The probability of sampling within-rotamer.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("prob_random_pert_current", xsct_real, "The probability of making a random perturbation to the current chi value.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("change_chi_without_replacing_residue", xsct_real, "The probability of changing chi but not replacing the residue.", "0.0" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Moves the side chain for a set of residues identified by a task operation in a manner that can be totally independent of rotamer assignments", attlist );
}

std::string SidechainMoverCreator::keyname() const {
	return SidechainMover::mover_name();
}

protocols::moves::MoverOP
SidechainMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SidechainMover );
}

void SidechainMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SidechainMover::provide_xml_schema( xsd );
}



/// APL NOTE:
/// The following block of comments was originally in make_move
/// but have been moved here

//Size chinum(1);
// use Dunbrack statistics to bias chi angle sampling
//if (numeric::random::rg().uniform() > prob_uniform_ && residue_dunbrack_library) {
//
// // get indices of the phi and psi bins, plus extra information we don't need
// Size phibin, psibin, phibin_next, psibin_next;
// Real phi_err, psi_err;
// residue_dunbrack_library->get_phi_interpolate_bin(pose.phi(resnum), pose.psi(resnum), phibin, psibin,
//  phibin_next, psibin_next, phi_err, psi_err);

// // choose a random rotamer according to its probability
// Real accprob = numeric::random::rg().uniform();
// for (Size rotnum = 1; rotnum <= residue_dunbrack_library->nrots_per_bin(); rotnum++) {
//
//  core::pack::dunbrack::DunbrackRotamer< core::pack::dunbrack::FOUR > rotamer(
//   residue_dunbrack_library->retrieve_rotamer(phibin, psibin, rotnum)
//  );

//  // keep looking if we haven't decremented the accumulated probability below zero
//  accprob -= rotamer.rotamer_probability();
//  if (accprob > 0) continue;

//  // for chi angles which have statistics, sample according to rotamer SD and mean
//  TR.Debug << "Using rotamer " << rotnum << ": ";
//  for ( ; chinum <= residue_dunbrack_library->nchi_aa(); ++chinum) {
//   TR.Debug << " " << rotamer.chi_mean(chinum) << "(" << rotamer.chi_sd(chinum) << ")";
//   last_chi_angles_[chinum] = radians(numeric::random::rg().gaussian()*rotamer.chi_sd(chinum) + rotamer.chi_mean(chinum));
//   pose.set_chi(chinum, resnum, degrees(last_chi_angles_[chinum]));
//  }
//  TR.Debug << std::endl;

//  last_uniform_ = false;
//  break;
// }
//}

// any/all remaning chi angles are chosen uniformly
//for ( ; chinum <= residue_type->nchi(); ++chinum) {
// last_chi_angles_[chinum] = numeric::random::rg().uniform()*numeric::constants::d::pi_2 - numeric::constants::d::pi;
// pose.set_chi(chinum, resnum, degrees(last_chi_angles_[chinum]));
//}


//  TR << "new chi angles:";
//  utility::vector1< Real > const & nchis(pose.residue(resnum).chi());
//  for (core::Size i = 1; i <= nchis.size(); ++i)
//   TR << " " << nchis[i];
//  TR << std::endl;


// TR.Debug << "Set residue chi angles to:";
// for (Size i = 1; i <= residue_type->nchi(); ++i) TR.Debug << " " << degrees(last_chi_angles_[i]);
// TR.Debug << std::endl;


} // sidechain_moves
} // simple_moves
} // protocols
