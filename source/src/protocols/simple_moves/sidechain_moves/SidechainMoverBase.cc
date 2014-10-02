// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/sidechain_moves/SidechainMoverBase.cc
/// @brief implementation of SidechainMoverBase class and functions
/// @author Oliver Lange ( oliver.lange@tum.de )


#include <protocols/simple_moves/sidechain_moves/SidechainMoverBase.hh>

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
#include <core/scoring/ScoringManager.hh>
#include <core/types.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <basic/prof.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

// Utility
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>

// for debug ( temporary )
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>


// C++ Headers
#include <sstream>
#include <fstream>
#include <utility/excn/Exceptions.hh>
#include <utility/fixedsizearray1.hh>
using namespace core;
using namespace core::pose;

static thread_local basic::Tracer tr( "protocols.simple_moves.sidechain_moves.SidechainMoverBase" );

namespace protocols {
namespace simple_moves {
namespace sidechain_moves {

using namespace chemical;
using namespace conformation;

SidechainMoverBase::SidechainMoverBase():
	rotamer_library_( * pack::dunbrack::RotamerLibrary::get_instance() )
{
	set_defaults();
	init_from_options();
}

SidechainMoverBase::SidechainMoverBase(
	pack::dunbrack::RotamerLibrary const & rotamer_library
):
	rotamer_library_(rotamer_library)
{
	set_defaults();
	init_from_options();
}

SidechainMoverBase::SidechainMoverBase(
	SidechainMoverBase const & mover
) :
	//utility::pointer::ReferenceCount(),
	protocols::canonical_sampling::ThermodynamicMover(mover),
	rotamer_library_(mover.rotamer_library_),
	packed_residues_(mover.packed_residues_),
	residue_packed_(mover.residue_packed_),
	preserve_detailed_balance_(mover.preserve_detailed_balance_),
	change_chi_without_replacing_residue_(mover.change_chi_without_replacing_residue_),
	last_chi_angles_(mover.last_chi_angles_),
	last_nchi_(mover.last_nchi_),
	last_proposal_density_ratio_(mover.last_proposal_density_ratio_)
{
	if (mover.task_factory_) task_factory_ = core::pack::task::TaskFactoryCOP( new pack::task::TaskFactory(*mover.task_factory_) );
	if (mover.task_) task_ = mover.task_->clone();
}

SidechainMoverBase::~SidechainMoverBase() {}

void
SidechainMoverBase::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	pose::Pose const & /*pose*/
)
{
	using core::pack::task::operation::TaskOperationCOP;
	pack::task::TaskFactoryOP new_task_factory( new pack::task::TaskFactory );

	if ( tag->hasOption("task_operations") ) {
		std::string const t_o_val( tag->getOption<std::string>("task_operations") );
		typedef utility::vector1< std::string > StringVec;
		StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );
		for ( StringVec::const_iterator t_o_key( t_o_keys.begin() ), end( t_o_keys.end() );
					t_o_key != end; ++t_o_key ) {
			if ( data.has( "task_operations", *t_o_key ) ) {
				new_task_factory->push_back( data.get_ptr< pack::task::operation::TaskOperation >( "task_operations", *t_o_key ) );
			} else {
				throw utility::excn::EXCN_RosettaScriptsOption("TaskOperation " + *t_o_key + " not found in basic::datacache::DataMap.");
			}
		}
	} else {
		new_task_factory->push_back( TaskOperationCOP( new pack::task::operation::RestrictToRepacking ) );
	}

	task_factory_ = new_task_factory;

	set_preserve_detailed_balance( tag->getOption<bool>( "preserve_detailed_balance", preserve_detailed_balance() ) );
	set_change_chi_without_replacing_residue( tag->getOption<Real>( "change_chi_without_replacing_residue", change_chi_without_replacing_residue() ) );
}

void
SidechainMoverBase::set_defaults() {
	change_chi_without_replacing_residue_=false;
	last_proposal_density_ratio_=1;
}

void
SidechainMoverBase::init_from_options() {
	using namespace basic::options;
	using core::pack::task::operation::TaskOperationCOP;
	core::pack::task::TaskFactoryOP new_task_factory( new pack::task::TaskFactory );
	new_task_factory->push_back( TaskOperationCOP( new pack::task::operation::InitializeFromCommandline ) );
	// design is not supported yet !!!
	new_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );
	if ( option[ OptionKeys::packing::resfile ].user() ) {
		new_task_factory->push_back( TaskOperationCOP( new pack::task::operation::ReadResfile ) );
	}
	set_task_factory( new_task_factory );
}


/// @detailed
/// Check to make sure that a packer task exists and matches the numer of residues in
/// the given pose. If that isn't the case, create a new one with the task factory.
/// Exits with an error if no task factory exists.
void
SidechainMoverBase::init_task( pose::Pose const &pose ) {
	// does a valid task exist ?
	if ( !task_ || task_->total_residue() != pose.total_residue() ) {
		if ( task_factory_ ) { //no - create a new one
			set_task(task_factory_->create_task_and_apply_taskoperations( pose ));
		} else {
			utility_exit_with_message("Cannot create task because no task factory is set");
		}
	}
}

void
SidechainMoverBase::initialize_simulation(
	pose::Pose & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const&, /*metropolis_hastings_mover*/
	core::Size /*cycle*/
) {
	idealize_sidechains( pose );
	init_task( pose );
}

core::Size SidechainMoverBase::suggest_residue_number( pose::Pose const& pose ) const {
	core::Size resid;
	do { //pick residue, respecting underlying packer task
		resid = packed_residues_[ numeric::random::rg().random_range(1, packed_residues_.size()) ];
	} while ( pose.residue( resid ).name1() == 'P'); //SidechainMover cannot sample proline rotamers? (is this true?)
	return resid;
}

/// @detailed
void
SidechainMoverBase::apply( Pose& pose ) {
	using numeric::conversions::degrees;
	using numeric::conversions::radians;

	init_task(pose);

	//select_resnum
	Size const resnum( suggest_residue_number( pose ) );

	ResidueOP newresidue( new Residue( pose.residue( resnum ) ) );
	ResidueOP final = make_move( newresidue );

	if ( !change_chi_without_replacing_residue() || have_mutated_residue() ) {
		pose.replace_residue( resnum, *final, true );
	} else {
		Size const nchi( final->nchi() );
		for ( Size i(1); i<=nchi; ++i ) {
			pose.set_chi(i, resnum, final->chi(i));
		} //for
	} //else
} //apply



conformation::ResidueOP
SidechainMoverBase::make_move( conformation::ResidueOP old_res )
{

	using numeric::conversions::degrees;
	using numeric::conversions::radians;
	using namespace ObjexxFCL;
	chemical::ResidueType  const& old_res_type( old_res->type() );
	chemical::ResidueTypeCOP new_res_type( old_res->type().get_self_ptr() ); //for now until we fix the design stuff back in...
	utility::vector1<Real> const old_chi( old_res->chi() );
	Size resnum = old_res->seqpos();
	Size nchi( old_chi.size() );

	utility::vector1< Real > new_chi;
	new_chi.reserve( 4 );
	new_chi.resize( nchi );

	make_chi_move( *old_res, old_chi, new_chi );
	tr << "chi1 (old/new): " << format::F(8,4,old_chi[ 1 ]) << " " << format::F( 8,4,new_chi[ 1 ]) << std::endl;
	Real proposal_density_reverse(1);
	Real proposal_density_forward(1);

	if (preserve_detailed_balance_) {
		proposal_density_reverse = compute_proposal_density(*old_res, resnum, *new_res_type, new_chi);
	}
	tr << "chi1 (old/new): " << format::F(8,4,old_chi[ 1 ]) << " " << format::F( 8,4,new_chi[ 1 ]) << std::endl;
	/// set the chi
	conformation::ResidueOP new_residue = old_res;
	conformation::ResidueOP previous_residue = old_res;
	if (  false ) { /* replace by	if (residue_type() != & previous_residue->type() ) */
		/* do design stuff */
		/* new_residue = < a new residue > */
		/*		new_residue =
			conformation::ResidueFactory::create_residue(*residue_type,
			*old_res,
			pose_->conformation(),
			residue_task.preserve_c_beta());
			for( Size ii = 1; ii <= residue_type->nchi(); ii++ ){
			new_residue->set_chi( ii, numeric::principal_angle_degrees(last_chi_angles_[ ii ] ));
		*/
	} else {
		for (Size i = 1; i <= old_res->nchi(); ++i) {
			new_residue->set_chi( i, new_chi[ i ] );
		}
	}
	tr << "chi1 (old/new): " << format::F(8,4,old_chi[ 1 ]) << " " << format::F( 8,4,new_chi[ 1 ]) << std::endl;
	if ( preserve_detailed_balance_ ) {
		proposal_density_forward = compute_proposal_density(*new_residue, resnum, old_res_type, old_chi);
	}

	tr << "chi1 (old/new): " << format::F(8,4,old_chi[ 1 ]) << " " << format::F( 8,4,new_chi[ 1 ]) << " reverse/forward" << format::E( 10,5, proposal_density_reverse) << " " << format::E(10,5,proposal_density_forward ) << std::endl;
	last_proposal_density_ratio_ = proposal_density_reverse / proposal_density_forward;
	return new_residue;

}



std::string
SidechainMoverBase::get_name() const {
	return "SidechainMoverBase";
}

/// @detailed
/// all sidechains that might be changed are replaced with ideal coordinates that have
/// the original chi angles
void
SidechainMoverBase::idealize_sidechains(
	pose::Pose & pose
)
{
	utility::vector1<Real> chi;
	init_task(pose);
	for (Size i = 1; i <= packed_residues_.size(); ++i) {
		Size const resnum(packed_residues_[i]);

		// get the residue type and residue level task
		ResidueType const & residue_type(pose.residue_type(resnum));
		pack::task::ResidueLevelTask const & residue_task(task_->residue_task(resnum));

		// disable proline idealization for now
		if (residue_type.aa() == aa_pro) continue;

		// save original chi angles
		if (residue_type.nchi() > chi.size()) chi.resize(residue_type.nchi());
		for (Size chinum = 1; chinum <= residue_type.nchi(); ++chinum) {
			chi[chinum] = pose.chi(chinum, resnum);
		}

		// create a new residue and replace the old one with it
		ResidueOP new_residue(
			ResidueFactory::create_residue(residue_type, pose.residue(resnum), pose.conformation(),
			                                             residue_task.preserve_c_beta())
		);
		pose.replace_residue(resnum, *new_residue, false);

		// put original chi angles back
		for (Size chinum = 1; chinum <= residue_type.nchi(); ++chinum) {
			pose.set_chi(chinum, resnum, chi[chinum]);
		}
	}
}

pack::dunbrack::RotamerLibrary const &
SidechainMoverBase::rotamer_library() const {
	return rotamer_library_;
}

pack::task::TaskFactoryCOP
SidechainMoverBase::task_factory() const {
	return task_factory_;
}

void
SidechainMoverBase::set_task_factory(	pack::task::TaskFactoryCOP task_factory ) {
	task_factory_ = task_factory;
}

pack::task::PackerTaskCOP
SidechainMoverBase::task() const {
	return task_;
}

void
SidechainMoverBase::set_task(	pack::task::PackerTaskCOP task ) {
	using pack::task::ResidueLevelTask;

	task_ = task;

	// update the list of residues being packed
	packed_residues_.clear();
	residue_packed_.resize(task_->total_residue());
	for (Size i = 1; i <= task_->total_residue(); ++i) {

		// loop over all possible residue types and check if any have chi angles
		residue_packed_[i] = false;
		pack::task::ResidueLevelTask const & residue_task(task->residue_task(i));
		for (ResidueLevelTask::ResidueTypeCOPListConstIter iter(residue_task.allowed_residue_types_begin());
		     iter != residue_task.allowed_residue_types_end(); ++iter) {

			// temporarily exclude proline residues from side chain sampling
			if ((*iter)->nchi() > 0 && (*iter)->aa() != aa_pro) {
				packed_residues_.push_back(i);
				residue_packed_[i] = true;
				break;
			}
		}
	}
}

bool
SidechainMoverBase::preserve_detailed_balance() const {
	return preserve_detailed_balance_;
}

void
SidechainMoverBase::set_preserve_detailed_balance(
	bool preserve_detailed_balance
) {
	preserve_detailed_balance_ = preserve_detailed_balance;
}

bool
SidechainMoverBase::change_chi_without_replacing_residue() const {
	return change_chi_without_replacing_residue_;
}

void
SidechainMoverBase::set_change_chi_without_replacing_residue(
	bool const change_chi_without_replacing_residue ) {
	change_chi_without_replacing_residue_ = change_chi_without_replacing_residue;
}

utility::vector1<id::TorsionID_Range>
SidechainMoverBase::torsion_id_ranges( pose::Pose & /*pose*/ ) {
	return utility::vector1<id::TorsionID_Range>();
}

utility::vector1<id::DOF_ID_Range>
SidechainMoverBase::dof_id_ranges( pose::Pose & pose ) {
	Real static const pi(numeric::NumericTraits<Real>::pi());

	init_task(pose);

	utility::vector1<id::DOF_ID_Range> range_vector;

	for (Size i = 1; i <= packed_residues_.size(); ++i) {

		Size const resnum(packed_residues_[i]);

		bool all_names_match = true;

		// we only monitor DOFs of residues that won't change primary type
		// check to see if pose residue type has the same name as all allowed residue types
		pack::task::ResidueLevelTask const & residueleveltask(task_->residue_task(resnum));
		for (pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter iter(residueleveltask.allowed_residue_types_begin());
				 iter != residueleveltask.allowed_residue_types_end(); ++iter) {

			if ((*iter)->name3() != pose.residue_type(resnum).name3()) {
				all_names_match = false;
				break;
			}
		}

		if (!all_names_match) break;

		for (Size j = 1; j <= pose.residue_type(resnum).nchi(); ++j) {

			id::TorsionID chi_torsion(resnum, id::CHI, j);
			id::DOF_ID chi_dof(pose.conformation().dof_id_from_torsion_id(chi_torsion));
			if (chi_dof.valid()) range_vector.push_back(id::DOF_ID_Range(chi_dof, -pi, pi));
		}
	}

	return range_vector;
}

utility::vector1<Size> const &
SidechainMoverBase::packed_residues() const {
	return packed_residues_;
}

utility::vector1<bool> const &
SidechainMoverBase::residue_packed() const {
	return residue_packed_;
}

Real SidechainMoverBase::last_proposal_density_ratio() {
	return last_proposal_density_ratio_;
}


} // sidechain_moves
} // simple_moves
} // protocols
