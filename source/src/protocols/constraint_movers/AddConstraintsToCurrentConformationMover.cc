// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#include <protocols/constraint_movers/AddConstraintsToCurrentConformationMover.hh>
#include <protocols/constraint_movers/AddConstraintsToCurrentConformationMoverCreator.hh>

// protocol headers
#include <protocols/constraint_generator/util.hh>
#include <protocols/residue_selectors/TaskSelector.hh>
#include <protocols/rosetta_scripts/util.hh>

// core headers
#include <core/pose/Pose.hh>

#include <core/chemical/ResidueType.hh>

#include <core/conformation/Residue.hh>



#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/SOGFunc.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>

#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// task operation
#include <core/pack/task/TaskFactory.fwd.hh>
#include <numeric/xyzVector.hh>

// utility
#include <utility/tag/Tag.hh>
#include <utility/fixedsizearray1.hh>
#include <basic/Tracer.hh>

// option
#include <basic/options/option.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.constraint_movers.AddConstraintsToCurrentConformationMover" );

namespace protocols {
namespace constraint_movers {

using namespace core;
using namespace basic::options;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;
using namespace constraints;

//Debugging utility, defined at the bottom of the page
#ifndef NDEBUG
void assert_that_all_atoms_are_moved_by_chi_angle( conformation::Residue const & res, Size chino, Size atomno_begin, Size atomno_end );
#endif


AddConstraintsToCurrentConformationMover::AddConstraintsToCurrentConformationMover():
	protocols::moves::Mover( AddConstraintsToCurrentConformationMover::mover_name() ),
	use_distance_cst_( false ),
	use_harmonic_func_( false ),
	inter_chain_( true ),
	cst_weight_( 1.0 ),
	max_distance_( 12.0 ),
	coord_dev_( 1.0 ),
	bound_width_( 0. ),
	min_seq_sep_( 8 ),
	selector_( new core::select::residue_selector::TrueResidueSelector ),
	use_bounded_func_( false )
{}

AddConstraintsToCurrentConformationMover::~AddConstraintsToCurrentConformationMover() = default;

void
AddConstraintsToCurrentConformationMover::apply( core::pose::Pose & pose )
{
	// Check that boolean logic options make sense
	runtime_assert( ! ( use_harmonic_func_ && use_bounded_func_ ) );
	if ( use_harmonic_func_ ) {
		TR << "Using harmonic function" << std::endl;
	} else if ( use_bounded_func_ ) {
		TR << "Using bounded function" << std::endl;
	} else {
		TR << "Using SOG function" << std::endl;
	}

	pose.add_constraints( generate_constraints( pose ) );
}

core::scoring::constraints::ConstraintCOPs
AddConstraintsToCurrentConformationMover::generate_constraints( core::pose::Pose const & pose )
{
	using core::id::AtomID;
	// generate residue subset
	if ( !selector_ ) {
		debug_assert( selector_ );
		throw CREATE_EXCEPTION(utility::excn::BadInput,  "Selector not set in AddConstraintsToCurrentConformationMover::generate_constraints()\n" );
	}
	core::select::residue_selector::ResidueSubset const subset = selector_->apply( pose );

	if ( !use_distance_cst_ ) {
		// this is not quite right without adding a virtual residue
		/*
		if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
		core::pose::addVirtualResAsRoot(pose);
		}
		*/
		return generate_coordinate_constraints( pose, subset );
	} else {
		return generate_atom_pair_constraints( pose, subset );
	}
}

core::Size
AddConstraintsToCurrentConformationMover::find_best_anchor( core::pose::Pose const & pose ) const
{
	core::Size const nres = pose.size();

	// find anchor residue
	numeric::xyzVector< core::Real > sum_xyz(0.0);
	numeric::xyzVector< core::Real > anchor_xyz(0.0);
	core::Real natom = 0.0;
	for ( core::Size ires=1; ires<=nres; ++ires ) {
		if ( pose.residue_type(ires).has("CA") ) {
			core::Size const iatom = pose.residue_type(ires).atom_index("CA");
			sum_xyz += pose.residue(ires).xyz(iatom);
			natom += 1.;
		}
		if ( natom > 1e-3 ) {
			anchor_xyz = sum_xyz / natom;
		}
	}
	core::Real min_dist2 = 1e9;
	core::Size best_anchor = 0;
	for ( core::Size ires=1; ires<=nres; ++ires ) {
		if ( pose.residue_type(ires).has("CA") ) {
			core::Size const iatom = pose.residue_type(ires).atom_index("CA");
			core::Real const dist2 = pose.residue(ires).xyz(iatom).distance_squared(anchor_xyz);
			if ( dist2 < min_dist2 ) {
				min_dist2 = dist2;
				best_anchor = ires;
			}
		}
	}
	return best_anchor;
}

core::scoring::constraints::ConstraintCOPs
AddConstraintsToCurrentConformationMover::generate_coordinate_constraints(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & subset ) const
{
	core::scoring::constraints::ConstraintCOPs csts;

	// TR << pose.fold_tree() << std::endl;
	core::Size const best_anchor_resid = find_best_anchor( pose );

	if ( best_anchor_resid == 0 ) {
		TR << "Best anchor resid == 0, not generating any constraints" << std::endl;
		return core::scoring::constraints::ConstraintCOPs();
	}
	core::Size const best_anchor_atom = pose.residue_type( best_anchor_resid ).atom_index("CA");
	core::id::AtomID const best_anchor_id( best_anchor_atom, best_anchor_resid );

	for ( core::Size ires=1; ires<=pose.size(); ++ires ) {
		if ( !subset[ ires ] ) continue;

		// find atom start, stop indices
		core::Size iatom_start=1, iatom_stop=pose.residue(ires).nheavyatoms();
		if ( pose.residue_type(ires).is_DNA() ) {
			iatom_stop = 0;
		} else if ( pose.residue_type(ires).is_protein() ) {
			if ( CA_only() && pose.residue_type(ires).has("CA") ) {
				iatom_start = iatom_stop = pose.residue_type(ires).atom_index("CA");
			} else if ( bb_only() ) {
				iatom_stop = pose.residue_type(ires).last_backbone_atom();
			} else if ( sc_tip_only() ) {
				utility::vector1< core::Size > const & final_chi_atoms =
					pose.residue( ires ).chi_atoms().back();;
#ifndef NDEBUG
				//core::Size const nchi = pose.residue( ires ).nchi();
				//core::Size const final_chi_atom = final_chi_atoms.back()
				//assert_that_all_atoms_are_moved_by_chi_angle( pose.residue( ires ), nchi, final_chi_atom, iatom_stop );
				//This assert is bad because it doesn't handle residues like TYR where there are atoms that shouldn't move
#endif
				iatom_start = final_chi_atoms[2]; //[1] moves relative to the tip, [2]-[4] do not
			}
		} else {
			continue;
		}
		TR.Debug << "Registering constraints for resid " << ires << " from atom " << iatom_start << " through " << iatom_stop << std::endl;

		// add constraints
		core::scoring::func::FuncOP cc_func; //NULL
		if ( !use_bounded_func_ || use_harmonic_func_ ) {
			cc_func = utility::pointer::make_shared< core::scoring::func::HarmonicFunc >( 0.0, coord_dev_ );
		} else {
			cc_func = utility::pointer::make_shared< BoundFunc >( 0, bound_width_, coord_dev_, "xyz" );
		}
		runtime_assert(cc_func);


		for ( core::Size iatom=iatom_start; iatom<=iatom_stop; ++iatom ) {
			csts.push_back( utility::pointer::make_shared< CoordinateConstraint >(
				core::id::AtomID(iatom,ires), best_anchor_id, pose.residue(ires).xyz(iatom), cc_func ) );
			TR.Debug << "coordinate constraint generated for residue " << ires << ", atom " << iatom << ", using func" << std::endl << *cc_func << std::endl;
		}
	} //loop through residues
	return csts;
}

core::scoring::constraints::ConstraintCOPs
AddConstraintsToCurrentConformationMover::generate_atom_pair_constraints(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & subset ) const
{
	if ( sc_tip_only() ) {
		utility_exit_with_message( "generate_atom_pair_constraints is not yet compatible with sc_tip_only" );
	}

	core::scoring::constraints::ConstraintCOPs csts;

	core::Size const nres =
		protocols::constraint_generator::compute_nres_in_asymmetric_unit( pose );

	for ( core::Size ires=1; ires<=nres; ++ires ) {
		if ( !subset[ ires ] ) continue;
		if ( pose.residue(ires).aa() == core::chemical::aa_vrt ) continue;

		utility::fixedsizearray1<core::Size,2> iatoms(0); // both are 0
		if ( pose.residue_type(ires).has("CA")  ) {
			iatoms[1] = pose.residue_type(ires).atom_index("CA");
		}
		if ( !CA_only() ) {
			iatoms[2] = pose.residue_type(ires).nbr_atom();
		}
		///I think this fails across chains
		///probably residues on different chains should always have sufficient
		///sequence separation?? --- SML dec 21 16
		for ( core::Size jres=ires+min_seq_sep_; jres<=pose.size(); ++jres ) {
			if ( !subset[ jres ] ) continue;
			if ( pose.residue(jres).aa() == core::chemical::aa_vrt ) continue;
			if ( !inter_chain_ && pose.chain(ires)!=pose.chain(jres) ) continue;
			utility::fixedsizearray1<core::Size,2> jatoms(0);
			if ( pose.residue_type(jres).has("CA") ) {
				jatoms[1] = pose.residue_type(jres).atom_index("CA");
			}
			if ( !CA_only() ) {
				jatoms[2] = pose.residue_type(jres).nbr_atom();
			}

			for ( utility::fixedsizearray1<core::Size,2>::const_iterator iiatom=iatoms.begin(); iiatom!=iatoms.end(); ++iiatom ) {
				for ( utility::fixedsizearray1<core::Size,2>::const_iterator jjatom=jatoms.begin(); jjatom!=jatoms.end(); ++jjatom ) {
					core::Size const &iatom(*iiatom), &jatom(*jjatom);
					if ( iatom==0 || jatom==0 ) continue;

					core::Real const dist = pose.residue(ires).xyz(iatom).distance( pose.residue(jres).xyz(jatom) );
					if ( dist > max_distance_ ) continue;

					core::scoring::func::FuncOP apc_func; //NULL!
					if ( use_harmonic_func_ ) {
						apc_func = utility::pointer::make_shared< core::scoring::func::HarmonicFunc >( dist, coord_dev_ );
						TR.Debug << "Using harmonic func with AtomPairConstraint: " << *apc_func << std::endl;
					} else if ( use_bounded_func_ ) {
						//establish bound_func limits based on basin around distance
						core::Real lower_limit( dist - ( 0.5  * bound_width_) );
						if ( lower_limit < 0 ) lower_limit = 0;
						core::Real const upper_limit( dist + ( 0.5 * bound_width_ ) );
						apc_func = utility::pointer::make_shared< BoundFunc >( lower_limit, upper_limit, coord_dev_, "xyz" );
						TR.Debug << "Using Bounded func with AtomPairConstraint: " << *apc_func << std::endl;
					} else {
						core::scoring::func::FuncOP weighted_func;
						core::scoring::func::FuncOP sog_func( new core::scoring::func::SOGFunc( dist, coord_dev_ ) );
						apc_func = utility::pointer::make_shared< core::scoring::func::ScalarWeightedFunc >( cst_weight_, sog_func );
						TR.Debug << "Using SOG func with AtomPairConstraint: " << *apc_func << std::endl;
					}
					runtime_assert(apc_func); //make sure that func got populated; should have been forced down one of the two if branches above

					core::scoring::constraints::ConstraintCOP newcst(
						new core::scoring::constraints::AtomPairConstraint( core::id::AtomID( iatom, ires ), core::id::AtomID( jatom, jres ), apc_func ) );
					csts.push_back( newcst );
					if ( use_harmonic_func_ ) {
						TR.Debug << "harmonic atom_pair_constraint generated for residue " << ires << ", atom " << iatom << " and residue " << jres << ", atom " << jatom  << std::endl;
					} else if ( use_bounded_func_ ) {
						TR.Debug << "bounded atom_pair_constraint generated for residue " << ires << ", atom " << iatom << " and residue " << jres << ", atom " << jatom << " with weight " << cst_weight_ << std::endl;
					} else {
						TR.Debug << "SOG atom_pair_constraint generated for residue " << ires << ", atom " << iatom << " and residue " << jres << ", atom " << jatom << " with weight " << cst_weight_ << std::endl;
					}
				}
			}
		} // jres loop
	} // ires loop
	return csts;
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
AddConstraintsToCurrentConformationMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap
) {
	use_distance_cst_ = tag->getOption< bool >( "use_distance_cst", use_distance_cst_ );
	max_distance_ = tag->getOption< core::Real >( "max_distance", max_distance_ );
	coord_dev_ = tag->getOption< core::Real >( "coord_dev", coord_dev_ );
	bound_width_ = tag->getOption< core::Real >( "bound_width", bound_width_ );
	min_seq_sep_ = tag->getOption< core::Size>( "min_seq_sep", min_seq_sep_ );
	cst_weight_ = tag->getOption< core::Real >( "cst_weight", cst_weight_ );
	inter_chain_ = tag->getOption< bool >( "inter_chain", inter_chain_ );

	AtomSelector sele = AtomSelector::ALL;

	if ( tag->hasOption( "bb_only" ) ) {
		bool const bb_only = tag->getOption< bool >( "bb_only" );
		if ( bb_only ) {
			runtime_assert_msg( sele == AtomSelector::ALL, "Please use only one of 'bb_only', 'ca_only', and 'sc_tip_only'" );
			sele = AtomSelector::BB_ONLY;
		}
	}

	if ( tag->hasOption( "CA_only" ) ) {
		bool const CA_only = tag->getOption< bool >( "CA_only" );
		if ( CA_only ) {
			runtime_assert_msg( sele == AtomSelector::ALL, "Please use only one of 'bb_only', 'ca_only', and 'sc_tip_only'" );
			sele = AtomSelector::CA_ONLY;
		}
	}

	if ( tag->hasOption( "sc_tip_only" ) ) {
		bool const sc_tip_only = tag->getOption< bool >( "sc_tip_only" );
		if ( sc_tip_only ) {
			runtime_assert_msg( sele == AtomSelector::ALL, "Please use only one of 'bb_only', 'ca_only', and 'sc_tip_only'" );
			sele = AtomSelector::SC_TIP_ONLY;
		}
	}

	set_atom_selector( sele );

	if ( bound_width_ < 1e-3 && cst_weight_ < 1e-3 ) {
		use_harmonic_func_ = true;
	} else if ( bound_width_ < 1e-3 ) {
		use_bounded_func_ = false;
	} else {
		use_bounded_func_ = true;
	}

	if ( tag->hasOption( "task_operations" ) ) {
		TR.Warning << "task_operations only active for proteins" << std::endl;
		parse_task_operations( tag, datamap );
	}

	core::select::residue_selector::ResidueSelectorCOP selector = protocols::rosetta_scripts::parse_residue_selector( tag, datamap );
	if ( selector ) residue_selector( selector );

	if ( tag->hasOption( "task_operations" ) && tag->hasOption( "residue_selector" ) ) {
		std::stringstream msg;
		msg << "AddConstraintsToCurrentConformationMover::parse_my_tag(): 'task_operations' and 'residue_selector' cannot both be specified." << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
	}
}

void
AddConstraintsToCurrentConformationMover::parse_task_operations(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap
) {
	TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == nullptr ) return;
	task_factory( new_task_factory );
}

void
AddConstraintsToCurrentConformationMover::residue_selector( core::select::residue_selector::ResidueSelectorCOP selector )
{
	selector_ = selector;
}

void AddConstraintsToCurrentConformationMover::task_factory( TaskFactoryOP tf )
{
	runtime_assert( tf != nullptr );
	protocols::residue_selectors::TaskSelectorOP task_selector( new protocols::residue_selectors::TaskSelector( tf, false, false, true ) );
	residue_selector( task_selector );
}

moves::MoverOP AddConstraintsToCurrentConformationMover::clone() const {
	return utility::pointer::make_shared< AddConstraintsToCurrentConformationMover >( *this );
}
moves::MoverOP AddConstraintsToCurrentConformationMover::fresh_instance() const {
	return utility::pointer::make_shared< AddConstraintsToCurrentConformationMover >();
}





std::string AddConstraintsToCurrentConformationMover::get_name() const {
	return mover_name();
}

std::string AddConstraintsToCurrentConformationMover::mover_name() {
	return "AddConstraintsToCurrentConformationMover";
}

void AddConstraintsToCurrentConformationMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "use_distance_cst", xsct_rosetta_bool, "use distance constraints instead of CoordinateConstraints. Probable default false." )
		+ XMLSchemaAttribute( "max_distance", xsct_real, "do not generate distance constraints beyond this distance.  Only active with use_distance_cst." )
		+ XMLSchemaAttribute( "coord_dev", xsct_real, "width (sd) for HarmonicFunc or BoundFunc." )
		+ XMLSchemaAttribute( "bound_width", xsct_real, "BoundFunc zero basin width BoundFunc; also activates use of BoundFunc (if non-zero)" )
		+ XMLSchemaAttribute( "min_seq_sep", xsct_non_negative_integer, "Do not generate distance constraints between residues within this sequence separation.  Only active with use_distance_cst." )
		+ XMLSchemaAttribute( "cst_weight", xsct_real, "use ScalarWeightedFunc to reweight constraints by this; also activates use of HarmonicFunc (if this and bound_width are both zero)" )
		+ XMLSchemaAttribute( "CA_only", xsct_rosetta_bool, "constrain only CA atoms. Sets bb_only and sc_tip_only to false." )
		+ XMLSchemaAttribute( "bb_only", xsct_rosetta_bool, "constrain only backbone atoms. Sets CA_only and sc_tip_only to false." )
		+ XMLSchemaAttribute( "sc_tip_only", xsct_rosetta_bool, "constrain only atoms affected by the final chi angle. Sets CA_only and bb_only to false." )
		+ XMLSchemaAttribute( "inter_chain", xsct_rosetta_bool, "Generate distance constraints between residues on different chains if true.  (Does not appear to generate ONLY interchain constraints.)  If false, skips constraints that would go between chains.  Only active with use_distance_cst." );

	rosetta_scripts::attributes_for_parse_task_operations(attlist);
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector", "make constraints for these residues; mutually exclisuve with task_operations" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string AddConstraintsToCurrentConformationMoverCreator::keyname() const {
	return AddConstraintsToCurrentConformationMover::mover_name();
}

protocols::moves::MoverOP
AddConstraintsToCurrentConformationMoverCreator::create_mover() const {
	return utility::pointer::make_shared< AddConstraintsToCurrentConformationMover >();
}

void AddConstraintsToCurrentConformationMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddConstraintsToCurrentConformationMover::provide_xml_schema( xsd );
}

//Debugging utilities
#ifndef NDEBUG
void
assert_that_all_atoms_are_moved_by_chi_angle(
	core::conformation::Residue const & res,
	core::Size const chino,
	core::Size const atomno_begin,
	core::Size const atomno_end
) {
	core::conformation::Residue res0( res );
	core::conformation::Residue res180( res );
	res0.set_chi(   chino, 0   );
	res180.set_chi( chino, 180 );
	for ( core::Size atom = atomno_begin; atom <= atomno_end; ++atom ) {
		auto const & xyz0 = res0.xyz( atom );
		auto const & xyz180 = res180.xyz( atom );
		//This isn't perfect because it doesn't handle atoms that are colinear with the rotation axis.
		debug_assert( xyz0.distance( xyz180 ) != 0 );
	}
}
#endif

} // moves
} // protocols
