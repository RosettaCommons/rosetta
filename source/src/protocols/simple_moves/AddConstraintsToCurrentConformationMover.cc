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

#include <protocols/simple_moves/AddConstraintsToCurrentConformationMover.hh>
#include <protocols/simple_moves/AddConstraintsToCurrentConformationMoverCreator.hh>

// protocol headers
#include <protocols/constraint_generator/util.hh>
#include <protocols/residue_selectors/TaskSelector.hh>
#include <protocols/rosetta_scripts/util.hh>

// core headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/SOGFunc.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>

#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// task operation
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <numeric/xyzVector.hh>

// utility
#include <utility/tag/Tag.hh>
#include <utility/fixedsizearray1.hh>
#include <basic/Tracer.hh>

// option
#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.AddConstraintsToCurrentConformationMover" );

namespace protocols {
namespace simple_moves {

using namespace core;
using namespace basic::options;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;
using namespace constraints;

AddConstraintsToCurrentConformationMover::AddConstraintsToCurrentConformationMover():
	protocols::moves::Mover( AddConstraintsToCurrentConformationMover::mover_name() ),
	use_distance_cst_( false ),
	CA_only_( true ),
	bb_only_( false ),
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
	pose.add_constraints( generate_constraints( pose ) );
}

core::scoring::constraints::ConstraintCOPs
AddConstraintsToCurrentConformationMover::generate_constraints( core::pose::Pose const & pose )
{
	using core::id::AtomID;
	// generate residue subset
	if ( !selector_ ) {
		debug_assert( selector_ );
		throw utility::excn::EXCN_BadInput( "Selector not set in AddConstraintsToCurrentConformationMover::generate_constraints()\n" );
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
	Size const nres = pose.size();

	// find anchor residue
	numeric::xyzVector< core::Real > sum_xyz(0.0);
	numeric::xyzVector< core::Real > anchor_xyz(0.0);
	core::Real natom = 0.0;
	for ( Size ires=1; ires<=nres; ++ires ) {
		if ( pose.residue_type(ires).has("CA") ) {
			Size const iatom = pose.residue_type(ires).atom_index("CA");
			sum_xyz += pose.residue(ires).xyz(iatom);
			natom += 1.;
		}
		if ( natom > 1e-3 ) {
			anchor_xyz = sum_xyz / natom;
		}
	}
	core::Real min_dist2 = 1e9;
	Size best_anchor = 0;
	for ( Size ires=1; ires<=nres; ++ires ) {
		if ( pose.residue_type(ires).has("CA") ) {
			Size const iatom = pose.residue_type(ires).atom_index("CA");
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
	Size const best_anchor_atom = pose.residue_type( best_anchor_resid ).atom_index("CA");
	core::id::AtomID const best_anchor_id( best_anchor_atom, best_anchor_resid );

	for ( Size ires=1; ires<=pose.size(); ++ires ) {
		if ( !subset[ ires ] ) continue;

		// find atom start, stop indices
		Size iatom_start=1, iatom_stop=pose.residue(ires).nheavyatoms();
		if ( pose.residue_type(ires).is_DNA() ) {
			iatom_stop = 0;
		} else if ( pose.residue_type(ires).is_protein() ) {
			if ( CA_only_ && pose.residue_type(ires).has("CA") ) {
				iatom_start = iatom_stop = pose.residue_type(ires).atom_index("CA");
			} else if ( bb_only_ ) {
				iatom_stop = pose.residue_type(ires).last_backbone_atom();
			}
		} else {
			continue;
		}

		// add constraints
		core::scoring::func::FuncOP cc_func; //NULL
		if ( !use_bounded_func_ ) {
			cc_func = core::scoring::func::FuncOP( new core::scoring::func::HarmonicFunc( 0.0, coord_dev_ ) );
		} else {
			cc_func = core::scoring::func::FuncOP( new BoundFunc( 0, bound_width_, coord_dev_, "xyz" ) );
		}
		runtime_assert(cc_func);


		for ( Size iatom=iatom_start; iatom<=iatom_stop; ++iatom ) {
			csts.push_back( scoring::constraints::ConstraintCOP( new CoordinateConstraint(
				core::id::AtomID(iatom,ires), best_anchor_id, pose.residue(ires).xyz(iatom), cc_func ) ) );
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
	core::scoring::constraints::ConstraintCOPs csts;

	Size const nres =
		protocols::constraint_generator::compute_nres_in_asymmetric_unit( pose );

	for ( Size ires=1; ires<=nres; ++ires ) {
		if ( !subset[ ires ] ) continue;
		if ( pose.residue(ires).aa() == core::chemical::aa_vrt ) continue;
		utility::fixedsizearray1<core::Size,2> iatoms(0); // both are 0
		if ( pose.residue_type(ires).has("CA")  ) {
			iatoms[1] = pose.residue_type(ires).atom_index("CA");
		}
		if ( !CA_only_ ) {
			iatoms[2] = pose.residue_type(ires).nbr_atom();
		}
		///I think this fails across chains
		///probably residues on different chains should always have sufficient
		///sequence separation?? --- SML dec 21 16
		for ( Size jres=ires+min_seq_sep_; jres<=pose.size(); ++jres ) {
			if ( !subset[ jres ] ) continue;
			if ( pose.residue(jres).aa() == core::chemical::aa_vrt ) continue;
			if ( !inter_chain_ && pose.chain(ires)!=pose.chain(jres) ) continue;
			utility::fixedsizearray1<core::Size,2> jatoms(0);
			if ( pose.residue_type(jres).has("CA") ) {
				jatoms[1] = pose.residue_type(jres).atom_index("CA");
			}
			if ( !CA_only_ ) {
				jatoms[2] = pose.residue_type(jres).nbr_atom();
			}

			for ( utility::fixedsizearray1<core::Size,2>::const_iterator iiatom=iatoms.begin(); iiatom!=iatoms.end(); ++iiatom ) {
				for ( utility::fixedsizearray1<core::Size,2>::const_iterator jjatom=jatoms.begin(); jjatom!=jatoms.end(); ++jjatom ) {
					Size const &iatom(*iiatom), &jatom(*jjatom);
					if ( iatom==0 || jatom==0 ) continue;

					core::Real const dist = pose.residue(ires).xyz(iatom).distance( pose.residue(jres).xyz(jatom) );
					if ( dist > max_distance_ ) continue;

					core::scoring::func::FuncOP apc_func; //NULL!
					if ( !use_bounded_func_ ) {
						core::scoring::func::FuncOP sog_func( new core::scoring::func::SOGFunc( dist, coord_dev_ ) );
						core::scoring::func::FuncOP weighted_func( new core::scoring::func::ScalarWeightedFunc( cst_weight_, sog_func ) );
						apc_func = weighted_func;
						TR.Debug << "Using SOG func with AtomPairConstraint: " << *apc_func << std::endl;
					} else {
						//establish bound_func limits based on basin around distance
						core::Real lower_limit( dist - ( 0.5  * bound_width_) );
						if ( lower_limit < 0 ) lower_limit = 0;
						core::Real const upper_limit( dist + ( 0.5 * bound_width_ ) );
						apc_func = core::scoring::func::FuncOP( new BoundFunc( lower_limit, upper_limit, coord_dev_, "xyz" ) );
						TR.Debug << "Using Bounded func with AtomPairConstraint: " << *apc_func << std::endl;
					}
					runtime_assert(apc_func); //make sure that func got populated; should have been forced down one of the two if branches above

					core::scoring::constraints::ConstraintCOP newcst(
						new core::scoring::constraints::AtomPairConstraint( core::id::AtomID( iatom, ires ), core::id::AtomID( jatom, jres ), apc_func ) );
					csts.push_back( newcst );
					TR.Debug << "atom_pair_constraint generated for residue " << ires << ", atom " << iatom << " and residue " << jres << ", atom " << jatom << " with weight " << cst_weight_ << std::endl;
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
	basic::datacache::DataMap & datamap,
	Filters_map const & filters,
	moves::Movers_map const & movers,
	Pose const & pose
) {
	use_distance_cst_ = tag->getOption< bool >( "use_distance_cst", use_distance_cst_ );
	max_distance_ = tag->getOption< core::Real >( "max_distance", max_distance_ );
	coord_dev_ = tag->getOption< core::Real >( "coord_dev", coord_dev_ );
	bound_width_ = tag->getOption< core::Real >( "bound_width", bound_width_ );
	min_seq_sep_ = tag->getOption< core::Size>( "min_seq_sep", min_seq_sep_ );
	cst_weight_ = tag->getOption< core::Real >( "cst_weight", cst_weight_ );
	CA_only_ = tag->getOption< bool >( "CA_only", CA_only_ );
	bb_only_ = tag->getOption< bool >( "bb_only", bb_only_ );
	inter_chain_ = tag->getOption< bool >( "inter_chain", inter_chain_ );

	if ( bound_width_ < 1e-3 ) {
		use_bounded_func_ = false;
	} else {
		use_bounded_func_ = true;
	}

	if ( tag->hasOption( "task_operations" ) ) {
		TR.Warning << "WARNING: task_operations only active for proteins" << std::endl;
		parse_task_operations( tag, datamap, filters, movers, pose );
	}

	core::select::residue_selector::ResidueSelectorCOP selector = protocols::rosetta_scripts::parse_residue_selector( tag, datamap );
	if ( selector ) residue_selector( selector );

	if ( tag->hasOption( "task_operations" ) && tag->hasOption( "residue_selector" ) ) {
		std::stringstream msg;
		msg << "AddConstraintsToCurrentConformationMover::parse_my_tag(): 'task_operations' and 'residue_selector' cannot both be specified." << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}
}

void
AddConstraintsToCurrentConformationMover::parse_task_operations(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const &
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
	return moves::MoverOP( new AddConstraintsToCurrentConformationMover( *this ) );
}
moves::MoverOP AddConstraintsToCurrentConformationMover::fresh_instance() const {
	return moves::MoverOP( new AddConstraintsToCurrentConformationMover );
}

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AddConstraintsToCurrentConformationMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new AddConstraintsToCurrentConformationMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AddConstraintsToCurrentConformationMoverCreator::keyname() const {
// XRW TEMP  return AddConstraintsToCurrentConformationMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AddConstraintsToCurrentConformationMover::mover_name() {
// XRW TEMP  return "AddConstraintsToCurrentConformationMover";
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AddConstraintsToCurrentConformationMover::get_name() const {
// XRW TEMP  return "AddConstraintsToCurrentConformationMover";
// XRW TEMP }

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
		+ XMLSchemaAttribute( "bound_width", xsct_real, "BoundFunc zero basin width BoundFunc; also activates use of BoundFunc over HarmonicFunc" )
		+ XMLSchemaAttribute( "min_seq_sep", xsct_non_negative_integer, "Do not generate distance constraints between residues within this sequence separation.  Only active with use_distance_cst." )
		+ XMLSchemaAttribute( "cst_weight", xsct_real, "use ScalarWeightedFunc to reweight constraints by this" )
		+ XMLSchemaAttribute( "CA_only", xsct_rosetta_bool, "constrain only CA atoms." )
		+ XMLSchemaAttribute( "bb_only", xsct_rosetta_bool, "constrain only backbone atoms." )
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
	return protocols::moves::MoverOP( new AddConstraintsToCurrentConformationMover );
}

void AddConstraintsToCurrentConformationMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddConstraintsToCurrentConformationMover::provide_xml_schema( xsd );
}


} // moves
} // protocols
