// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//////////////////////////////////////////////////////////////////////
/// @file
///
/// @brief
/// A general mover class for performing BCL mutates on Rosetta small molecule residues
///
/// @details
/// This class makes a mover for BCL small molecule drug design mutates.
/// The design leverages the BCL util::Implementation< template> class, which allows
/// the templated argument type to be constructed from an object data label, which
/// itself can be generated from a string. This is essentially how users interact
/// with the BCL from the command-line: they pass serializable arguments to the BCL
/// that are interpreted as object data labels for construction of objects within the
/// the specified application.
///
/// @author
/// Benjamin P. Brown (benjamin.p.brown17@gmail.com)
////////////////////////////////////////////////////////////////////////

// unit include
#include <protocols/drug_design/bcl/BCLFragmentMutateMover.hh>
#include <protocols/drug_design/bcl/BCLFragmentMutateMoverCreator.hh>

// rosetta scripts
#include <protocols/rosetta_scripts/util.hh>

// protocols includes
#include <protocols/drug_design/bcl/BCLSampleConfsManager.hh>

// projects includes
#include <core/chemical/Atom.hh>
#include <core/chemical/Atom.fwd.hh>
#include <core/chemical/AtomProperty.hh>
#include <core/chemical/AtomProperties.fwd.hh>
#include <core/chemical/bcl/util.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/Bond.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>

// utility includes
#include <utility/numbers.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ includes
#include <string>

// debug
#include <utility/io/ozstream.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

// BCL includes
#ifdef USEBCL
#include <bcl/include/chemistry/bcl_chemistry_atom_complete.h>
#include <bcl/include/chemistry/bcl_chemistry_atom_vector.h>
#include <bcl/include/chemistry/bcl_chemistry_fragment_track_mutable_atoms.h>
#include <bcl/include/chemistry/bcl_chemistry_rotamer_library_file.h>
#include <bcl/include/io/bcl_io_ofstream.h>
#include <bcl/include/io/bcl_io_file.h>
#include <bcl/include/linal/bcl_linal_vector_3d.h>
#include <bcl/include/storage/bcl_storage_vector.h>
#include <bcl/include/util/bcl_util_format.h>
#endif

namespace protocols
{
namespace drug_design
{
namespace bcl
{
// tracer
static basic::Tracer TR("protocols.drug_design.bcl.BCLFragmentMutateMover");

/////////////////
// static data //
/////////////////

//////////////////////////////////
// construction and destruction //
//////////////////////////////////

//! @brief default constructor
BCLFragmentMutateMover::BCLFragmentMutateMover() :
	BCLFragmentBaseMover(),
	max_mutate_retry_( 10)
{
#ifdef USEBCL
	mutate_ = ::bcl::util::Implementation< ::bcl::chemistry::FragmentMutateInterface>();
#endif
}

//! @brief default sample confs constructor
#ifdef USEBCL
BCLFragmentMutateMover::BCLFragmentMutateMover
(
	const core::Size &n_max_mutate_attempts,
	const ::bcl::util::Implementation< ::bcl::chemistry::FragmentMutateInterface> &mutate
) :
	BCLFragmentBaseMover(),
	max_mutate_retry_( n_max_mutate_attempts),
	mutate_( mutate)
{
	// extras=bcl required for construction of this class object
	core::chemical::bcl::require_bcl();
}
#endif

//! @brief clone this object
protocols::moves::MoverOP
BCLFragmentMutateMover::clone() const
{
	return utility::pointer::make_shared< BCLFragmentMutateMover >( *this );
}

//! @brief create this type of object
protocols::moves::MoverOP
BCLFragmentMutateMover::fresh_instance() const
{
	return utility::pointer::make_shared< BCLFragmentMutateMover >();
}

/////////////////
// data access //
/////////////////

std::string BCLFragmentMutateMover::mover_name()
{
	return "BCLFragmentMutateMover";
}

std::string BCLFragmentMutateMover::get_name() const
{
	return mover_name();
}

std::string BCLFragmentMutateMover::get_tag( Pose & pose, core::Size resid ) const
{
	core::pose::PDBInfoCOP pdb_info( pose.pdb_info() );
	std::stringstream tag;
	tag << pdb_info->chain( resid) << pdb_info->number( resid);
	if ( pdb_info->icode( resid) && pdb_info->icode( resid) != ' ' ) {
		tag << pdb_info->icode( resid);
	}
	return tag.str();
}

#ifdef USEBCL
//! @brief return the type of mutate being performed
std::string const &BCLFragmentMutateMover::get_mutate_label() const
{
	return mutate_->GetAlias();
}
#else
void BCLFragmentMutateMover::get_mutate_label() const
{
	core::chemical::bcl::require_bcl();
}
#endif

////////////////
// operations //
////////////////

//! @brief apply BCLFragmentMutateMover
//! @details extract ligand residue from pose, apply mutate,
//! add mutated ligand back to the pose; additional protein-ligand
//! interface optimization NOT included; I recommend repacking/minimization
void BCLFragmentMutateMover::apply(
#ifdef USEBCL
	core::pose::Pose & pose )
#else
	core::pose::Pose & /* pose */ )
#endif
{
	// extras=bcl required for construction of this class object
	core::chemical::bcl::require_bcl();

#ifdef USEBCL

	// restrict chain  ID to single letter
	if ( get_base_ligand_chain().size() != 1 ) {
		utility_exit_with_message
			(
			"BCLFragmentMutateMover chain ID '" + get_base_ligand_chain() + "' must be a single letter."
		);
	}

	// get ligand identifiers
	utility::vector1< core::Size> ligand_resnums( core::pose::get_resnums_for_chain( pose, get_base_ligand_chain()[0] ));

	// require that we actually have a ligand specified by the selected resnum
	if ( !ligand_resnums.size() ) {
		utility_exit_with_message( "The provided chain ID does not correspond to a residue!");
	}
	// restrict chain to 1 residue
	if ( ligand_resnums.size() > 1 ) {
		utility_exit_with_message( "BCLFragmentMutateMover can only operate on chains of size 1");
	}

	// output for user
	core::Size ligand_resnum( ligand_resnums[ 1]);
	TR << "Ligand chain ID: " + get_base_ligand_chain() << std::endl;
	TR << "Ligand resnum: " + ::bcl::util::Format()( ligand_resnum) << std::endl;

	// create a BCL fragment from the ligand in the pose
	::bcl::chemistry::FragmentComplete fragment( this->pose_residue_to_fragment( pose, ligand_resnum));

	// debug
	TR << this->get_mutate_label() << std::endl;

	// mutate
	// not all mutates work - molecules can fail the drug likeness filter, not have a 3D conformer
	// that can be easily generated, or not find a valid mutate to perform given the constraints
	// on mutable atoms and whatnot. depending on which mutate object you construct and
	// how you construct it, the mutate may also be highly stochastic. so, it is worth
	// trying a few times perhaps, but you also want to bail eventually.
	bool success( false);
	core::Size num_tries( 1);
	while ( !success ) {
		TR << "Mutate trial: " << num_tries << std::endl;
		::bcl::math::MutateResult< ::bcl::chemistry::FragmentComplete> mutated_fragment( ( *mutate_)( fragment));
		if ( num_tries == max_mutate_retry_ ) {
			TR << "Maximum mutate trials attempted. Returning input fragment!" << std::endl;
			break;
		} else if ( !mutated_fragment.GetArgument().IsDefined() ) {
			TR << "Mutate fragment returned null!" << std::endl;
			++num_tries;
			continue;
		} else {
			success = true;
			fragment = *( mutated_fragment.GetArgument());
		}
	}

	// can skip this last bit if we do not actually make a structural perturbation to our molecule
	if ( success ) {
		// generate conformers of the mutated fragment
		TR << "SampleByParts indices: " + ::bcl::util::Format()( fragment.GetMDLProperty("SampleByParts")) << std::endl;
		TR << "Generating conformers of the mutated fragment" << std::endl;
		::bcl::chemistry::FragmentEnsemble fragment_confs( BCLSampleConfsManager::get_instance()->get_sample_confs()( fragment).First());

		// revert to Rosetta type
		core::chemical::MutableResidueTypeOP reverted_res_ptr( get_base_handler().fragment_to_restype( fragment, fragment_confs));
		core::chemical::ResidueTypeCOP base_rosetta_type_ptr( core::chemical::ResidueType::make(( *reverted_res_ptr)));
		core::conformation::Residue new_ligand( base_rosetta_type_ptr, true);

		// set new residue in pose
		pose.replace_residue( ligand_resnum, new_ligand, false); // orient backbone bool NA for ligands
		pose.pdb_info()->set_resinfo( pose.total_residue(), get_base_ligand_chain()[0] , pose.total_residue(), ' ');
	}
#else
	core::chemical::bcl::require_bcl();
#endif
}

//! @brief set the maximum number of mutate attempts
void BCLFragmentMutateMover::set_n_max_mutate_attempts( core::Size const &n_max_mutate_attempts)
{
	max_mutate_retry_ = n_max_mutate_attempts;
}

//! @brief set the mutate using a BCL object data label string
void BCLFragmentMutateMover::set_mutate(
#ifdef USEBCL
	std::string const & object_data_label )
#else
	std::string const & /* object_data_label */ )
#endif
{
#ifdef USEBCL
	// make the initial object
	mutate_ = ::bcl::util::Implementation<::bcl::chemistry::FragmentMutateInterface >( object_data_label);
#else
	core::chemical::bcl::require_bcl();
#endif
}

//////////////////////
// helper functions //
//////////////////////

//! @brief parse xml file
void BCLFragmentMutateMover::parse_my_tag( TagCOP tag, basic::datacache::DataMap &data )
{
	// extras=bcl required for construction of this class object
	core::chemical::bcl::require_bcl();

	// parse the parent class tag
	BCLFragmentBaseMover::parse_my_tag( tag, data);

	// construct the BCL mutate implementation with the object data label string
	this->set_mutate( tag->getOption< std::string>( "object_data_label"));

	// set max number of mutate attempts before returning the input ligand
	this->set_n_max_mutate_attempts( tag->getOption< core::Size>( "n_max_mutate_attempts", max_mutate_retry_));
}

//! @brief set xml options
void BCLFragmentMutateMover::provide_xml_schema( utility::tag::XMLSchemaDefinition &xsd )
{
	// initialize attribute list
	using namespace utility::tag;
	AttributeList attlist;

	// tack on parent class attributes
	BCLFragmentBaseMover::base_attributes( attlist);

	// load child class attributes
	attlist + XMLSchemaAttribute::required_attribute
		(
		"object_data_label", xs_string,
		"The BCL Object Data Label specifying the mutate to use and the options "
		"with which that mutate will be constructed. To view all of the options "
		"available to a mutate, pass 'help' as the argument to the mutate. For "
		"example, if you are making an ExtendWithLinker flavor of this mover, "
		"in your XML script Mover definition for BCLFragmentMutateMover, the "
		"'object_data_label' attribute would read 'ExtendWithLinker(help)'."
	);

	attlist + XMLSchemaAttribute::attribute_w_default
		(
		"n_max_mutate_attempts", xsct_non_negative_integer,
		"The maximum number of retry attempts that can occur for any single "
		"mutate move. Retry attempts are initiated when a mutated fragment "
		"fails the internal druglikeness filter or when 3D coordinates cannot "
		"be generated for the mutated fragment.",
		"10"
	);

	// Mover description and final attlist
	protocols::moves::xsd_type_definition_w_attributes
		(
		xsd, mover_name(),
		"Perform small molecule drug design mutations via the BCL FragmentMutateInterface",
		attlist
	);
}

//////////////////////////////
// mover creator operations //
//////////////////////////////

std::string
BCLFragmentMutateMoverCreator::keyname() const
{
	return BCLFragmentMutateMover::mover_name();
}

protocols::moves::MoverOP
BCLFragmentMutateMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< BCLFragmentMutateMover >();
}

void
BCLFragmentMutateMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition &xsd ) const
{
	BCLFragmentMutateMover::provide_xml_schema( xsd );
}

//! @brief Provide citations to the passed CitationCollectionList.
//! @details This mover is unpublished. Hopefully we can publish it soon and I can change
//! this to a citation for an actual manuscript
void
BCLFragmentMutateMover::provide_citation_info
(
	basic::citation_manager::CitationCollectionList &citations
) const
{
	using namespace basic::citation_manager;
	UnpublishedModuleInfoOP mover_info( utility::pointer::make_shared< UnpublishedModuleInfo>
		(
		"BCLFragmentMutateMover",
		CitedModuleType::Mover
		));

	// add me
	mover_info->add_author
		(
		"Benjamin P. Brown", "Vanderbilt University", "benjamin.p.brown17@gmail.com"
	);

	// add Jeff
	mover_info->add_author
		(
		"Jeffrey Mendenhall", "Vanderbilt University", "jeffreymendenhall@gmail.com"
	);

	// add Jens the bossman
	mover_info->add_author
		(
		"Jens Meiler", "Vanderbilt University", "jens@meilerlab.org"
	);

	citations.add( mover_info);
}

} // bcl
} // drug_design
} // protocols
