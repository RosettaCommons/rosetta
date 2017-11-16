// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/abinitio/DomainAssembly.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/abinitio/DomainAssembly.hh>
#include <protocols/abinitio/DomainAssemblyCreator.hh>

// Package headers

// Project headers
#include <core/fragment/FragSet.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <core/pack/pack_rotamers.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/abinitio/AbrelaxApplication.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <utility>
#include <utility/tag/Tag.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

//Auto Headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Residue.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace abinitio {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.DomainAssembly" );

// XRW TEMP std::string
// XRW TEMP DomainAssemblyCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return DomainAssembly::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP DomainAssemblyCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new DomainAssembly );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP DomainAssembly::mover_name()
// XRW TEMP {
// XRW TEMP  return "DomainAssembly";
// XRW TEMP }

DomainAssembly::DomainAssembly() :
	protocols::moves::Mover( DomainAssembly::mover_name() )
{}

DomainAssembly::DomainAssembly(
	std::string const & linker_start,
	std::string const & linker_end,
	FragSetOP fragset_large,
	FragSetOP fragset_small
) :
	protocols::moves::Mover( DomainAssembly::mover_name() ),
	linker_start_( linker_start ),
	linker_end_( linker_end ),
	fragset_large_(std::move( fragset_large )),
	fragset_small_(std::move( fragset_small )),
	fragments_set_( true )
{}

DomainAssembly::~DomainAssembly() = default;

protocols::moves::MoverOP
DomainAssembly::clone() const {
	return( protocols::moves::MoverOP( new DomainAssembly( *this ) ) );
}

protocols::moves::MoverOP
DomainAssembly::fresh_instance() const
{
	return protocols::moves::MoverOP( new DomainAssembly );
}

void
DomainAssembly::apply( core::pose::Pose & pose )
{
	protocols::simple_moves::ReturnSidechainMover recover_sidechains( pose );

	if ( !fragments_set_ ) {
		TR<<"*******WARNING WARNING********: fragments not set, skipping domain assembly"<<std::endl;
		runtime_assert( fragments_set_ );
		return;
	}

	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap() );
	mm->set_bb( false );
	core::Size linker_start( core::pose::parse_resnum( linker_start_, pose ) );
	core::Size linker_end( core::pose::parse_resnum( linker_end_, pose ) );
	runtime_assert( linker_end > linker_start );
	runtime_assert( linker_start > 0 );
	runtime_assert( linker_end < pose.size() );

	for ( core::Size i=linker_start; i<=linker_end; ++i ) mm->set_bb( i, true );
	/* //The following code only makes sense if we do 'blind' prediction stuff. In all cases where we have a starting structure,
	//it's probably a good idea to start near that.
	for(i = 1; i <= flexible_regions_.size(); i++ ){
	ir = std::max(0, std::min( (int)pose_full_centroid.size() ,  (int) flexible_regions_[i] ) );

	std::cout << "Linker: " << ir << std::endl;

	Real const init_phi  ( -150.0 );
	Real const init_psi  (  150.0 );
	Real const init_omega(  180.0 );
	pose_full_centroid.set_phi( ir ,  init_phi  );
	pose_full_centroid.set_psi( ir ,  init_psi  );
	pose_full_centroid.set_omega( ir, init_omega);
	movemap->set_bb(ir, true);
	}
	*/
	core::scoring::ScoreFunctionCOP scorefxn( get_score_function() );
	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( core::chemical::CENTROID );
	protocols::simple_moves::SwitchResidueTypeSetMover to_fullatom( core::chemical::FA_STANDARD );
	to_centroid.apply( pose );
	protocols::abinitio::AbrelaxApplication abrelax_app;
	do {
		protocols::abinitio::ClassicAbinitio abinit( fragset_large_, fragset_small_, mm );
		abinit.init( pose );
		abinit.apply( pose );
	} while( !abrelax_app.check_filters( pose ) );

	//recover sidechains from starting structures
	to_fullatom.apply( pose );
	recover_sidechains.apply( pose );
	// pose.update_residue_neighbors(); // o/w fails assertion `graph_state_ == GOOD`
	// (*scorefxn)( pose );
	// scorefxn->accumulate_residue_total_energies( pose );

	//Repack regions around the linker
	core::pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().or_include_current( true );
	task->restrict_to_repacking();

	for ( core::Size i = linker_end+1; i <= pose.size(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;
		if ( pose.residue(i).type().is_disulfide_bonded() ) {
			task->nonconst_residue_task( i ).prevent_repacking();
			continue;
		}

		core::conformation::Residue const resi( pose.residue( i ) );
		core::Size j;
		for ( j = 1; j<=linker_end; ++j ) {
			core::conformation::Residue const resj( pose.residue( j ) );

			core::Real const distance( resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) );
			if ( distance <= 8.0 ) break;
		}
		if ( j>linker_end ) task->nonconst_residue_task( i ).prevent_repacking();
	}
	for ( core::Size i = 1; i <= linker_start - 1; ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;
		if ( pose.residue(i).type().is_disulfide_bonded() ) {
			task->nonconst_residue_task( i ).prevent_repacking();
			continue;
		}

		core::conformation::Residue const resi( pose.residue( i ) );
		core::Size j;
		for ( j = linker_start; j<=pose.size(); ++j ) {
			core::conformation::Residue const resj( pose.residue( j ) );

			core::Real const distance( resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) );
			if ( distance <= 8.0 ) break;
		}
		if ( j>pose.size() ) task->nonconst_residue_task( i ).prevent_repacking();
	}
	//in case there is a resfile, information in this resfile overrides the computed task
	if ( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) {
		core::pack::task::parse_resfile(pose, *task);
	}

	pack::pack_rotamers( pose, *scorefxn, task);
	(*scorefxn)( pose );
	/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( pose );
}

// XRW TEMP std::string
// XRW TEMP DomainAssembly::get_name() const {
// XRW TEMP  return DomainAssembly::mover_name();
// XRW TEMP }

void
DomainAssembly::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	linker_start_ = core::pose::get_resnum_string( tag, "linker_start_" );
	linker_end_   = core::pose::get_resnum_string( tag, "linker_end_" );

	std::string const frag_large_fname( tag->getOption< std::string >( "frag9", "frag9" ) );
	std::string const frag_small_fname( tag->getOption< std::string >( "frag3", "frag3" ) );

	using namespace core::fragment;
	using namespace basic::options;
	fragset_large_ = FragmentIO(option[ OptionKeys::abinitio::number_9mer_frags ] ).read_data( frag_large_fname );
	fragset_small_ = FragmentIO(option[ OptionKeys::abinitio::number_3mer_frags ] ).read_data( frag_small_fname );
	fragments_set_ = true;
}

std::string DomainAssembly::get_name() const {
	return mover_name();
}

std::string DomainAssembly::mover_name() {
	return "DomainAssembly";
}

void DomainAssembly::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "frag9", xs_string, "Path to fragment file containing 9mers", "frag9" )
		+ XMLSchemaAttribute::attribute_w_default( "frag3", xs_string, "Path to fragment file containing 3mers", "frag3" );

	core::pose::attributes_for_get_resnum_string( attlist, "linker_start_" );
	core::pose::attributes_for_get_resnum_string( attlist, "linker_end_" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Do domain assembly sampling by fragment insertion in a linker region", attlist );
}

std::string DomainAssemblyCreator::keyname() const {
	return DomainAssembly::mover_name();
}

protocols::moves::MoverOP
DomainAssemblyCreator::create_mover() const {
	return protocols::moves::MoverOP( new DomainAssembly );
}

void DomainAssemblyCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DomainAssembly::provide_xml_schema( xsd );
}


} // abinitio
} // protocols
