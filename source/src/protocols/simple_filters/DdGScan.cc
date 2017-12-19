// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/DdGScan.cc
/// @brief
/// @author Neil King (neilking@uw.edu)
/// @author Kyle Barlow (kb@kylebarlow.com)

// Project Headers
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/types.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <map>
#include <numeric/random/random.hh>
#include <protocols/jd2/util.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/simple_filters/DdgFilter.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_moves/ddG.hh>
//#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <string>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <protocols/simple_filters/DdGScan.hh>
#include <protocols/simple_filters/DdGScanCreator.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

static basic::Tracer TR( "protocols.simple_filters.DdGScan" );

namespace protocols {
namespace simple_filters {

/// @brief default constructor
DdGScan::DdGScan():
	task_factory_( /* NULL */ ),
	repeats_( 3 ),
	scorefxn_( /* NULL */ ),
	report_diffs_( true ),
	write2pdb_( false )
{
	initialize();
}

/// @brief constructor with arguments
DdGScan::DdGScan(
	core::pack::task::TaskFactoryOP task_factory,
	core::Size const repeats,
	core::scoring::ScoreFunctionCOP /*scorefxn*/,
	bool report_diffs,
	bool write2pdb
):
	Filter( "DdGScan" ),
	task_factory_(std::move( task_factory )),
	repeats_( repeats ),
	report_diffs_( report_diffs ),
	write2pdb_( write2pdb )
{
	initialize();
}

/// @brief copy constructor
DdGScan::DdGScan( DdGScan const & )= default;

/// @brief destructor
DdGScan::~DdGScan() = default;

void DdGScan::initialize() {
	ddG_mover( protocols::simple_moves::ddGOP( new protocols::simple_moves::ddG() ) );
}

protocols::filters::FilterOP DdGScan::clone() const { return protocols::filters::FilterOP( new DdGScan( *this ) ); }
protocols::filters::FilterOP DdGScan::fresh_instance() const { return protocols::filters::FilterOP( new DdGScan() ); }

// setters
void DdGScan::task_factory( core::pack::task::TaskFactoryOP task_factory ) { task_factory_ = task_factory; }
void DdGScan::repeats( core::Size const r ) { repeats_ = r; }
void DdGScan::scorefxn( core::scoring::ScoreFunctionOP const scorefxn ) { scorefxn_ = scorefxn; }
void DdGScan::report_diffs( bool const report_diffs ) { report_diffs_ = report_diffs; }
void DdGScan::write2pdb( bool const write ) { write2pdb_ = write; }
void DdGScan::ddG_mover( protocols::simple_moves::ddGOP const ddG_mover_op ) {
	ddG_mover_ = ddG_mover_op ;
}

// getters
core::pack::task::TaskFactoryOP DdGScan::task_factory() const { return task_factory_; }
core::Size DdGScan::repeats() const { return repeats_; }
bool DdGScan::report_diffs() const { return report_diffs_; }
bool DdGScan::write2pdb() const { return write2pdb_; }
protocols::simple_moves::ddGOP DdGScan::ddG_mover() const { return ddG_mover_; }

/// @brief Dummy Filter apply function
bool DdGScan::apply( core::pose::Pose const & ) const { return true; }

/// @brief parse xml
void
DdGScan::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & /*pose*/
)
{
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	repeats( tag->getOption< core::Size >( "repeats", 1 ) );
	report_diffs( tag->getOption< bool >( "report_diffs", true ) );
	write2pdb( tag->getOption< bool >( "write2pdb", false ) );
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );

	// Handle definition of a special ddG mover
	if ( tag->hasOption("ddG_mover") ) {
		moves::MoverOP mover = protocols::rosetta_scripts::parse_mover( tag->getOption< std::string >( "ddG_mover" ), movers );
		ddG_mover( utility::pointer::dynamic_pointer_cast < protocols::simple_moves::ddG > (mover));
	}
}

/// @brief Calculate the binding ddG for a mutation at the specified position (or wildtype no-mutant if resi is 0)
core::Real
DdGScan::ddG_for_single_residue( core::pose::Pose const & const_pose, core::Size const resi, core::pack::task::PackerTaskOP const general_task, core::pose::Pose & pose_to_mutate ) const
{
	// Setup packer task to only mutate single residue (and not pack or design anything else)
	using namespace core::pack::task;
	PackerTaskOP task = general_task->clone();
	task->initialize_from_command_line().or_include_current( true );
	for ( core::Size resj=1; resj<=pose_to_mutate.size(); ++resj ) {
		if ( resi != resj ) {
			task->nonconst_residue_task( resj ).prevent_repacking();
		}
	}

	core::scoring::ScoreFunctionOP symmetry_ready_scorefxn;
	if ( core::pose::symmetry::is_symmetric( pose_to_mutate ) ) {
		// Except for score function conversion, symmetry is now handled inline (pack_rotamers does autodispatch).
		symmetry_ready_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *scorefxn_ );
	} else {
		symmetry_ready_scorefxn = scorefxn_;
	}

	core::pack::pack_rotamers( pose_to_mutate, *symmetry_ready_scorefxn, task );

	// Set options for ddG mover, calculate, and return value
	ddG_mover_->scorefxn( symmetry_ready_scorefxn );

	core::Real average( 0.0 );
	for ( core::Size i=1; i<=repeats(); ++i ) {
		ddG_mover_->calculate( pose_to_mutate );
		average += ddG_mover_->sum_ddG();
		ddG_mover_->report_ddG( TR );
	}
	core::Real const calculated_ddG( average / (core::Real)repeats() );

	if ( resi == 0 ) {
		TR << protocols::jd2::current_output_name() << " wild-type binding ddG " <<
			" = " << calculated_ddG << std::endl;
	} else {
		TR << protocols::jd2::current_output_name() << " DdGScan binding ddG for mutation " <<
			const_pose.residue( resi ).name3() << resi << pose_to_mutate.residue( resi ).name3() << " = " << calculated_ddG << std::endl;
	}

	TR.flush();
	return( calculated_ddG );
}

/// @brief Implementation of virtual report function
/// @details Calls report function with return value
void
DdGScan::report( std::ostream & out, core::pose::Pose const & const_pose ) const
{
	calculate( out, const_pose );
}

/// @brief calculate and report the per-residue ddGs
utility::vector1< ddG_data_tuple >
DdGScan::calculate( std::ostream & out, core::pose::Pose const & const_pose ) const
{

	// Set up, get data
	core::pose::Pose pose( const_pose );
	// Used to save data for return function
	utility::vector1< ddG_data_tuple > ddg_saved_data;

	// Apply TaskOperations from xml
	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
	if ( task_factory_ != nullptr ) {
		task = task_factory_->create_task_and_apply_taskoperations( pose );
	} else {
		TR.Warning << "You have not provided any TaskOperations. A default will be used." << std::endl;
	}

	// Calculate wt ddG score of binding (by using special case of residue "0", handled by function called)
	core::Real const wt_ddG = ddG_for_single_residue( const_pose, 0, task, pose );

	// *** Entire loop goal: Calculate the binding ddG score upon mutation ***
	// Loop through residues
	for ( core::Size resi = 1; resi <= pose.size(); resi++ ) {
		// Check for canonical protein residue and that position is designable to something
		// Here we check only to see if the residue is set to at least be packable
		if ( pose.residue( resi ).is_protein() && task->pack_residue( resi ) ) {

			std::string const res_type( const_pose.residue( resi ).name3() );
			core::pose::PDBInfoCOP pose_info( const_pose.pdb_info() );
			char const chain( pose_info->chain( resi ) );
			core::Size output_resi = resi;
			if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
				output_resi = pose.pdb_info()->number( resi );
			}

			// Loop through designable residues at this position (AKA make packer do each one one at a time instead of try all at once)
			// allowed_residue_types only contain residues defined through PIKAA - so a NATRO-only residue that made it this far won't have nything else happen
			for (
					auto aa_iter = task->nonconst_residue_task( resi ).allowed_residue_types_begin();
					aa_iter != task->nonconst_residue_task( resi ).allowed_residue_types_end();
					++aa_iter
					) {
				pose = core::pose::Pose( const_pose );
				core::pack::task::PackerTaskOP point_mutant_task( task->clone() );
				utility::vector1<bool> restrict_to_aa(20, false);
				restrict_to_aa[ (*aa_iter)->aa() ] = true;
				point_mutant_task->nonconst_residue_task(resi).restrict_absent_canonical_aas(restrict_to_aa);

				core::Real const mut_ddG( ddG_for_single_residue( const_pose, resi, point_mutant_task, pose ) );
				core::Real const diff_ddG( mut_ddG - wt_ddG );
				std::string const mutant_res_type( pose.residue( resi ).name3() );
				core::Real output_ddG = ( report_diffs() == 1 ) ? diff_ddG : mut_ddG;

				// Output to pdb file
				if ( write2pdb() ) {
					write_to_pdb( pose, resi, output_resi, res_type, output_ddG );
				}

				// Save output data for return value (for features database, among other things)
				ddg_saved_data.push_back( ddG_data_tuple( resi, mutant_res_type, output_ddG ) );

				// Output to tracer
				out << " " << "Residue " << chain << output_resi << res_type << "->" <<  mutant_res_type << " : " << ObjexxFCL::format::F (9,4,output_ddG) << std::endl;

			} // End loop through all designable residues
		} // Close if statement
	} // Close outer for loop

	out << std::endl;
	return ddg_saved_data;
}

void DdGScan::write_to_pdb(
	core::pose::Pose const & pose,
	core::Size const & residue,
	core::Size const & output_resi,
	std::string const & residue_name,
	core::Real const & ddG
) const {

	std::string filter_name = this->name();
	std::string user_name = this->get_user_defined_name();

	std::string output_string = filter_name + " " + user_name + ": " + residue_name + ObjexxFCL::string_of(output_resi) + pose.residue( residue ).name3() + " = " + ObjexxFCL::format::F (9,4,ddG);
	protocols::jd2::add_string_to_current_job(output_string);

}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP DdGScanCreator::create_filter() const { return protocols::filters::FilterOP( new DdGScan ); }

// XRW TEMP std::string
// XRW TEMP DdGScanCreator::keyname() const { return "DdGScan"; }

std::string DdGScan::name() const {
	return class_name();
}

std::string DdGScan::class_name() {
	return "DdGScan";
}

void DdGScan::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "repeats" , xsct_non_negative_integer , "How many times to repeat the ddg calculations; the average of all the repeats is returned." , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "report_diffs" , xsct_rosetta_bool , "Whether to report the changes in binding energy upon mutation (pass true), or the total binding energy for the mutated structure (pass false)." , "1" )
		+ XMLSchemaAttribute( "ddG_mover" , xs_string, "Handle definition of a special ddG mover." ) //XRW TO DO: should this be an attributes_for_parse_mover??
		+ XMLSchemaAttribute::attribute_w_default( "write2pdb" , xsct_rosetta_bool , "Whether to write the residue-specific ddG information to the output .pdb file." , "0" ) ;

	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist ) ;
	//The score function used for the calculations. If a ddG mover is defined, this score function will be used in that mover as well, overriding any different score function defined in that mover.
	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist ) ;
	//The task operations to use to identify which residues to scan. Designable or packable residues are scanned.

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Takes a set of task operations from the user in order to more precisely specify a set of residues to analyze via ddG scanning. Individually mutates each of the residues to alanine (or whatever other residue is defined in the task operations) and calculates the change in binding energy (ddG).", attlist );
}

std::string DdGScanCreator::keyname() const {
	return DdGScan::class_name();
}

protocols::filters::FilterOP
DdGScanCreator::create_filter() const {
	return protocols::filters::FilterOP( new DdGScan );
}

void DdGScanCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DdGScan::provide_xml_schema( xsd );
}


} // simple_filters
} // protocols
