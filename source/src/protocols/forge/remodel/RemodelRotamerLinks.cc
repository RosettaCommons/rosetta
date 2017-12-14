// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <protocols/forge/remodel/RemodelRotamerLinks.hh>
#include <protocols/forge/remodel/RemodelRotamerLinksCreator.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>

#include <core/chemical/ResidueType.hh>
#include <basic/options/option.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// option key includes
#include <basic/options/keys/remodel.OptionKeys.gen.hh>

using namespace basic::options;

namespace protocols {
namespace forge {
namespace remodel {

using namespace core;
using namespace chemical;
using namespace conformation;
using namespace basic::options;
using namespace pack;
using namespace rotamer_set;
using namespace task;
using namespace operation;
using namespace utility::tag;

using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static basic::Tracer TR( "protocols.forge.remodel.RemodelRotamerLinks", t_info );

TaskOperationOP RemodelRotamerLinksCreator::create_task_operation() const
{
	return TaskOperationOP( new RemodelRotamerLinks );
}

void RemodelRotamerLinksCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RemodelRotamerLinks::provide_xml_schema( xsd );
}

std::string RemodelRotamerLinksCreator::keyname() const
{
	return RemodelRotamerLinks::keyname();
}

RemodelRotamerLinks::~RemodelRotamerLinks() = default;

TaskOperationOP RemodelRotamerLinks::clone() const
{
	return TaskOperationOP( new RemodelRotamerLinks( *this ) );
}

RemodelRotamerLinks::RemodelRotamerLinks() = default;

void
RemodelRotamerLinks::parse_tag( TagCOP /*tag*/ , DataMap & )
{}

void RemodelRotamerLinks::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname(), "XRW TO DO: no documentation online, setup RemodelRotamerLinks." );
}

void
RemodelRotamerLinks::apply(
	Pose const & pose,
	PackerTask & ptask
) const
{
	Size nres;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
		nres = symm_info->num_independent_residues();
	} else {
		nres = pose.size();
	}

	if ( nres == 0 ) { utility_exit_with_message("Unable to setup RemodelRotamerLinks on an empty pose."); }

	// setup residue couplings
	RotamerLinksOP links( new RotamerLinks );
	links->resize( nres );

	Size repeat_number =option[OptionKeys::remodel::repeat_structure];
	Size segment_length = nres / repeat_number;
	if ( segment_length == 0 ) { utility_exit_with_message("Unable to setup RemodelRotamerLinks with zero-length segments."); }

	utility::vector1< utility::vector1< Size > > equiv_pos;

	//find all the equivalent positions, first pass iterate over the base
	for ( Size res = 1; res<= segment_length ; res++ ) {
		utility::vector1< Size> list;

		for ( Size rep = 0; rep < repeat_number; rep++ ) {
			list.push_back(res+(segment_length*rep));
		}
		equiv_pos.push_back(list);
	}


	//second pass, iterate over to populate the entire chain

	for ( Size i = 1; i <= nres ; i++ ) {
		Size subcounter = (i%segment_length);
		if ( subcounter == 0 ) {
			subcounter = segment_length;
		}

		links->set_equiv(i, equiv_pos[subcounter]);

		//std::cout << "linking " << i << " with " << subcounter << "array with ";
		for ( Size k=1; k<= equiv_pos[subcounter].size(); k++ ) {
			//std::cout << " " << equiv_pos[subcounter][k];
		}
		//std::cout << std::endl;
		// check for similarities
		//std::cout << ptask.task_string(pose);

	}

	//std::cout << ptask << std::endl;


	ptask.rotamer_links( links );
}

} // namespace remodel
} // namespace forge
} // namespace protocols

