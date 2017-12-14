// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/SingleFragmentMover.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/nonlocal/SingleFragmentMover.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/pose/Pose.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>

// Package headers
#include <protocols/nonlocal/Chunk.hh>
#include <protocols/nonlocal/Policy.hh>
#include <protocols/nonlocal/PolicyFactory.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/nonlocal/SingleFragmentMoverCreator.hh>

using Parent = protocols::moves::Mover;

namespace protocols {
namespace nonlocal {

static basic::Tracer TR( "protocols.nonlocal.SingleFragmentMover" );


SingleFragmentMover::SingleFragmentMover() : Parent("SingleFragmentMover") {}

SingleFragmentMover::SingleFragmentMover(const core::fragment::FragSetOP& fragments,
	const core::kinematics::MoveMapOP& movable)
: Parent("SingleFragmentMover") {
	initialize(fragments, movable, PolicyFactory::get_policy("uniform", fragments));
}

SingleFragmentMover::SingleFragmentMover(const core::fragment::FragSetOP& fragments,
	const core::kinematics::MoveMapOP& movable,
	const PolicyOP& policy)
: Parent("SingleFragmentMover") {
	initialize(fragments, movable, policy);
}

void SingleFragmentMover::initialize(const core::fragment::FragSetOP& fragments,
	const core::kinematics::MoveMapOP& movable,
	const PolicyOP& policy) {
	debug_assert(fragments);
	//debug_assert(movable); // Can be nullptr
	debug_assert(policy);

	// Initialize member variables
	fragments_ = fragments;
	movable_ = movable;
	policy_ = policy;

	// Create a position-indexable lookup for core::fragment::Frame's
	initialize_library();
}

void SingleFragmentMover::apply(core::pose::Pose& pose) {
	using core::Size;
	using core::fragment::FragDataCOP;
	using core::fragment::Frame;
	using core::kinematics::FoldTree;

	// ensure that preconditions on the input have been met
	debug_assert(pose.size() > 0);
	debug_assert(pose.fold_tree().check_fold_tree());
	debug_assert(valid());

	// determine whether the pose is in fullatom mode. if so, warn the user and
	// convert it to centroid mode automatically.
	// bool was_fullatom = // Unused variable causes warning.
	to_centroid(&pose);

	MoveMapOP movable( movemap( pose ) );

	// reuse <chunks_> when possible
	const FoldTree& current_tree = pose.fold_tree();
	if ( !previous_tree_ || *previous_tree_ != current_tree ) {
		chunks_.clear();
		probs_.clear();
		initialize_chunks(current_tree, movable);
		previous_tree_ = FoldTreeOP( new core::kinematics::FoldTree(current_tree) );
	}

	// randomly select the insertion position
	const Chunk* chunk = random_chunk();
	if ( !chunk ) {
		TR.Warning << "No move possible-- 0 chunks" << std::endl;
		return;
	}

	Size insertion_pos = chunk->choose();

	while ( !movable->get_bb(insertion_pos) )
			insertion_pos = chunk->choose();

	// delegate responsibility for choosing the fragment to the policy
	const Frame& frame = library_[insertion_pos];
	Size fragment_pos = policy_->choose(frame, pose);
	frame.apply(*movable, fragment_pos, pose);
	TR.Debug << "Inserted fragment at position " << insertion_pos << std::endl;
}

// XRW TEMP std::string SingleFragmentMover::get_name() const {
// XRW TEMP  return "SingleFragmentMover";
// XRW TEMP }

protocols::moves::MoverOP SingleFragmentMover::clone() const {
	return protocols::moves::MoverOP( new SingleFragmentMover(*this) );
}

protocols::moves::MoverOP SingleFragmentMover::fresh_instance() const {
	return protocols::moves::MoverOP( new SingleFragmentMover() );
}

void SingleFragmentMover::parse_my_tag(const utility::tag::TagCOP tag,
	basic::datacache::DataMap& data,
	const protocols::filters::Filters_map&,
	const protocols::moves::Movers_map&,
	const core::pose::Pose& ) {
	using core::fragment::FragmentIO;
	using core::fragment::FragSetOP;
	using core::kinematics::MoveMap;
	using core::kinematics::MoveMapOP;
	using std::string;

	// required flags
	// fragments, movemap
	string fragments_file = tag->getOption<string>("fragments");
	FragSetOP fragments = FragmentIO().read_data(fragments_file);

	// initially, all backbone torsions are movable
	movemap_factory( protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data ) );

	// optional flags
	// policy -- default => uniform
	string policy_type = tag->getOption<string>("policy", "uniform");
	PolicyOP policy = PolicyFactory::get_policy(policy_type, fragments);

	// Movable is set by the movemap_factory
	initialize(fragments, /*movable=*/nullptr, policy);
}

bool SingleFragmentMover::valid() const {
	return fragments_ && ( movable_ || movemap_factory_ ) && policy_;
}

void SingleFragmentMover::initialize_library() {
	core::fragment::ConstFrameIterator i;
	for ( i = fragments_->begin(); i != fragments_->end(); ++i ) {
		Size position = (*i)->start();
		library_[position] = **i;
	}
}

core::kinematics::MoveMapOP
SingleFragmentMover::movemap( core::pose::Pose const & pose ) {
	if ( movable_ ) {
		return movable_;
	} else if ( movemap_factory_ ) {
		return movemap_factory_->create_movemap_from_pose( pose );
	} else {
		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
		return movemap;
	}
}

void SingleFragmentMover::initialize_chunks(const core::kinematics::FoldTree& tree, MoveMapOP const & movable) {
	using core::Size;
	using core::Real;

	// Assumption: constant-length fragments at each insertion position
	Size fragment_len = fragments_->max_frag_length();

	// Last position to examine
	Size last_pos = fragments_->max_pos() - fragment_len + 1;

	Size p = fragments_->min_pos();
	do {
		Size q = p + 1;
		while ( !tree.is_cutpoint(q++) ) {}

		// During fold tree construction, it may happen that the distance between
		// the <p> and the next cutpoint is less than fragment length. In this case,
		// the Region::start() exceeds Region::stop(). It's impossible to perform
		// fragment insertions on this region given the current fragment library.
		RegionOP region( new Region(p, q - fragment_len) );
		if ( region->increasing() ) {
			// Ensure that the chunk is valid before adding it to the list. Mainly, this
			// means that there must be at least 1 movable residue.
			Chunk chunk(region, movable);

			if ( chunk.is_movable() ) {
				TR.Debug << "Added chunk: " << region->start() << "-" << region->stop() << std::endl;
				chunks_.push_back(chunk);
			} else {
				TR.Debug << "Skipped chunk: " << region->start() << "-" << region->stop() << ": no movable positions" << std::endl;
			}
		}

		p = q - 1;
	} while (p < last_pos);

	// Assign probabilities of selection to each Chunk according to its length
	Real N = fragments_->max_pos();
	Real sum = 0;
	for ( Size i = 1; i <= chunks_.size(); ++i ) {
		Real prob = chunks_[i].length() / N;
		probs_.push_back(prob);
		sum += prob;
	}

	// Normalize probabilities
	for ( Size i = 1; i <= probs_.size(); ++i ) {
		probs_[i] /= sum;
		TR << "P(c_" << i << ")=" << probs_[i] << std::endl;
	}
}

const Chunk* SingleFragmentMover::random_chunk() const {
	using core::Real;
	using core::Size;
	using utility::vector1;

	debug_assert(chunks_.size() > 0);
	vector1<Real> fitnesses(probs_);

	// convert to a CDF
	Size n = fitnesses.size();
	for ( Size i = 2; i <= n; ++i ) {
		fitnesses[i] += fitnesses[i-1];
	}

	// random number from 0 to f[n] (inclusive)
	Real x = numeric::random::uniform();

	// this could be done more efficiently with binary search
	for ( Size i = 2; i <= n; ++i ) {
		if ( fitnesses[i-1] < x && x <= fitnesses[i] ) {
			return &chunks_[i];
		}
	}
	return &chunks_[1];
}

bool SingleFragmentMover::to_centroid(core::pose::Pose* pose) const {
	debug_assert(pose);
	if ( pose->is_fullatom() ) {
		TR.Warning << "Input pose is full atom (centroid required)" << std::endl;
		TR.Warning << "Performing implicit conversion..." << std::endl;
		core::util::switch_to_residue_type_set(*pose, core::chemical::CENTROID_t );
		return true;
	}
	return false;
}

std::string SingleFragmentMover::get_name() const {
	return mover_name();
}

std::string SingleFragmentMover::mover_name() {
	return "SingleFragmentMover";
}

void SingleFragmentMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( [] (std::string const& name) {
		return "SingleFragmentMover_subelement_" + name + "Type";
	});
	rosetta_scripts::append_subelement_for_parse_movemap_factory_legacy(xsd, subelements);

	XMLSchemaRestriction policy_type_enum;
	policy_type_enum.name("policy_type");
	policy_type_enum.base_type( xs_string );
	policy_type_enum.add_restriction( xsr_enumeration, "uniform" );
	policy_type_enum.add_restriction( xsr_enumeration, "smooth" );
	xsd.add_top_level_element(policy_type_enum);

	AttributeList attlist;

	attlist + XMLSchemaAttribute::required_attribute(
		"fragments", xs_string,
		"Fragment file. Fragments are used in the assembly of proteins whether for structure prediction or design, to cut down on the size of the protein-folding search space. They are a core part of the Rosetta design. : Fragment libraries are used by many protocols but are a core part of ab initio.");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"policy", "policy_type",
		"Policy object is responsible for choosing from among the possible fragments contained in the fragment file. Currently, two policies are supported-- 'uniform' and 'smooth.' The former chooses uniformly amongst the set of possibilities. The latter chooses the fragment that, if applied, causes minimal distortion to the pose.",
		"uniform");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, mover_name(),
		"Performs a single fragment insertion move on the pose. Respects the restrictions imposed by the user-supplied MoveMap and underlying kinematics of the pose -i.e. FoldTree . By default, all backbone torsions are movable. ",
		attlist, subelements );
}

std::string SingleFragmentMoverCreator::keyname() const {
	return SingleFragmentMover::mover_name();
}

protocols::moves::MoverOP
SingleFragmentMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SingleFragmentMover );
}

void SingleFragmentMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SingleFragmentMover::provide_xml_schema( xsd );
}


}  // namespace nonlocal
}  // namespace protocols
