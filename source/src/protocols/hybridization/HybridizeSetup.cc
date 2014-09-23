// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Yifan Song

#include <protocols/hybridization/HybridizeSetup.hh>
#include <protocols/hybridization/HybridizeSetupMoverCreator.hh>

#include <protocols/hybridization/FoldTreeHybridize.hh>
#include <protocols/hybridization/CartesianHybridize.hh>
#include <protocols/hybridization/TemplateHistory.hh>
#include <protocols/hybridization/util.hh>
#include <protocols/hybridization/DomainAssembly.hh>
#include <protocols/hybridization/DDomainParse.hh>
#include <protocols/hybridization/TMalign.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/util.hh>

#include <protocols/loops/util.hh>
#include <protocols/loops/loops_main.hh>

#include <protocols/electron_density/SetupForDensityScoringMover.hh>

// dynamic fragpick
#include <protocols/moves/DsspMover.hh>
#include <core/fragment/picking_old/vall/util.hh>

#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.fwd.hh>

#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/util.hh>

#include <core/sequence/util.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// task operation
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

//#include <protocols/moves/DataMap.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// utility
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <numeric/random/WeightedSampler.hh>
#include <ObjexxFCL/format.hh>
#include <boost/foreach.hpp>

// evaluation
#include <core/scoring/rms_util.hh>
#include <protocols/comparative_modeling/coord_util.hh>

// option
#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh> // strand pairings
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>

#include <string>

#define foreach BOOST_FOREACH

static thread_local basic::Tracer TR( "protocols.hybridization.HybridizeSetup" );

namespace protocols {
namespace hybridization {

using namespace core;
using namespace kinematics;
using namespace sequence;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;
using namespace constraints;


HybridizeSetup::HybridizeSetup() :
template_weights_sum_(0)
{
	init();
}

void
HybridizeSetup::init() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	domain_assembly_ = false;
	add_hetatm_ = false;

	add_non_init_chunks_ = option[cm::hybridize::add_non_init_chunks]();
	frag_weight_aligned_ = option[cm::hybridize::frag_weight_aligned]();
	auto_frag_insertion_weight_ = option[cm::hybridize::auto_frag_insertion_weight]();
	frag_1mer_insertion_weight_ = option[cm::hybridize::frag_1mer_insertion_weight]();
	small_frag_insertion_weight_ = option[cm::hybridize::small_frag_insertion_weight]();
	big_frag_insertion_weight_ = option[cm::hybridize::big_frag_insertion_weight]();

	hetatm_self_cst_weight_ = 10.;
	hetatm_prot_cst_weight_ = 0.;
}

/////////////
// creator
std::string
HybridizeSetupMoverCreator::keyname() const {
	return HybridizeSetupMoverCreator::mover_name();
}

protocols::moves::MoverOP
HybridizeSetupMoverCreator::create_mover() const {
	return new HybridizeSetupMover;
}

std::string
HybridizeSetupMoverCreator::mover_name() {
	return "HybridizeSetupMover";
}

/////////////
// mover
HybridizeSetupMover::HybridizeSetupMover()
{
	hybridize_setup_ = new protocols::hybridization::HybridizeSetup();
}

// sets default options
void
HybridizeSetup::check_and_create_fragments( core::pose::Pose const & pose ) {
	core::pose::Pose pose_copy(pose);
	if (fragments_big_.size() > 0 && fragments_small_.size() > 0) return;

	if (!fragments_big_.size()) {
		core::fragment::FragSetOP frags = new core::fragment::ConstantLengthFragSet( 9 );

		// number of residues
		core::Size nres_tgt = pose.total_residue();
		core::conformation::symmetry::SymmetryInfoCOP symm_info;
		if ( core::pose::symmetry::is_symmetric(pose) ) {
			core::conformation::symmetry::SymmetricConformation & SymmConf (
				dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose_copy.conformation() ) );
			symm_info = SymmConf.Symmetry_Info();
			nres_tgt = symm_info->num_independent_residues();
		}
		if (pose.residue(nres_tgt).aa() == core::chemical::aa_vrt) nres_tgt--;
		while (!pose.residue(nres_tgt).is_protein()) nres_tgt--;

		// sequence
		std::string tgt_seq = pose.sequence();
		std::string tgt_ss(nres_tgt, '0');

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		if (option[ OptionKeys::in::file::psipred_ss2 ].user()) {
			utility::vector1< char > psipred = read_psipred_ss2_file( pose );
			for ( core::Size j=1; j<=nres_tgt; ++j ) {
				tgt_ss[j-1] = psipred[j];
			}
		} else {
			// templates vote on secstruct
			for (core::Size i=1; i<=templates_.size(); ++i) {
				for (core::Size j=1; j<=templates_[i]->total_residue(); ++j ) {
					if (!templates_[i]->residue(j).is_protein()) continue;
					core::Size tgt_pos = templates_[i]->pdb_info()->number(j);

					runtime_assert( tgt_pos<=nres_tgt );
					char tgt_ss_j = templates_[i]->secstruct(j);

					if (tgt_ss[tgt_pos-1] == '0') {
						tgt_ss[tgt_pos-1] = tgt_ss_j;
					} else if (tgt_ss[tgt_pos-1] != tgt_ss_j) {
						tgt_ss[tgt_pos-1] = 'D'; // templates disagree
					}
				}
			}
			for ( core::Size j=1; j<=nres_tgt; ++j ) {
				if (tgt_ss[j-1] == '0') tgt_ss[j-1] = 'D';
			}
		}

		// pick from vall based on template SS + target sequence
		for ( core::Size j=1; j<=nres_tgt-8; ++j ) {
			core::fragment::FrameOP frame = new core::fragment::Frame( j, 9 );
			frame->add_fragment(
				core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa( tgt_ss.substr( j-1, 9 ), tgt_seq.substr( j-1, 9 ), 25, true, core::fragment::IndependentBBTorsionSRFD() ) );
			frags->add( frame );
		}
		fragments_big_.push_back( frags );
	}
	if (!fragments_small_.size()) {
		core::fragment::FragSetOP frags = new core::fragment::ConstantLengthFragSet( 3 );

		// make them from big fragments
		core::fragment::chop_fragments( *fragments_big_[1], *frags );
		fragments_small_.push_back( frags );
	}
}

void HybridizeSetup::add_template(
	std::string template_fn,
	std::string cst_fn,
	std::string symm_file,
	core::Real weight,
	core::Real domain_assembly_weight,
	utility::vector1<core::Size> cst_reses)
{
	core::chemical::ResidueTypeSetCOP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	core::pose::PoseOP template_pose = new core::pose::Pose();
	core::import_pose::pose_from_pdb( *template_pose, *residue_set, template_fn );

	add_template( template_pose, cst_fn, symm_file, weight, domain_assembly_weight, cst_reses, template_fn );
}


void HybridizeSetup::add_template(
	core::pose::PoseOP template_pose,
	std::string cst_fn,
	std::string symm_file,
	core::Real weight,
	core::Real domain_assembly_weight,
	utility::vector1<core::Size> cst_reses,
	std::string name )
{
	// add secondary structure information to the template pose
	core::scoring::dssp::Dssp dssp_obj( *template_pose );
	dssp_obj.insert_ss_into_pose( *template_pose );

	// if pdbinfo doesn't exist (not read from pdb file)
	if (!template_pose->pdb_info()) {
		core::pose::PDBInfoOP new_pdb_info = new core::pose::PDBInfo( *template_pose );

		utility::vector1< int > pdb_numbering;
		utility::vector1< char > pdb_chains;
		for (Size ires=1; ires<=template_pose->total_residue(); ++ires) {
			pdb_numbering.push_back( ires );
			pdb_chains.push_back( 'A' );
		}

		new_pdb_info->set_numbering( pdb_numbering );
		new_pdb_info->set_chains( pdb_chains );
		template_pose->pdb_info( new_pdb_info );
		template_pose->pdb_info()->obsolete( false );
	}

	// find ss chunks in template
	protocols::loops::Loops contigs = protocols::loops::extract_continuous_chunks(*template_pose);
	protocols::loops::Loops chunks = protocols::loops::extract_secondary_structure_chunks(*template_pose, "HE", 3, 6, 3, 4);

	if (chunks.num_loop() == 0)
		chunks = contigs;

	TR.Debug << "Chunks from template\n" << chunks << std::endl;
	TR.Debug << "Contigs from template\n" << contigs << std::endl;

	template_fn_.push_back(name);
	templates_.push_back(template_pose);
	template_cst_fn_.push_back(cst_fn);
	symmdef_files_.push_back(symm_file);
	template_weights_.push_back(weight);
	domain_assembly_weights_.push_back(domain_assembly_weight);
	template_chunks_.push_back(chunks);
	template_contigs_.push_back(contigs);
	template_cst_reses_.push_back(cst_reses);
}


void HybridizeSetup::pick_starting_template()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if (starting_templates_.size() > 0) {
		numeric::random::WeightedSampler weighted_sampler;
		for (Size i=1; i<=starting_templates_.size(); ++i) {
			weighted_sampler.add_weight(template_weights_[starting_templates_[i]]);
		}
		Size k = weighted_sampler.random_sample(numeric::random::rg());
		initial_template_index_ = starting_templates_[k];
	}
	else {
		numeric::random::WeightedSampler weighted_sampler;
		weighted_sampler.weights(template_weights_);
		initial_template_index_ = weighted_sampler.random_sample(numeric::random::rg());
	}
}

void HybridizeSetup::realign_templates(core::pose::PoseCOP ref_pose)
{
	// realign each template to the starting template by domain
	// does not to domain realignment if in domain assembly mode
	// domain parsing
	DDomainParse ddom(pcut_,hcut_,length_);
	utility::vector1< utility::vector1< loops::Loops > > domains_all_templ;
	domains_all_templ.resize( templates_.size() );
	for (Size i_template=1; i_template<=templates_.size(); ++i_template) {
		domains_all_templ[i_template] = ddom.split( *templates_[i_template], templates_[i_template]->total_residue() );

		// protocols::loops::Loops my_chunks(template_chunks_[initial_template_index_]);
		// convert domain numbering to target pose numbering
		for (Size iloops=1; iloops<=domains_all_templ[i_template].size(); ++iloops) {
			for (Size iloop=1; iloop<=domains_all_templ[i_template][iloops].num_loop(); ++iloop) {
				Size seqpos_start_pose = templates_[i_template]->pdb_info()->number(domains_all_templ[i_template][iloops][iloop].start());
				domains_all_templ[i_template][iloops][iloop].set_start( seqpos_start_pose );

				Size seqpos_stop_pose = templates_[i_template]->pdb_info()->number(domains_all_templ[i_template][iloops][iloop].stop());
				domains_all_templ[i_template][iloops][iloop].set_stop( seqpos_stop_pose );
			}
		}

		TR << "Found " << domains_all_templ[i_template].size() << " domains for template " << template_fn_[i_template] << std::endl;
		for (Size i=1; i<=domains_all_templ[i_template].size(); ++i) {
			TR << "domain " << i << ": " << domains_all_templ[i_template][i] << std::endl;
		}
	}

	// combine domains that are not in the initial template
	domains_ = expand_domains_to_full_length(domains_all_templ, initial_template_index(), nres_tgt_);
	TR << "Final decision: " << domains_.size() << " domains" << std::endl;
	for (Size i=1; i<= domains_.size(); ++i) {
		TR << "domain " << i << ": " << domains_[i] << std::endl;
	}

	// local align
	align_by_domain(templates_, domains_, ref_pose);

	// update chunk, contig informations
	for (Size i_template=1; i_template<=templates_.size(); ++i_template) {
		// default minimum length is 3 and CA distance is 4
		template_contigs_[i_template] = protocols::loops::extract_continuous_chunks(*templates_[i_template]); // for chunk insertions
		template_chunks_[i_template] = protocols::loops::extract_secondary_structure_chunks(*templates_[i_template], "HE", 3, 6, 3, 4); // for fold tree setup
		if (template_chunks_[i_template].num_loop() == 0)
			template_chunks_[i_template] = template_contigs_[i_template];
	}
}

void HybridizeSetupMover::apply( core::pose::Pose & pose )
{
	using namespace protocols::moves;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose::datacache;
	using namespace core::io::silent;
	using namespace ObjexxFCL::format;

	// starting structure pool
	std::vector < SilentStructOP > post_centroid_structs;

	// Three opts for template preprocessing:
	//  -> adding hetero residues
	//  -> domain assembly
	//  -> local realignment
	// Currently they are mutually exclusive
	utility::vector1< std::pair< core::Size,core::Size > > hetatms;
	if (hybridize_setup_->add_hetatm()) {
		for ( Size ires=1; ires <= hybridize_setup_->template_poses()[hybridize_setup_->initial_template_index()]->total_residue(); ++ires ) {
			if (hybridize_setup_->template_poses()[hybridize_setup_->initial_template_index()]->pdb_info()->number(ires) > (int)hybridize_setup_->nres_tgt_asu()) {
				TR.Debug << "Insert hetero residue: " << hybridize_setup_->template_poses()[hybridize_setup_->initial_template_index()]->residue(ires).name3() << std::endl;
				if ( hybridize_setup_->template_poses()[hybridize_setup_->initial_template_index()]->residue(ires).is_polymer()
					&& !hybridize_setup_->template_poses()[hybridize_setup_->initial_template_index()]->residue(ires).is_lower_terminus() ) {
					pose.append_residue_by_bond(hybridize_setup_->template_poses()[hybridize_setup_->initial_template_index()]->residue(ires));
				} else {
					pose.append_residue_by_jump(hybridize_setup_->template_poses()[hybridize_setup_->initial_template_index()]->residue(ires), 1);
				}
				hetatms.push_back( std::make_pair( ires, pose.total_residue() ) );
			}
		}
	}


	// symmetrize
	std::string symmdef_file = hybridize_setup_->symmdef_files()[hybridize_setup_->initial_template_index()];
	if (!symmdef_file.empty() && symmdef_file != "NULL") {
		protocols::simple_moves::symmetry::SetupForSymmetryMover makeSymm( symmdef_file );
		makeSymm.apply(pose);

        TR << "After applying symmetry" << std::endl;

		//fpd   to get the right rotamer set we need to do this
		basic::options::option[basic::options::OptionKeys::symmetry::symmetry_definition].value( "dummy" );

		// xyz copy hetatms -- this makes helical symmetries a little easier to setup
		if (hybridize_setup_->add_hetatm()) {
			for ( Size ihet=1; ihet <= hetatms.size(); ++ihet ) {
				core::conformation::Residue const &res_in = hybridize_setup_->template_poses()[hybridize_setup_->initial_template_index()]->residue(hetatms[ihet].first);

				for (Size iatm=1; iatm<=res_in.natoms(); ++iatm) {
					core::id::AtomID tgt(iatm,hetatms[ihet].second);
					pose.set_xyz( tgt, res_in.xyz( iatm ) );
				}
			}
		}
	}

    // initialize template history
    // >> keep this after symmetry
    TemplateHistoryOP history = new TemplateHistory(pose);
    history->setall( hybridize_setup_->initial_template_index() );
    pose.data().set( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY, history );

	// set pose for density scoring if a map was input
	// keep this after symmetry
	if ( option[ OptionKeys::edensity::mapfile ].user() ) {
		MoverOP dens( new protocols::electron_density::SetupForDensityScoringMover );
		dens->apply( pose );
	}
}

utility::vector1 <Loops>
HybridizeSetup::expand_domains_to_full_length(utility::vector1 < utility::vector1 < Loops > > all_domains, Size ref_domains_index, Size n_residues)
{
	utility::vector1 < Loops > domains = all_domains[ref_domains_index];
	if (all_domains[ref_domains_index].size() == 0) return domains;

	utility::vector1 < bool > residue_mask(n_residues, false);
	residue_mask.resize(n_residues);

	//mask residues from all domains not to be cut
	for (Size i_pose=1; i_pose <= all_domains.size(); ++i_pose) {
		for (Size idomain=1; idomain <= all_domains[i_pose].size(); ++idomain) {
			for (core::Size iloop=1; iloop<=all_domains[i_pose][idomain].num_loop(); ++iloop) {
				for (core::Size ires=all_domains[i_pose][idomain][iloop].start()+1; ires<=all_domains[i_pose][idomain][iloop].stop(); ++ires) {
					residue_mask[ires] = true;
				}
			}
		}
	}

	// assuming loops in domains are sequential
	for (Size idomain=1; idomain <= all_domains[ref_domains_index].size(); ++idomain) {
		for (core::Size iloop=1; iloop<=all_domains[ref_domains_index][idomain].num_loop(); ++iloop) {
			if (idomain == 1 && iloop == 1) {
				domains[idomain][iloop].set_start(1);
			}
			if (idomain == all_domains[ref_domains_index].size() && iloop == all_domains[ref_domains_index][idomain].num_loop()) {
				domains[idomain][iloop].set_stop(n_residues);
			}

			// the next loop
			Size jdomain = idomain;
			Size jloop = iloop+1;
			if (jloop > all_domains[ref_domains_index][idomain].num_loop()) {
				++jdomain;
				jloop = 1;
				if (jdomain > all_domains[ref_domains_index].size()) continue;

				Size gap_start = all_domains[ref_domains_index][idomain][iloop].stop()+1;
				Size gap_stop  = all_domains[ref_domains_index][jdomain][jloop].start();
				utility::vector1<Size> cut_options;
				for (Size ires=gap_start; ires<=gap_stop; ++ires) {
					if (residue_mask[ires]) continue;
					cut_options.push_back(ires);
				}
				if (cut_options.size() == 0) {
					for (Size ires=gap_start; ires<=gap_stop; ++ires) {
						cut_options.push_back(ires);
					}
				}
				Size cut = cut_options[numeric::random::rg().random_range(1,cut_options.size())];

				domains[idomain][iloop].set_stop(cut-1);
				domains[jdomain][jloop].set_start(cut);
			}
		}
	}
	return domains;
}

void
HybridizeSetup::align_by_domain(utility::vector1<core::pose::PoseOP> & poses, utility::vector1 < Loops > domains, core::pose::PoseCOP & ref_pose)
{
	for (Size i_pose=1; i_pose <= poses.size(); ++i_pose) {
		if (poses[i_pose] == ref_pose) continue; // same pose
		align_by_domain(*poses[i_pose], *ref_pose, domains);

		//std::string out_fn = template_fn_[i_pose] + "_realigned.pdb";
		//poses[i_pose]->dump_pdb(out_fn);
	}
}


void
HybridizeSetup::align_by_domain(core::pose::Pose & pose, core::pose::Pose const & ref_pose, utility::vector1 <Loops> domains)
{
	for (Size i_domain = 1; i_domain <= domains.size() ; ++i_domain) {
		// To avoid doing TMalign from a small fragment to a large fragment, adjust residue_list to have a large overlap in residue numbers between the two lists
		Loop pose_domain(domains[i_domain][1]);
		for (core::Size ires=1; ires<=pose.total_residue(); ++ires) {
			if (!pose.residue_type(ires).is_protein()) continue;
			int pose_res = (pose.pdb_info()) ? pose.pdb_info()->number(ires) : ires;
			for (core::Size iloop=1; iloop<=domains[i_domain].num_loop(); ++iloop) {
				if ( pose_res < (int)domains[i_domain][iloop].start() || pose_res > (int)domains[i_domain][iloop].stop() ) continue;
				if (pose_res < (int)pose_domain.start()) {
					pose_domain.set_start(pose_res);
				}
				if (pose_res > (int)pose_domain.stop()) {
					pose_domain.set_stop(pose_res);
				}
			}
		}

		Loop ref_domain(domains[i_domain][1]);
		for (core::Size jres=1; jres<=ref_pose.total_residue(); ++jres) {
			if (!ref_pose.residue_type(jres).is_protein()) continue;
			int ref_pose_res = (ref_pose.pdb_info()) ? ref_pose.pdb_info()->number(jres) : jres;
			for (core::Size iloop=1; iloop<=domains[i_domain].num_loop(); ++iloop) {
				if ( ref_pose_res < (int)domains[i_domain][iloop].start() || ref_pose_res > (int)domains[i_domain][iloop].stop() ) continue;
				if (ref_pose_res < (int)ref_domain.start()) {
					ref_domain.set_start(ref_pose_res);
				}
				if (ref_pose_res > (int)ref_domain.stop()) {
					ref_domain.set_stop(ref_pose_res);
				}
			}
		}

		Loop overlapped_domain;
		overlapped_domain.set_start(pose_domain.start() > ref_domain.start() ? pose_domain.start() : ref_domain.start());
		overlapped_domain.set_stop(pose_domain.stop() < ref_domain.stop() ? pose_domain.stop() : ref_domain.stop());
		// done looking for overlapped domain


		core::id::AtomID_Map< core::id::AtomID > atom_map;
		core::pose::initialize_atomid_map( atom_map, pose, core::id::BOGUS_ATOM_ID );
		std::list <Size> residue_list;
		core::Size n_mapped_residues=0;
		for (core::Size ires=1; ires<=pose.total_residue(); ++ires) {
			if (!pose.residue_type(ires).is_protein()) continue;
			int pose_res = (pose.pdb_info()) ? pose.pdb_info()->number(ires) : ires;
			if ( pose_res < (int)overlapped_domain.start() || pose_res > (int)overlapped_domain.stop() ) continue;
			for (core::Size iloop=1; iloop<=domains[i_domain].num_loop(); ++iloop) {
				if ( pose_res < (int)domains[i_domain][iloop].start() || pose_res > (int)domains[i_domain][iloop].stop() ) continue;
				residue_list.push_back(ires);
			}
		}
		std::list <Size> ref_residue_list;
		for (core::Size jres=1; jres<=ref_pose.total_residue(); ++jres) {
			if (!ref_pose.residue_type(jres).is_protein()) continue;
			int ref_pose_res = (ref_pose.pdb_info()) ? ref_pose.pdb_info()->number(jres) : jres;
			if ( ref_pose_res < (int)overlapped_domain.start() || ref_pose_res > (int)overlapped_domain.stop() ) continue;
			for (core::Size iloop=1; iloop<=domains[i_domain].num_loop(); ++iloop) {
				if ( ref_pose_res < (int)domains[i_domain][iloop].start() || ref_pose_res > (int)domains[i_domain][iloop].stop() ) continue;
				ref_residue_list.push_back(jres);
			}
		}

		TMalign tm_align;
		std::string seq_pose, seq_ref, aligned;
		int reval = tm_align.apply(pose, ref_pose, residue_list, ref_residue_list);
		if (reval == 0) {
			tm_align.alignment2AtomMap(pose, ref_pose, residue_list, ref_residue_list, n_mapped_residues, atom_map);
			tm_align.alignment2strings(seq_pose, seq_ref, aligned);

			using namespace ObjexxFCL::format;
			Size norm_length = residue_list.size() < ref_residue_list.size() ? residue_list.size():ref_residue_list.size();
			TR << "Align domain with TMscore of " << F(8,3,tm_align.TMscore(norm_length)) << std::endl;
			TR << seq_pose << std::endl;
			TR << aligned << std::endl;
			TR << seq_ref << std::endl;

			if (n_mapped_residues >= 6) {
				utility::vector1< core::Real > aln_cutoffs;
				aln_cutoffs.push_back(6);
				aln_cutoffs.push_back(4);
				aln_cutoffs.push_back(3);
				aln_cutoffs.push_back(2);
				aln_cutoffs.push_back(1.5);
				aln_cutoffs.push_back(1);
				core::Real min_coverage = 0.2;
				partial_align(pose, ref_pose, atom_map, residue_list, true, aln_cutoffs, min_coverage); // iterate_convergence = true
			}
		}
		//else {
		//	TR << "This domain cannot be aligned: " << n_mapped_residues<< std::endl;
		//	TR << domains[i_domain];
		//}
	}
}

void
HybridizeSetup::run_domain_assembly(){
	DomainAssembly domain_assembly(templates_, domain_assembly_weights_);
	domain_assembly.run();
}

HybridizeSetupOP HybridizeSetup::clone() const { return new HybridizeSetup( *this ); }


protocols::moves::MoverOP HybridizeSetupMover::clone() const { return new HybridizeSetupMover( *this ); }
protocols::moves::MoverOP HybridizeSetupMover::fresh_instance() const { return new HybridizeSetupMover; }

std::string
HybridizeSetupMover::get_name() const {
	return "HybridizeSetupMover";
}

void
HybridizeSetupMover::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		filters::Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const & pose )
{
	std::string setup_data_name = "hybridize_setup";
	if( tag->hasOption( "name" ) ) {
		setup_data_name = tag->getOption< std::string >( "name" );
	}

	// force starting template
	if( tag->hasOption( "starting_template" ) ) {
		utility::vector1<std::string> buff = utility::string_split( tag->getOption<std::string>( "starting_template" ), ',' );
		foreach(std::string field, buff){
			Size const value = std::atoi( field.c_str() );
			hybridize_setup_->add_starting_templates(value);
		}
	}

	// tons of ab initio options
	if( tag->hasOption( "add_hetatm" ) )
		hybridize_setup_->set_add_hetatm( tag->getOption< bool >( "add_hetatm" ) );
	if( tag->hasOption( "hetatm_cst_weight" ) )
		hybridize_setup_->set_hetatm_self_cst_weight( tag->getOption< core::Real >( "hetatm_cst_weight" ) );
	if( tag->hasOption( "hetatm_to_protein_cst_weight" ) )
		hybridize_setup_->set_hetatm_prot_cst_weight( tag->getOption< core::Real >( "hetatm_to_protein_cst_weight" ) );

	if( tag->hasOption( "domain_assembly" ) )
		hybridize_setup_->set_domain_assembly( tag->getOption< bool >( "domain_assembly" ) );

	if( tag->hasOption( "align_domains_to_template" ) )
		hybridize_setup_->set_align_domains_to_template( tag->getOption< bool >( "align_domains_to_template" ) );
	if( tag->hasOption( "align_domains_to_pose" ) )
		hybridize_setup_->set_align_domains_to_pose( tag->getOption< bool >( "align_domains_to_pose" ) );

	//if( tag->hasOption( "disulf_file" ) )
	//	disulf_file_ = tag->getOption< std::string >( "disulf_file" );

	// ddomain options
	hybridize_setup_->set_domain_hcut( tag->getOption< core::Real >( "domain_hcut" , 0.18) );
	hybridize_setup_->set_domain_pcut( tag->getOption< core::Real >( "domain_pcut" , 0.81) );
	hybridize_setup_->set_domain_length( tag->getOption< core::Size >( "domain_length" , 38) );

	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
	for (tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it) {
		// fragments
		if ( (*tag_it)->getName() == "Fragments" ) {
			using namespace core::fragment;
			if ( (*tag_it)->hasOption( "3mers" ) ) {
				core::fragment::FragSetOP frags = core::fragment::FragmentIO().read_data( (*tag_it)->getOption<std::string>( "3mers" )  );
				hybridize_setup_->add_fragments_small(frags);
			} else if ( (*tag_it)->hasOption( "small" ) ) {
				utility::vector1<std::string> frag_files = utility::string_split((*tag_it)->getOption<std::string>( "small" ), ',');
				for (core::Size i=1; i<= frag_files.size(); ++i )
					hybridize_setup_->add_fragments_small(core::fragment::FragmentIO().read_data( frag_files[i] ));
			}

			if ( (*tag_it)->hasOption( "9mers" ) ) {
				core::fragment::FragSetOP frags = core::fragment::FragmentIO().read_data( (*tag_it)->getOption<std::string>( "9mers" ) );
				hybridize_setup_->add_fragments_big(frags);
			} else if ( (*tag_it)->hasOption( "big" ) ) {
				utility::vector1<std::string> frag_files = utility::string_split((*tag_it)->getOption<std::string>( "big" ), ',');
				for (core::Size i=1; i<= frag_files.size(); ++i )
					hybridize_setup_->add_fragments_big(core::fragment::FragmentIO().read_data( frag_files[i] ));
			}
		}

		// templates
		if ( (*tag_it)->getName() == "Template" ) {
			std::string template_fn = (*tag_it)->getOption<std::string>( "pdb" );
			std::string cst_fn = (*tag_it)->getOption<std::string>( "cst_file", "AUTO" );
			core::Real weight = (*tag_it)->getOption<core::Real>( "weight", 1 );
			core::Real domain_assembly_weight = (*tag_it)->getOption<core::Real>( "domain_assembly_weight", 0. );
			std::string symm_file = (*tag_it)->getOption<std::string>( "symmdef", "" );
			utility::vector1<core::Size> cst_reses;
			hybridize_setup_->add_template(template_fn, cst_fn, symm_file, weight, domain_assembly_weight, cst_reses);
		}

	}

	// Additional setup
	// number of residues in asu without VRTs
	core::Size nres_tgt = pose.total_residue();
	core::Size nres_protein_tgt = pose.total_residue();
	core::pose::Pose pose_copy(pose);
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose_copy) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
				dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose_copy.conformation()) );

		symm_info = SymmConf.Symmetry_Info();
		nres_tgt = symm_info->num_independent_residues();
		nres_protein_tgt = symm_info->num_independent_residues();
	}
	if (pose.residue(nres_tgt).aa() == core::chemical::aa_vrt) nres_tgt--;
	while (!pose.residue(nres_protein_tgt).is_protein()) nres_protein_tgt--;
	hybridize_setup_->set_nres_tgt(nres_tgt);
	hybridize_setup_->set_nres_protein(nres_protein_tgt);

	// make fragments if we don't have them at this point
	hybridize_setup_->check_and_create_fragments( pose );

	// pick starting template
	hybridize_setup_->pick_starting_template();

	if (hybridize_setup_->domain_assembly()) {
		hybridize_setup_->run_domain_assembly();
	}

	using namespace ObjexxFCL::format;
	TR << "Using initial template: " << I(4,hybridize_setup_->initial_template_index()) << " " << hybridize_setup_->template_file_names()[hybridize_setup_->initial_template_index()] << std::endl;

	data.add( "HybridizeSetup", setup_data_name, hybridize_setup_);

}

} // hybridization
} // protocols
