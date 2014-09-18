// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SewingDesigner.cc
///
/// @brief
/// @author Tim Jacobs

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SewingDesigner.cc
///
/// @brief
/// @author Tim Jacobs

// Unit Headers
#include <devel/sewing/SewingDesigner.hh>
#include <devel/sewing/SewingDesignerCreator.hh>
#include <devel/sewing/util.hh>

//Basic
#include <basic/MetricValue.hh>
#include <basic/datacache/DiagnosticData.hh>
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>

//Core
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/util.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/rotamer_set/AddResiduesRotamerSetOperation.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//Protocols
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/Loop.hh>
#include <devel/loop_creation/LoopCreationMover.hh>
#include <devel/loop_creation/LoopCloser.hh>
#include <devel/loop_creation/LoopInserter.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>

//Utility
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/sort_predicates.hh>
#include <utility/LexicographicalIterator.hh>
#include <utility/tag/Tag.hh>
#include <utility/file/FileName.hh>

//Basic
#include <basic/Tracer.hh>

//Devel
#include <devel/sewing/NativeResidueReader.hh>

namespace devel {
namespace sewing {

using namespace core;
using namespace io;
using namespace pdb;
using namespace chemical;
using namespace conformation;
using namespace utility;

using core::pose::PoseOP;
using core::pose::PoseAP;
using core::pose::Pose;
using utility::vector1;
using core::Size;
using std::string;

static thread_local basic::Tracer TR( "protocols.loophash.SewingDesigner" );

/****CREATOR FUNCTIONS*****/
protocols::moves::MoverOP
SewingDesignerCreator::create_mover() const
{
	return new SewingDesigner;
}

std::string
SewingDesignerCreator::keyname() const
{
	return SewingDesignerCreator::mover_name();
}

std::string
SewingDesignerCreator::mover_name()
{
	return "SewingDesigner";
}
/****END CREATOR FUNCTIONS*****/

SewingDesigner::SewingDesigner():
Mover("SewingDesigner")
{
	init();
}

SewingDesigner::SewingDesigner(
	NativeRotamersMap native_helix_residues
):
Mover("SewingDesigner"),
native_helix_residues_(native_helix_residues)
{
	init();
}

SewingDesigner::SewingDesigner(
	devel::loop_creation::LoopCreationMoverOP loop_creation_mover,
	NativeRotamersMap native_helix_residues,
	utility::vector1<core::Size> loop_anchors
):
Mover("SewingDesigner"),
loop_creation_mover_(loop_creation_mover),
native_helix_residues_(native_helix_residues),
loop_anchors_(loop_anchors)
{
	init();
}

protocols::moves::MoverOP
SewingDesigner::clone() const {
	return( protocols::moves::MoverOP( new SewingDesigner( *this ) ) );
}
protocols::moves::MoverOP
SewingDesigner::fresh_instance() const {
	return protocols::moves::MoverOP( new SewingDesigner );
}

string
SewingDesigner::get_name() const {
	return "SewingDesigner";
}

void
SewingDesigner::init()
{
	using namespace basic::resource_manager;
	using namespace basic::options;

	if ( ResourceManager::get_instance()->has_option( OptionKeys::sewing::nat_ro_file ) ||  option[ OptionKeys::sewing::nat_ro_file ].user() ) {
		utility::file::FileName native_res_filename = ResourceManager::get_instance()->get_option( OptionKeys::sewing::nat_ro_file );
		NativeResidueReader native_res_reader;
		native_helix_residues_  =
				native_res_reader.generateResiduesFromFile(native_res_filename.name());
	}
	else{
		utility_exit_with_message("You must provide a nat ro file to the sewing designer mover.");
	}
}

void
SewingDesigner::apply( pose::Pose & pose )
{
	//Setup scorefunction and movers needed for design
	scorefxn_ = scoring::get_score_function();
	scorefxn_->set_weight(core::scoring::res_type_constraint, 1.0);

	repack_mover_ = new protocols::simple_moves::PackRotamersMover;
	repack_mover_->score_function(scorefxn_);

	//Look for breaks in the pose and keep residues as anchors for loops
	if(loop_anchors_.size()==0){
		detect_loop_anchors(pose);
	}

	//optimally rearrange helix segments to minimize loop distance
	rearrange_pose(pose, native_helix_residues_, loop_anchors_);

	//Do an initial design of the entire protein
	TR << "Statistic before design: " << std::endl;
	std::map<core::Size, char> empty;
	record_statistics(pose, native_helix_residues_, empty);
	design(pose, native_helix_residues_);
	pose.dump_pdb("initial_design.pdb");

	TR << "Statistic after initial design: " << std::endl;
	record_statistics(pose, native_helix_residues_, empty);

	//Try to design a loop after each loop_anchor
	std::map<Size,char> native_loop_residues;
	for(core::Size loop_index=1; loop_index<=loop_anchors_.size(); ++loop_index)
	{
		core::Size cur_anchor = loop_anchors_[loop_index];
		TR << "Trying to insert loop after residue: " << cur_anchor << std::endl;

		//update the loop_inserter's loop_anchor and apply
		loop_creation_mover_->loop_inserter()->loop_anchor(cur_anchor);
		loop_creation_mover_->apply(pose);

		protocols::loops::Loop new_loop=loop_creation_mover_->get_last_created_loop();

		std::set<Size> loop_residues;
		for(core::Size i=new_loop.start(); i<=new_loop.stop(); ++i)
		{
			loop_residues.insert(i);
		}

		std::string const nb_calc("loop_neighboorhood_calculator");
		pose::metrics::CalculatorFactory::Instance().register_calculator( nb_calc,
			new protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( loop_residues ) );

		basic::MetricValue< std::set< Size > > neighbor_mv;
		pose.metric( nb_calc, "neighbors", neighbor_mv);
		std::set<Size> const neighbor_set ( neighbor_mv.value() );

		//Update native rotamers map for the given insert
		update_native_helix_residues(native_helix_residues_, new_loop);
		update_anchors(loop_anchors_, new_loop, loop_index);

		//record all the loop residues before design so we can generate recovery metrics
		native_loop_residues.clear();
		for(core::Size i=1; i<=pose.total_residue(); ++i)
		{
			//If this position isn't a helix, it's a loop
			if(native_helix_residues_.find(i)==native_helix_residues_.end())
			{
				native_loop_residues.insert(std::make_pair(i, pose.residue(i).type().name1()));
			}
		}

		asym_size_+=new_loop.size();
		core::Size dup_anchor=cur_anchor+asym_size_;
		TR << "asym size is: " << asym_size_ << std::endl;
		TR << "size of new loop: " << new_loop.size() << std::endl;
		while(dup_anchor<pose.total_residue())//strictly less than so we don't try to build a loop that connects to nothing
		{
			TR << "Duplication loop build on anchor " << cur_anchor << " to new anchor: " << dup_anchor << std::endl;
			pose.dump_pdb("pre_dup_anchor_"+utility::to_string(cur_anchor)+".pdb");
			loop_creation_mover_->copy_last_loop_to_new_anchor(pose, dup_anchor);
			dup_anchor+=asym_size_;
		}

		//Design/repack & score
		design_neighborhood(pose, native_helix_residues_, neighbor_set);

		//remove calculator so we don't get name conflicts
		pose::metrics::CalculatorFactory::Instance().remove_calculator( nb_calc );
	}
	//final design of everything
//	design(pose, native_helix_residues_);

	//record final statistics
	record_statistics(pose, native_helix_residues_, native_loop_residues);
	if( !pose.data().has( core::pose::datacache::CacheableDataType::SCORE_MAP ) )
	{
		utility_exit_with_message("No native info was added to the pose.");
	}
	std::map<std::string, core::Real> const & native_retention_map =
		( static_cast< basic::datacache::DiagnosticData const &>( pose.data().get( core::pose::datacache::CacheableDataType::SCORE_MAP ))).data();

	protocols::jd2::JobOP const job_me = protocols::jd2::JobDistributor::get_instance()->current_job();
	job_me->add_string_real_pair("native_helix_percent", native_retention_map.find("native_helix_percent")->second);
	job_me->add_string_real_pair("native_helix_percent_buried", native_retention_map.find("native_helix_percent_buried")->second);
	job_me->add_string_real_pair("native_helix_percent_exposed", native_retention_map.find("native_helix_percent_exposed")->second);
	job_me->add_string_real_pair("native_loop_percent", native_retention_map.find("native_loop_percent")->second);
	job_me->add_string_real_pair("native_loop_pro_percent", native_retention_map.find("native_loop_pro_percent")->second);
	job_me->add_string_real_pair("native_loop_gly_percent", native_retention_map.find("native_loop_gly_percent")->second);
}

///@brief Keep track of the percentages of native residues that
///were retained throughout the design process. Add these statistics
///to the job data.
void
SewingDesigner::record_statistics(
	core::pose::Pose & pose,
	NativeRotamersMap const & native_helix_residues,
	std::map<core::Size, char> const & native_loop_residues
){
	//Calculate number of "native" helix and loop residues are preserved
	std::string sequence = pose.sequence();
	Size total_helix_exposed=0;
	Size total_helix_buried=0;

	Size recovered_helix_total=0;
	Size recovered_helix_exposed=0;
	Size recovered_helix_buried=0;

	Size recovered_loop=0;

	Size total_loop_pro=0;
	Size total_loop_gly=0;

	Size recovered_loop_pro=0;
	Size recovered_loop_gly=0;

	std::set<core::Size> retained_helix_residues;
	std::set<core::Size> retained_loop_residues;

	//Score to make sure we have a good tenA neighbor graph
	scorefxn_->score(pose);
	for ( Size cur_res = 1 ; cur_res <= pose.total_residue(); cur_res++)
	{
		NativeRotamersMap::const_iterator it = native_helix_residues.find(cur_res);
		if(it != native_helix_residues.end())
		{
			bool buried=true;
			if(pose.energies().tenA_neighbor_graph().get_node(cur_res)->num_neighbors_counting_self() <= 16)
			{
				buried=false;
				total_helix_exposed++;
			}
			else
			{
				total_helix_buried++;
			}
			std::set<char> res_type_names;//prevents double counting
			utility::vector1<core::conformation::ResidueOP> const & position_residues =
			it->second;
			for(Size pos_index=1; pos_index<=position_residues.size(); ++pos_index)
			{
				char const & aa = position_residues[pos_index]->type().name1();
				if(sequence[cur_res-1]==aa)
				{
					if(res_type_names.find(aa) == res_type_names.end())
					{
						recovered_helix_total++;
						retained_helix_residues.insert(cur_res);
						if(buried)
						{
							recovered_helix_buried++;
						}
						else
						{
							recovered_helix_exposed++;
						}
						res_type_names.insert(aa);
					}
				}
			}
		}
		else
		{
			bool pro=false;
			bool gly=false;
			if(native_loop_residues.find(cur_res)->second == 'P')
			{
				pro=true;
				total_loop_pro++;
			}
			else if(native_loop_residues.find(cur_res)->second == 'G')
			{
				gly=true;
				total_loop_gly++;
			}

			if(sequence[cur_res-1]==native_loop_residues.find(cur_res)->second)
			{
				recovered_loop++;
				retained_loop_residues.insert(cur_res);
				if(pro) recovered_loop_pro++;
				if(gly) recovered_loop_gly++;
			}
		}
	}
	std::map<std::string, core::Real> native_retention_map;

	core::Real native_helix_percent = (core::Real)recovered_helix_total/native_helix_residues.size();
	native_retention_map.insert(std::make_pair("native_helix_percent", native_helix_percent));

	core::Real native_helix_percent_buried = (core::Real)recovered_helix_buried/total_helix_buried;
	native_retention_map.insert(std::make_pair("native_helix_percent_buried", native_helix_percent_buried));

	core::Real native_helix_percent_exposed = (core::Real)recovered_helix_exposed/total_helix_exposed;
	native_retention_map.insert(std::make_pair("native_helix_percent_exposed", native_helix_percent_exposed));

	core::Real native_loop_percent = 0;
	if(native_loop_residues.size() != 0)
	{
		native_loop_percent = (core::Real)recovered_loop/native_loop_residues.size();
		native_retention_map.insert(std::make_pair("native_loop_percent", native_loop_percent));
	}

	core::Real pro_loop_percent = (core::Real)recovered_loop_pro/total_loop_pro;
	native_retention_map.insert(std::make_pair("native_loop_pro_percent", pro_loop_percent));

	core::Real gly_loop_percent = (core::Real)recovered_loop_gly/total_loop_gly;
	native_retention_map.insert(std::make_pair("native_loop_gly_percent", gly_loop_percent));

	std::string pymol_helix_select = "select resi ";
	for(std::set<Size>::const_iterator it=retained_helix_residues.begin(); it!= retained_helix_residues.end(); ++it)
	{
		pymol_helix_select += "+"+utility::to_string(*it);
	}
	pymol_helix_select += "/helix_residues";

	std::string pymol_loop_select = "select resi ";
	for(std::set<Size>::const_iterator it=retained_loop_residues.begin(); it!= retained_loop_residues.end(); ++it)
	{
		pymol_loop_select += "+"+utility::to_string(*it);
	}
	pymol_loop_select += "/loop_residues";


	protocols::jd2::JobOP const job_me ( protocols::jd2::JobDistributor::get_instance()->current_job() );
	std::string const job_name ( protocols::jd2::JobDistributor::get_instance()->job_outputter()->output_name(job_me) );

//	dump_native_residue_file(native_helix_residues, job_name+".rots");

	TR << "Helix residue select: " << pymol_helix_select << std::endl;
	TR << "Loop residue select: " << pymol_loop_select << std::endl;

	pose.data().set(core::pose::datacache::CacheableDataType::SCORE_MAP,
					new basic::datacache::DiagnosticData(native_retention_map));
	TR << "Percent of native helix residues (buried): " << native_helix_percent_buried << std::endl;
	TR << "Percent of native helix residues (exposed): " << native_helix_percent_exposed << std::endl;
	TR << "Percent of native loop residues: " << native_loop_percent << std::endl;
}

void
SewingDesigner::design(
	core::pose::Pose & pose,
	NativeRotamersMap const & native_helix_residues
){
	std::set<core::Size> all_residues;
	for(core::Size i=1; i<=pose.total_residue(); ++i){
		all_residues.insert(i);
	}
	design_neighborhood(pose, native_helix_residues, all_residues);
}

///@brief Setup the packer with the appropriate
///constraints and design the residues contained
///within the neighbor_residues list
void
SewingDesigner::design_neighborhood(
	core::pose::Pose & pose,
	NativeRotamersMap const & native_helix_residues,
	std::set<core::Size> const & neighbor_residues
){
	//Update loops for each of the new structures and continue closing
	pack::task::TaskFactoryOP task_factory = new pack::task::TaskFactory;
	//	task_factory->push_back( new pack::task::operation::InitializeFromCommandline );

	//Add native rotamers for each position in the updated native rotamers map
	for(NativeRotamersMap::const_iterator map_it = native_helix_residues.begin();
		map_it != native_helix_residues.end(); ++map_it)
	{
		//Create rotamer set operation from a list of residues
		pack::rotamer_set::AddResiduesRotamerSetOperation const & nat_ro_set(map_it->second);

//		TR.Debug << "Adding the following native rotamers for position: " << map_it->first << std::endl;
		for(Size temp=1; temp<=map_it->second.size(); ++temp)
		{
			conformation::ResidueOP cur_res = map_it->second[temp];
//			TR.Debug << "Residue Type: " << cur_res->type().name() << std::endl;
		}

		//Add AppendResidueRotamerSet task operation to the task factory. This task operation
		//adds the rotamer set to the residue-level task for the given residue
		task_factory->push_back(
			new core::pack::task::operation::AppendResidueRotamerSet(map_it->first, nat_ro_set.clone()) );
	}


	//Add constraints for each of the native residues at each position that is a neighbor
	//and just repack at non-neighbor residues
	core::Real native_bonus=0.5;
	pack::task::operation::RestrictResidueToRepackingOP repack_res =
	new pack::task::operation::RestrictResidueToRepacking();
	for(Size resnum=1; resnum<=pose.total_residue(); ++resnum)
	{
		//if this residue is a a nieghbor of the loop residues than favor all natives and allow design
		if(neighbor_residues.find(resnum) != neighbor_residues.end())
		{
			//helical position, so favor all "natives"
			if(native_helix_residues.find(resnum) != native_helix_residues.end())
			{
				utility::vector1<core::conformation::ResidueOP> const & position_residues =
					native_helix_residues.find(resnum)->second;
				for(Size pos_index=1; pos_index<=position_residues.size(); ++pos_index)
				{
					core::scoring::constraints::ResidueTypeConstraintOP nat_res_constraint =
						new core::scoring::constraints::ResidueTypeConstraint(pose, resnum,
							position_residues[pos_index]->name3(), native_bonus);

					pose.add_constraint(nat_res_constraint);

//					TR.Debug << "Favored residue for helix position " << resnum << ": "
//						<< position_residues[pos_index]->name() << std::endl;
				}
			}
			//newly built loop position, favor single native residue
			else
			{
				core::scoring::constraints::ResidueTypeConstraintOP nat_res_constraint =
					new core::scoring::constraints::ResidueTypeConstraint(pose, resnum, native_bonus);
				pose.add_constraint(nat_res_constraint);
				TR << "Favored residue for loop position " << resnum << ": "
					<< pose.residue(resnum).type().name() << std::endl;
			}
		}
		//not a neighbor or a loop residue, repack only
		else
		{
			repack_res->include_residue(resnum);
		}
	}

	task_factory->push_back(repack_res);
	repack_mover_->score_function(scorefxn_);
	repack_mover_->task_factory(task_factory);
	repack_mover_->apply(pose);
}

///@brief update the native rotamers map to reflect
///the position changes after an insertion
void
SewingDesigner::update_native_helix_residues(
	NativeRotamersMap & native_helix_residues,
	protocols::loops::Loop const & new_loop
){
	core::Size start_of_residue_change=new_loop.stop();
	core::Size num_residues_inserted=new_loop.size();

	NativeRotamersMap new_native_helix_residues;
	for(NativeRotamersMap::const_iterator map_it = native_helix_residues.begin();
		map_it != native_helix_residues.end(); ++map_it)
	{
		Size resnum=map_it->first;
		if(resnum >= start_of_residue_change)
		{
			resnum+=num_residues_inserted;
		}
		new_native_helix_residues.insert(
			std::pair<core::Size, utility::vector1<core::conformation::ResidueOP> >(resnum, map_it->second));
	}
	native_helix_residues = new_native_helix_residues;
}


///@brief update the loops to reflect
///the position changes after an insertion
void
SewingDesigner::update_anchors(
	utility::vector1<core::Size> & loop_anchors,
	protocols::loops::Loop const & new_loop,
	core::Size index_of_new_loop
){
	for(core::Size i=index_of_new_loop+1; i<=loop_anchors.size(); ++i)
	{
		loop_anchors_[i]+=new_loop.size();
	}
}

struct protein_segment
{
	Size start;
	Size end;

	bool operator<(const protein_segment & rhs) const
	{
		return start < rhs.start;
	}
};

///@brief rearrange the segments before and after the given
///loops file in order to minimize the sum of all loop distances
void
SewingDesigner::rearrange_pose(
	Pose & pose,
	NativeRotamersMap & native_helix_residues,
	utility::vector1<core::Size> & loop_anchors
){
	utility::vector1<protein_segment> protein_segments;
	for(core::Size i=1; i<=loop_anchors.size()+1; ++i)
	{
		protein_segment segment;
		if(i==1)
		{
			segment.start=1;
			segment.end=loop_anchors_[i];
		}
		else if(i==loop_anchors.size()+1)
		{
			segment.start=loop_anchors_[i-1]+1;
			segment.end=pose.total_residue();
		}
		else
		{
			segment.start=loop_anchors_[i-1]+1;
			segment.end=loop_anchors_[i];
		}
		protein_segments.push_back(segment);
	}

	utility::vector1<core::Size> dim_sizes;
//	if(num_helices_in_repeat_ > 0)
//	{
//		dim_sizes = utility::vector1<core::Size>(num_helices_in_repeat_, num_helices_in_repeat_);
//	}
//	else
//	{
		dim_sizes = utility::vector1<core::Size>(protein_segments.size(), protein_segments.size());
//	}
	core::Size n_dims = dim_sizes.size();
	Real min_distance_sum = 1000000000;
	utility::vector1<core::Size> rearrangement_order;
	for ( utility::LexicographicalIterator lex( dim_sizes ); ! lex.at_end(); ++lex )
	{
		bool valid=true;

		//Force either the first or second helix to be the first in the new bundle
		if(lex[1] != 1 && lex[1] != 2)
		{
			valid=false;
		}

		//Don't allow duplicate helix segments
		for(Size i=1; i<=n_dims; ++i)
		{
			for(Size j=i+1; j<=n_dims; ++j)
			{
				if(protein_segments[lex[i]].start == protein_segments[lex[j]].start)
				{
					valid=false;
				}
			}
		}

		if(valid)
		{
			Real dist_sq=0;
			utility::vector1<core::Size> new_order(n_dims);
			for(Size i=1; i<=n_dims-1; ++i)
			{
				dist_sq += pose.residue(protein_segments[lex[i]].end).atom("CA").xyz().distance_squared(
					pose.residue(protein_segments[lex[i+1]].start).atom("CA").xyz());
				new_order[i]=lex[i];
			}
			new_order[n_dims]=lex[n_dims];
//			if(num_helices_in_repeat_ != 0)//if this is a repeat then make sure to minize distance between repeating units
//			{
//				dist_sq += pose.residue(protein_segments[lex[n_dims]].end).atom("CA").xyz().distance_squared(
//					pose.residue(protein_segments[num_helices_in_repeat_+lex[1]].start).atom("CA").xyz());
//			}
			if(dist_sq < min_distance_sum)
			{
				min_distance_sum=dist_sq;
				rearrangement_order = new_order;
			}
		}
	}

	//Use the rearrangment order of the unique (non-repeat) helices to create order for the entire pose
	utility::vector1<protein_segment> best_order;
//	if(num_helices_in_repeat_ > 0)
//	{
//		//assert(protein_segments.size() % num_helices_in_repeat_ == 0);
//		core::Size num_repeats = protein_segments.size()/num_helices_in_repeat_;
//		for(core::Size i=0; i<num_repeats; ++i)
//		{
//			for(core::Size j=1; j<=rearrangement_order.size(); ++j)
//			{
//				core::Size adjusted_index = (num_helices_in_repeat_*i)+rearrangement_order[j];
//				best_order.push_back(protein_segments[adjusted_index]);
//			}
//		}
//	}
//	else
//	{
	for(core::Size j=1; j<=rearrangement_order.size(); ++j)
	{
		best_order.push_back(protein_segments[rearrangement_order[j]]);
	}
//	}

	TR << "Best order for pose:\n " << std::endl;
	for(Size i=1; i<=best_order.size(); ++i)
	{
		TR << best_order[i].start << " " << best_order[i].end << std::endl;
	}

	core::pose::Pose reordered_pose;
	NativeRotamersMap reordered_native_helix_residues;
//	protocols::loops::Loops reordered_loop;
	Size residue_counter=1;
	for(Size i=1; i<=best_order.size(); ++i)
	{
		protein_segment const & cur_frag = best_order[i];
		for(Size j=cur_frag.start; j<=cur_frag.end; ++j)
		{
			if(i==1 && j==cur_frag.start)
			{
				reordered_pose.conformation().append_residue_by_jump(pose.residue(j), 0);
			}
			else{
				conformation::remove_lower_terminus_type_from_conformation_residue(pose.conformation(), j);
				conformation::remove_upper_terminus_type_from_conformation_residue(pose.conformation(), j);
				reordered_pose.conformation().append_residue_by_bond(pose.residue(j));
			}

//			if(i!=1 && j==cur_frag.start)
//			{
//				reordered_loops.add_loop(reordered_pose.total_residue()-1, reordered_pose.total_residue());
//			}

			reordered_native_helix_residues[residue_counter] = native_helix_residues[j];
			++residue_counter;
		}
	}
	//add back final terminus type
	conformation::add_upper_terminus_type_to_conformation_residue(pose.conformation(), 1);
	conformation::add_lower_terminus_type_to_conformation_residue(pose.conformation(), pose.total_residue());
	pose = reordered_pose;
	native_helix_residues = reordered_native_helix_residues;
//	loops = reordered_loops;

	protocols::jd2::JobOP const job_me ( protocols::jd2::JobDistributor::get_instance()->current_job() );
	std::string const job_name ( protocols::jd2::JobDistributor::get_instance()->job_outputter()->output_name(job_me) );
	reordered_pose.dump_pdb(job_name+"_reordered.pdb");
	//dump_native_residue_file(native_helix_residues, job_name+"_reordered.rots");
}

void
SewingDesigner::parse_my_tag(
		TagCOP const tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const & movers,
		Pose const &)
{
	if(tag->hasOption("loop_creation_mover"))
	{
		string const lcm_mover_name = tag->getOption< string >( "loop_creation_mover" );
		protocols::moves::Movers_map::const_iterator lcm_it = movers.find( lcm_mover_name );
		if(lcm_it == movers.end())
		{
			utility_exit_with_message( "Mover " + lcm_mover_name + " not found" );
		}
		loop_creation_mover_ = dynamic_cast< devel::loop_creation::LoopCreationMover * >(lcm_it->second());
	}
	else
	{
		utility_exit_with_message("You must provide a 'loop_creation_mover' option to the SewingDesigner tag");
	}

	if(tag->hasOption("loop_anchors"))
	{
		string const loop_anchors_string = tag->getOption<string>("loop_anchors");
		utility::vector1<string> loop_anchor_strings=utility::string_split(loop_anchors_string, ',');
		for(core::Size i=1; i<=loop_anchor_strings.size(); ++i)
		{
			loop_anchors_.push_back(utility::string2int(loop_anchor_strings[i]));
		}
	}
	if(tag->hasOption("asym_size"))
	{
		asym_size_ = tag->getOption<core::Size>("asym_size");
	}
	else
	{
		utility_exit_with_message("You must provide a asym_size option to the SewingDesigner tag");
	}
}

void
SewingDesigner::detect_loop_anchors(
	Pose const & pose
){
	// squared distance at which bond is considered discontinuous
	Real const chain_break_cutoff = { 4.0 };

	for ( Size i = 1; i < pose.total_residue(); ++i )
	{
		Size j = i+1;
		conformation::Residue const & rsd = pose.residue(i);
		conformation::Residue const & next_rsd = pose.residue(j);
		if (rsd.is_polymer() && next_rsd.is_polymer())
		{
			Real dist_squared = rsd.atom( rsd.upper_connect_atom() ).xyz().distance_squared(next_rsd.atom( next_rsd.lower_connect_atom() ).xyz());
			if (dist_squared > chain_break_cutoff || dist_squared < 0.1)
			{
				loop_anchors_.push_back(i);
			}
		}
	}
}

} //sewing
} //devel
