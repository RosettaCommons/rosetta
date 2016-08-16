// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/delucasl/roc_optimizer.cc
/// @author Sam DeLuca

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

//basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

#include <basic/database/open.hh>


#include <core/optimization/ParticleSwarmMinimizer.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/qsar/scoring_grid/AtrGrid.hh>
#include <protocols/qsar/scoring_grid/RepGrid.hh>
#include <protocols/qsar/scoring_grid/VdwGrid.hh>
#include <protocols/qsar/scoring_grid/HbaGrid.hh>
#include <protocols/qsar/scoring_grid/HbdGrid.hh>
#include <protocols/qsar/qsarOptFunc.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <devel/init.hh>
#include <numeric/roc_curve.hh>

#include <protocols/ligand_docking/Transform.hh>
#include <protocols/ligand_docking/SlideTogether.hh>
#include <protocols/ligand_docking/InterfaceScoreCalculator.hh>
#include <protocols/features/ProteinSilentReport.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//protocol headers
#include <protocols/ligand_docking/StartFrom.hh>
#include <protocols/qsar/scoring_grid/GridManager.cc>

#include <basic/database/sql_utils.hh>
#include <basic/options/option_macros.hh>
#include <utility/io/izstream.hh>

//external headers
#include <cppdb/frontend.h>


//These things really sketch me out but all the cool kids are using them sooooo
OPT_1GRP_KEY(String,roc_opt,active_list)
OPT_1GRP_KEY(String,roc_opt,inactive_list)
//OPT_1GRP_KEY(Integer,roc_opt,outer_cycles)


static THREAD_LOCAL basic::Tracer roc_tracer( "ROC_optimizer" );


class GridWeights {

public:

	GridWeights(
		core::Real const & atr,
		core::Real const & rep,
		core::Real const & vdw,
		core::Real const & hba,
		core::Real const & hbd) :
		atr_(atr),
		rep_(rep),
		vdw_(vdw),
		hba_(hba),
		hbd_(hbd)
	{
		//
	}

	core::optimization::Multivec get_multivec()
	{
		core::optimization::Multivec current_weights;
		current_weights.push_back(atr_);
		current_weights.push_back(rep_);
		current_weights.push_back(vdw_);
		current_weights.push_back(hba_);
		current_weights.push_back(hbd_);
		return current_weights;
	}

	std::map<std::string,core::Size> get_indices()
	{
		std::map<std::string,core::Size> current_indices;
		current_indices.insert(std::make_pair("AtrGrid",1));
		current_indices.insert(std::make_pair("RepGrid",2));
		current_indices.insert(std::make_pair("VdwGrid",3));
		current_indices.insert(std::make_pair("HbaGrid",4));
		current_indices.insert(std::make_pair("HbdGrid",5));
		return current_indices;
	}

	core::Real atr() const
	{
		return atr_;
	}

	core::Real rep() const
	{
		return rep_;
	}

	core::Real vdw() const
	{
		return vdw_;
	}

	core::Real hba() const
	{
		return hba_;
	}

	core::Real hbd() const
	{
		return hbd_;
	}

private:

	core::Real atr_;
	core::Real rep_;
	core::Real vdw_;
	core::Real hba_;
	core::Real hbd_;
};


void setup_grid_manager(GridWeights const & /*weights*/)
{
	using namespace protocols::qsar::scoring_grid;
	GridManager* grid_manager( GridManager::get_instance() );
	grid_manager->reset();

	SingleGridOP atr( new AtrGrid() );
	SingleGridOP rep( new RepGrid() );
	SingleGridOP vdw( new VdwGrid() );
	SingleGridOP hba( new HbaGrid() );
	SingleGridOP hbd( new HbdGrid() );

	grid_manager->insert_grid("atr",atr);
	grid_manager->insert_grid("rep",rep);
	grid_manager->insert_grid("vdw",vdw);
	grid_manager->insert_grid("hba",hba);
	grid_manager->insert_grid("hbd",hbd);
	roc_tracer << "setup grids" <<std::endl;

}


protocols::moves::MoverOP setup_start_coords()
{
	utility::vector1<core::Real> coord_vect = basic::options::option[basic::options::OptionKeys::docking::ligand::start_from]();
	core::Vector start_coords(coord_vect[1],coord_vect[2],coord_vect[3]);

	protocols::ligand_docking::StartFromOP start_from( new protocols::ligand_docking::StartFrom() );
	start_from->add_coords(start_coords,"default");
	start_from->chain("X");

	return start_from;

}

protocols::moves::MoverOP setup_transform_mover()
{
	//TODO boo hard coding bad
	std::string chain("X");
	core::Real box_size(5.0);
	core::Real move_distance(1.0);
	core::Real angle(45.0);
	core::Size cycles(8000);
	core::Real temp(100.0);

	protocols::moves::MoverOP transform( new protocols::ligand_docking::Transform(chain,box_size,move_distance,angle,cycles,temp) );
	return transform;


}

protocols::moves::MoverOP setup_slide_mover()
{
	std::string chain("X");
	protocols::moves::MoverOP slide_together( new protocols::ligand_docking::SlideTogether(chain) );
	return slide_together;
}

protocols::moves::MoverOP setup_score_mover()
{
	std::vector<std::string> chains;
	chains.push_back("X");

	core::scoring::ScoreFunctionOP score_fxn(core::scoring::get_score_function_legacy( "score12prime.wts" ));

	protocols::ligand_docking::InterfaceScoreCalculatorOP score_mover( new protocols::ligand_docking::InterfaceScoreCalculator() );

	score_mover->chains(chains);
	score_mover->score_fxn(score_fxn);
	return score_mover;

}

protocols::moves::MoverOP setup_lowres_protocol()
{
	protocols::moves::SequenceMoverOP lowres_protocol_mover( new protocols::moves::SequenceMover() );

	lowres_protocol_mover->add_mover(setup_start_coords());
	lowres_protocol_mover->add_mover(setup_transform_mover());
	lowres_protocol_mover->add_mover(setup_slide_mover());
	lowres_protocol_mover->add_mover(setup_score_mover());

	return lowres_protocol_mover;

}

void setup_activity_table(
	utility::sql_database::sessionOP & db_session,
	std::string const & active_list_filename,
	std::string const & inactive_list_filename
)
{
	std::string schema =
		"CREATE TABLE IF NOT EXISTS structure_activity (\n"
		"\tinput_tag TEXT,\n"
		"\tactivity BOOLEAN);";

	cppdb::statement schema_statement( basic::database::safely_prepare_statement(schema,db_session));
	basic::database::safely_write_to_database(schema_statement);

	std::string insert_string = "INSERT INTO structure_activity VALUES(?,?);";
	cppdb::statement insert_statement(basic::database::safely_prepare_statement(insert_string,db_session));

	utility::io::izstream active_file;
	active_file.open(active_list_filename, std::ios_base::in);
	while ( !active_file.eof() )
			{
		std::string line;
		getline(active_file,line);
		if ( line.size() > 0 ) {
			insert_statement.bind(1,line);
			insert_statement.bind(2,true);
			basic::database::safely_write_to_database(insert_statement);
		}
	}
	active_file.close();

	utility::io::izstream inactive_file;
	inactive_file.open(inactive_list_filename, std::ios_base::in);
	while ( !inactive_file.eof() )
			{
		std::string line;
		getline(inactive_file,line);
		if ( line.size() >0 ) {
			insert_statement.bind(1,line);
			insert_statement.bind(2,false);
			basic::database::safely_write_to_database(insert_statement);
		}
	}
	inactive_file.close();

}


numeric::RocCurve setup_roc_curve(utility::sql_database::sessionOP & db_session, core::Real cutoff)
{
	numeric::RocCurve roc_curve;
	std::string select_string =
		"SELECT structures.tag, job_string_real_data.data_value, structure_activity.activity\n"
		"\tFROM structures\n"
		"\t\tINNER JOIN structure_activity ON structures.input_tag = structure_activity.input_tag\n"
		"\t\tINNER JOIN job_string_real_data ON structures.struct_id = job_string_real_data.struct_id\n"
		"\tWHERE\n"
		"\t\tjob_string_real_data.data_key= ?";

	cppdb::statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));

	select_statement.bind(1,"total_score_X");

	cppdb::result result(basic::database::safely_read_from_database(select_statement));
	while ( result.next() )
			{
		std::string tag;
		core::Real value = 0.0;
		int activity = 0; //

		result >> tag >> value >> activity;
		bool predicted(value < cutoff);
		roc_curve.insert_point(predicted,static_cast<bool>(activity),tag,value);
	}

	return roc_curve;

}

void dump_curve_to_db(
	utility::sql_database::sessionOP & db_session,
	core::Size const & curve_id,
	utility::vector1<std::pair<platform::Real, platform::Real> > const & curve )
{
	std::string schema =
		"CREATE TABLE IF NOT EXISTS roc_curve (\n"
		"\tpoint_id INTEGER PRIMARY KEY AUTO_INCREMENT,\n"
		"\tcurve_id INTEGER,\n"
		"\tTPR REAL,\n"
		"\tFPR REAL);";

	cppdb::statement schema_statement( basic::database::safely_prepare_statement(schema,db_session));
	basic::database::safely_write_to_database(schema_statement);

	std::string insert_string = "INSERT INTO roc_curve VALUES(NULL,?,?,?);";
	cppdb::statement insert_statement(basic::database::safely_prepare_statement(insert_string,db_session));
	insert_statement.bind(1,curve_id);
	for ( utility::vector1<std::pair<platform::Real, platform::Real> >::const_iterator curve_it = curve.begin(); curve_it != curve.end(); ++curve_it ) {
		insert_statement.bind(2,curve_it->first);
		insert_statement.bind(3,curve_it->second);
		basic::database::safely_write_to_database(insert_statement);
	}
}

GridWeights optimize_weights(core::optimization::Multivec min, core::optimization::Multivec max, core::optimization::Multifunc & fitness,core::optimization::Multivec initial_values)
{
	core::optimization::ParticleSwarmMinimizer roc_pso(min,max);
	core::optimization::ParticleOPs optimized_particles(roc_pso.run(100,fitness,100,initial_values));

	core::optimization::ParticleOP best_particle(optimized_particles[1]);
	core::optimization::Multivec best_scores = best_particle->pbest();
	return GridWeights(best_scores[1],best_scores[2],best_scores[3],best_scores[4],best_scores[5]);
}

void clean_up_database(utility::sql_database::sessionOP & db_session)
{

	std::string struct_select =
		"SELECT struct_id FROM structures";
	cppdb::statement select_statement(basic::database::safely_prepare_statement(struct_select,db_session));
	cppdb::result res(basic::database::safely_read_from_database(select_statement));

	protocols::features::ProteinSilentReport reporter;

	while ( res.next() )
			{
		protocols::features::StructureID struct_id;
		res >> struct_id;

		reporter.delete_pose(db_session,struct_id);

	}

}

void write_weights_to_db(utility::sql_database::sessionOP & db_session, GridWeights const & weights,core::Size const & cycle)
{
	std::string schema =
		"CREATE TABLE IF NOT EXISTS grid_weights (\n"
		"\tcycle INTEGER PRIMARY KEY,\n"
		"\thba REAL,\n"
		"\thbd REAL,\n"
		"\tvdw REAL,\n"
		"\tatr REAL,\n"
		"\trep REAL);";
	cppdb::statement schema_statement( basic::database::safely_prepare_statement(schema,db_session));
	basic::database::safely_write_to_database(schema_statement);

	std::string insert_string = "INSERT INTO grid_weights VALUES(?,?,?,?,?,?);";
	cppdb::statement insert_statement(basic::database::safely_prepare_statement(insert_string,db_session));
	insert_statement.bind(1,cycle);
	insert_statement.bind(2,weights.hba());
	insert_statement.bind(3,weights.hbd());
	insert_statement.bind(4,weights.vdw());
	insert_statement.bind(5,weights.atr());
	insert_statement.bind(6,weights.rep());
	basic::database::safely_write_to_database(insert_statement);

}

int main(int argc, char* argv[])
{
	try {

		NEW_OPT(roc_opt::active_list,"a list of active protein_ligand complexes","");
		NEW_OPT(roc_opt::inactive_list,"a list of inactive protein_ligand complexes","");
		//NEW_OPT(roc_opt::outer_cycles,"the number of dock/optimize cycles to perfrom","");

		devel::init(argc,argv);

		GridWeights current_weights(1.0,1.0,1.0,1.0,1.0);

		core::Size mpi_rank = 0;

		core::optimization::Multivec min_weights(5,-3.0);
		core::optimization::Multivec max_weights(5,3.0);

#ifdef USEMPI
	MPI_Comm_rank( MPI_COMM_WORLD, ( int* )( &mpi_rank ) );
#endif

		std::string active_list_filename = basic::options::option[basic::options::OptionKeys::roc_opt::active_list]();
		std::string inactive_list_filename = basic::options::option[basic::options::OptionKeys::roc_opt::inactive_list]();

		utility::sql_database::sessionOP db_session(basic::database::get_db_session());

		if ( mpi_rank == 0 ) {
			roc_tracer <<"setting up activity table"<<std::endl;
			setup_activity_table(db_session,active_list_filename,inactive_list_filename);
			roc_tracer << "done setting up activity table" <<std::endl;
		}


		//resync processes after table setup

#ifdef USEMPI
		MPI_Barrier( MPI_COMM_WORLD );
#endif

		for ( core::Size cycle = 1; cycle <= 10; ++cycle ) {

			setup_grid_manager(current_weights);
			protocols::moves::MoverOP mover(setup_lowres_protocol());
			protocols::jd2::JobDistributor::get_instance()->mpi_finalize(false);
			roc_tracer << "starting docking for cycle " << cycle << std::endl;
			protocols::jd2::JobDistributor::get_instance()->go(mover);

#ifdef USEMPI
		MPI_Barrier( MPI_COMM_WORLD );
#endif


			if ( mpi_rank == 0 ) {
				roc_tracer << "starting roc optimization for cycle " << cycle << std::endl;
				protocols::qsar::qsarOptFunc qsar_func(db_session,current_weights.get_multivec(),current_weights.get_indices());
				numeric::RocCurve curve(setup_roc_curve(db_session,0.0));
				curve.generate_roc_curve();
				dump_curve_to_db(db_session,cycle,curve.roc_curve());
				roc_tracer << curve.calculate_auc() <<std::endl;

				//query the database and process score information
				qsar_func.setup_data_map();
				current_weights = optimize_weights(min_weights,max_weights,qsar_func,current_weights.get_multivec());
				roc_tracer <<"New weights: " << current_weights.hba() << " "
					<<current_weights.hbd() << " "
					<<current_weights.atr() << " "
					<<current_weights.rep() << " "
					<<current_weights.vdw() << std::endl;
				write_weights_to_db(db_session,current_weights,cycle);
				clean_up_database(db_session);
			}

			protocols::jd2::JobDistributor::get_instance()->restart();

#ifdef USEMPI
		MPI_Barrier( MPI_COMM_WORLD );
#endif
		}

		//one last sync and finalize to shut everything down
#ifdef USEMPI
		MPI_Barrier( MPI_COMM_WORLD );
		MPI_Finalize();
#endif
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
