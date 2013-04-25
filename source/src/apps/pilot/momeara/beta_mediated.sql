
--find all hydrogen bonds at beta sheet mediated protein interfaces
--where the structures are taken from the richardson's top5200 dataset
--run with the command
--$sqlite3 structure_data_db_091002.sqlite < beta_mediated.sql 
select
	hbonds.energy,
	hbonds.HBEvalType,
	structures.fileName
--	geom.AHdist,
--	geom.cosBAH,
--	geom.cosAHD,
--	don_sites.chain,
--	acc_sites.chain,
--	don_env.dssp,
--	acc_env.dssp
from
	hbond_geom_coords as geom,
	hbonds,
	structures,
	sites as don_sites,
	sites as acc_sites,
	site_environment as don_env,
	site_environment as acc_env
where
	geom.hbond_id=hbonds.hbond_id and
	hbonds.struct_id=structures.struct_id and	
	hbonds.don_id=don_sites.site_id and
	hbonds.acc_id=acc_sites.site_id and
	hbonds.don_id=don_env.site_id and
	hbonds.acc_id=acc_env.site_id and
	don_env.dssp="E" and
	acc_env.dssp="E" and
	don_sites.chain!=acc_sites.chain;