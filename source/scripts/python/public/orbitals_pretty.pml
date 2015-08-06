cmd.color(52, "*")
util.cnc("*")
select orbitals, name LP*
unbond LP*, *
set sphere_scale, 0.1
cmd.show("spheres"   ,"orbitals")
unbond orbitals, *
cmd.color(5263,"orbitals")
hide everything
show cartoon
