<loop_design>

     <io infile="bCD_UR1.pdb" output_prefix=ttt nstruct=1/>

     <param max_cycles=1 max_cycles_start=1 max_cycles_start_3mer=1 max_cycles_start_3mer_1mer=1 max_cycles_start_3mer_anchor=1 max_cycles_start_1mer=1 max_cycles_rand=1 max_cycles_1mer=1 max_cycles_rep=1 max_cycles_ramp=1 max_cycles_ramp_sm=1 max_cycles_ramp_sm_min=1 max_cycles_design=2 max_cycles_design_sm=1 max_cycles_design_sm_min=1 />

     <anchored_loop name=x nres_pre=3 nres_post=3 vary_pre=2 vary_post=2 ss=LLLLLLLL>
	  <begin res_id=A:150/>
	  <end res_id=A:157/>
	  <from res_id=L:1 atom=O1/>
	  <to res_type=GLN atom=1HE2/>
	  <template filename="bCD_UR1.pdb" from=L:1 to=A:153 to_atom=1HE2/>
     </anchored_loop>

</loop_design>
