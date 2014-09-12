__kernel void get_ray_score(
	   __global float4 *rays,	
	      __global float4 *atomcoords,
	          __global float4 *particles_translation_offset,
		       __global float4 *particles_rotation_offset,
		              __global float4 *ligCoM,
			               __global float *weights,
				                __global float4 *RAYorigins,
   							  __global float *raytest_scores,
							     unsigned int NUM_RAYS
 )
{
  int rayID = get_global_id(0);
  if(rayID >= NUM_RAYS) return;

  //scoring variables
  float const missing_point_weight = weights[0];
  float const steric_weight = weights[1];
  float const extra_point_weight = weights[2];
int const NUM_PARTICLES = (int)weights[3];

  float4 ray = rays[rayID];

    int origin_index = (int)ray.w;
    float4 ray_origin = RAYorigins[origin_index];

  const int large_dist = 999.;
  float best_rho_sq, min_intersect_SQ, dist_deviation;
  //quadratic setup starts for ray
  float dirX = large_dist*sin(ray.x)*cos(ray.y);
  float dirY = large_dist*sin(ray.x)*sin(ray.y);
  float dirZ = large_dist*cos(ray.x);
  // setup our quadratic equation
  float a = (dirX*dirX) + (dirY*dirY) + (dirZ*dirZ);

  for(int particleID = 0; particleID < NUM_PARTICLES; particleID++){
    int ray_score_ID = (particleID*NUM_RAYS+rayID);
    best_rho_sq = 9999.;

float4 rot_ang = particles_rotation_offset[particleID];
float4 trans_offset = particles_translation_offset[particleID];

int ligconfID = (int)trans_offset.w;
float4 ligand_CoM = ligCoM[ligconfID];
int const NUM_ATOMS = (int)ligand_CoM.w;

//tot_rot_mat is the total rotation matrix which is the multiplication of z*y*x rotation matrix
float tot_rot_mat_0 = cos(rot_ang.y)*cos(rot_ang.z);
float tot_rot_mat_1 = (cos(rot_ang.z)*sin(rot_ang.x)*sin(rot_ang.y)) + (-1*cos(rot_ang.x)*sin(rot_ang.z));
float tot_rot_mat_2 = (cos(rot_ang.x)*cos(rot_ang.z)*sin(rot_ang.y)) + (sin(rot_ang.x)*sin(rot_ang.z));
float tot_rot_mat_3 = cos(rot_ang.y)*sin(rot_ang.z);
float tot_rot_mat_4 = (cos(rot_ang.x)*cos(rot_ang.z)) + (sin(rot_ang.x)*sin(rot_ang.y)*sin(rot_ang.z));
float tot_rot_mat_5 = (cos(rot_ang.x)*sin(rot_ang.y)*sin(rot_ang.z)) + (-1*cos(rot_ang.z)*sin(rot_ang.x));
float tot_rot_mat_6 = sin(rot_ang.y)*(-1);
float tot_rot_mat_7 = cos(rot_ang.y)*sin(rot_ang.x);
float tot_rot_mat_8 = cos(rot_ang.x)*cos(rot_ang.y);

//rotate and translate the ligand pose
    int atom_ID_start = ligconfID*NUM_ATOMS;
    int atom_ID_stop = (ligconfID+1)*NUM_ATOMS;
    for(int atomID = atom_ID_start; atomID < atom_ID_stop; atomID ++) {

float4 ligatom = atomcoords[atomID];
//move pose such that ligand CoM is at the origin
ligatom.x -= ligand_CoM.x;
ligatom.y -= ligand_CoM.y;
ligatom.z -= ligand_CoM.z;
//apply_rotation offset
float ligatom_x = tot_rot_mat_0*ligatom.x + tot_rot_mat_1*ligatom.y + tot_rot_mat_2*ligatom.z;
float ligatom_y = tot_rot_mat_3*ligatom.x + tot_rot_mat_4*ligatom.y + tot_rot_mat_5*ligatom.z;
float ligatom_z = tot_rot_mat_6*ligatom.x + tot_rot_mat_7*ligatom.y + tot_rot_mat_8*ligatom.z;
//move pose back to original ligand CoM
ligatom.x = ligatom_x + ligand_CoM.x;
ligatom.y = ligatom_y + ligand_CoM.y;
ligatom.z = ligatom_z + ligand_CoM.z;
//Apply translation offset
ligatom.x -= (ray_origin.x + trans_offset.x);
ligatom.y -= (ray_origin.y + trans_offset.y);
ligatom.z -= (ray_origin.z + trans_offset.z);

      //quadratic setup starts for atom
      min_intersect_SQ = 9999.;
      float b = 2.0 * ( (dirX*(-ligatom.x)) + (dirY*(-ligatom.y)) + (dirZ*(-ligatom.z)) );
      float c = ligatom.x*ligatom.x + ligatom.y*ligatom.y + ligatom.z*ligatom.z - (ligatom.w * ligatom.w);
      // test for intersection
      float inside_sq = ( b * b ) - ( 4.0 * a * c );
if(inside_sq > 0.0) {
        float inside = sqrt(inside_sq);
        float mu1 = -(b-inside) / ( 2.0 * a);
        float x1 =  mu1 * dirX;
        float y1 =  mu1 * dirY;
        float z1 =  mu1 * dirZ;
        float dist1_sq = x1*x1 + y1*y1 + z1*z1;
        float mu2 = -(b+inside) / ( 2.0 * a);
        float x2 = mu2 * dirX;
        float y2 = mu2 * dirY;
        float z2 = mu2 * dirZ;
        float dist2_sq = x2*x2 + y2*y2 + z2*z2;
	if(dist2_sq < dist1_sq){
	  min_intersect_SQ = dist2_sq;
	  }else{
	    min_intersect_SQ = dist1_sq;
	    }
      }
      if(min_intersect_SQ < best_rho_sq ){
        best_rho_sq  = min_intersect_SQ;
      }
    }

    float plaid_rho = 9999.;
    if(best_rho_sq < 9998.0){
      plaid_rho = sqrt(best_rho_sq );
    }
    if( (plaid_rho > 9998.0) && (ray.z > 0.001) ) {
      raytest_scores[ray_score_ID] = missing_point_weight;
    }else if( (plaid_rho < 9998.0) && (ray.z < 0.001) ) {
      raytest_scores[ray_score_ID] = extra_point_weight;
    }else if( (plaid_rho < 9998.0) && (ray.z > 0.001) ) {
      dist_deviation = plaid_rho - ray.z;
      if(dist_deviation < 0.0){
        dist_deviation = (ray.z - plaid_rho)*steric_weight;
      }
      raytest_scores[ray_score_ID] = dist_deviation;
    }else{
      raytest_scores[ray_score_ID] = 0.0;
    }
 }
}

__kernel void calc_electrostatics(
	__global float4 *eatoms,
	__global float *elec_scores,
        __global float4 *gpu_espGrid,
	__global float4 *gpu_typGrid,
	__global float *grid_dim,
        __global float *grid_mid,
        const unsigned int NUM_EATOMS
    	       	   )         
  	       {

int espID = get_global_id(0);																		       
 if(espID >= NUM_EATOMS) return;

  unsigned int dim_x = grid_dim[0];
  unsigned int dim_y = grid_dim[1];
  unsigned int dim_z = grid_dim[2];
float mid_x = grid_mid[0];
float mid_y = grid_mid[1];
float mid_z = grid_mid[2];
float spacing = grid_mid[3];

 float atm_esp_energy = 0;
 float4 atm = eatoms[espID];      
      float grid_coord_x = (atm.x - ( mid_x - ( (float)((dim_x-1)/2) * spacing) ))/spacing;
      float grid_coord_y = (atm.y - ( mid_y - ( (float)((dim_y-1)/2) * spacing) ))/spacing;
      float grid_coord_z = (atm.z - ( mid_z - ( (float)((dim_z-1)/2) * spacing) ))/spacing;

      float X = grid_coord_x;
      float Y = grid_coord_y;
      float Z = grid_coord_z;;

      float X1 = floor(X+1);
      float Y1 = floor(Y+1);
      float Z1 = floor(Z+1);

      float X0 = ceil(X-1);
      float Y0 = ceil(Y-1);
      float Z0 = ceil(Z-1);
      if( ((X0 >= dim_x)||(X0 < 0)) || ((Y0 > dim_y)||(Y0 < 0)) || ((Z0 >= dim_z)||(Z0 < 0)) || ((X1 >= dim_x)||(X1 < 0)) || ((Y1 >= dim_y)||(Y1 < 0)) || ((Z1 >= dim_z)||(Z1 < 0)) ) {
      atm_esp_energy = 100;
      } else {

      float Xd = (X-X0)/(X1-X0);
      float Yd = (Y-Y0)/(Y1-Y0);
      float Zd = (Z-Z0)/(Z1-Z0);
      int inx_000 = (X0*dim_y*dim_z) + (Y0*dim_z) + Z0;
      int inx_010 = (X0*dim_y*dim_z) + (Y1*dim_z) + Z0;
      int inx_001 = (X0*dim_y*dim_z) + (Y0*dim_z) + Z1;
      int inx_011 = (X0*dim_y*dim_z) + (Y1*dim_z) + Z1;
      int inx_100 = (X1*dim_y*dim_z) + (Y0*dim_z) + Z0;
      int inx_110 = (X1*dim_y*dim_z) + (Y1*dim_z) + Z0;
      int inx_101 = (X1*dim_y*dim_z) + (Y0*dim_z) + Z1;
      int inx_111 = (X1*dim_y*dim_z) + (Y1*dim_z) + Z1;

      float4 espGrid_inx_000 = gpu_espGrid[inx_000];
      float4 espGrid_inx_010 = gpu_espGrid[inx_010];
      float4 espGrid_inx_001 = gpu_espGrid[inx_001];
      float4 espGrid_inx_011 = gpu_espGrid[inx_011];
      float4 espGrid_inx_100 = gpu_espGrid[inx_100];
      float4 espGrid_inx_110 = gpu_espGrid[inx_110];
      float4 espGrid_inx_101 = gpu_espGrid[inx_101];
      float4 espGrid_inx_111 = gpu_espGrid[inx_111];

      float4 typGrid_inx_000 = gpu_typGrid[inx_000];
      float4 typGrid_inx_010 = gpu_typGrid[inx_010];
      float4 typGrid_inx_001 = gpu_typGrid[inx_001];
      float4 typGrid_inx_011 = gpu_typGrid[inx_011];
      float4 typGrid_inx_100 = gpu_typGrid[inx_100];
      float4 typGrid_inx_110 = gpu_typGrid[inx_110];
      float4 typGrid_inx_101 = gpu_typGrid[inx_101];
      float4 typGrid_inx_111 = gpu_typGrid[inx_111];

      float C00 = (espGrid_inx_000.x * (1-Xd)) + (espGrid_inx_100.x * Xd);
      float C10 = (espGrid_inx_010.x * (1-Xd)) + (espGrid_inx_110.x * Xd);
      float C01 = (espGrid_inx_001.x * (1-Xd)) + (espGrid_inx_101.x * Xd);
      float C11 = (espGrid_inx_011.x * (1-Xd)) + (espGrid_inx_111.x * Xd);

      float C0 = C00*(1-Yd) + C10*Yd;
      float C1 = C01*(1-Yd) + C11*Yd;
      float C = C0*(1-Zd) + C1*Zd;

atm_esp_energy = C * atm.w;

	              if ( ( gpu_typGrid[inx_000].x < 0.1 ) ||
   		       ( typGrid_inx_010.x < 0.1 ) ||
		        ( typGrid_inx_001.x < 0.1 ) ||
     			 ( typGrid_inx_011.x < 0.1 ) ||
     			  ( typGrid_inx_100.x < 0.1 ) ||
 			   ( typGrid_inx_110.x < 0.1 ) ||
			    ( typGrid_inx_101.x < 0.1 ) ||
   			     ( typGrid_inx_111.x < 0.1 ) ) {
			      if (atm_esp_energy < 0.) atm_esp_energy = 0.;
			       }
			        }
				 elec_scores[espID] = atm_esp_energy;
}

__kernel void Get_scores(
                         __global float *raytest_scores,
                         __global float *elec_scores,
                         __global float *elec_weights,
                         __global float *particletest_scores,
                  	  const unsigned int NUM_ATOMS,
                         const unsigned int NUM_RAYS,
         		 const unsigned int NUM_PARTICLES
                         )
{
  int particleID = get_global_id(0);
  if(get_global_id(0) >= NUM_PARTICLES) return;

//get darc score
  int start = particleID*NUM_RAYS;
  int stop = start + NUM_RAYS;
  int ray_score_ID;
  float score = 0;
  int num = 0;
  for(ray_score_ID=start; ray_score_ID<stop; ray_score_ID++){
    if(raytest_scores[ray_score_ID]>0.0){
      score += raytest_scores[ray_score_ID];
      num++;
    }
  }

//get electrostatics score
  int esp_start = particleID*NUM_ATOMS;
  int esp_stop = esp_start + NUM_ATOMS;
  int esp_score_ID;
  float const esp_weight = elec_weights[0];
  float esp_score = 0;
  for(esp_score_ID=esp_start; esp_score_ID<esp_stop; esp_score_ID++){
      esp_score += elec_scores[esp_score_ID];
      }
//add darc score and electrostatics score
//particletest_scores[particleID] = esp_score * esp_weight;
//particletest_scores[particleID] = score/num;
  particletest_scores[particleID] = (score/num) + (esp_score * esp_weight);
  
}
