__kernel void Check_for_intersection(
	__global const float4 *rays,
	__global const float4 *atoms,
	__global float *distances,
	const unsigned int NUM_ATOMS
)
{
	const int large_dist = 9999.0;
	register int rayID = get_global_id(0);
	register int atomID, i, j, l;
	register float dirX, dirY, dirZ;
	register float distance = 9999.0;
	register float4 ray;
	
	ray = rays[rayID];
	dirX = large_dist*sin(ray.x)*cos(ray.y);
	dirY = large_dist*sin(ray.x)*sin(ray.y);
	dirZ = large_dist*cos(ray.x);
	
	// setup our quadratic equation
	float a = (dirX*dirX) + (dirY*dirY) + (dirZ*dirZ);

	for(atomID = 0; atomID < NUM_ATOMS; atomID ++) {

			float4 atom = atoms[atomID];
			float b = 2.0 * ( (dirX*(-atom.x)) + (dirY*(-atom.y)) + (dirZ*(-atom.z)) );
			float c = atom.x*atom.x + atom.y*atom.y + atom.z*atom.z - (atom.w * atom.w);
	
			// test for intersection
			float inside_sq = ( b * b ) - ( 4.0 * a * c );

			if (inside_sq > 0) {
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

				if(dist1_sq < distance)
					distance = dist1_sq;
				if(dist2_sq < distance)
					distance = dist2_sq;
			}				
	}
	if(distance < 9999.0) {
		distance = sqrt(distance) - ray.z;
		if(distance < 0)
			distance = -distance;
		else
			distance *= 3.12; // check this
	} else
		distance = 0;

	distances[rayID] = distance;
}
