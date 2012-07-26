#define MIN(a,b) ((a) < (b) ? (a): (b))
#define MAX(a,b) ((a) > (b) ? (a): (b))

__kernel void TrimPeripheralBand(
	__global const float4 *dAccDotCoords, 
	const uint nAcc, 
	__global const float4 *dBurDotCoords, 
	__global char *dDotColl, 
	const float r2)
{
	int i, j, l;
	float4 dot1;
	__local char sColl[512];
	__local float4 sCoords[512];

	sColl[get_local_id(0)] = 0;
 	dot1 = dBurDotCoords[get_global_id(0)];

	for(i = 0; i < nAcc; i += get_local_size(0)) {
		barrier(CLK_LOCAL_MEM_FENCE);
		sCoords[get_local_id(0)] = dAccDotCoords[i + get_local_id(0)];
		barrier(CLK_LOCAL_MEM_FENCE);

		l = MIN(nAcc - i, get_local_size(0));
		for(j = 0; j < l; j++) {
			float4 dot2 = sCoords[j];
			dot2.x -= dot1.x;
			dot2.y -= dot1.y;
			dot2.z -= dot1.z;
			sColl[get_local_id(0)] |= (dot2.x*dot2.x + dot2.y*dot2.y + dot2.z*dot2.z) <= r2;
		}
	}

	dDotColl[get_global_id(0)] = sColl[get_local_id(0)];
}

__kernel void FindClosestNeighbor(
	__global float4 *dMyDotCoords, 
	__global float4 *dTheirDotCoords, 
	uint nTheirDots, 
	__global uint *dNeighbors)
{
	int i, j, l;
	float4 dot1;
	__local uint sNeighbors[512];
	__local float4 sCoords[512];
	float distmin = 99999.0, d2;

 	dot1 = dMyDotCoords[get_global_id(0)];
	sNeighbors[get_local_id(0)] = 42;

	for(i = 0; i < nTheirDots; i += get_local_size(0)) {
		barrier(CLK_LOCAL_MEM_FENCE);
		sCoords[get_local_id(0)] = dTheirDotCoords[i + get_local_id(0)];
		barrier(CLK_LOCAL_MEM_FENCE);

		l = MIN(nTheirDots - i, get_local_size(0));
		for(j = 0; j < l; j++) {
			float4 dot2 = sCoords[j];
			dot2.x -= dot1.x;
			dot2.y -= dot1.y;
			dot2.z -= dot1.z;
			d2 = dot2.x*dot2.x + dot2.y*dot2.y + dot2.z*dot2.z;
			if(d2 <= distmin) {
				distmin = d2;
				sNeighbors[get_local_id(0)] = i+j;
			}
		}
	}
	dNeighbors[get_global_id(0)] = sNeighbors[get_local_id(0)];
}
