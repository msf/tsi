#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>

#define VER "3.0"

using namespace std;


void usage(char *name)
{
	printf("tsi-convert %s\n Program to convert tsi-binary grid format into ascii and vice-versa.\n\n",VER);
	printf("usage:\t %s grid.tsi grid.out\n\t %s grid.out grid.tsi\n",name,name);
}



float * gslib_grid_read(char *in_file, int *nx, int *ny, int *nz)
{
	char *name;
	float val;
	float *data;
	int tmp, i;
	int xlen, ylen, zlen;

	ifstream file( in_file, ifstream::in );

	file >> name; // discard name
	file >> xlen >> ylen >> zlen;
	file >> tmp;
	if(tmp != 1) //only 1 property please
		return NULL;
	file >> name; // discard property name
	
	cout << "Reading grid with size: " 
		<< xlen << ", " << ylen << ", " << zlen << endl;

	int size = xlen * ylen * zlen;

	data = new float[size];

	while(file.good() && i < size) {
		file >> val;
		data[i++] = val;
	}

	file.close();

	*nx = xlen;
	*ny = ylen;
	*nz = zlen;

	return data;
}

int gslib_grid_write(char *dst, float *values, int x, int y, int z)
{
	return 0;
}

float * sgems_grid_read(char *src, int *x, int *y, int *z)
{
	return NULL;
}


int	sgems_grid_write(char *out_file, float *values, int x, int y, int z)
{
	int magic = 0x11b25d17;
	int version = 100;
	float p1 = 1.0, p0 = 0.0;

	ofstream file(out_file, ofstream::binary);

	/* write sgems "magic header" */
	file << 0x5d17 << 0x11b2;
	file << "Cgrid"; // cartesian grid
	file << "tsi"; // grid name
	file << version;
	file << x << y << z; //size of grid
	file << p1 << p1 << p1; //size for cell 
	file << p0 << p0 << p0; //origin coords
	file << 1; // only one property per cell
	file << "data"; // property name
	file.write( (char *) values, sizeof(float) * x*y*z);

	file.close();

	return 0;
}
	
	


int main(int argc, char *argv[])
{
	int x,y,z;
	float *values = NULL;
	char *src, *dst;

	if(argc != 4){
		usage(argv[0]);
		return 1;
	}

	src = argv[2];
	dst = argv[3];

	if( argv[1][0] == 'a' ) {
		/* ascii gslib to sgems binary format */
		values = gslib_grid_read(src, &x, &y, &z);

		sgems_grid_write(dst, values, x, y, z );
	} else {
		values = sgems_grid_read(src, &x, &y, &z);

		gslib_grid_write(dst, values, x, y, z);
	}
	delete[] values;
	

	return 0;
}



