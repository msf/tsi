#include <iostream>
#include <fstream>
#include <vector>

#include <QFile>
#include <QDataStream>


#include <stdio.h>
#include <stdlib.h>

#define VER "3.0"

using namespace std;


void usage(char *name)
{
	printf("tsi-convert %s\n Program to convert tsi-binary grid format into ascii and vice-versa.\n\n",VER);
	printf("usage:\t %s grid.tsi grid.out\n\t %s grid.out grid.tsi\n",name,name);
}

int swap_bytes(int in)
{
	int out;
	char *a = (char *) &in;
	char *b = (char *) &out;
	b[3] = a[0];
	b[0] = a[3];

	b[2] = a[1];
	b[1] = a[2];
	return out;
}



float * gslib_grid_read(char *in_file, int *nx, int *ny, int *nz)
{
	char *name = NULL;
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
	char *type;
	char *grid_name;
	uint32_t magic_nb;
	uint32_t nx, ny, nz;
	float xsiz, ysiz, zsiz;
	float ox, oy, oz;
	int version;

	cout << "IN: " << src << endl;
//	ifstream file(src, ifstream::binary);
//
//
//	file >> magic_nb;
//	printf("magic_nb: %d: sw:%d\n", magic_nb, ntohl(magic_nb));
//
//	file >> type;
//	printf("type: %s\n", type);
//	//delete[] type;
//
//	file >> grid_name;
//	printf("grid_name: %s\n", grid_name);
//
//	file >> version;
//	printf("version: %d\n", version);
//
//
//	file >> nx >> ny >> nz;
//	file >> ox >> oy >> oz;
//	file >> xsiz >> ysiz >> zsiz;
//
//	printf("dimensions: %d.%d.%d orig: %f.%f.%f cell: %f.%f.%f\n",
//			nx, ny, nz,
//			ox, oy, oz,
//			xsiz, ysiz, zsiz);
//
//	printf("QT CENAS\n");

	QFile qfile(src);
	if( !qfile.open(QIODevice::ReadOnly) )
		return NULL;

	QDataStream file(&qfile);

	file >> magic_nb;
	printf("magic_nb: %d: sw:%d\n", magic_nb, swap_bytes(magic_nb));

	file >> type;
	printf("type: %s\n", type);

	file >> grid_name;
	printf("grid: %s\n", grid_name);

	file >> version;
	printf("version: %d\n", version);


	file >> nx >> ny >> nz;
	file >> ox >> oy >> oz;
	file >> xsiz >> ysiz >> zsiz;

	printf("dimensions: %d.%d.%d orig: %f.%f.%f cell: %f.%f.%f\n",
			nx, ny, nz,
			ox, oy, oz,
			xsiz, ysiz, zsiz);
	int prop_count;
	file >> prop_count;
	printf("prop_count: %d\n", prop_count);

	vector<char *> prop_names( prop_count );
	vector<float *> properties( prop_count );

	for(int i = 0; i < prop_count; i++) {
		file >> prop_names[i];

		int grid_size = nx * ny * nz;
		properties[i] = new float[grid_size];
		float val;
		for(int j = 0; i < grid_size; j++) {
			file >> val;
			properties[i][j];
		}

	}

	float *grid = (float *) properties[0];
	*x = nx;
	*y = ny;
	*z = nz;
	
	return grid;
}


int	sgems_grid_write(char *out_file, float *values, int x, int y, int z)
{
	int magic = 2987464541;
	int version = 100;
	float p1 = 1.0, p0 = 0.0;

	ofstream file(out_file, ofstream::binary);

	/* write sgems "magic header" */
	file << magic;
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

	setvbuf(stdout, NULL, _IONBF, 0);

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



