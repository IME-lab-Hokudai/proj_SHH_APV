#include "cubeMap_SH.h"
#include "sphericalHarmonics.h"
//-----------------------------------------------------------------------
float	*deltaFormFactor[6];
float	*cubeSHTable[6];
int		 shOrder = -1;
//-----------------------------------------------------------------------
//void sphericalHarmonics16( const float x, const float y, const float z, float sh[ 16 ] )
//{
//	// A is normalization factor: A = sqrt( ( 2 * l + 1 ) / ( 4 * PI ) ) )
//	float c[ 16 ];
//
//	c[ 0] =  0.282095f;
//	c[ 1] = -0.488603f;
//	c[ 2] =  0.488603f;
//	c[ 3] = -0.488603f;
//
//	c[ 4] =  1.092548f;
//	c[ 5] = -1.092548f;
//	c[ 6] =  0.315392f;
//	c[ 7] = -1.092548f;
//	c[ 8] =  0.546274f;
//
//	c[ 9] = -0.590044f;			// Y{3, -3}
//	c[10] =  2.89061f;			// Y{3, -2}
//	c[11] = -0.457046f;			// Y{3, -1}
//	c[12] =  0.373176f;			// Y{3,  0}
//	c[13] = -0.457046f;			// Y{3,  1}
//	c[14] =  1.44531f;			// Y{3,  2}
//	c[15] = -0.590044f;			// Y{3,  3}
//
//	//for(int i=0; i<16; i++) c[i] = 1.0;
//
//	// Y_{0,0} 
//	sh[ 0] = c[ 0];
//
//	// Y_{1,-1} Y_{1,0}, Y_{1,1}
//	sh[ 1] = c[ 1] * y;
//	sh[ 2] = c[ 2] * z;
//	sh[ 3] = c[ 3] * x;
//
//	// Y_{2, -2} Y_{2,-1}, Y_{2,1}
//	sh[ 4] = c[ 4] * x * y;
//	sh[ 5] = c[ 5] * y * z;
//	sh[ 7] = c[ 7] * x * z;
//
//	// Y_{2,0} 
//	sh[ 6] = c[ 6] * (3.0f * z * z - 1.0f);
//
//	// Y_{2,2} 
//	sh[ 8] = c[ 8] * (x * x - y * y);
//
//	// Y_{3, -3} = A * sqrt(5/8) * (3 * x^2 * y - y^3)
//	sh[ 9] = c[ 9] * (3.0f * x * x * y - y * y * y); 
//
//	// Y_{3, -2} = A * sqrt(15) * x * y * z 
//	sh[10] = c[10] * x * y * z;
//
//	// Y_{3, -1} = A * sqrt(3/8) * y * (5 * z^2 - 1)
//	sh[11] = c[11] * y * (5.0f * z * z - 1.0f);
//
//	// Y_{3,  0} = A * (1/2) * (5 * z^3 - 3 *z)	
//	sh[12] = c[12] * (5.0f * z * z * z - 3 * z);
//
//	// Y_{3,  1} = A * sqrt(3/8) * x * (5 * z^2 - 1)
//	sh[13] = c[13] * x * (5.0f * z * z  - 1.0f);
//
//	// Y_{3,  2} = A * sqrt(15/4) * z *(x^2 - y^2)
//	sh[14] = c[14] * z * (x * x - y * y);
//
//	// Y_{3,  3} = A * sqrt(5/8) * (x^3 - 3 * x * y^2)
//	sh[15] = c[15] * (x * x * x - 3.0f * x * y * y);
//}

//-----------------------------------------------------------------------
void getPositiveXPoint(int s, int t, int width, float vec[3]) {
	vec[ 0 ] = 1.0f;	/* x */
	vec[ 1 ] = -( 2.0f * ( ( (float)s + 0.5f ) / (float)width) - 1.0f );
	vec[ 2 ] =  ( 2.0f * ( ( (float)t + 0.5f ) / (float)width) - 1.0f );
}
void getNegativeXPoint(int s, int t, int width, float vec[3]) {
	vec[ 0 ] =   -1.0f;	/* x */
	vec[ 1 ] =  ( 2.0f * ( ( (float)s + 0.5f ) / (float)width ) - 1.0f );
	vec[ 2 ] =  ( 2.0f * ( ( (float)t + 0.5f ) / (float)width ) - 1.0f );
}
void getPositiveYPoint(int s, int t, int width, float vec[3]) {
	vec[ 1 ] = 1.0f;	/* y */
	vec[ 0 ] =  ( 2.0f * ( ( (float)s + 0.5f ) / (float)width ) - 1.0f );
	vec[ 2 ] =  ( 2.0f * ( ( (float)t + 0.5f ) / (float)width ) - 1.0f );
}
void getNegativeYPoint(int s, int t, int width, float vec[3]) {
	vec[ 1 ] = -1.0f;	/* y */
	vec[ 0 ] = -( 2.0f * ( ( (float)s + 0.5f ) / (float)width ) - 1.0f );
	vec[ 2 ] =  ( 2.0f * ( ( (float)t + 0.5f ) / (float)width ) - 1.0f );
}
void getPositiveZPoint(int s, int t, int width, float vec[3]) {
	vec[ 2 ] = 1.0f;	/* z */
	vec[ 0 ] = -( 2.0f * ( ( (float)t + 0.5f ) / (float)width ) - 1.0f );
	vec[ 1 ] = -( 2.0f * ( ( (float)s + 0.5f ) / (float)width ) - 1.0f );
}
void getNegativeZPoint(int s, int t, int width, float vec[3]) {
	vec[ 2 ] = -1.0f;	/* -z */
	vec[ 0 ] =  ( 2.0f * ( ( (float)t + 0.5f ) / (float)width ) - 1.0f );
	vec[ 1 ] = -( 2.0f * ( ( (float)s + 0.5f ) / (float)width ) - 1.0f );
}
//-----------------------------------------------------------------------

void calDeltaFormFactor(int size) {
	float vec[ 3 ];
	const float deltaA = 2.0f * 2.0f / (float)(size * size);

	for(int i=0; i<6; i++) {
		if(deltaFormFactor[i] != NULL) delete [] deltaFormFactor[i];
		deltaFormFactor[i] = new float[size*size];
	}

	for ( int t = 0; t <size; t++ ) {
		for ( int s = 0; s < size; s++ ) {
		
			// +X•ûŒü‚Ì–Ê
			getPositiveXPoint(s, t, size, vec);

			float d = ( vec[ 0 ] * vec[ 0 ] + vec[ 1 ] * vec[ 1 ] + vec[ 2 ] * vec[ 2 ] );
			d = 1.0f / ( d * sqrt( d ) );
			deltaFormFactor[ 0 ][ t * size + s ] = d * deltaA;

			// -X•ûŒü‚Ì–Ê
			getNegativeXPoint(s, t, size, vec);

			d = ( vec[ 0 ] * vec[ 0 ] + vec[ 1 ] * vec[ 1 ] + vec[ 2 ] * vec[ 2 ] );
			d = 1.0f / ( d * sqrt( d ) );
			deltaFormFactor[ 1 ][ t * size + s ] = d * deltaA;

			// +Y•ûŒü‚Ì–Ê
			getPositiveYPoint(s, t, size, vec);

			d = ( vec[ 0 ] * vec[ 0 ] + vec[ 1 ] * vec[ 1 ] + vec[ 2 ] * vec[ 2 ] );
			d = 1.0f / ( d * sqrt( d ) );
			deltaFormFactor[ 2 ][ t * size + s ] = d * deltaA;

			// -Y•ûŒü‚Ì–Ê
			getNegativeYPoint(s, t, size, vec);

			d = ( vec[ 0 ] * vec[ 0 ] + vec[ 1 ] * vec[ 1 ] + vec[ 2 ] * vec[ 2 ] );
			d = 1.0f / ( d * sqrt( d ) );
			deltaFormFactor[ 3 ][ t * size + s ] = d * deltaA;

			// +Z•ûŒü‚Ì–Ê
			getPositiveZPoint(s, t, size, vec);

			d = ( vec[ 0 ] * vec[ 0 ] + vec[ 1 ] * vec[ 1 ] + vec[ 2 ] * vec[ 2 ] );
			d = 1.0f / ( d * sqrt( d ) );
			deltaFormFactor[ 4 ][ t * size + s ] = d * deltaA;

			// -Z•ûŒü‚Ì–Ê
			getNegativeZPoint(s, t, size, vec);

			d = ( vec[ 0 ] * vec[ 0 ] + vec[ 1 ] * vec[ 1 ] + vec[ 2 ] * vec[ 2 ] );
			d = 1.0f / ( d * sqrt( d ) );	
			deltaFormFactor[ 5 ][ t * size + s ] = d * deltaA;
		}
	}
}
//-----------------------------------------------------------------------
void initSHTable(int sh_order, int size)
{
	calDeltaFormFactor(size);

	shOrder = sh_order;
	SphericalHarmonics sh(sh_order);
	int num_basis = sh.getNumBasis();

	for(int i=0; i<6; i++) {
		if(cubeSHTable[i] != NULL) delete [] cubeSHTable[i];
		cubeSHTable[i] = new float[num_basis*size*size];
	}

	// cubeSHTable‚Ì‰Šú‰»
	float vec[ 3 ];
	double	ra, th, ph;
	vector<double> y(num_basis);

	for ( int t = 0; t < size; t++ ) {
		for ( int s = 0; s < size; s++ ) {
		
			// +X•ûŒü‚Ì–Ê
			getPositiveXPoint(s, t, size, vec);
			XYZ2RTP(vec[0], vec[1], vec[2], ra, th, ph);

			sh.calc(y, cos(th), ph);
			for ( int i = 0; i < num_basis; ++i )
			{
				cubeSHTable[ 0 ][IX3(i, s, t, num_basis, size)] = y[ i ]; 
			}

			// -X•ûŒü‚Ì–Ê
			getNegativeXPoint(s, t, size, vec);
			XYZ2RTP(vec[0], vec[1], vec[2], ra, th, ph);

			sh.calc(y, cos(th), ph);
			for ( int i = 0; i < num_basis; ++i )
			{
				cubeSHTable[ 1 ][IX3(i, s, t, num_basis, size)] = y[ i ]; 
			}

			// +Y•ûŒü‚Ì–Ê
			getPositiveYPoint(s, t, size, vec);
			XYZ2RTP(vec[0], vec[1], vec[2], ra, th, ph);

			sh.calc(y, cos(th), ph);
			for ( int i = 0; i < num_basis; ++i )
			{
				cubeSHTable[ 2 ][IX3(i, s, t, num_basis, size)] = y[ i ]; 
			}

			// -Y•ûŒü‚Ì–Ê
			getNegativeYPoint(s, t, size, vec);
			XYZ2RTP(vec[0], vec[1], vec[2], ra, th, ph);

			sh.calc(y, cos(th), ph);
			for ( int i = 0; i < num_basis; ++i )
			{
				cubeSHTable[ 3 ][IX3(i, s, t, num_basis, size)] = y[ i ]; 
			}

			// +Z•ûŒü‚Ì–Ê
			getPositiveZPoint(s, t, size, vec);
//			getNegativeZPoint(s, t, size, vec);
			XYZ2RTP(vec[0], vec[1], vec[2], ra, th, ph);

			sh.calc(y, cos(th), ph);
			for ( int i = 0; i < num_basis; ++i )
			{
				cubeSHTable[ 4 ][IX3(i, s, t, num_basis, size)] = y[ i ]; 
			}

			// -Z•ûŒü‚Ì–Ê
			getNegativeZPoint(s, t, size, vec);
//			getPositiveZPoint(s, t, size, vec);
			XYZ2RTP(vec[0], vec[1], vec[2], ra, th, ph);

			sh.calc(y, cos(th), ph);
			for ( int i = 0; i < num_basis; ++i )
			{
				cubeSHTable[ 5 ][IX3(i, s, t, num_basis, size)] = y[ i ]; 
			}
		}
	}

}
//-----------------------------------------------------------------------
void decomposeSH(vector<Tcolor4<float>>& out, TcubeMapTexturef& cubemap)
{

	if(shOrder == -1) {
		fprintf(stderr, "call initSHTable before decompositionSH!\n");
		return;
	}

	int	num_basis = (shOrder+1)*(shOrder+1);
	int size = cubemap.getSize();
	out.resize(num_basis);

	for ( int l = 0; l < num_basis; ++l )
	{
		double r = 0.0;
		double g = 0.0;
		double b = 0.0;
		double a = 0.0;

		for ( int s = 0; s < 6; ++s )
		{
			int t = 0;
			for( int iy = 0; iy < size; iy++) {
				for( int ix = 0; ix < size; ix++) {
					float *pix = cubemap.getFaceTex(s).getPixelAddress(ix, iy);
					float yd = cubeSHTable[ s ][IX3(l, ix, iy, num_basis, size)] * deltaFormFactor[ s ][ t ];

					r += (double)pix[0] * (double)yd;
					g += (double)pix[1] * (double)yd;
					b += (double)pix[2] * (double)yd;
					a += (double)pix[3] * (double)yd;

					t++;
				}
			}
		}

		out[l].assign((float)r, (float)g, (float)b, (float)a);
	}
}
//-----------------------------------------------------------------------
void reconstructSH(vector<Tcolor4<float>>& sh_coeff, TcubeMapTexturef& cubemap) {
	if(shOrder == -1) {
		fprintf(stderr, "call initSHTable before reconstructionSH!\n");
		return;
	}

	int num_basis = (shOrder+1)*(shOrder+1);
	int size = cubemap.getSize();

	for(int s = 0; s<6; s++) {
		cubemap.getFaceTex(s).clear();
		for(int j=0; j<size; j++) {
			for(int i=0; i<size; i++) {
				Tcolor4<float> c(0., 0., 0., 0.);
				for(int l=0; l<num_basis; l++) {
					c += sh_coeff[l] * cubeSHTable[s][IX3(l, i, j, num_basis, size)];
				}
				cubemap.getFaceTex(s).set(i, j, c);

			}
		}
	}
}
//-----------------------------------------------------------------------
