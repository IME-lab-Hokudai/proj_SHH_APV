#include "sphericalHarmonics.h"
#include "CommonDefines.h"
#include "envMap_SH.h"

float*  dOmega;
float*  envSHTable;
int shOrder = -1;
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image.h"
#include "stb_image_write.h"
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

//  precomputed pixel solid angle weight sin(theta)*dtheta*dphi
void calcDeltaFormFactorEquirect(int width, int height)
{
    dOmega = new float[width * height];
    double dTheta = M_PI / height;
    double dPhi = 2.0 * M_PI / width;

    for (int y = 0; y < height; ++y)
    {
        // theta from 0 to pi
        double theta = M_PI * (y + 0.5) / height;
        double sinTheta = sin(theta);

        for (int x = 0; x < width; ++x)
        {
            dOmega[TWO_D_TO_ONE_D(x,y,width)] = (float)sinTheta * dTheta * dPhi;
        }
    }

    //double sum = 0.0;
    //for (int y = 0; y < height; ++y)
    //    for (int x = 0; x < width; ++x)
    //        sum += dOmega[y * width + x];

    //std::cout << "Total Omega = " << sum << std::endl; // Should be ≈ 4π
}

//this function precalculate value of all SH basis in all direction and store in SHTableLookup 
//then when we precompute coefficients we look up from SHTableLookup to calculate Coeff = env_map dot Y
void initSHTable(int sh_order, int width, int height)
{
    calcDeltaFormFactorEquirect(width, height);
    shOrder = sh_order;
    SphericalHarmonics sh(sh_order);
    int num_basis = sh.getNumBasis();
    // Allocate SH basis table for each pixel (flattened 3D array: [basis][width][height])
    envSHTable = new float[num_basis * width * height];

    // Preallocate vector for SH basis
    vector<double> y(num_basis);
    for (int y_idx = 0; y_idx < height; ++y_idx)
    {
        // θ = latitude angle from 0 (top) to π (bottom)
        double v = (y_idx + 0.5) / height;
        double theta = v * M_PI;

        for (int x_idx = 0; x_idx < width; ++x_idx)
        {
            // φ = longitude angle from 0 to 2π
            double u = (x_idx + 0.5) / width;
            double phi = u * 2.0 * M_PI;

            // Compute SH basis at (θ, φ)
            sh.calcSHBasis(y, cos(theta), phi); // ct = cos(θ)

            //storing all Ylm for this direction (i.e 9 Y if l = 2)
            for (int i = 0; i < num_basis; ++i)
            {
                envSHTable[THREE_D_TO_ONE_D(i, x_idx, y_idx, width , height)] = y[i];
            }
        }
    }
    //float Y00 = envSHTable[THREE_D_TO_ONE_D(0, width / 2, height / 2, width, height)];
    //std::cout << "Y00 = " << Y00 << std::endl;
}

void decomposeSH(std::vector<float4>& out, const Falcor::ref<EnvMap>& envMap)
{
    if (shOrder == -1)
    {
        logError( "call initSHTable before decompositionSHEnvMap!");
        return;
    }

    int num_basis = (shOrder + 1) * (shOrder + 1);
    out.resize(num_basis);
    const Falcor::ref<Texture> envMapTex = envMap->getEnvMap();
    int width = envMapTex->getWidth();
    int height = envMapTex->getHeight();

    //read texture data into an array
    float4* data = new float4[width*height];
    envMapTex->getSubresourceBlob(0, &data[0], sizeof(float4) * width * height);

   // float4* tranposedData = TranposeData(data, width, height);

    for (int l = 0; l < num_basis; ++l)
    {
        double r = 0.0;
        double g = 0.0;
        double b = 0.0;
        double a = 0.0;

        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                //read value from env map
                float4 envMapValue = data[TWO_D_TO_ONE_D(x,y,width)]; 

                // Lookup SH basis value at (l, x, y)
                float yd = envSHTable[THREE_D_TO_ONE_D(l, x, y, width, height)]; 

               //  precomputed pixel solid angle weight sin(theta)*dtheta*dphi
                float delta = dOmega[TWO_D_TO_ONE_D(x, y, width)]; 

                double weight = (double)yd * (double)delta;

                r += (double)envMapValue[0] * weight;
                g += (double)envMapValue[1] * weight;
                b += (double)envMapValue[2] * weight;
                a += (double)envMapValue[3] * weight;
            }
        }
        float4 tmp(r, g, b, a);
        out[l] = tmp;
    }
    delete [] data;
}

void reconstructSH(std::vector<float4>& sh_coeff, const Falcor::ref<EnvMap>& envMap, Falcor::ref<Device> pDevice)
{
    int num_basis = (shOrder + 1) * (shOrder + 1);
    const Falcor::ref<Texture> envMapTex = envMap->getEnvMap();
    int width = envMapTex->getWidth();
    int height = envMapTex->getHeight();

  // float3* reconstructedData = new float3[width * height];

    //for (int y = 0; y < height; ++y)
    //{
    //    for (int x = 0; x < width; ++x)
    //    {
    //        float4 tmp = float4 (0,0,0,0);
    //        for (int i = 0; i < num_basis; ++i)
    //        {
    //            tmp += sh_coeff[i] * envSHTable[THREE_D_TO_ONE_D(i, x, y, width, height)]; 
    //        }
    //        float3 color = float3(tmp[0], tmp[1], tmp[2]);
    //        reconstructedData[TWO_D_TO_ONE_D(x, y, width)] = color;
    //    }
    //}

   float* reconstructedData = new float[width * height*3];

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            float4 tmp = float4 (0,0,0,0);
            for (int i = 0; i < num_basis; ++i)
            {
                tmp += sh_coeff[i] * envSHTable[THREE_D_TO_ONE_D(i, x, y, width, height)]; 
            }
            float3 color = float3(tmp[0], tmp[1], tmp[2]);
            reconstructedData[TWO_D_TO_ONE_D(x, y, width)*3 + 0] = color[0];
            reconstructedData[TWO_D_TO_ONE_D(x, y, width)*3 + 1] = color[1];
            reconstructedData[TWO_D_TO_ONE_D(x, y, width)*3 + 2] = color[2];
        }
    }

    //Falcor::ref<Texture> outTex =
    //    pDevice->createTexture2D(width, height, ResourceFormat::RGB32Float, 1, 1, reconstructedData, ResourceBindFlags::ShaderResource);

    //outTex->captureToFile(0, 0, "reconstructed.pfm", Bitmap::FileFormat::PfmFile, Bitmap::ExportFlags::None, false);
    //cout << "First 3" << endl;
    //for (int i = 0; i < 9 ; i+=3)
    //{
    //    cout << reconstructedData[i] << " " << reconstructedData[i + 1] << " " << reconstructedData[i + 2] << endl;
    //}
    //cout << "end first 3" << endl;

    stbi_write_hdr("reconstructed_env.hdr", width, height, 3, reconstructedData);

    delete[] reconstructedData;
}

float4* TranposeData(float4* data, int width, int height)
{
    // tranpose data because input hdr file is column major
    float4* tranposedData = new float4[width * height];

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int srcIdx = (x * height + y); // column-major
            int dstIdx = (y * width + x);  // row-major

            tranposedData[dstIdx + 0] = data[srcIdx + 0];
        }
    }
    return tranposedData;
}


