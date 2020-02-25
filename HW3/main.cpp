#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <cfloat>

#include "SETTINGS.h"

using namespace std;

#define PI_ 3.14159265358979323846

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
VEC3 truncate(const VEC4& v) { return VEC3(v[0], v[1], v[2]); }
VEC4 extend(const VEC3& v) { return VEC4(v[0], v[1], v[2], 1.0); }
// clear buffers
void clean(int s, float *a, vector<VEC3>& v, vector<VEC3I>& ind, vector<VEC3>& c) {
  for (int i = 0; i < s; i++) { a[i] = 0; }
  v.clear(); ind.clear(); c.clear();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void readPPM(const string& filename, int& xRes, int& yRes, float*& values)
{
  // try to open the file
  FILE *fp;
  fp = fopen(filename.c_str(), "rb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for reading." << endl;
    cout << " Make sure you're not trying to read from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  // get the dimensions
  fscanf(fp, "P6\n%d %d\n255\n", &xRes, &yRes);
  int totalCells = xRes * yRes;

  // grab the pixel values
  unsigned char* pixels = new unsigned char[3 * totalCells];
  fread(pixels, 1, totalCells * 3, fp);

  // copy to a nicer data type
  values = new float[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = pixels[i];

  // clean up
  delete[] pixels;
  fclose(fp);
  cout << " Read in file " << filename.c_str() << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void writePPM(const string& filename, int& xRes, int& yRes, const float* values)
{
  int totalCells = xRes * yRes;
  unsigned char* pixels = new unsigned char[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    pixels[i] = values[i];

  FILE *fp;
  fp = fopen(filename.c_str(), "wb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for writing." << endl;
    cout << " Make sure you're not trying to write from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  fprintf(fp, "P6\n%d %d\n255\n", xRes, yRes);
  fwrite(pixels, 1, totalCells * 3, fp);
  fclose(fp);
  delete[] pixels;
}

//////////////////////////////////////////////////////////////////////////////////
// build out a single square
//////////////////////////////////////////////////////////////////////////////////
void buildBigSquare(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
{
  vertices.push_back(VEC3(1.0,  1.0,  1.5));
  vertices.push_back(VEC3(11.0, 1.0,  1.5));
  vertices.push_back(VEC3(1.0,  11.0, 1.5));
  vertices.push_back(VEC3(11.0, 11.0, 1.5));

  // front face
  indices.push_back(VEC3I(0, 1, 2));
  indices.push_back(VEC3I(2, 1, 3));
  colors.push_back(VEC3(1.0, 0.0, 0.0));
  colors.push_back(VEC3(0.0, 0.0, 1.0));
}

//////////////////////////////////////////////////////////////////////////////////
// build out a single square
//////////////////////////////////////////////////////////////////////////////////
void buildSquare(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
{
  vertices.push_back(VEC3(-0.5, -0.5,  0.5));
  vertices.push_back(VEC3( 0.5, -0.5,  0.5));
  vertices.push_back(VEC3(-0.5,  0.5,  0.5));
  vertices.push_back(VEC3( 0.5,  0.5,  0.5));

  // front face
  // use these indices to loop thru vertices
  indices.push_back(VEC3I(0, 1, 2));
  indices.push_back(VEC3I(2, 1, 3));
  colors.push_back(VEC3(1.0, 0.0, 0.0));
  colors.push_back(VEC3(0.0, 1.0, 0.0));
}

//////////////////////////////////////////////////////////////////////////////////
// build out a cube
//////////////////////////////////////////////////////////////////////////////////

void buildCube(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
{
  vertices.push_back(VEC3(-0.5, -0.5,  0.5));
  vertices.push_back(VEC3( 0.5, -0.5,  0.5));
  vertices.push_back(VEC3(-0.5,  0.5,  0.5));
  vertices.push_back(VEC3( 0.5,  0.5,  0.5));
  vertices.push_back(VEC3(-0.5, -0.5, -0.5));
  vertices.push_back(VEC3( 0.5, -0.5, -0.5));
  vertices.push_back(VEC3(-0.5,  0.5, -0.5));
  vertices.push_back(VEC3( 0.5,  0.5, -0.5));

  // front face
  indices.push_back(VEC3I(0, 1, 2));
  indices.push_back(VEC3I(2, 1, 3));
  colors.push_back(VEC3(1.0, 0.0, 0.0));
  colors.push_back(VEC3(1.0, 0.0, 0.0));

  // back face
  indices.push_back(VEC3I(5, 4, 7));
  indices.push_back(VEC3I(7, 4, 6));
  colors.push_back(VEC3(0.0, 1.0, 0.0));
  colors.push_back(VEC3(0.0, 1.0, 0.0));

  // left face
  indices.push_back(VEC3I(4, 0, 6));
  indices.push_back(VEC3I(6, 0, 2));
  colors.push_back(VEC3(0.0, 0.0, 1.0));
  colors.push_back(VEC3(0.0, 0.0, 1.0));

  // right face
  indices.push_back(VEC3I(1, 5, 3));
  indices.push_back(VEC3I(3, 5, 7));
  colors.push_back(VEC3(0.0, 1.0, 1.0));
  colors.push_back(VEC3(0.0, 1.0, 1.0));

  // top face
  indices.push_back(VEC3I(2, 3, 6));
  indices.push_back(VEC3I(6, 3, 7));
  colors.push_back(VEC3(1.0, 1.0, 0.0));
  colors.push_back(VEC3(1.0, 1.0, 0.0));

  // bottom face
  indices.push_back(VEC3I(4, 5, 0));
  indices.push_back(VEC3I(0, 5, 1));
  colors.push_back(VEC3(1.0, 0.0, 1.0));
  colors.push_back(VEC3(1.0, 0.0, 1.0));
}

//////////////////////////////////////////////////////////////////////////////////
// build out a cube
//////////////////////////////////////////////////////////////////////////////////
void buildCubePerVertexColors(vector<VEC3>& vertices, vector<VEC3I>& indices, vector<VEC3>& colors)
{
  vertices.push_back(VEC3(-0.5, -0.5,  0.5));
  vertices.push_back(VEC3( 0.5, -0.5,  0.5));
  vertices.push_back(VEC3(-0.5,  0.5,  0.5));
  vertices.push_back(VEC3( 0.5,  0.5,  0.5));
  vertices.push_back(VEC3(-0.5, -0.5, -0.5));
  vertices.push_back(VEC3( 0.5, -0.5, -0.5));
  vertices.push_back(VEC3(-0.5,  0.5, -0.5));
  vertices.push_back(VEC3( 0.5,  0.5, -0.5));
  colors.push_back(VEC3(1.0, 0.0, 0.0));
  colors.push_back(VEC3(1.0, 0.0, 0.0));
  colors.push_back(VEC3(0.0, 1.0, 0.0));
  colors.push_back(VEC3(0.0, 1.0, 0.0));
  colors.push_back(VEC3(0.0, 0.0, 1.0));
  colors.push_back(VEC3(0.0, 0.0, 1.0));
  colors.push_back(VEC3(1.0, 1.0, 0.0));
  colors.push_back(VEC3(1.0, 1.0, 0.0));

  // front face
  indices.push_back(VEC3I(0, 1, 2));
  indices.push_back(VEC3I(2, 1, 3));

  // back face
  indices.push_back(VEC3I(5, 4, 7));
  indices.push_back(VEC3I(7, 4, 6));

  // left face
  indices.push_back(VEC3I(4, 0, 6));
  indices.push_back(VEC3I(6, 0, 2));

  // right face
  indices.push_back(VEC3I(1, 5, 3));
  indices.push_back(VEC3I(3, 5, 7));

  // top face
  indices.push_back(VEC3I(2, 3, 6));
  indices.push_back(VEC3I(6, 3, 7));

  // bottom face
  indices.push_back(VEC3I(4, 5, 0));
  indices.push_back(VEC3I(0, 5, 1));
}

// to be used in rasterization
int f01(int x0, int y0, int x1, int y1, int x, int y) {
  return ((y0 - y1) * x) + ((x1 - x0) * y) + (x0 * y1) - (x1 * y0); }

int f12(int x1, int y1, int x2, int y2, int x, int y) {
  return ((y1 - y2) * x) + ((x2 - x1) * y) + (x1 * y2) - (x2 * y1); }

int f20(int x2, int y2, int x0, int y0, int x, int y) {
  return ((y2 - y0) * x) + ((x0 - x2) * y) + (x2 * y0) - (x0 * y2); }

// viewport matrix
void make1(vector<VEC3>& v, vector<VEC3I>& ind, float *values, vector<VEC3>& c, bool flag) {
  int xRes = 800;
  int yRes = 600;

  MATRIX4 Mvp = Mvp.setZero();
  Mvp(0, 0) = xRes/2; Mvp(0, 3) = xRes/2;
  Mvp(1, 1) = yRes/2; Mvp(1, 3) = yRes/2;
  Mvp(2, 2) = 1;
  Mvp(3, 3) = 1;

  for (int i = 0; i < v.size(); i++) {
    VEC4 ve = extend(v[i]);
    v[i] = truncate(Mvp * ve);
  }

  // only do this for 1.ppm
  if (flag == true) {
    for (int i = 0; i < ind.size(); i++) {
      for (int j = 0; j < 3; j++) {
        int index = ind[i][j];
        int x = v[index][0];
        int y = v[index][1];
        int in = 3 * (x + (600 - y) * 800);

        values[in    ] = c[i][0] * 255;
        values[in + 1] = c[i][1] * 255;
        values[in + 2] = c[i][2] * 255;
      }
    }
  }
}

// rasterization
void make2(vector<VEC3>& v, vector<VEC3I>& ind, float *values, vector<VEC3>& c, float *zVal, bool buffer, bool interp) {
  for (int i = 0; i < ind.size(); i++) {
    VEC3I index = ind[i];

    float x0 = v[index[0]][0];
    float y0 = v[index[0]][1];
    float z0 = v[index[0]][2];

    float x1 = v[index[1]][0];
    float y1 = v[index[1]][1];
    float z1 = v[index[1]][2];

    float x2 = v[index[2]][0];
    float y2 = v[index[2]][1];
    float z2 = v[index[2]][2];

    // average depth coords
    float depth = (z0 + z1 + z2) / 3.0;

    // bounds for rasterization
    int xmin = floor(min(min(x0, x1), x2));
    int ymin = floor(min(min(y0, y1), y2));
    int xmax = ceil(max(max(x0, x1), x2));
    int ymax = ceil(max(max(y0, y1), y2));

    int falpha = f12(x1, y1, x2, y2, x0, y0);
    int fbeta  = f20(x2, y2, x0, y0, x1, y1);
    int fgamma = f01(x0, y0, x1, y1, x2, y2);

    for (int y = ymin; y < ymax; y++) {
      for (int x = xmin; x < xmax; x++) {
        float alpha = (float)f12(x1, y1, x2, y2, x, y) / falpha;
        float beta = (float)f20(x2, y2, x0, y0, x, y) / fbeta;
        float gamma = (float)f01(x0, y0, x1, y1, x, y) / fgamma;

        if (alpha >= 0 && beta >= 0 && gamma >= 0) {
          // if ((alpha > 0 || (falpha * f12(x1, y1, x2, y2, -1, -1)) > 0) &&
          //     (beta  > 0 || (fbeta  * f20(x2, y2, x0, y0, -1, -1)) > 0) &&
          //     (gamma > 0 || (fgamma * f01(x0, y0, x1, y1, -1, -1)) > 0)) {
                int in = 3 * (x + (600 - y) * 800);
                int zin = (x + (600 - y) * 800);
                VEC3 color = (c[index[0]] * alpha) + (c[index[1]] * beta) + (c[index[2]] * gamma);

                // z-buffering && color interpolation
                if (buffer == true && interp == true) {
                  if (depth < zVal[zin]) {
                    values[in    ] = color[0] * 255;
                    values[in + 1] = color[1] * 255;
                    values[in + 2] = color[2] * 255;

                    zVal[zin] = depth;
                  }
                } // just z-buffering
                else if (buffer == true) {
                  if (depth < zVal[zin]) {
                    values[in    ] = c[i][0] * 255;
                    values[in + 1] = c[i][1] * 255;
                    values[in + 2] = c[i][2] * 255;

                    zVal[zin] = depth;
                  }
                } // just interpolation
                else if (interp == true) {
                  values[in    ] = color[0] * 255;
                  values[in + 1] = color[1] * 255;
                  values[in + 2] = color[2] * 255;
                }
                else {
                  values[in    ] = c[i][0] * 255;
                  values[in + 1] = c[i][1] * 255;
                  values[in + 2] = c[i][2] * 255;
                }
              }
          //}
        }
      }
    }
}

// orthographic matrix
void make3(vector<VEC3>& v, float l, float r, float b, float t, float f, float n) {
  MATRIX4 Mortho = Mortho.setZero();
  Mortho(0, 0) = 2.0 / (r - l); Mortho(0, 3) = -1.0 * (((2 * l) / (r - l)) + 1);
  Mortho(1, 1) = 2.0 / (t -b); Mortho(1, 3) = -1.0 * (((2 * b) / (t - b)) + 1);
  Mortho(2, 2) = 2.0 / (n - f); Mortho(2, 3) = -1.0 * (((2 * f) / (n - f)) + 1);
  Mortho(3, 3) = 1.0;

  for (int i = 0; i < v.size(); i++) {
    VEC4 ve = extend(v[i]);
    v[i] = truncate(Mortho * ve);
  }
}

// camera matrix
void make4(vector<VEC3>& v, float ex, float ey, float ez, float lx, float ly, float lz, float tx, float ty, float tz) {
  VEC3 e(ex, ey, ez); // eye
  VEC3 l(lx, ly, lz); // lookat
  VEC3 t(tx, ty, tz); // up

  VEC3 w = (e - l).normalized();
  VEC3 u = t.cross(w).normalized();
  VEC3 vv = w.cross(u);

  MATRIX4 Mcb = Mcb.setZero();
  Mcb(0,0) = u[0];  Mcb(0,1) = u[1];  Mcb(0,2) = u[2];
  Mcb(1,0) = vv[0]; Mcb(1,1) = vv[1]; Mcb(1,2) = vv[2];
  Mcb(2,0) = w[0];  Mcb(2,1) = w[1];  Mcb(2,2) = w[2];
  Mcb(3,3) = 1;

  MATRIX4 Mcam = Mcam.setZero();
  Mcam(0, 0) = 1; Mcam(0, 3) = -1.0 * e[0];
  Mcam(1, 1) = 1; Mcam(1, 3) = -1.0 * e[1];
  Mcam(2, 2) = 1; Mcam(2, 3) = -1.0 * e[2];
  Mcam(3, 3) = 1;

  for (int i = 0; i < v.size(); i++) {
    VEC4 ve = extend(v[i]);
    v[i] = truncate(Mcb * Mcam * ve);
  }
}

// perspective projection matrix
void make5(vector<VEC3>& v, float fovy, float aspect, float near, float far, bool flag) {
  float f = 3.0 * (1.0 / tan(fovy/2.0));

  // https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluPerspective.xml
  MATRIX4 Mproject = Mproject.setZero();
  Mproject(0,0) = f / aspect;
  Mproject(1,1) = f;
  Mproject(2,2) = (far + near) / (near - far);
  Mproject(2,3) = (2.0 * far * near) / (near - far);
  Mproject(3,2) = -1.0;

  for (int i = 0; i < v.size(); i++) {
    VEC4 ve = Mproject * extend(v[i]);
    // perspective divide
    if (flag == true) { ve /= ve[3]; }
    v[i] = truncate(ve);
  }
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
  int xRes = 800;
  int yRes = 600;
  int size = 800 * 600 * 3;

  vector<VEC3>  v;
  vector<VEC3I> ind;
  vector<VEC3>  c;

  float values[800 * 600 * 3];
  float zVal[800 * 600];

  if (argc == 1) {
    buildSquare(v, ind, c);
    make1(v, ind, values, c, true);
    writePPM("1.ppm", xRes, yRes, values);

    clean(size, values, v, ind, c);

    buildSquare(v, ind, c);
    make1(v, ind, values, c, false);
    make2(v, ind, values, c, zVal, false, false);
    writePPM("2.ppm", xRes, yRes, values);

    clean(size, values, v, ind, c);

    buildBigSquare(v, ind, c);
    make3(v, 0.0, 12.0, 0.0, 12.0, 0.0, 12.0);
    make1(v, ind, values, c, false);
    make2(v, ind, values, c, zVal, false, false);
    writePPM("3.ppm", xRes, yRes, values);

    clean(size, values, v, ind, c);

    buildBigSquare(v, ind, c);
    make4(v, 0.2, 0.2, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    make3(v, 0.0, 12.0, 0.0, 12.0, 0.0, 12.0);
    make1(v, ind, values, c, false);
    make2(v, ind, values, c, zVal, false, false);
    writePPM("4.ppm", xRes, yRes, values);

    clean(size, values, v, ind, c);

    buildCube(v, ind, c);
    for (int i = 0; i < v.size(); i++) { v[i] *= 0.5; }
    make4(v, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    make5(v, 65.0, (4.0 / 3.0), 1.0, 100.0, false);
    make1(v, ind, values, c, false);
    make2(v, ind, values, c, zVal, false, false);
    writePPM("5.ppm", xRes, yRes, values);

    clean(size, values, v, ind, c);

    buildCube(v, ind, c);
    make4(v, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    make5(v, 65.0, (4.0 / 3.0), 1.0, 100.0, true);
    make1(v, ind, values, c, false);
    make2(v, ind, values, c, zVal, false, false);
    writePPM("6.ppm", xRes, yRes, values);

    clean(size, values, v, ind, c);
    for (int i = 0; i < 800 * 600; i++) { zVal[i] = FLT_MAX; }

    buildCube(v, ind, c);
    make4(v, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    make5(v, 65.0, (4.0 / 3.0), 1.0, 100.0, true);
    make1(v, ind, values, c, false);
    make2(v, ind, values, c, zVal, true, false);
    writePPM("7.ppm", xRes, yRes, values);

    clean(size, values, v, ind, c);
    for (int i = 0; i < 800 * 600; i++) { zVal[i] = FLT_MAX; }

    buildCubePerVertexColors(v, ind, c);
    make4(v, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    make5(v, 65.0, (4.0 / 3.0), 1.0, 100.0, true);
    make1(v, ind, values, c, false);
    make2(v, ind, values, c, zVal, true, true);
    writePPM("8.ppm", xRes, yRes, values);
  }
  else if (argc == 10) {
    float eyex = atof(argv[1]);
    float eyey = atof(argv[2]);
    float eyez = atof(argv[3]);

    float lookx = atof(argv[4]);
    float looky = atof(argv[5]);
    float lookz = atof(argv[6]);

    float upx = atof(argv[7]);
    float upy = atof(argv[8]);
    float upz = atof(argv[9]);

    for (int i = 0; i < 800 * 600; i++) { zVal[i] = FLT_MAX; }

    buildCubePerVertexColors(v, ind, c);
    make4(v, eyex, eyey, eyez, lookx, looky, lookz, upx, upy, upz);
    make5(v, 65.0, (4.0 / 3.0), 1.0, 100.0, true);
    make1(v, ind, values, c, false);
    make2(v, ind, values, c, zVal, true, true);
    writePPM("custom.ppm", xRes, yRes, values);
  }
  else {
    cout << "usage: ./run eye.x eye.y eye.z lookat.x lookat.y lookat.z up.x up.y up.z" << endl;
    return 1;
  }
  return 0;
}
