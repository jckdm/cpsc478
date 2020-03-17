#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <cfloat>

#include "SETTINGS.h"

#define PI_ 3.14159265358979323846
#define XRES 800
#define YRES 600
#define SIZE 480000

using namespace std;

// clear buffers
void clean(float *a) { for (int i = 0; i < SIZE*3; i++) { a[i] = 0.0; } }

VEC3 truncate(const VEC4& v) { return VEC3(v[0], v[1], v[2]); }
VEC4 extend(const VEC3& v) { return VEC4(v[0], v[1], v[2], 4.0); }

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

// VEC3 rayColor(VEC3 ray) {
//   VEC3 color(0.0, 0.0, 0.0);
//   if (intersectScene(ray)) {
//     // for each light {
//     //   color += 3 term lighting
//     // }
//   }
//   return color;
// }

VEC3 fix(VEC3 c) {
  for (int i = 0; i < c.size(); i++) {
    float t = c[i] * 255.0;
    if (t < 0.0) { t = 0.0; }
    if (t > 255.0) { t = 255.0; }
    c[i] = t;
  }
  return c;
}

void setUp(vector<VEC3>& cen, vector<VEC3>& col, vector<VEC3>& lpos, vector<VEC3>& lcol) {
  cen.push_back(VEC3(0.0, -1000.0, 10.0));
  cen.push_back(VEC3(-3.5, 0.0, 10.0));
  cen.push_back(VEC3(3.5, 0.0, 10.0));

  col.push_back(VEC3(0.5, 0.5, 0.5));
  col.push_back(VEC3(1.0, 0.25, 0.25));
  col.push_back(VEC3(0.25, 0.25, 1.0));

  lpos.push_back(VEC3(10.0, 3.0, 5.0));
  lpos.push_back(VEC3(-10.0, 3.0, 7.5));

  lcol.push_back(VEC3(1.0, 1.0, 1.0));
  lcol.push_back(VEC3(0.5, 0.0, 0.0));
}

VEC3 intersectSphere(VEC3 ray, VEC3 cen, float r) {
  VEC3 bad(FLT_MAX, FLT_MAX, 0.0);
  VEC3 e(0.0, 0.0, 0.0); // eye
  VEC3 l = e - cen;
  float a = ray.dot(ray);
  float b = 2 * ray.dot(l);
  float c = l.dot(l) - pow(r, 2.0);

  float discr = (b * b) - 4 * a * c;
  if (discr < 0.0) { return bad; }

  float root1 = (-1.0 * b - sqrt(discr)) / (2.0 * a);
  float root2 = (-1.0 * b + sqrt(discr)) / (2.0 * a);

  if (root1 > 0.0) { return root1 * ray; }
  else if (root2 > 0.0) { return root2 * ray; }
  return bad;
}

VEC4 intersectScene(VEC3 ray, vector<VEC3>& cen, vector<VEC3>& col, VEC3 r) {
  VEC4 closest = VEC4(FLT_MAX, FLT_MAX, 3.0, 3.0);

  for (int i = 0; i < cen.size(); i++) {
    VEC3 n = intersectSphere(ray, cen[i], r[i]);
    if (n != VEC3(FLT_MAX, FLT_MAX, 0.0)) {
      closest = extend(n);
      closest[3] = i;
    }
  }
  return closest;
}

VEC3 genRay(int x, int y) {
  VEC3 e(0.0, 0.0, 0.0);
  VEC3 l(0.0, 0.0, 1.0);
  VEC3 t(0.0, 1.0, 0.0);

  VEC3 ww = (e - l).normalized();
  VEC3 uu = t.cross(ww).normalized();
  VEC3 vv = ww.cross(uu);

  float h = tan(PI_ * 32.5 / 180) * 2.0;
  float w = (h * 800.0) / 600.0;

  float right = w/2.0;
  float left = -1.0 * right;
  float top = h/2.0;
  float bottom = -1.0 * top;

  float u = left + (right - left) * (x + 0.5) / XRES;
  float v = bottom + (top - bottom) * (y + 0.5) / YRES;

  VEC3 dir = uu*u*-1.0 + vv*v + ww*-1.0;
  dir.normalize();
  return dir;
}

VEC3 rayColor(VEC3 ray, vector<VEC3>& col, int s, vector<VEC3>& lpos, vector<VEC3>& lcol, VEC3 n, bool flag, bool phong) {
  VEC3 color(0.0, 0.0, 0.0);
  n.normalize();

  // 3.ppm
  if (flag == false) {
    VEC3 l = lpos[0] - n;
    l.normalize();
    float m = max(0.0, n.dot(l));
    for (int i = 0; i < 3; i++) { color[i] = col[s][i] * lcol[0][i] * m; }
  }
  else {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 2; j++) {
        VEC3 l = lpos[j] - n;
        l.normalize();

        float m = max(0.0, n.dot(l));
        // 4.ppm
        if (phong == false) {
          color[i] += col[s][i] * lcol[j][i] * m;
        }
        // 5.ppm
        else {
          VEC3 h = l + VEC3(0.0, 0.0, 1.0);
          h.normalize();
          color[i] += (col[s][i] * lcol[j][i] * m) + (col[s][i] * lcol[j][i] * pow(max(0.0, n.dot(h)), 10.0));
        }
      }
    }
  }
  return color;
}

void make1(float *values, bool xy, bool ab) {
  for (int y = 0; y < YRES; y++) {
    for (int x = 0; x < XRES; x++) {

      VEC3 ray = genRay(x, y);
      int in = 3 * (y * XRES + x);

      if (xy == true) {
        if (ab == true) { ray[0] = abs(ray[0]); }
        if (ray[0] < 0.0) { ray[0] = 0.0; }
        values[in    ] = ray[0] * 255.0;
        values[in + 1] = 0.0;
        values[in + 2] = 0.0;
      }
      if (xy == false) {
        if (ab == true) { ray[1] = abs(ray[1]); }
        if (ray[1] < 0.0) { ray[1] = 0.0; }
        values[in    ] = 0.0;
        values[in + 1] = ray[1] * 255.0;
        values[in + 2] = 0.0;
      }
    }
  }
}

//ray-sphere intersection
void make2(float *values, vector<VEC3>& cen, vector<VEC3>& col, VEC3 r) {
  for (int y = 0; y < YRES; y++) {
    for (int x = 0; x < XRES; x++) {

      VEC3 ray = genRay(x, (YRES - y));
      VEC4 result = intersectScene(ray, cen, col, r);
      int s = result[3];
      VEC3 c = fix(col[s]);
      int in = 3 * (y * XRES + x);

      if (s < 3) {
        values[in    ] = c[0];
        values[in + 1] = c[1];
        values[in + 2] = c[2];
      }
    }
  }
}

// difuse shading
void make345(float *values, vector<VEC3>& cen, vector<VEC3>& col, VEC3 r, vector<VEC3>& lpos, vector<VEC3>& lcol, bool flag, bool phong) {
  for (int y = 0; y < YRES; y++) {
    for (int x = 0; x < XRES; x++) {

      VEC3 ray = genRay(x, (YRES - y));
      VEC4 n = intersectScene(ray, cen, col, r);
      VEC3 c = fix(rayColor(ray, col, n[3], lpos, lcol, truncate(n), flag, phong));
      int in = 3 * (y * XRES + x);

      if (n[3] < 3) {
        values[in    ] = c[0];
        values[in + 1] = c[1];
        values[in + 2] = c[2];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
  int xRes = 800;
  int yRes = 600;
  float values[SIZE * 3];

  VEC3 e(0.0, 0.0, 0.0); // eye
  VEC3 l(0.0, 0.0, 1.0); // lookat
  VEC3 t(0.0, 1.0, 0.0); // up

  VEC3 r(997.0, 3.0, 3.0);

  vector<VEC3> cen;
  vector<VEC3> col;
  vector<VEC3> lpos;
  vector<VEC3> lcol;
  setUp(cen, col, lpos, lcol);

  make1(values, true, false);
  writePPM("1x.ppm", xRes, yRes, values);
  make1(values, true, true);
  writePPM("1xabs.ppm", xRes, yRes, values);
  make1(values, false, false);
  writePPM("1y.ppm", xRes, yRes, values);
  make1(values, false, true);
  writePPM("1yabs.ppm", xRes, yRes, values);

  clean(values);

  make2(values, cen, col, r);
  writePPM("2.ppm", xRes, yRes, values);

  clean(values);

  make345(values, cen, col, r, lcol, lpos, false, false);
  writePPM("3.ppm", xRes, yRes, values);

  clean(values);

  make345(values, cen, col, r, lcol, lpos, true, false);
  writePPM("4.ppm", xRes, yRes, values);

  make345(values, cen, col, r, lcol, lpos, true, true);
  writePPM("5.ppm", xRes, yRes, values);

  return 0;
}
