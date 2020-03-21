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

VEC3 intersectSphere(VEC3 ray, VEC3 cen, float r, bool shadows) {
  VEC3 bad(FLT_MIN, FLT_MIN, 0.0);
  VEC3 e(0.0, 0.0, 0.0); // eye
  VEC3 l = e - cen;
  float a = ray.dot(ray);
  float b = 2 * ray.dot(l);
  float c = l.dot(l) - pow(r, 2.0);

  float discr = (b * b) - 4 * a * c;
  if (discr < 0.0) { return bad; }

  float root1 = (-b - sqrt(discr)) / (2.0 * a);
  float root2 = (-b + sqrt(discr)) / (2.0 * a);

  // point of intersection
  float limit = (shadows == false) ? 0.0 : 0.001;
  if (root1 > limit) { return root1 * ray; }
  else if (root2 > limit) { return root2 * ray; }
  return bad;
}

VEC3 intersectScene(VEC3 ray, vector<VEC3>& cen, vector<VEC3>& col, VEC3 r, int &s, bool flag) {
  VEC3 closest(FLT_MIN, FLT_MIN, 0.0);

  for (int i = 0; i < cen.size(); i++) {
    VEC3 p = intersectSphere(ray, cen[i], r[i], flag);
    if (p != VEC3(FLT_MIN, FLT_MIN, 0.0)) {
      closest = p;
      s = i;
    }
  }
  return closest;
}

VEC3 genRay(int x, int y, vector<VEC3>& cen, vector<VEC3>& col, VEC3 r, int &s, bool flag, bool shadows) {
  VEC3 e(0.0, 0.0, 0.0);
  VEC3 l(0.0, 0.0, 1.0);
  VEC3 t(0.0, 1.0, 0.0);

  VEC3 ww = (e - l).normalized();
  VEC3 uu = t.cross(ww).normalized();
  VEC3 vv = ww.cross(uu);

  float h = tan(PI_ * 32.5 / 180) * 2.0;
  float w = (h * 800.0) / 600.0;

  float right = w/2.0;
  float left = -right;
  float top = h/2.0;
  float bottom = -top;

  float u = left + (right - left) * (x + 0.5) / XRES;
  float v = bottom + (top - bottom) * (y + 0.5) / YRES;

  VEC3 ray = -ww + (-u*uu) + v*vv;
  ray.normalize();

  return (flag == true) ? ray : intersectScene(ray, cen, col, r, s, shadows);
}

VEC3 rayColor(VEC3 p, vector<VEC3>& col, vector<VEC3>& cen, int s, VEC3 r, vector<VEC3>& lpos, vector<VEC3>& lcol, bool flag, bool phong, bool shadows) {
  VEC3 color(0.0, 0.0, 0.0);
  VEC3 v(0.0, 0.0, 1.0);

  VEC3 n = p - cen[s]; n.normalize();
  VEC3 Kd = col[s];

  if (flag == false) {
    VEC3 l = lpos[0] - p; l.normalize();
    float nl = max(0.0, n.dot(l));
    VEC3 Ii = lcol[0];
    color = Kd.cwiseProduct(Ii) * nl;
  }
  else {
    for (int i = 0; i < lpos.size(); i++) {
      VEC3 l = lpos[i] - p; l.normalize();
      float nl = max(0.0, n.dot(l));
      VEC3 Ii = lcol[i];
      if (phong == false) { color += Kd.cwiseProduct(Ii) * nl; }
      if (phong == true) {
        VEC3 h = l - v; h.normalize();
        float nh = pow(max(0.0, n.dot(h)), 10.0);
        color += ((Kd.cwiseProduct(Ii) * nl) + (Kd.cwiseProduct(Ii) * nh));
        if (shadows == true) {
          int ss = 3;
          VEC3 test = intersectScene((p + lpos[i]), cen, col, r, ss, true);
          if (ss >= 3) {
            n = test - cen[ss]; n.normalize();
            Kd = col[ss];
            l = lpos[i] - test; l.normalize();
            nl = max(0.0, n.dot(l));
            h = l - v; h.normalize();
            nh = pow(max(0.0, n.dot(h)), 10.0);
            color += (Kd.cwiseProduct(Ii) * nl) + (Kd.cwiseProduct(Ii) * nh);
          }
          else { color = VEC3(0.0, 0.0, 0.0); }
        }
      }
    }
  }
  return fix(color);
}

void make1(float *values, bool xy, bool ab, vector<VEC3>& cen, vector<VEC3>& col, VEC3 r) {
  for (int y = 0; y < YRES; y++) {
    for (int x = 0; x < XRES; x++) {

      int s;
      VEC3 ray = genRay(x, y, cen, col, r, s, true, false);
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

// ray-sphere intersection
void make2(float *values, vector<VEC3>& cen, vector<VEC3>& col, VEC3 r) {
  for (int y = 0; y < YRES; y++) {
    for (int x = 0; x < XRES; x++) {

      int s = 3;
      genRay(x, (YRES - y), cen, col, r, s, false, false);
      VEC3 c = fix(col[s]);
      int in = 3 * (y * XRES + x);

      values[in    ] = c[0];
      values[in + 1] = c[1];
      values[in + 2] = c[2];
    }
  }
}

// difuse shading
void make3456(float *values, vector<VEC3>& cen, vector<VEC3>& col, VEC3 r, vector<VEC3>& lpos, vector<VEC3>& lcol, bool flag, bool phong, bool shadows) {
  for (int y = 0; y < YRES; y++) {
    for (int x = 0; x < XRES; x++) {

      int s = 3;
      VEC3 ray = genRay(x, (YRES - y), cen, col, r, s, false, false);
      if (s < 3) {
        VEC3 c = rayColor(ray, col, cen, s, r, lpos, lcol, flag, phong, shadows);
        int in = 3 * (y * XRES + x);

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
  VEC3 r(997.0, 3.0, 3.0);
  vector<VEC3> cen;
  vector<VEC3> col;
  vector<VEC3> lpos;
  vector<VEC3> lcol;
  setUp(cen, col, lpos, lcol);

  make1(values, true, false, cen, col, r);
  writePPM("1x.ppm", xRes, yRes, values);
  make1(values, true, true, cen, col, r);
  writePPM("1xabs.ppm", xRes, yRes, values);
  make1(values, false, false, cen, col, r);
  writePPM("1y.ppm", xRes, yRes, values);
  make1(values, false, true, cen, col, r);
  writePPM("1yabs.ppm", xRes, yRes, values);

  make2(values, cen, col, r);
  writePPM("2.ppm", xRes, yRes, values);

  make3456(values, cen, col, r, lpos, lcol, false, false, false);
  writePPM("3.ppm", xRes, yRes, values);

  make3456(values, cen, col, r, lpos, lcol, true, false, false);
  writePPM("4.ppm", xRes, yRes, values);

  make3456(values, cen, col, r, lpos, lcol, true, true, false);
  writePPM("5.ppm", xRes, yRes, values);

  make3456(values, cen, col, r, lpos, lcol, true, true, true);
  writePPM("6.ppm", xRes, yRes, values);

  return 0;
}
