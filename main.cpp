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
void readPPM(const string& filename, int& xRes, int& yRes, float*& values) {
  // try to open the file
  FILE *fp;
  fp = fopen(filename.c_str(), "rb");
  if (fp == NULL) {
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
void writePPM(const string& filename, int& xRes, int& yRes, const float* values) {
  int totalCells = xRes * yRes;
  unsigned char* pixels = new unsigned char[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    pixels[i] = values[i];

  FILE *fp;
  fp = fopen(filename.c_str(), "wb");
  if (fp == NULL) {
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

VEC3 fix(VEC3 c) {
  for (int i = 0; i < c.size(); i++) {
    float x = c[i] * 255.0;
    c[i] = (x < 0.0) ? 0.0 : x;
    c[i] = (x > 255.0) ? 255.0 : x;
  }
  return c;
}

float dist(VEC3 a, VEC3 b) {
  float d = 0.0;
  for (int i = 0; i < 3; i++) { d += pow((b[i] - a[i]), 2.0); }
  return sqrt(d);
}

void setUp(vector<VEC3>& cen, vector<VEC3>& col, vector<VEC3>& lpos, vector<VEC3>& lcol) {
  cen.push_back(VEC3(-3.5, 0.0, 10.0));
  cen.push_back(VEC3(3.5, 0.0, 10.0));
  cen.push_back(VEC3(0.0, -1000.0, 10.0));

  col.push_back(VEC3(1.0, 0.25, 0.25));
  col.push_back(VEC3(0.25, 0.25, 1.0));
  col.push_back(VEC3(0.5, 0.5, 0.5));

  lpos.push_back(VEC3(10.0, 3.0, 5.0));
  lpos.push_back(VEC3(-10.0, 3.0, 7.5));

  lcol.push_back(VEC3(1.0, 1.0, 1.0));
  lcol.push_back(VEC3(0.5, 0.0, 0.0));
}

void setUp2(vector<VEC3>& cen2, vector<VEC3>& lpos2, vector<VEC3>& lcol2) {
  cen2.push_back(VEC3(-3.5, 0.0, 10.0));
  cen2.push_back(VEC3(3.5, 0.0, 10.0));
  cen2.push_back(VEC3(0.0, -1000.0, 10.0));

  for (int i = -20; i <= 18; i += 2) {
    for (int j = -2; j <= 16; j += 2) {
      cen2.push_back(VEC3(i, j, 20.0));
    }
  }

  lpos2.push_back(VEC3(10.0, 10.0, 5.0));
  lpos2.push_back(VEC3(-10.0, 10.0, 7.5));

  lcol2.push_back(VEC3(1.0, 1.0, 1.0));
  lcol2.push_back(VEC3(0.5, 0.25, 0.25));
}

VEC3 genRay(int x, int y) {
  VEC3 e(0.0, 0.0, 0.0);
  VEC3 l(0.0, 0.0, 1.0);
  VEC3 tt(0.0, 1.0, 0.0);

  VEC3 ww = (e - l).normalized();
  VEC3 uu = tt.cross(ww).normalized();
  VEC3 vv = ww.cross(uu);

  float h = tan(PI_ * 32.5 / 180) * 2.0;
  float w = (h * 800.0) / 600.0;

  float right = w/2.0;
  float left = -right;
  float top = h/2.0;
  float bottom = -top;

  float u = left + (right - left) * (x + 0.5) / XRES;
  float v = bottom + (top - bottom) * (y + 0.5) / YRES;

  return (-ww + (-u*uu) + v*vv).normalized();
}

VEC3 intersectSphere(VEC3 dir, VEC3 e, VEC3 cen, float r, bool shadows) {
  VEC3 bad(FLT_MIN, FLT_MIN, 0.0);
  VEC3 l = e - cen;
  float a = dir.dot(dir);
  float b = 2 * dir.dot(l);
  float c = l.dot(l) - pow(r, 2.0);

  float discr = (b * b) - 4 * a * c;
  if (discr < 0.0) { return bad; }

  double root1 = (-b - sqrt(discr)) / (2.0 * a);
  double root2 = (-b + sqrt(discr)) / (2.0 * a);

  float limit = (shadows == false) ? 0.0 : 0.001;
  if (root1 > limit) { return root1 * dir; }
  else if (root2 > limit) { return root2 * dir; }
  return bad;
}

VEC3 intersectScene(VEC3 ray, VEC3 e, vector<VEC3>& cen, float *r, int &s, bool shadows) {
  VEC3 closest(FLT_MIN, FLT_MIN, 0.0);

  for (int i = 0; i < cen.size(); i++) {
    float rad = (i >= 3) ? 1.0 : r[i];
    VEC3 p = intersectSphere(ray, e, cen[i], rad, shadows);
    if (p != VEC3(FLT_MIN, FLT_MIN, 0.0)) {
      if (dist(ray, e) > dist(closest, e)) {
        closest = p;
        s = i;
      }
    }
  }
  return closest;
}

VEC3 rayColor(VEC3 dir, vector<VEC3>& col, vector<VEC3>& cen, int s, float *r, vector<VEC3>& lpos, vector<VEC3>& lcol, int recurse, bool flag, bool phong, bool shadows, bool reflect) {
  VEC3 color(0.0, 0.0, 0.0);
  VEC3 e(0.0, 0.0, 0.0);
  VEC3 w(0.0, 0.0, 1.0);

  VEC3 p = intersectScene(dir, e, cen, r, s, shadows);
  if (s == 203) { return color; }

  VEC3 n = (p - cen[s]).normalized();
  VEC3 Kd;
  if (s >= 3) { Kd = VEC3(1.0, 1.0, 1.0); }
  else { Kd = col[s]; }

  if (reflect == true && s == 0 && recurse < 10) {
    VEC3 ref = (w + 2 * -w.dot(n) * n).normalized();
    return rayColor(ref, col, cen, s, r, lpos, lcol, recurse+1, true, true, true, false);
  }
  if (flag == false) {
    VEC3 l = (lpos[0] - p).normalized();
    float nl = max(0.0, n.dot(l));
    VEC3 Ii = lcol[0];
    color = Kd.cwiseProduct(Ii) * nl;
  }
  else {
    for (int i = 0; i < lpos.size(); i++) {
      VEC3 l = (lpos[i] - p).normalized();
      float nl = max(0.0, n.dot(l));
      VEC3 Ii = lcol[i];
      if (phong == false) { color += Kd.cwiseProduct(Ii) * nl; }
      if (phong == true) {
        VEC3 v = (-p).normalized();
        VEC3 rr = (-l + 2 * (l.dot(n)) * n).normalized();
        float nh = pow(max(0.0, rr.dot(v)), 10.0);
        VEC3 c = (Kd.cwiseProduct(Ii) * nl) + (Kd.cwiseProduct(Ii) * nh);
        if (shadows == true) {
          int ss = 3;
          intersectScene(l, p, cen, r, ss, true);
          if (ss == 3) { color += c; }
        }
        else { color += c; }
      }
    }
  }
  return fix(color);
}

void make1(float *values, bool xy, bool ab) {
  for (int y = 0; y < YRES; y++) {
    for (int x = 0; x < XRES; x++) {

      int s; double t;
      VEC3 dir = genRay(x, (YRES - y));
      int in = 3 * (y * XRES + x);

      if (xy == true) {
        if (ab == true) { dir[0] = abs(dir[0]); }
        if (dir[0] < 0.0) { dir[0] = 0.0; }
        values[in    ] = dir[0] * 255.0;
        values[in + 1] = 0.0;
        values[in + 2] = 0.0;
      }
      if (xy == false) {
        if (ab == true) { dir[1] = abs(dir[1]); }
        if (dir[1] < 0.0) { dir[1] = 0.0; }
        values[in    ] = 0.0;
        values[in + 1] = dir[1] * 255.0;
        values[in + 2] = 0.0;
      }
    }
  }
}
// ray-sphere intersection
void make2(float *values, vector<VEC3>& cen, vector<VEC3>& col, float *r) {
  for (int y = 0; y < YRES; y++) {
    for (int x = 0; x < XRES; x++) {

      int s = 3;
      intersectScene(genRay(x, (YRES - y)), VEC3(0.0, 0.0, 0.0), cen, r, s, false);
      VEC3 c = fix(col[s]);
      int in = 3 * (y * XRES + x);

      values[in    ] = c[0];
      values[in + 1] = c[1];
      values[in + 2] = c[2];
    }
  }
}
// difuse shading, shadows, reflection
void make34567(float *values, vector<VEC3>& cen, vector<VEC3>& col, float *r, vector<VEC3>& lpos, vector<VEC3>& lcol, bool flag, bool phong, bool shadows, bool reflect) {
  for (int y = 0; y < YRES; y++) {
    for (int x = 0; x < XRES; x++) {

      int s = 203;
      VEC3 c = rayColor(genRay(x, (YRES - y)), col, cen, s, r, lpos, lcol, 0, flag, phong, shadows, reflect);
      int in = 3 * (y * XRES + x);

      values[in    ] = c[0];
      values[in + 1] = c[1];
      values[in + 2] = c[2];

    }
  }
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
  int xRes = 800, yRes = 600;
  float values[SIZE * 3];

  make1(values, true, false);
  writePPM("1x.ppm", xRes, yRes, values);
  make1(values, true, true);
  writePPM("1xabs.ppm", xRes, yRes, values);
  make1(values, false, false);
  writePPM("1y.ppm", xRes, yRes, values);
  make1(values, false, true);
  writePPM("1yabs.ppm", xRes, yRes, values);

  float r[3] = {3.0, 3.0, 997.0};
  vector<VEC3> cen;
  vector<VEC3> col;
  vector<VEC3> lpos;
  vector<VEC3> lcol;
  setUp(cen, col, lpos, lcol);

  make2(values, cen, col, r);
  writePPM("2.ppm", xRes, yRes, values);

  make34567(values, cen, col, r, lpos, lcol, false, false, false, false);
  writePPM("3.ppm", xRes, yRes, values);

  make34567(values, cen, col, r, lpos, lcol, true, false, false, false);
  writePPM("4.ppm", xRes, yRes, values);

  make34567(values, cen, col, r, lpos, lcol, true, true, false, false);
  writePPM("5.ppm", xRes, yRes, values);

  make34567(values, cen, col, r, lpos, lcol, true, true, true, false);
  writePPM("6.ppm", xRes, yRes, values);

  cen.clear(); lpos.clear(); lcol.clear();
  setUp2(cen, lpos, lcol);

  make34567(values, cen, col, r, lpos, lcol, true, true, true, true);
  writePPM("7.ppm", xRes, yRes, values);

  return 0;
}
