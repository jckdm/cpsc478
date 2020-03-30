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
#define EPSILON 0.001

bool triangle = false;

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

// clamp, but works with clang
VEC3 fix(VEC3 c) {
  for (int i = 0; i < c.size(); i++) {
    float x = c[i] * 255.0;
    c[i] = (x < 0.0) ? 0.0 : x;
    c[i] = (x > 255.0) ? 255.0 : x;
  }
  return c;
}

// 2 - 6.ppm
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

// 7 - 8.ppm
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

// 10 - 11.ppm
void setUp3(vector<VEC3>& cen3, vector<VEC3>& col, vector<VEC3>& tri) {
  cen3.push_back(VEC3(-3.5, 0.0, 10.0));

  cen3.push_back(VEC3(0.0, -1000.0, 10.0));

  for (int i = -20; i <= 18; i += 2) {
    for (int j = -2; j <= 16; j += 2) {
      cen3.push_back(VEC3(i, j, 20.0));
    }
  }

  col.push_back(VEC3(1.0, 0.25, 0.25));
  col.push_back(VEC3(0.5, 0.5, 0.5));

  tri.push_back(VEC3(0.5, -3.0, 10.0));
  tri.push_back(VEC3(0.5, 3.0, 10.0));
  tri.push_back(VEC3(6.5, 3.0, 10.0));

  tri.push_back(VEC3(0.5, -3.0, 10.0));
  tri.push_back(VEC3(6.5, 3.0, 10.0));
  tri.push_back(VEC3(6.5, -3.0, 10.0));

  float v = sin(45.0 * (PI_ / 180.0));

  // rotate by 45 degrees
  MATRIX3 Rot = Rot.setZero();
  Rot(0, 0) = Rot(0, 2) = Rot(2, 2) = v;
  Rot(1, 1) = 1.0;
  Rot(2, 0) = -v;

  VEC3 center(3.5, 0.0, 10.0);

  for (int i = 0; i < tri.size(); i++) {
    tri[i] -= center;
    tri[i] = Rot * tri[i];
    tri[i] += center;
  }
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

// Möller–Trumbore intersection algorithm, via Wikipedia
VEC3 intersectTriangle(VEC3 v0, VEC3 v1, VEC3 v2, VEC3 dir, VEC3 o, double &t) {
  VEC3 bad(FLT_MIN, FLT_MIN, 0.0);

  VEC3 e1, e2, h, s, q;
  float a, f, u, v;

  e1 = v1 - v0;
  e2 = v2 - v0;

  h = dir.cross(e2);
  a = e1.dot(h);
  if (a > -EPSILON && a < EPSILON) { return bad; }
  f = 1.0 / a;
  s = o - v0;
  u = f * s.dot(h);
  if (u < 0.0 || (u + v) > 1.0) { return bad; }
  q = s.cross(e1);
  v = f * dir.dot(q);
  if (v < 0.0 || (u + v) > 1.0) { return bad; }
  t = f * e2.dot(q);

  if (t > EPSILON) { return o + (t * dir); }
  else { return bad; }
}

// calculates the Normal of a triangle
VEC3 triN(VEC3 v0, VEC3 v1, VEC3 v2) { return ((v1 - v0).cross((v2 - v0))).normalized(); }

VEC3 intersectSphere(VEC3 dir, VEC3 o, double &t, VEC3 cen, float r, bool shadows, bool& inout) {
  VEC3 bad(FLT_MIN, FLT_MIN, 0.0);

  VEC3 l = o - cen;
  float a = dir.dot(dir);
  float b = 2 * dir.dot(l);
  float c = l.dot(l) - pow(r, 2.0);

  float discr = (b * b) - 4 * a * c;
  if (discr < 0.0) { return bad; }

  double root1 = (-b - sqrt(discr)) / (2.0 * a);
  double root2 = (-b + sqrt(discr)) / (2.0 * a);

  float epsi = (shadows == false) ? 0.0 : EPSILON;

  if (root1 > epsi) { inout = true; t = root1; return o + (root1 * dir); }
  else if (root2 > epsi) { inout = false; t = root2; return o + (root2 * dir); }
  return bad;
}

VEC3 intersectScene(VEC3 ray, VEC3 o, vector<VEC3>& cen, float *r, int &s, bool shadows, bool& inout, vector<VEC3>& tri) {
  VEC3 closest(FLT_MIN, FLT_MIN, 0.0);
  VEC3 bad(FLT_MIN, FLT_MIN, 0.0);
  double t, close = FLT_MAX;

  for (int i = 0; i < cen.size(); i++) {
    int y = (triangle == true) ? 2 : 3;
    float rad = (i >= y) ? 1.0 : r[i];
    VEC3 p = intersectSphere(ray, o, t, cen[i], rad, shadows, inout);
    if (p != bad) {
      if (t < close) {
        close = t;
        closest = p;
        s = i;
      }
    }
  }
  if (triangle == true) {
    for (int i = 0; i < 2; i++) {
      VEC3 v0 = tri[(3*i)    ];
      VEC3 v1 = tri[(3*i) + 1];
      VEC3 v2 = tri[(3*i) + 2];
      VEC3 tp = intersectTriangle(v0, v1, v2, ray, o, t);
      if (tp != bad) {
        if (t < close) {
          close = t;
          closest = tp;
          s = i - 2;
          return closest;
        }
      }
    }
  }
  return closest;
}

VEC3 refraction(VEC3 dir, VEC3 n, bool& inout) {
  float In, Out;
  if (inout) { In = 1.0; Out = 1.5; }
  else { In = 1.5; Out = 1.0; }

  if (n.dot(dir) > 0.0) { n = -n; }

  float coss = n.dot(-dir);
  float sinn = sqrt(1.0 - pow(coss, 2.0));
  VEC3 b = ((dir + (coss * dir)) / sinn).normalized();
  float cosPP = 1.0 - pow(In/Out, 2.0) * (1 - pow(n.dot(-dir), 2.0));

  if (cosPP < EPSILON) { return VEC3(FLT_MIN, FLT_MIN, 0.0); }
  return (((In / Out) * sinn) * b - sqrt(cosPP) * n).normalized();
}

VEC3 reflection(VEC3 dir, VEC3 n) { return (dir - 2 * dir.dot(n) * n).normalized(); }

VEC3 rayColor(VEC3 dir, VEC3 o, vector<VEC3>& tri, vector<VEC3>& col, vector<VEC3>& cen, int &s, float *r, vector<VEC3>& lpos, vector<VEC3>& lcol, bool flag, bool phong, bool shadows, int Lrec, int Rrec, bool& inout, bool triref) {
  VEC3 color(0.0, 0.0, 0.0);
  VEC3 p = intersectScene(dir, o, cen, r, s, shadows, inout, tri);
  if (s == 203) { return color; }

  VEC3 n, Kd, l, Ii, e1, e2;
  float nl;
  VEC3 c = color;

  if (triangle == true) {
    if (s == -2) {
      Kd = VEC3(0.0, 1.0, 0.0);
      n = triN(tri[0], tri[1], tri[2]);
    }
    else if (s == -1) {
      Kd = VEC3(1.0, 0.0, 0.0);
      n = triN(tri[3], tri[4], tri[5]);
    }
    else {
      n = (p - cen[s]).normalized();
      if (s > 1) { Kd = VEC3(1.0, 1.0, 1.0); }
      else if (s >= 0) { Kd = col[s]; }
    }
  }
  if (triangle == false) {
    n = (p - cen[s]).normalized();
    if (s > 2) { Kd = VEC3(1.0, 1.0, 1.0); }
    else { Kd = col[s]; }
  }

  VEC3 bias = (0.001 * n);

  if (flag == false) {
    l = (lpos[0] - p).normalized();
    nl = max(0.0, n.dot(l));
    return Kd.cwiseProduct(lcol[0]) * nl;
  }
  if (s == 0 && Lrec != 0 && triref == false) {
    col[0] = VEC3(0.0, 0.0, 0.0);
    VEC3 ref = reflection(dir, n);
    VEC3 orig;
    if (inout) { orig = p + bias; }
    else { orig = p - bias; }
    color += rayColor(ref, orig, tri, col, cen, s, r, lpos, lcol, true, true, true, Lrec-1, Rrec, inout, triref);
  }
  if (triref == true && Lrec != 0 && ((s == -1) || (s == -2))) {
    Kd = VEC3(0.0, 0.0, 0.0);
    VEC3 ref = reflection(dir, n);
    VEC3 orig;
    if (inout) { orig = p + bias; }
    else { orig = p - bias; }
    color += rayColor(ref, orig, tri, col, cen, s, r, lpos, lcol, true, true, true, Lrec-1, Rrec, inout, triref);
  }
  if (s == 0 && Rrec != 0) {
    col[0] = VEC3(0.0, 0.0, 0.0);
    VEC3 refrac = refraction(dir, n, inout);
    if (refrac != VEC3(FLT_MIN, FLT_MIN, 0.0)) {
      VEC3 orig;
      if (inout) { orig = p + bias; }
      else { orig = p - bias; }
      color += rayColor(refrac, orig, tri, col, cen, s, r, lpos, lcol, true, true, true, Lrec, Rrec-1, inout, triref);
    }
  }
  for (int i = 0; i < lpos.size(); i++) {
    l = (lpos[i] - p).normalized();
    nl = max(0.0, n.dot(l));
    Ii = lcol[i];
    if (phong == false) { color += Kd.cwiseProduct(Ii) * nl; }
    if (phong == true) {
      VEC3 rr = (-l + 2 * (l.dot(n)) * n).normalized();
      float nh = pow(max(0.0, rr.dot(-p.normalized())), 10.0);
      c = (Kd.cwiseProduct(Ii) * nl) + (Kd.cwiseProduct(Ii) * nh);
    }
    if (shadows == true) {
      int ss = 3;
      intersectScene(l, p, cen, r, ss, true, inout, tri);
      if (ss == 3) { color += c; }
    }
    else { color += c; }
  }
  return color;
}

void make1(float *values, bool xy, bool ab) {
  for (int y = 0; y < YRES; y++) {
    for (int x = 0; x < XRES; x++) {

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

void make2(float *values, vector<VEC3>& cen, vector<VEC3>& col, float *r, vector<VEC3>& tri) {
  for (int y = 0; y < YRES; y++) {
    for (int x = 0; x < XRES; x++) {

      int s = 3; bool inout;
      intersectScene(genRay(x, (YRES - y)), VEC3(0.0, 0.0, 0.0), cen, r, s, false, inout, tri);
      VEC3 c = fix(col[s]);
      int in = 3 * (y * XRES + x);

      values[in    ] = c[0];
      values[in + 1] = c[1];
      values[in + 2] = c[2];
    }
  }
}

void makerest(float *values, vector<VEC3>& cen, vector<VEC3>& col, float *r, vector<VEC3>& lpos, vector<VEC3>& lcol, bool flag, bool phong, bool shadows, int Lrec, int Rrec, vector<VEC3>& tri, bool triref) {
  if (tri.size() > 0) { triangle = true; }

  for (int y = 0; y < YRES; y++) {
    for (int x = 0; x < XRES; x++) {

      int s = 203; bool inout;
      VEC3 c = fix(rayColor(genRay(x, (YRES - y)), VEC3(0.0, 0.0, 0.0), tri, col, cen, s, r, lpos, lcol, flag, phong, shadows, Lrec, Rrec, inout, triref));
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
  vector<VEC3> tri;
  setUp(cen, col, lpos, lcol);

  make2(values, cen, col, r, tri);
  writePPM("2.ppm", xRes, yRes, values);

  makerest(values, cen, col, r, lpos, lcol, false, false, false, 0, 0, tri, false);
  writePPM("3.ppm", xRes, yRes, values);

  makerest(values, cen, col, r, lpos, lcol, true, false, false, 0, 0, tri, false);
  writePPM("4.ppm", xRes, yRes, values);

  makerest(values, cen, col, r, lpos, lcol, true, true, false, 0, 0, tri, false);
  writePPM("5.ppm", xRes, yRes, values);

  makerest(values, cen, col, r, lpos, lcol, true, true, true, 0, 0, tri, false);
  writePPM("6.ppm", xRes, yRes, values);

  cen.clear(); lpos.clear(); lcol.clear();
  setUp2(cen, lpos, lcol);

  makerest(values, cen, col, r, lpos, lcol, true, true, true, 1, 0, tri, false);
  writePPM("7.ppm", xRes, yRes, values);

  makerest(values, cen, col, r, lpos, lcol, true, true, true, 0, 5, tri, false);
  writePPM("8.ppm", xRes, yRes, values);

  // This comment in memoriam of the Fresnel Effects that could've been, but never were.

  cen.clear(); col.clear();
  float r2[2] = {3.0, 997.0};
  setUp3(cen, col, tri);

  makerest(values, cen, col, r2, lpos, lcol, true, true, true, 0, 5, tri, false);
  writePPM("10.ppm", xRes, yRes, values);

  makerest(values, cen, col, r2, lpos, lcol, true, true, true, 1, 5, tri, true);
  writePPM("11.ppm", xRes, yRes, values);

  return 0;
}
