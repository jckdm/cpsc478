#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <float.h>
#include "SETTINGS.h"
#include "skeleton.h"
#include "displaySkeleton.h"
#include "motion.h"

#define EPSILON 0.001

using namespace std;

// Stick-man classes
DisplaySkeleton displayer;
Skeleton* skeleton;
Motion* motion;

int windowWidth = 640;
int windowHeight = 480;

VEC3 eye(-6, 0.5, 1);
VEC3 lookingAt(5, 0.5, 1);
VEC3 up(0,1,0);

// scene geometry
vector<VEC3> sphereCenters;
vector<float> sphereRadii;
vector<VEC3> sphereColors;

vector<VEC3> triVerts;

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
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

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
bool raySphereIntersect(const VEC3& center,
                        const float radius,
                        const VEC3& rayPos,
                        const VEC3& rayDir,
                        float& t)
{
  const VEC3 op = center - rayPos;
  const float eps = 1e-8;
  const float b = op.dot(rayDir);
  float det = b * b - op.dot(op) + radius * radius;

  // determinant check
  if (det < 0) return false;

  det = sqrt(det);
  t = b - det;
  if (t <= eps) {
    t = b + det;
    if (t <= eps)
      t = -1;
  }
  if (t < 0) return false;
  return true;
}

// Möller–Trumbore intersection algorithm, via Wikipedia
bool rayTriangleIntersect(VEC3 v0, VEC3 v1, VEC3 v2, VEC3 o, VEC3 dir, float &t) {
  VEC3 e1, e2, h, s, q;
  float a, f, u;
  float v = 0.0;

  e1 = v1 - v0;
  e2 = v2 - v0;

  h = dir.cross(e2);
  a = e1.dot(h);
  if (a > -EPSILON && a < EPSILON) { return false; }
  f = 1.0 / a;
  s = o - v0;
  u = f * s.dot(h);
  if (u < 0.0 || (u + v) > 1.0) { return false; }
  q = s.cross(e1);
  v = f * dir.dot(q);
  if (v < 0.0 || (u + v) > 1.0) { return false; }
  t = f * e2.dot(q);

  if (t > EPSILON) { return true; }
  else { return false; }
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void rayColor(const VEC3& rayPos, const VEC3& rayDir, VEC3& pixelColor) {
  pixelColor = VEC3(1,1,1);

  // look for intersections
  int hitID = -1;
  float tMinFound = FLT_MAX;
  for (int y = 0; y < sphereCenters.size(); y++) {
    float tMin = FLT_MAX;
    if (raySphereIntersect(sphereCenters[y], sphereRadii[y], rayPos, rayDir, tMin)) {
      // is the closest so far?
      if (tMin < tMinFound) {
        tMinFound = tMin;
        hitID = y;
      }
    }
  }
  for (int i = 0; i < triVerts.size()/3; i++) {
    float tMin = FLT_MAX;
    VEC3 v0 = triVerts[(3*i)    ];
    VEC3 v1 = triVerts[(3*i) + 1];
    VEC3 v2 = triVerts[(3*i) + 2];

    if (rayTriangleIntersect(v0, v1, v2, rayPos, rayDir, tMin)) {
      if (tMin < tMinFound) {
        tMinFound = tMin;
        if (i % 4 == 0 || i % 4 == 1) { hitID = -2; }
        if (i % 4 == 2 || i % 4 == 3) { hitID = -3; }
      }
    }
  }
  // No intersection, return white
  if (hitID == -1) return;
  // triangles are all blue
  if (hitID == -2) { pixelColor = VEC3(0.35, 0.35, 0.35); }
  else if (hitID == -3) { pixelColor = VEC3(0.5, 0.5, 0.5); }
  else {
    // set to the sphere color
    pixelColor = sphereColors[hitID];
  }
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
float clamp(float value) {
  if (value < 0.0)      return 0.0;
  else if (value > 1.0) return 1.0;
  return value;
}

VEC3 dist(VEC3 one, VEC3 two) {
  return VEC3(two[0] - one[0], two[1] - one[1], two[2] - one[2]);
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void renderImage(int& xRes, int& yRes, const string& filename) {
  // allocate the final image
  const int totalCells = xRes * yRes;
  float* ppmOut = new float[3 * totalCells];

  vector<VEC4>& translations = displayer.translations();
  const VEC4 pT = translations[1];
  VEC3 pelvis(pT[0], pT[1], pT[2]);

  lookingAt -= dist(pelvis, lookingAt);

  // compute image plane
  const float halfY = (lookingAt - eye).norm() * tan(45.0f / 360.0f * M_PI);
  const float halfX = halfY * 4.0f / 3.0f;

  const VEC3 cameraZ = (lookingAt - eye).normalized();
  const VEC3 cameraX = up.cross(cameraZ).normalized();
  const VEC3 cameraY = cameraZ.cross(cameraX).normalized();

  for (int y = 0; y < yRes; y++)
    for (int x = 0; x < xRes; x++) {
      // generate the ray, making x-axis go left to right
      const float ratioX = 1.0f - ((xRes - 1) - x) / float(xRes) * 2.0f;
      const float ratioY = 1.0f - y / float(yRes) * 2.0f;
      const VEC3 rayHitImage = lookingAt +
                               ratioX * halfX * cameraX +
                               ratioY * halfY * cameraY;
      const VEC3 rayDir = (rayHitImage - eye).normalized();

      // get the color
      VEC3 color;
      rayColor(eye, rayDir, color);

      int in = 3 * (y * xRes + x);

      // set, in final image
      ppmOut[in] = clamp(color[0]) * 255.0f;
      ppmOut[in + 1] = clamp(color[1]) * 255.0f;
      ppmOut[in + 2] = clamp(color[2]) * 255.0f;
    }
  writePPM(filename, xRes, yRes, ppmOut);

  delete[] ppmOut;
}

//////////////////////////////////////////////////////////////////////////////////
// Load up a new motion captured frame
//////////////////////////////////////////////////////////////////////////////////
void setSkeletonsToSpecifiedFrame(int frameIndex) {
  if (frameIndex < 0) {
    printf("Error in SetSkeletonsToSpecifiedFrame: frameIndex %d is illegal.\n", frameIndex);
    exit(0);
  }
  if (displayer.GetSkeletonMotion(0) != NULL) {
    int postureID;
    if (frameIndex >= displayer.GetSkeletonMotion(0)->GetNumFrames()) {
      cout << " We hit the last frame! You might want to pick a different sequence. " << endl;
      postureID = displayer.GetSkeletonMotion(0)->GetNumFrames() - 1;
    }
    else
    postureID = frameIndex;
    displayer.GetSkeleton(0)->setPosture(* (displayer.GetSkeletonMotion(0)->GetPosture(postureID)));
  }
}

//////////////////////////////////////////////////////////////////////////////////
// Build a list of spheres in the scene
//////////////////////////////////////////////////////////////////////////////////
void buildScene() {
  sphereCenters.clear();
  sphereRadii.clear();
  sphereColors.clear();
  triVerts.clear();
  displayer.ComputeBonePositions(DisplaySkeleton::BONES_AND_LOCAL_FRAMES);

  // retrieve all the bones of the skeleton
  vector<MATRIX4>& rotations = displayer.rotations();
  vector<MATRIX4>& scalings  = displayer.scalings();
  vector<VEC4>& translations = displayer.translations();
  const VEC4 pelvisTranslation = translations[1];
  vector<float>& lengths     = displayer.lengths();

  // build a sphere list, but skip the first bone,
  // it's just the origin
  int totalBones = rotations.size();
  for (int x = 1; x < totalBones; x++) {
    MATRIX4& rotation = rotations[x];
    MATRIX4& scaling = scalings[x];
    VEC4& translation = translations[x];

    // get the endpoints of the cylinder
    VEC4 leftVertex(0,0,0,1);
    VEC4 rightVertex(0,0,lengths[x],1);

    leftVertex = rotation * scaling * leftVertex + translation;
    rightVertex = rotation * scaling * rightVertex + translation;

    // get the direction vector
    VEC3 direction = (rightVertex - leftVertex).head<3>();
    const float magnitude = direction.norm();
    direction *= 1.0 / magnitude;

    // how many spheres?
    const float sphereRadius = 0.05;
    const int totalSpheres = magnitude / (2.0 * sphereRadius);
    const float rayIncrement = magnitude / (float)totalSpheres;

    // store the spheres
    sphereCenters.push_back(leftVertex.head<3>());
    sphereRadii.push_back(0.05);
    sphereColors.push_back(VEC3(1,0,0));

    sphereCenters.push_back(rightVertex.head<3>());
    sphereRadii.push_back(0.05);
    sphereColors.push_back(VEC3(1,0,0));
    for (int y = 0; y < totalSpheres; y++) {
      VEC3 center = ((float)y + 0.5) * rayIncrement * direction + leftVertex.head<3>();
      sphereCenters.push_back(center);
      sphereRadii.push_back(0.05);
      sphereColors.push_back(VEC3(1,0,0));
    }
  }
  // make the floor
  for (int x = -5; x < 10; x++) {
    for (int z = -5; z < 10; z++) {
      triVerts.push_back(VEC3(x, 0, z));
      triVerts.push_back(VEC3(x, 0, z-1));
      triVerts.push_back(VEC3(x+1, 0, z-1));

      triVerts.push_back(VEC3(x, 0, z));
      triVerts.push_back(VEC3(x+1, 0, z-1));
      triVerts.push_back(VEC3(x+1, 0, z));
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
  string skeletonFilename("14.asf");
  string motionFilename("14_30.amc");

  // load up skeleton stuff
  skeleton = new Skeleton(skeletonFilename.c_str(), MOCAP_SCALE);
  skeleton->setBasePosture();
  displayer.LoadSkeleton(skeleton);

  // load up the motion
  motion = new Motion(motionFilename.c_str(), MOCAP_SCALE, skeleton);
  displayer.LoadMotion(motion);
  skeleton->setPosture(*(displayer.GetSkeletonMotion(0)->GetPosture(0)));

  // Note we're going 8 frames at a time, otherwise the animation
  // is really slow.
  for (int x = 0; x < 2400; x += 8) {
    setSkeletonsToSpecifiedFrame(x);
    buildScene();

    char buffer[256];
    sprintf(buffer, "./frames/frame.%04i.ppm", x / 8);
    renderImage(windowWidth, windowHeight, buffer);
    cout << "Rendered " + to_string(x / 8) + " frames" << endl;
  }
  return 0;
}
