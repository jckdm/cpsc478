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
#define PI_ 3.14159265358979323846

using namespace std;

// Stick-man classes
DisplaySkeleton displayer;
Skeleton* skeleton;
Motion* motion;

int windowWidth = 640;
int windowHeight = 480;

VEC3 eye(-6, 1.0, 1);
VEC3 lookingAt(5, 0.5, 1);
VEC3 up(0,1,0);

vector<VEC3> cylCenters;
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

VEC3 dist(VEC3 one, VEC3 two) { return VEC3(two[0] - one[0], two[1] - one[1], two[2] - one[2]); }

bool rayCylinderIntersect(const VEC3& rayPos,
                          const VEC3& rayDir,
                          float &t,
                          MATRIX4& rotation,
                          MATRIX4& scaling,
                          VEC4& translation,
                          float& lengths) {
  MATRIX4 Rot = scaling.inverse() * rotation.inverse();

  VEC4 pt = Rot * (VEC4(rayPos[0], rayPos[1], rayPos[2], 0.0) - translation);
  VEC4 Dr = Rot * VEC4(rayDir[0], rayDir[1], rayDir[2], 0.0);

  double a = pow(Dr[0], 2.0) + pow(Dr[1], 2.0);
	double b = (Dr[0] * pt[0]) + (Dr[1] * pt[1]);
  double c = pow(pt[0], 2.0) + pow(pt[1], 2.0) - pow(0.15, 2.0);

	double delta = pow(b, 2.0) - (a * c);
	if (delta < EPSILON) { return false; }
	t = (-b - sqrt(delta)) / a;

  VEC4 ray = pt + (t * Dr);
  if ((ray[2] > lengths) || (ray[2] < 0.0)) { return false; }

	if (t <= EPSILON) { return false; }
	return true;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
bool raySphereIntersect(const VEC3& rayPos,
                        const VEC3& rayDir,
                        float& t) {
  VEC3 center(0, 0.85, 1.5);
  const VEC3 op = center - rayPos;
  const float eps = 1e-8;
  const float b = op.dot(rayDir);
  float det = b * b - op.dot(op) + pow(0.35, 2.0);

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
  float v = 0.0;

  VEC3 e1 = v1 - v0;
  VEC3 e2 = v2 - v0;

  VEC3 h = dir.cross(e2);
  float a = e1.dot(h);
  if (a > -EPSILON && a < EPSILON) { return false; }
  float f = 1.0 / a;
  VEC3 s = o - v0;
  float u = f * s.dot(h);
  if (u < 0.0 || (u + v) > 1.0) { return false; }
  VEC3 q = s.cross(e1);
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

  vector<MATRIX4>& rotations = displayer.rotations();
  vector<MATRIX4>& scalings  = displayer.scalings();
  vector<VEC4>& translations = displayer.translations();
  vector<float>& lengths     = displayer.lengths();

  // look for intersections
  int hitID = -1;
  float tMinFound = FLT_MAX;

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
        if (i >= 450) { hitID = -4; }
        if (i >= 452) { hitID = -5; }
        if (i >= 454) { hitID = -6; }
        if (i >= 456) { hitID = -7; }
        if (i >= 458 || i >= 462) { hitID = -10; }
        if (i >= 464) { hitID = -9; }
      }
    }
  }
  for (int y = 1; y < cylCenters.size(); y++) {
    float tMin = FLT_MAX;
    MATRIX4& rotation = rotations[y];
    MATRIX4& scaling = scalings[y];
    VEC4& translation = translations[y];
    if (rayCylinderIntersect(rayPos, rayDir, tMin, rotation, scaling, translation, lengths[y])) {
      // is the closest so far?
      if (tMin < tMinFound) {
        tMinFound = tMin;
        hitID = y;
      }
    }
    if (raySphereIntersect(rayPos, rayDir, tMin)) {
      if (tMin < tMinFound) {
        tMinFound = tMin;
        hitID = -8;
      }
    }
  }
  // No intersection, return white
  if (hitID == -1) return;
  // triangles are all blue
  if (hitID == -2) { pixelColor = VEC3(0.25, 0.25, 0.25); }
  else if (hitID == -3) { pixelColor = VEC3(0.4, 0.4, 0.4); }
  else if (hitID == -4) { pixelColor = VEC3(0.35, 0.35, 0.5); }
  else if (hitID == -5) { pixelColor = VEC3(0.8, 0.8, 1.0); }
  else if (hitID == -6) { pixelColor = VEC3(0.8, 1.0, 0.8); }
  else if (hitID == -7) { pixelColor = VEC3(1.0, 0.8, 0.8); }
  else if (hitID == -8) { pixelColor = VEC3(0.1, 0.3, 0.7); }
  else if (hitID == -9) { pixelColor = VEC3(0.1, 0.1, 0.1); }
  else if (hitID == -10) { pixelColor = VEC3(0.9, 0.9, 0.9); }
  else { pixelColor = VEC3(0.85,0.65,0.5); }
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
float clamp(float value) {
  if (value < 0.0)      return 0.0;
  else if (value > 1.0) return 1.0;
  return value;
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
  cylCenters.clear();
  triVerts.clear();
  displayer.ComputeBonePositions(DisplaySkeleton::BONES_AND_LOCAL_FRAMES);

  // retrieve all the bones of the skeleton
  vector<MATRIX4>& rotations = displayer.rotations();
  vector<MATRIX4>& scalings  = displayer.scalings();
  vector<VEC4>& translations = displayer.translations();
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
    const int totalSpheres = magnitude / (2.0 * 0.05);
    const float rayIncrement = magnitude / (float)totalSpheres;

    // store the spheres
    cylCenters.push_back(leftVertex.head<3>());
    cylCenters.push_back(rightVertex.head<3>());

    for (int y = 0; y < totalSpheres; y++) {
      VEC3 center = ((float)y + 0.5) * rayIncrement * direction + leftVertex.head<3>();
      cylCenters.push_back(center);
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

  // back wall
  triVerts.push_back(VEC3(10, 0, -5));
  triVerts.push_back(VEC3(10, 0, 10));
  triVerts.push_back(VEC3(10, 5, -5));

  triVerts.push_back(VEC3(10, 5, -5));
  triVerts.push_back(VEC3(10, 5, 10));
  triVerts.push_back(VEC3(10, 0, 10));

  // box side
  triVerts.push_back(VEC3(0, 0.00, -0.32));
  triVerts.push_back(VEC3(1, 0.00, -0.32));
  triVerts.push_back(VEC3(1, 0.45, -0.32));

  triVerts.push_back(VEC3(0, 0.00, -0.32));
  triVerts.push_back(VEC3(0, 0.45, -0.32));
  triVerts.push_back(VEC3(1, 0.45, -0.32));

  // box front
  triVerts.push_back(VEC3(0, 0.00, -0.32));
  triVerts.push_back(VEC3(0, 0.45, -0.32));
  triVerts.push_back(VEC3(0, 0.45, -1.00));

  triVerts.push_back(VEC3(0, 0.00, -0.32));
  triVerts.push_back(VEC3(0, 0.00, -1.00));
  triVerts.push_back(VEC3(0, 0.45, -1.00));

  // box top
  triVerts.push_back(VEC3(0, 0.45, -0.32));
  triVerts.push_back(VEC3(1, 0.45, -0.32));
  triVerts.push_back(VEC3(0, 0.45, -1.00));

  triVerts.push_back(VEC3(0, 0.45, -1.00));
  triVerts.push_back(VEC3(1, 0.45, -0.32));
  triVerts.push_back(VEC3(1, 0.45, -1.00));

  // pedestal side?
  triVerts.push_back(VEC3(-0.5, 0, 2.0));
  triVerts.push_back(VEC3(0.5, 0, 2.0));
  triVerts.push_back(VEC3(0.5, 0.5, 2.0));

  triVerts.push_back(VEC3(-0.5, 0, 2.0));
  triVerts.push_back(VEC3(-0.5, 0.5, 2.0));
  triVerts.push_back(VEC3(0.5, 0.5, 2.0));

  // othr side
  triVerts.push_back(VEC3(-0.5, 0, 1.0));
  triVerts.push_back(VEC3(0.5, 0, 1.0));
  triVerts.push_back(VEC3(0.5, 0.5, 1.0));

  triVerts.push_back(VEC3(-0.5, 0, 1.0));
  triVerts.push_back(VEC3(-0.5, 0.5, 1.0));
  triVerts.push_back(VEC3(0.5, 0.5, 1.0));

  // pedestal front
  triVerts.push_back(VEC3(-0.5, 0, 2));
  triVerts.push_back(VEC3(-0.5, 0.5, 2));
  triVerts.push_back(VEC3(-0.5, 0.5, 1));

  triVerts.push_back(VEC3(-0.5, 0, 2));
  triVerts.push_back(VEC3(-0.5, 0, 1));
  triVerts.push_back(VEC3(-0.5, 0.5, 1));

  // pedestal top
  triVerts.push_back(VEC3(-0.5, 0.5, 2));
  triVerts.push_back(VEC3(0.5, 0.5, 2));
  triVerts.push_back(VEC3(-0.5, 0.5, 1));

  triVerts.push_back(VEC3(-0.5, 0.5, 1));
  triVerts.push_back(VEC3(0.5, 0.5, 2));
  triVerts.push_back(VEC3(0.5, 0.5, 1));
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
  for (int x = 0; x < 8; x += 8) {
    setSkeletonsToSpecifiedFrame(x);
    buildScene();

    char buffer[256];
    sprintf(buffer, "./frames/frame.%04i.ppm", x / 8);
    renderImage(windowWidth, windowHeight, buffer);
    cout << "Rendered " + to_string(x / 8) + " frames" << endl;
  }
  return 0;
}
