#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void writePPM(string filename, int xRes, int yRes, float* values) {
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

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  if (argc != 5) {
    cout << "Usage: ./run x1 y1 x2 y2" << endl;
    return 1;
  }

  // Print out the command line args
  cout << " Command line arguments were: " << endl;
  for (int x = 0; x < argc; x++)
    cout << argv[x] << " ";
  cout << endl;

  // try to convert them into ints
  vector<int> commandInts(argc);
  for (int x = 0; x < argc; x++)
    commandInts[x] = atoi(argv[x]);

  cout << " I tried to convert each into an integer and got: " << endl;
  for (int x = 0; x < argc; x++)
    cout << commandInts[x] << " ";
  cout << endl;

  int size = 500 * 500 * 3;
  float* values = NULL;

  values = new float[size];
  for (int i = 0; i < size; i++)
    values[i] = 255.0;

  int i0 = commandInts[1];
  int j0 = commandInts[2];
  int i1 = commandInts[3];
  int j1 = commandInts[4];

  int dY = j1 - j0;
  int dX = i1 - i0;
  int d = dX - 2 * dY;

  while (i0 < i1) {
    if (d >= 0) {
      i0 += 1;
      j0 += 1;
      d = d + (dX - dY);
    }
    else {
      i0 += 1;
      d -= dY;
    }
    int index = 3 * (j0 * 500) + i0 * 3;
    values[index  ] = 0.0;
    values[index+1] = 0.0;
    values[index+2] = 0.0;
  }

  // THE BELOW HORIZONTALLY FLIPS THE IMAGE
  // so that the origin is that of a cartesian plane
  //
  // int w = 500;
  //
  // float temp;
  // for (int i = 0; i < size/2; i++) {
  //   temp = values[i];
  //   values[i] = values[(w * 499) + (i % w) - int(floor(i / w)) * w];
  //   values[(w * 499) + (i % w) - int(floor(i / w)) * w] = temp;
  // }

  writePPM("line.ppm", 500, 500, values);
}
