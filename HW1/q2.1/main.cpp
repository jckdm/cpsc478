#include <cmath>
#include <iostream>
#include <vector>
#include <MERSENNE_TWISTER.h>

using namespace std;

// a random number generator
MERSENNE_TWISTER twister(123456);

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void writePPM(string filename, int xRes, int yRes, float* values)
{
  int totalCells = xRes * yRes;
  unsigned char* pixels = new unsigned char[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++) {
    pixels[i] = values[i];
  }

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

int* prob(int p) {

  static int arr[10];

  for (int i = 0; i < 10; i++) {
    if (i < p) { arr[i] = 1; }
    else { arr[i] = 2; }
  }

  return arr;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {

  if (argc != 8) {
    cout << "Usage: ./run probability R1 G1 B1 R2 G2 B2" << endl;
    return 1;
  }

  // Print out the command line args
  cout << " Command line arguments were: " << endl;
  for (int x = 0; x < argc; x++) {
    cout << argv[x] << " ";
  }
    cout << endl;

  // try to convert them into ints
  vector<int> commandInts(argc);
  for (int x = 0; x < argc; x++) {
    int next = atoi(argv[x]);
    if ((next < 0) || (next > 255)) {
      cout << "RGB values must be in range [0, 255]" << endl;
      return 2;
    }
    else { commandInts[x] = next; }
  }

  float probability = stof(argv[1]);

  if ((probability < 0.0) || (probability > 1.0)) {
    cout << "Probability must be in range [0.0, 1.0]" << endl;
    return 3;
  }

  cout << " I tried to convert each into an integer and got: " << endl;
  for (int x = 0; x < argc; x++) {
    cout << commandInts[x] << " ";
  }
    cout << endl;

  int *choices;
  choices = prob(probability * 10);

  float* values = new float[500 * 500 * 3];

  for (int x = 0; x < 500 * 500; x++) {

    int c = int(ceil(twister.rand() * 10)) % 10;

    int choice = choices[c];

    if (choice == 1) {
      values[3 * x    ] = commandInts[2];
      values[3 * x + 1] = commandInts[3];
      values[3 * x + 2] = commandInts[4];
    }
    else {
      values[3 * x    ] = commandInts[5];
      values[3 * x + 1] = commandInts[6];
      values[3 * x + 2] = commandInts[7];
    }
  }

  writePPM("noise.ppm", 500, 500, values);
}
