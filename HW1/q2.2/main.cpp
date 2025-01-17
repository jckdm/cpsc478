#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

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
void writePPM(const string& filename, const int xRes, const int yRes, const float* values)
{
  // copy to the data type PPM expects
  int totalCells = xRes * yRes;
  unsigned char* pixels = new unsigned char[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    pixels[i] = values[i];

  // try to open the file
  FILE *fp;
  fp = fopen(filename.c_str(), "wb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for writing." << endl;
    cout << " Make sure you're not trying to write from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  // write what PPM expects
  fprintf(fp, "P6\n%d %d\n255\n", xRes, yRes);
  fwrite(pixels, 1, totalCells * 3, fp);

  // clean up
  fclose(fp);
  delete[] pixels;
  cout << " Wrote out file " << filename.c_str() << endl;
}

///////////////////////////////////////////////////////////////////////
// Example of parsing command line arguments
///////////////////////////////////////////////////////////////////////
void parseArgs(int& argc, char**& argv)
{
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

  // try to convert them into doubles
  vector<double> commandDoubles(argc);
  for (int x = 0; x < argc; x++)
    commandDoubles[x] = atof(argv[x]);

  cout << " I tried to convert each into double and got: " << endl;
  for (int x = 0; x < argc; x++)
    cout << commandDoubles[x] << " ";
  cout << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  int xRes, yRes;
  float* values = NULL;

  if (argc != 2) {
    cout << " Expecting args " << argv[0] << " <input filename> " << endl;
    return 0;
  }

  readPPM(argv[1], xRes, yRes, values);

  // ADAPTED FROM STAFF SOLUTION OF CS50 FALL 19 PSET "filter/blur"
  for (int y = 0; y < yRes; y++) {
    for (int x = 0; x < xRes; x++) {

      int red = 0;
      int green = 0;
      int blue = 0;
      int num = 0;

      for (int r = MAX(0, y - 2); r <= MIN(yRes - 1, y + 2); r++) {
        for (int c = MAX(0, x - 2); c <= MIN(xRes - 1, x + 2); c++) {

          int index = 3 * (r * xRes) + c * 3;

          red += values[index];
          green += values[index+1];
          blue += values[index+2];
          num++;
        }
      }

      int in = 3 * (y * xRes) + x;

      values[in]     = round(red/num);
      values[in + 1] = round(green/num);
      values[in + 2] = round(blue/num);
    }
  }

  writePPM("filtered.ppm", xRes, yRes, values);

  return 0;
}
