#include "header.h"

double io::column_increment(const string& filename, uint32_t col) {
  /* determine the first increment of the given column in a numeric text file */

  ifstream file(filename);
  if (!file.is_open()) {
    cerr << "could not open file " + filename << "\n";
    exit(1);
  }

  double increment = 0;
  
  uint32_t lines_read = 0;
  double old = 0, now = 0;
  while (lines_read <= 2) {
    string skip_string; double skip_value;
    
    if (file.peek() == '#') {
      getline(file, skip_string);
      continue;
    }
    for (uint32_t i=0; i<col; ++i) file >> skip_value;
    
    file >> now;
    increment = now-old;
    old = now;
    lines_read++;
    getline(file, skip_string);
  }
  
  if (lines_read < 2) {
    cerr << "could not find increments of column " << col << " in file " + filename << "\n";
    exit(1);
  }

  return increment;
};