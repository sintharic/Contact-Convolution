#ifndef HEADER_H
#define HEADER_H
//endif at EOF

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

namespace io {
  template <typename T>
  vector<T> read_column(const string& filename, uint32_t col) {
    /* return data of the given column in a numeric text file */

    ifstream file(filename);
    if (!file.is_open()) {
      cerr << "could not open file " + filename << "\n";
      exit(1);
    }

    vector<T> result(0);
    double val = 0;
    while (!file.eof()) {
      // skip comment lines and previous columns
      double skip_value; string skip_string;
      if (file.peek() == '#') getline(file, skip_string);
      for (uint32_t i=0; i<col; ++i) file >> skip_value;
      
      if (file >> val) result.push_back(val);
      getline(file, skip_string);
    }
    
    file.close();
    return result;
  };

  double column_increment(const string&, uint32_t);
}

#endif