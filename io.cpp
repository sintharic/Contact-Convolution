#include "header.h"

template <typename T>
vector<T> io::read_column(const string& filename, uint32_t col) {
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

double io::column_increment(const string& filename, uint32_t col) {
  /* 
  determines the first increment of the given column in a numeric text file 
  */

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

void io::write_vectors(const string& filename, vector< vector<double>* > list, const string& header) {
  /*
  writes numeric data to text file, where each vector in <list> represents a column.
  this assumes that all vectors in <list> have the same size.
  */

  uint32_t nData = (*list[0]).size();
  for (auto array : list) {
    if ((*array).size() != nData) {
      cerr << "io::write_vectors() cannot write vectors of different sizes!\n";
      exit(1);
    }
  }

  ofstream output(filename);
  output << "# "+header << "\n";
  if (list.size() == 1) {
    for (double val : *list[0]) output << val << "\n";
  }
  else {
    for (int iData = 0; iData < nData; ++iData) {
      output << (*list[0])[iData];
      for (int iArray = 1; iArray < list.size(); ++iArray) {
        output << "\t" << (*list[iArray])[iData];
      }
      output << "\n";
    }
  }
  output.close();
};