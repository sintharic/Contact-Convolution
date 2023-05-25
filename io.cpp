#include "header.h"

template <typename T>
vector<T> io::read_column(const string& filename, uint32_t col) {
/* 
  return data of the given column in a numeric text file 
*/

  ifstream file(filename);
  if (!file.is_open()) {
    terminate("io::read_column() could not open file " + filename);
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
    terminate("io::column_increment() could not open "+filename);
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
    terminate("could not find increments of column "+to_string(col)+" in "+filename);
  }

  return increment;
};



void io::write_vectors(const string& filename, vector< vector<double>* > list, const string& header) {
/*
  writes numeric data to text file, where each vector in <list> represents a column.
  this assumes that all vectors in <list> have the same size.
*/

  uint32_t nVec = list.size(), nData = (*list[0]).size();
  for (uint32_t iVec = 1; iVec < nVec; ++iVec) {
    if ((*list[iVec]).size() != nData) {
      terminate("io::write_vectors() cannot write vectors of different sizes!");
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
      for (int iVec = 1; iVec < list.size(); ++iVec) {
        output << "\t" << (*list[iVec])[iData];
      }
      output << "\n";
    }
  }
  output.close();
};



void io::write_array(const string& filename, vector< vector<double> >& array, const string& header) {
/*
  writes numeric data to text file, where each vector in <array> represents a column.
  this assumes that all vectors in <array> have the same size.
*/

  uint32_t nVec = array.size(), nData = array[0].size();
  for (uint32_t iVec = 1; iVec < nVec; ++iVec) {
    if (array[iVec].size() != nData) {
      terminate("io::write_array() cannot write vectors of different sizes!");
    }
  }

  ofstream output(filename);
  output << "# "+header << "\n";
  for (int iVec = 0; iVec < nVec; ++iVec) {
    output << array[iVec][0];
    for (int iData = 1; iData < nData; ++iData) {
      output << "\t" << array[iVec][iData];
    }
    output << "\n";
  }
  output.close();
};



const char DELIMITER = '=', COMMENT = '#';

vector< map<string,string> > io::read_pairs(const string& filename) {
  // inspired by https://www.walletfox.com/course/parseconfigfile.php
  map<string, string> global, elastic, indenter, inter;
  map<string, string> *active = &global;
  ifstream input(filename); 
  
  if (!input.is_open()) terminate("io::read_pairs() could not open "+filename+".");
  
  string line; 
  while (getline(input, line)) {
    // erase leading and trailing whitespaces
    line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end()); 
    
    // ignore empty lines and comments, i.e. lines starting with '#'
    if(line[0] == COMMENT || line.empty()) continue;

    if (line.find("ElasticBody:") != string::npos) {active = &elastic; continue;}
    if (line.find("Indenter:") != string::npos) {active = &indenter; continue;}
    if (line.find("Interaction:") != string::npos) {active = &inter; continue;}

    // find {name,value} pairs separated by delimiter
    size_t delim = line.find(DELIMITER);
    string name = line.substr(0, delim); string value = line.substr(delim + 1); 
    for (const auto& [key, subs] : SUBS) {
      if (key == value) value = subs;
    }
    (*active)[name] = value;
  }
  input.close();

  // add dTime to elastic params for Maxwell parameter checking
  if (global.find("dTime") != global.end()) elastic["dTime"] = global["dTime"];

  return {global, elastic, indenter, inter};
}