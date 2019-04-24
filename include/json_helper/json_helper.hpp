#ifndef JSONHELPER_HPP_
#define JSONHELPER_HPP_

#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <armadillo>
#include <iostream>
#include "json.hpp"
#include <stack>
#include <ctime>
#include <vector>
#include <iomanip>

using json = nlohmann::json;
using namespace arma;
using namespace std;

namespace json_helper {

  class Jsonhelper{

  public:

    Jsonhelper();
    Jsonhelper(string path);
    Jsonhelper(json jsn);

    json parse(string path);
    void save(string path);
    bool is_empty(json jsonObj);

    void get_dimension(json jsonObj, int& dim);
    void get_size(json jsonObj, vector<unsigned int>& sze);

    template<typename T>
    void json2mat(string key,Mat<T>& m);

    template<typename T>
    void json2mat(Mat<T>& m);

    template<typename T>
    void json2vector(string key,vector<T>& v);

    template<typename T>
    void json2vector2D(string key,vector<vector<T>>& v);

    template<typename T>
    void json2vector3D(string key,vector<vector<vector<T>>>& v);

    template<typename T>
    void mat2json(string key,Mat<T> m);

  public:
    json j;
  
  };

}

#endif
