#ifndef JSONHELPER_HPP_
#define JSONHELPER_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <armadillo>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <stack>
#include <vector>
#include "json.hpp"

using json = nlohmann::json;
using namespace arma;
using namespace std;

namespace json_helper
{

class Jsonhelper
{
public:
    // Jsonhelper();
    // Jsonhelper(string path);
    // Jsonhelper(json jsn);
    //
    // json parse(string path);
    // void save(string path);
    // bool is_empty(json jsonObj);
    //
    // void get_dimension(json jsonObj, int& dim);
    // void get_size(json jsonObj, vector<unsigned int>& sze);
    //
    // template<typename T>
    // void json2mat(string key,Mat<T>& m);
    //
    // template<typename T>
    // void json2mat(Mat<T>& m);
    //
    // template<typename T>
    // void json2vector(string key,vector<T>& v);
    //
    // template<typename T>
    // void json2vector2D(string key,vector<vector<T>>& v);
    //
    // template<typename T>
    // void json2vector3D(string key,vector<vector<vector<T>>>& v);
    //
    // template<typename T>
    // void mat2json(string key,Mat<T> m);

    Jsonhelper()
    {
    }

    Jsonhelper(string path)
    {
        j = parse(path);
    }

    Jsonhelper(json jsn)
    {
        j = jsn;
    }

    json parse(string path)
    {
        ifstream i(path);
        j = json::parse(i);
        return j;
    }

    void save(string path)
    {
        ofstream o(path, std::ios::out | std::ios::trunc);
        o << std::setw(4) << j << endl;  // print pretty-> readable
        o.close();
    }

    bool is_empty(json jsonObj)
    {
        return jsonObj.empty();
    }

    void get_dimension(json jsonObj, int& dim)
    {
        /* saves the dimension of a multidimensional array contained n jsonObj.
         * jsonObj must have same amount of elements in each dimension.*/
        if (jsonObj.is_primitive())
        {
            dim = 1;
            //			cout << "primitive"<<endl;
            //			cout << setw(2)<<jsonObj<<endl;
        }
        else
        {
            dim = 0;
            while (!jsonObj.empty() && !jsonObj.is_primitive())
            {  // jsonObj.size()>1){
                try
                {
                    jsonObj = jsonObj.at(0);
                    //					cout << "trying"<<endl;
                    //					cout <<
                    // setw(2)<<jsonObj<<endl;
                }
                catch (exception& e)
                {
                    cerr << e.what() << '\n';
                    cerr << std::setw(4) << jsonObj << std::endl;
                    dim--;
                }
                dim++;
            }
        }
    }

    void get_size(json jsonObj, vector<unsigned int>& sze)
    {
        /* saves the number of elements of a multidimensional array contained n
         * jsonObj. jsonObj must have same amount of elements in each
         * dimension.*/
        // check if is empty
        if (!this->is_empty(jsonObj))
        {
            int dim;
            this->get_dimension(jsonObj, dim);
            // print for testing
            //			cout <<endl<<endl<< "dim="<<dim<< "   ";
            //			cout<<jsonObj.dump();
            //			cout <<"   ["<<endl;
            for (int i_dim = 1; i_dim <= dim; i_dim++)
            {
                sze.push_back(jsonObj.size());
                if (i_dim + 1 <= dim) jsonObj = jsonObj.at(0);
                //				cout << sze.at(i_dim-1) <<
                //"/"<<endl;
            }
        }
    }

    /*---*/
    template <typename T>
    void json2mat(string key, Mat<T>& m)
    {
        if (j.is_object())
        {
            if (!this->is_empty(j[key]))
            {
                // find size
                vector<unsigned int> sze;
                this->get_size(j[key], sze);
                // print for testing
                //								cout
                //<<"sze:
                //"; for(auto const& value: sze) cout << value<<" "; cout
                //<<"//"<<endl;
                if (sze.size() > 2)
                    std::cerr << "json2mat: json obj " << key
                              << "  must be <=2 but is " << sze.size()
                              << " dimensional!\n"
                              << endl;
                else
                {
                    if (sze.size() == 1) sze.push_back(1);
                    m.zeros(sze.at(0), sze.at(1));
                    // fill
                    for (unsigned int itr = 0; itr < sze.at(0); itr++)
                    {
                        for (unsigned int itc = 0; itc < sze.at(1); itc++)
                        {
                            if (sze.at(0) == 1 && sze.at(1) == 1)
                                m(itr, itc) = j[key];
                            if (sze.at(0) > 1 && sze.at(1) == 1)
                                m(itr, itc) = j[key].at(itr);
                            if (sze.at(0) > 1 && sze.at(1) > 1)
                                m(itr, itc) = j[key].at(itr).at(itc);
                        }
                    }
                }
            }
            else
                m.zeros(0, 0);
            //						m.print(key);
        }
        else
        {
            std::cerr << "Jsonhelper: json not initialized!\n" << endl;
        }
    }
    /*---*/

    template <typename T>
    void json2mat(Mat<T>& m)
    {
        if (j.is_object())
        {
            if (!this->is_empty(j))
            {
                // find size
                vector<unsigned int> sze;
                this->get_size(j, sze);
                // print for testing
                //				cout <<"sze: ";
                //				for(auto const& value: sze) cout
                //<< value<<" "; 				cout
                // <<"//"<<endl;
                if (sze.size() > 2)
                    std::cerr << "json2mat: json obj " << j.dump()
                              << "  must be <=2 but is " << sze.size()
                              << " dimensional!\n"
                              << endl;
                else
                {
                    if (sze.size() == 1) sze.push_back(1);
                    m.zeros(sze.at(0), sze.at(1));
                    // fill
                    for (unsigned int itr = 0; itr < sze.at(0); itr++)
                    {
                        for (unsigned int itc = 0; itc < sze.at(1); itc++)
                        {
                            if (sze.at(0) == 1 && sze.at(1) == 1)
                                m(itr, itc) = j;
                            if (sze.at(0) > 1 && sze.at(1) == 1)
                                m(itr, itc) = j.at(itr);
                            if (sze.at(0) > 1 && sze.at(1) > 1)
                                m(itr, itc) = j.at(itr).at(itc);
                        }
                    }
                }
            }
            else
                m.zeros(0, 0);
            //			m.print(key);
        }
        else
        {
            std::cerr << "Jsonhelper: json not initialized!\n" << endl;
        }
    }

    /*---*/
    template <typename T>
    void json2vector(string key, vector<T>& v)
    {
        v.clear();
        if (j.is_object())
        {
            // find size
            vector<unsigned int> sze;
            this->get_size(j[key], sze);
            // print for testing
            // cout <<"sze: ";
            // for(auto const& value: sze1) cout << value<<" ";
            // cout <<"//"<<endl;
            if (sze.size() > 1)
            {
                std::cerr << "json2vec: json obj " << key
                          << "  must be 1 but is " << sze.size()
                          << " dimensional!\n"
                          << endl;
            }
            else
            {
                // fill
                for (unsigned int itr = 0; itr < sze.at(0); itr++)
                {
                    if (sze.at(0) == 1 && !j[key].is_array())
                        v.push_back(j[key]);
                    else
                        v.push_back(j[key].at(itr));
                }
                // print for testing
                //				cout <<key <<": ";
                //				for(auto const& value: v) cout
                //<< value<<" "; 				cout
                // <<"//"<<endl;
            }
        }
        else
        {
            std::cerr << "Jsonhelper: json not initialized!\n" << endl;
        }
    }

    /*---*/
    template <typename T>
    void json2vector2D(string key, vector<vector<T>>& v)
    {
        v.clear();
        if (j.is_object())
        {
            // find size
            vector<unsigned int> sze;
            this->get_size(j[key], sze);
            // print for testing
            //									cout <<"sze
            //of
            //"<<key<<":
            //"; 									for(auto const& value: sze) cout <<
            //value<<"
            //";
            // cout
            //<<"//"<<endl;
            if (!(sze.size() == 2))
            {
                std::cerr << "json2vec2D: json obj " << key
                          << "  must be 2 but is " << sze.size()
                          << " dimensional!\n"
                          << endl;
            }
            else
            {
                // fill
                //				cout << setw(4)<<j[key]<<endl;
                for (unsigned int itr = 0; itr < sze.at(0); itr++)
                {
                    v.push_back(vector<T>());
                    for (unsigned int itc = 0; itc < sze.at(1); itc++)
                    {
                        if (sze.at(1) == 1 && !j[key].at(itr).is_array())
                            v[itr].push_back(j[key].at(itr));
                        else
                            v[itr].push_back(j[key].at(itr).at(itc));
                    }
                    // print for testing
                    //					cout <<key <<itr<<": ";
                    //					for(auto const& value: v[itr]) cout
                    //<< value<<"
                    //"; 					cout
                    //<<"//"<<endl;
                }
                // print for testing
                //								cout
                //<<key
                //<<": "<<endl;
                // for(unsigned int i_v1=0;i_v1<v.size();i_v1++) {
                // cout<<"|
                //";
                // for(unsigned int i_v2=0;i_v2<v.at(i_v1).size();i_v2++){
                // cout
                //<< v[i_v1][i_v2]<<" ";
                //									}
                //									cout<<"|"<<endl;
                //								}
            }
        }
        else
        {
            std::cerr << "Jsonhelper: json not initialized!\n" << endl;
        }
    }

    /*---*/
    template <typename T>
    void json2vector3D(string key, vector<vector<vector<T>>>& v)
    {
        v.clear();
        if (j.is_object())
        {
            // find size
            vector<unsigned int> sze;
            this->get_size(j[key], sze);
            // print for testing
            //			cout <<"sze: ";
            //			for(auto const& value: sze) cout << value<<" ";
            //			cout <<"//"<<endl;
            if (!(sze.size() == 3))
            {
                std::cerr << "json2vec3D: json obj " << key
                          << "  must be 3 but is " << sze.size()
                          << " dimensional!\n"
                          << endl;
            }
            else
            {
                // fill
                for (unsigned int it1 = 0; it1 < sze.at(0); it1++)
                {
                    v.push_back(vector<vector<T>>());
                    for (unsigned int it2 = 0; it2 < sze.at(1); it2++)
                    {
                        v.at(it1).push_back(vector<T>());
                        for (unsigned int it3 = 0; it3 < sze.at(2); it3++)
                        {
                            if (sze.at(2) == 1 &&
                                !j[key].at(it1).at(it2).is_array())
                                v[it1][it2].push_back(j[key].at(it1).at(it2));
                            else
                                v[it1][it2].push_back(
                                    j[key].at(it1).at(it2).at(it3));
                        }
                        // print for testing
                        //						cout
                        //<<key
                        //<<it1<<
                        //"/"<<it2<<": ";
                        // for(auto const& value: v[it1][it2]) cout << value<<"
                        // "; 						cout
                        // <<"//"<<endl;
                    }
                }
                // print for testing
                //								cout
                //<<key
                //<<": "<<endl;
                // for(unsigned int i_v1=0;i_v1<v.size();i_v1++) {
                // cout<<i_v1<< " "; for(unsigned int
                // i_v2=0;i_v2<v.at(i_v1).size();i_v2++){ 										cout<<i_v2<< " |";
                // for(unsigned int
                // i_v3=0;i_v3<v.at(i_v1).at(i_v2).size();i_v3++){
                // cout << v[i_v1][i_v2][i_v3]<<" ";
                //										}
                //										cout<<"|"<<endl;
                //									}
                //									cout<<"|"<<endl;
                //								}
            }
        }
        else
        {
            std::cerr << "Jsonhelper: json not initialized!\n" << endl;
        }
    }

    /*---*/
    template <typename T>
    void mat2json(string key, Mat<T> m)
    {
        if (m.n_rows == 1 && m.n_cols > 1)
        {
            m = m.t();
            j[key] = m;
            m = m.t();
        }
        if (m.n_rows > 1 && m.n_cols == 1)
        {
            j[key] = m;
        }
        if (m.n_rows > 1 && m.n_cols > 1)
        {
            j[key] = {(Mat<T>)m(0, span::all), (Mat<T>)m(1, span::all)};
            for (unsigned int i = 2; i < m.n_rows; i++)
            {
                j[key] += (Mat<T>)m(i, span::all);
            }
        }
        if (m.n_rows == 1 && m.n_cols == 1) j[key] = m;
    }

public:
    json j;
};

}  // namespace json_helper

#endif
