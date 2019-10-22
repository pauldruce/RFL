#include <sstream>
#include <iostream>
#include <algorithm>
#include <ctime>
#include "utils.hpp"
#include "geometry.hpp"
#include "clifford.hpp"

using namespace std;


string filename_from_data(const int& p, const int& q, const int& dim, const double& g2, const string& prefix)
{
    string sp = to_string(p);
    string sq = to_string(q);
    string sdim = to_string(dim);
    
    ostringstream osg2;
    osg2 << g2;
    string sg2 = osg2.str();
    replace(sg2.begin(), sg2.end(), '.', 'd');
    
    return prefix + "p" + sp + "q" + sq + "dim" + sdim + "g" + sg2;
}

string basename_from_data(const int& p, const int& q, const int& dim, const string& prefix)
{
    string sp = to_string(p);
    string sq = to_string(q);
    string sdim = to_string(dim);
    
    return prefix + "p" + sp + "q" + sq + "dim" + sdim;
}

void data_from_filename(const string& s, int& p, int& q, int& dim, double& g2, const string& prefix)
{
    // generate a clean string s1 from s that
    // ignores everything before the first appearence
    // of prefix (if found)
    size_t start = s.find(prefix);
    string s1;
    if(start != string::npos)
        s1 = s.substr(start);
    else
        s1 = s;

    // extract data from string.
    // data is formatted as follows:
    //
    // GEOMp[int]q[int]dim[int]g[int]d[int][garbage].txt
    // 
    // where [int] denotes an integer numerical value.
    // the last part g[int]d[int] is a double written as
    // integer part and decimal part separated by a letter d
    
    string s_p = s1.substr(s1.find("p")+1, s1.find("q")-1);
    string s_q = s1.substr(s1.find("q")+1, s1.find("d")-s1.find("q")-1);
    string s_dim = s1.substr(s1.find("m")+1, s1.find("g")-s1.find("m")-1);
    string s_g = s1.substr(s1.find("g")+1, s1.find(".")-s1.find("g")-1);
    replace(s_g.begin(), s_g.end(), 'd', '.');

    p = stoi(s_p);
    q = stoi(s_q);
    dim = stoi(s_dim);
    g2 = stod(s_g);
}


string foldername_from_time(const time_t& t)
{
    tm* timePtr = localtime(&t);

    string day;
    if(timePtr->tm_mday < 10)
        day = "0" + to_string(timePtr->tm_mday);
    else
        day = to_string(timePtr->tm_mday);
    
    string month;
    if(timePtr->tm_mon < 9)
        month = "0" + to_string(timePtr->tm_mon + 1);
    else
        month =to_string(timePtr->tm_mon + 1);

    string year = to_string(timePtr->tm_year + 1900);

    return year+month+day+"/";
}


int n_meas(const int& n_tot, const int& gap)
{
    int res = 0;

    for(int i=0; i<n_tot; ++i)
    {
        if( !(i%gap) )
            ++res;
    }

    return res;
}



void dofs_analysis(const string& path)
{
    // extract geometric data
    int p, q, dim;
    double g2;
    data_from_filename(path, p, q, dim, g2, "PRELIM");
        
    
    // read quadratic and quartic part from file
    // and store separately in vectors
    ifstream in_s;
    in_s.open(path + "_S.txt");
    
    vector<double> vec2;
    vector<double> vec4;
    double temp2, temp4;

    while(in_s >> temp2 >> temp4)
    {
        vec2.push_back(temp2);
        vec4.push_back(temp4);
    }
    size_t len = vec2.size();

    in_s.close();

    
    // calculate dofs
    Geom24 G(p, q, 1, 1);
    int c = dim*dim*G.get_nHL();

    // calculate average of 2gTrD2 + 4TrD4 based on the last
    // 1/10th of samples and store it in path_dofs.txt
    ofstream out_dofs;
    out_dofs.open(path + "_dofs.txt");

    for(size_t i=0; i<(len-(len/10)); ++i)
    {
        double res = 0;
        for(size_t j=0; j<len/10; ++j)
            res += 2*g2*vec2[i+j] + 4*vec4[i+j];

        out_dofs << 10*res/len << " " << c << endl;
    }
    out_dofs.close();
}
