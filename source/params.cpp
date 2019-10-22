#include <iostream>
#include <string>
#include "params.hpp"

using namespace std;

bool read_init_stream(istream& in, struct Simul_params& sm)
{
    bool success = false;
    if(in)
    {
        string temp;
        
        while(in >> temp)
        {
            if(temp == "p:")
            {
                sm.control += temp;
                in >> sm.p;
            }
            else if(temp == "q:")
            {
                sm.control += temp;
                in >> sm.q;
            }
            else if(temp == "dim:")
            {
                sm.control += temp;
                in >> sm.dim;
            }
            else if(temp == "L:")
            {
                sm.control += temp;
                in >> sm.L;
            }
            else if(temp == "dL:")
            {
                sm.control += temp;
                in >> sm.dL;
            }
            else if(temp == "dt:")
            {
                sm.control += temp;
                in >> sm.dt;
            }
            else if(temp == "ddt:")
            {
                sm.control += temp;
                in >> sm.ddt;
            }
            else if(temp == "AR:")
            {
                sm.control += temp;
                in >> sm.AR;
            }
            else if(temp == "dAR:")
            {
                sm.control += temp;
                in >> sm.dAR;
            }
            else if(temp == "M:")
            {
                sm.control += temp;
                in >> sm.M;
            }
            else if(temp == "scale:")
            {
                sm.control += temp;
                in >> sm.scale;
            }
            else if(temp == "iter_therm:")
            {
                sm.control += temp;
                in >> sm.iter_therm;
            }
            else if(temp == "iter_simul:")
            {
                sm.control += temp;
                in >> sm.iter_simul;
            }
            else if(temp == "gap:")
            {
                sm.control += temp;
                in >> sm.gap;
            }
            else if(temp == "adj:")
            {
                sm.control += temp;
                in >> sm.adj;
            }
            else if(temp == "g2_i:")
            {
                sm.control += temp;
                in >> sm.g2_i;
            }
            else if(temp == "g2_f:")
            {
                sm.control += temp;
                in >> sm.g2_f;
            }
            else if(temp == "g2_step:")
            {
                sm.control += temp;
                in >> sm.g2_step;
            }
            else if(temp == "mode:")
            {
                sm.control += temp;
                in >> sm.mode;
            }
        }

        success = true;

    }

    return success;
}

bool params_validity(const struct Simul_params& sm)
{
    if(sm.control.find("p:") == std::string::npos)
    {
        cerr << "p not found" << endl;
        return 0;
    }
    if(sm.control.find("q:") == std::string::npos)
    {
        cerr << "q not found" << endl;
        return 0;
    }
    if(sm.control.find("dim:") == std::string::npos)
    {
        cerr << "dim not found" << endl;
        return 0;
    }
    if(sm.control.find("g2_i:") == std::string::npos)
    {
        cerr << "g2_i not found" << endl;
        return 0;
    }
    
    if(sm.control.find("mode:") == std::string::npos)
    {
        cerr << "mode not found" << endl;
        return 0;
    }
    if(sm.control.find("iter_therm:") == std::string::npos)
    {
        cerr << "iter_therm not found" << endl;
        return 0;
    }
    if(sm.control.find("iter_simul:") == std::string::npos)
    {
        cerr << "iter_simul not found" << endl;
        return 0;
    }

    if(sm.mode == "fix_nosplit")
    {
        if(sm.control.find("L:") == std::string::npos)
        {
            cerr << "L not found" << endl;
            return 0;
        }
        if((sm.control.find("dt:") == std::string::npos) && (sm.control.find("AR:") == std::string::npos))
        {
            cerr << "dt and AR not found" << endl;
            return 0;
        }
    }
    
    if(sm.mode == "fix_split")
    {
        if(sm.control.find("L:") == std::string::npos)
        {
            cerr << "L not found" << endl;
            return 0;
        }
        if((sm.control.find("dt:") == std::string::npos) && (sm.control.find("AR:") == std::string::npos))
        {
            cerr << "dt and AR not found" << endl;
            return 0;
        }
    }

    if(sm.mode == "rand_nosplit")
    {
        if(sm.control.find("L:") == std::string::npos)
        {
            cerr << "L not found" << endl;
            return 0;
        }
        if((sm.control.find("dt:") == std::string::npos) && (sm.control.find("AR:") == std::string::npos))
        {
            cerr << "dt and AR not found" << endl;
            return 0;
        }
    }

    if(sm.mode == "rand_split")
    {
        if(sm.control.find("L:") == std::string::npos)
        {
            cerr << "L not found" << endl;
            return 0;
        }
        if((sm.control.find("dt:") == std::string::npos) && (sm.control.find("AR:") == std::string::npos))
        {
            cerr << "dt and AR not found" << endl;
            return 0;
        }
    }

    if(sm.mode == "mmc")
    {
        if((sm.control.find("scale:") == std::string::npos) && (sm.control.find("AR:") == std::string::npos))
        {
            cerr << "scale and AR not found" << endl;
            return 0;
        }
    }

    return 1;
}




