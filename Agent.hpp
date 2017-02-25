//
//  Agent.hpp
//  Multi_objective_Project
//
//  Created by Scott S Forer on 1/19/17.
//  Copyright © 2017 Scott S Forer. All rights reserved.
//

#ifndef Agent_hpp
#define Agent_hpp

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <algorithm>
#include <fstream>
#include <iomanip>


using namespace std;


class Agent
{
    friend class Parameters;
    friend class EA;
    friend class F_val;
    
protected:
    
    
public:
    vector <double> x;
    double GXm;
    vector<double> Xm;
    vector<double> F;
    double fitness;
    int neg;
    double length;
private:
};

#endif /* Agent_hpp */
