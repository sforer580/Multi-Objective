//
//  Test_Point.hpp
//  Multi-Objective_Project
//
//  Created by Scott S Forer on 5/3/17.
//  Copyright © 2017 Scott S Forer. All rights reserved.
//

#ifndef Test_Point_hpp
#define Test_Point_hpp

#include <stdio.h>

using namespace std;


class Test_Point
{
    friend class Parameters;
    friend class Agent;
    friend class F_val;
    friend class PaCcET;
    friend class EA;
    
protected:
    
    
public:
    vector<double> point;
    int dom;
    
private:
};

#endif /* Test_Point_hpp */
