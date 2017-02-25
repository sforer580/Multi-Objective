//
//  F_value.hpp
//  Multi-Objective_Project
//
//  Created by Scott S Forer on 2/24/17.
//  Copyright © 2017 Scott S Forer. All rights reserved.
//

#ifndef F_value_hpp
#define F_value_hpp

#include <stdio.h>

using namespace std;


class F_value
{
    friend class Parameters;
    friend class EA;
    friend class Agent;
    
protected:
    
    
public:
    vector<double> F;
    int optimal = 1;            //1=optimal, 0=dominated
    
private:
};

#endif /* F_value_hpp */
