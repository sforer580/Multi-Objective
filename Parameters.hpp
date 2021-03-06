//
//  Parameters.hpp
//  Multi_objective_Project
//
//  Created by Scott S Forer on 1/19/17.
//  Copyright © 2017 Scott S Forer. All rights reserved.
//

#ifndef Parameters_hpp
#define Parameters_hpp

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


class Parameters
{
    friend class Agent;
    friend class EA;
    friend class F_val;
    
protected:
    
    
public:
    int num_agents = 100;
    int num_x_val = 1;
    int m = 20;
    double x_val_upper_limit = 1;
    double x_val_lower_limit = 0;
    double Xm_val_upper_limit = 1;
    double Xm_val_lower_limit = 0;
    int num_F = 2;
    int linear_combination = 0;                 //0=off, 1=on
    int volumetirc = 0;                         //0=off, 1=on
    int use_PaCcet = 0;                         //0=off, 1=on
    int use_quartet = 1;                        //0=off, 1=on
    vector<double> set_point;
    int num_tp = 100000;
    double set_point_val_0 = 1;
    double set_point_val_1 = 1;
    double set_point_val_2 = 1;
    double W1 = 1;
    double W2 = 1;
    double W3 = 1;
    int gen_max = 300;
    int to_kill = num_agents/2;
    double prob_mutate = 50;
    double mutate_range = 0.15;
    int num_sr = 1;

private:
};

#endif /* Parameters_hpp */
