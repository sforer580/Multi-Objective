//
//  main.cpp
//  Multi_objective_Project
//
//  Created by Scott S Forer on 1/19/17.
//  Copyright © 2017 Scott S Forer. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <ctime>

#include "Agent.hpp"
#include "F_value.hpp"
#include "Parameters.hpp"
#include "PaCcET.h"
#include "Quartet.hpp"
#include "Test_Point.hpp"
#include "EA.hpp"


void test_A(EA E){
    Quartet Q;
    Quartet* pQ=&Q;
    
    E.tp.clear();
    E.Build_Hyper_Volume();
    
    vector<double> a = {2.5,2.5};
    Q.T.Pareto_Check(a);
    vector<double> b = {3.75,0};
    Q.T.Pareto_Check(b);
    vector<double> c = {0,3.75};
    Q.T.Pareto_Check(c);
    
    E.hyper_dom.clear();
    E.Run_Hyper_Volume_Check_Quartet(pQ);
    E.Write_Hyper_Dom_To_File();
}


int main()
{
    srand(time(NULL));
    Parameters P;
    EA E;
    E.pP = &P;
    E.Run_Multi_Objective();
    
    //test_A(E);
}
