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
#include "EA.hpp"

int main()
{
    srand(time(NULL));
    Parameters P;
    EA E;
    E.pP = &P;
    E.Run_Multi_Objective();
}
