//
//  Quartet.hpp
//  Multi-Objective_Project
//
//  Created by Scott S Forer on 4/26/17.
//  Copyright © 2017 Scott S Forer. All rights reserved.
//

#ifndef Quartet_hpp
#define Quartet_hpp

#include <stdio.h>

class Quartet
{
    friend class Parameters;
    friend class Agent;
    friend class F_val;
    friend class PaCcET;
    friend class EA;
    friend class Test_Point;
    
protected:
    
    
public:
    PaCcET T;
    PaCcET* pT = &T;
    vector<vector<double>> to_be_updated;
    vector<vector<int>> quartet_counter;
    
    int counter_1 = 0;
    int counter_2 = 0;
    int counter_P = 0;
    vector<int> Q_1;
    vector<int> Q_2;
    vector<int> Q_P;
    
    double Quartet_main(vector<double> coord);
    void End_Generation();
    
private:
};

//Runs the end of generation functions for the Quartet
void Quartet::End_Generation()
{
    //cout << "TBU PROCESSING" << endl;
    for(int i=0; i<to_be_updated.size(); i++)
    {
        T.Pareto_Check(to_be_updated.at(i));
    }
    to_be_updated.clear();
    assert(to_be_updated.size() == 0);
    vector<int> V;
    V.push_back(counter_1);
    V.push_back(counter_2);
    V.push_back(counter_P);
    quartet_counter.push_back(V);
    Q_1.push_back(counter_1);
    Q_2.push_back(counter_2);
    Q_P.push_back(counter_P);
    counter_1 = 0;
    counter_2 = 0;
    counter_P = 0;
}

//Runs the Quartet main
double Quartet::Quartet_main(vector<double> coord)
{
    double fit = -1991;
    
    double quad1value = -100000000;
    double quad2value = -10000;
    double quadPValue = 10000;
    
    to_be_updated.push_back(coord);
    
    static int AAA;
    if(AAA=0){
    vector<double> BBB;
    BBB.push_back(5);
    BBB.push_back(5);
    BBB.push_back(5);
    T.Pareto_Check(BBB);
        AAA = 1;
    }
    
    if(T.get_PFront_size() == 0){
        counter_1 +=1;
        return quad1value;
    }
    
    /// coords vs. utopia. 3 cases.
    vector<double> p1 = coord;
    vector<double> p2 = T.get_ut();
    //T.does_v1_dominate_v2(v1, v2);
    // quad 1; return quad1value
    if (p1.at(0)<p2.at(0) && p1.at(1)<p2.at(1))
    {
        fit = quad1value + (double)rand()/RAND_MAX;
        counter_1 +=1;
    }
    // quad2 return quad3value
    else if (p1.at(0)>p2.at(0) && p1.at(1)<p2.at(1))
    {
        fit = quad2value + (double)rand()/RAND_MAX;
        counter_2 +=1;
    }
    else if (p1.at(0)<p2.at(0) && p1.at(1)>p2.at(1))
    {
        fit = quad2value + (double)rand()/RAND_MAX;
        counter_2 +=1;
    }
    else if (p1.at(0)<p2.at(0) && p1.at(1)<p2.at(1))
    {
        fit = quad2value+ (double)rand()/RAND_MAX;
        counter_2 +=1;
    }
    // quad3 return paccet linear combination of thingy.
    else if (p1.at(0)>p2.at(0) && p1.at(1)>p2.at(1))
    {
        vector<double> MO;
        vector<double>* pMO = &MO;
        MO = coord;
        vector<double> OMO = MO;
        T.execute_N_transform(pMO);
        fit = MO.at(0)+MO.at(1);
        counter_P +=1;
    }
    else{
        fit = quad2value + (double)rand()/RAND_MAX;
        counter_2 ++;
    }
    // put coords into "to be updated" list.
    
    return fit;
}


#endif /* Quartet_hpp */
