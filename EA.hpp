//
//  EA.hpp
//  Multi-Objective_Project
//
//  Created by Scott S Forer on 1/19/17.
//  Copyright © 2017 Scott S Forer. All rights reserved.
//

#ifndef EA_hpp
#define EA_hpp

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


class EA
{
    friend class Parameters;
    friend class Agent;
    
protected:
    
    
public:
    Parameters* pP;
    vector <Agent> indv;
    
    //Data structure
    void Build_Pop();
    void Create_set_point();
    
    //Fitness functions
    void Function_0(int a);
    void Function_1(int a);
    void Function_2(int a);
    void Sub_Function(int a);
    void Get_Fitness();
    void Linear_Combination_Fitness(int a);
    void Volumetric_fitness(int a);
    
    //EA functions
    int Down_Select();
    void Mutation(Agent &M);
    int kill;
    void Natural_Selection();
    void Sort_indivduals_fitness();
    
    //Statistical functions
    struct Greater_than_agent_fitness;
    void Write_final_pop_to_file();
    
    //EA main
    void Run_Multi_Objective();
    
private:
};


//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------
//Builds population of individuals
void EA::Build_Pop()
{
    for (int a=0; a<pP->num_agents; a++)
    {
        Agent A;
        indv.push_back(A);
        indv.at(a).neg = 0;
        indv.at(a).x.resize(pP->num_x_val);
        for (int i=0; i<pP->num_x_val; i++)
        {
            double max = pP->x_val_upper_limit;
            double min = pP->x_val_lower_limit;
            double range = max-min;
            //randomly gives an x value between the upper and lower limit
            indv.at(a).x.at(i) = range*((double)rand()/RAND_MAX)+min;
        }
        indv.at(a).Xm.resize(pP->m);
        for (int i=0; i<pP->m; i++)
        {
            double max = pP->Xm_val_upper_limit;
            double min = pP->Xm_val_lower_limit;
            double range = max-min;
            //randomly gives an Xm value between the upper and lower limit
            indv.at(a).Xm.at(i) = range*((double)rand()/RAND_MAX)+min;
        }
        indv.at(a).F.resize(pP->num_F);
        for (int i=0; i<pP->num_F; i++)
        {
            indv.at(a).F.at(i) = 0;
        }
        
    }
}


//-------------------------------------------------------------------------
//Creates the set point for the volumetirc fitness evaluation
void EA::Create_set_point()
{
    pP->set_point.resize(pP->num_F);
    pP->set_point.at(0) = pP->set_point_val_0;
    pP->set_point.at(1) = pP->set_point_val_1;
    pP->set_point.at(2) = pP->set_point_val_2;
}


//-------------------------------------------------------------------------
//Function 0
void EA::Function_0(int a)
{
    indv.at(a).F.at(0) = 0;
    indv.at(a).F.at(0) = (1 + indv.at(a).GXm);
    for (int i=0; i<pP->num_x_val; i++)
    {
        if (i < pP->num_x_val-1)
        {
            indv.at(a).F.at(0) = indv.at(a).F.at(0)*cos((indv.at(a).x.at(i)*3.14159)/2);
        }
        if (i == pP->num_x_val-1)
        {
            indv.at(a).F.at(0) = indv.at(a).F.at(0)*cos((indv.at(a).x.at(i)*3.14159)/2);
        }
        assert(indv.at(a).F.at(0) >= 0);
    }
}


//-------------------------------------------------------------------------
//Function 1
void EA::Function_1(int a)
{
    indv.at(a).F.at(1) = 0;
    indv.at(a).F.at(1)= (1 + indv.at(a).GXm);
    for (int i=0; i<pP->num_x_val; i++)
    {
        if (i < pP->num_x_val-1)
        {
            indv.at(a).F.at(1) = indv.at(a).F.at(1)*cos((indv.at(a).x.at(i)*3.14159)/2);
        }
        if (i == pP->num_x_val-1)
        {
            indv.at(a).F.at(1) = indv.at(a).F.at(1)*sin((indv.at(a).x.at(i)*3.14159)/2);
        }
        assert(indv.at(a).F.at(1) >= 0);
    }
}


//-------------------------------------------------------------------------
//Function 2
void EA::Function_2(int a)
{
    indv.at(a).F.at(2) = 0;
    indv.at(a).F.at(2) = (1 + indv.at(a).GXm);
    indv.at(a).F.at(2) = indv.at(a).F.at(2)*sin((indv.at(a).x.at(0)*3.14159)/2);
    assert(indv.at(a).F.at(2) >= 0);
}


//-------------------------------------------------------------------------
//Sub function that calculates Xm
void EA::Sub_Function(int a)
{
    indv.at(a).GXm = 0;
    for (int i=0; i<pP->m; i++)
    {
        indv.at(a).GXm += (indv.at(a).Xm.at(i) - 0.5)*(indv.at(a).Xm.at(i) - 0.5);
    }
}


//-------------------------------------------------------------------------
//Gets fitness for each individual
void EA::Get_Fitness()
{
    for (int a=0; a<pP->num_agents; a++)
    {
        indv.at(a).fitness = 0;
        Sub_Function(a);
        Function_0(a);
        Function_1(a);
        Function_2(a);
        //linear combination
        if (pP->linear_combination == 1)
        {
            Linear_Combination_Fitness(a);
        }
        if (pP->volumetirc == 1)
        {
            Volumetric_fitness(a);
        }

    }
}


//-------------------------------------------------------------------------
//Volumetric fitness
void EA::Volumetric_fitness(int a)
{
    indv.at(a).fitness = 0;
    for (int i=0; i<pP->num_F; i++)
    {
        indv.at(a).neg = 0;         //resets negative identifier
        if (i == 0)
        {
            if (indv.at(a).F.at(i) > pP->set_point.at(i))
            {
                indv.at(a).neg = 1;
                indv.at(a).fitness = abs((pP->set_point.at(i) - indv.at(a).F.at(i)));
            }
            else
            {
                indv.at(a).fitness = abs((pP->set_point.at(i) - indv.at(a).F.at(i)));
            }
        }
        else
        {
            if (indv.at(a).F.at(i) > pP->set_point.at(i))
            {
                indv.at(a).neg = 1;
                indv.at(a).fitness = abs((pP->set_point.at(i) - indv.at(a).F.at(i)));
            }
            else
            {
                 indv.at(a).fitness = abs((pP->set_point.at(i) - indv.at(a).F.at(i)));
            }
        }
        if (indv.at(a).neg == 1)
        {
            indv.at(a).fitness = indv.at(a).fitness*-1;
        }
    }
}


//-------------------------------------------------------------------------
//Linear combintation with weights
void EA::Linear_Combination_Fitness(int a)
{
    indv.at(a).fitness = 0;
    indv.at(a).fitness = indv.at(a).F.at(0)*pP->W1 + indv.at(a).F.at(1)*pP->W2 + indv.at(a).F.at(2)*pP->W3;
    //swithces to a maximization problem
    indv.at(a).fitness = 1/indv.at(a).fitness;
}


//-------------------------------------------------------------------------
//Down selects by binarly selecting two indivduals, comapring their fitness and returns the loser
int EA::Down_Select()
{
    int loser;
    int index_1 = rand() % indv.size();
    int index_2 = rand() % indv.size();
    while (index_1 == index_2)
    {
        index_2 = rand() % indv.size();
    }
    //indvidual with higher fitness wins
    if(indv.at(index_1).fitness > indv.at(index_2).fitness)
    {
        loser = index_2;
    }
    else
    {
        loser = index_1;
    }
    return loser;
}


//-------------------------------------------------------------------------
//Mutates the copies of the winning individuals
void EA::Mutation(Agent &M)
{
    //mutates x values
    for (int i=0; i<pP->num_x_val; i++)
    {
        double R = 0;
        R = ((double) rand() / (RAND_MAX));
        if (R <= pP->prob_mutate/100)
        {
            //triangular distribution
            double R1 = ((double)rand()/RAND_MAX) * pP->mutate_range;
            double R2 = ((double)rand()/RAND_MAX) * pP->mutate_range;
            M.x.at(i) = M.x.at(i) + (R1-R2);
            if (M.x.at(i) < pP->x_val_lower_limit)
            {
                M.x.at(i) = pP->x_val_lower_limit;
            }
            if (M.x.at(i) > pP->x_val_upper_limit)
            {
                M.x.at(i) = pP->x_val_upper_limit;
            }
        }
    }
    //muttaes Xm values
    for (int i=0; i<pP->m; i++)
    {
        double R = 0;
        R = ((double) rand() / (RAND_MAX));
        if (R <= pP->prob_mutate/100)
        {
            //triangular distribution
            double R1 = ((double)rand()/RAND_MAX) * pP->mutate_range;
            double R2 = ((double)rand()/RAND_MAX) * pP->mutate_range;
            M.Xm.at(i) = M.Xm.at(i) + (R1-R2);
            if (M.Xm.at(i) < pP->Xm_val_lower_limit)
            {
                M.Xm.at(i) = pP->Xm_val_lower_limit;
            }
            if (M.Xm.at(i) > pP->Xm_val_upper_limit)
            {
                M.Xm.at(i) = pP->Xm_val_upper_limit;
            }
        }
    }
}


//-------------------------------------------------------------------------
//Runs the entire natural selection process
void EA::Natural_Selection()
{
    for(int k=0; k<pP->to_kill; k++)
    {
        kill = Down_Select();
        //kills off the loserr from the down select
        indv.erase(indv.begin() + kill);
    }
    int to_replicate = pP->to_kill;
    for (int r=0; r<to_replicate; r++)
    {
        Agent M;
        int spot = rand() % indv.size();
        M = indv.at(spot);
        Mutation(M);
        indv.push_back(M);
    }
}

//-------------------------------------------------------------------------
//Sorts the population based on their fitness from lowest to highest
struct EA::Greater_than_agent_fitness
{
    inline bool operator() (const Agent& struct1, const Agent& struct2)
    {
        return (struct1.fitness < struct2.fitness);
    }
};


//-------------------------------------------------------------------------
//Runs sort function for entire population
void EA::Sort_indivduals_fitness()
{
    for (int a=0; a<pP->num_agents; a++)
    {
        sort(indv.begin(), indv.end(), Greater_than_agent_fitness());
    }
}


//-------------------------------------------------------------------------
//Writes final population to a txt file
void EA::Write_final_pop_to_file()
{
    ofstream File1;
    File1.open("Final_Population_Fitness.txt");
    ofstream File2;
    File2.open("Final_Population_Function_Values.txt");
    for (int a=0; a<pP->num_agents; a++)
    {
        File1 << indv.at(a).fitness << endl;
        for (int i=0; i<pP->num_F; i++)
        {
            File2 << indv.at(a).F.at(i) << "\t";
        }
        File2 << endl;
    }
    File1.close();
    File2.close();
}


//-------------------------------------------------------------------------
//Runs entire multi-objective problem
void EA::Run_Multi_Objective()
{
    Build_Pop();
    Create_set_point();
    for (int gen=0; gen<pP->gen_max; gen++)
    {
        if (gen < pP->gen_max-1)
        {
            cout << "Generation" << "\t" << gen << endl;
            Get_Fitness();
            Sort_indivduals_fitness();
            cout << "Best Individual" << endl;
            cout << "Fitness" << "\t" << indv.at(0).fitness << endl;
            cout << "Xm values" << "\t";
            for (int i=0; i<pP->m; i++)
            {
                cout << indv.at(0).Xm.at(i) << "\t";
            }
            cout << endl;
            cout << "Function values" << "\t";
            for (int i=0; i<pP->num_F; i++)
            {
                cout << indv.at(0).F.at(i) << "\t";
            }
            cout << endl;
            Natural_Selection();
            cout << endl;
            cout << "--------------------------------------------------------------------" << endl;
        }
        if (gen == pP->gen_max-1)
        {
            cout << "Final Generation" << "\t" << gen << endl;
            Get_Fitness();
            Sort_indivduals_fitness();
            cout << "Best Individual" << endl;
            cout << "Fitness" << "\t" << indv.at(0).fitness << endl;
            cout << "Xm values" << "\t";
            for (int i=0; i<pP->m; i++)
            {
                cout << indv.at(0).Xm.at(i) << "\t";
            }
            cout << endl;
            cout << "Function values" << "\t";
            for (int i=0; i<pP->num_F; i++)
            {
                cout << indv.at(0).F.at(i) << "\t";
            }
            cout << endl;
            Write_final_pop_to_file();
        }
    }
}

#endif /* EA_hpp */
