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
    friend class F_val;
    friend class PaCcET;
    friend class Test_Point;
    
protected:
    
    
public:
    Parameters* pP;
    vector<Agent> indv;
    vector<F_value> point;
    vector<Test_Point> tp;
    vector<vector<double> > hyper_val;
    int num_hyper_dom;
    vector<int> hyper_dom;
    
    //Data structure
    void Build_Pop();
    void Create_set_point();
    
    //Fitness functions
    void Function_0(int a);
    void Function_1(int a);
    void Function_2(int a);
    void Sub_Function(int a);
    void Get_Fitness(Quartet* pQ);
    void Get_Quartet_Fitness(Quartet *pQ);
    void Get_PaCcET_Fitness(PaCcET *pT);
    void PaCcET_Fitness(PaCcET* pT, int a);
    void Run_Quartet(Quartet *pQ, int a);
    void Linear_Combination_Fitness(int a);
    void Volumetric_fitness(int a);
    void Build_Hyper_Volume();
    void Run_Hyper_Volume_Check_Quartet(Quartet *pQ);
    void Run_Hyper_Volume_Check_PaCcET(PaCcET *pT);
    
    //EA functions
    int Down_Select();
    void Mutation(Agent &M);
    int kill;
    void Natural_Selection();
    void Sort_indivduals_fitness();
    void Output_Best_Individual_Info();
    
    //Statistical functions
    struct Greater_than_agent_fitness;
    void Write_final_pop_to_file();
    void Store_f_values(int gen);
    void Compare_Points(int p, int pp);
    void Find_Pareto_Optimal_Points();
    void Write_Pareto_Optimal_Points_To_File();
    void Write_Counter_File(Quartet* pQ);
    void Write_Hyper_Dom_To_File();
    void Delete_Files();
    
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
            //double max = pP->x_val_upper_limit;
            //double min = pP->x_val_lower_limit;
            double max = 0.51;
            double min = 0.49;
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
        assert(indv.at(a).x.size()==pP->num_x_val);
        assert(indv.at(a).Xm.size()==pP->m);
        assert(indv.at(a).F.size()==pP->num_F);
    }
    assert(indv.size()==pP->num_agents);
}


//-------------------------------------------------------------------------
//Creates the set point for the volumetirc fitness evaluation
void EA::Create_set_point()
{
    pP->set_point.resize(pP->num_F);
    pP->set_point.at(0) = pP->set_point_val_0;
    pP->set_point.at(1) = pP->set_point_val_1;
    //pP->set_point.at(2) = pP->set_point_val_2;
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
//Gets PaCcET fitness for each individual
void EA::PaCcET_Fitness(PaCcET* pT, int a)
{
    vector<double> MO;
    vector<double>* pMO = &MO;
    MO = indv.at(a).F;
    vector<double> OMO = MO;
    pT->execute_N_transform(pMO);
    indv.at(a).fitness = MO.at(0) + MO.at(1);
}



//-------------------------------------------------------------------------
//Runs the Quartet process
void EA::Run_Quartet(Quartet *pQ, int a)
{
    vector<double> coord;
    coord = indv.at(a).F;
    indv.at(a).fitness = pQ->Quartet_main(coord);
}


//-------------------------------------------------------------------------
//Build the hyper volume space
void EA::Build_Hyper_Volume()
{
    
    ifstream ifile("test_points.txt", ios::in);
    vector<double> sw;
    
    //check to see that the file was opened correctly:
    if (!ifile.is_open())
    {
        cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }
    
    double num = 0.0;
    //keep storing values from the text file so long as data exists:
    while (ifile >> num)
    {
        sw.push_back(num);
    }
    
    //verify that the scores were stored correctly:
    for (int i = 0; i < sw.size(); ++i)
    {
        //cout << sw[i] << endl;
    }
    if (sw.size()==pP->num_tp*2)
    {
        int count = 0;
        for (int i=0; i<pP->num_tp; i++)
        {
            Test_Point TP;
            tp.push_back(TP);
            for (int j=0; j<2; j++)
            {
                tp.at(i).point.push_back(sw.at(count));
                count += 1;
            }
        }
        assert(tp.size()==pP->num_tp);
    }
    
    
    if (sw.size()!=pP->num_tp*2)
    {
        for (int i=0; i<pP->num_tp; i++)
        {
            Test_Point TP;
            tp.push_back(TP);
            vector<double> r;
            for (int j=0; j<2; j++)
            {
                r.push_back(((double)rand()/RAND_MAX)*5);
            }
            tp.at(i).point = r;
            //hyper_val.push_back(r);
        }
        assert(tp.size()==pP->num_tp);
    }
    
    ofstream File13;
    File13.open("test_points.txt");
    for (int i=0; i<pP->num_tp; i++)
    {
        for (int j=0; j<2; j++)
        {
            File13 << tp.at(i).point.at(j) << "\t";
        }
        File13 << endl;
    }
    File13.close();
}


//-------------------------------------------------------------------------
//Checks how much of the hyper volume is dominated
void EA::Run_Hyper_Volume_Check_Quartet(Quartet *pQ)
{
    num_hyper_dom = 0;
    vector<vector<double> > PFront_ph;
    PFront_ph = pQ->pT->get_PFront();
    for (int i=0; i<PFront_ph.size(); i++)
    {
        vector<double> v1  = PFront_ph.at(i);
        //v1.at(0) = 0;
        //v1.at(1) = 0;
        for (int j=0; j<pP->num_tp; j++)
        {
            if(tp.at(j).dom == 0)
            {
                vector<double> v2 = tp.at(j).point;
                if (pQ->pT->does_v1_dominate_v2(v1, v2) == true)
                {
                    //hyper volume point is dominated by pareto front
                    //num_hyper_dom += 1;
                    tp.at(j).dom = 1;
                }
            }
        }
    }
    for (int j=0; j<pP->num_tp; j++)
    {
        if(tp.at(j).dom == 1)
        {
            num_hyper_dom += 1;
        }
    }
    assert (num_hyper_dom <= pP->num_tp);
    hyper_dom.push_back(num_hyper_dom);
}



//-------------------------------------------------------------------------
//Checks how much of the hyper volume is dominated
void EA::Run_Hyper_Volume_Check_PaCcET(PaCcET *pT)
{
    num_hyper_dom = 0;
    vector<vector<double> > PFront_ph;
    PFront_ph = pT->get_PFront();
    for (int i=0; i<PFront_ph.size(); i++)
    {
        vector<double> v1  = PFront_ph.at(i);
        //v1.at(0) = 0;
        //v1.at(1) = 0;
        for (int j=0; j<pP->num_tp; j++)
        {
            if(tp.at(j).dom == 0)
            {
                vector<double> v2 = tp.at(j).point;
                if (pT->does_v1_dominate_v2(v1, v2) == true)
                {
                    //hyper volume point is dominated by pareto front
                    //num_hyper_dom += 1;
                    tp.at(j).dom = 1;
                }
            }
        }
    }
    for (int j=0; j<pP->num_tp; j++)
    {
        if(tp.at(j).dom == 1)
        {
            num_hyper_dom += 1;
        }
    }
    assert (num_hyper_dom <= pP->num_tp);
    hyper_dom.push_back(num_hyper_dom);
}


//-------------------------------------------------------------------------
//Gets fitness for each individual
void EA::Get_PaCcET_Fitness(PaCcET *pT)
{
    //cout << "XXXX"  << endl;
    for (int a=0; a<pP->num_agents; a++)
    {
        indv.at(a).fitness = 0;
        Sub_Function(a);
        Function_0(a);
        Function_1(a);
        //Function_2(a);
        if (pP->use_PaCcet == 1)
        {
            pT->Pareto_Check(indv.at(a).F);
        }
    }
    
    //cout << "XXyX"  << endl;

    for (int a=0; a<pP->num_agents; a++)
    {
        if (pP->use_PaCcet == 1)
        {
            PaCcET_Fitness(pT, a);
        }
    }
       //cout << "XyyX"  << endl;
    if (pP->use_PaCcet==1)
    {
        Run_Hyper_Volume_Check_PaCcET(pT);
    }
}



//-------------------------------------------------------------------------
//Gets fitness for each individual
void EA::Get_Quartet_Fitness(Quartet *pQ)
{
    for (int a=0; a<pP->num_agents; a++)
    {
        indv.at(a).fitness = 0;
        Sub_Function(a);
        Function_0(a);
        Function_1(a);
        //Function_2(a);
    }
    
    
    for (int a=0; a<pP->num_agents; a++)
    {
        if (pP->use_quartet ==1)
        {
            Run_Quartet(pQ, a);
        }
    }
    if (pP->use_quartet ==1)
    {
        pQ->End_Generation();
    }
    if (pP->use_quartet==1)
    {
        Run_Hyper_Volume_Check_Quartet(pQ);
    }
}


//-------------------------------------------------------------------------
//Volumetric fitness
void EA::Volumetric_fitness(int a)
{
    indv.at(a).fitness = 0;
    indv.at(a).length = 0;
    indv.at(a).neg = 0;         //resets negative identifier
    for (int i=0; i<pP->num_F; i++)
    {
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
                indv.at(a).fitness = indv.at(a).fitness*abs((pP->set_point.at(i) - indv.at(a).F.at(i)));
            }
            else
            {
                 indv.at(a).fitness = indv.at(a).fitness*abs((pP->set_point.at(i) - indv.at(a).F.at(i)));
            }
        }
        indv.at(a).length += indv.at(a).F.at(i)*indv.at(a).F.at(i);
    }
    indv.at(a).length = sqrt(indv.at(a).length);
    assert (indv.at(a).length >= 1);
    //applies the negatvie volume if needed
    if (indv.at(a).neg == 1)
    {
        indv.at(a).fitness = indv.at(a).fitness*-1;
    }
}


//-------------------------------------------------------------------------
//Linear combintation with weights
void EA::Linear_Combination_Fitness(int a)
{
    indv.at(a).fitness = 0;
    indv.at(a).fitness = indv.at(a).F.at(0)*pP->W1 + indv.at(a).F.at(1)*pP->W2 + indv.at(a).F.at(2)*pP->W3;
    //swithces to a maximization problem
    //indv.at(a).fitness = 1/indv.at(a).fitness;
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
    if(indv.at(index_1).fitness < indv.at(index_2).fitness)
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
//Writes Pareto optimal points to a txt file
void EA::Write_Pareto_Optimal_Points_To_File()
{
    ofstream File3;
    File3.open("Pareto_Optimal_Points.txt");
    for (int p=0; p<point.size(); p++)
    {
        for (int i=0; i<pP->num_F; i++)
        {
            File3 << point.at(p).F.at(i) << "\t";
        }
        File3 << endl;
    }
    File3.close();
}


//-------------------------------------------------------------------------
//Stores all the f values that ever existed
void EA::Store_f_values(int gen)
{
    for (int i=0; i<pP->num_agents; i++)
    {
        F_value F;
        point.push_back(F);
        for (int v=0; v<pP->num_F; v++)
        {
            point.at(pP->num_agents*gen+i).F.push_back(indv.at(i).F.at(v));
        }
    }
}


//-------------------------------------------------------------------------
//Compares two points
void EA::Compare_Points(int p, int pp)
{
    if (point.at(p).F.at(0) > point.at(pp).F.at(0))
    {
        if (point.at(p).F.at(1) > point.at(pp).F.at(1))
        {
            if (point.at(p).F.at(2) > point.at(pp).F.at(2))
            {
                point.erase(point.begin() + pp);
            }
        }
    }
}


//-------------------------------------------------------------------------
//Finds the Pareto optimal points out of all the points ever
void EA::Find_Pareto_Optimal_Points()
{
    for (int p=0; p<point.size(); p++)
    {
        for (int pp=0; pp<point.size(); pp++)
        {
            if (p != pp)
            {
                Compare_Points(p, pp);
            }
        }
    }
}


//-------------------------------------------------------------------------
//Writes the quratet counter to a txt file
void EA::Write_Counter_File(Quartet* pQ)
{
    /*
    ofstream File11;
    File11.open("Quartet_Counter.txt");
    for (int i=0; i<pQ->quartet_counter.size(); i++)
    {
        for (int j=0; j<pQ->quartet_counter.at(i).size(); j++)
        {
            File11 << pQ->quartet_counter.at(i).at(j) << "\t";
        }
        File11 << endl;
    }
     */
    assert(pQ->Q_1.size()==pP->gen_max);
    assert(pQ->Q_2.size()==pP->gen_max);
    assert(pQ->Q_P.size()==pP->gen_max);
    
    ofstream File14;
    File14.open("Quartet_Counter_1.txt", ios_base::app);
    for (int i=0; i<pQ->Q_1.size(); i++)
    {
        File14 << pQ->Q_1.at(i) << "\t";
    }
    File14 << endl;
    
    ofstream File15;
    File15.open("Quartet_Counter_2.txt", ios_base::app);
    for (int i=0; i<pQ->Q_2.size(); i++)
    {
        File15 << pQ->Q_2.at(i) << "\t";
    }
    File15 << endl;
    
    ofstream File16;
    File16.open("Quartet_Counter_P.txt", ios_base::app);
    for (int i=0; i<pQ->Q_P.size(); i++)
    {
        File16 << pQ->Q_P.at(i) << "\t";
    }
    File16 << endl;
    
    pQ->quartet_counter.clear();
    pQ->Q_1.clear();
    pQ->Q_2.clear();
    pQ->Q_P.clear();
    
    //File11.close();
    File14.close();
    File15.close();
    File16.close();
    
    
}


//-------------------------------------------------------------------------
//Outputs the best individual info to the console
void EA::Output_Best_Individual_Info()
{
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
    cout << endl;
    cout << "//////////////////////////////////////////////" << endl;
}


//-------------------------------------------------------------------------
//Writes the amount of hyper volume points that were dominated each generation to a txt file
void EA::Write_Hyper_Dom_To_File()
{
    ofstream File12;
    File12.open("Hyper_Dom.txt", ios_base::app);
    assert(hyper_dom.size()==pP->gen_max);
    for (int i=0; i<hyper_dom.size(); i++)
    {
        File12 << hyper_dom.at(i) << "\t";
    }
    hyper_dom.clear();
    File12 << endl;
    File12.close();
}


//-------------------------------------------------------------------------
//Deltes the txt files
void EA::Delete_Files()
{
    if( remove( "Hyper_Dom.txt" ) != 0 )
        perror( "ERROR DELETING FILE Hyper_Dom" );
    else
        puts( "Hyper_Dom FILE SUCCEDDFULLY DELETED" );
    cout << endl;
    
    if( remove( "Quartet_Counter_1.txt" ) != 0 )
        perror( "ERROR DELETING FILE Quartet_Counter_1" );
    else
        puts( "Quartet_Counter_1 FILE SUCCEDDFULLY DELETED" );
    cout << endl;
    
    if( remove( "Quartet_Counter_2.txt" ) != 0 )
        perror( "ERROR DELETING FILE Quartet_Counter_2" );
    else
        puts( "Quartet_Counter_2 FILE SUCCEDDFULLY DELETED" );
    cout << endl;
    
    if( remove( "Quartet_Counter_P.txt" ) != 0 )
        perror( "ERROR DELETING FILE Quartet_Counter_P" );
    else
        puts( "Quartet_Counter_P FILE SUCCEDDFULLY DELETED" );
    cout << endl;
}


//-------------------------------------------------------------------------
//Runs entire multi-objective problem
void EA::Run_Multi_Objective()
{
    Delete_Files();
    Build_Hyper_Volume();
    for (int sr=0; sr<pP->num_sr; sr++)
    {
        for (int i=0; i<pP->num_tp; i++)
        {
            tp.at(i).dom = 0;
        }
        cout << endl;
        cout << "--------------------------------------------------------------------" << endl;
        PaCcET* pT;
        if (pP->use_PaCcet==1)
        {
            PaCcET T;
            pT = &T;
        }
        Quartet* pQ;
        if (pP->use_quartet==1)
        {
            Quartet Q;
            pQ = &Q;
        }
        Build_Pop();
        Create_set_point();
        for (int gen=0; gen<pP->gen_max; gen++)
        {
            if (gen < pP->gen_max-1)
            {
                if (gen%25 == 0)
                {
                    cout << sr << "::" << gen << endl;
                }
                if (pP->use_quartet==1)
                {
                    Get_Quartet_Fitness(pQ);
                }
                if (pP->use_PaCcet==1)
                {
                    Get_PaCcET_Fitness(pT);
                }
                
                //Output_Best_Individual_Info();
                
                Sort_indivduals_fitness();
                
                //Store_f_values(gen);
                Natural_Selection();
            }
            
            if (gen == pP->gen_max-1)
            {
                cout << sr << "::" << gen << endl;
                if (pP->use_quartet==1)
                {
                    Get_Quartet_Fitness(pQ);
                }
                if (pP->use_PaCcet==1)
                {
                    Get_PaCcET_Fitness(pT);
                }
                
                //Output_Best_Individual_Info();
                
                Sort_indivduals_fitness();
                
                //Store_f_values(gen);
                Write_final_pop_to_file();
                indv.clear();
            }
        }
        //Find_Pareto_Optimal_Points();
        //Write_Pareto_Optimal_Points_To_File();
        //T.exhaustive_to_file();
        //T.PFront_to_file();
        Write_Hyper_Dom_To_File();
        if (pP->use_quartet==1)
        {
            Write_Counter_File(pQ);
        }
        cout << "END STAT RUN" << endl;
    }
}

#endif /* EA_hpp */
