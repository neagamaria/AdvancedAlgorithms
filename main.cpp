#include <iostream>
#include <fstream>
#include <vector>
#include <bits/stdc++.h>
#include <random>
#include <cstdlib>
#include <ctime>

using namespace std;

ofstream out("Evolutie.txt");

int n; //dimensiune pop.
int a, b; //intervalul
int c1, c2, c3; //coeficientii
int p; //precizia
double pc, pm; //probab. de crossover si mutatie
int m, l; //numar pasi
int iMax;
double Max;

vector<string> chrom, selectedChrom;
vector<double> x, f, q;

int find_l()
{
    int arg = (b-a) * pow(10, p);
    double log = log2(arg);
    int l = ceil(log);
    return l;
}

double get_u()
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1
    double u = dis(gen);

    return u;
}

double get_x(int l, string c)
{
    double val = (b - a) / (pow(2, l) - 1);
    double x = 0;

    for(int i=c.size()-1; i>=0; i--)
    {
        if(c[i] == '1')
            x += val;

        val *= 2;
    }
    return x + a;
}

double apply_f(double x)
{
    return c1*x*x + c2*x + c3;
}

string generate_chrom(int l)
{
    string chrom = "";
    for(int i=0; i<l; i++)
    {
        int bit = rand() % 2;
        chrom += to_string(bit);
    }
    return chrom;
}

int find_interval(double x, int st, int dr)
{
    while(st <= dr)
    {
        int m = (st + dr) / 2;
        double i1 = q[m];
        double i2 = q[m+1];

        if(x >= i1 && x <= i2)
            return m;

        if(x < i1)
            dr = m - 1;

        else
            st = m + 1;
    }

    return -1;
}

pair<string, string> cross(string chrom1, string chrom2, int i)
{
    string r1 = chrom1, r2 = chrom2;

    for(int j=i; j<chrom1.size(); j++)
    {
        r1[j] = chrom2[j];
        r2[j] = chrom1[j];
    }

    return make_pair(r1, r2);
}


///**********************

void initial_population(int step)
{
    l = find_l();

    for(int i=0; i<n; i++)
    {
        string currentChrom = generate_chrom(l);
        chrom[i] = currentChrom;
        x[i] = get_x(l, currentChrom);
        f[i] = apply_f(x[i]);
    }

    if(step == 1)
    {
        out<<"Populatia initiala:\n";
        for(int i=0; i<n; i++)
            out<<"b"<<i<<"="<<chrom[i]<<" x"<<i<<"="<<x[i]<<" f"<<i<<"="<<f[i]<<"\n";
    }


    for(int i=0; i<n; i++)
    {
        if(f[i] > Max)
        {
            Max = f[i];
            iMax = i;
        }
    }
}

void calculate_probabilities(int step)
{
    double F = 0;
    for(int i=0; i<n; i++)
        F += f[i];

    if(step == 1)
    {
        out<<"\nF="<<F<<"\n";
        out<<"\n\nProbabilitati selectie:\n";
    }


    double sum = 0;
    double p;

    for(int i=0; i<n; i++)
    {
        p = f[i] / F;

        if(step == 1)
            out<<"cromozom "<<i<<" p"<<i<<"="<<p<<"\n";

        q[i] = sum;
        sum += p;
    }
    q[n] = sum;

    if(step == 1)
    {
        out<<"\n\nIntervale probabilitati selectie :\n";

        for(int i=0; i<n; i++)
        {
            out<<"q"<<i<<"="<<q[i]<<" ";
        }
        out<<"q"<<n<<"="<<q[n]<<"\n";
    }
}

void selection_process(int step)
{
    random_device rd;
    unsigned seed = rd() ^ static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::mt19937 gen(seed);
    uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1

    if(step == 1)
        out<<"\n\nProces selectie\n";

    int number = 0;

    for (int i=0; i<n; i++)
    {
        double u = dis(gen);
        if(step == 1)
            out<<"u="<<u<<"\n";

        int intervalIndex = find_interval(u, 0, n);

        if(intervalIndex != -1)
        {
            if(step == 1)
                out<<"cromozom selectat pentru intervalul "<<intervalIndex<<": "<<chrom[intervalIndex]<<"\n";

            f[i] = f[intervalIndex];
            x[i] = x[intervalIndex];
            selectedChrom[number++] = chrom[intervalIndex];
        }

            else
                if(step == 1)
                    out<<-1<<"\n";
    }

    if(step == 1)
    {
        out<<"\n\nDupa selectie:\n";
        for(int i=0; i<selectedChrom.size(); i++)
            out<<selectedChrom[i]<<" x="<<x[i]<<" f="<<f[i]<<"\n";
    }
}


void crossover(int step)
{
    vector<int> crossoverIndexes;

    random_device rd;
    unsigned seed = rd() ^ static_cast<unsigned>(chrono::high_resolution_clock::now().time_since_epoch().count());
    std::mt19937 gen(seed);
    uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1

    srand(static_cast<unsigned int>(std::time(0)));

    int Max = 0;
    int iMax;

    if(step == 1)
        out<<"\n\nProbabilitatea de incrucisare: "<<pc<<"\n";

    for(int i=0; i<selectedChrom.size(); i++)
    {
        double u = dis(gen);

        if(step == 1)
            out<<selectedChrom[i]<<" u="<<u;

        if(u < pc)
        {
            crossoverIndexes.push_back(i);
            if(step == 1)
                out<<" participa\n";
        }

        else
            if(step == 1)
                out<<"\n";
    }

    int number = crossoverIndexes.size();
    int start = 0;
    if(number > 1)
    {
        if(number % 2 == 1)
        {
            string chrom1 = selectedChrom[crossoverIndexes[0]], chrom2 = selectedChrom[crossoverIndexes[1]], chrom3 = selectedChrom[crossoverIndexes[2]];
            int index  = rand() % (chrom1.size() - 1);

            if(step == 1)
            {
                out<<"\n\nIncrucisare "<<crossoverIndexes[0]<<" "<<crossoverIndexes[1]<<" "<<" "<<crossoverIndexes[2]<<"\n";
                out<<chrom1<<" "<<chrom2<<" "<<chrom3<<" index:"<<index;
            }

            pair<string, string> children1 = cross(chrom1, chrom2, index);
            pair<string, string> children2 = cross(children1.second, chrom3, index);

            selectedChrom[crossoverIndexes[0]] = children1.first;
            selectedChrom[crossoverIndexes[1]] = children1.second;
            selectedChrom[crossoverIndexes[2]] = children2.second;

            if(step == 1)
                out<<"\nDupa incrucisare: "<<children1.first<<" "<<children1.second<<" "<<children2.second<<"\n";

            start = 3;

        }

        for(int i=start; i<crossoverIndexes.size()-1; i+=2)
        {
            string chrom1 = selectedChrom[crossoverIndexes[i]];
            string chrom2 = selectedChrom[crossoverIndexes[i+1]];
            int index  = rand() % (chrom1.size() - 1);

            if(step == 1)
            {
                out<<"\n\nIncrucisare "<<crossoverIndexes[i]<<" "<<crossoverIndexes[i+1]<<" ";
                out<<chrom1<<" "<<chrom2<<" index:"<<index;
            }


            pair<string, string> children = cross(chrom1, chrom2, index);
            selectedChrom[crossoverIndexes[i]] = children.first;
            selectedChrom[crossoverIndexes[i+1]] = children.second;

            if(step == 1)
                out<<"\nDupa incrucisare: "<<children.first<<" "<<children.second<<"\n";

        }
    }
}

void mutation(int step)
{
    random_device rd;
    unsigned seed = rd() ^ static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::mt19937 gen(seed);
    uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1
    srand(static_cast<unsigned int>(std::time(0)));

    if(step == 1)
    {
        out<<"\n\nProbabilitate de mutatie: "<<pm<<"\n";
        out<<"\n\nAu fost modificati cromozomii:\n";
    }


    for(int i=0; i<selectedChrom.size(); i++)
    {
        double u = dis(gen);
        //out<<u<<" ";
        if(u < pm)
        {
            if(step == 1)
                out<<i<<" ";

            int pos = rand() % (selectedChrom[i].size() - 1);
            selectedChrom[i][pos] = 1 - (selectedChrom[i][pos] == '1');
        }
    }

    if(step == 1)
    {
        out<<"\n\nPopulatia dupa mutatii aleatoare:\n";
        for(int i=0; i<selectedChrom.size(); i++)
        {
            out<<selectedChrom[i]<<"\n";
        }
    }
}

void f_max_mean(double &fMax, double &fMean)
{
    selectedChrom.pop_back();
    selectedChrom.push_back(chrom[iMax]);
    fMax = 0;
    fMean = 0;
    double sum = 0;
    for(int i=0; i<selectedChrom.size(); i++)
    {
        double fMutation = apply_f(get_x(l, selectedChrom[i]));

        if(fMutation > fMax)
            fMax = fMutation;

        sum += fMutation;
    }

    fMean = sum / n;

}


int main()
{
    mt19937 rng(time(0));

    cin>>n>>a>>b>>c1>>c2>>c3>>p>>pc>>pm>>m;

    pc /= 100;
    pm /= 100;

    chrom.resize(n);
    f.resize(n);
    x.resize(n);
    selectedChrom.resize(n);
    q.resize(n);

    double fMax, fMean;
    int step = 1;

    initial_population(1);
    calculate_probabilities(1);
    selection_process(1);
    crossover(1);
    mutation(1);
    f_max_mean(fMax, fMean);

    do
    {

        out<<setprecision(16)<<fMax<<" "<<fMean<<"\n";
        step++;
        initial_population(step);
        calculate_probabilities(step);
        selection_process(step);
        crossover(step);
        mutation(step);
        f_max_mean(fMax, fMean);
    }
    while(step < m);

    return 0;
}
