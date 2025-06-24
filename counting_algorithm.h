#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <chrono>
#include <limits>
#include <numbers>
#include <iomanip>
#include <random>


bool is_valid(std::vector<long double> state, std::vector<long double> a, long double b);
class KnapsackMCMC
{
    private:
        std::vector<long double> m_a;
        long double m_b;
        std::vector<long double> m_state;
    public:
        long double current_sum = 0;
        KnapsackMCMC(std::vector <long double> a, long double b);
        void single_step();
        void run_chain(long long n_steps);
};
long double estimate_ratio(std::vector<long double> a, long double b_current, long double b_prev, long double epsilon, long double delta, long long mixing_steps);
std::vector<long double> fpras_knapscak(std::vector<long double> a, long double b, long double epsilon, long double delta);
long long exact_counting(std::vector<long double> a, long double b);