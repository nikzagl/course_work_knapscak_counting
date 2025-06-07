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

std::independent_bits_engine<std::mt19937,64,uint_fast64_t> generator (0);


bool is_valid(std::vector<long double> state, std::vector<long double> a, long double b)
{
    long double dot_product = 0;
    for (int i = 0; i < a.size(); ++i)
    {
        dot_product += a[i]*state[i];
    }
    return dot_product <= b;
}
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

KnapsackMCMC::KnapsackMCMC(std::vector<long double> a, long double b)
{
    m_a = a;
    m_b = b;
    m_state = std::vector<long double>(a.size(), 0);
}
void KnapsackMCMC::single_step()
{
    int choice = generator()%2;
    if (!(choice))
    {
        return;
    }
    int i = generator() % (m_a.size());
    if (m_state[i] == 1)
    {
        m_state[i] = 0;
        current_sum -= m_a[i];
    }
    else if (current_sum + m_a[i] <= m_b)
    {
        m_state[i] = 1;
        current_sum += m_a[i];
    }
}
void KnapsackMCMC::run_chain(long long n_steps)
{
    for (long long i = 0; i < n_steps; ++i)
    {
        KnapsackMCMC::single_step();
    }
}

long double estimate_ratio(std::vector<long double> a, long double b_current, long double b_prev, long double epsilon, long double delta, long long mixing_steps)
{
     long long num_samples = size(a)*ceil(1/(pow(epsilon, 2))*log(2/delta));

     long long count = 0;
     for (long long i = 0; i <= num_samples; ++i)
     {
        KnapsackMCMC chain = KnapsackMCMC(a, b_current);
        chain.run_chain(mixing_steps);
        if (chain.current_sum <= b_prev)
        {
            count += 1;
        }
     }
     return (long double)num_samples/(long double)count;
}
std::vector<long double> fpras_knapscak(std::vector<long double> a, long double b, long double epsilon=0.2, long double delta=0.2)
{
    std::vector<long double> a_sorted = a;
    std::sort(a_sorted.begin(), a_sorted.end());
    std::vector<long double> b_sequence = {0};
    long double current_sum = 0;
    for (int i = 0; i < a.size(); ++i)
    {
        current_sum += a_sorted[i];
        if (current_sum > b)
            break;
        b_sequence.push_back(current_sum);
    }
    int k = b_sequence.size() - 1;
    if (k == 0)
        return {1, epsilon, delta};
    
    long long mixing_steps = 2*ceil(a.size() + log(1/epsilon));
    long long num_samples = size(a)*ceil(1/pow(epsilon/k, 2)*log(2/delta));
    long long count_omega = 0;
    
    for (int i = 0; i < num_samples; ++i)
    {
        KnapsackMCMC chain_omega = KnapsackMCMC(a_sorted, b);
        chain_omega.run_chain(mixing_steps);
        if (chain_omega.current_sum <= b_sequence[b_sequence.size()-1])
            count_omega += 1;
    }
    long double fpras_result = (long double)(num_samples)/(long double)(count_omega);
    for (int i = 1; i < k + 1; ++i)
    {
        long double ratio = estimate_ratio(a_sorted, b_sequence[i], b_sequence[i - 1], epsilon/k, delta, mixing_steps);
        fpras_result *= ratio;
    }
    std::vector<long double> result;
    result.push_back(fpras_result);
    result.push_back(fpras_result*epsilon);
    result.push_back(delta);
    return result;
}
long long exact_counting(std::vector<long double> a, long double b)
{
    long long result = 0;
    for(long long i = 0; i < (long long)ceil(pow(2, a.size())); ++i)
    {
        std::vector<long double> m_state;
        long long i_t = i;
        for (int j = 0; j < a.size(); ++j)
        {
            m_state.push_back(i_t%2);
            i_t = i_t/2;
        }
        if (is_valid(m_state, a, b))
        {
            result += 1;
        }
    }
    return result;
}
int main()
{
    std::ofstream file;
    file.open("dep_times_n.txt");
    std::vector<long double> a = {1, 2};
    long double b = 2;
    for (int i = 3; i <= 22; ++i)
    {
        auto beg = std::chrono::high_resolution_clock::now();
        fpras_knapscak(a, b, 0.1, 0.1);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - beg);
        std::cout<<a.size()<<" "<<duration.count()<<"\n";
        file << a.size() << " "<<duration.count() << "\n";
        b = i;
        a.push_back(i);
    }
    file.close();
    file.open("dep_times_eps.txt");
    a = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    b = 12;
    for (long double eps = 0.5; eps < 0.98; eps += 0.02)
    {
         auto beg = std::chrono::high_resolution_clock::now();
        fpras_knapscak(a, b, 1-eps, 0.1);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - beg);
        std::cout<<eps<<" "<<duration.count()<<"\n";
        file <<eps << " "<<duration.count() << "\n";
    }
    file.close();
    file.open("dep_times_delta.txt");
    for (long double delta = 0.8;delta < 0.99; delta += 0.01)
    {
         auto beg = std::chrono::high_resolution_clock::now();
        fpras_knapscak(a, b, 0.1, 1-delta);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - beg);
        std::cout<<delta<<" "<<duration.count()<<"\n";
        file << delta<< " "<<duration.count() << "\n";
    }
    file.close();
   /*std::vector<long double> a = {1, 2, 3, 4, 5, 6, 7, 8, 9};
   long double b = 10;
   long long exact_result = exact_counting(a, b);
   int n_successes = 0;
   std::vector<long double> fpras_result;
   for (int i = 0; i < 200; ++i)
   {
    std::cout << i + 1;
    fpras_result = fpras_knapscak(a, b, 0.05, 0.05);
    long double fpras_answer = fpras_result[0];
    long double error = fpras_result[1];
    if (abs(fpras_answer - exact_result) <= error)
    {
        n_successes += 1;
    }
   }
   long double ratio = (long double)n_successes/200;
   std::cout<<"\n"<< ratio<<"\n"<<n_successes;
   */
   /*
   std::vector<long double> a(16,1);
   long double b = 1;
   long long exact_result = exact_counting(a, b);
   std::cout<<exact_result<<"\n";
   int n_successes = 0;
   std::vector<long double> fpras_result;
   for (int i = 0; i < 200; ++i)
   {
    fpras_result = fpras_knapscak(a, b, 0.05, 0.05);
    long double fpras_answer = fpras_result[0];
    long double error = fpras_result[1];
    if (abs(fpras_answer - exact_result) <= error)
    {
        n_successes += 1;
    }
    std::cout<<fpras_result[0] << " "<<fpras_result[1] <<"\n";
    std::cout<<i<<"...";
   }
   long double ratio = (long double)n_successes/200;
   std::cout<<"\n"<< ratio<<"\n"<<n_successes;
   */
}