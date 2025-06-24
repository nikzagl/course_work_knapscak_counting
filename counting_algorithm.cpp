
 #include "counting_algorithm.h"

std::random_device dev;
std::mt19937 rng(dev());
std::uniform_int_distribution<uint32_t> uint_2(0, 1);         // by default range [0, MAX]



bool is_valid(std::vector<long double> state, std::vector<long double> a, long double b)
{
    long double dot_product = 0;
    for (int i = 0; i < a.size(); ++i)
    {
        dot_product += a[i]*state[i];
    }
    return dot_product <= b;
}

KnapsackMCMC::KnapsackMCMC(std::vector<long double> a, long double b)
{
    m_a = a;
    m_b = b;
    m_state = std::vector<long double>(a.size(), 0);
}
void KnapsackMCMC::single_step()
{
    int choice = uint_2(rng);
    std::uniform_int_distribution<uint32_t> uint_dist_masize(0, m_a.size() - 1);
    int i = uint_dist_masize(rng);
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
     long long num_samples = a.size()*ceil(1/(pow(epsilon, 2))*log(2/delta));

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
std::vector<long double> fpras_knapscak(std::vector<long double> a, long double b, long double epsilon, long double delta)
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
    
    long long mixing_steps = pow(a.size(), 2)*log(a.size()/epsilon);
    long long num_samples = a.size()*ceil(1/pow(epsilon/k, 2)*log(2/delta));
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
