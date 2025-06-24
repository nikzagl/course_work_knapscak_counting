#include "counting_algorithm.h"
#include "number_gen.h"
#include <vector>
#include <windows.h>
int main()
{
 /*   
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
    */
   /*
   std::vector<long double> a = {1, 2, 3, 4, 5, 6, 7, 8, 9};
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
   std::vector<long double> a(16, 1);
   long long b = 1;
   long long exact_result = exact_counting(a, b);
   std::cout<<exact_result<<"\n";
   int n_successes = 0;
   std::vector<long double> fpras_result;
   for (int i = 0; i < 200; ++i)
   {
    fpras_result = fpras_knapscak(a, b, 0.1, 0.05);
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

  //std::string s = "sequences.txt";
  //generate_file(1000, 65536,s, 9, 200);
  /*
  std::ifstream f;
  f.open("sequences.txt");
  int N;
  int n_successes = 0;
  f >> N;
  for (int i = 0; i < N; ++i)
  {
    std::vector<long double> a;
    long double a_i, b;
    for (int j = 0; j < 9; ++j)
    {
       f >> a_i;
       a.push_back(a_i); 
    }
    f >> b;
    for (auto elem: a)
    {
        std::cout << elem <<" ";
    }
    std::cout<<"\n";
  std::cout<<b<<"\n";
long long exact_result = exact_counting(a, b);
std::cout<< exact_result<<"\n";
 std::vector<long double> fpras_result = fpras_knapscak(a, b, 0.1, 0.5);
 long double fpras_answer = fpras_result[0];
 long double error = fpras_result[1];
  if (abs(fpras_answer - exact_result) <= error)
    {
        n_successes += 1;
    }
    std::cout<<fpras_result[0] << " "<<fpras_result[1] <<"\n";
}
std::cout << long double(n_successes)/(long double)N;
*/
long double b;
std::vector<long double> a;
std::ifstream f;
f.open("p07_c.txt");
f >> b;
f.close();
f.clear();
f.open("p07_w.txt");
long double a_i;
while(f >> a_i)
{
    a.push_back(a_i);
}
long long exact_result = exact_counting(a, b);
double N_successes = 0;
for (int i = 0; i < 200; ++i)
{
    std::vector<long double> fpras_result = fpras_knapscak(a, b, 0.1, 0.2);
    long double fpras_value = fpras_result[0];
    long double error = fpras_result[1];
    std::cout<<fpras_value<<" "<<error<<"\n";
    std::cout<<exact_result<<"\n";
    if (abs(fpras_value - exact_result) <= error)
    {
        N_successes += 1;
    }
    std::cout << i << "...";
}
std::cout<<N_successes/200;
}