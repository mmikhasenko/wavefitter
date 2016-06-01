
int f(double *x) {
  return x[0]+x[1]+x[2];
}

#include <vector>
#include <iostream>
int main() {
  
  double factor = 5;
  auto func0 = [](double *x) -> int { return f(x); };
  auto func1 = [&factor](double *x) -> int { return factor*f(x); };

  auto fl2 = func1;
  int (*lf)(double *x) = func0;
  
  double arr[3] = {1,2,3};
  
  std::vector<double> v1; v1.push_back(0.1); v1.push_back(3.); v1.push_back(4.);
  std::vector<double> v2(1);
  std::cout << v1.size() << std::endl;
  v1 = v2;
  std::cout << v1.size() << std::endl;

  return lf(arr);
}


