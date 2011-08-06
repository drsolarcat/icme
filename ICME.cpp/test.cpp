
#include <iostream>

using namespace std;

class Base {
  public:
    int f(int i, int k) {return 3*i;}
    int f(int i) {return i;}
};

class Derived : public Base {
  public:
    using Base::f;
    int f(int i, int k) {return 2*i;}
};

int main() {
  Derived obj;

  cout << obj.f(2) << endl;

  return 0;
}

