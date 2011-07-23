
#include <iostream>

using namespace std;

class A {
  public:
    int test() {return 2;};
};

class B : public A {
  //public:
    //int test() {return 1;}
};

int B::test() {return 1;}

int main() {
  B b;
  cout << b.test() << endl;
  return 0;
}

