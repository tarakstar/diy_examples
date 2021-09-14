#include <iostream>
#include <diy/serialization.hpp>

class Test{
  int a,b;
  public:
  Test(int a_,int b_) : a(a_), b(b_) {};
  void print()
  { std::cout<<a<<"\t"<<b<<std::endl;}

};


// Point does not need any special handling for serialization; it will be saved
// as plain binary data
struct Point
{
  int&  operator[](int i)       { return coords[i]; }
  int   coords[3];
};

struct PointVec
{
        PointVec()              {}          // need this for de-serialization
        PointVec(size_t dim) :
          coords(dim)           {}
  int&  operator[](size_t i)       { return coords[i]; }
  std::vector<int>  coords;
};

namespace diy
{
  template<>
  struct Serialization<PointVec>
  {
    static void save(BinaryBuffer& bb, const PointVec& p)       { diy::save(bb, p.coords); }
    static void load(BinaryBuffer& bb, PointVec& p)             { diy::load(bb, p.coords); }
  };
}

namespace diy
{
  template<>
  struct Serialization<Test>
  {
    static void save(BinaryBuffer& bb, const Test& p)       { diy::save(bb, p); }
    static void load(BinaryBuffer& bb, Test& p)             { diy::load(bb, p); }
  };
}


int main()
{
  diy::MemoryBuffer bb;
  std::cout << "Position: " << bb.position << std::endl;

  {
    std::vector<Point>    points;
    Point p; p[0] = p[1] = p[2] = 5;
    points.push_back(p);
    p[0] = 1;
    points.push_back(p);
    p[2] = 1;
    points.push_back(p);

    save(bb, points);
  }
  std::cout << "Position: " << bb.position << std::endl;

  {
    std::vector<PointVec> points;
    PointVec p(3);
    points.push_back(p);
    p[0] = 2;
    points.push_back(p);
    p[2] = 2;
    points.push_back(p);

    save(bb, points);
  }
  std::cout << "Position: " << bb.position << std::endl;

  bb.reset();
  std::cout << "Position: " << bb.position << std::endl;

  {
    std::vector<Point>    points;
    load(bb, points);
    for (unsigned i = 0; i < points.size(); ++i)
    {
      for (unsigned j = 0; j < 3; ++j)
        std::cout << points[i][j] << ' ';
      std::cout << std::endl;
    }
  }
  std::cout << "Position: " << bb.position << std::endl;

  {
    std::vector<PointVec> points;
    load(bb, points);
    for (unsigned i = 0; i < points.size(); ++i)
    {
      for (unsigned j = 0; j < 3; ++j)
        std::cout << points[i][j] << ' ';
      std::cout << std::endl;
    }
  }
  std::cout << "Position: " << bb.position << std::endl;


  //--------- Serialize Test class-----------------------

  diy::MemoryBuffer bb1;
  Test t1(1,2);
  save(bb1, &t1);

  Test t2(2,3);
  save(bb1,&t2);

  bb1.reset();

  Test* r1;
  load(bb1,r1);
  r1->print();

   load(bb1,r1);
  r1->print();

 



}
