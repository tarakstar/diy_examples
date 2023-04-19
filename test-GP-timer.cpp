
// mpic++ test-GP-timer.cpp -I${DIY_INC} -lpthread
// mpirun -np 2 ./a.out -d 1 -b 4 -v false

#include <cmath>
#include <vector>
#include <thread>
#include <chrono>
#include <mutex>
#include <condition_variable>
#include <atomic>


#include <diy/master.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>
#include <diy/partners/broadcast.hpp>
#include <diy/partners/all-reduce.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/master.hpp>
#include <diy/serialization.hpp>

#include "opts.h"


using namespace std;

typedef     diy::ContinuousBounds       Bounds;
typedef     diy::RegularContinuousLink  RCLink;

// A class to test serialized data exchange
class TestData{
  public:
    int x=0;
    int y=0;
    TestData() {}
    TestData(int a) : x(a){}
    void compute() {y = x*x;}
    void print(){
      std::cout<<"TestData "<<x<<"\t"<<y<<std::endl;
    }
};

class GlobalData{
  
  vector<int> data;
  public:
    GlobalData() {}

    void add(int i){
      data.push_back(i);
    }

    void print(){
      cout<<"data : ";
      for(auto i:data)
        cout<<i<<"\t";
      cout<<endl;
    }

};

class myglobal{
  int timeout_sec;
  public:
    std::mutex m;
    std::atomic_bool timer;
    std::condition_variable flag;
    myglobal(int t):timeout_sec(t){}

    void operator()(){
       std::this_thread::sleep_for(std::chrono::seconds(timeout_sec));
       timer.store(false);
    }

};


namespace diy
{
  template<>
    struct Serialization<TestData>
    {
      static void save(BinaryBuffer& bb, const TestData& p)       
      { diy::save(bb, p.x); 
        diy::save(bb, p.y);
      }
      static void load(BinaryBuffer& bb, TestData& p)   
      { diy::load(bb, p.x); 
        diy::load(bb, p.y);
      }
    };
}



// --- block structure ---//

// the contents of a block are completely user-defined
// however, a block must have functions defined to create, destroy, save, and load it
// create and destroy allocate and free the block, while
// save and load serialize and deserialize the block
// these four functions are called when blocks are cycled in- and out-of-core
// they can be member functions of the block, as below, or separate standalone functions
struct Block
{
  Block(const Bounds& bounds_):
    bounds(bounds_)                         {}

  static void*    create()                                   // allocate a new block
  { return new Block; }
  static void     destroy(void* b)                           // free a block
  { delete static_cast<Block*>(b); }
  static void     save(const void* b, diy::BinaryBuffer& bb) // serialize the block and write it
  {
    diy::save(bb, static_cast<const Block*>(b)->bounds);
    diy::save(bb, static_cast<const Block*>(b)->data);
  }
  static void     load(void* b, diy::BinaryBuffer& bb)       // read the block and deserialize it
  {
    diy::load(bb, static_cast<Block*>(b)->bounds);
    diy::load(bb, static_cast<Block*>(b)->data);
  }

    // block data
  Bounds          bounds{1};
  int     data;
  int myresult;
  vector<TestData> myresult2;
  TestData td;

  private:
  Block() {}
};

// diy::decompose needs to have a function defined to create a block
// here, it is wrapped in an object to add blocks with an overloaded () operator
// it could have also been written as a standalone function
struct AddBlock
{
  AddBlock(diy::Master& master_):
    master(master_)
  {}

  // this is the function that is needed for diy::decompose
  void  operator()(int gid,                // block global id
      const Bounds& core,     // block bounds without any ghost added
      const Bounds& bounds,   // block bounds including any ghost region added
      const Bounds& domain,   // global data bounds
      const RCLink& link)     // neighborhood
    const
    {
      Block*          b   = new Block(core);
      RCLink*         l   = new RCLink(link);
      diy::Master&    m   = const_cast<diy::Master&>(master);

      m.add(gid, b, l); // add block to the master (mandatory)

    }

  diy::Master&  master;
};




//////////////////////////////////////////////////////////

// Perform aggregate calculations
// In this example : form a vector with results of TestData
void collect_results(Block* b,                                  // local block
    const diy::ReduceProxy& rp,                // communication proxy
    const diy::RegularMergePartners& partners)// partners of the current block 

{
  unsigned   round    = rp.round();               // current round number

  //cout<<"rplink size "<<rp.in_link().size()<<"\t"<<rp.out_link().size()<<endl;

  // step 1: dequeue and merge
  for (int i = 0; i < rp.in_link().size(); ++i)
  {
    int nbr_gid = rp.in_link().target(i).gid;
    if (nbr_gid == rp.gid())
     continue;

    vector<TestData> vtd;
    rp.dequeue(nbr_gid,vtd);
    for(auto &td:vtd){
      //td.print();
      b->myresult2.push_back(td);
    }

  }

  // step 2: enqueue
  for (int i = 0; i < rp.out_link().size(); ++i)    // redundant since size should equal to 1
  {

    // only send to root of group, but not self
    if (rp.out_link().target(i).gid != rp.gid())
    {
      //cout<<"round out gid "<<rp.round()<<"\t"<<rp.gid()<<"\t"<<rp.out_link().target(i).gid<<endl;
      rp.enqueue(rp.out_link().target(i), b->myresult2);

    }
  }

}


//////////////////////////////////////////////////////////

void run_serial_code(Block* b,                             // local block
    const diy::Master::ProxyWithLink& cp, // communication proxy
    bool verbose)                         //
{

  //cout<<"gid b->td.y myresult "<<cp.gid()<<"\t"<<b->td.y<<"\t"<<b->myresult<<endl;
  std::chrono::seconds one_sec(1);
  //std::chrono::time_point starting_time = std::chrono::steady_clock::now();


  if(cp.gid()==0){
    cout<<"-----------------------------"<<endl;
    cout<<"The root block..."<<endl;
    cout<<"This should be executed only once for serial computing."<<endl;
    cout<<"myresult = "<<b->myresult<<endl;
    //cout<<"b->td.y "<<b->td.y<<endl;
    cout<<"Result vector size "<<b->myresult2.size()<<endl;
    for(auto &td:b->myresult2)
      td.print();
    cout<<"-----------------------------"<<endl;
    //std::this_thread::sleep_until(starting_time + one_sec);
    std::this_thread::sleep_for(std::chrono::milliseconds(2000));
    cout<<"Woke up...."<<endl;
  }

}


//////////////////////////////////////////////////////////
// needs to be called only by the root block
void initialize_GP(Block* b,                                  // local block
    const diy::ReduceProxy& rp,                // communication proxy
    const diy::RegularBroadcastPartners& partners, // partners of the current block
    int iteration) 
{

  

}
//////////////////////////////////////////////////////////
// This function involves computation on block data
void run_computations(Block* b,                             // local block
    const diy::Master::ProxyWithLink& cp, // communication proxy
    bool verbose)                         //
{

  b->td.compute();
  b->myresult = b->td.y;
  b->myresult2.push_back(b->td);


}

//////////////////////////////////////////////////////////
// This function broadcasts a vector<int> to all blocks and also receive it.
void assign_data(Block* b,                                  // local block
    const diy::ReduceProxy& rp,                // communication proxy
    const diy::RegularBroadcastPartners& partners, // partners of the current block
    int iteration) 
{
  unsigned   round    = rp.round();               // current round number
  //cout<<"round = "<<round<<endl;
  std::vector<int> mydata; // data to be broadcasted
  std::vector<TestData> mydata2;
  TestData tdobj(15);

  for(int i=0;i<10;i++){
    //TestData td(i);
    TestData td(i+iteration);
    mydata2.push_back(td);
    mydata.push_back(i);
  }

  //cout<<"Broadcasting data...."<<rp.gid()<<endl;
  //cout<<"in out : "<<rp.in_link().size()<<"\t"<<rp.out_link().size()<<endl;

  // send data to others but not to the self
  for (int i = 0; i < rp.out_link().size(); ++i)  
  {
    // The following lines try eneueing different data types
    //rp.enqueue(rp.out_link().target(i), mydata);
    rp.enqueue(rp.out_link().target(i), mydata2);
    //rp.enqueue(rp.out_link().target(i), tdobj);
  }

  // now receive data which was sent and set block's vector<int> values from it.
  //cout<<"Receiving bcast..."<<rp.gid()<<endl;
  for (int i = 0; i < rp.in_link().size(); ++i)
  {
    int nbr_gid = rp.in_link().target(i).gid;

    std::vector<TestData>   datar;
    TestData tdobjr;

    //rp.dequeue(nbr_gid, in_vals); // vector<int>
    rp.dequeue(nbr_gid, datar); // vector<TestData>
    //rp.dequeue(nbr_gid, tdobjr); // TestData

    b->td = datar.at(rp.gid()); // pick-up msg for at current block gid
    //b->td =  tdobjr;
    //b->mynum = in_vals[rp.gid()];

  }


  //cout<<"Broadcast complete."<<endl;
}
//////////////////////////////////////////////////////////

//
// prints the block values
//
void print_block(Block* b,                             // local block
    const diy::Master::ProxyWithLink& cp, // communication proxy
    bool verbose,int iteration)                         // user-defined additional arguments
{

  //if (verbose && cp.gid() == 0)
  if (verbose)
  {
    //cout<<"myresult : "<<b->myresult<<endl;
	cout<<"iteration = "<<iteration<<"\t"<<cp.gid()<<endl;
    //std::cout<<cp.gid()<<"\t";
    //b->td.print();
    //cout<<"Result = "<<b->myresult<<endl;
  }
}

//////////////////////////////////////////////////////////

// --- main program ---//

int main(int argc, char* argv[])
{
  diy::mpi::environment     env(argc, argv); // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
  diy::mpi::communicator    world;           // equivalent of MPI_COMM_WORLD

  int                       nblocks     = world.size(); // global number of blocks
  // in this example, one per process
  size_t                    num_points  = 10;           // points per block
  int                       mem_blocks  = -1;           // all blocks in memory
  int                       threads     = 1;            // no multithreading
  int                       dim         = 3;            // 3-d blocks

  // get command line arguments
  using namespace opts;
  Options ops;

  bool verbose, contiguous, help;
  ops
    >> Option('v', "verbose",    verbose,    "verbose output")
    >> Option('c', "contiguous", contiguous, "use contiguous partners")
    >> Option('h', "help",       help,       "show help")
    ;

  ops
    >> Option('d', "dim",     dim,            "dimension")
    >> Option('b', "blocks",  nblocks,        "number of blocks")
    >> Option('t', "thread",  threads,        "number of threads")
    ;

  if (!ops.parse(argc, argv) || help)
  {
    std::cout << "Usage: " << argv[0] << " [OPTIONS]\n";
    std::cout << ops;
    return 1;
  }

  // diy initialization
  diy::FileStorage       storage("./DIY.XXXXXX"); // used for blocks moved out of core
  diy::Master            master(world,            // master is the top-level diy object
      threads,
      mem_blocks,
      &Block::create,
      &Block::destroy,
      &storage,
      &Block::save,
      &Block::load);
  AddBlock               create(master); // an object for adding new blocks to master

  // set some global data bounds
  Bounds domain {dim};
  for (int i = 0; i < dim; ++i)
  {
    domain.min[i] = 0;
    domain.max[i] = 128.;
  }

  // choice of contiguous or round robin assigner
  diy::ContiguousAssigner   assigner(world.size(), nblocks);
  //diy::RoundRobinAssigner   assigner(world.size(), nblocks);

  // decompose the domain into blocks
  diy::RegularDecomposer<Bounds> decomposer(dim, domain, nblocks);
  decomposer.decompose(world.rank(), assigner, create);

  // merge-based reduction: create the partners that determine how groups are formed
  // in each round and then execute the reduction

  int k = 2;                               // the radix of the k-ary reduction tree

  // partners for merge over regular block grid
  diy::RegularMergePartners  partners(decomposer,  // domain decomposition
      k,           // radix of k-ary reduction
      contiguous); // contiguous = true: distance doubling
  // contiguous = false: distance halving

  diy::RegularBroadcastPartners bpartners(decomposer,k,contiguous);

  GlobalData gb;


  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //cout<<"Rank = "<<myrank<<endl;

  for(int it=0;it<2;++it){


	  cout<<"iteration : "<<it<<"\t"<<myrank<<endl;

	  // broadcast with RegularBroadcastPartners
	  // Assign points on which to perform computations
	  /*diy::reduce(master,                              // Master object
		  assigner,                            // Assigner object
		  bpartners,                            // RegularMergePartners object
		  &assign_data);                               // merge operator callback function
		*/


	  diy::reduce(master,                              // Master object
		  assigner,                            // Assigner object
		  bpartners,                            // RegularMergePartners object
		  [&] (Block* b, const diy::ReduceProxy& rp, const diy::RegularBroadcastPartners& partners )  
		  { assign_data (b,rp,bpartners,(myrank+1)*100+(it+1)*10); });                               // merge operator callback function
	  

	  //cout<<"Now we have TestData."<<endl;
	  master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
		  { print_block(b, cp, verbose,(myrank+1)*100+(it+1)*10); });  // callback function for each local block
	  

	  // Run computations in blocks
	  master.foreach([verbose](Block* b, const diy::Master::ProxyWithLink& cp)
		  { run_computations(b, cp, verbose); });  // callback function for each local block

	  int rounds = partners.rounds();
	  //cout<<"round = "<<rounds<<endl;

	 // reduction
	  diy::reduce(master,                              // Master object
		  assigner,                            // Assigner object
		  partners,                            // RegularMergePartners object
		  [&](Block* b, const diy::ReduceProxy& rp, const diy::RegularMergePartners& partners )
		  { collect_results(b, rp, partners); });


	  master.foreach([verbose](Block* b, const diy::Master::ProxyWithLink& cp)
	  	  { run_serial_code(b, cp, verbose); });  

	  cout<<"Printing block TestData contents and results after computations."<<endl;
	  master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
		  { print_block(b, cp, verbose,it); });  // callback function for each local block
	  


  }

}
