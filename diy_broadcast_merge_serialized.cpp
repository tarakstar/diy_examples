//
// Merge reduction merges blocks together, computing a sum of their values. At
// each round, one block of a group of k blocks is the root of the group. The
// other blocks send their data to the root, which computes the sum, and the
// root block (only) proceeds to the next round. After log_k(numblocks) rounds,
// one block contains the global sum of the values.
//

// mpic++ diy_broadcast_merge.cpp -I${DIY_INC} -lpthread
// mpirun -np 2 ./a.out -d 1 -b 2 -v true

#include <cmath>
#include <vector>

#include <diy/master.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>
#include <diy/partners/broadcast.hpp>
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
      TestData(int a) : x(a) {y=x*x;}
      void print(){
        std::cout<<"TestData "<<x<<"\t"<<y<<std::endl;
      }
};

namespace diy
{
  template<>
  struct Serialization<TestData>
  {
    static void save(BinaryBuffer& bb, const TestData& p)       { diy::save(bb, p); }
    static void load(BinaryBuffer& bb, TestData& p)             { diy::load(bb, p); }
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

    void            generate_data(size_t n)                    // initialize block values
    {
        data.resize(n);
        for (size_t i = 0; i < n; ++i)
            data[i] = i;
    }
    // block data
    Bounds          bounds;
    vector<int>     data;
    TestData td;

private:
    Block()                                     {}
};

// diy::decompose needs to have a function defined to create a block
// here, it is wrapped in an object to add blocks with an overloaded () operator
// it could have also been written as a standalone function
struct AddBlock
{
    AddBlock(diy::Master& master_, size_t num_points_):
        master(master_),
        num_points(num_points_)
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

            b->generate_data(num_points);          // initialize block data (typical)
        }

    diy::Master&  master;
    size_t        num_points;
};

// --- callback functions ---//

//
// callback function for merge operator, called in each round of the reduction
// one block is the root of the group
// link is the neighborhood of blocks in the group
// root block of the group receives data from other blocks in the group and reduces the data
// nonroot blocks send data to the root
//
void sum(Block* b,                                  // local block
         const diy::ReduceProxy& rp,                // communication proxy
         const diy::RegularMergePartners& partners) // partners of the current block
{
    unsigned   round    = rp.round();               // current round number

    // step 1: dequeue and merge
    for (int i = 0; i < rp.in_link().size(); ++i)
    {
        int nbr_gid = rp.in_link().target(i).gid;
        if (nbr_gid == rp.gid())
        {
            fmt::print(stderr, "[{}:{}] Skipping receiving from self\n", rp.gid(), round);
            continue;
        }

        std::vector<int>    in_vals;
        rp.dequeue(nbr_gid, in_vals);
        fmt::print(stderr, "[{}:{}] Received {} values from [{}]\n",
                   rp.gid(), round, (int)in_vals.size(), nbr_gid);
        for (size_t j = 0; j < in_vals.size(); ++j)
            (b->data)[j] += in_vals[j];
    }

    // step 2: enqueue
    for (int i = 0; i < rp.out_link().size(); ++i)    // redundant since size should equal to 1
    {
        // only send to root of group, but not self
        if (rp.out_link().target(i).gid != rp.gid())
        {
            rp.enqueue(rp.out_link().target(i), b->data);
            fmt::print(stderr, "[{}:{}] Sent {} valuess to [{}]\n",
                       rp.gid(), round, (int)b->data.size(), rp.out_link().target(i).gid);
        } else
            fmt::print(stderr, "[{}:{}] Skipping sending to self\n", rp.gid(), round);

    }
}

//////////////////////////////////////////////////////////
// This function broadcasts a vector<int> to all blocks and also receive it.
void assign_data(Block* b,                                  // local block
         const diy::ReduceProxy& rp,                // communication proxy
         const diy::RegularBroadcastPartners& partners) // partners of the current block
{
    unsigned   round    = rp.round();               // current round number
    //cout<<"round = "<<round<<endl;
    std::vector<int> mydata (10,1); // data to be broadcasted
    std::vector<TestData> mydata2(10);
    //int mydata=5;
  
    for(int i=0;i<10;i++){
      TestData td(i);
      mydata2.push_back(td);
    }

    diy::MemoryBuffer bb;
    save(bb,&mydata2);
    
    cout<<"Broadcasting data...."<<rp.gid()<<endl;
    cout<<"in out : "<<rp.in_link().size()<<"\t"<<rp.out_link().size()<<endl;

    // send data to others but not to the self
    for (int i = 0; i < rp.out_link().size(); ++i)    // redundant since size should equal to 1
    {
        // send to everyone without the self
        //if (rp.out_link().target(i).gid != rp.gid())
        //{
            //rp.enqueue(rp.out_link().target(i), mydata);
            //rp.enqueue(rp.out_link().target(i), mydata2);
            rp.enqueue(rp.out_link().target(i), bb);
            //rp.enqueue(rp.out_link().target(i), bb);
            //fmt::print(stderr, "[{}:{}] Sent {} valuess to [{}]\n",
            //           rp.gid(), round, (int)mydata.size(), rp.out_link().target(i).gid);
        //} else
        //    fmt::print(stderr, "[{}:{}] Skipping sending to self\n", rp.gid(), round);
    }

    //diy::MemoryBuffer bb2;
    //std::vector<TestData>*   datar;
    //load(bb2,&mydatar);

    // now receive data which was sent and set block's vector<int> values from it.
    //cout<<"Receiving bcast..."<<rp.gid()<<endl;
    for (int i = 0; i < rp.in_link().size(); ++i)
    {
        int nbr_gid = rp.in_link().target(i).gid;

        std::vector<int>    in_vals;
        std::vector<TestData>*   datar;

        diy::MemoryBuffer bb2;
        //int in_vals;
        //cout<<"dequeue start "<<nbr_gid<<endl;
        //rp.dequeue(nbr_gid, in_vals);
        rp.dequeue(nbr_gid, bb2);
        load(bb2,datar);
        std::cout << "Position: " << bb2.position << std::endl;
        //cout<<"dequeue success "<<endl;
        std::cout<<datar->size()<<std::endl;
        //fmt::print(stderr, "[{}:{}] Received {} values from [{}]\n",
        //           rp.gid(), round, (int)in_vals.size(), nbr_gid);
        b->td = datar->at(nbr_gid); // pick-up msg for at current block gid
        for (size_t j = 0; j < datar->size(); ++j)
            //(b->data)[j] = in_vals[j];
            datar->at(i).print();
    }
   

    cout<<"Broadcast complete."<<endl;
}
//////////////////////////////////////////////////////////

//
// prints the block values
//
void print_block(Block* b,                             // local block
                 const diy::Master::ProxyWithLink& cp, // communication proxy
                 bool verbose)                         // user-defined additional arguments
{
    /*fmt::print(stderr, "[{}] Bounds: {} {} {} -- {} {} {}\n",
               cp.gid(),
               b->bounds.min[0], b->bounds.min[1], b->bounds.min[2],
               b->bounds.max[0], b->bounds.max[1], b->bounds.max[2]);*/

    //if (verbose && cp.gid() == 0)
    if (verbose)
    {
        /*fmt::print(stderr, "[{}] {} vals: ", cp.gid(), b->data.size());
        for (size_t i = 0; i < b->data.size(); ++i)
            fmt::print(stderr, "{}  ", b->data[i]);
        fmt::print(stderr, "\n");*/
        std::cout<<cp.gid()<<"\t";
        b->td.print();
    }
}

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
    AddBlock               create(master, num_points); // an object for adding new blocks to master

    // set some global data bounds
    Bounds domain;
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

    master.foreach([verbose](Block* b, const diy::Master::ProxyWithLink& cp)
             { print_block(b, cp, verbose); });  // callback function for each local block


    // broadcast with RegularBroadcastPartners
    diy::reduce(master,                              // Master object
                assigner,                            // Assigner object
                bpartners,                            // RegularMergePartners object
                &assign_data);                               // merge operator callback function
   

    master.foreach([verbose](Block* b, const diy::Master::ProxyWithLink& cp)
             { print_block(b, cp, verbose); });  // callback function for each local block


   // reduction
   /*diy::reduce(master,                              // Master object
                assigner,                            // Assigner object
                partners,                            // RegularMergePartners object
                &sum);                               // merge operator callback function
   

    master.foreach([verbose](Block* b, const diy::Master::ProxyWithLink& cp)
             { print_block(b, cp, verbose); });  // callback function for each local block
    */

   
}
