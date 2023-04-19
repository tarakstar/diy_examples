#include <boost/thread.hpp>
#include <boost/chrono.hpp>
#include <boost/atomic.hpp>
//#include <boost/bind/bind.hpp>  
#include <iostream>

using namespace std;

// g++ bthread_02.cpp -lboost_thread -lboost_system -lboost_chrono


class th_timer{
  int timeout_sec;
  public:
    //std::mutex m;
    boost::atomic_bool timer{true};
    boost::thread& th;
    //std::condition_variable flag;
    th_timer(int timeout,boost::thread& t):timeout_sec(timeout),th(t){}

    void operator()(){
       boost::this_thread::sleep_for(boost::chrono::seconds(timeout_sec));
       timer.store(false);
       cout<<"Timed out...."<<endl;
       th.interrupt();
       
    }

};

void ThreadFunction()
{
    int counter = 0;

    for(;;)
    {
        cout << "thread iteration " << ++counter << " Press Enter to stop" << endl;

        try
        {
            // Sleep and check for interrupt.
            // To check for interrupt without sleep,
            // use boost::this_thread::interruption_point()
            // which also throws boost::thread_interrupted
            boost::this_thread::interruption_point();
            boost::this_thread::sleep(boost::posix_time::milliseconds(500));
            boost::this_thread::interruption_point();
          }
        catch(boost::thread_interrupted&)
        {
            cout << "Thread is interrupted." << endl;
            return;
        }
    }
}

int main()
{

    std::cout << "Boost version: " << BOOST_LIB_VERSION << std::endl;
    // Start thread
    boost::thread t(&ThreadFunction);
    //t.join();

    th_timer gb(2,boost::ref(t));
    gb.timer.store(true);
    boost::thread timer(boost::ref(gb));


    // Join - wait when thread actually exits
    t.join();
    cout << "main: thread ended" << endl;

    return 0;
}
