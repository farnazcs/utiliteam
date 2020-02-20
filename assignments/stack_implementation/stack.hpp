#include <array>
#include <vector>
#include <iostream>

//#define Max_size 300000000
using namespace std;
class IntStack
{
    public:

        /// Returns true if there are no elements on the stack
        bool empty() const
        {
            if (size_amount==0)
                return true;
            else
                return false;
        };

        /// Returns the number of elements on the stack
        int size() const
        {
            return size_amount;

        };

        /// Return the top element of the stack, i.e. the last element that was pushed
        int top() const
        {if (size_amount==0)
            {cout<<"Error: stack empty\n"<<endl;
                return -1;
            }
            return num[size_amount-1];
        };

        /// Adds a new element to the top of the stack
        void push(const int& new_value)
        {
            if (size_amount<300000000)
                num[size_amount++] = new_value;
            else
                cout<<"Error: stack full\n"<<endl;
        };

        /// Remove the topmost element of the stack
        void pop()
        {
            if (size_amount == 0)
               cout<<"Error: stack empty"<<endl;
            else 
                size_amount--;

        };

    private:
        //int num[300000000];
//        std::array<int,300000000> num;
       vector<int> num = vector<int>(300000000);
        int size_amount = 0;
};
