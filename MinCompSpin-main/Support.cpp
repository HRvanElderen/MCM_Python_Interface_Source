#include <sstream>
#include <algorithm>

using namespace std;

// create a string of binary values from an int.
string int_to_bstring(unsigned int x, unsigned int n)
{
    string s;
    do
    {
        s.push_back('0' + (x & 1));
    } while (x >>= 1);
    reverse(s.begin(), s.end());
    s = string(n - s.length(), '0') + s;

    return s;
}

// Counts all the set bits of an int.
unsigned int countSetBits(unsigned int n)
{
    unsigned int count = 0;
    while (n) {
        count += n & 1;
        n >>= 1;
    }
    return count;
}