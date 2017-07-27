#include <iostream>
#include <vector>
#include "interpolation.h"

using namespace std;

int main()
{
    SineInterp sine;
    CubicHermiteInterp cubic;
    SplineInterp spline;

    vector<double> y = {23 93 33 46 80 14 96 42 16};

    auto y1 = sine.interp(y, 10);
    auto y2 = cubic.interp(y, 10);
    auto y3 = spline.interp(y, 10);

    auto print = [](const vector<double> &y)
    {
        for(auto d : y)
        {
            cout << d << "\t";
        }
        cout << endl;
    };

    cout << "sine interp data is: " << y1.size() << endl;
    print(y1);

    cout << "cubic interp data is: " << y2.size() << endl;
    print(y2);

    cout << "spline interp data is: " << y3.size() << endl;
    print(y3);

    return 0;
}

