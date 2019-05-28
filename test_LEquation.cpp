#include <iostream>
#include "LEquation.h"
#include "Test.h"

void test ()
{
        double x1[4];
        one_equation One_E(1,2);
        if (One_E.SolveP1(x1,1,2)==0) cout<<"Tested_1\n";
        two_equation Two_E(2,9,9);
        if( Two_E.SolveP2(x1,9/2,9/2)==2 ) cout<<"Tested_2\n";
        three_equation Three_E (5,15,-4,-12);
        if (Three_E.SolveP3(x1,15/5,-4/5,-12/5)==3) cout<<"Tested_3\n";
        four_equation Four_E (1,1,1,1,1);
        if (Four_E.SolveP4(x1,1,1,1,1)==0) cout<<"Tested_4\n";

}
