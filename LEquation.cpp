#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include "LEquation.h"
#include "Test.h"
#define	TwoPi  6.28318530717958648
#define SWAP(a,b) { t=b; b=a; a=t; }
const double eps=1e-10;
using namespace std;

        one_equation::one_equation ( double a, double b ) {
            this->a=a;
            this->b=b;
            }
        one_equation:: SolveP1 (double *x, double a, double b) { // solve equation a*x + b = 0

            if (a==0&&b==0 ) return 1;
            if (a==0&&b!=0) return -1;
            if (a!=0&&b==0) {x[0] = 0; return 0;}
            if (a!=0&&b!=0) {x[0] = -b/a;    return 0;}
        }
        void one_equation:: print_P1 (double *x, double a, double b){
                cout<<"Linear equation (a*x+b=0): ";
                if (one_equation::SolveP1(x,a,b)==0){
                cout<<"x_1 = "<<x[0]<<endl;
                }
                if (one_equation::SolveP1(x,a,b)==1){
                cout<<"x_1 = inf"<<endl;
                }
                if (one_equation::SolveP1(x,a,b)==-1){
                cout<<"x_1 = no root"<<endl;
                }
                }
        void one_equation:: print_file_P1 (double *x, double a, double b){
                ofstream out("Output_LEquation.txt");
                out<<"Linear equation (a*x+b=0): \n";
                if (one_equation::SolveP1(x,a,b)==0){
                out<<"x_1 = "<<x[0];
                }
                if (one_equation::SolveP1(x,a,b)==1){
                out<<"x_1 = inf";
                }
                if (one_equation::SolveP1(x,a,b)==-1){
                out<<"x_1 = no root";
                }
                out.close();
                }
        one_equation::~one_equation(){}

        two_equation::two_equation ( double a, double b, double c ) {
            this->a=a;
            this->b=b;
            this->c=c;
            }
        two_equation:: SolveP2(double *x, double a, double b ) { // solve equation x^2 + a*x + b = 0
            double D;
            D = 0.25*a*a - b;
            if (D >= 0) {
                D = sqrt(D);
                x[0] = -0.5*a + D;
                x[1] = -0.5*a - D;
                return 2;
            }
            x[0] = -0.5*a;
            x[1] = sqrt(-D);
            return 0;
        }

        void two_equation ::print_P2 (double *x, double a, double b, double c){
            cout<<"Quadratic equation (a*x^2+b*x+c=0): ";
            if (!a==0) {
        if( two_equation ::SolveP2(x,b/a,c/a)==2 )
        cout<<"x_1 = "<<x[0]<<"; x_2 = "<<x[1]<<endl;
        if( two_equation ::SolveP2(x,b/a,c/a)==0 ){
        cout<<"x_1 = "<<x[0]<<" + i*"<<x[1]<<"; ";
        cout<<"x_2 = "<<x[0]<<" - i*"<<x[1]<<endl;
        }
        }
        else cout<<"It's not Quadratic equation\n";
        }

        void two_equation ::print_file_P2 (double *x, double a, double b, double c){

        ofstream out("Output_LEquation.txt", ios::app);
        out<<"\nQuadratic equation (a*x^2+b*x+c=0): ";
        if (!a==0) {

        if( two_equation ::SolveP2(x,b/a,c/a)==2 )
        out<<"x_1 = "<<x[0]<<"; x_2 = "<<x[1];
        if( two_equation ::SolveP2(x,b/a,c/a)==0 ){
        out<<"x_1 = "<<x[0]<<" + i*"<<x[1]<<"; ";
        out<<"x_2 = "<<x[0]<<" - i*"<<x[1];
        }
        }
        else out<<"\nIt's not Quadratic equation\n";
        out.close();
        }

       two_equation :: ~two_equation(){}


    three_equation::three_equation ( double a, double b, double c, double d ) {
            this->a=a;
            this->b=b;
            this->c=c;
            this->d=d;
            }

           double three_equation:: _root3 ( double x )
            {
                double s = 1.;
                while ( x < 1. )
                {
                    x *= 8.;
                    s *= 0.5;
                }
                while ( x > 8. )
                {
                    x *= 0.125;
                    s *= 2.;
                }
                double r = 1.5;
                r -= 1./3. * ( r - x / ( r * r ) );
                r -= 1./3. * ( r - x / ( r * r ) );
                r -= 1./3. * ( r - x / ( r * r ) );
                r -= 1./3. * ( r - x / ( r * r ) );
                r -= 1./3. * ( r - x / ( r * r ) );
                r -= 1./3. * ( r - x / ( r * r ) );
                return r * s;
            }

           double three_equation:: root3 ( double x )
            {
                if ( x > 0 ) return three_equation::_root3 ( x ); else
                if ( x < 0 ) return -three_equation::_root3 (-x ); else
                return 0.;
            }

                three_equation::SolveP3(double *x,double a,double b,double c) { // solve equation x^3 + a*x^2 + b*x + c = 0
                double a2 = a*a;
                double q  = (a2 - 3*b)/9;
                double r  = (a*(2*a2-9*b) + 27*c)/54;
                double r2 = r*r;
                double q3 = q*q*q;
                double A,B;
                if (r2 <= (q3 + eps)) {
                    double t=r/sqrt(q3);
                    if( t<-1) t=-1;
                    if( t> 1) t= 1;
                    t=acos(t);
                    a/=3; q=-2*sqrt(q);
                    x[0]=q*cos(t/3)-a;
                    x[1]=q*cos((t+TwoPi)/3)-a;
                    x[2]=q*cos((t-TwoPi)/3)-a;
                    return(3);
                } else {
                    A =-three_equation::root3(fabs(r)+sqrt(r2-q3));
                    if( r<0 ) A=-A;
                    B = (A==0? 0 : B=q/A);

                    a/=3;
                    x[0] =(A+B)-a;
                    x[1] =-0.5*(A+B)-a;
                    x[2] = 0.5*sqrt(3.)*(A-B);
                    if(fabs(x[2])<eps) { x[2]=x[1]; return(2); }
                    return(1);
                }
            }

        void three_equation:: print_P3 (double *x, double a, double b, double c, double d){
            cout<<"Cubic equation (a*x^3+b*x^2+c*x+d=0): ";
            if (!a==0) {


        if( three_equation::SolveP3(x,b/a,c/a,d/a)==3 )
        cout<<"x_1 = "<<x[0]<<"; x_2 = "<<x[1]<<"; x_3 = "<<x[2]<<endl;
        if( three_equation::SolveP3(x,b/a,c/a,d/a)==2 )
        cout<<"x_1 = "<<x[0]<<"; x_2 = "<<x[1]<<endl;
        if( three_equation::SolveP3(x,b/a,c/a,d/a)==1 ){
        cout<<"x_1 = "<<x[0]<<"; ";
        cout<<"x_2 = "<<x[1]<<" + i*"<<x[2]<<"; ";
        cout<<"x_3 = "<<x[1]<<" - i*"<<x[2]<<endl;
        }
        }
        else cout<<"It's not Cubic equation\n";
        }
        void three_equation:: print_file_P3 (double *x, double a, double b, double c, double d){
            ofstream out("Output_LEquation.txt", ios::app);
            out<<"Cubic equation (a*x^3+b*x^2+c*x+d=0): \n";
            if (!a==0) {
        if( three_equation::SolveP3(x,b/a,c/a,d/a)==3 )
        out<<"x_1 = "<<x[0]<<"; x_2 = "<<x[1]<<"; x_3 = "<<x[2];
        if( three_equation::SolveP3(x,b/a,c/a,d/a)==2 )
        out<<"x_1 = "<<x[0]<<"; x_2 = "<<x[1];
        if( three_equation::SolveP3(x,b/a,c/a,d/a)==1 ){
        out<<"x_1 = "<<x[0]<<"; ";
        out<<"x_2 = "<<x[1]<<" + i*"<<x[2]<<"; ";
        out<<"x_3 = "<<x[1]<<" - i*"<<x[2];

        }
        }
         else out<<"It's not Cubic equation";
         out.close();
        }

       three_equation:: ~three_equation(){}


    four_equation::four_equation ( double a, double b, double c, double d, double e ) {
            this->a=a;
            this->b=b;
            this->c=c;
            this->d=d;
            this->e=e;
            }

            void four_equation:: CSqrt( double x, double y, double &a, double &b) // returns:  a+i*s = sqrt(x+i*y)
            {
                double r  = sqrt(x*x+y*y);
                if( y==0 ) {
                    r = sqrt(r);
                    if(x>=0) { a=r; b=0; } else { a=0; b=r; }
                } else {
                    a = sqrt(0.5*(x+r));
                    b = 0.5*y/a;
                }
            }

            four_equation :: SolveP4Bi(double *x, double b, double d) // solve equation x^4 + b*x^2 + d = 0
            {
                double D = b*b-4*d;
                if( D>=0 )
                {
                    double sD = sqrt(D);
                    double x1 = (-b+sD)/2;
                    double x2 = (-b-sD)/2;
                    if( x2>=0 )                         // 0 <= x2 <= x1, 4 real roots
                    {
                        double sx1 = sqrt(x1);
                        double sx2 = sqrt(x2);
                        x[0] = -sx1;
                        x[1] =  sx1;
                        x[2] = -sx2;
                        x[3] =  sx2;
                        return 4;
                    }
                    if( x1 < 0 )                    // x2 <= x1 < 0, two pair of imaginary roots
                    {
                        double sx1 = sqrt(-x1);
                        double sx2 = sqrt(-x2);
                        x[0] =    0;
                        x[1] =  sx1;
                        x[2] =    0;
                        x[3] =  sx2;
                        return 0;
                    }
                                                    // now x2 < 0 <= x1 , two real roots and one pair of imginary root
                        double sx1 = sqrt( x1);
                        double sx2 = sqrt(-x2);
                        x[0] = -sx1;
                        x[1] =  sx1;
                        x[2] =  0;
                        x[3] =  sx2;
                        return 2;
                } else {                            // if( D < 0 ), two pair of compex roots
                    double sD2 = 0.5*sqrt(-D);
                    four_equation::CSqrt(-0.5*b, sD2, x[0],x[1]);
                    four_equation::CSqrt(-0.5*b,-sD2, x[2],x[3]);
                    return 0;
                }
            }

             void four_equation:: dblSort3( double &a, double &b, double &c)        // make: a <= b <= c
            {
                double t;
                if( a>b ) SWAP(a,b);
                if( c<b ) {
                    SWAP(b,c);
                    if( a>b ) SWAP(a,b);
                }
            }

            four_equation:: SolveP4De(double *x, double b, double c, double d)      // solve equation x^4 + b*x^2 + c*x + d
            {
               if( fabs(c)<1e-14*(fabs(b)+fabs(d)) ) return four_equation::SolveP4Bi(x,b,d);

                three_equation Object(*x, a, b, c);

                int res3 = Object.SolveP3(x, 2*b, b*b-4*d, -c*c); // by Viet theorem:  x1*x2*x3=-c*c not equals to 0, so x1!=0, x2!=0, x3!=0
                if( res3>1 )                                        // 3 real roots,
                {
                    four_equation::dblSort3(x[0], x[1], x[2]);
                    if( x[0] > 0)                   // all roots are positive
                    {
                        double sz1 = sqrt(x[0]);
                        double sz2 = sqrt(x[1]);
                        double sz3 = sqrt(x[2]);
                        if( c>0 )
                        {
                            x[0] = (-sz1 -sz2 -sz3)/2;
                            x[1] = (-sz1 +sz2 +sz3)/2;
                            x[2] = (+sz1 -sz2 +sz3)/2;
                            x[3] = (+sz1 +sz2 -sz3)/2;
                            return 4;
                        }

                        x[0] = (-sz1 -sz2 +sz3)/2;
                        x[1] = (-sz1 +sz2 -sz3)/2;
                        x[2] = (+sz1 -sz2 -sz3)/2;
                        x[3] = (+sz1 +sz2 +sz3)/2;
                        return 4;
                    }
                    double sz1 = sqrt(-x[0]);
                    double sz2 = sqrt(-x[1]);
                    double sz3 = sqrt( x[2]);

                    if( c>0 )
                    {
                        x[0] = -sz3/2;
                        x[1] = ( sz1 -sz2)/2;
                        x[2] =  sz3/2;
                        x[3] = (-sz1 -sz2)/2;
                        return 0;
                    }
                    x[0] =   sz3/2;
                    x[1] = (-sz1 +sz2)/2;
                    x[2] =  -sz3/2;
                    x[3] = ( sz1 +sz2)/2;
                    return 0;
                }
                if (x[0] < 0) x[0] = 0;
                double sz1 = sqrt(x[0]);
                double szr, szi;
               four_equation:: CSqrt(x[1], x[2], szr, szi);
                if( c>0 )
                {
                    x[0] = -sz1/2-szr;
                    x[1] = -sz1/2+szr;
                    x[2] = sz1/2;
                    x[3] = szi;
                    return 2;
                }
                x[0] = sz1/2-szr;
                x[1] = sz1/2+szr;
                x[2] = -sz1/2;
                x[3] = szi;
                return 2;
            }

            double four_equation::N4Step(double x, double a,double b,double c,double d) // for x^4 + a*x^3 + b*x^2 + c*x + d
            {
                double fxs= ((4*x+3*a)*x+2*b)*x+c;	// f'(x)
                if (fxs == 0) return x;
                double fx = (((x+a)*x+b)*x+c)*x+d;	// f(x)
                return x - fx/fxs;
            }

            four_equation::  SolveP4(double *x,double a,double b,double c,double d) { // solve equation x^4 + a*x^3 + b*x^2 + c*x + d by Dekart-Euler method
                double d1 = d + 0.25*a*( 0.25*b*a - 3./64*a*a*a - c);
                double c1 = c + 0.5*a*(0.25*a*a - b);
                double b1 = b - 0.375*a*a;
                int res = four_equation::SolveP4De( x, b1, c1, d1);
                if( res==4) { x[0]-= a/4; x[1]-= a/4; x[2]-= a/4; x[3]-= a/4; }
                else if (res==2) { x[0]-= a/4; x[1]-= a/4; x[2]-= a/4; }
                else             { x[0]-= a/4; x[2]-= a/4; }

                if( res>0 )
                {
                    x[0] = N4Step(x[0], a,b,c,d);
                    x[1] = N4Step(x[1], a,b,c,d);
                }
                if( res>2 )
                {
                    x[2] = N4Step(x[2], a,b,c,d);
                    x[3] = N4Step(x[3], a,b,c,d);
                }
                return res;
            }

        void four_equation:: print_P4 (double *x, double a, double b, double c, double d, double e){
            cout<<"Quatre equation (a*x^4+b*x^3+c*x^2+d*x+e=0): ";
            if (!a==0) {

        if( SolveP4(x,b/a,c/a,d/a,e/a)==4 )
        cout<<"x_1 = "<<x[0]<<"; x_2 = "<<x[1]<<"; x_3 = "<<x[2]<<"; x_4 = "<<x[3]<<endl;

    if( SolveP4(x,b/a,c/a,d/a,e/a)==2 )
        cout<<"x_1 = "<<x[0]<<"; x_2 = "<<x[1]<<"; x_3 = "<<x[2]<<" + i*"<<x[3]<<"; x_4 = "<<x[2]<<" - i*"<<x[3]<<endl;
    if( SolveP4(x,b/a,c/a,d/a,e/a)==0 )
        cout<<"x_1 = "<<x[0]<<" + i*"<<x[1]<<"; x_2 = "<<x[0]<<" - i*"<<x[1]<<"; x_3 = "<<x[2]<<" + i*"<<x[3]<<"; x_4 = "<<x[2]<<" - i*"<<x[3]<<endl;

        }
        else cout<<"It's not Quatre equation\n";
        }

        void four_equation:: print__file_P4 (double *x, double a, double b, double c, double d, double e){
            ofstream out("Output_LEquation.txt", ios::app);
            out<<"\nQuatre equation (a*x^4+b*x^3+c*x^2+d*x+e=0): \n";
            if (!a==0) {
        if( SolveP4(x,b/a,c/a,d/a,e/a)==4 )
        out<<"x_1 = "<<x[0]<<"; x_2 = "<<x[1]<<"; x_3 = "<<x[2]<<"; x_4 = "<<x[3];
    if( SolveP4(x,b/a,c/a,d/a,e/a)==2 )
        out<<"x_1 = "<<x[0]<<"; x_2 = "<<x[1]<<"; x_3 = "<<x[2]<<" + i*"<<x[3]<<"; x_4 = "<<x[2]<<" - i*"<<x[3];
    if( SolveP4(x,b/a,c/a,d/a,e/a)==0 )
        out<<"x_1 = "<<x[0]<<" + i*"<<x[1]<<"; x_2 = "<<x[0]<<" - i*"<<x[1]<<"; x_3 = "<<x[2]<<" + i*"<<x[3]<<"; x_4 = "<<x[2]<<" - i*"<<x[3];
        out.close();
        }
        else out<<"It's not Quatre equation\n";
        }
        four_equation::~four_equation(){}
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
        if (Four_E.SolveP4(x1,1,1,1,1)==0) cout<<"Tested_4\n\n";

}

int main()
{
    int a,b,c,d,e;
    double x[4];
    char R;
    test();
    cout<<"Enter from file/keyboard press f/k ";
    cin>>R;
    if (R=='k')
    {
    cout<<endl<<"Enter koeficients: a, b, c, d, e \n";
    cin>>a>>b>>c>>d>>e;

        one_equation One_E(a,b);
        One_E.print_P1(x,a,b);

        two_equation Two_E(a,b,c);
        Two_E.print_P2(x,a,b,c);

        three_equation Three_E (a,b,c,d);
        Three_E.print_P3(x,a,b,c,d);

        four_equation Four_E (a,b,c,d,e);
        Four_E.print_P4(x,a,b,c,d,e);
    }

    if (R=='f') {
    cout<<"Koeficients from Input_LEquation\n";
    ifstream in("Input_LEquation.txt");
    in>>a;
    in>>b;
    in>>c;
    in>>d;
    in>>e;
    {
    one_equation One_E(a,b);
    One_E.print_file_P1(x,a,b);

    two_equation Two_E(a,b,c);
    Two_E.print_file_P2(x,a,b,c);

    three_equation Three_E (a,b,c,d);
    Three_E.print_file_P3(x,a,b,c,d);

    four_equation Four_E (a,b,c,d,e);
    Four_E.print__file_P4(x,a,b,c,d,e);
    }

    cout<<"Writed to Output_LEquation";
    in.close();
    }
    return 0;
}
