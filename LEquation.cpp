#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include "LEquation.h"
#define	TwoPi  6.28318530717958648
#define SWAP(a,b) { t=b; b=a; a=t; }
const double eps=1e-10;
using namespace std;

        one_equation::one_equation ( double a, double b ) {
            this->a=a;
            this->b=b;
            }
        one_equation:: SolveP1 (double *x, double a, double b) {
            if ( a==0&&b!=0) return -1;
            if ( a==0&&b==0) return 255;
            if ( b==0 )  x[0] = 0;
            else x[0] = -b/a;
            return 0;
        }
        void one_equation:: print_P1 (double *x, double a, double b){
                cout<<"Linear equation: ";
                if (one_equation::SolveP1(x,a,b)==0){
                cout<<"x_1 = "<<x[0]<<endl;
                }
                if (one_equation::SolveP1(x,a,b)==-1){
                cout<<"No root"<<endl;
                }
                if (one_equation::SolveP1(x,a,b)==255){
                cout<<"Inf root"<<endl;
                }
                }
        one_equation::~one_equation(){}

        two_equation::two_equation ( double a, double b, double c ) {
            this->a=a;
            this->b=b;
            this->c=c;
            }
        two_equation:: SolveP2(double *x, double a, double b ) {
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
            cout<<"Quadratic equation: ";
            if (!a==0) {
        if( two_equation ::SolveP2(x,b/a,c/a)==2 )
        cout<<"x_1 = "<<x[0]<<"; x_2 = "<<x[1]<<endl;
        if( two_equation ::SolveP2(x,b/a,c/a)==0 ){
        cout<<"x_1 = "<<x[0]<<" + i*"<<x[1]<<"; ";
        cout<<"x_2 = "<<x[0]<<" - i*"<<x[1]<<endl;
        }
        }
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

                three_equation::SolveP3(double *x,double a,double b,double c) {
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
            cout<<"Cubic equation: ";
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
        }
       three_equation:: ~three_equation(){}


    four_equation::four_equation ( double a, double b, double c, double d, double e ) {
            this->a=a;
            this->b=b;
            this->c=c;
            this->d=d;
            this->e=e;
            }

            void four_equation:: CSqrt( double x, double y, double &a, double &b)
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

            four_equation::   SolveP4Bi(double *x, double b, double d)
            {
                double D = b*b-4*d;
                if( D>=0 )
                {
                    double sD = sqrt(D);
                    double x1 = (-b+sD)/2;
                    double x2 = (-b-sD)/2;
                    if( x2>=0 )
                    {
                        double sx1 = sqrt(x1);
                        double sx2 = sqrt(x2);
                        x[0] = -sx1;
                        x[1] =  sx1;
                        x[2] = -sx2;
                        x[3] =  sx2;
                        return 4;
                    }
                    if( x1 < 0 )
                    {
                        double sx1 = sqrt(-x1);
                        double sx2 = sqrt(-x2);
                        x[0] =    0;
                        x[1] =  sx1;
                        x[2] =    0;
                        x[3] =  sx2;
                        return 0;
                    }

                        double sx1 = sqrt( x1);
                        double sx2 = sqrt(-x2);
                        x[0] = -sx1;
                        x[1] =  sx1;
                        x[2] =  0;
                        x[3] =  sx2;
                        return 2;
                } else {
                    double sD2 = 0.5*sqrt(-D);
                    four_equation::CSqrt(-0.5*b, sD2, x[0],x[1]);
                    four_equation::CSqrt(-0.5*b,-sD2, x[2],x[3]);
                    return 0;
                }
            }

             void four_equation:: dblSort3( double &a, double &b, double &c)
            {
                double t;
                if( a>b ) SWAP(a,b);
                if( c<b ) {
                    SWAP(b,c);
                    if( a>b ) SWAP(a,b);
                }
            }
           double four_equation:: _root3 ( double x )
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

           double four_equation:: root3 ( double x )
            {
                if ( x > 0 ) return four_equation::_root3 ( x ); else
                if ( x < 0 ) return -four_equation::_root3 (-x ); else
                return 0.;
            }

                four_equation::SolveP3(double *x,double a,double b,double c) {
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

                    A =-four_equation::root3(fabs(r)+sqrt(r2-q3));
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

            four_equation:: SolveP4De(double *x, double b, double c, double d)
            {
               if( fabs(c)<1e-14*(fabs(b)+fabs(d)) ) return four_equation::SolveP4Bi(x,b,d);

                int res3 = four_equation::SolveP3(x, 2*b, b*b-4*d, -c*c);
                if( res3>1 )
                {
                    four_equation::dblSort3(x[0], x[1], x[2]);
                    if( x[0] > 0)
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

            double four_equation::N4Step(double x, double a,double b,double c,double d)
            {
                double fxs= ((4*x+3*a)*x+2*b)*x+c;	// f'(x)
                if (fxs == 0) return x;
                double fx = (((x+a)*x+b)*x+c)*x+d;	// f(x)
                return x - fx/fxs;
            }

            four_equation::  SolveP4(double *x,double a,double b,double c,double d) {
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
            cout<<"Quatre equation: ";
            if (!a==0) {
        if( SolveP4(x,b/a,c/a,d/a,e/a)==4 )
        cout<<"x_1 = "<<x[0]<<"; x_2 = "<<x[1]<<"; x_3 = "<<x[2]<<"; x_4 = "<<x[3]<<endl;
    if( SolveP4(x,b/a,c/a,d/a,e/a)==2 )
        cout<<"x_1 = "<<x[0]<<"; x_2 = "<<x[1]<<"; x_3 = "<<x[2]<<" + i*"<<x[3]<<"; x_4 = "<<x[2]<<" - i*"<<x[3]<<endl;
    if( SolveP4(x,b/a,c/a,d/a,e/a)==0 )
        cout<<"x_1 = "<<x[0]<<" + i*"<<x[1]<<"; x_2 = "<<x[0]<<" - i*"<<x[1]<<"; x_3 = "<<x[2]<<" + i*"<<x[3]<<"; x_4 = "<<x[2]<<" - i*"<<x[3]<<endl;
        }
        }

        four_equation::~four_equation(){}

int main()
{
    int a,b,c,d,e;
    double x[4];
    char R;
    cout<<"Enter from file or keyboard press f/k ";
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
    ifstream in("Input_LEquation.txt");
    in>>a;
    in>>b;
    in>>c;
    in>>d;
    in>>e;
    one_equation One_E(a,b);
    One_E.print_P1(x,a,b);

    two_equation Two_E(a,b,c);
    Two_E.print_P2(x,a,b,c);

    three_equation Three_E (a,b,c,d);
    Three_E.print_P3(x,a,b,c,d);

    four_equation Four_E (a,b,c,d,e);
    Four_E.print_P4(x,a,b,c,d,e);
    in.close();
    }
    return 0;
}
