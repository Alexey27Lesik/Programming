class one_equation
{
    protected:
        double a;
        double b;
    public:

        one_equation ( double a, double b );
        int SolveP1 (double *x, double a, double b);
        void print_P1 (double *x, double a, double b);
        ~one_equation();
};
class two_equation
{
    protected:
        double a;
        double b;
        double c;

    public:
        two_equation ( double a, double b, double c );
        int   SolveP2(double *x, double a, double b );
        void print_P2 (double *x, double a, double b, double c);
        ~two_equation();
};

class three_equation
{
    protected:
        double a;
        double b;
        double c;
        double d;

    public:
        three_equation ( double a, double b, double c, double d );
        static double _root3 ( double x );
        double root3 ( double x );
         int SolveP3(double *x,double a,double b,double c);
        void print_P3 (double *x, double a, double b, double c, double d);
        ~three_equation();
};

class four_equation
{

    protected:
        double a;
        double b;
        double c;
        double d;
        double e;

    public:
        four_equation ( double a, double b, double c, double d, double e );
        void  CSqrt( double x, double y, double &a, double &b);
        int   SolveP4Bi(double *x, double b, double d);
        static void  dblSort3( double &a, double &b, double &c);
        int   SolveP4De(double *x, double b, double c, double d);
        double N4Step(double x, double a,double b,double c,double d);
        int   SolveP4(double *x,double a,double b,double c,double d);
        void  print_P4 (double *x, double a, double b, double c, double d, double e);
        ~four_equation();
};
