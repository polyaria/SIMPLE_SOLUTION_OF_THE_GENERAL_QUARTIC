#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

using namespace std;

template<typename fp_t>
fp_t discriminant(fp_t a, fp_t b, fp_t c)
{
    fp_t D = fma(b,b, -4*a*c);
    return D;
}
template <typename fp_t>
int sgn(fp_t val) {
    return (fp_t(0) < val) - (val < fp_t(0));
}

//РЕШЕНИЕ КУБИЧЕСКОГО УРАВНЕНИЯ
template<typename fp_t>
fp_t solve_cubic(std::vector<fp_t> coefficients, std::vector<fp_t> &roots) {
    fp_t a, b, c, d, t, A, B, Q, R, D;

    //Coefficients - делим на коэффициент, стоящий перед z^3
    a = coefficients[1]/coefficients[0];
    b = coefficients[2]/coefficients[0];
    c = coefficients[3]/coefficients[0];

    cout << "\t " << "z^3 + " << a << "*z^2 + " << b <<"*z + " << c << " = 0 "<< endl;
    Q = ((long double)pow(a,2) - 3*b)/9;
    R = ((long double)2*pow(a,3) - 9*a*b + 27*c)/54;

    if (pow(R,2) < pow(Q,3)) //то уравнение имеет три действительных корня
    {
        t = acos(R/(sqrt((long double)pow(Q,3))))/3;
        roots[0] = -2*sqrt(Q)*cos(t) - a/3;
        roots[1] = -2*sqrt(Q)*cos(t+(2*(numbers::pi/3))) - a/3;
        roots[2] = -2*sqrt(Q)*cos(t - (2*(numbers::pi/3))) - a/3;
        cout<<"\tThree real roots" << endl;
        cout << "\t\tRoots: " << endl;
        cout << "\t z1 = " << roots[0] << "; z2 = " << roots[1] << "; z3 = " << roots[2] << endl;

    }
    if (pow(R,2) >= pow(Q,3)) //то 1 действительный корень(общий случай) или 2 (вырожденные случаи)
    {
        A = -sgn(R)*pow((long double)fabs(R)+sqrt(pow(R,2)-pow(Q,3)),(long double)1/3);
        B = sgn(A) == 0 ? 0 : Q/A;
        roots[0] = (A+B) - a/3; //действительный корень
        cout<< "\t\t1 real root"<<endl;
        cout << "\t\tz1 = " << roots[0]<<endl;
        //Комплексно-сопряженные корни
        if (A!=B) {
            std::complex<fp_t> root1(-(A + B) / 2 - a / 3, sqrt(3) * (A - B) / 2);
            std::complex<fp_t> root2(-(A + B) / 2 - a / 3, -sqrt(3) * (A - B) / 2);
            cout << "\t\tComplex roots: " << endl;
            cout << "\tz2 = " << root1 << "; z3 = " << root2 << endl;
        }
        else { //комплекнсые корни вырождаются в действительный
            roots[1] = -A - a/3;
            cout<< "\t\t2 real root"<<endl;
            cout << "\t\tz2 = " << roots[1] << endl;
        }
     }
    return roots[0]; //возвращаем один действительный корень - он нужен для дальнейшего решения

}

//РЕШЕНИЕ УРАВНЕНИЯ 4ОЙ СТЕПЕНИ
template<typename fp_t>
void solution(std::vector<fp_t> coefficients,std::vector<fp_t> coefficients1, std::vector<fp_t> &roots)
{
    fp_t a, b, c, d, e, a1, b1, c1, d1, P, Q, R,K, L, KL;
    std::vector<fp_t> roots_y(4);

    //Coefficients
    a = coefficients[0];b = coefficients[1];c = coefficients[2];d = coefficients[3];e = coefficients[4];
    // a*x^4 + b*x^3 + c*x^2 + d*x + c = 0
    cout << endl << "\t " << a << "*x^4 + " << b << "*x^3 + " << c << "*x^2 + " << d <<"*x + " << e << " = 0 "<< endl;

    // y^4 + P*y^2 + Q*y +R = 0
    //Замена
    //y = x + b/(4*a);
    P = (-3*pow(b,2) + 8*a*c)/(8*pow(a,2));
    Q = (pow(b,3) - 4*a*b*c +8*pow(a,2)*d)/((8*pow(a,3)));
    R = (-3*pow(b,4) + 16*a*pow(b,2)*c - 64*pow(a,2)*b*d + 256*pow(a,3)*e)/(256*pow(a,4));
    cout << "\t " << "y = x + " << b/(4*a) << endl;
    cout << "\t "  << "y^4 + " << P << "*y^2 + " << Q <<"*y + " << R << " = 0 "<< endl;
    cout << "\t" << "(y^2 + z)^2 = y^4 + 2y^2*z + z^2" <<endl;
    cout << "\t "<< "(y^2 + z)^2 = " << "(2*z" << - P<< ") *y^2 + "<<- Q <<"*y +(z^2 + " << - R << ")" << endl;
    cout << "\t "<< "Discriminant = 0" << endl;
    cout << "\t "<< "8*z^3 - 4*"<<P<<"*z^2 - 8*" << R<<"*z +4*"<< (P*R - pow(Q,2)) <<" = 0 " << endl;

    a1 = coefficients1[0];b1 = coefficients1[1];c1 = coefficients1[2];d1 = coefficients1[3];

    a1 = 8; b1 = -4*P; c1 = -8*R; d1 = 4*P*R - pow(Q,2);
    vector<fp_t> koef = {a1, b1, c1, d1};
    vector<fp_t> r(3);
    //Решаем кубическое уравнение
    // 8*z^3 - 4*P*z^2 - 8*R*z +4*P*R - Q^2 = 0
    fp_t z1 =solve_cubic(koef, r); // возвращаем один действительный корень

    //(y^2 + z1)^2 = K^2*y^2-2*K*L*y +L^2, где
    //K^2 =2*z1 - P,  L^2 = z1^2 - R, 2*K*L = Q
    //+-(y^2 + z1) = K*y - L
    //Раскрывается на 2 уравнения:
    // * y^2 - K*y +z1 +L = 0
    // * y^2 + K*y +z1 -L = 0
    // Корни которых:
    // * y = (K+- sqrt(K^2 - 4*(z1+L))/2
    // * y = (-K+- sqrt(K^2 - 4*(z1-L))/2

    //КАК ВЫБИРАТЬ ЗНАК?
    //тут получаются случаи, когда K, L  - КОМПЛЕКСНЫЕ - как быть?

    if (2*z1 - P >=0 ) K = sqrt(2*z1 -P);
    else cout << "\t "<< "K - complex" << endl ;
    if (pow(z1,2) -R >= 0)  L = -sqrt(pow(z1,2) - R);
    else cout << "\t "<< "L - complex" << endl ;
    KL = Q/2;

    cout << "\t "<< "K = "<< K<<endl;
    cout<< "\t "<< "L = "<< L << endl;
    cout << "\t "<< "KL = "<< KL<<endl;
    if (pow(K,2) - 4*(z1 + L ) >= 0) { // вещественные корни

        roots_y[0] = (K + sqrt(pow(K, 2) - 4 * (z1 + L))) / 2;
        roots_y[1] = (K - sqrt(pow(K, 2) - 4 * (z1 + L))) / 2;

        //cout << "\t y0 = " <<   roots_y[0]<<endl;
       // cout << "\t y1 = " <<   roots_y[1]<<endl;
        //обратно возвращаемся к x
        cout << "\t" << "x = y - "<< + b/(4*a) << endl;
        roots[0] = roots_y[0] - b/(4*a);
        roots[1] = roots_y[1] - b/(4*a);
        cout << "\t"<< "Real roots: "<< endl;
        cout << "\t x0 = " <<  roots[0]<<endl;
        cout << "\t x1 = " <<  roots[1]<<endl;
    }
    else { //комплексные корни
        std::complex<fp_t> sqrt_y = std::sqrt(std::complex<fp_t>(pow(K, 2) - 4 * (z1 + L)));
        std::complex<fp_t> root_y0(K/2, sqrt_y.imag()/2);
        std::complex<fp_t> root_y1(K/2, -sqrt_y.imag()/2);

       // cout << "\t y0 =  "<< root_y0<<endl;
      //  cout << "\t y1 =  "<< root_y1<<endl;
        std::complex<fp_t> root_x0(K/2 - b/(4*a), sqrt_y.imag()/2);
        std::complex<fp_t> root_x1(K/2- b/(4*a), -sqrt_y.imag()/2);

        cout << "\t" << "x = y - "<< + b/(4*a) << endl;
        cout << "\t "<< "Complex roots: "<< endl;
        cout << "\t x0 =  "<< root_x0<<endl;
        cout << "\t x1 =  "<< root_x1<<endl;
    }
    if (pow(K,2) - 4*(z1 - L ) >= 0)
    {
        roots_y[2] = (-K+sqrt(pow(K,2) - 4*(z1 - L )))/2;
        roots_y[3] = (-K-sqrt(pow(K,2) - 4*(z1 - L )))/2;

       // cout<< "\t y2 =  "<< roots_y[2]<<endl;
       // cout<< "\t y3 =  "<< roots_y[3]<<endl;
        //обратно возвращаемся к x
        cout << "\t" << "x = y - "<< + b/(4*a) << endl;
        roots[2] = roots_y[2] - b/(4*a);
        roots[3] = roots_y[3] - b/(4*a);
        cout << "\t "<< "Real roots: "<< endl;
        cout<< "\t x2 =  "<< roots[2]<<endl;
        cout<< "\t x3 =  "<< roots[3]<<endl;
    }
    else {
        std::complex<fp_t> sqrt_y1 = std::sqrt(std::complex<fp_t>(pow(K, 2) - 4 * (z1 - L)));
        std::complex<fp_t> root_y3(-K/2, sqrt_y1.imag()/2 );
        std::complex<fp_t> root_y4(-K/2, -sqrt_y1.imag()/2 );

       // cout << "\t y2 = "<< root_y3<<endl;
       // cout << "\t y3 = "<< root_y4<<endl;
        std::complex<fp_t> root_x3(-K/2 - b/(4*a), sqrt_y1.imag()/2);
        std::complex<fp_t> root_x4(-K/2 - b/(4*a), -sqrt_y1.imag()/2);
        cout << "\t" << "x = y - "<< + b/(4*a) << endl;
        cout << "\t "<< "Complex roots: "<< endl;
        cout << "\t x2 = "<< root_x3<<endl;
        cout << "\t x3 = "<< root_x4<<endl;
    }

}

int main() {

    cout << "\t\tEXAMPLE 1";
    vector<float> koef1 = {1, -4, 7, -10, 3};
    vector<float> koef2 = {0, 0, 0, 0};
    vector<float> r1(4);
    solution(koef1,koef2, r1);


    cout << "\n\t\tEXAMPLE 2";
    vector<float> koef3 = {2,0,5,0,-3};
    vector<float> koef4 = {0, 0, 0, 0};
    vector<float> r2(4);
    solution(koef3,koef4, r2);

    cout << "\n\t\tEXAMPLE 3";
    vector<float> koef5 = {16,0,145,0,9};
    vector<float> koef6 = {0, 0, 0, 0};
    vector<float> r3(4);
    solution(koef5,koef6, r3);

    cout << "\n\t\tEXAMPLE 4";
    vector<float> koef7 = {1,3,3,-1,-6};
    vector<float> koef8 = {0, 0, 0, 0};
    vector<float> r4(4);
    solution(koef7,koef8, r4);

    cout << "\n\t\tEXAMPLE 5";
    vector<float> koef9 = {1,4,-4,-20,-5};
    vector<float> koef10 = {0, 0, 0, 0};
    vector<float> r5(4);
    solution(koef9,koef10, r5);


    return 0;
}
