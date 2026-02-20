#include <bits/stdc++.h>
#include <cryptopp/des.h>
#include <cryptopp/osrng.h>
#include <cryptopp/modes.h>
#include <cryptopp/integer.h>
#include <cryptopp/nbtheory.h>
#include <fstream>

// g++ lab01.cpp -o lab01.out -std=c++17 -lcryptopp
using namespace std;
using bigint = CryptoPP::Integer;
using point = tuple<bigint, bigint, bigint>;
using ll = long long;
const point INF = {0, 1, 0};

// Overload operators
string operator+(string a, bigint b) {
    return a + CryptoPP::IntToString(b);
}
ostream& operator<<(ostream& os, const point& a) {
    os << "{" 
       << get<0>(a) << ", " 
       << get<1>(a) << ", " 
       << get<2>(a) 
       << "}";
    return os;
}

// Helper functions
bigint getPrime(int bits) {
    CryptoPP::AutoSeededRandomPool prng;
    return CryptoPP::MihailescuProvablePrime(prng, bits);
}

bool checkDiscriminant(bigint a, bigint b, bigint p) {
    return ((4*a*a*a + 27*b*b) % p) != 0;
}

bigint evaluateEC(bigint a, bigint b, bigint p, bigint x) {
    return (x*x*x + a*x + b) % p;
}

string formatEC(bigint a, bigint b, bigint p) {
    string y2 = "y^2 = x^3";

    if(a == 1) y2 += " + x";
    else if(a != 0) y2 += " + " + a + "x";
    if(b != 0) y2 += " + " + b;
    y2 += " (mod " + p + ")";

    return y2;
}

ll validEllipticCurves(bigint p) {
    string path = "valid_curves.txt";
    ofstream file(path);
    bool f = 0;
    ll cnt = 0;

    if(file.is_open()) f = 1;
    else {
        cout << "Failed to open file. Continue without saving!" << endl;
    }
    for(bigint a=0; a<p; a++) {
        for(bigint b=0; b<p; b++) {
            if(checkDiscriminant(a, b, p)) {
                string curve = formatEC(a, b, p);
                if(f) file << curve << endl;
                cnt++;
            }
        }
    }
    if(f) file.close();
    return cnt;
}

ll validEllipticCurvesRandomPrime(int bits) {
    string path = "valid_curves_random.txt";
    ofstream file(path);
    ll cnt = 0;
    bool f = 0;
    bigint p = getPrime(bits);

    cout << "Using p: " << p << endl;
    if(file.is_open()) f = 1;
    else {
        cout << "Failed to open file. Continue without saving!" << endl;
    }
    for(bigint a=0; a<p; a++) {
        for(bigint b=0; b<p; b++) {
            if(checkDiscriminant(a, b, p)) {
                string curve = formatEC(a, b, p);
                if(f) file << curve << endl;
                cnt++;
            }
        }
    }
    if(f) file.close();
    return cnt;
}

array<bigint, 2> hassesTheorem(bigint p) {
    bigint sr = (4*p).SquareRoot();   // |_2*√p_| => |√4*p_|
    cout << sr << endl;
    return {p + 1 - sr, p + 1 + sr};
}

ll rationalPoints(bigint a, bigint b, bigint p) {
    // First compute Quadratic Residues (QR)
    // and map valid solutions
    map<bigint, bigint> sr;
    for(bigint i=0; i<=(p-1)/2; i++) {
        bigint r = i*i % p;
        if(!sr.count(r))
            sr[r] = i;
    }

    // For every value of x, evaluate the EC
    vector<array<bigint, 3>> points;
    points.push_back({0, 1, 0});
    for(bigint x=0; x<p; x++) {
        bigint r = evaluateEC(a, b, p, x);
        if(sr.count(r)) {
            bigint y1 = sr[r], y2 = (p - y1) % p;
            points.push_back({x, y1, 1});
            if(y1 != y2) points.push_back({x, y2, 1});
        }   
    }
    sort(points.begin(), points.end());
    points.resize(unique(points.begin(), points.end()) - points.begin());
    for(auto &[x, y, z] : points){
        cout << "(" << x << ", " << y <<", " << z << ")" << endl;
    }
    return points.size();
}

bigint add(bigint x, bigint p) {
    return (x % p + p) % p;
}

point addPoints(point p1,  point q1, bigint p) {
    if(p1 == INF) return q1;
    if(q1 == INF) return p1;
    auto &[x1, y1, z1] = p1;
    auto &[x2, y2, z2] = q1;
    
    bigint num = (y2 - y1), den = (x2 - x1);
    if(!den) return INF;

    bigint lambda = num * den.InverseMod(p) % p;
    bigint x3 = add(lambda*lambda - x1 - x2, p);
    bigint y3 = add(lambda * (x1 - x3) - y1, p);
    return {x3, y3, 1};
}

point doublePoint(point p1, bigint a, bigint p) {
    if(p1 == INF) return INF;
    auto &[x, y, z] = p1;
    if(!y) return INF;
    bigint lambda = (3*x*x + a) * (2*y).InverseMod(p) % p;
    bigint x3 = add(lambda*lambda - 2*x, p);
    bigint y3 = add(lambda * (x - x3) - y, p);
    return {x3, y3, 1};
}

void count() {
    bigint p;
    cout << "Enter a prime number: ";
    cin >> p;
    auto [l, r] = hassesTheorem(p);
    cout << "[" << l << ", " << r << "]" << endl;
}

void menu() {
    cout << "* * * Lab01 - Elliptic Curves * * *" << endl;
    cout << "[1] Find valid curves." << endl;
    cout << "[2] Finv valid curves for a random prime." << endl;
    cout << "[3] Find rational points." << endl;
    cout << "[4] Add two points." << endl;
    cout << "[5] Double a point." << endl;
    cout << "[0] Exit." << endl;
    cout << "Select an option: ";
}

int main() {
    ios::sync_with_stdio(0);
    int opc, bits;
    ll res;
    bool f = true;
    bigint a, b, p;
    point p1, q1;
    auto &[x1, y1, z1] = p1;
    auto &[x2, y2, z2] = q1;

    while(f) {
        menu();
        cin >> opc;
        cout << endl;

        switch(opc) {
            case 1:
                cin >> p;
                res = validEllipticCurves(p);
                cout << "\nfound " << res << " non-singular curves!" << endl;
                break;

            case 2:
                cout << "Enter number of bits: ";
                cin >> bits;
                while(bits < 3) {
                    cout << "Please, enter a number greater than 2!: ";
                    cin >> bits;
                }
                res = validEllipticCurvesRandomPrime(bits);
                cout << "\nfound " << res << " non-singular curves!" << endl;
                break;

            case 3:
                cout << "Enter value of a: ";
                cin >> a;
                
                cout << "Enter value of b: ";
                cin >> b;
                cout << "Enter a prime number: ";
                cin >> p;
                cout << rationalPoints(a, b, p) << endl;
                break;

            case 4:
                cout << "Enter Point P1" << endl;
                cout << "Enter value of x1: ";
                cin >> x1;
                cout << "Enter value of y1: ";
                cin >> y1;
                cout << "Enter value of z1: ";
                cin >> z1;

                cout << "Enter Point Q1" << endl;
                cout << "Enter value of x2: ";
                cin >> x2;
                cout << "Enter value of y2: ";
                cin >> y2;
                cout << "Enter value of z2: ";
                cin >> z2;

                cout << "Enter a prime number: ";
                cin >> p;

                cout << addPoints(p1, q1, p);
                break;
            
            case 5:
                cout << "Enter Point P1" << endl;
                cout << "Enter value of x: ";
                cin >> x1;
                cout << "Enter value of y: ";
                cin >> y1;
                cout << "Enter value of a: ";
                cin >> a;
                cout << "Enter a prime number: ";
                cin >> p;
                cout << doublePoint(p1, a, p) << endl;
                break;

            case 0:
                f = false;
                break;
            default:
                cout << "Ingresa una opción válida" << endl;
        }
    }

    return 0;
}

