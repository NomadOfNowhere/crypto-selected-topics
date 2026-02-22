#include "crow_all.h"
#include <cryptopp/des.h>
#include <cryptopp/osrng.h>
#include <cryptopp/modes.h>
#include <cryptopp/integer.h>
#include <cryptopp/nbtheory.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <map>

// g++ lab01.cpp -o lab01 -std=c++17 -pthread -lcryptopp -O3
using namespace std;
using bigint = CryptoPP::Integer;
using point = tuple<bigint, bigint, bigint>;
using ll = long long;
const point INF = {0, 1, 0};
const int LIMIT = 1e6;

// Overload operators
string to_string(const bigint &a) { return CryptoPP::IntToString(a); }
string operator+(const string &a, const bigint &b) { return a + to_string(b); }
string to_string(const point &a) { 
    return "{" + get<0>(a) + ", " + 
                 get<1>(a) + ", " +
                 get<2>(a) + "}";
}
ostream& operator<<(ostream &os, const point &a) {
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

bigint add(bigint x, bigint p) {
    return (x % p + p) % p;
}

vector<string> getValidCurves(bigint p) {
    vector<string> curves;
    for(bigint a=0; a<p; a++) {
        for(bigint b=0; b<p; b++) {
            if(checkDiscriminant(a, b, p)) {
                if(curves.size() > 1e8) break;
                    curves.push_back(formatEC(a, b, p));
            }
        }
    }
    return curves;
}

vector<point> getRationalPoints(bigint a, bigint b, bigint p) {
    // First compute Quadratic Residues (QR)
    // and map valid solutions
    map<bigint, bigint> sr;
    for(bigint i=0; i<=(p-1)/2; i++) {
        bigint r = i*i % p;
        if(!sr.count(r))
            sr[r] = i;
    }

    // For every value of x, evaluate the EC
    vector<point> points;
    points.push_back({0, 1, 0});
    for(bigint x=0; x<p; x++) {
        bigint r = evaluateEC(a, b, p, x);
        if(sr.count(r)) {
            bigint y1 = sr[r], y2 = (p - y1) % p;
            points.push_back({x, y1, 1});
            if(y1 != y2) points.push_back({x, y2, 1});
        }   
    }

    // Sort and get rid of duplicates
    sort(points.begin(), points.end());
    points.resize(unique(points.begin(), points.end()) - points.begin());

    return points;
}

point getAddPoints(point p1,  point q1, bigint p) {
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

point getDoublePoint(point p1, bigint a, bigint p) {
    if(p1 == INF) return INF;
    auto &[x, y, z] = p1;
    if(!y) return INF;
    bigint lambda = (3*x*x + a) * (2*y).InverseMod(p) % p;
    bigint x3 = add(lambda*lambda - 2*x, p);
    bigint y3 = add(lambda * (x - x3) - y, p);
    return {x3, y3, 1};
}

pair<bigint, bigint> getHassesTheorem(bigint p) {
    bigint sr = (4*p).SquareRoot();   // |_2*√p_| => |√4*p_|
    return {p + 1 - sr, p + 1 + sr};
}

// Parse functions
bigint parse_bigint(const crow::json::rvalue &body, const string &key) {
    string value = body[key].s();
    return bigint(value.c_str());
}
point parse_point(const crow::json::rvalue &body, const string &key) {
    crow::json::rvalue arr = body[key];
    // if (!arr || arr.size() != 3)
    //     throw std::runtime_error("Invalid point");

    bigint x = parse_bigint(arr, "x");
    bigint y = parse_bigint(arr, "y");
    bigint z = parse_bigint(arr, "z");

    return {x, y, z};
}

// Send response as Headers CORS
crow::response response(crow::json::wvalue res) {
    auto response = crow::response(200, res);
    response.add_header("Access-Control-Allow-Origin", "*");
    response.add_header("Access-Control-Allow-Headers", "Content-Type");
    return response;
}

// Convert to JSON
crow::json::wvalue serialize(const string &s) { return s; }
crow::json::wvalue serialize(const point &p) {
    // return to_string(p);
    crow::json::wvalue obj;
    obj["x"] = to_string(get<0>(p));
    obj["y"] = to_string(get<1>(p));
    obj["z"] = to_string(get<2>(p));
    return obj;
}
crow::json::wvalue serialize(const pair<bigint, bigint> &p) {
    crow::json::wvalue obj;
    obj["l"] = to_string(p.first);
    obj["r"] = to_string(p.second);
    return obj;
}
// Convert vector to JSON list
template <typename T>
crow::json::wvalue vector_to_json(const vector<T> &values) {
    crow::json::wvalue obj;
    for(int i=0; i<values.size(); i++) {
        obj[i] = serialize(values[i]);
    }
    return obj;
}

int main() {
    crow::SimpleApp app;

    CROW_ROUTE(app, "/")
    ([](){
        return "Server running!";
    });

    CROW_ROUTE(app, "/api/valid_curves").methods(crow::HTTPMethod::POST)
    ([](const crow::request &req) {
        // Parse JSON input
        std::cout << "CROW RECIBIÓ ESTO: [" << req.body << "]" << std::endl;
        auto body = crow::json::load(req.body);
        bool valid = body.has("p");
        if(!body || !valid) return crow::response(400, "Invalid JSON");

        // Extract prime number as string
        bigint p = parse_bigint(body, "p");

        vector<string> curves = getValidCurves(p);

        // Pack response as JSON
        crow::json::wvalue res;
        res["count"] = curves.size();
        res["curves"] = vector_to_json(curves);
        
        return response(res);
    });

    CROW_ROUTE(app, "/api/valid_random").methods(crow::HTTPMethod::POST)
    ([](const crow::request &req) {
        // Parse JSON input
        auto body = crow::json::load(req.body);
        bool valid = body.has("bits");
        if(!body || !valid) return crow::response(400, "Invalid JSON");

        int bits = body["bits"].i();
        if(bits < 2) return crow::response(400, "Number of bits must be greater than 2");

        bigint p = getPrime(bits);
        vector<string> curves = getValidCurves(p);

        // Pack response as JSON
        crow::json::wvalue res;
        res["count"] = curves.size();
        res["prime"] = to_string(p); 
        res["curves"] = vector_to_json(curves);

        return response(res);
    });
    
    
    CROW_ROUTE(app, "/api/rational_points").methods(crow::HTTPMethod::POST)
    ([](const crow::request &req) {
        // Parse JSON input
        auto body = crow::json::load(req.body);
        bool valid = body.has("a") & body.has("b") & body.has("p");
        if(!body || !valid) return crow::response(400, "Invalid JSON");

        bigint a = parse_bigint(body, "a");
        bigint b = parse_bigint(body, "b");
        bigint p = parse_bigint(body, "p");

        vector<point> points = getRationalPoints(a, b, p);
        
        // Pack response as JSON
        crow::json::wvalue res;
        res["count"] = points.size();
        res["points"] = vector_to_json(points);
        
        return response(res);
    });
    
    CROW_ROUTE(app, "/api/add_points").methods(crow::HTTPMethod::POST)
    ([](const crow::request &req) {
         // Parse JSON input
        auto body = crow::json::load(req.body);
        bool valid = body.has("p1") & body.has("q1") & body.has("p");
        if(!body || !valid) return crow::response(400, "Invalid JSON");

        point p1 = parse_point(body, "p1");
        point q1 = parse_point(body, "q1");
        bigint p = parse_bigint(body, "p");

        point r1 = getAddPoints(p1, q1, p);

        // Pack response as JSON
        crow::json::wvalue res;
        res["point"] = serialize(r1);

        return response(res);
    });

    CROW_ROUTE(app, "/api/double_point").methods(crow::HTTPMethod::POST)
    ([](const crow::request &req) {
         // Parse JSON input
        auto body = crow::json::load(req.body);
        bool valid = body.has("p1") & body.has("a") & body.has("p");
        if(!body || !valid) return crow::response(400, "Invalid JSON");

        // Parse input
        point p1 = parse_point(body, "p1");
        bigint a = parse_bigint(body, "a");
        bigint p = parse_bigint(body, "p");

        point p2 = getDoublePoint(p1, a, p);

        // Pack response as JSON
        crow::json::wvalue res;
        res["point"] = serialize(p2); 

        return response(res);
    });

    CROW_ROUTE(app, "/api/hasse").methods(crow::HTTPMethod::POST)
    ([](const crow::request &req) {
         // Parse JSON input
        auto body = crow::json::load(req.body);
        bool valid = body.has("p");
        if(!body || !valid) return crow::response(400, "Invalid JSON");

        // Parse input
        bigint p = parse_bigint(body, "p");

        pair<bigint, bigint> lr = getHassesTheorem(p);

        // Pack response as JSON
        crow::json::wvalue res;
        res["bound"] = serialize(lr);

        return response(res);
    });

    // app.port(18080).multithreaded().run();
    app.port(18080).concurrency(4).run();
}

/*
API test
curl -X POST http://localhost:18080/api/valid_curves -H "Content-Type: application/json" -d '{"p": "3"}'
curl -X POST http://localhost:18080/api/valid_random -H "Content-Type: application/json" -d '{"bits": "3"}'
curl -X POST http://localhost:18080/api/rational_points -H "Content-Type: application/json" -d '{"a": "3", "b": "1", "p": "7"}'
curl -X POST http://localhost:18080/api/add_points -H "Content-Type: application/json" -d '{"p1": {"x": "1", "y": "2", "z": "1"}, "q1": {"x": "3","y": "4","z": "1"}, "p":"7"}'
curl -X POST http://localhost:18080/api/double_point -H "Content-Type: application/json" -d '{"p1": {"x": "1", "y": "2", "z": "1"}, "a": "2", "p":"7"}'
curl -X POST http://localhost:18080/api/hasse -H "Content-Type: application/json" -d '{"p": "7"}'
*/