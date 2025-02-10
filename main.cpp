#include <iostream>
#include <cmath>
#include <string>
#include <random>
#include <fstream>
#include <ctype.h>
#include <cstring>
#include <filesystem>

using namespace std;
using namespace std::filesystem;

#define ARGS 6
#define SOLAR_N 3

const double G = 6.674*pow(10,-11); // G for grav force
const double SOFTENING_FACTOR = 0.0001; // this is just an arbitary float (I wasn't sure what "softening factor" actuall was...)

void particleFieldInit(int flag, struct particleNode* p, const int size);
double calculateForceComponent(double mass1, double coordinate1, double mass2, double coordinate2);
double calculateAcceleration(double force, double mass);
double calculatePosition(double coordinate, double velocity, double delta_t);
double calculateVelocity(double velocity, double acceleration, double delta_t);
double calculateLinearDistance(double coordinate1, double coordinate2);
double calculateEuclideanDistance(double x1, double y1, double z1, double x2, double y2, double z2);
double calculateForce(struct particleNode* p1, struct particleNode* p2);

struct acceleration {
    double ax;
    double ay;
    double az;
};

struct position {
    double x;
    double y;
    double z;
};

struct velocity {
    double vx;
    double vy;
    double vz;
};

struct force {
    double fx;
    double fy;
    double fz;
};

struct particleNode {
    double mass;
    struct position position;
    struct velocity velocity;
    struct force force;
    struct acceleration acceleration;
};


int main (int argc, char* argv[]) {
    if (argc != ARGS) {
        cerr << "Invalid number of arguments <" << argc << ">" << ". Expected " << ARGS << endl;
        return 1;
    } else {
        for (int i = 1; i < ARGS-1; i++) {
            for (int j = 0; j < strlen(argv[i]); j++) {
                if (!isdigit(argv[i][j])) {
                    cerr << "Invalid numerical argument type passed: " << argv[i] << endl;
                    return 1;
                }
            }
        }
    }
    size_t FLAG = stoi(argv[1]);
    size_t DELTA_T = stod(argv[3]); // delta t
    size_t TIME_STEPS = stol(argv[4]); // iterations
    size_t DUMP_RATE = stol(argv[5]); // how often we should dump the stats
    size_t N;
    if (FLAG == 2) {
        N = SOLAR_N;
    } else {
        N = stol(argv[2]);
    }
    struct particleNode* particleField = (struct particleNode*)malloc(sizeof(struct particleNode) * N);
    if (!particleField) {
        cerr << "Cannot dynamically allocate mem for particle field" << endl;
        return -1;
    }

    particleFieldInit(FLAG, particleField, N);

    ofstream dataFile;
    dataFile.open("results.tsv", ios::trunc);
    for (int i = 0; i < TIME_STEPS; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << "Particle " << j << endl;
            for (int k = 0; k < N; ++k) {
                if (k == j) {
                    cout << "Skipping particle " << k << endl;
                    continue;
                }
    
                double totalForce = calculateForce((particleField+j), (particleField+k));
                cout << "Time Step " << j << " total force on particle " << j << " from particle " << j <<  " is: " << totalForce << endl;
    
    
                double distance = calculateEuclideanDistance(
                    (particleField+j)->position.x,
                    (particleField+j)->position.y,
                    (particleField+j)->position.z,
                    (particleField+k)->position.x,
                    (particleField+k)->position.y,
                    (particleField+k)->position.z
                );
                cout << "Time Step " << j << " distance in euclidean space: " << distance << endl;
    
                (particleField+j)->force.fx = calculateForceComponent(
                    totalForce,
                    (particleField+j)->position.x,
                    (particleField+k)->position.x,
                    distance
                );
                (particleField+j)->force.fy = calculateForceComponent(
                    totalForce, 
                    (particleField+j)->position.y,
                    (particleField+k)->position.y,
                    distance
                );
                (particleField+j)->force.fz = calculateForceComponent(
                    totalForce, 
                    (particleField+j)->position.z,
                    (particleField+k)->position.z,
                    distance
                );
                cout << "Time Step " << j << " force component: (" << (particleField+j)->force.fx << "," << (particleField+j)->force.fy << "," << (particleField+j)->force.fz << ")" << endl;
    
                (particleField+j)->acceleration.ax = calculateAcceleration((particleField+j)->force.fx, (particleField+j)->mass);
                (particleField+j)->acceleration.ay = calculateAcceleration((particleField+j)->force.fy, (particleField+j)->mass);
                (particleField+j)->acceleration.az = calculateAcceleration((particleField+j)->force.fz, (particleField+j)->mass);
                cout << "Time Step " << j << " acceleration: (" << (particleField+j)->acceleration.ax << "," << (particleField+j)->acceleration.ay << "," << (particleField+j)->acceleration.az << ")" << endl;
    
                (particleField+j)->velocity.vx = calculateVelocity((particleField+j)->velocity.vx, (particleField+j)->acceleration.ax, DELTA_T);
                (particleField+j)->velocity.vy = calculateVelocity((particleField+j)->velocity.vy, (particleField+j)->acceleration.ay, DELTA_T);
                (particleField+j)->velocity.vz = calculateVelocity((particleField+j)->velocity.vz, (particleField+j)->acceleration.az, DELTA_T);
                cout << "Time Step " << j << " velocity: (" << (particleField+j)->velocity.vx << "," << (particleField+j)->velocity.vy << "," << (particleField+j)->velocity.vz << ")" << endl;
    
                (particleField+j)->position.x = calculatePosition((particleField+j)->position.x, (particleField+j)->velocity.vx, DELTA_T);
                (particleField+j)->position.y = calculatePosition((particleField+j)->position.y, (particleField+j)->velocity.vy, DELTA_T);
                (particleField+j)->position.z = calculatePosition((particleField+j)->position.z, (particleField+j)->velocity.vz, DELTA_T);
                cout << "Time Step " << j << " position: (" << (particleField+j)->position.x << "," << (particleField+j)->position.y << "," << (particleField+j)->position.z << ")" << endl;
                
                if (!(k % DUMP_RATE)) {
                    dataFile << N << "\t";
                    dataFile << (particleField+j)->mass << "\t";
                    dataFile << (particleField+j)->position.x << "\t";
                    dataFile << (particleField+j)->position.y << "\t";
                    dataFile << (particleField+j)->position.z << "\t";
                    dataFile << (particleField+j)->velocity.vx << "\t";
                    dataFile << (particleField+j)->velocity.vx << "\t";
                    dataFile << (particleField+j)->velocity.vy << "\t";
                    dataFile << (particleField+j)->force.fx << "\t";
                    dataFile << (particleField+j)->force.fy << "\t";
                    dataFile << (particleField+j)->force.fz << "\t";
                    dataFile << "\n";
                }
            }
        }
    }
    dataFile.close();
    return 0;
}

void particleFieldInit(int flag, struct particleNode* p, const int size) {

    switch (flag) {
        case 1:
            for (int i = 0; i < size; ++i) {
                (p+i)->velocity.vx = rand() % 100000;
                (p+i)->velocity.vy = rand() % 100000;
                (p+i)->velocity.vx = rand() % 100000;
                (p+i)->force.fx = rand() % 100000;
                (p+i)->force.fy = rand() % 100000;
                (p+i)->force.fz = rand() % 100000;
                (p+i)->position.x = rand() % 100000;
                (p+i)->position.y = rand() % 100000;
                (p+i)->position.z = rand() % 100000;
                (p+i)->mass = rand() % 10000000000;
            }
            break;
        case 2: // 3 body pre-defined model
            (p+0)->velocity.vx = 225;
            (p+0)->velocity.vy = 0.0; // sun (aka the origin)
            (p+0)->velocity.vz = 0.0;
            (p+0)->force.fx = 0;
            (p+0)->force.fy = 0;
            (p+0)->force.fz = 0;
            (p+0)->position.x = 0;
            (p+0)->position.y = 0;
            (p+0)->position.z = 0;
            (p+0)->mass = 1.9891*10e30;

            (p+1)->velocity.vx = 29.8;
            (p+1)->velocity.vy = 0.0;
            (p+1)->velocity.vz = 0.0; // earth
            (p+1)->force.fx = 0;
            (p+1)->force.fy = 0;
            (p+1)->force.fz = 0;
            (p+1)->position.x = 149.6*10e6;
            (p+1)->position.y = 0;
            (p+1)->position.z = 0;
            (p+1)->mass = 5.97219*1e24;

            (p+2)->velocity.vx = 0.0549; // moon
            (p+2)->velocity.vy = 0.0;
            (p+2)->velocity.vz = 0.0;
            (p+2)->force.fx = 0;
            (p+2)->force.fy = 0;
            (p+2)->force.fz = 0;
            (p+2)->position.x = 0.384*10e6;
            (p+2)->position.y = 0;
            (p+2)->position.z = 0;
            (p+2)->mass = 7.34767309*10e22;
            break;
        case 3:
            // ifstream loadedFile("solar.tsv");
            break;
        default:
            cout << "Invalid init type: " << flag << endl;
            exit(1);
    }
}

double calculateEuclideanDistance(double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt((pow((x2-x1), 2) + pow((y2-y1), 2) + pow((z2-z1), 2)));
}

double calculateLinearDistance(double coordinate1, double coordinate2) {
    return coordinate2-coordinate1;
}

double calculateForce(struct particleNode* p1, struct particleNode* p2) {
    double totalMass = p1->mass*p2->mass;
    cout << "Total Mass: " << totalMass << endl;
    double distance = calculateEuclideanDistance(
        p1->position.x,
        p1->position.y,
        p1->position.z,
        p2->position.x,
        p2->position.y,
        p2->position.z
    );
    cout << "Distance: " << distance << endl;
    return (distance/abs(distance))*G*(totalMass/(pow(distance,2)+SOFTENING_FACTOR));
}

double calculatePosition(double coordinate, double velocity, double delta_t) {
    return coordinate + velocity*delta_t;
}

double calculateVelocity(double velocity, double acceleration, double delta_t) {
    return velocity + acceleration*delta_t;
}

double calculateAcceleration(double force, double mass) {
    return force/mass;
}

double calculateForceComponent(double force, double coordinate1, double coordinate2, double distance) {
    return force*((coordinate2-coordinate1)/(distance));
}