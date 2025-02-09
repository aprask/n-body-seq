#include <iostream>
#include <cmath>
#include <string>
#include <random>
#include <chrono>
#include <fstream>
#include <ctype.h>
#include <cstring>

using namespace std;
using namespace std::chrono;

#define ARGS 5 // num of particles, time step size, num of iterations, how often to dump state
const double G = 6.674*pow(10,-11); // G for grav force
const double SOFTENING_FACTOR = 0.0001; // this is just an arbitary float (I wasn't sure what "softening factor" actuall was...)

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
        for (int i = 0; i < ARGS-1; i++) {
            for (int j = 0; j < strlen(argv[i]); j++) {
                if (!isdigit(argv[i][j])) {
                    cerr << "Invalid numerical argument type passed: " << argv[i] << endl;
                    return 1;
                }
            }
        }
    }
    const size_t N = stol(argv[1]); // num of particles
    const size_t DELTA_T = stod(argv[2]); // delta t
    const size_t TIME_STEPS = stol(argv[3]); // iterations
    const size_t DUMP_RATE = stol(argv[4]); // how often we should dump the stats

    struct particleNode* particleField = (struct particleNode*)malloc(sizeof(struct particleNode) * N);
    if (!particleField) {
        cerr << "Cannot dynamically allocate mem for particle field" << endl;
        return -1;
    }

    ofstream dataFile;
    dataFile.open("data/results.tsv", ios::app);
    dataFile << N << "\t"; // num of particles for tsv

    for (int i = 0; i < N; ++i) {
        (particleField+i)->velocity.vx = rand() % 100000;
        (particleField+i)->velocity.vx = rand() % 100000;
        (particleField+i)->velocity.vx = rand() % 100000;
        (particleField+i)->force.fx = rand() % 100000;
        (particleField+i)->force.fy = rand() % 100000;
        (particleField+i)->force.fz = rand() % 100000;
        (particleField+i)->position.x = rand() % 100000;
        (particleField+i)->position.y = rand() % 100000;
        (particleField+i)->position.z = rand() % 100000;
        (particleField+i)->mass = rand() % 10000000000;
    }

    // (particleField+0)->velocity.vx = 3;
    // (particleField+0)->velocity.vy = 4;
    // (particleField+0)->velocity.vz = 5;
    // (particleField+0)->force.fx = 6;
    // (particleField+0)->force.fy = 7;
    // (particleField+0)->force.fz = 8;
    // (particleField+0)->position.x = 1;
    // (particleField+0)->position.y = 2;
    // (particleField+0)->position.z = 3;
    // (particleField+0)->mass = 5000;

    // (particleField+1)->velocity.vx = -2;
    // (particleField+1)->velocity.vy = 3;
    // (particleField+1)->velocity.vz = -4;
    // (particleField+1)->force.fx = -5;
    // (particleField+1)->force.fy = 6;
    // (particleField+1)->force.fz = -7;
    // (particleField+1)->position.x = 4;
    // (particleField+1)->position.y = 5;
    // (particleField+1)->position.z = 6;
    // (particleField+1)->mass = 6000;
    
    // auto t_initial = steady_clock::now();
    for (int i = 0; i < N; ++i) {
        cout << "Particle " << i << endl;
        dataFile << (particleField)->mass << "\t";
        for (int j = 0; j < N; ++j) {
            if (i == j) {
                cout << "Skipping particle " << j << endl;
                continue;
            }
            // auto delta_t = duration_cast<seconds>(steady_clock::now() - t_initial).count();

            double totalForce = calculateForce((particleField+i), (particleField+j));
            cout << "Time Step " << i << " total force on particle " << i << " from particle " << j <<  " is: " << totalForce << endl;


            double distance = calculateEuclideanDistance(
                (particleField+i)->position.x,
                (particleField+i)->position.y,
                (particleField+i)->position.z,
                (particleField+j)->position.x,
                (particleField+j)->position.y,
                (particleField+j)->position.z
            );
            cout << "Time Step " << i << " distance in euclidean space: " << distance << endl;

            (particleField+i)->force.fx = calculateForceComponent(
                totalForce,
                (particleField+i)->position.x,
                (particleField+j)->position.x,
                distance
            );
            (particleField+i)->force.fy = calculateForceComponent(
                totalForce, 
                (particleField+i)->position.y,
                (particleField+j)->position.y,
                distance
            );
            (particleField+i)->force.fz = calculateForceComponent(
                totalForce, 
                (particleField+i)->position.z,
                (particleField+j)->position.z,
                distance
            );
            cout << "Time Step " << i << " force component: (" << (particleField+i)->force.fx << "," << (particleField+i)->force.fy << "," << (particleField+i)->force.fz << ")" << endl;

            (particleField+i)->acceleration.ax = calculateAcceleration((particleField+i)->force.fx, (particleField+i)->mass);
            (particleField+i)->acceleration.ay = calculateAcceleration((particleField+i)->force.fy, (particleField+i)->mass);
            (particleField+i)->acceleration.az = calculateAcceleration((particleField+i)->force.fz, (particleField+i)->mass);
            cout << "Time Step " << i << " acceleration: (" << (particleField+i)->acceleration.ax << "," << (particleField+i)->acceleration.ay << "," << (particleField+i)->acceleration.az << ")" << endl;

            (particleField+i)->velocity.vx = calculateVelocity((particleField+i)->velocity.vx, (particleField+i)->acceleration.ax, DELTA_T);
            (particleField+i)->velocity.vy = calculateVelocity((particleField+i)->velocity.vy, (particleField+i)->acceleration.ay, DELTA_T);
            (particleField+i)->velocity.vz = calculateVelocity((particleField+i)->velocity.vz, (particleField+i)->acceleration.az, DELTA_T);
            cout << "Time Step " << i << " velocity: (" << (particleField+i)->velocity.vx << "," << (particleField+i)->velocity.vy << "," << (particleField+i)->velocity.vz << ")" << endl;

            (particleField+i)->position.x = calculatePosition((particleField+i)->position.x, (particleField+i)->velocity.vx, DELTA_T);
            (particleField+i)->position.y = calculatePosition((particleField+i)->position.y, (particleField+i)->velocity.vy, DELTA_T);
            (particleField+i)->position.z = calculatePosition((particleField+i)->position.z, (particleField+i)->velocity.vz, DELTA_T);
            cout << "Time Step " << i << " position: (" << (particleField+i)->position.x << "," << (particleField+i)->position.y << "," << (particleField+i)->position.z << ")" << endl;
            
            if (!(j % DUMP_RATE)) {
                dataFile << (particleField)->position.x << "\t";
                dataFile << (particleField)->position.y << "\t";
                dataFile << (particleField)->position.z << "\t";
                dataFile << (particleField)->velocity.vx << "\t";
                dataFile << (particleField)->velocity.vx << "\t";
                dataFile << (particleField)->velocity.vy << "\t";
                dataFile << (particleField)->force.fx << "\t";
                dataFile << (particleField)->force.fy << "\t";
                dataFile << (particleField)->force.fz << "\t";    
            }
        }
    }
    dataFile << "\n";
    dataFile.close();
    return 0;
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