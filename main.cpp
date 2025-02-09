#include <iostream>
#include <cmath>
#include <string>
#include <random>
#include <fstream>
#include <ctype.h>
#include <cstring>

using namespace std;

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
    // if (argc != ARGS) {
    //     cerr << "Invalid number of arguments <" << argc << ">" << ". Expected " << ARGS << endl;
    //     return 1;
    // } else {
    //     for (int i = 0; i < ARGS-1; i++) {
    //         for (int j = 0; j < strlen(argv[i]); j++) {
    //             if (!isdigit(argv[i][j])) {
    //                 cerr << "Invalid numerical argument type passed: " << argv[i] << endl;
    //                 return 1;
    //             }
    //         }
    //     }
    // }
    // const size_t N = stol(argv[1]); // num of particles
    // const size_t DELTA_T = stod(argv[2]); // delta t
    // const size_t TIME_STEPS = stol(argv[3]); // iterations
    // const size_t DUMP_RATE = stol(argv[4]); // how often we should dump the stats

    const size_t N = 2; // Number of particles
    const double DELTA_T = 15; // Time step size
    const size_t TIME_STEPS = 1; // Number of iterations
    const size_t DUMP_RATE = 1; // How often to dump the state


    struct particleNode* particleField = (struct particleNode*)malloc(sizeof(struct particleNode) * N);
    if (!particleField) {
        cerr << "Cannot dynamically allocate mem for particle field" << endl;
        return -1;
    }

    ofstream dataFile;
    dataFile.open("data/results.tsv", ios::app);
    dataFile << N << "\t"; // num of particles for tsv

    // for (int i = 0; i < N; ++i) {
    //     (particleField+i)->velocity.vx = rand() % 100000;
    //     (particleField+i)->velocity.vx = rand() % 100000;
    //     (particleField+i)->velocity.vx = rand() % 100000;
    //     (particleField+i)->force.fx = rand() % 100000;
    //     (particleField+i)->force.fy = rand() % 100000;
    //     (particleField+i)->force.fz = rand() % 100000;
    //     (particleField+i)->position.x = rand() % 100000;
    //     (particleField+i)->position.y = rand() % 100000;
    //     (particleField+i)->position.z = rand() % 100000;
    //     (particleField+i)->mass = rand() % 10000000000;
    // }

    (particleField+0)->velocity.vx = 3;
    (particleField+0)->velocity.vy = 4;
    (particleField+0)->velocity.vz = 5;
    (particleField+0)->force.fx = 6;
    (particleField+0)->force.fy = 7;
    (particleField+0)->force.fz = 8;
    (particleField+0)->position.x = 1;
    (particleField+0)->position.y = 2;
    (particleField+0)->position.z = 3;
    (particleField+0)->mass = 5000;

    (particleField+1)->velocity.vx = -2;
    (particleField+1)->velocity.vy = 3;
    (particleField+1)->velocity.vz = -4;
    (particleField+1)->force.fx = -5;
    (particleField+1)->force.fy = 6;
    (particleField+1)->force.fz = -7;
    (particleField+1)->position.x = 4;
    (particleField+1)->position.y = 5;
    (particleField+1)->position.z = 6;
    (particleField+1)->mass = 6000;

    
    for (int i = 0; i < TIME_STEPS; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << "Particle " << j << endl;
            dataFile << (particleField)->mass << "\t";
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
                    (particleField+j)->position.x,
                    distance
                );
                (particleField+j)->force.fy = calculateForceComponent(
                    totalForce, 
                    (particleField+j)->position.y,
                    (particleField+j)->position.y,
                    distance
                );
                (particleField+j)->force.fz = calculateForceComponent(
                    totalForce, 
                    (particleField+j)->position.z,
                    (particleField+j)->position.z,
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
                    dataFile << (particleField+j)->position.x << "\t";
                    dataFile << (particleField+j)->position.y << "\t";
                    dataFile << (particleField+j)->position.z << "\t";
                    dataFile << (particleField+j)->velocity.vx << "\t";
                    dataFile << (particleField+j)->velocity.vx << "\t";
                    dataFile << (particleField+j)->velocity.vy << "\t";
                    dataFile << (particleField+j)->force.fx << "\t";
                    dataFile << (particleField+j)->force.fy << "\t";
                    dataFile << (particleField+j)->force.fz << "\t";    
                }
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