#include <iostream>
#include <cmath>
#include <string>
#include <random>
#include <fstream>
#include <ctype.h>
#include <cstring>
#include <filesystem>
#include <sstream>
#include <chrono>

#define ARGS 5

const double G = 6.674*pow(10,-11); // G for grav force
const double SOFTENING_FACTOR = 0.001; // this is just an arbitary float (I wasn't sure what "softening factor" actuall was...)

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
    bool readFromFile = false;
    size_t N = 0;
    struct particleNode* particleField;
    if (argc != ARGS) {
        std::cerr << "Invalid number of arguments <" << argc << ">" << ". Expected " << ARGS << std::endl;
        return 1;
    } else {
        for (int i = 1; i < ARGS; i++) {
            for (int j = 0; j < strlen(argv[i]); j++) {
                if (!isdigit(argv[i][j])) {
                    if (i == 1) {
                        std::string file = argv[i];
                        file = file += ".tsv";
                        std::ifstream sourceFile(file);
                        if (!sourceFile) {
                            std::cerr << "Failed to open file: " << argv[i] << std::endl;
                            return 1;
                        }
                        if (sourceFile.is_open()) {
                            readFromFile = true;
                            std::string line;
                            while(getline(sourceFile, line)) {
                                N++;
                            }
                            particleField = (struct particleNode*)malloc(sizeof(struct particleNode) * N);
                            if (!particleField) {
                                std::cerr << "Cannot dynamically allocate mem for particle field" << std::endl;
                                return -1;
                            }
                            sourceFile.close();
                            std::ifstream sourceFile(file);
                            int i = 0;
                            while (getline(sourceFile, line)) {
                                std::stringstream s(line);
                                std::string token;
                                int attrib = 10; // num of attrib (mass, pos, vel, force)
                                for (int tokenIdx = 1; getline(s, token, '\t'); ++tokenIdx) {
                                    switch (tokenIdx) {
                                        case 1:
                                            (particleField+i)->mass = stod(token);
                                            break; // mass
                                        case 2:
                                            (particleField+i)->position.x = stod(token);
                                            break; // x
                                        case 3:
                                            (particleField+i)->position.y = stod(token);
                                            break; // y
                                        case 4:
                                            (particleField+i)->position.z = stod(token);
                                            break; // z
                                        case 5:
                                            (particleField+i)->velocity.vx = stod(token);
                                            break; // vx
                                        case 6:
                                            (particleField+i)->velocity.vy = stod(token);
                                            break; // vy
                                        case 7:
                                            (particleField+i)->velocity.vz = stod(token);
                                            break; // vz
                                        case 8:
                                            (particleField+i)->force.fx = stod(token);
                                            break; // fx
                                        case 9:
                                            (particleField+i)->force.fy = stod(token);
                                            break; // fy
                                        case 10:
                                            (particleField+i)->force.fz = stod(token);
                                            break; // fz
                                        default:
                                            break;
                                            
                                    }
                                }
                                i++; 
                        }
                        sourceFile.close();
                        break;
                    }
                    } else {
                        std::cerr << "Invalid numerical argument type passed: " << argv[i] << std::endl;
                        return 1;
                    }
                }
            }
        }
    }
    if (!readFromFile) {
        N = std::stol(argv[1]);
        particleField = (struct particleNode*)malloc(sizeof(struct particleNode) * N);
        for (int i = 0; i < N; ++i) {
            (particleField+i)->velocity.vx = rand() % 100000;
            (particleField+i)->velocity.vy = rand() % 100000;
            (particleField+i)->velocity.vz = rand() % 100000;
            (particleField+i)->position.x = rand() % 100000;
            (particleField+i)->position.y = rand() % 100000;
            (particleField+i)->position.z = rand() % 100000;
            (particleField+i)->mass = rand() % 10000000000;
        }
    }
    if (!particleField) {
        std::cerr << "Cannot dynamically allocate mem for particle field" << std::endl;
        return -1;
    }
    size_t DELTA_T = std::stod(argv[2]);
    size_t TIME_STEPS = std::stol(argv[3]);
    size_t DUMP_RATE = std::stol(argv[4]);
    
    std::ofstream dataFile;
    dataFile.open("results.tsv", std::ios::trunc);
    auto global_time = std::chrono::steady_clock::now();
    for (int i = 0; i < TIME_STEPS; ++i) {
        for (int k = 0; k < N; ++k) {
            (particleField+k)->force.fx = 0;
            (particleField+k)->force.fy = 0;
            (particleField+k)->force.fz = 0;
        }
        for (int j = 0; j < N; ++j) {
           for (int k = 0; k < N; ++k) {
                if (k == j) continue;
                double totalForce = calculateForce((particleField+j), (particleField+k));
                double distance = calculateEuclideanDistance(
                    (particleField+j)->position.x,
                    (particleField+j)->position.y,
                    (particleField+j)->position.z,
                    (particleField+k)->position.x,
                    (particleField+k)->position.y,
                    (particleField+k)->position.z
                );    
                (particleField+j)->force.fx += calculateForceComponent(
                    totalForce,
                    (particleField+j)->position.x,
                    (particleField+k)->position.x,
                    distance
                );
                (particleField+j)->force.fy += calculateForceComponent(
                    totalForce, 
                    (particleField+j)->position.y,
                    (particleField+k)->position.y,
                    distance
                );
                (particleField+j)->force.fz += calculateForceComponent(
                    totalForce, 
                    (particleField+j)->position.z,
                    (particleField+k)->position.z,
                    distance
                );
    
                (particleField+j)->acceleration.ax = calculateAcceleration((particleField+j)->force.fx, (particleField+j)->mass);
                (particleField+j)->acceleration.ay = calculateAcceleration((particleField+j)->force.fy, (particleField+j)->mass);
                (particleField+j)->acceleration.az = calculateAcceleration((particleField+j)->force.fz, (particleField+j)->mass);
    
                (particleField+j)->velocity.vx = calculateVelocity((particleField+j)->velocity.vx, (particleField+j)->acceleration.ax, DELTA_T);
                (particleField+j)->velocity.vy = calculateVelocity((particleField+j)->velocity.vy, (particleField+j)->acceleration.ay, DELTA_T);
                (particleField+j)->velocity.vz = calculateVelocity((particleField+j)->velocity.vz, (particleField+j)->acceleration.az, DELTA_T);
    
                (particleField+j)->position.x = calculatePosition((particleField+j)->position.x, (particleField+j)->velocity.vx, DELTA_T);
                (particleField+j)->position.y = calculatePosition((particleField+j)->position.y, (particleField+j)->velocity.vy, DELTA_T);
                (particleField+j)->position.z = calculatePosition((particleField+j)->position.z, (particleField+j)->velocity.vz, DELTA_T);
                
                if (!(i % DUMP_RATE)) {
                    dataFile << N << "\t";
                    for (int l = 0; l < N; ++l) {
                        dataFile << (particleField+l)->mass << "\t";
                        dataFile << (particleField+l)->position.x << "\t";
                        dataFile << (particleField+l)->position.y << "\t";
                        dataFile << (particleField+l)->position.z << "\t";
                        dataFile << (particleField+l)->velocity.vx << "\t";
                        dataFile << (particleField+l)->velocity.vy << "\t";
                        dataFile << (particleField+l)->velocity.vz << "\t";
                        dataFile << (particleField+l)->force.fx << "\t";
                        dataFile << (particleField+l)->force.fy << "\t";
                        dataFile << (particleField+l)->force.fz << "\t";
                    }
                    dataFile << std::endl;
                }
            }
        }
    }
    dataFile.close();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - global_time).count();
    std::ofstream logFile;
    logFile.open("data.csv", std::ios::app);
    logFile << N << "," << DELTA_T << "," << TIME_STEPS << "," << duration << std::endl;
    free(particleField);
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
    double distance = calculateEuclideanDistance(
        p1->position.x,
        p1->position.y,
        p1->position.z,
        p2->position.x,
        p2->position.y,
        p2->position.z
    );
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