#include <iostream>
#include "LatticeModel.h"
#include "Simulation.h"
#include "OutputFile.h"
#include <mpi.h>
#include "Subdomain.h"
#include "Input.h"
#include "ComplexGeometry.h"
#include "disclaimer.h"
#include <iostream>
#include <chrono>
#include <thread>
#include <iomanip>
#include <omp.h>
#include "Utility.h"
#include "GeometryReader.h"
#include "Input.h"
#include <cmath>  

int main(int argc, char* argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    int myrank, num_procs;
    Utility Utility;
    Reader Reader;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (myrank == 0) {
        Print_disclaimer();

        omp_set_num_threads(omp_get_max_threads());
        int numThreads = omp_get_max_threads();
        omp_set_num_threads(numThreads);
        bool stopChrono = false;

        std::thread chronoThread([&]() {
            while (!stopChrono) {
                auto now = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = now - start;
                auto hours = std::chrono::duration_cast<std::chrono::hours>(elapsed).count();
                auto minutes = std::chrono::duration_cast<std::chrono::minutes>(elapsed).count() % 60;
                auto seconds = std::chrono::duration_cast<std::chrono::seconds>(elapsed).count() % 60;
                auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() % 1000;

                std::cout << "\rElapsed time: "
                    << std::setw(2) << std::setfill('0') << hours << ":"
                    << std::setw(2) << std::setfill('0') << minutes << ":"
                    << std::setw(2) << std::setfill('0') << seconds << "."
                    << std::setw(3) << std::setfill('0') << milliseconds << std::flush;

                std::this_thread::sleep_for(std::chrono::milliseconds(100));  // Update every 100 ms
            }
            });
        chronoThread.detach();
    }

    double nu_blood = 3.5e-6;           
    double dt = 8e-6;        

    double cs = spacing*1e-3 / (sqrt(3) * dt);
  //  double tau = nu_blood / (cs * cs * dt) + 0.5;
    double tau=1.5;

    if (myrank == 0) {
        std::cout << "\nRelaxation Time (tau): " << tau << std::endl;
}

    std::shared_ptr<LatticeModel> velSet = std::make_shared<LatticeModel>();
    Geometry Geometry("PV.stl", velSet, spacing, 'z');
    SubDomain SubDomain(Geometry);
    Simulation(SubDomain, tau, pressure, dt);  

    MPI_Finalize();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    if (myrank == 0) {
        std::cout << "\nTotal Elapsed time: " << std::fixed << std::setprecision(2)
            << std::setw(8) << elapsed.count() << " seconds." << std::endl;
    }

    return 0;
}
