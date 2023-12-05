#include <iostream>
#include <chrono>
#include "construct_coeff_constraint.h"
#include "create2D2DExperiment.h"
#include "essential_matrix.h"
#include "npt_pose.h"

int main()
{
    // parameters for synthetic scene
    int pts_num = 100;
    double noise_level = 1.0;
    double outlier_rate = 0;

    std::vector<Eigen::Vector3d> pts_1(pts_num);
    std::vector<Eigen::Vector3d> pts_2(pts_num);
    Eigen::Matrix3d R_gt, E_gt;
    Eigen::Vector3d T_gt;
    double* P1 = new double[pt_number*3];
    double* P2 = new double[pt_number*3];

    std::vector<Eigen::Matrix<double, 12, 12>> A;
    std::vector<double> b;
    double* C = new double[81];
    Eigen::Matrix<double, 12, 12> X_sol;
    Eigen::Matrix3d E_est;

    int iter_times = 100;
    for (int i = 0; i < iter_times; i++)
    {
        std::cout << "---------------------------------------" << std::endl;
        std::cout << "i: " << i << std::endl;
        // prepare data
        generateRandomPoints(pts_num, noise_level, outlier_rate, pts_1, pts_2, R_gt, T_gt);
        construct_essential_matrix(R_gt, T_gt, E_gt);
        for (int i = 0; i < pts_num; i++)
        {
            P1[i*3] = pts_1[i](0);
            P1[i*3+1] = pts_1[i](1);
            P1[i*3+2] = pts_1[i](2);
            P2[i*3] = pts_2[i](0);
            P2[i*3+1] = pts_2[i](1);
            P2[i*3+2] = pts_2[i](2);
        }

        // do the job
        // timer begin
        std::chrono::time_point<std::chrono::system_clock> time_start = std::chrono::system_clock::now();
        npt_pose(P1, P2, C, pts_num, X_sol, E_est, true);
        //timer end
        std::chrono::time_point<std::chrono::system_clock> time_end = std::chrono::system_clock::now();
        std::chrono::duration<double> time_spend = time_end - time_start;
        std::cout << "Runtime: " << time_spend.count()*1000.0 << " ms" << std::endl;

        // show the results
        if (1)
        {
            std::cout << "ground truth" << std::endl;
//            std::cout << "R = " << std::endl << R_gt << std::endl;
//            std::cout << "T = " << std::endl << T_gt.transpose() << std::endl;
            std::cout << "E = " << std::endl << E_gt << std::endl;
            std::cout << "optimization" << std::endl;
            std::cout << "E = " << std::endl << E_est << endl;
        }
    }
    delete[] P1;
    delete[] P2;
    delete[] C;

    return 1;
}

