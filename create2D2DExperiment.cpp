#include <random>
#include <vector>
#include "create2D2DExperiment.h"

// converted from some functions in openGV
// https://github.com/laurentkneip/opengv/tree/master/matlab/helpers

void generateRandomR(Eigen::Matrix3d& R) {
    // random generate rotation axis
    Eigen::Vector3d rot_axis = Eigen::Vector3d::Random(); // [-1, 1]
    rot_axis.normalize();
    
    // random generate rotation angle
    double rot_angle = M_PI * 2.0 * ((double)rand() / (RAND_MAX) - 0.5); // [-pi, pi]
    
    // generate rotation matrix
    Eigen::AngleAxisd rotation(rot_angle, rot_axis);
    R = rotation.toRotationMatrix();
    
    return;
}


void generateRandomPoints(
    int pts_num, double noise_level, double outlier_rate,
    std::vector<Eigen::Vector3d>& pts_1, // points in the first image, size = pts_num
    std::vector<Eigen::Vector3d>& pts_2, // points in the second image, size = pts_num
    Eigen::Matrix3d& R_gt, // rotation matrix, ground truth
    Eigen::Vector3d& t_gt // translation matrix, ground truth
    ) {
    
    // generate R
    generateRandomR(R_gt);
    
    // generate pts_1
    int sum0 = 0, sum1 = 0, sum2 = 0;
    for (int i = 0; i < pts_num; i++) {
        pts_1[i] = Eigen::Vector3d::Random(); // [-1, 1]
        pts_1[i](2) += 2; // [1, 3]
        
        sum0 += pts_1[i](0);
        sum1 += pts_1[i](1);
        sum2 += pts_1[i](2);
     }
    
     // generate t
     Eigen::Vector3d pts_1_c(sum0/pts_num, sum1/pts_num, sum2/pts_num);
     Eigen::Vector3d random_vec = Eigen::Vector3d::Random(); // [-1, 1]
     Eigen::Vector3d scaled_vec = 0.5 * (random_vec.array() + 1); // [0, 1]
     t_gt = -R_gt * pts_1_c + pts_1_c + scaled_vec;
     
     // generate pts_2
     for (int i = 0; i < pts_num; i++) {
        pts_2[i] = R_gt * pts_1[i] + t_gt;
     }
     
     int outlier_num = pts_num * outlier_rate;
     int inlier_num = pts_num - outlier_num;
     
     double focal_length = 1000;
     
     // add gaussian noise to pts_2
     std::random_device rd;
     std::mt19937 gen(rd());
     std::normal_distribution<double> distribution_inlier(0, noise_level); // inlier
     std::normal_distribution<double> distribution_outlier(0, 20); // outlier
     
     for (int i = 0; i < inlier_num; i++) {
        pts_2[i](0) = focal_length * pts_2[i](0) / pts_2[i](2) + distribution_inlier(gen);
        pts_2[i](1) = focal_length * pts_2[i](1) / pts_2[i](2) + distribution_inlier(gen);
        pts_2[i](2) = focal_length;
     }
     
     for (int i = inlier_num; i < pts_num; i++) {
        pts_2[i](0) = focal_length * pts_2[i](0) / pts_2[i](2) + distribution_outlier(gen);
        pts_2[i](1) = focal_length * pts_2[i](1) / pts_2[i](2) + distribution_outlier(gen);
        pts_2[i](2) = focal_length;
     }
     
     // normalization
     t_gt.normalize();
     
     for (int i = 0; i < pts_num; i++) {
        pts_1[i].normalize();
        pts_2[i].normalize();
     }
}

