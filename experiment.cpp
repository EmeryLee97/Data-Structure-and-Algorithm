#include <iostream>
#include <ofstream<
#include <cmath>
#include <eigen>
#include <opencv2/dense>
#include <opencv2/core/eigen.hpp>


Eigen::Matrix3d skewSymmetricMatrix(const Eigen::Vector3d vec) {
    Eigen::Matrix3d mat;
    mat << 0, -vec(2), vec(1),
           vec(2), 0, -vec(0),
           -vec(1), vec(0), 0;
    return mat;
}

// This method requires opencv version > 4.5.2
void outlierRateExperiment(
    const int num_pts=100,
    const std::vector<double> outlier_rate_list={},
    const double noise_level=1,
    const int iter_times=100,
    const double epsilon=0.001,
    const double rho=0.99
    ) {
    
    for (const auto outlier_rate : outlier_rate_list) {
    
        // generate points correspondences
        std::vector<Eigen::Vector3d> pts1_eigen(num_pts); // (3, num_pts), norm of each column is 1
        std::vector<Eigen::Vector3d> pts2_eigen(num_pts); // (3, num_pts), norm of each column is 1
    
        Eigen::Matrix3d R_gt;
        Eigen::Vector3d t_gt;
    
        generateRandomPoints(num_pts, noise_level, outlier_rate, pts1_eigen, pts2_eigen, R_gt, t_gt);
        
        Eigen::Matrix3d t_gt_mat = skewSymmetricMatrix(t_gt);
        Eigen::Matrix3d E_gt = t_gt_mat * R_gt;
        
        // experiment of npt-pose
        
    
        // experiment of gc-ransac
        std::vector<cv::Point2d> pts1_cv(num_pts);
        std::vector<cv::Point2d> pts2_cv(num_pts);
        
        for (int i = 0; i < num_pts; i++) {
            pts1_cv[i] = (cv::Point2d(pts1_eigen[i].x()/pts1_eigen[i].z(), pts1_eigen[i].y()/pts1_eigen[i].z()));
            pts2_cv[i] = (cv::Point2d(pts2_eigen[i].x()/pts2_eigen[i].z(), pts2_eigen[i].y()/pts2_eigen[i].z()));
        }
        cv::Mat K = cv::Mat::eye(3, 3, CV_64F);
        cv::Mat E_gcransac = cv::findEssentialMat(pts1_cv, pts2_cv, K, USAC_ACCURATE, rho, epsilon);
        
        // recover pose from both algorithms
        cv::Mat R_gcransac_cv;
        cv::Mat t_gcransac_cv;
        int ignore = cv::recoverPose(E_gcransac, pts1_cv, pts2_cv, K, R_gcransac_cv, t_gcransac_cv);
        
        Eigen::Marix3d R_gcransac_eigen;
        Eigen::Vector3d t_gcransac_eigen;
        cv::cv2eigen(R_gcransac_cv, R_gcransac_eigen);
        cv::cv2eigen(t_gcransac_cv, t_gcransac_eigen);
        
        double R_error_gcransac = acos(((R_gt.transpose()*R_gcransac_eigen).trace()-1)/2 * 180 / M_PI;
        double t_error_gcransac = acos(abs(t_gt.transpose()*t_gcransac_eigen)) * 180 / M_PI;
        
        // write to txt file
        std::ofstream output_gcransac("gcransac.txt", std::ios::app);
        if (output_gcransac.is_open()) {
            output_gcransac << R_error_gcransac << ' ' << t_error_gcransac << ' ';
        } else {
            std::cerr << "Unable to open the file for writing.\n";
        }
    }
}

void noiseLevelExperiment(
    const int num_pts=100,
    const double outlier_rate=0.2,
    const std::vecor<double> noise_level_list=[],
    const int iter_times=100,
    const double epsilon=0.001,
    const double rho=0.99
    ) {
    
    for (const auto noise_level : noise_level_list) {
    
        // generate points correspondences
        std::vector<Eigen::Vector3d> pts1_eigen(num_pts); // (3, num_pts), norm of each column is 1
        std::vector<Eigen::Vector3d> pts2_eigen(num_pts); // (3, num_pts), norm of each column is 1
    
        Eigen::Matrix3d R_gt;
        Eigen::Vector3d t_gt;
    
        generate2d2dExperiments(num_pts, outlier_rate, noise_level, pts1_eigen, pts2_eigen, R_gt, t_gt);
        
        Eigen::Matrix3d t_gt_mat = skewSymmetricMatrix(t_gt);
        Eigen::Matrix3d E_gt = t_gt_mat * R_gt;
        
        // experiment of npt-pose
        
    
        // experiment of gc-ransac
        std::vector<cv::Point2d> pts1_cv(num_pts);
        std::vector<cv::Point2d> pts2_cv(num_pts);
        
        for (int i = 0; i < num_pts; i++) {
            pts1_cv[i] = (cv::Point2d(pts1_eigen[i].x()/pts1_eigen[i].z(), pts1_eigen[i].y()/pts1_eigen[i].z()));
            pts2_cv[i] = (cv::Point2d(pts2_eigen[i].x()/pts2_eigen[i].z(), pts2_eigen[i].y()/pts2_eigen[i].z()));
        }
        cv::Mat K = cv::Mat::eye(3, 3, CV_64F);
        cv::Mat E_gcransac = cv::findEssentialMat(pts1_cv, pts2_cv, K, USAC_ACCURATE, rho, epsilon);
        
        // recover pose from both algorithms
        cv::Mat R_gcransac_cv;
        cv::Mat t_gcransac_cv;
        int ignore = cv::recoverPose(E_gcransac, pts1_opencv, pts2_opencv, K, R_gcransac_cv, t_gcransac_cv);
        
        Eigen::Marix3d R_gcransac_eigen;
        Eigen::Vector3d t_gcransac_eigen;
        cv::cv2eigen(R_gcransac_cv, R_gcransac_eigen);
        cv::cv2eigen(t_gcransac_cv, t_gcransac_eigen);
        
        double R_error_gcransac = acos(((R_gt.transpose()*R_gcransac_eigen).trace()-1)/2 * 180 / M_PI;
        double t_error_gcransac = acos(abs(t_gt.transpose()*t_gcransac_eigen)) * 180 / M_PI;
        
        // write to txt file
        std::ofstream output_gcransac("gcransac.txt", std::ios::app);
        if (output_gcransac.is_open()) {
            output_gcransac << R_error_gcransac << ' ' << t_error_gcransac << ' ';
        } else {
            std::cerr << "Unable to open the file for writing.\n";
        }
    }
}

void pointsExperiment(
    const std::vector<int> num_pts_list=[],
    const double outlier_rate=0.2,
    const double noise_level_list=0.5,
    const int iter_times=100,
    const double epsilon=0.001,
    const double rho=0.99
    ) {
    
    for (const auto num_pts : num_pts_list) {
    
        // generate points correspondences
        std::vector<Eigen::Vector3d> pts1_eigen(num_pts); // (3, num_pts), norm of each column is 1
        std::vector<Eigen::Vector3d> pts2_eigen(num_pts); // (3, num_pts), norm of each column is 1
    
        Eigen::Matrix3d R_gt;
        Eigen::Vector3d t_gt;
    
        generate2d2dExperiments(num_pts, outlier_rate, noise_level, pts1_eigen, pts2_eigen, R_gt, t_gt);
        
        Eigen::Matrix3d t_gt_mat = skewSymmetricMatrix(t_gt);
        Eigen::Matrix3d E_gt = t_gt_mat * R_gt;
        
        // experiment of npt-pose
        
    
        // experiment of gc-ransac
        std::vector<cv::Point2d> pts1_cv(num_pts);
        std::vector<cv::Point2d> pts2_cv(num_pts);
        
        for (int i = 0; i < num_pts; i++) {
            pts1_cv[i] = (cv::Point2d(pts1_eigen[i].x()/pts1_eigen[i].z(), pts1_eigen[i].y()/pts1_eigen[i].z()));
            pts2_cv[i] = (cv::Point2d(pts2_eigen[i].x()/pts2_eigen[i].z(), pts2_eigen[i].y()/pts2_eigen[i].z()));
        }
        cv::Mat K = cv::Mat::eye(3, 3, CV_64F);
        cv::Mat E_gcransac = cv::findEssentialMat(pts1_cv, pts2_cv, K, USAC_ACCURATE, rho, epsilon);
        
        // recover pose from both algorithms
        cv::Mat R_gcransac_cv;
        cv::Mat t_gcransac_cv;
        int ignore = cv::recoverPose(E_gcransac, pts1_opencv, pts2_opencv, K, R_gcransac_cv, t_gcransac_cv);
        
        Eigen::Marix3d R_gcransac_eigen;
        Eigen::Vector3d t_gcransac_eigen;
        cv::cv2eigen(R_gcransac_cv, R_gcransac_eigen);
        cv::cv2eigen(t_gcransac_cv, t_gcransac_eigen);
        
        double R_error_gcransac = acos(((R_gt.transpose()*R_gcransac_eigen).trace()-1)/2 * 180 / M_PI;
        double t_error_gcransac = acos(abs(t_gt.transpose()*t_gcransac_eigen)) * 180 / M_PI;
        
        // write to txt file
        std::ofstream output_gcransac("gcransac.txt", std::ios::app);
        if (output_gcransac.is_open()) {
            output_gcransac << R_error_gcransac << ' ' << t_error_gcransac << ' ';
        } else {
            std::cerr << "Unable to open the file for writing.\n";
        }
    }
}
