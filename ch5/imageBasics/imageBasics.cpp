#include <iostream>
#include <chrono>

// OpenCV 4 headers
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgcodecs.hpp>

int main(int argc, char **argv) {
  // Read the image specified by argv[1]
  cv::Mat image = cv::imread(argv[1]);

  // Check if the image file has been read correctly
  if (image.empty()) {
    std::cerr << "The file " << argv[1] << " does not exist." << std::endl;
    return 0;
  }

  // Output basic information about the image
  std::cout << "Image width: " << image.cols << ", height: " << image.rows 
            << ", number of channels: " << image.channels() << std::endl;
  cv::imshow("image", image);
  cv::waitKey(0);

  // Check the type of image
  if (image.type() != CV_8UC1 && image.type() != CV_8UC3) {
    std::cout << "Please enter a color or grayscale image." << std::endl;
    return 0;
  }

  // Time the pixel access operation
  auto t1 = std::chrono::steady_clock::now();
  for (int y = 0; y < image.rows; ++y) {
    auto *row_ptr = image.ptr<unsigned char>(y);
    for (int x = 0; x < image.cols; ++x) {
      auto *data_ptr = &row_ptr[x * image.channels()];
      for (int c = 0; c < image.channels(); ++c) {
        unsigned char data = data_ptr[c];
      }
    }
  }
  auto t2 = std::chrono::steady_clock::now();
  auto time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  std::cout << "Time taken to traverse the image: " << time_used.count() << " seconds." << std::endl;

  // Direct assignment does not copy data
  cv::Mat image_another = image;
  image_another(cv::Rect(0, 0, 100, 100)).setTo(0);
  cv::imshow("image", image);
  cv::waitKey(0);

  // Use clone to make a deep copy
  cv::Mat image_clone = image.clone();
  image_clone(cv::Rect(0, 0, 100, 100)).setTo(255);
  cv::imshow("image", image);
  cv::imshow("image_clone", image_clone);
  cv::waitKey(0);

  cv::destroyAllWindows();
  return 0;
}
