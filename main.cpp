#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <numeric>
#include <stdexcept>

// Constants for BMP file
const int BYTES_PER_PIXEL = 3; // red, green, & blue
const int FILE_HEADER_SIZE = 14;
const int INFO_HEADER_SIZE = 40;

// Constants for data processing
const int SIGNAL_THRESHOLD = 2048;
const int BETA_THORNE_ROWS_COUNT = 15;
const int MEDIAN_DETECTOR_COUNT = 50;
const double THICKNESS_CALIBRATION_FACTOR = 0.247;

void generateBitmapImage(const unsigned char* image, int height, int width, const char* imageFileName);
unsigned char* createBitmapFileHeader(int height, int stride);
unsigned char* createBitmapInfoHeader(int height, int width);

// Struct to hold pixel data and calibration status
struct PixelData {
    double value;
    bool is_calibrated;
};
void generateBitmapImage(const unsigned char* image, int height, int width, const char* imageFileName)
{
    int widthInBytes = width * BYTES_PER_PIXEL;

    unsigned char padding[3] = {0, 0, 0};
    int paddingSize = (4 - (widthInBytes) % 4) % 4;

    int stride = (widthInBytes) + paddingSize;

    FILE* imageFile = fopen(imageFileName, "wb");

    unsigned char* fileHeader = createBitmapFileHeader(height, stride);
    fwrite(fileHeader, 1, FILE_HEADER_SIZE, imageFile);

    unsigned char* infoHeader = createBitmapInfoHeader(height, width);
    fwrite(infoHeader, 1, INFO_HEADER_SIZE, imageFile);

    int i;
    for (i = 0; i < height; i++) {
        fwrite(image + (i*widthInBytes), BYTES_PER_PIXEL, width, imageFile);
        fwrite(padding, 1, paddingSize, imageFile);
    }

    fclose(imageFile);
}

unsigned char* createBitmapFileHeader(int height, int stride)
{
    int fileSize = FILE_HEADER_SIZE + INFO_HEADER_SIZE + (stride * height);

    static unsigned char fileHeader[] = {
        0,0,     /// signature
        0,0,0,0, /// image file size in bytes
        0,0,0,0, /// reserved
        0,0,0,0, /// start of pixel array
    };

    fileHeader[ 0] = (unsigned char)('B');
    fileHeader[ 1] = (unsigned char)('M');
    fileHeader[ 2] = (unsigned char)(fileSize      );
    fileHeader[ 3] = (unsigned char)(fileSize >>  8);
    fileHeader[ 4] = (unsigned char)(fileSize >> 16);
    fileHeader[ 5] = (unsigned char)(fileSize >> 24);
    fileHeader[10] = (unsigned char)(FILE_HEADER_SIZE + INFO_HEADER_SIZE);

    return fileHeader;
}

unsigned char* createBitmapInfoHeader(int height, int width)
{
    static unsigned char infoHeader[] = {
        0,0,0,0, /// header size
        0,0,0,0, /// image width
        0,0,0,0, /// image height
        0,0,     /// number of color planes
        0,0,     /// bits per pixel
        0,0,0,0, /// compression
        0,0,0,0, /// image size
        0,0,0,0, /// horizontal resolution
        0,0,0,0, /// vertical resolution
        0,0,0,0, /// colors in color table
        0,0,0,0, /// important color count
    };

    infoHeader[ 0] = (unsigned char)(INFO_HEADER_SIZE);
    infoHeader[ 4] = (unsigned char)(width      );
    infoHeader[ 5] = (unsigned char)(width >>  8);
    infoHeader[ 6] = (unsigned char)(width >> 16);
    infoHeader[ 7] = (unsigned char)(width >> 24);
    infoHeader[ 8] = (unsigned char)(height      );
    infoHeader[ 9] = (unsigned char)(height >>  8);
    infoHeader[10] = (unsigned char)(height >> 16);
    infoHeader[11] = (unsigned char)(height >> 24);
    infoHeader[12] = (unsigned char)(1);
    infoHeader[14] = (unsigned char)(BYTES_PER_PIXEL*8);

    return infoHeader;
}

// Function to read data from file
std::vector<std::vector<int>> read_data_from_file(const std::string& filename, unsigned& height, unsigned& width) {
    std::ifstream inf(filename, std::fstream::in | std::fstream::binary);
    if (!inf.is_open()) {
        throw std::runtime_error("Error: could not open file " + filename);
    }

    inf.read(reinterpret_cast<char*>(&width), sizeof(unsigned));
    inf.read(reinterpret_cast<char*>(&height), sizeof(unsigned));
    
    // Skip the next 14 unsigned integers (header info)
    inf.seekg(static_cast<unsigned>(inf.tellg()) + sizeof(unsigned) * 14);

    if (height == 0 || width == 0) {
        throw std::runtime_error("Error: image dimensions cannot be zero.");
    }
    
    std::vector<std::vector<int>> data(height, std::vector<int>(width));
    unsigned number;
    for (unsigned i = 0; i < height; ++i) {
        for (unsigned j = 0; j < width; ++j) {
            inf.read(reinterpret_cast<char*>(&number), 4);
            data[i][j] = static_cast<int>(number);
        }
    }
    inf.close();
    return data;
}

// function to process data
std::vector<std::vector<PixelData>> process_data(const std::vector<std::vector<int>>& data) {
    unsigned m = data.size();
    unsigned n = data[0].size();
    std::vector<std::vector<PixelData>> processed_data(m, std::vector<PixelData>(n));
    // Background normalization
    for (unsigned i = 0; i < m; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            if (data[i][j] > SIGNAL_THRESHOLD) {
                processed_data[i][j].value = data[i][j] - SIGNAL_THRESHOLD;
            } else {
                processed_data[i][j].value = 0;
            }
            processed_data[i][j].is_calibrated = false;
        }
    }

    // Calibration by beta-thorne (last 15 rows)
    if (m < BETA_THORNE_ROWS_COUNT) {
        throw std::runtime_error("Error: Not enough rows for beta-thorne calibration.");
    }
    std::vector<double> median_betathrone(n, 0.0);
    for (unsigned j = 0; j < n; ++j) {
        double sum = 0.0;
        for (unsigned i = m - BETA_THORNE_ROWS_COUNT; i < m; ++i) {
            sum += processed_data[i][j].value;
            processed_data[i][j].is_calibrated = true;
        }
        median_betathrone[j] = sum / BETA_THORNE_ROWS_COUNT;
    }
    double overall_median = std::accumulate(median_betathrone.begin(), median_betathrone.end(), 0.0) / n;
    
    for (unsigned i = 0; i < m; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            if (!processed_data[i][j].is_calibrated) {
                if (median_betathrone[j] != 0) {
                    double proportion = overall_median / median_betathrone[j];
                    processed_data[i][j].value *= proportion;
                }
                else{
                    processed_data[i][j].value = 0.0;
                }
            }
        }
    }

    // Calibration by detectors (last 50 columns)
    if (n < MEDIAN_DETECTOR_COUNT) {
         throw std::runtime_error("Error: Not enough columns for detector calibration.");
    }
    std::vector<double> median_detector(m, 0.0);
    for (unsigned i = 0; i < m; ++i) {
        double sum = 0.0;
        for (unsigned j = n - MEDIAN_DETECTOR_COUNT; j < n; ++j) {
            if (!processed_data[i][j].is_calibrated) {
                sum += processed_data[i][j].value;
                processed_data[i][j].is_calibrated = true;
            }
        }
        median_detector[i] = sum / MEDIAN_DETECTOR_COUNT;
    }

    for (unsigned i = 0; i < m; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            if (!processed_data[i][j].is_calibrated) {
                processed_data[i][j].value /= median_detector[i];
            }
            if (processed_data[i][j].value > 1.0) {
                processed_data[i][j].value = 1.0;
            }
        }
    }
    return processed_data;
}

void create_and_save_image(const std::vector<std::vector<PixelData>>& data, const std::string& filename) {
    unsigned m = data.size();
    unsigned n = data[0].size();
    std::vector<unsigned char> image(m * n * BYTES_PER_PIXEL);
    
    for (unsigned i = 0; i < m; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            int pixel_index = (i * n + j) * BYTES_PER_PIXEL;
            if (data[i][j].is_calibrated) {
                image[pixel_index + 2] = 255; // Red
                image[pixel_index + 1] = 0;
                image[pixel_index + 0] = 0;
            } else {
                unsigned char color_value = static_cast<unsigned char>(data[i][j].value * 255);
                image[pixel_index + 2] = color_value; // Red
                image[pixel_index + 1] = color_value; // Green
                image[pixel_index + 0] = color_value; // Blue
            }
        }
    }
    generateBitmapImage(image.data(), m, n, filename.c_str());
}

// Function to calculate and save thickness image
void calculate_and_save_thickness(const std::vector<std::vector<PixelData>>& data, const std::string& filename) {
    unsigned m = data.size();
    unsigned n = data[0].size();
    std::vector<unsigned char> image(m * n * BYTES_PER_PIXEL);
    double min_thickness = std::numeric_limits<double>::max();

    std::vector<std::vector<double>> thickness_data(m, std::vector<double>(n));

    for (unsigned i = 0; i < m; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            int pixel_index = (i * n + j) * BYTES_PER_PIXEL;
            double v = data[i][j].value;
            double t = (v > 0.0) ? -std::log(v) : 10.0;
            int iv = static_cast<int>(std::round(t * 25.0));
            if (iv < 0) iv = 0;
            if (iv > 255) iv = 255;
            unsigned char color_value = static_cast<unsigned char>(iv);
            image[pixel_index + 2] = color_value;
            image[pixel_index + 1] = color_value;
            image[pixel_index + 0] = color_value;
        }
    }
    generateBitmapImage(image.data(), m, n, filename.c_str());
}

int main() {
    try {
        unsigned m, n;
        auto data = read_data_from_file("block.int", m, n);
        auto processed_data = process_data(data);
        create_and_save_image(processed_data, "normalized_image.bmp");
        
        std::cout << "Image 'normalized_image.bmp' generated successfully." << std::endl;

        int choice;
        std::cout << "Input 1 to check thickness: ";
        std::cin >> choice;
        if (choice == 1) {
            calculate_and_save_thickness(processed_data, "thickness_image.bmp");
            std::cout << "Image 'thickness_image.bmp' generated successfully." << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}




