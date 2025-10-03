#include <chrono>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <algorithm>

std::vector<cv::Point2f> control_points;

void mouse_handler(int event, int x, int y, int flags, void* userdata)
{
	if (event == cv::EVENT_LBUTTONDOWN && control_points.size() < 4)
	{
		std::cout << "Left button of the mouse is clicked - position (" << x << ", "
			<< y << ")" << '\n';
		control_points.emplace_back(x, y);
	}
}

void naive_bezier(const std::vector<cv::Point2f>& points, cv::Mat& window)
{
	auto& p_0 = points[0];
	auto& p_1 = points[1];
	auto& p_2 = points[2];
	auto& p_3 = points[3];

	for (double t = 0.0; t <= 1.0; t += 0.001)
	{
		auto point = std::pow(1 - t, 3) * p_0 + 3 * t * std::pow(1 - t, 2) * p_1 +
			3 * std::pow(t, 2) * (1 - t) * p_2 + std::pow(t, 3) * p_3;

		window.at<cv::Vec3b>(point.y, point.x)[2] = 255;
	}
}

cv::Point2f recursive_bezier(const std::vector<cv::Point2f>& control_points, float t)
{
	// TODO: Implement de Casteljau's algorithm
	// return cv::Point2f();

	if (control_points.empty())
	{
		throw std::invalid_argument("control_points cannot be empty");
	}

	if (control_points.size() == 1)
	{
		return control_points[0];
	};

	std::vector<cv::Point2f> next_control_points;
	for (int i = 0; i < control_points.size() - 1; i++)
	{
		const cv::Point2f& a = control_points[i];
		const cv::Point2f& b = control_points[i + 1];
		cv::Point2f point = (1 - t) * a + t * b;
		next_control_points.emplace_back(point);
	}

	return recursive_bezier(next_control_points, t);
}


void bezier(const std::vector<cv::Point2f>& control_points, cv::Mat& window)
{
	// TODO: Iterate through all t = 0 to t = 1 with small steps, and call de Casteljau's 
	// recursive Bezier algorithm.

	for (double t = 0.0; t <= 1.0; t += 0.001)
	{
		cv::Point2f point = recursive_bezier(control_points, t);
		//window.at<cv::Vec3b>(point.y, point.x)[2] = 255;

		// Anti-aliasing for the point
		int start_x = static_cast<int>(std::floor(point.x));
		int start_y = static_cast<int>(std::floor(point.y));
		int end_x = static_cast<int>(std::ceil(point.x));
		int end_y = static_cast<int>(std::ceil(point.y));
		if (point.x - start_x < 0.5) start_x -= 1;
		if (point.y - start_y < 0.5) start_y -= 1;

		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				int subpixel_x = start_x + i;
				int subpixel_y = start_y + j;
				if (subpixel_x < 0 || subpixel_y < 0 || subpixel_x >= window.cols || subpixel_y >= window.rows) continue;

				float dx = point.x - (subpixel_x + 0.5f);
				float dy = point.y - (subpixel_y + 0.5f);
				float distance = std::sqrt(dx * dx + dy * dy);
				float weight = std::max(0.0f, 1.0f - distance);

				float old_value = window.at<cv::Vec3b>(subpixel_y, subpixel_x)[1];
				// Lerp between original color and target 
				float mixed = std::min(255.0f, (1.0f - weight) * old_value + 255 * weight);

				window.at<cv::Vec3b>(subpixel_y, subpixel_x)[1] = mixed;
			}
		}
	}
}

int main()
{
	cv::Mat window = cv::Mat(700, 700, CV_8UC3, cv::Scalar(0));
	cv::cvtColor(window, window, cv::COLOR_BGR2RGB);
	cv::namedWindow("Bezier Curve", cv::WINDOW_AUTOSIZE);

	cv::setMouseCallback("Bezier Curve", mouse_handler, nullptr);

	int key = -1;

	// test for bezier function
	control_points.push_back({ 100, 300 });
	control_points.push_back({ 200, 200 });
	control_points.push_back({ 300, 200 });
	control_points.push_back({ 400, 300 });
	//control_points.push_back({ 500, 300 });

	while (key != 27)
	{
		for (auto& point : control_points)
		{
			cv::circle(window, point, 3, { 255, 255, 255 }, 3);
		}


		if (control_points.size() == 4)
		{
			naive_bezier(control_points, window);
			bezier(control_points, window);

			cv::imshow("Bezier Curve", window);
			cv::imwrite("my_bezier_curve.png", window);
			key = cv::waitKey(0);

			return 0;
		}

		cv::imshow("Bezier Curve", window);
		key = cv::waitKey(20);
	}

	return 0;
}
