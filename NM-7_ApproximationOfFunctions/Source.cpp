#include <iostream>
#include <vector>
#include <gnuplot-iostream.h>

#define ANSI_COLOR_BLUE "\033[34m"
#define ANSI_COLOR_RESET "\033[0m"
#define ANSI_COLOR_GREEN "\033[32m"
#define ANSI_COLOR_RED "\033[31m"
#define ANSI_COLOR_PURPLE "\033[35m"
#define ANSI_COLOR_LIGHTBLUE "\033[94m"

using namespace std;

/**
 * @brief Функція знаходження мінору матриці
 * @param matrix Матриця коефіцієнтів Х
 * @param row Рядок для видалення
 * @param column Стовпець для видалення
 * @return Нової матриці після видалення
 */
vector<vector<double>> GetMinor(const vector<vector<double>>& matrix, const int& row, const int& column) {

	vector<vector<double>> minor;
	for (int i = 0; i < matrix.size(); i++) {
		if (i == row) continue;
		vector<double> rowMinor;
		for (int j = 0; j < matrix[i].size(); j++) {
			if (j == column) continue;
			rowMinor.push_back(matrix[i][j]);
		}
		minor.push_back(rowMinor);
	}
	return minor;
}

/**
 * @brief Функція знаходження визначника матриці
 * @param matrix Матриця
 * @return Визначник матриці
 */
double GetDeterminate(const vector<vector<double>>& matrix) {

	if (matrix.size() == 1) return matrix[0][0];
	if (matrix.size() == 2) {
		double a = matrix[0][0]; double b = matrix[0][1];
		double c = matrix[1][0]; double d = matrix[1][1];
		return a * d - c * b;
	}
	double det = 0;
	for (int j = 0; j < matrix.size(); j++) {
		double coef = matrix[0][j];
		double sign;
		if (j % 2 == 0) sign = 1;
		else sign = -1;
		det += sign * coef * GetDeterminate(GetMinor(matrix, 0, j));
	}
	return det;
}

/**
 * @brief Функція заміни стовпця на вектор вільних членів
 * @param matrix Матриця коефіцієнтів Х
 * @param B Вектор вільних членів
 * @param column Стовпець заміни
 * @return Матриця з заміним стовпцем на вільні члени
 */
vector<vector<double>> ReplaceMatrixColumn(const vector<vector<double>>& matrix, const vector<double>& B, const int& column) {

	vector<vector<double>> newMatrix(matrix.size(), vector<double>(matrix[0].size(), 0));
	int index = 0;
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix[i].size(); j++) {
			if (j == column) { newMatrix[i][j] = B[index]; index++; }
			else newMatrix[i][j] = matrix[i][j];
		}
	}
	return newMatrix;
}

/**
 * @brief Функція, яка реалізує метод Крамера для вирішення СЛАР
 * @param matrix Матриця А
 * @param B Вектор вільних членів
 */
vector<double> CramerMethod(const vector<vector<double>>& matrix, const vector<double>& B) {

	cout << "-----Welcome to Cramer Method-----" << endl;
	if (matrix.empty() || matrix[0].size() != matrix.size() || B.size() != matrix.size()) {
		throw "Error: <Invalid matrix size>";
	}
	double det = 0;
	vector<double> detX;
	vector<double> x;
	det = GetDeterminate(matrix);
	if (abs(det) < 1e-9) throw "Error: <System has no unique solution (det(A) = 0)";
	else cout << "det A = " << ANSI_COLOR_GREEN << det << ANSI_COLOR_RESET << endl;
	for (int i = 0; i < matrix.size(); i++) {
		detX.push_back(GetDeterminate(ReplaceMatrixColumn(matrix, B, i)));
	}
	cout << "detXn : ";
	for (auto i : detX) {
		cout << ANSI_COLOR_GREEN << i << " " << ANSI_COLOR_RESET;
	}cout << endl;
	for (int i = 0; i < matrix.size(); i++) {
		x.push_back(detX[i] / det);
	}
	return x;
}

struct PointsXY {
	double x, y;

	friend ostream& operator<<(ostream& os, const PointsXY& point) {
		os << point.x << " " << point.y;
		return os;
	}
};

/**
 * @brief Функція побудови графіка
 * @param points Система точок
 * @param coefLinear Параметри лінійної апроксимації
 * @param coefQudratic Параметри квадратичної апроксимації
 */
void MakeGraph(const vector<PointsXY>& points, const vector<double>& coefLinear, const vector<double>& coefQudratic) {

	double minX = points.front().x, maxX = points.back().x;
	vector<pair<double, double>> linearApprox;
	vector<pair<double, double>> quadraticApprox;
	for (double x = minX - 1; x <= maxX + 1; x += 0.1) {
		linearApprox.emplace_back(x, coefLinear[0] * x + coefLinear[1]);
		quadraticApprox.emplace_back(x, coefQudratic[0] * x * x + coefQudratic[1] * x + coefQudratic[2]);
	}
	Gnuplot gr;
	gr << "set title 'Linear vs Quadratic Approximations'\n";
	gr << "set xlabel 'X'\n";
	gr << "set ylabel 'Y'\n";
	gr << "set grid\n";
	gr << "plot '-' with points pointtype 7 title 'Points', "
		"'-' with lines title 'y = -2x + 3,98', "
		"'-' with lines title 'y = 0,9714x^2 - 2x + 2,03714'\n";
	gr.send1d(points);
	gr.send1d(linearApprox);
	gr.send1d(quadraticApprox);
}

/**
 * @brief Функція отримання коефіцієнта детермінації R^2
 * @param points Система точок
 * @param coefs Параметри апроксимації
 * @return Коефіцієнт детермінації
 */
double GetR2(const vector<PointsXY>& points, const vector<double>& coefs) {

	double ssRes = 0, ssTot = 0;
	double yMean = 0;

	for (const auto& p : points) {
		yMean += p.y;
	}
	yMean /= points.size();

	for (const auto& p : points) {
		double preddictedY = 0;
		if (coefs.size() == 2) preddictedY = coefs[0] * p.x + coefs[1];
		else if (coefs.size() == 3) preddictedY = coefs[0] * p.x * p.x + coefs[1] * p.x + coefs[2];
		ssRes += pow(p.y - preddictedY, 2);
		ssTot += pow(p.y - yMean, 2);
	}

	return 1 - (ssRes / ssTot);
}

/**
 * @brief Функція реалізації лінійної апроксимації
 * @param points Система точок
 * @return Параметри лінійної апроксимації
 */
vector<double> LinearApproximation(const vector<PointsXY>& points) {

	double sumX = 0, sumY = 0, sumXX = 0, sumXY = 0;
	const int n = points.size();

	for (const auto& p : points) {
		sumX += p.x;
		sumY += p.y;
		sumXX += pow(p.x, 2);
		sumXY += p.x * p.y;
	}

	double a = (n * sumXY - sumX * sumY) / (n * sumXX - pow(sumX, 2));
	double b = (sumY - a * sumX) / n;
	vector<double> coefLinear = { a, b };
	cout << "Linear approximation: y = " << ANSI_COLOR_PURPLE << a << "x + " << b << ANSI_COLOR_RESET << endl;
	cout << "R^2 : " << ANSI_COLOR_LIGHTBLUE << GetR2(points, coefLinear) << ANSI_COLOR_RESET << endl;

	return coefLinear;
}

/**
 * @brief Функція реалізації квадратичної апроксимації
 * @param points Система точок
 * @return Параметри квадратичної апроксимації
 */
vector<double> QuadraticApproximation(const vector<PointsXY>& points) {

	double sumX = 0, sumY = 0, sumXX = 0, sumXY = 0, sumXXX = 0, sumXXXX = 0, sumXXY = 0;
	const int n = points.size();

	for (const auto& p : points) {
		sumX += p.x;
		sumY += p.y;
		sumXX += pow(p.x, 2);
		sumXY += p.x * p.y;
		sumXXX += pow(p.x, 3);
		sumXXXX += pow(p.x, 4);
		sumXXY += pow(p.x, 2) * p.y;
	}

	vector<vector<double>> matrix = {
		{sumXXXX, sumXXX, sumXX},
		{sumXXX, sumXX, sumX},
		{sumXX, sumX, static_cast<double>(n)}
	};
	vector<double> B = { sumXXY, sumXY, sumY };
	vector<double> coefQudratic = CramerMethod(matrix, B);
	cout << "Quadratic approximation: y = " << ANSI_COLOR_PURPLE << coefQudratic[0] << "x^2 " << coefQudratic[1] << "x + " << coefQudratic[coefQudratic.size() - 1] << ANSI_COLOR_RESET << endl;
	cout << "R^2 : " << ANSI_COLOR_LIGHTBLUE << GetR2(points, coefQudratic) << ANSI_COLOR_RESET << endl;

	return coefQudratic;
}

int main() {

	vector<PointsXY> points = {
		{-2, 9.9}, {-1, 5.1}, {0, 1.9}, {1, 1.1}, {2, 1.9}
	};

	try {
		vector<double> coefLinear = LinearApproximation(points);
		vector<double> coefQudratic = QuadraticApproximation(points);
		MakeGraph(points, coefLinear, coefQudratic);
	}
	catch (const char* err) {
		cerr << ANSI_COLOR_RED << err << ANSI_COLOR_RESET << endl;
	}	

	return 0;
}