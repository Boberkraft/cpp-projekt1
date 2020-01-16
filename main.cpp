//#include "stdafx.h"
#include <iostream>
#include <bitset>
#include <climits>
#include <stdio.h>  
#include <intrin.h>
#include <vector>
#include <iomanip>
#include <assert.h>
#include <math.h>
#include <complex>

using namespace std;

struct Point {
	double x;
	double y;
};

void loadData(vector<int> &container, int size) {
	int number;
	while (container.size() < size_t(size)) {
		cin >> number;
		container.push_back(number);
	}
}

void print(string data, bool sep = true) {
	cout << data;
	if (sep) {
		cout << " ";
	}
}

void print(unsigned long long int data, bool sep = true) {
	cout << setprecision(20) << fixed << data;
	if (sep) {
		cout << " ";
	}
}

void print(int data, bool sep = true) {
	cout << setprecision(0) << fixed << floor(data);
	if (sep) {
		cout << " ";
	}
}


void print(double data, bool sep = true) {
	print(int(floor(data)), sep);
}

void print(float data, bool sep = true) {
	print(int(floor(data)), sep);
}

void inkrementujTablice(vector<int> &data) {
	for (auto &num : data) {
		num += 1;
	}
}

void pokazTablice(vector<string> data) {
	for (size_t i = 0; i < data.size(); i++) {
		print(data[i]);
	}
	cout << endl;
}

void pokazTablice(vector<int> data) {
	for (size_t i = 0; i < data.size(); i++) {
		print(data[i]);
	}
	cout << endl;
}

void pokazTablice(vector<unsigned long long int> data) {
	for (size_t i = 0; i < data.size(); i++) {
		print(data[i]);
	}
	cout << endl;
}

void pokazTablice(vector<complex<double>> data) {
	for (size_t i = 0; i < data.size(); i++) {
		int real = int(floor(data[i].real()));
		int imag = int(floor(data[i].imag()));
		if (real != 0) {
			print(real, false);
		} else if (real == 0 and imag == 0) {
			print("0");
			continue;
		}

		if (imag != 0) {
			if (imag == 1) {
				if (real != 0) {
					print("+", false);
				}
				print("i", false);
			}
			else if (imag == -1) {
				print("-i", false);
			}
			else {
				print(imag, false);
				print(imag);
			}
		}

		print(" ", false);

	}
	cout << endl;
}

bool isDiv(int val, int mod) {
	return val % mod == 0;
}

// So, a more efficient method is to test if n is divisible by 2 or 3, 
// then to check through all the numbers of the form 6k +- 1 <= sqrt(n)
bool isPrime(int val) {
	if (val <= 3) {
		return true;
	}
	if (isDiv(val, 2)) {
		return false;
	}
	if (isDiv(val, 3)) {
		return false;
	}
	int acc = 5;
	while (acc <= sqrt(val)) {
		if (isDiv(val, acc) or isDiv(val, acc + 2)) {
			return false;
		}
	}
	return true;
}


// Znajdź pozycje zawierające najmniejszą wartość dla wejściowego ciągu liczb.
// wejscie
// 8 3 2 5
// wyjscie
// 3

void f0(vector<int> &data, vector<int> &najmniejsze, int from = NULL, int to = NULL) {
	if (from == NULL) { from = 0; }
	if (to == NULL) { to = int(data.size()); }
	int min_i = from;
	int min = data[from];

	for (int i = from; i < to; i++) {
		if (data[i] < min) {
			min = data[i];
			min_i = i;
		}
	}

	for (int i = from; i < to; i++) {
		if (data[i] == min) {
			najmniejsze.push_back(i);
		}
	}
}

// Posortuj wejściowy ciąg liczb. Musisz zastosować funkcję z podproblemu 0).
// wejscie
// 4  2 2 5 1 12
// wyjscie
// 12 5 4 2 2  1
void f1(vector<int> &data, size_t from = 0, size_t to = NULL) {
	if (to == NULL) { to = data.size(); }

	vector<int> _najmniejsza = {};
	f0(data, _najmniejsza, from, int(to));
	int najmniejsza = _najmniejsza[0];
	swap(data[najmniejsza], data[to - 1]);

	if (from < to - 1 and to > 0) {
		f1(data, from, to - 1);
	}
}


// Długość wektora
// wejscie
// 3 4
// wyjscie
// 5
double f2(vector <double> &data) {
	double pod_pierwiastkiem = 0;
	for (auto el : data) {
		pod_pierwiastkiem += pow(el, 2);
	}
	return sqrt(pod_pierwiastkiem);
}

// Znajdź odchylenie standardowe* dla wejściowego ciągu liczb. Musisz zastosować funkcję z podproblemu 2)
// wejscie
// 1 2 3
// wyjscie
// 1
double f3(vector<int> &old_data) {
	vector<double> data(old_data.begin(), old_data.end());

	double srednia = 0;
	for (auto el : data) {
		srednia += el;
	}
	srednia = double(srednia / data.size());

	for (size_t i = 0; i < data.size(); i++) {
		data[i] -= srednia;
	}
	double dlugosc = pow(f2(data), 2);

	// https://wikimedia.org/api/rest_v1/media/math/render/svg/41416d8501c3c4fd63240c29e7fe98fae70cac0e
	return sqrt(dlugosc / double(data.size()));
}

// Zapisz podany wejściowy ciąg liczbowy w tablicy i odwróć go w miejscu (używając tylko wspominanej tablicy).
// wejscie
// 6 4 2 2 5 1 12
// wyjscie
// 12 1 5 2 2 4 6
void f4(vector<int> &data) {
	size_t size = data.size();
	for (size_t i = 0; i < size / 2; i++) {
		swap(data[i], data[size - i - 1]);
	}
}

void orderPoints(vector<Point> &points) {
	int sum_x = 0;
	int sum_y = 0;

	for (auto point : points) {
		sum_x += point.x;
		sum_y += point.y;
	}

	double middle_x = double(sum_x / points.size());
	double middle_y = double(sum_y / points.size());

	vector<int> oryginal_values;
	vector<int> sorted_values;
	for (auto point : points) {
		int angle = int(atan2(point.x - middle_x,
			point.y - middle_y) * 1000);
		oryginal_values.push_back(angle);
		// razy 1000 bo chce miec duzo miejsc po przecinku w incie do porównania
		sorted_values.push_back(angle);
	}

	f1(sorted_values);

	size_t i = 0;
	for (; i < points.size(); i++) {
		size_t j = 0;
		for (; sorted_values[j] != oryginal_values[i]; j++) {}
		swap(points[i], points[j]);
	}

}

void cordsToPoints(vector<double> &data, vector<Point> &points) {
	//auto y = data.back();
	//data.pop_back();
	//auto x = data.back();
	//data.pop_back();
	for (size_t i = 0; i, i + 1 < data.size(); i += 2) {
		Point p;
		p.x = data[i];
		p.y = data[i + 1];
		points.push_back(p);
	}
}

void pokazTabliceRealZComplex(vector<complex<double>> &wejscie) {
	
	vector<Point> punkty;

	vector<double> punkty_x_y;
	for (auto wej : wejscie) {
		punkty_x_y.push_back(wej.real());
		punkty_x_y.push_back(wej.imag());

	}
	cordsToPoints(punkty_x_y, punkty);
	orderPoints(punkty);

	vector<complex<double>> posortowane_wejscie;
	for (auto pkt : punkty) {
		complex<double> zzz(double(pkt.x), double(pkt.y));
		posortowane_wejscie.push_back(zzz);
	}
	pokazTablice(posortowane_wejscie);
}
// Odpowiedz, czy liczba jest pierwsza.
// https://en.wikipedia.org/wiki/Primality_test
// wejscie
// 7 4 5
// wyjscie
// 1 0 1
void f5(vector<int> &data, vector<int> &wyniki) {
	for (size_t i = 0; i < data.size(); i++) {
		wyniki.push_back(isPrime(data[i]));
	}
	//	for (auto el : data) {
	//		el = isPrime(el);
	//	}
}



// Znajdź pole wielokąta wypukłego. Zapisz punkty jako tablicę struktur.
// wejscie
// 0 0 0 2 2 2 2 0
// wyjscie
// 4
double f6(vector<Point> &points) {
	orderPoints(points);
	double area = 0;
	int j;
	for (size_t i = 0; i < points.size(); i++) {
		j = int((i + 1) % points.size());
		area += points[i].x * points[j].y;
		area -= points[j].x * points[i].y;
	}
	area = abs(area) / 2.0;
	return area;
}


// Rozwiąż równanie kwadratowe oraz sześcienne (równanie kwadratowe jest za 50%).
// wejscie
// 0 1 2 1
// wyjscie
// -1
// http://web.cs.iastate.edu/~cs577/handouts/polyroots.pdf
void f7(vector<double> &data, vector<complex<double>> &wyniki) {


	double zzz = data[0], p = data[1], q = data[2], r = data[3];
	if (zzz == 0) {
		double a = p;
		double b = q;
		double c = r;
		double delta = b * b - 4 * a * c;


		complex<double> y1(real(-q / (2 * p) - sqrt(delta) / (2 * p)), imag(-q / (2 * p) + sqrt(delta) / (2 * p)));
		wyniki.push_back(y1);

		if (delta > 0.000001 or delta < -0.000001) {
			complex<double> y2(real(-q / (2 * p) + sqrt(delta) / (2 * p)), imag(-q / (2 * p) - sqrt(delta) / (2 * p)));
			wyniki.push_back(y2);
		}
	}
	else {
		p /= zzz;
		q /= zzz;
		r /= zzz;

		double a = (1 / 3.) * (3 * q - pow(p, 2));
		double b = (1 / 27.) * (2 * pow(p, 3) - 9 * p*q + 27 * r);


		double A = cbrt((-b / 2.) + sqrt(((b*b) / 4.) + pow(a, 3) / 24.));
		double B = cbrt((-b / 2.) - sqrt(((b*b) / 4.) + pow(a, 3) / 24.));

		complex<double> y1((A + B) - (p / 3.), 0);
		complex<double> y2(-1 / 2.*(A + B) - (p / 3.), sqrt(3) / 2.*(A - B));
		complex<double> y3(-1 / 2.*(A + B) - (p / 3.), -sqrt(3) / 2.*(A - B));

		wyniki.push_back(y1);
		wyniki.push_back(y2);
		wyniki.push_back(y3);
	}
}

// Wyznacz wartość wyrażenia 1 * 2 ^ 2 + 2 * 3 ^ 2 + ... + n(n + 1) ^ 2 dla zadanego n.
// wejscie
// 2
// wyjscie
// 22
unsigned long long int f8(unsigned long long int n) {
	// 1 * 2 ^ 2 + 2 * 3 ^ 2 + ... + n(n + 1) ^ 2
	// n(n + 1)(n + 1)
	// i * i+n * i+n
	// i^3 + 2*i^2 + i
	return ((unsigned long long) (n * (n + 1) * (n + 2) * (3 * n + 5))) / 12;
}

// Zlicz liczbę ustawionych bitów w liczbie naturalnej (unsigned).
// wejscie
// 5
// wyjscie
// 2

int p(char a) {
	return int(a - '0');
}
string dziel(string a) {
	int val = 0;
	string bb = "";
	for (int i = 0; i < a.size(); i++) {
		int dup = 0;
		val = val * 10 + p(a[i]);
		if ((val & 1) == 1) {
			dup = 1;
		}
		val = val / 2;
		bb = bb + to_string(val);
		val = dup;
	}
	if (bb[0] == '0') {
		return bb.substr(1);
	}

	return bb;
}

int prznies(string a, long b) {
	if (a.size() == 0) {
		return b;
	}
	if (a.size() == 1 and a[0] == '1') {
		return ++b;
	}
	int n = p(a[a.size() - 1]);
	if ((n & 1) == 1) {
		++b;
	}
	a = dziel(a);
	return prznies(a, b);
}

void f9(vector<string> data, vector<int> &wynik) {
	for (auto napis : data) {
		wynik.push_back(prznies(napis, 0));
	}
}


void test();

int main() {
	//test();
	//system("PAUSE");

	int subprogram, n;
	vector<int> data;
	vector<int> wyniki;
	vector<double> dane_dobule;
	vector<string> dane_string;
	vector<complex<double>> wyniki_complex;
	vector<Point> points;
	while (cin >> subprogram >> n) {
		data.clear();
		wyniki.clear();
		dane_dobule.clear();
		wyniki_complex.clear();
		dane_string.clear();
		points.clear();
		if (subprogram != 9 and subprogram != 7 and subprogram != 6) {
			loadData(data, n);
		}
		switch (subprogram) {
		case 0:
			f0(data, wyniki);
			inkrementujTablice(wyniki);
			pokazTablice(wyniki);
			break;
		case 1:
			f1(data, 0);
			pokazTablice(data);
			break;
		case 2:
		{
			vector<double> double_data(data.begin(), data.end());
			print(f2(double_data));
			cout << endl;
			break;
		}
		case 3:
			print(f3(data));
			cout << endl;
			break;
		case 4:
			f4(data);
			pokazTablice(data);
			break;
		case 5:
			f5(data, wyniki);
			pokazTablice(wyniki);
			break;
		case 6:
			points.clear();
			for (int i = 0; i < 4; i++) {
				double x;
				cin >> x;
				dane_dobule.push_back(x);
			}
			cordsToPoints(dane_dobule, points);
			print(f6(points));
			cout << endl;
			break;
		case 7:
			for (int i = 0; i < 4; i++) {
				double x;
				cin >> x;
				dane_dobule.push_back(x);
			}
			f7(dane_dobule, wyniki_complex);
			if (dane_dobule[0] == 1) {
				pokazTabliceRealZComplex(wyniki_complex);
			}
			else {
				pokazTablice(wyniki_complex);
			}

			break;
		case 8:
			print(f8(data[0]));
			cout << endl;
			break;
		case 9:
			for (int i = 0; i < n; i++) {
				string wejscie_9;
				cin >> wejscie_9;
				dane_string.push_back(wejscie_9);
			}
			f9(dane_string, wyniki);
			pokazTablice(wyniki);
			break;
		}
	}
	return 0;
}

void test0() {
	vector<int> data = { 8,3,2,5 };
	vector<int> data_p;
	vector<int> data_out = { 3 };
	f0(data, data_p);
	inkrementujTablice(data_p);
	pokazTablice(data_p);
	pokazTablice(data_out);
}
void test1() {
	vector<int> data1 = { 20, 14, 6, 4, 10, 5, 7, 7, 12, 5, 14, 3, 9, 13, 15, 13, 2, 16, 6, 5, 19, 3, 9 };
	vector<int> data1_out = { 20, 19, 16, 15, 14, 14, 13, 13, 12, 10, 9, 9, 7, 7, 6, 6, 5, 5, 5, 4, 3, 3, 2 };
	f1(data1);
	pokazTablice(data1);
	pokazTablice(data1_out);
}
void test2() {
	vector<double> data2 = { 3., 4. };
	vector<int> data2_out = { 5 };
	cout << f2(data2) << endl;
	pokazTablice(data2_out);
}
void test3() {
	vector<int> data3 = { 1, 2, 3 };
	vector<int> data3_out = { 1 };
	cout << f3(data3) << endl;
	pokazTablice(data3_out);
}
void test4() {
	vector<int> data4 = { 1, 2, 2, 4, 5, 12 };
	vector<int> data4_out = { 12, 5, 4, 2, 2, 1 };
	f4(data4);
	pokazTablice(data4);
	pokazTablice(data4_out);
}
void test5() {
	vector<int> data5 = { 7, 4, 5 };
	vector<int> data5_p;
	vector<int> data5_out = { 1, 0, 1 };
	f5(data5, data5_p);
	pokazTablice(data5);
	pokazTablice(data5_p);
	pokazTablice(data5_out);
}
void test6() {
	vector<double> data6 = { 91, 84, -75, 76, -88, 67, -33, -49, 54, 17 };
	vector<Point> data6_points;
	vector<int> data6_out = { 13238 };
	cordsToPoints(data6, data6_points);
	cout << f6(data6_points) << endl;
	pokazTablice(data6_out);
}
void test7() {
	vector<double> data7 = { 1, 1.08243, -0.785653 , 1.52854 };
	vector<complex<double>> data7_p;
	vector<string> data7_out = { "-2", "-i", "0"};
	f7(data7, data7_p);
	pokazTabliceRealZComplex(data7_p);
	pokazTablice(data7_out);
}
void test8() {
	int data8 = 938;
	vector<unsigned long long int> data8_out = { 194495749210 };
	cout << f8(data8) << endl;
	pokazTablice(data8_out);

}
void test9() {
	vector<string> data9 = { "31123123" };
	vector<int> data5_p;
	vector<int> data9_out = { 16 };
	f9(data9, data5_p);
	pokazTablice(data5_p);
	pokazTablice(data9_out);
}

void test() {
	cout << " --- Test 0" << endl; test0(); cout << endl;
	cout << " --- Test 1" << endl; test1(); cout << endl;
	cout << " --- Test 2" << endl; test2(); cout << endl;
	cout << " --- Test 3" << endl; test3(); cout << endl;
	cout << " --- Test 4" << endl; test4(); cout << endl;
	cout << " --- Test 5" << endl; test5(); cout << endl;
	cout << " --- Test 6" << endl; test6(); cout << endl;
	cout << " --- Test 7" << endl; test7(); cout << endl;
	cout << " --- Test 8" << endl; test8(); cout << endl;
	cout << " --- Test 9" << endl; test9(); cout << endl;
}
