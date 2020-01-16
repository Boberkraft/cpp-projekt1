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
#include <cstdlib>
#include <numeric>

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

void print(long double data, bool sep = true) {
	print(int(floor(data)), sep);
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
		}
		else if (real == 0 and imag == 0) {
			print("0");
			continue;
		}

		if (imag > 0) {
			if (real != 0) {
				print("+", false);
			}
		}
		if (imag < 0) {
			print("-", false);
		}
		if (imag != 0) {
			if (abs(imag) > 1) {
				print(abs(imag), false);
			}
			print("i", false);
		}

		print(" ", false);

	}
	cout << endl;
}

inline bool isDiv(unsigned int val, unsigned int mod) {
	return val % mod == 0;
}





// So, a more efficient method is to test if n is divisible by 2 or 3, 
// then to check through all the numbers of the form 6k +- 1 <= sqrt(n)
inline bool isPrime(unsigned val) {
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
	while (acc * acc <=val) {
		if (isDiv(val, acc) or isDiv(val, acc + 2)) {
			return false;
		}
		acc += 6;
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
long double f2(vector <double> &data) {
	long double pod_pierwiastkiem = 0;
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
long double f3(vector<int> &old_data) {
	vector<double> data(old_data.begin(), old_data.end());

	double srednia = 0;
	for (auto el : data) {
		srednia += el;
	}
	srednia = double(srednia / data.size());
	for (size_t i = 0; i < data.size(); i++) {
		data[i] -= srednia;
	}

	long double dlugosc = pow(f2(data), 2);

	
	// https://wikimedia.org/api/rest_v1/media/math/render/svg/41416d8501c3c4fd63240c29e7fe98fae70cac0e
	return sqrt((dlugosc + 0.04) / ((long double) data.size()));
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

void orderPoints(vector<Point> &points, vector<Point> &out_points) {
	int sum_x = 0;
	int sum_y = 0;

	for (auto point : points) {
		sum_x += point.x;
		sum_y += point.y;
	}

	double middle_x = double(sum_x / (double)points.size());
	double middle_y = double(sum_y / (double)points.size());

	vector<int> oryginal_values;
	vector<int> sorted_values;
	for (auto point : points) {
		int angle = int(atan2(point.y - middle_y,
			point.x - middle_x) * 1000);
		//cout << setprecision(3) << atan2(point.y - middle_y,
		//	point.x - middle_x);
		//cout << endl;
		oryginal_values.push_back(angle);
		// razy 1000 bo chce miec duzo miejsc po przecinku w incie do porównania
		sorted_values.push_back(angle);
	}

	f1(sorted_values);


	for (size_t i = 0; i < points.size(); i++) {
		size_t j = 0;
		for (;  sorted_values[i] != oryginal_values[j]; j++) {}
		out_points.push_back(points[j]);
	}

}

void cordsToPoints(vector<double> &data, vector<Point> &points) {
	//auto y = data.back();
	//data.pop_back();
	//auto x = data.back();
	//data.pop_back();
	for (size_t i = 0; i, i < data.size(); i += 2) {
		Point p;
		p.x = data[i];
		p.y = data[i + 1];
		points.push_back(p);
	}
}

void pokazTabliceRealZComplex(vector<complex<double>> &wejscie,  int &start_logging) {
	

	vector<int> oryginalne_real;
	vector<int> oryginalne_imag;
	vector<int> posortowane_real;
	vector<int> posortowane_imag;

	vector<complex<double>> wejscie_kopia(wejscie.begin(), wejscie.end());
	vector<complex<double>> wejscie_kopia1(wejscie.begin(), wejscie.end());
	vector<complex<double>> posortowane_real_complex;
	vector<complex<double>> posortowane_imag_complex;

	vector<complex<double>> posortowane_wejscie;

	for (auto pkt : wejscie) {
		posortowane_real.push_back(-int(floor(pkt.real())));
		posortowane_imag.push_back(-int(floor(pkt.imag())));
	}
	
	f1(posortowane_real);
	f1(posortowane_imag);

	for (auto pkt : posortowane_real) {
		for (auto &org_pkt : wejscie_kopia) {
			if (int(floor(org_pkt.real())) == -pkt) {
				posortowane_real_complex.push_back(org_pkt);
				complex<double> y(-99999, -999999);
				org_pkt = y;
				break;
			}
		}
	}
	for (auto pkt : posortowane_imag) {
		for (auto &org_pkt : wejscie_kopia1) {
			if (int(floor(org_pkt.imag())) == -pkt) {
				
				posortowane_imag_complex.push_back(org_pkt);
				complex<double> y(-99999, -999999);
				org_pkt = y;
				break;
			}
		}
	}

	for (auto real : posortowane_real_complex) {
		// 5-2i, 5-3i
		for (auto &imag : posortowane_imag_complex) {
			if (int(floor(imag.real()))== int(floor(real.real()))) {
				posortowane_wejscie.push_back(imag);
				complex<double> y(-99999, -999999);
				imag = y;
				break;
			}
		}
	}

	if (posortowane_wejscie.size() == 3) {
		if (int(floor(posortowane_wejscie[0].real())) == -1 and
			int(floor(posortowane_wejscie[0].imag())) == 0 and
			int(floor(posortowane_wejscie[1].real())) == 0 and
			int(floor(posortowane_wejscie[1].imag())) == 0 and
			int(floor(posortowane_wejscie[2].real())) == 2 and
			int(floor(posortowane_wejscie[2].imag())) == 0) {
			//start_logging = 10;
		}
	}

	//cout << setprecision(5) << posortowane_wejscie[0] << " ";
	//cout << setprecision(5) << posortowane_wejscie[1] << " ";
	//cout << setprecision(5) << posortowane_wejscie[2] << " " << endl;
	//pokazTablice(posortowane_wejscie);
	pokazTablice(posortowane_wejscie);
}
// Odpowiedz, czy liczba jest pierwsza.
// https://en.wikipedia.org/wiki/Primality_test
// wejscie
// 7 4 5
// wyjscie
// 1 0 1
void f5(vector<int> &data, int &x) {
	for (auto el : data) {

		if (isPrime(el)) {
			cout << "1 ";
		} else {
			cout << "0 ";
		}
	}
	cout << "\n";
}

// Znajdź pole wielokąta wypukłego. Zapisz punkty jako tablicę struktur.
// wejscie
// 0 0 0 2 2 2 2 0
// wyjscie
// 4
double f6(vector<Point> &up) {

	vector<Point> p;
 	orderPoints(up, p);

	double area = 0.0;
	int j = p.size() - 1;
	for (size_t i = 0; i < p.size(); i++) {
		area += (p[j].x + p[i].x) * (p[j].y - p[i].y);
		j = i;
	}
	area = abs(area) / 2.0;
	return area;
}

complex<double> cube_root(complex<double> x) {
	if (x.real() >= 0) {
		return pow(x, 1. / 3.);
	} else {
		return -pow(-x, 1. / 3.);
	}
}

// Rozwiąż równanie kwadratowe oraz sześcienne (równanie kwadratowe jest za 50%).
// wejscie
// 0 1 2 1
// wyjscie
// -1
// http://web.cs.iastate.edu/~cs577/handouts/polyroots.pdf
void f7(vector<double> &data, vector<complex<double>> &wyniki, int start_logging) {
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
		if (start_logging == 10) {
			cout << zzz << " " << p << " " << q << " " << r << " " << endl;
		}
		//
		p /= zzz;
		q /= zzz;
		r /= zzz;

		double a = (1 / 3.) * (3 * q - pow(p, 2));
		double b = (1 / 27.) * (2 * pow(p, 3) - 9 * p*q + 27 * r);

		complex<double> wnetrze_A = ((b*b) / 4.) + pow(a, 3) / 27.;
		complex<double> wnetrze_B = ((b*b) / 4.) + pow(a, 3) / 27.;
		complex<double> A = cube_root((-b / 2.) + sqrt(wnetrze_A));
		complex<double> B = cube_root((-b / 2.) - sqrt(wnetrze_B));

		complex<double> y1 = ((A + B) - (p / 3.));
		complex<double> y2 = (-1. / 2.)*(A + B) + (complex<double>(0, 1) * sqrt(3)) / 2.*(A - B) - (p / 3.);
		complex<double> y3 = (-1. / 2.)*(A + B) - (complex<double>(0, 1) * sqrt(3)) / 2.*(A - B) - (p / 3.);
		//cout << setprecision(5) << (wnetrze_A) << endl;
		//cout << setprecision(5) << " " <<y1 << " " << y2 << " " << y3 << " " << endl;
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
	int how_many = 4;
	int start_logging = 9;
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
			f5(data,  how_many);
			break;
		case 6:
			for (int i = 0; i < n; i++) {
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
			if (n != 4) {
				cout << "XDD";
			}
			f7(dane_dobule, wyniki_complex, start_logging);
			pokazTabliceRealZComplex(wyniki_complex, start_logging);

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
	vector<int> data3 = { 11, 17, 17, 20};
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
	// 
	
	vector<int> data5 = { 1073958442, 1175845667, 2037467347, 355373657 };
	vector<int> data5_out = { 0, 1, 1, 1};
	pokazTablice(data5);
	int x = 2;
	f5(data5, x);
	pokazTablice(data5_out);
}
void test6() {
	vector<double> data6 = { 91, 84, -75, 76, -88, 67, -33, -49, 54, 17 };
	vector<Point> data6_points;
	vector<int> data6_out = { 13238 };
	cordsToPoints(data6, data6_points);
	cout << f6(data6_points) << endl;
	pokazTablice(data6_out);


	vector<double> adata6 = { 94, -21, 94, 50, 81, 56, 26, 70, -26, 54, -88, - 35, -79, -80, 4, -71 };
	vector<Point> adata6_points;
	vector<int> adata6_out = { 18891 };
	cordsToPoints(adata6, adata6_points);
	cout << f6(adata6_points) << endl;
	pokazTablice(adata6_out);
}
void test7() {
	int start_logging = 9;

	  

	vector<double> data7 = { -1686, -1056,  -5363, 8066 };
	vector<complex<double>> data7_p;
	vector<string> data7_out = { "-1-3i",  "-1+2i", "0" };
	f7(data7, data7_p, start_logging);
	pokazTabliceRealZComplex(data7_p, start_logging);
	pokazTablice(data7_out);
	
	return;

	vector<double> adata7 = { 0, 9, 126, 441 };
	vector<complex<double>> adata7_p;
	vector<int> adata7_out = { -7 };
	f7(adata7, adata7_p, start_logging);
	pokazTabliceRealZComplex(adata7_p, start_logging);
	pokazTablice(adata7_out);
	
	vector<double> bdata7 = { -3342, -5139, 2139, 1690 };
	vector<complex<double>> bdata7_p;
	vector<string> bdata7_out = { "-2", "-i", "0" };
	f7(bdata7, bdata7_p, start_logging);
	pokazTabliceRealZComplex(bdata7_p, start_logging);
	pokazTablice(bdata7_out);

	vector<double> cdata7 = { -5799, -6277, 4556, -8864 };
	vector<complex<double>> cdata7_p;
	vector<string> cdata7_out = { "-2", "-i", "0" };
	f7(cdata7, cdata7_p, start_logging);
	pokazTabliceRealZComplex(cdata7_p, start_logging);
	pokazTablice(cdata7_out);

	vector<double> ddata7 = { -5393 , 7073, -5087, 3555 };
	vector<complex<double>> ddata7_p;
	vector<string> ddata7_out = { "-i", "0", "1" };
	f7(ddata7, ddata7_p, start_logging);
	pokazTabliceRealZComplex(ddata7_p, start_logging);
	pokazTablice(ddata7_out);
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
