#include "stdafx.h"
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

void loadData(vector<int> &container, int size) {
  int number;
  while (container.size() < size_t(size)) {
    cin >> number;
    container.push_back(number);
  }
}

void print(int data) {
  cout << setprecision(0) << fixed << floor(data) << ' ';
}

void print(double data) {
  print(int(data));
}

void print(float data) {
  print(int(data));
}


void pokazTablice(vector<int> data) {
  for (size_t i = 0; i < data.size(); i++) {
    print(data[i]);
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

void f0(vector<int> &data, vector<int> &najmniejsze, int from=NULL, int to=NULL) {
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

  vector<int> najmniejsze_indeksy;

  for (int i = from, i < to; i++) {
    if (data[i] == min) {
      najmniejsze_indeksy.push_back(data[i]);
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

  vector<int> _najmniejsza;
  f0(data, _najmniejsza, from, int(data.size()));
	int najmniejsza = int(f0(data, from, to));
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
double f2(vector <int> &data) {
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
double f3(vector<int> &data) {
  double srednia = 0;
  for (auto el : data) {
    srednia += el;
  }
  srednia = srednia / data.size();
  pokazTablice(data);
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

struct Point {
  int x;
  int y;
};

// Znajdź pole wielokąta wypukłego. Zapisz punkty jako tablicę struktur.
// http://en.wikipedia.org/wiki/Shoelace_formula
// wejscie
// 0 0 0 2 2 2 2 0
// wyjscie
// 4

void orderPoints(vector<Point> &points) {
  int sum_x = 0;
  int sum_y = 0;

  for (auto point : points) {
    sum_x += point.x;
    sum_y += point.y;
  }

  double middle_x = sum_x / points.size();
  double middle_y = sum_y / points.size();

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

void cordsToPoints(vector<int> &data, vector<Point> &points) {
  //auto y = data.back();
  //data.pop_back();
  //auto x = data.back();
  //data.pop_back();
  pokazTablice(data);
  for (size_t i = 0; i, i+1 < data.size(); i += 2) {
    Point p;
    p.x = data[i];
    p.y = data[i+1];
    points.push_back(p);
  }
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
    j = (i + 1) % points.size();
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
void f7(vector<int> &data, vector<complex<double>> &wyniki) {
	double zzz = data[0], p = data[1], q = data[2], r = data[3];
	double a = 1 / 3. * (3 * q - p * p);
	double b = 1 / 27. * (2 * pow(p, 3) - 9 * p*q + 27 * r);
	
	double A = pow(-b / 2. + sqrt(b*b / 4. + pow(a, 3) / 24.), 1 / 3.);
	double B = pow(-b / 2. - sqrt(b*b / 4. + pow(a, 3) / 24.), 1 / 3.);

	double y1 = A + B;
	
	complex<double> y1(A + B, 0);
	complex<double> y2(-1 / 2.*(A + B), sqrt(3) / 2.*(A - B));
	complex<double> y3(-1 / 2.*(A + B), -sqrt(3) / 2.*(A - B));


	wyniki.push_back(y1);
	wyniki.push_back(y1);
	wyniki.push_back(y1);

}

// Wyznacz wartość wyrażenia 1 * 2 ^ 2 + 2 * 3 ^ 2 + ... + n(n + 1) ^ 2 dla zadanego n.
// wejscie
// 2
// wyjscie
// 22
double f8(int n) {
  // 1 * 2 ^ 2 + 2 * 3 ^ 2 + ... + n(n + 1) ^ 2
  // n(n + 1)(n + 1)
  // i * i+n * i+n
  // i^3 + 2*i^2 + i
  double sum_natural = n * (n + 1) / 2.;
  double sum_cubes = pow(sum_natural, 2);
  double sum_squeres = sum_natural * (2 * n + 1) / 3.0;

  return sum_cubes + 2 * sum_squeres + sum_natural;
}

// Zlicz liczbę ustawionych bitów w liczbie naturalnej (unsigned).
// wejscie
// 5
// wyjscie
// 2

int f9(unsigned long long int n) {
	int counter = 0;
	while (n != 0) {
		counter += n & 1;
		n >>= 1;
	}

	return counter;
}

void test();

int main() {
	int subprogram, n;
	vector<int> data;
	vector<int> wyniki;
	vector<complex<double>> wyniki_complex;
	vector<Point> points;
  	while (cin >> subprogram >> n) {
		data.clear();
		wyniki.clear();
		wyniki_complex.clear();
    if (n != 9) {
      loadData(data, n);
    }
    switch (subprogram) {
    case 0:
      f0(data, wyniki, 0, data.size());
      pokazTablice(wyniki);
      break;
    case 1:
			f1(data, 0);
			f4(data);
			pokazTablice(data);
			break;
    case 2:
      print(f2(data));
      break;
    case 3:
      print(f3(data));
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
      cordsToPoints(data, points);
      print(f6(points));
      break;
    case 7:
      f7(data, wyniki);
      pokazTablice(wyniki);
      break;
    case 8:
      print(f8(data[0]));
      break;
    case 9:
      vector<unsigned long long> wejscie_9;
      for (int i = 0; i < n; i++) {
        unsigned long long x;
        cin >> x;
        wejscie_9.push_back(x);
      }
      print(f9(wejscie_9);
      break;
    }
  }
  return 0;

}

void test0() {
  vector<int> data = { 8,3,2,5 };
  vector<int> data_p;
  vector<int> data_out = { 3 };
  f0(data, data_p, 0, data.size()) + 1 << endl;
  pokazTablice(data_p)
  pokazTablice(data_out);
}

void test1() {
  vector<int> data1 = { 4, 2, 2, 5, 1, 12 };
  vector<int> data1_out = { 12, 5, 4, 2, 2, 1 };
  f1(data1, 0, data1.size());
  pokazTablice(data1);
  pokazTablice(data1_out);
}
void test2() {
  vector<int> data2 = { 3, 4 };
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
  vector<int> data6 = { 0, 0, 0, 2, 2, 2, 2, 0 };
  vector<Point> data6_points;
  vector<int> data6_out = { 4 };
  cordsToPoints(data6, data6_points);
  cout << f6(data6_points) << endl;
  pokazTablice(data6_out);
}
void test7() {
  vector<int> data7 = { -8,12,-6,1 };
  vector<int> data7_p;
  vector<int> data7_out = { 2 };
  f7(data7, data7_p);
  pokazTablice(data7_p);
  pokazTablice(data7_out);
}
void test8() {
  int data8 = 307;
  vector<int> data8_out = { 307 };
  cout << f8(data8) << endl;
  pokazTablice(data8_out);

}
void test9() {
  int data9 = 31123123;
  vector<int> data9_out = { 16 };
  cout << f9(data9) << endl;
  pokazTablice(data9_out);
}

void test() {
  test0(); count << endl;
  test1(); count << endl;
  test2(); count << endl;
  test3(); count << endl;
  test4(); count << endl;
  test5(); count << endl;
  test6(); count << endl;
  test7(); count << endl;
  test8(); count << endl;
  test9(); count << endl;
}
