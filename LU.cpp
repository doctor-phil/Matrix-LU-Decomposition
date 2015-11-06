#include <iostream>
#include <vector>
#include <complex>

using namespace std;

typedef complex<double> com;
typedef vector<com> vec;				//define vec type as a vector filled with double type variables
typedef vector<vec> mat;				//define mat type as a matrix aka a vector filled with vectors

void LU(mat&, mat&);					//LU Decomp & forward-backward substitution function prototypes
vec FB(mat&, mat&, vec);

int main() {
	int n;
	char cont;
	do {
		cout << "LU DECOMPOSITION OF nxn MATRIX\n";
		cout << "n = ";									
		cin >> n;
		vec row(n);							//new double vector with n entries
		mat A, I;							//matrixes to hold A & I
		mat L, U;							//these will hold LU after running the LU() function
		vec b(n);							//new vector for b with n entries
		vec x(n);							//this vector will hold the solution to Ax=b after FB()
	
		for (int i = 0; i<n; i++) {
			for (int i2=0; i2<n; i2++)		//initialize nxn identity matrix
				if (i2==i)
					row[i2]=1;
				else
					row[i2]=0;
			I.push_back(row);
		}
	
		cout << "Enter A by row, with each digit followed by the enter key: \n";
		for (int i = 0; i<(n); i++) {
			for (int i2 = 0; i2<n; i2++)	//user enters A
				cin >> row[i2];				//row is a temp variable that holds rows before they're entered in A
			A.push_back(row);				//append A with this row vector
		}
		cout << "Enter b with each digit followed by the enter key: \n";
		for (int i = 0; i<n; i++) {
			cin >> b[i];					//user enters b
		}
	
		cout << "\n\nA =\n";				//print A & b to the screen
		for (int i = 0; i<(n); i++) {
			cout << "\t[\t";
			for (int i2 = 0; i2<n; i2++)
				cout << A[i][i2] << "\t";
			cout << "\t]" << endl;
		}
		cout << "\nb =\n";
		for (int i = 0; i<n; i++)
			cout << "\t[\t" << b[i] << "\t]" << "\n";
		
		L = I;								//copy I into L, and A into U
		U = A;								//function can change these without changing I or L in main()
		LU(U,L);							//decompose A into L and U
	
		cout << endl << endl;
		cout << "RESULT:\nL = " << endl;	//print results of LU decomposition
		for (int i = 0; i<n; i++) {
			cout << "\t[\t";
			for (int i2=0;i2<n;i2++)
				cout << L[i][i2] << "\t";
			cout << "\t]\n";
		}
		cout << endl << endl;
		cout << "U = " << endl;
		for (int i = 0; i<n; i++) {
			cout << "\t[\t";
			for (int i2=0;i2<n;i2++)
				cout << U[i][i2] << "\t";
			cout << "\t]\n";
		}
	
		x = FB(L,U,b);						//perform forward & backward substitution to find x: Ax=b
	
		cout << "\nSOLUTION (Ax=b):\nx =\n";//print x
		for (int i = 0; i<n; i++)
			cout << "\t[\t" << x[i] << "\t]\n";
		cout << endl;
		cout << "Enter Q to quit, or any other character to continue: ";
		cin >> cont;
		cout << endl;
	} while (cont!='Q'&&cont!='q');
	return 0;
}

void LU(mat &U, mat &L) {
	int n = U.size();					//set n based on # of columns in A
										//array indexes in C++ are 0:n-1 rather than 1:n
										//loop bounds were changed to reflect this
	for (int k = 0; k<n-1; k++) {
		for (int k2 = k+1; k2<n; k2++) {
			L[k2][k] = U[k2][k]/U[k][k];//translated from prof. online matlab code for LU decomposition:
		}								//C++ has no : operator, so an addt'l loop is req'd for i2
		for (int j = k+1; j<n; j++) {
			for (int j2=k; j2<n; j2++) {
				U[j][j2] = U[j][j2] - L[j][k]*U[k][j2];
			}
		}								//this function uses pass by address instead of returning
	}									//a new value, allowing main() to get both L and U from LU()
}

vec FB(mat &L, mat &U, vec b) {
	int n = b.size();
	vec x(n), y(n);						//LUx=b, Ux=y, Ly=b
										
										//begin solving for y:Ly=b
	for (int i=0;i<n;i++) {
		y[i] = b[i];					//using derived algorithm yk=(bk-sum[i=1,k-1](lki*yi))/lkk
										//runs "forward" from y1 to yn
		for (int ii=0;ii<i;ii++) {
			y[i] = y[i] - (L[i][ii] * y[ii]);	//y[i]=b[i]-sum(L[i][0:i-1]*y[0:i-1])
		}
		y[i] = y[i] / L[i][i];			//y[i]=(b[i]-sum(L[i][0:i-1]*y[0:i-1]))/L[i][i]
	}

										//begin solving for x:Ux=y
	for (int i=n-1;i>=0;i--) {
		x[i]=y[i];						//derived algorithm xk=(yk-sum[i=k+1,n](Uki*xi))/Ukk
										//runs "backwards" from xn to x1
		for (int ii=i+1; ii<n; ii++) {
			x[i] = x[i] - U[i][ii] * x[ii];		//x[i]=(y[i]-sum(U[i][i+1:n-1]*x[i+1:n-1]))
		}
		x[i] = x[i] / U[i][i];			//x[i]=(y[i]-sum(U[i][i+1:n-1]*x[i+1:n-1]))/U[i][i]
	}
	return x;							//return solution as vector
}
