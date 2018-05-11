
#include <iostream>
#include <fstream>
#include <cassert>
#include <set>
#include <vector>
#include <cmath>
#include <queue>
#include <string>
#include <complex>

using namespace std;

const double eps = 1e-12;
const double sqrt2 = sqrt(2);
const double alpha = sqrt(2. + sqrt2) / 2.;
const double beta = sqrt(2. - sqrt2) / 2.;
const int numClifford = 24;
const double pi = acos(0.) * 2;
const string gatestr[] = {"I","X","(I+iX)","(I-iX)","(I+iY)","(I-iY)","Z","(I+iZ)","(I-iZ)","T"};
const int dagId[] =		  { 0, 1, 3, 2, 5, 4, 6, 8, 7, -1 };
const int TID = 9;

// allowed number of T-gate in epsilon net
const int maxHierarchy = 3;

// maximum recursive number of Solovay-Kitaev
const int maxRecursiveSK = 10;

// If ture, normalize when unitary operator is defined
const bool normalize = true;

// If ture, print construction of gates with Clifford gates and T-gates.
const bool showConstruction = false;

struct Uop {
public:
	double I;
	double X;
	double Y;
	double Z;

	// required T-gate
	const int hierarchy = 0;

	// construction of gate with cliffords and T-gates
	vector<int> construction;

	// True if construction is known
	const bool buildFlag;

	Uop(double _I, double _X, double _Y, double _Z, int _hierarchy, vector<int> _construction,bool _buildFlag=true) 
		: hierarchy(_hierarchy),construction(_construction),buildFlag(_buildFlag) {

		// to omit double-counting in SU(2), every variable is inverted if the first non-zero value is negative
		double v[4] = { _I,_X,_Y,_Z };
		int mi = 0;
		for (int i = 0; i < 4; i++) if (abs(v[i]) > eps) { mi = i; break; }
		if (v[mi]>0) {
			I = _I;
			X = _X;
			Y = _Y;
			Z = _Z;
		}
		else {
			I = -_I;
			X = -_X;
			Y = -_Y;
			Z = -_Z;
		}
		if (normalize) {
			double norm = sqrt(pow(_I, 2) + pow(_X, 2) + pow(_Y, 2) + pow(_Z, 2));
			I /= norm;
			X /= norm;
			Y /= norm;
			Z /= norm;
		}
#ifdef _DEBUG
		if(!normalize) {
			double norm = sqrt(pow(_I, 2) + pow(_X, 2) + pow(_Y, 2) + pow(_Z, 2));
			if (fabs(norm - 1.) > eps) {
				cout << toConstStr() << endl;
				assert(fabs(norm - 1.) < eps);
			}
		}
#endif _DEBUG
	}

	// compare operator with allowing eps
	bool operator<(const Uop& rhs) const {
		if (I < rhs.I-eps) return true;
		else if (fabs(I - rhs.I)<eps) {
			if (X < rhs.X-eps) return true;
			else if (fabs(X - rhs.X)<eps) {
				if (Y < rhs.Y-eps) return true;
				else if (fabs(Y-  rhs.Y)<eps) {
					if (Z < rhs.Z-eps) return true;
				}
			}
		}
		return false;
	}

	/*
	 convert to matrix form:
	 [ [mat[0], mat[1]]
	   [mat[2], mat[3]] ]
	*/
	void matrixForm(complex<double>* mat) const {
		complex<double> imag = complex<double>(0, 1);
		mat[0] = complex<double>(I, Z);
		mat[1] = complex<double>(Y, X);
		mat[2] = complex<double>(-Y, X);
		mat[3] = complex<double>(I, -Z);
	}
	Uop operator*(const Uop& rhs) const {
		double ni, nx, ny, nz;
		double a, b, c, d;
		double p, q, r, s;
		a = I;
		b = X;
		c = Y;
		d = Z;
		p = rhs.I;
		q = rhs.X;
		r = rhs.Y;
		s = rhs.Z;
		ni = a*p - b*q - c*r - d*s;
		nx = a*q + b*p - c*s + d*r;
		ny = a*r + b*s + c*p - d*q;
		nz = a*s - b*r + c*q + d*p;

		int newHierarchy = hierarchy + rhs.hierarchy;
		vector<int> newConstruction(construction);
		newConstruction.insert(newConstruction.end(), rhs.construction.begin(), rhs.construction.end());
		bool newBuildFlag = buildFlag & rhs.buildFlag;
		Uop newOp(ni, nx, ny, nz, newHierarchy, newConstruction, newBuildFlag);

#ifdef _DEBUG
		// check product processing
		complex<double> m0[4];
		complex<double> m1[4];
		complex<double> m2[4];
		complex<double> m3[4];
		this->matrixForm(m0);
		rhs.matrixForm(m1);
		newOp.matrixForm(m2);
		m3[0] = m0[0] * m1[0] + m0[1] * m1[2];
		m3[1] = m0[0] * m1[1] + m0[1] * m1[3];
		m3[2] = m0[2] * m1[0] + m0[3] * m1[2];
		m3[3] = m0[2] * m1[1] + m0[3] * m1[3];
		bool f1, f2;
		f1 = f2 = true;
		for (int i = 0; i < 4; i++) {
			f1 &= abs(m3[i] - m2[i]) < 1e-10;
			f2 &= abs(m3[i] + m2[i]) < 1e-10;
		}
		if ((!f1) && (!f2)) {
			this->print();
			rhs.print();
			newOp.print();
			assert(false);
		}
#endif
		return newOp;
	}

	// create Hermitian conjugate
	virtual Uop dag() const{
		vector<int> ncon;
		for (auto ite = this->construction.rbegin(); ite != this->construction.rend(); ++ite) {
			if ((*ite) == TID) {
				// T^dag = XTX
				ncon.push_back(1);
				ncon.push_back(TID);
				ncon.push_back(1);
			}
			else {
				ncon.push_back(dagId[*ite]);
			}
		}
		Uop newOp(I, -X, -Y, -Z, hierarchy, ncon, buildFlag);
#ifdef _DEBUG
		newOp.verifyConstruction();
#endif
		return newOp;
	}
	string toConstStr() const {
		string s = "";
		for (auto ite = construction.begin(); ite != construction.end(); ++ite) {
			s += gatestr[*ite]+" ";
		}
		if (s.length() == 0) s = "I";
		return s;
	}
	virtual void print() const{
		cout << " ***** " << endl;
		cout << "#T gate : " << hierarchy << endl;
		cout << I << " I + " << X << " iX + " << Y << " iY + " << Z << " iZ" << endl;
		double nn;
		nn = sqrt(pow(X,2)+pow(Y,2)+pow(Z,2));
		cout << "rot = " << 2*acos(I)/pi <<  " * pi, axis = (" << X/nn << "," << Y/nn << "," << Z/nn << ")" << endl;
		if (showConstruction) {
			if (buildFlag) cout << "Seq: " << toConstStr() << endl;
			else cout << "No build insturction" << endl;
		}
		cout << " ***** " << endl;
	}
#ifdef _DEBUG
	// verify its construction regenerates operator
	virtual void verifyConstruction() const {
		if (!this->buildFlag) return;
		Uop oI(1, 0, 0, 0, 0, {}, true);
		Uop oX(0, 1, 0, 0, 0, {1}, true);
		Uop osX(1 / sqrt2, 1 / sqrt2, 0, 0, 0, {2}, true);
		Uop osXd(1 / sqrt2, -1 / sqrt2, 0, 0, 0, {3}, true);
		Uop osY(1 / sqrt2, 0, 1 / sqrt2, 0, 0, {4}, true);
		Uop osYd(1 / sqrt2, 0, -1 / sqrt2, 0, 0, {5}, true);
		Uop oZ(0, 0, 0, 1, 0, {6}, true);
		Uop osZ(1 / sqrt2, 0, 0, 1 / sqrt2, 0, {7}, true);
		Uop osZd(1 / sqrt2, 0, 0, -1 / sqrt2, 0, {8}, true);
		Uop T(cos(pi/8.), 0, 0, sin(pi/8.), 0, {9}, true);
		Uop gates[] = {oI, oX, osX, osXd, osY, osYd, oZ, osZ, osZd, T};

		vector<Uop> ops;
		vector<int> cons;
		ops.push_back(Uop(1, 0, 0, 0, 0, {}, true));
		for (int i = 0; i < this->construction.size(); i++) {
			Uop op = ops[ops.size()-1]*gates[this->construction[i]];
			ops.push_back(op);
		}
		if (!(  (ops[ops.size() - 1] < (*this)) == false 
			&& ((*this) < ops[ops.size() - 1])==false  )) {
			this->print();
			ops[ops.size() - 1].print();
			assert(false);
		}
		return;
	}
#endif
virtual ~Uop() {};
};

// create length24 list of [(u*C) for C in Clifford]
// Clifford can be constructed with product((I,iX,I+iX,I-iX,I+iY,I-iY) , (I,iZ,I+iZ,I-iZ))
void setClifford(const Uop& u,vector<Uop>& list) {
	double I, X, Y, Z;
	I = u.I;
	X = u.X;
	Y = u.Y;
	Z = u.Z;
	int hie = u.hierarchy;
	vector<int> con(u.construction);

	list.clear();
	list.push_back(Uop(I, X, Y, Z, hie,con)); // I
	con.push_back(1);	list.push_back(Uop(-X, I, -Z, Y, hie, con)); con.pop_back();  // iX
	con.push_back(2);	list.push_back(Uop((I - X) / sqrt2, (X + I) / sqrt2, (Y - Z) / sqrt2, (Z + Y) / sqrt2, hie, con)); con.pop_back();   // I+iX
	con.push_back(3);	list.push_back(Uop((I + X) / sqrt2, (X - I) / sqrt2, (Y + Z) / sqrt2, (Z - Y) / sqrt2, hie, con)); con.pop_back();   // I-iX
	con.push_back(4);	list.push_back(Uop((I - Y) / sqrt2, (X + Z) / sqrt2, (Y + I) / sqrt2, (Z - X) / sqrt2, hie, con)); con.pop_back();   // I+iY
	con.push_back(5);	list.push_back(Uop((I + Y) / sqrt2, (X - Z) / sqrt2, (Y - I) / sqrt2, (Z + X) / sqrt2, hie, con)); con.pop_back();   // I-iY
	assert(list.size() == 6);
	for (int i = 0; i < 6; i++) {
		double a, b, c, d;
		a = list[i].I;
		b = list[i].X;
		c = list[i].Y;
		d = list[i].Z;
		if (i != 0) con.push_back(i);
		con.push_back(6);	list.push_back(Uop(-d, -c, b, a, hie, con)); con.pop_back(); // iZ
		con.push_back(7);	list.push_back(Uop((a - d) / sqrt2, (b - c) / sqrt2, (c + b) / sqrt2, (d + a) / sqrt2, hie, con)); con.pop_back(); // I+iZ
		con.push_back(8);	list.push_back(Uop((a + d) / sqrt2, (b + c) / sqrt2, (c - b) / sqrt2, (d - a) / sqrt2, hie, con)); con.pop_back(); // I-iZ
		if (i != 0) con.pop_back();
	}
	assert(list.size() == numClifford);
}

// enumerate every unitary operator constructed with T-gates and Clifford gates where #T-gate < maxHierarchy.
void genreateEpsilonNetwork(set<Uop>& s) {
	queue<Uop> q;
	vector<Uop> list;
	int currentHierarchy = 0;
	Uop I(1, 0, 0, 0, 0, {});

	setClifford(I, list);
	for (unsigned int i = 0; i < list.size(); i++) {
		s.insert(list[i]);
		q.push(list[i]);
	}
	while (!q.empty()) {
		double a, b, c, d;
		double ti, tx, ty, tz;
		a = q.front().I;
		b = q.front().X;
		c = q.front().Y;
		d = q.front().Z;
		ti = a*alpha - d*beta;
		tx = b*alpha - c*beta;
		ty = c*alpha + b*beta;
		tz = d*alpha + a*beta;

		vector<int> ncon(q.front().construction);
		ncon.push_back(TID);
		setClifford(Uop(ti, tx, ty, tz, q.front().hierarchy + 1, ncon), list);

		if (q.front().hierarchy + 1 > currentHierarchy) {
			cout << currentHierarchy << " " << s.size() << endl;
			currentHierarchy++;
		}
		for (unsigned int i = 0; i < list.size(); i++) {
			if (s.insert(list[i]).second) {
#ifdef _DEBUG
				list[i].verifyConstruction();
#endif
				if (list[i].hierarchy < maxHierarchy) {
					q.push(list[i]);
				}
			}
		}
		q.pop();
	}
}

/*
	for given matrix A = [ [ mat[0],mat[1] ] , [mat[2],mat[3]]],
	obtain max eigenvalue of A^dag A
*/
double maxeigen(complex<double>* mat) {
	double p, q, r2;
	p = norm(mat[0]) + norm(mat[2]);
	q = norm(mat[1]) + norm(mat[3]);
	r2 = norm(mat[0] * conj(mat[1]) + conj(mat[2])*mat[3]);
	return (p + q + sqrt(pow(p - q, 2) + 4 * r2))/2;
}

double operatorNorm(const Uop& u) {
	complex<double> um[4];
	u.matrixForm(um);
	return sqrt(maxeigen(um));
}

// for given u and v, operatorDist(u,v) = operatorNorm(u-v)
double operatorDist(const Uop& u, const Uop& v) {
	complex<double> um[4],vm[4];
	u.matrixForm(um);
	v.matrixForm(vm);
	for (int i = 0; i < 4; i++) um[i] -= vm[i];
	return sqrt(maxeigen(um));
}

Uop getSimilar(const set<Uop> U0set, const Uop& u) {
	double mindist = 100;
	double i, x, y, z;
	int hie;
	vector<int> con;
	for (auto ite = U0set.begin(); ite != U0set.end(); ite++) {
		double dist = operatorDist(u, *ite);
		if (mindist > dist) {
			mindist = dist;
			i = (*ite).I;
			x = (*ite).X;
			y = (*ite).Y;
			z = (*ite).Z;
			hie = (*ite).hierarchy;
			con = (*ite).construction;
		}
	}
	return Uop(i,x,y,z,hie,con,true);
}

Uop solovayKitaev(const set<Uop> U0set, const Uop& u, int rec) {
	if (rec == 0) {
		return getSimilar(U0set, u);
	}
	else {
		Uop ud = solovayKitaev(U0set, u, rec - 1);
		Uop udd = u*ud.dag();
		/*

		SK

		for given u, rec

		ud = SK(u,rec-1)
		udd = u ud
		decompose udd = S V W V.dag W.dag S.dag
			where	
					udd	= a I - i b X - i c Y - i d Z
					a	= cos(theta/2)
					V	= cos(phi/2) I - i sin(phi/2) X
					W	= cos(phi/2) I - i sin(phi/2) Y

		V W V.dag W.dag = (1-2s^4, 2s^3c,-2s^3c, -2c^2s^2)

		a = 1-2s^4 <==>	s = pow((1-a)/2,1/4)

		Vt = S V S.dag
		Wt = S W S.dag

		udd = Vt Wt Vt.dag Wt.dag
		Vd = SK(Vt,rec-1)
		Wd = SK(Wt,rec-1)

		return Vd Wd Vd.ag Wd.dag ud
		
		*/
		double s = pow((1. - udd.I) / 2, 0.25);
		double c = sqrt(1. - pow(s, 2));
		Uop V(c, s, 0, 0, 0, {}, false);
		Uop W(c, 0, s, 0, 0, {}, false);

		double mx, my, mz,mn;
		double nx, ny, nz,nn;
		double x, y, z,n;
		nn = sqrt(1 - pow(udd.I,2));
		nx = udd.X/nn;
		ny = udd.Y/nn;
		nz = udd.Z/nn;
		mn = sqrt(1.-pow(udd.I,2));
		mx = 2 * pow(s, 3)*c/mn;
		my = -2 * pow(s, 3)*c/mn;
		mz = -2 * pow(s, 2)*pow(c, 2)/mn;
		x = (nx + mx) / 2;
		y = (ny + my) / 2;
		z = (nz + mz) / 2;
		n = sqrt(x*x+y*y+z*z);
		x /= n;
		y /= n;
		z /= n;

		Uop S(0, x, y, z, 0, {}, false);
		Uop Vt = S*V*S.dag();
		Uop Wt = S*W*S.dag();

#ifdef _DEBUG
		double dist = operatorDist(udd, Vt*Wt*Vt.dag()*Wt.dag());
		if (dist > 1e-10) {
			udd.print();
			(Vt*Wt*Vt.dag()*Wt.dag()).print();
			cout << dist << endl;
			assert(false);
		}
#endif

		Uop vd = solovayKitaev(U0set, Vt, rec - 1);
		Uop wd = solovayKitaev(U0set, Wt, rec - 1);
		return vd*wd*vd.dag()*wd.dag()*ud;
	}
}

int main() {

	set<Uop> U0set;
	genreateEpsilonNetwork(U0set);
	Uop sqT(cos(pi / 16), 0, 0, sin(pi / 16), 0, {},false);
	sqT.print();

	for (int i = 0; i <= maxRecursiveSK; i++) {
		Uop appr = solovayKitaev(U0set, sqT, i);
		cout << "***********************" << endl << "recursive: " << i << endl;
		cout << "dist:" << operatorDist(sqT, appr) << endl;
		appr.print();
#ifdef _DEBUG
		appr.verifyConstruction();
#endif
	}
}
