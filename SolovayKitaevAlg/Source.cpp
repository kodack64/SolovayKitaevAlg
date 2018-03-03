
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
const double alpha = sqrt(2. + sqrt(2.)) / 2.;
const double beta = sqrt(2. - sqrt(2.)) / 2.;
const int maxHierarchy = 5;
const int numClifford = 24;
const double pi = acos(0.) * 2;
const string gatestr[numClifford+1] = {
	"I","X","sX","sXd","sY","sYd",
	"Z","Y","ZsX","ZsXd","ZsY","ZsYd",
	"sZ","sZX","sZsX","sZsXd","sZsY","sZsYd",
	"sZd","sZdX","sZdsX","sZdsXd","sZdsY","sZdsYd","T"
};
const int dagId[numClifford] = { 0,1,3,2,5,4,6,7,8,9,10,11,18,13,22,23,21,20,12,19,17,16,14,15 };


/*
	instance uniquely represents singnle qubit unitary operator (SU(2))
	U = a*I - b*iX - c*iY - d*iZ = exp(i * 2 * arccos(a) * (b,c,d)/sqrt(b^2,c^2,d^2) \dot (X,Y,Z)  )
	where first non-zero element of (a,b,c,d) is positive

	@param I,X,Y,Z : represents coefficients of I,X,Y,Z
	@param hierarchy : required number of T gates
	@param str : construction string
	@param dagFlag : is conjugated or not
	@param buildFlag : is construction known or not
*/
struct Uop {
public:
	double I;
	double X;
	double Y;
	double Z;
	const int hierarchy = 0;
	vector<int> con;
	const bool buildFlag;

	Uop(double i, double x, double y, double z, int hie,vector<int> _con,bool _build=true) : hierarchy(hie),con(_con),buildFlag(_build) {
		I = i;
		X = x;
		Y = y;
		Z = z;
		if (i < 0) {
			I = -i;
			X = -x;
			Y = -y;
			Z = -z;
		}
		else if (i == 0.) {
			if (x < 0) {
				X = -x;
				Y = -y;
				Z = -z;
			}
			else if (x == 0.) {
				if (y < 0) {
					Y = -y;
					Z = -z;
				}
				else if (y == 0.) {
					if (z < 0) {
						Z = -z;
					}
				}
			}
		}
#ifdef _DEBUG
		double norm = sqrt(pow(i, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2));
		if (fabs(norm - 1.) > eps) {
			cout << toConstStr() << endl;
			assert(fabs(norm - 1.) < eps);
		}
#endif _DEBUG
	}

	/*
		compare operator with allowing eps
	*/
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

	void matrixForm(complex<double>* mat) const {
		complex<double> imag = complex<double>(0, 1);
		mat[0] = complex<double>(I, -Z);
		mat[1] = complex<double>(-Y, -X);
		mat[2] = complex<double>(Y, -X);
		mat[3] = complex<double>(I, Z);
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

		int nhie = hierarchy + rhs.hierarchy;
		vector<int> ncon(con);
		ncon.insert(ncon.end(), rhs.con.begin(), rhs.con.end());
		bool nbuild = buildFlag & rhs.buildFlag;
		return Uop(ni, nx, ny, nz, nhie, ncon, nbuild);
	}
	virtual Uop dag() const{
		vector<int> ncon;
		for (auto ite = ncon.rbegin(); ite != ncon.rend(); ite++) {
			if ((*ite) == numClifford) {
				// Tdag = XTX
				ncon.push_back(1);
				ncon.push_back(24);
				ncon.push_back(1);
			}
			else {
				ncon.push_back(dagId[*ite]);
			}
		}
		return Uop(I, -X, -Y, -Z, hierarchy, ncon, buildFlag);
	}
	string toConstStr() const {
		string s = "";
		for (auto ite = con.begin(); ite != con.end(); ite++) {
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
		if(buildFlag) cout << "Seq: " << toConstStr() << endl;
		else cout << "No build insturction" << endl;
		cout << " ***** " << endl;
	}
	virtual ~Uop() {};
};

void addClifford(const Uop& u,vector<Uop>& list) {
	double I, X, Y, Z;
	I = u.I;
	X = u.X;
	Y = u.Y;
	Z = u.Z;
	int hie = u.hierarchy;
	vector<int> con(u.con);

	list.push_back(Uop(I, X, Y, Z, hie,con)); // I
	con.push_back(1);	list.push_back(Uop(X, -I, Z, -Y, hie, con)); con.pop_back();  // iX
	con.push_back(2);	list.push_back(Uop((I + X) / sqrt2, (X - I) / sqrt2, (Y + Z) / sqrt2, (Z - Y) / sqrt2, hie,con)); con.pop_back();   // I+iX
	con.push_back(3);	list.push_back(Uop((I - X) / sqrt2, (X + I) / sqrt2, (Y - Z) / sqrt2, (Z + Y) / sqrt2, hie, con)); con.pop_back();   // I-iX
	con.push_back(4);	list.push_back(Uop((I + Y) / sqrt2, (X - Z) / sqrt2, (Y - I) / sqrt2, (Z + X) / sqrt2, hie, con)); con.pop_back();   // I+iY
	con.push_back(5);	list.push_back(Uop((I - Y) / sqrt2, (X + Z) / sqrt2, (Y + I) / sqrt2, (Z - X) / sqrt2, hie, con)); con.pop_back();   // I-iY
	assert(list.size() == 6);
	for (int i = 0; i < 6; i++) {
		double a, b, c, d;
		a = list[i].I;
		b = list[i].X;
		c = list[i].Y;
		d = list[i].Z;
		con.push_back(6);	list.push_back(Uop(d, c, -b, -a, hie, con)); con.pop_back(); // iZ
		con.push_back(12);	list.push_back(Uop((a + d) / sqrt2, (b + c) / sqrt2, (c - b) / sqrt2, (d - a) / sqrt2, hie, con)); con.pop_back(); // I+iZ
		con.push_back(18);	list.push_back(Uop((a - d) / sqrt2, (b - c) / sqrt2, (c + b) / sqrt2, (d + a) / sqrt2, hie, con)); con.pop_back(); // I-iZ
	}
	assert(list.size() == 24);
}

void genreateEpsilonNetwork(set<Uop>& s) {
	queue<Uop> q;
	vector<Uop> list;
	int currentHierarchy = 0;
	Uop I(1, 0, 0, 0, 0, {});
	addClifford(I, list);
	for (unsigned int i = 0; i < list.size(); i++) {
		if (s.find(list[i]) == s.end()) {
			s.insert(list[i]);
			if (list[i].hierarchy<maxHierarchy) q.push(list[i]);
		}
	}
	while (!q.empty()) {
		double a, b, c, d;
		list.clear();
		double ti, tx, ty, tz;
		a = q.front().I;
		b = q.front().X;
		c = q.front().Y;
		d = q.front().Z;
		ti = a*alpha + d*beta;
		tx = b*alpha + c*beta;
		ty = c*alpha - b*beta;
		tz = d*alpha - a*beta;
		vector<int> ncon(q.front().con);
		ncon.push_back(24);
		addClifford(Uop(ti, tx, ty, tz, q.front().hierarchy + 1, ncon), list);

		if (q.front().hierarchy + 1 > currentHierarchy) {
			cout << currentHierarchy << " " << s.size() << endl;
			currentHierarchy++;
		}
		for (unsigned int i = 0; i < list.size(); i++) {
			if (s.find(list[i]) == s.end()) {
				s.insert(list[i]);
				if (list[i].hierarchy<maxHierarchy) q.push(list[i]);
			}
		}
		q.pop();
	}

	/*
	// output all operators in epsilon net
	ofstream ofs[maxHierarchy + 1];
	for (int hie = 0; hie <= maxHierarchy; hie++) {
		ofs[hie].open("result_" + to_string(hie) + ".txt");
	}
	for (auto ite = s.begin(); ite != s.end(); ite++) {
		ofs[ite->hierarchy] << ite->I << " " << ite->X << " " << ite->Y << " " << ite->Z << endl;
	}
	for (int hie = 0; hie <= maxHierarchy; hie++) {
		ofs[hie].close();
	}
	*/
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
			con = (*ite).con;
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
		decompose udd = S V W V.dag W.dag S.dag
			where	
					udd	= a I - i b X - i c Y - i d Z
					a	= cos(theta/2)
					V	= cos(phi/2) I - i sin(phi/2) X
					W	= cos(phi/2) I - i sin(phi/2) Y

		VWV.dagW.dag = (1-2s^4, 2s^3c,-2s^3c, 2c^2s^2)
		VWV.dagW.dag = (1-2s^4, 2s^3c,-2s^3c, -2c^2s^2) ?

		a = 1-2s^4
		s = pow((1-a)/2,1/4)
		
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

		// udd.print();
		//cout << nx << " " << ny << " " << nz << endl;


		mn = sqrt(1.-pow(udd.I,2));
		mx = 2 * pow(s, 3)*c/mn;
		my = -2 * pow(s, 3)*c/mn;
		mz = -2 * pow(s, 2)*pow(c, 2)/mn;

		//(V*W*V.dag()*W.dag()).print();
		//cout << mx << " " << my << " " << mz << endl;

		x = (nx + mx) / 2;
		y = (ny + my) / 2;
		z = (nz + mz) / 2;
		n = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
		x /= n;
		y /= n;
		z /= n;

		Uop S(0, x, y, z, 0, {}, false);
		Uop Vt = S*V*S.dag();
		Uop Wt = S*W*S.dag();
		//cout << x << " " << y << " " << z << endl;
		//Uop u = V*W*V.dag()*W.dag();

		//udd.print();
		//(Vt*Wt*Vt.dag()*Wt.dag()).print();

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
	Uop appr0 = solovayKitaev(U0set, sqT, 0);
	cout << operatorDist(sqT, appr0) << " "; appr0.print();
	Uop appr1 = solovayKitaev(U0set, sqT, 1);
	cout << operatorDist(sqT, appr1) << " "; appr1.print();
	Uop appr2 = solovayKitaev(U0set, sqT, 2);
	cout << operatorDist(sqT, appr2) << " "; appr2.print();
	Uop appr3 = solovayKitaev(U0set, sqT, 3);
	cout << operatorDist(sqT, appr3) << " "; appr3.print();
	Uop appr4 = solovayKitaev(U0set, sqT, 4);
	cout << operatorDist(sqT, appr4) << " "; appr4.print();
	Uop appr5 = solovayKitaev(U0set, sqT, 5);
	cout << operatorDist(sqT, appr5) << " "; appr5.print();
}
