#ifndef __Random__
#define __Random__

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // methods
  void Set();
  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
  double Lorentz(double mean, double width);
  double Exponential(double lambda);
  void random_unif(double * r, int n);
  void random_exp(double * r, int n, double lambda);
  void random_lorentz(double * r, int n, double mean, double width);
  void random_gauss(double * r, int n, double mean, double sigma);
  double error(double *AV, double *AV2, int n);
  void angolo_rand(double ** omega, int dim);
  double psi110(double * vettore);
  double psi210(double * vettore);
};

#endif // __Random__
