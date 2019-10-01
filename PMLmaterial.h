
//material written by  Debasis mohapatra(M.Tech.), IITR,Email Id.-debasismpt@gmail.com
#ifndef PMLmaterial_h
#define PMLmaterial_h
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
class PMLmaterial:public NDMaterial
{
public:
  PMLmaterial (int tag, double K, double G, double rho, double jita, double base);
  PMLmaterial();
  ~PMLmaterial();
    double getRho(void);
    int setTrialStrain(const Vector &v);
    const Matrix &getTangent(void);
    const Matrix &getInitialTangent(void);
    const Vector &getStress(void);
    const Vector &getStrain(void);
    int  commitState(void);
    int  revertToLastCommit(void);
    int  revertToStart(void);
    const char *getType (void) const;
	NDMaterial *getCopy(void);
    NDMaterial *getCopy(const char *code) ;
	int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker); 
    void Print(OPS_Stream &s, int flag =0);
 
  protected:
  private:  
    static Vector sigma;         // Stress vector ... class-wide for returns
    static Matrix D;	        // Elastic constants
    Vector strain_nplus1;	   //strain at (n+1)th step
    Vector strain_n;  	      //strain at nth step
	double K;                      //Bulk modulus
	double G;                     //Shear Modulus
	double rho;                  //mass density
	double jita;                //damping ratio
	double base;               //width of footing
};

 #endif