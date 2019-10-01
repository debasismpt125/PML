//code written by Debasis Mohapatra,IITR(email id.-debasismpt@gmail.com)for PMLelement
#ifndef PMLelement_h
#define PMLelement_h
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
class Node;
class NDMaterial;
class Response;
class PMLelement: public Element
{
  public:
    PMLelement(int tag, int nd1, int nd2, int nd3, int nd4,
		 NDMaterial &m, const char *type,double t,double b,double G,double r,double L,double H,double jita,double Lp);
    PMLelement();
    ~PMLelement();
    const char *getClassType(void) const {return "PMLelement";};
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    int getNumDOF(void);
    void setDomain(Domain *theDomain);
	int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getDamp(void);     
    const Matrix &getMass(void);    
	void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);
    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);

    int getResponse(int responseID, Information &eleInformation);

  protected:
    
  private:
    // private attributes - a copy for each object of the class
    NDMaterial **theMaterial;                // pointer to the ND material objects
    ID connectedExternalNodes;               // Tags of quad nodes
    Node *theNodes[4];
	static double matrixData[64];             // array data for matrix
    static Matrix K;		            // Element stiffness, damping, and mass Matrix(ke+kt)
	static Matrix M;                     //element stiffness,damping,mass matrix       
    static Vector P;		              // Element resisting force vector
    Vector Q ;	                           // Applied nodal loads
	static  Vector eps;           //strain of current step
	static  Vector peps;          //strain of previous step
	static  Vector En;                       //summation of  eps*Delta(t)  
	static  Vector E_nminus1;             //summation of eps*Delta(t) till previous step
	static  Vector Fn;                  //summation of sigma*Delta(t)
	static  Vector F_nminus1;        //summation of sigma*Delta(t) upto previous step
	int applyLoad;      
    double thickness;	           // Element thickness
    double pressure;	          // Normal surface traction (pressure) over entire element
	double rho;                   //value of massdensity
	double G;      //value of shear modulus
	double L;          //value of horizontal distance from where pml started
	double Lp;
	double H;         //value of vertical distance from where pml started
	double base;      // base width of foundation
	double jita;
    static double shp[3][4];	// Stores shape functions and derivatives (overwritten)
    static double pts[4][2];	       // Stores quadrature points
    static double wts[4];		      // Stores quadrature weights
	static double FMCK[1][3];           //fm,fc,fk
	static double stretching[2][8];    //Fe,Fp,Fbe,Fbp
    // private member functions - only objects of this class can call these
    double shapeFunction(double xi, double eta);
	void   getpmlparameter(double xi,double eta);
    Matrix *Ki;
	Matrix *Mi;
	
};

#endif

