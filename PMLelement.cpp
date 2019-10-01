//PMLelement written by Mr. Debasis Mohapatra(M.Tech.)IITR,India.(Email Id.-debasismpt@gmail.com)
#include "PMLelement.h"
#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <math.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <G3Globals.h>
double PMLelement::matrixData[64];
Matrix PMLelement::M(matrixData,8,8);
Matrix PMLelement::K(matrixData,8,8);
Vector PMLelement::P(8);
Vector PMLelement::eps(3);
Vector PMLelement::peps(3);
Vector PMLelement::En(3);
Vector PMLelement::E_nminus1(3);
Vector PMLelement::Fn(3);
Vector PMLelement::F_nminus1(3);
double PMLelement::shp[3][4];
double PMLelement::pts[4][2];
double PMLelement::wts[4];
double PMLelement::FMCK[1][3];
double PMLelement::stretching[2][8];


#define OPS_Export

static int numPMLelement = 0;
OPS_Export void *
OPS_PMLelement(void)
{
  if (numPMLelement == 0) {
      opserr << "PMLelement element - Written by Debasis Mohapatra for project work\n \n";
    numPMLelement++;
   }
  // Pointer to an element that will be returned
  Element *theElement = 0;
  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();;
  
  if (numRemainingInputArgs < 15) {
    opserr << "Invalid #args, want: element element PMLelement\n";
    return 0;
  }

  int iData[6];
  char *theType;
  double dData[8];
  int numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element PMLelement\n";
    return 0;
  }
  
  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid thickness data: element PMLelement " << iData[0] << endln;
    return 0;
  }
  
  if (OPS_GetStringCopy(&theType) != 0) {
    opserr << "WARNING invalid type, want: ""PlaneStress"" or ""PlaneStrain""  element SSPquad " << iData[0] << endln;
    return 0;
  }
  
  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[5]) != 0) {
    opserr << "WARNING invalid integer data: element PMlelement\n";
    return 0;
  }
  int matID = iData[5];
  
  NDMaterial *theMaterial = OPS_GetNDMaterial(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element PMLelement " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }
  
  if (numRemainingInputArgs == 15) {
    numData = 4;
    if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
      opserr << "WARNING invalid optional data: element PMLelement " << iData[0] << endln;
      return 0;
    }
  }
  
  // parsing was successful, allocate the element
  theElement = new PMLelement(iData[0], iData[1], iData[2], iData[3],iData[4],
			 *theMaterial, theType, 
			 dData[0], dData[1], dData[2], dData[3], dData[4],dData[5] ,dData[6], dData[7]);
  
  if (theElement == 0) {
    opserr << "WARNING could not create element of type PMLelement\n";
    return 0;
  }
  
  return theElement;
}
PMLelement::PMLelement(int tag, int nd1, int nd2, int nd3, int nd4,
			   NDMaterial &m, const char *type, double t,double b,double g,double r,double l,double h,double j,double Lp)
			   :Element (tag,0), theMaterial(0), connectedExternalNodes(4), Q(8),thickness(t),applyLoad(0),rho(r), G(g),base(b),L(l),H(h),Ki(0),jita(j),Lp(Lp)
{  //Gauss points
	pts[0][0] = -0.5773502691896258;
	pts[0][1] = -0.5773502691896258;
	pts[1][0] =  0.5773502691896258;
	pts[1][1] = -0.5773502691896258;
	pts[2][0] =  0.5773502691896258;
	pts[2][1] =  0.5773502691896258;
	pts[3][0] = -0.5773502691896258;
	pts[3][1] =  0.5773502691896258;
	//weighting function
    wts[0] = 1.0;
	wts[1] = 1.0;
	wts[2] = 1.0;
	wts[3] = 1.0;
    if (strcmp(type,"PlaneStrain") != 0 && strcmp(type,"PlaneStress") != 0
	    && strcmp(type,"PlaneStrain2D") != 0 && strcmp(type,"PlaneStress2D") != 0) {
	  opserr << "PMLelement -- improper material type: " << type << "for FourNodeQuad\n";
	  exit(-1);
	}
    theMaterial = new NDMaterial *[4];
    if (theMaterial == 0) {
      opserr << "PMLelement - failed allocate material model pointer\n";
      exit(-1);
    }

	int i;
    for (i = 0; i < 4; i++) {
      // Get copies of the material model for each integration point
      theMaterial[i] = m.getCopy(type);
      // Check allocation
      if (theMaterial[i] == 0) {
	opserr << "PMLmaterial -- failed to get a copy of material model\n";
	exit(-1);
      }
    }

    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    connectedExternalNodes(2) = nd3;
    connectedExternalNodes(3) = nd4;
    
    for (i=0; i<4; i++)
      theNodes[i] = 0;
}

PMLelement::PMLelement()
:Element (0,0),
  theMaterial(0), connectedExternalNodes(4), 
 Q(8), thickness(0.0), applyLoad(0), Ki(0),G(0.0),base(0.0),L(0.0),H(0.0),Lp(0.0)
{
  pts[0][0] = -0.577350269189626;
  pts[0][1] = -0.577350269189626;
  pts[1][0] =  0.577350269189626;
  pts[1][1] = -0.577350269189626;
  pts[2][0] =  0.577350269189626;
  pts[2][1] =  0.577350269189626;
  pts[3][0] = -0.577350269189626;
  pts[3][1] =  0.577350269189626;
  
  wts[0] = 1.0;
  wts[1] = 1.0;
  wts[2] = 1.0;
  wts[3] = 1.0;

    for (int i=0; i<4; i++)
      theNodes[i] = 0;
}

PMLelement::~PMLelement()
{    
  for (int i = 0; i < 4; i++) {
    if (theMaterial[i])
      delete theMaterial[i];
  }
// Delete the array of pointers to NDMaterial pointer arrays
  if (theMaterial)
    delete [] theMaterial;
if (Ki != 0)
    delete Ki;
}

int
PMLelement::getNumExternalNodes() const
{
    return 4;
}

const ID&
PMLelement::getExternalNodes()
{
    return connectedExternalNodes;
}


Node **
PMLelement::getNodePtrs(void) 
{
  return theNodes;
}

int
PMLelement::getNumDOF()
{
    return 8;
}

void
PMLelement::setDomain(Domain *theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	theNodes[0] = 0;
	theNodes[1] = 0;
	theNodes[2] = 0;
	theNodes[3] = 0;
	return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    int Nd3 = connectedExternalNodes(2);
    int Nd4 = connectedExternalNodes(3);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
    theNodes[2] = theDomain->getNode(Nd3);
    theNodes[3] = theDomain->getNode(Nd4);
    if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 || theNodes[3] == 0) {
    return;
    }
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    int dofNd3 = theNodes[2]->getNumberDOF();
    int dofNd4 = theNodes[3]->getNumberDOF();
    
    if (dofNd1 != 2 || dofNd2 != 2 || dofNd3 != 2 || dofNd4 != 2) {
    return;
    }
    this->DomainComponent::setDomain(theDomain);
}

int
PMLelement::commitState()
{    peps=eps;
     E_nminus1=En;
	 F_nminus1=Fn;
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "PMLelement::commitState () - failed in base class";
    }    
    // Loop over the integration points and commit the material states
    for (int i = 0; i < 4; i++)
      retVal += theMaterial[i]->commitState();
	  return retVal;
}

int
PMLelement::revertToLastCommit()
{   eps=peps; 
    En=E_nminus1;
	Fn=F_nminus1;
    int retVal = 0;
    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < 4; i++)
		retVal += theMaterial[i]->revertToLastCommit();
     return retVal;
}

int
PMLelement::revertToStart()
{   eps.Zero();
    peps.Zero();
	En.Zero();
	E_nminus1.Zero();
	Fn.Zero();
	F_nminus1.Zero();
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < 4; i++)
		retVal += theMaterial[i]->revertToStart();

    return retVal;
}


int
PMLelement::update()
{
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	const Vector &disp4 = theNodes[3]->getTrialDisp();
	const Vector &vel1  = theNodes[0]->getTrialVel();
	const Vector &vel2  = theNodes[1]->getTrialVel();
	const Vector &vel3  = theNodes[2]->getTrialVel();
	const Vector &vel4  = theNodes[3]->getTrialVel();
	static double u[2][4];
	static double v[2][4];
    u[0][0] = disp1(0);
	u[1][0] = disp1(1);
	u[0][1] = disp2(0);
	u[1][1] = disp2(1);
	u[0][2] = disp3(0);
	u[1][2] = disp3(1);
	u[0][3] = disp4(0);
	u[1][3] = disp4(1);
    v[0][0] = vel1(0);
	v[1][0] = vel1(1);
	v[0][1] = vel2(0);
	v[1][1] = vel2(1);
	v[0][2] = vel3(0);
	v[1][2] = vel3(1);
	v[0][3] = vel4(0);
	v[1][3] = vel4(1);
	int ret = 0;
	// Loop over the integration points
for (int i = 0; i < 4; i++) {
     // Determine Jacobian for this integration point
	 this->shapeFunction(pts[i][0], pts[i][1]);
     double k1=stretching[0][0]/ops_Dt+stretching[0][2];
	 double k2=stretching[1][1]/ops_Dt+stretching[1][3];
	  En(0)+=E_nminus1(0)+peps(0)*ops_Dt;
	  En(1)+=E_nminus1(1)+peps(1)*ops_Dt;
	  En(2)+=E_nminus1(2)+peps(2)*ops_Dt;
	  for (int beta = 0; beta < 4; beta++) {
	  eps(0) += (stretching[0][0]*k1*k1*shp[0][beta]*v[0][beta]+stretching[0][2]*k1*k1*shp[0][beta]*u[0][beta]+stretching[0][0]*k1*stretching[0][0]*k1*peps(0)-stretching[0][2]*k1*stretching[0][2]*k1*En(0))/ops_Dt;
	  eps(1) += (stretching[1][1]*k2*k2*shp[1][beta]*v[1][beta]+stretching[1][3]*k2*k2*shp[1][beta]*u[1][beta]+stretching[1][1]*k2*stretching[1][1]*k2*peps(1)-stretching[1][3]*k2*stretching[1][3]*k2*En(1))/ops_Dt;
	  eps(2) += 2*(stretching[0][0]*k1*k2*shp[1][beta]*v[0][beta]+stretching[1][1]*k1*k2*shp[0][beta]*v[1][beta]+stretching[0][2]*k1*k2*shp[1][beta]*u[0][beta]+stretching[1][3]*k1*k2*shp[0][beta]*u[1][beta]+stretching[0][0]*stretching[1][1]*k1*k2*peps(2)-stretching[0][2]*stretching[1][3]*k1*k2*En(2))/ops_Dt;	}
      // Set the material strain
      ret += theMaterial[i]->setTrialStrain(eps);
	}

	return ret;
}
const Matrix&
PMLelement::getTangentStiff()
{

	M.Zero();
	K.Zero();
	double dvol;
	double DB[3][2];
    int k, l, i1, j1;
    double Nstiff;
    for (int m = 0; m < 4; m++){
    dvol = this->shapeFunction(pts[m][0], pts[m][1]);
    dvol *= (G*thickness*wts[m]);
    for (k = 0, i1 = 0; k < 8; k += 2, i1++) {
    for (l= 0, j1 = 0; l < 8;  l += 2, j1++) {
    Nstiff = dvol*shp[2][i1]*shp[2][j1]*FMCK[0][2];
    M(k,l) += Nstiff;
    M(k+1,l+1) += Nstiff;
    }
    }
  }
  for (int i = 0; i < 4; i++) {
     // Determine Jacobian for this integration point
        dvol = this->shapeFunction(pts[i][0], pts[i][1]);
        dvol *= (thickness*wts[i]);
        const Matrix &D = theMaterial[i]->getTangent();
		double cs=(sqrt(G/rho)*ops_Dt);
        double C=(1+2*jita*base/cs)/ops_Dt;
	    double D00 = C*D(0,0); double D01 = C*D(0,1); double D02 = C*D(0,2);
	    double D10 = C*D(1,0); double D11 = C*D(1,1); double D12 = C*D(1,2);
	    double D20 = C*D(2,0); double D21 = C*D(2,1); double D22 = C*D(2,2);
		double X1= FMCK[0][0]/ops_Dt+FMCK[0][2];
		double X2= FMCK[1][1]/ops_Dt+FMCK[1][3];
		double D1= FMCK[0][4]+ops_Dt*FMCK[0][6];
		double D2= FMCK[1][5]+ops_Dt*FMCK[0][7];
        for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
	    for (int beta = 0, ib = 0; beta < 4; beta++, ib += 2) {
	    DB[0][0] = dvol * (D00 *FMCK[0][2]*X1*X1*shp[0][beta] + D02 * FMCK[0][2]*X1*X2* shp[1][beta]);
	    DB[1][0] = dvol * (D10 *FMCK[0][2]*X1*X1*shp[0][beta] + D12 * FMCK[0][2]*X1*X2* shp[1][beta]);
	    DB[2][0] = dvol * (D20 *FMCK[0][2]*X1*X1*shp[0][beta] + D22 * FMCK[0][2]*X1*X2* shp[1][beta]);
	    DB[0][1] = dvol * (D01 *FMCK[1][3]*X2*X2*shp[1][beta] + D02 * FMCK[1][3]*X1*X2* shp[0][beta]);
	    DB[1][1] = dvol * (D11 *FMCK[1][3]*X2*X2*shp[1][beta] + D12 * FMCK[1][3] *X1*X2* shp[0][beta]);
	    DB[2][1] = dvol * (D21 *FMCK[1][3]*X2*X2*shp[1][beta] + D22 * FMCK[1][3]*X1*X2* shp[0][beta]);
	      

	    K(ia,ib) += D1*shp[0][alpha]*DB[0][0] + D2*shp[1][alpha]*DB[2][0];
	    K(ia,ib+1) += D1*shp[0][alpha]*DB[0][1] +D2*shp[1][alpha]*DB[2][1];
	    K(ia+1,ib) += D2*shp[1][alpha]*DB[1][0] + D1*shp[0][alpha]*DB[2][0];
	    K(ia+1,ib+1) +=D2*shp[1][alpha]*DB[1][1] +D1*shp[0][alpha]*DB[2][1];
	     }
	  }
	}
	K.addMatrix( 1.0,M, 1.0);
	return K;
}
const Matrix&
PMLelement::getInitialStiff()
{ K.Zero();
    double dvol;
    int k, l, i1, j1;
    double Nstiff;
    for (int m = 0; m < 4; m++){
    dvol = this->shapeFunction(pts[m][0], pts[m][1]);
    dvol *= (G*thickness*wts[m]);
    for (k = 0, i1 = 0; k < 8; k += 2, i1++) {
    for (l= 0, j1 = 0; l < 8;  l += 2, j1++) {
    Nstiff = dvol*shp[2][i1]*shp[2][j1]*FMCK[0][2];
    K(k,l) += Nstiff;
    K(k+1,l+1) += Nstiff;
    }
    }
  }
  return K;
}

const Matrix&
PMLelement::getMass()
{
  K.Zero();
  int m;
  static double rhoi[4];
  double sum = 0.0;
    for (m = 0; m < 4; m++) {
	  if (rho == 0)
	    rhoi[m] = theMaterial[m]->getRho();
	  else
	    rhoi[m] = rho;	    
	  sum += rhoi[m];
	}

    if (sum == 0.0)
	  return K;
    int k, l, i1, j1;
    double Nrho,rhodvol;
    for (int m = 0; m < 4; m++){
    rhodvol = this->shapeFunction(pts[m][0], pts[m][1]);
    rhodvol *= (rhoi[m]*thickness*wts[m]);
    for (k = 0, i1 = 0; k < 8; k += 2, i1++) {
    for (l= 0, j1 = 0; l < 8;  l += 2, j1++) {
    Nrho = rhodvol*shp[2][i1]*shp[2][j1]*FMCK[0][0];
    K(k,l) += Nrho;
    K(k+1,l+1) += Nrho;
    }
    }
  }
  return K;
}
const Matrix&
PMLelement::getDamp()
{

	M.Zero();
	K.Zero();
  int m;
  static double rhoi[4];
  double sum = 0.0;
  for (m = 0; m < 4; m++) {
	  if (rho == 0)
	    rhoi[m] = theMaterial[m]->getRho();
	  else
	    rhoi[m] = rho;	    
	    sum += rhoi[m];
	}

  if (sum == 0.0)
	return K;
    double dvol;
	double DB[3][2];
	int k, l, i1, j1;
    double Ndamp;
    for (int m = 0; m < 4; m++){
    dvol = this->shapeFunction(pts[m][0], pts[m][1]);
    dvol *= (rhoi[m]*sqrt(G/rho)*thickness*wts[m]*FMCK[0][1]);
    for (k = 0, i1 = 0; k < 8; k += 2, i1++) {
    for (l= 0, j1 = 0; l < 8;  l += 2, j1++) {
    Ndamp = dvol*shp[2][i1]*shp[2][j1]*FMCK[0][1];
    M(k,l) += Ndamp;
    M(k+1,l+1) += Ndamp;
    }
    }
  }
  for (int i = 0; i < 4; i++) {
     // Determine Jacobian for this integration point
        dvol = this->shapeFunction(pts[i][0], pts[i][1]);
        dvol *= (thickness*wts[i]);
        const Matrix &D = theMaterial[i]->getTangent();
        double cs=(sqrt(G/rho)*ops_Dt);
        double C=(1+2*jita*base/cs)/ops_Dt;
	    double D00 = C*D(0,0); double D01 = C*D(0,1); double D02 = C*D(0,2);
	    double D10 = C*D(1,0); double D11 = C*D(1,1); double D12 = C*D(1,2);
	    double D20 = C*D(2,0); double D21 = C*D(2,1); double D22 = C*D(2,2);
		double X1= FMCK[0][0]/ops_Dt+FMCK[0][2];
		double X2= FMCK[1][1]/ops_Dt+FMCK[1][3];
		double D1 = FMCK[0][4]+ops_Dt*FMCK[0][6];
		double D2 = FMCK[1][5]+ops_Dt*FMCK[0][7];
        for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
	    for (int beta = 0, ib = 0; beta < 4; beta++, ib += 2) {
	      DB[0][0] = dvol * (D00 *FMCK[0][0]*X1*X1*shp[0][beta] + D02 * FMCK[0][0]*X1*X2* shp[1][beta]);
	      DB[1][0] = dvol * (D10 *FMCK[0][0]*X1*X1*shp[0][beta] + D12 * FMCK[0][0]*X1*X2* shp[1][beta]);
	      DB[2][0] = dvol * (D20 *FMCK[0][0]*X1*X1*shp[0][beta] + D22 * FMCK[0][0]*X1*X2* shp[1][beta]);
	      DB[0][1] = dvol * (D01 *FMCK[1][1]*X2*X2*shp[1][beta] + D02 * FMCK[1][1]*X1*X2* shp[0][beta]);
	      DB[1][1] = dvol * (D11 *FMCK[1][1]*X2*X2*shp[1][beta] + D12 * FMCK[1][1] *X1*X2* shp[0][beta]);
	      DB[2][1] = dvol * (D21 *FMCK[1][1]*X2*X2*shp[1][beta] + D22 * FMCK[1][1]*X1*X2* shp[0][beta]);
	      

	      K(ia,ib) += D1*shp[0][alpha]*DB[0][0] + D2*shp[1][alpha]*DB[2][0];
	      K(ia,ib+1)+= D1*shp[0][alpha]*DB[0][1]+D2*shp[1][alpha]*DB[2][1];
	      K(ia+1,ib)+= D2*shp[1][alpha]*DB[1][0]+D1*shp[0][alpha]*DB[2][0];
	      K(ia+1,ib+1)+= D2*shp[1][alpha]*DB[1][1]+D1*shp[0][alpha]*DB[2][1];
	     }
	  }
	}
	K.addMatrix( 1.0,M, 1.0);
	return K;
}
void
PMLelement::zeroLoad(void)
{
	Q.Zero();
	applyLoad = 0;
    return;
}

int
PMLelement::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	int type;
	const Vector &data = theLoad->getData(type, loadFactor);

	if (type == LOAD_TAG_SelfWeight) {
		applyLoad = 1;
		return 0;
	} else {
		opserr << "FourNodeQuad::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
		return -1;
	} 

	return -1;
}
int 
PMLelement::addInertiaLoadToUnbalance(const Vector &accel)
{
  int i;
  static double rhoi[4];
  double sum = 0.0;
  for (i = 0; i < 4; i++) {
    rhoi[i] = theMaterial[i]->getRho();
    sum += rhoi[i];
  }
  
  if (sum == 0.0)
    return 0;
  
  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  const Vector &Raccel3 = theNodes[2]->getRV(accel);
  const Vector &Raccel4 = theNodes[3]->getRV(accel);
  
  if (2 != Raccel1.Size() || 2 != Raccel2.Size() || 2 != Raccel3.Size() ||
      2 != Raccel4.Size()) {
    opserr << "FourNodeQuad::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
    return -1;
  }
  
  static double ra[8];
  
  ra[0] = Raccel1(0);
  ra[1] = Raccel1(1);
  ra[2] = Raccel2(0);
  ra[3] = Raccel2(1);
  ra[4] = Raccel3(0);
  ra[5] = Raccel3(1);
  ra[6] = Raccel4(0);
  ra[7] = Raccel4(1);
  
  // Compute mass matrix
  this->getMass();
  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
    Q(i) += -K(i,j)*ra[j];
    }
    }
  
  return 0;
}

const Vector&
PMLelement::getResistingForce()
{
	P.Zero();

	double dvol;

	// Loop over the integration points
	for (int i = 0; i < 4; i++) {

		// Determine Jacobian for this integration point
		dvol = this->shapeFunction(pts[i][0], pts[i][1]);
		dvol *= (thickness*wts[i]);
        // Get material stress response
		const Vector &sigma = theMaterial[i]->getStress();
		 Fn(0)+=F_nminus1(0)+ops_Dt*sigma(0);
		 Fn(1)+=F_nminus1(1)+ops_Dt*sigma(1);
		 Fn(2)+=F_nminus1(2)+ops_Dt*sigma(2);

		for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
		    P(ia) += dvol*(FMCK[0][4]*shp[0][alpha]*sigma(0) + FMCK[1][5]*shp[1][alpha]*sigma(2))+(FMCK[0][5]*shp[0][alpha]*Fn(0) + FMCK[1][7]*shp[1][alpha]*Fn(2));
			
			P(ia+1) += dvol*(FMCK[1][5]*shp[1][alpha]*sigma(1) + FMCK[0][4]*shp[0][alpha]*sigma(2)+FMCK[1][7]*shp[1][alpha]*Fn(1) + FMCK[0][5]*shp[0][alpha]*Fn(2));

				}
	}
     // Subtract other external nodal loads ... P_res = P_int - P_ext
	//P = P - Q;
	P.addVector(1.0, Q, -1.0);

	return P;
}
const Vector&
PMLelement::getResistingForceIncInertia()
{
	
    const Vector &accel1 = theNodes[0]->getTrialAccel();
	const Vector &accel2 = theNodes[1]->getTrialAccel();
	const Vector &accel3 = theNodes[2]->getTrialAccel();
	const Vector &accel4 = theNodes[3]->getTrialAccel();
	static double a[8];
    a[0] = accel1(0);
	a[1] = accel1(1);
	a[2] = accel2(0);
	a[3] = accel2(1);
	a[4] = accel3(0);
	a[5] = accel3(1);
	a[6] = accel4(0);
	a[7] = accel4(1);
    // Compute the current resisting force
	this->getResistingForce();
	// Compute the mass matrix
    this->getMass();
    for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++){
      P(i) += K(i,j)*a[j];
     }
	 }
	const Vector &Vel1  = theNodes[0]->getTrialVel();
	const Vector &Vel2  = theNodes[1]->getTrialVel();
	const Vector &Vel3  = theNodes[2]->getTrialVel();
	const Vector &Vel4  = theNodes[3]->getTrialVel();
    a[0] = Vel1(0);
	a[1] = Vel1(1);
	a[2] = Vel2(0);
	a[3] = Vel2(1);
	a[4] = Vel3(0);
	a[5] = Vel3(1);
	a[6] = Vel4(0);
	a[7] = Vel4(1);
	this->getDamp();
    for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
    P(i) += M(i,j)*a[j];
    }
    }
  	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	const Vector &disp4 = theNodes[3]->getTrialDisp();
	a[0] = disp1(0);
	a[1] = disp1(1);
	a[2] = disp2(0);
	a[3] = disp2(1);
	a[4] = disp3(0);
	a[5] = disp3(1);
	a[6] = disp4(0);
	a[7] = disp4(1);
	this->getTangentStiff();
	for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
    P(i) += M(i,j)*a[j];
    }
    }
	return P;
	
}
int
PMLelement::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  int dataTag = this->getDbTag();
  
  // Quad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(8);
  data(0) = this->getTag();
  data(1) = thickness;
  data(2) = base;
  data(3) = G;
  data(4) = rho;
  data(5) = L;
  data(6) = H;
  data(7) =jita;
  
  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING FourNodeQuad::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }	      
  // Now quad sends the ids of its materials
  int matDbTag;
  
  static ID idData(12);
  
  int i;
  for (i = 0; i < 4; i++) {
    idData(i) = theMaterial[i]->getClassTag();
    matDbTag = theMaterial[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
			  theMaterial[i]->setDbTag(matDbTag);
    }
    idData(i+4) = matDbTag;
  }
  
  idData(8) = connectedExternalNodes(0);
  idData(9) = connectedExternalNodes(1);
  idData(10) = connectedExternalNodes(2);
  idData(11) = connectedExternalNodes(3);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING FourNodeQuad::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 4; i++) {
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING FourNodeQuad::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  return res;
}
int
PMLelement::recvSelf(int commitTag, Channel &theChannel,
						FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();
  // Quad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(8);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING FourNodeQuad::recvSelf() - failed to receive Vector\n";
    return res;
  }
  
  this->setTag((int)data(0));
  thickness=data(1);
  base = data(2);
  G =data(3);
  rho = data(4);
  L = data(5);
  H = data(6);
  jita = data(7);

  static ID idData(12);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING FourNodeQuad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  connectedExternalNodes(0) = idData(8);
  connectedExternalNodes(1) = idData(9);
  connectedExternalNodes(2) = idData(10);
  connectedExternalNodes(3) = idData(11);
  

  if (theMaterial == 0) {
    // Allocate new materials
    theMaterial = new NDMaterial *[4];
    if (theMaterial == 0) {
      opserr << "FourNodeQuad::recvSelf() - Could not allocate NDMaterial* array\n";
      return -1;
    }
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Allocate new material with the sent class tag
      theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
      if (theMaterial[i] == 0) {
	opserr << "FourNodeQuad::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
	return -1;
      }
      // Now receive materials into the newly allocated space
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
opserr << "NLBeamColumn3d::recvSelf() - material " << i << "failed to recv itself\n";
	return res;
      }
    }
  }

  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theMaterial[i]->getClassTag() != matClassTag) {
	delete theMaterial[i];
	theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
	if (theMaterial[i] == 0) {
opserr << "NLBeamColumn3d::recvSelf() - material " << i << "failed to create\n";
				
	  return -1;
	}
      }
      // Receive the material
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
opserr << "NLBeamColumn3d::recvSelf() - material " << i << "failed to recv itself\n";
	return res;
      }
    }
  }
  
  return res;
}
void
PMLelement::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {

    s << "#FourNodeQuad\n";
    
    int i;
    const int numNodes = 4;
    const int nstress = 3 ;
    
    for (i=0; i<numNodes; i++) {
      const Vector &nodeCrd = theNodes[i]->getCrds();
      const Vector &nodeDisp = theNodes[i]->getDisp();
      s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << endln;
     }
    
    // spit out the section location & invoke print on the scetion
    const int numMaterials = 4;

    static Vector avgStress(nstress);
    static Vector avgStrain(nstress);
    avgStress.Zero();
    avgStrain.Zero();
    for (i=0; i<numMaterials; i++) {
      avgStress += theMaterial[i]->getStress();
      avgStrain += theMaterial[i]->getStrain();
    }
    avgStress /= numMaterials;
    avgStrain /= numMaterials;

    s << "#AVERAGE_STRESS ";
    for (i=0; i<nstress; i++)
      s << avgStress(i) << " " ;
    s << endln;

    s << "#AVERAGE_STRAIN ";
    for (i=0; i<nstress; i++)
      s << avgStrain(i) << " " ;
    s << endln;

  } else {
	s << "\nFourNodeQuad, element id:  " << this->getTag() << endln;
	s << "\tConnected external nodes:  " << connectedExternalNodes;
	s << "\tthickness:  " << thickness << endln;
    theMaterial[0]->Print(s,flag);
	s << "\tStress (xx yy xy)" << endln;
	for (int i = 0; i < 4; i++)
		s << "\t\tGauss point " << i+1 << ": " << theMaterial[i]->getStress();
  }
}
int
PMLelement::displaySelf(Renderer &theViewer, int displayMode, float fact)
{   // first set the quantity to be displayed at the nodes;
    // if displayMode is 1 through 3 we will plot material stresses otherwise 0.0
	static Vector values(4);
    for (int j=0; j<4; j++)
	   values(j) = 0.0;
     if (displayMode < 4 && displayMode > 0) {
	for (int i=0; i<4; i++) {
	  const Vector &stress = theMaterial[i]->getStress();
	  values(i) = stress(displayMode-1);
	}
    }

    // now  determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    // store this information in 4 3d vectors v1 through v4
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    const Vector &end3Crd = theNodes[2]->getCrds();	
    const Vector &end4Crd = theNodes[3]->getCrds();	

    static Matrix coords(4,3);

    if (displayMode >= 0) {    
      
      const Vector &end1Disp = theNodes[0]->getDisp();
      const Vector &end2Disp = theNodes[1]->getDisp();
      const Vector &end3Disp = theNodes[2]->getDisp();
      const Vector &end4Disp = theNodes[3]->getDisp();

      for (int i = 0; i < 2; i++) {
	coords(0,i) = end1Crd(i) + end1Disp(i)*fact;
	coords(1,i) = end2Crd(i) + end2Disp(i)*fact;    
	coords(2,i) = end3Crd(i) + end3Disp(i)*fact;    
	coords(3,i) = end4Crd(i) + end4Disp(i)*fact;    
      }
    } else {
      int mode = displayMode * -1;
      const Matrix &eigen1 = theNodes[0]->getEigenvectors();
      const Matrix &eigen2 = theNodes[1]->getEigenvectors();
      const Matrix &eigen3 = theNodes[2]->getEigenvectors();
      const Matrix &eigen4 = theNodes[3]->getEigenvectors();
      if (eigen1.noCols() >= mode) {
	for (int i = 0; i < 2; i++) {
	  coords(0,i) = end1Crd(i) + eigen1(i,mode-1)*fact;
	  coords(1,i) = end2Crd(i) + eigen2(i,mode-1)*fact;
	  coords(2,i) = end3Crd(i) + eigen3(i,mode-1)*fact;
	  coords(3,i) = end4Crd(i) + eigen4(i,mode-1)*fact;
	}    
      } else {
	for (int i = 0; i < 2; i++) {
	  coords(0,i) = end1Crd(i);
	  coords(1,i) = end2Crd(i);
	  coords(2,i) = end3Crd(i);
	  coords(3,i) = end4Crd(i);
	}    
      }
    }
    
    int error = 0;

    // finally we draw the element using drawPolygon
    error += theViewer.drawPolygon (coords, values);

    return error;
}
Response*
PMLelement::setResponse(const char **argv, int argc, 
			  OPS_Stream &output)
{
  Response *theResponse =0;

  output.tag("ElementOutput");
  output.attr("eleType","PMLelement");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);
  output.attr("node3",connectedExternalNodes[2]);
  output.attr("node4",connectedExternalNodes[3]);

  char dataOut[10];
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {

    for (int i=1; i<=4; i++) {
      sprintf(dataOut,"P1_%d",i);
      output.tag("ResponseType",dataOut);
      sprintf(dataOut,"P2_%d",i);
      output.tag("ResponseType",dataOut);
    }
    
    theResponse =  new ElementResponse(this, 1, P);
  }   

  else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {
    
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 4) {
      
      output.tag("GaussPoint");
      output.attr("number",pointNum);
      output.attr("eta",pts[pointNum-1][0]);
      output.attr("neta",pts[pointNum-1][1]);

      theResponse =  theMaterial[pointNum-1]->setResponse(&argv[2], argc-2, output);
      
      output.endTag();

    }
  }
  else if ((strcmp(argv[0],"stresses") ==0) || (strcmp(argv[0],"stress") ==0)) {
    for (int i=0; i<4; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.attr("eta",pts[i][0]);
      output.attr("neta",pts[i][1]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());
      
      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma12");
      
      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 3, Vector(12));
  }

  else if ((strcmp(argv[0],"strain") ==0) || (strcmp(argv[0],"strains") ==0)) {
    for (int i=0; i<4; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.attr("eta",pts[i][0]);
      output.attr("neta",pts[i][1]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());
      
      output.tag("ResponseType","eta11");
      output.tag("ResponseType","eta22");
      output.tag("ResponseType","eta12");
      
      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 4, Vector(12));
  }

  output.endTag(); // ElementOutput

  return theResponse;
}
int 
PMLelement::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1) {

    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 3) {

    // Loop over the integration points
    static Vector stresses(12);
    int cnt = 0;
    for (int i = 0; i < 4; i++) {

      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStress();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      cnt += 3;
    }
    
    return eleInfo.setVector(stresses);
      
  } else if (responseID == 4) {

    // Loop over the integration points
    static Vector stresses(12);
    int cnt = 0;
    for (int i = 0; i < 4; i++) {

      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStrain();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      cnt += 3;
    }

    return eleInfo.setVector(stresses);
	
  } else

    return -1;
}

double 
PMLelement::shapeFunction(double xi, double eta)
{
	const Vector &nd1Crds = theNodes[0]->getCrds();
	const Vector &nd2Crds = theNodes[1]->getCrds();
	const Vector &nd3Crds = theNodes[2]->getCrds();
	const Vector &nd4Crds = theNodes[3]->getCrds();
    double x1,x2,x3,x4,y1,y2,y3,y4;
	double oneMinuseta = 1.0-eta;
	double onePluseta = 1.0+eta;
	double oneMinusxi = 1.0-xi;
	double onePlusxi = 1.0+xi;

	shp[2][0] = 0.25*oneMinusxi*oneMinuseta;	// N_1
	shp[2][1] = 0.25*onePlusxi*oneMinuseta;		// N_2
	shp[2][2] = 0.25*onePlusxi*onePluseta;		// N_3
	shp[2][3] = 0.25*oneMinusxi*onePluseta;		// N_4

	double J[2][2];

	J[0][0] = 0.25 * (-nd1Crds(0)*oneMinuseta + nd2Crds(0)*oneMinuseta +
				nd3Crds(0)*(onePluseta) - nd4Crds(0)*(onePluseta));

	J[0][1] = 0.25 * (-nd1Crds(0)*oneMinusxi - nd2Crds(0)*onePlusxi +
				nd3Crds(0)*onePlusxi + nd4Crds(0)*oneMinusxi);

	J[1][0] = 0.25 * (-nd1Crds(1)*oneMinuseta + nd2Crds(1)*oneMinuseta +
				nd3Crds(1)*onePluseta - nd4Crds(1)*onePluseta);

	J[1][1] = 0.25 * (-nd1Crds(1)*oneMinusxi - nd2Crds(1)*onePlusxi +
				nd3Crds(1)*onePlusxi + nd4Crds(1)*oneMinusxi);

	double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
	double oneOverdetJ = 1.0/detJ;
	double K[2][2];

	// K = inv(J)
	K[0][0] =  J[1][1]*oneOverdetJ;
	K[1][0] = -J[0][1]*oneOverdetJ;
	K[0][1] = -J[1][0]*oneOverdetJ;
	K[1][1] =  J[0][0]*oneOverdetJ;

    double L00 = 0.25*K[0][0];
    double L10 = 0.25*K[1][0];
    double L01 = 0.25*K[0][1];
    double L11 = 0.25*K[1][1];
	
	double L00oneMinuseta = L00*oneMinuseta;
	double L00onePluseta  = L00*onePluseta;
	double L01oneMinusxi  = L01*oneMinusxi;
	double L01onePlusxi   = L01*onePlusxi;

	double L10oneMinuseta = L10*oneMinuseta;
	double L10onePluseta  = L10*onePluseta;
	double L11oneMinusxi  = L11*oneMinusxi;
	double L11onePlusxi   = L11*onePlusxi;

	// See Cook, Malkus, Plesha p. 169 for the derivation of these terms
    shp[0][0] = -L00oneMinuseta - L01oneMinusxi;	// N_1,1
    shp[0][1] =  L00oneMinuseta - L01onePlusxi;		// N_2,1
    shp[0][2] =  L00onePluseta  + L01onePlusxi;		// N_3,1
    shp[0][3] = -L00onePluseta  + L01oneMinusxi;	// N_4,1
	
    shp[1][0] = -L10oneMinuseta - L11oneMinusxi;	// N_1,2
    shp[1][1] =  L10oneMinuseta - L11onePlusxi;		// N_2,2
    shp[1][2] =  L10onePluseta  + L11onePlusxi;		// N_3,2
    shp[1][3] = -L10onePluseta  + L11oneMinusxi;	// N_4,2
	double B,X,Y,X1,Y1,A;
	x1=nd1Crds(0); y1=nd1Crds(1);
	x2=nd2Crds(0); y2=nd2Crds(1);
	x3=nd3Crds(0) ; y3=nd3Crds(1);
	x4=nd4Crds(0) ; y4=nd4Crds(1);
	 X = shp[2][0]*x1+shp[2][1]*x2+shp[2][2]*x3+shp[2][3]*x4;
	 Y =  shp[2][0]*y1+shp[2][1]*y2+shp[2][2]*y3+shp[2][3]*y4;
	 B = abs(Y);
	 A = abs(X);
	 
	double F1ye,F1yp,F2xe,F2xp;
	Y1 =( B - H);
	if (Y1>0)
    { F1ye=10*Y1/Lp;
	  F1yp=10*Y1/Lp;}
	else 
	{F1ye=0;
	F1yp=0;}
   X1 = (A-L);
	if(X1>0)
	 { F2xe=10*X1/Lp;
	  F2xp=10*X1/Lp;}
	 else 
	 {F2xe=0;
	 F2xp=0;}
	 FMCK[0][0]=(1+F1ye)*(1+F2xe);
	 FMCK[0][1]=(1+F1ye)*F2xp+(1+F2xe)*F1yp;
	 FMCK[0][2]=F1yp*F2xp;
	 stretching[0][0]=1+F1ye;
	 stretching[0][1]=0;
	 stretching[0][2]=F1yp*sqrt(G/rho)/base;
	 stretching[0][3]=0;
	 stretching[0][4]=1+F2xe;
	 stretching[0][5]=0;
	 stretching[0][6]=F2xp*sqrt(G/rho)/base;
	 stretching[0][7]=0;
	 stretching[1][0]=0;
	 stretching[1][1]=1+F2xe;
	 stretching[1][2]=0;
	 stretching[1][3]=F2xp*sqrt(G/rho)/base;
	 stretching[1][4]=0;
	 stretching[1][5]=1+F1ye;
	 stretching[1][6]=0;
	 stretching[1][7]=F1yp*sqrt(G/rho)/base;
	 return detJ;
	}