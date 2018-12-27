#ifndef _FMM2DTree_HPP__
#define _FMM2DTree_HPP__

#include <vector>
#include <Eigen/Dense>
#include <fstream>

#define EIGEN_DONT_PARALLELIZE
#include <iostream>
#include <cmath>
#include <stdlib.h>     /* srand, rand */
#include <time.h>

#include "domain2D.hpp"
long double epsilon_f	=	1e-6;
long double epsilon	=	1e-9;

using namespace std;

const double PI	=	3.1415926535897932384;

struct pts2D {
	double x,y;
};

struct charge {
	double q,x,y;
};

struct level_box {
	int level,box;
};

class FMM2DBox {
public:
	bool exists;
	int boxNumber;
	int parentNumber;
	int childrenNumbers[4];
	int neighborNumbers[8];
	int innerNumbers[16];
	int outerNumbers[24];

	std::vector<level_box> list1;		//	Neighbors or boxes that share a boundary
	//std::vector<level_box> list2;		//	Interactions at same level
	std::vector<level_box> list3;		//	Descendants of neighbors at the same level, which are not neighbors
	std::vector<level_box> list4;		//	Opposite of List 3

	FMM2DBox () {
		boxNumber		=	-1;
		parentNumber	=	-1;
		for (int l=0; l<4; ++l) {
			childrenNumbers[l]	=	-1;
		}
		for (int l=0; l<8; ++l) {
			neighborNumbers[l]	=	-1;
		}
		for (int l=0; l<16; ++l) {
			innerNumbers[l]		=	-1;
		}
		for (int l=0; l<24; ++l) {
			outerNumbers[l]		=	-1;
		}		
	}
	
	Eigen::VectorXd multipoles;
	Eigen::VectorXd locals;

	pts2D center;

	//	The following will be stored only at the leaf nodes
	std::vector<pts2D> chebNodes;
};


class user_Greens_function: public Greens_function {
public:

	// inputs for integrand  chebnodes
	pts2D chebNode_i;
	pts2D chebNode_j;
	int nChebNodes;
	void initialise(pts2D chebNode_i1, pts2D chebNode_j1, int nChebNodes1, double xcenter1, double ycenter1, double L1, int nNodes1) {
		xcenter	=	xcenter1;
		ycenter	=	ycenter1;
		Lx	=	L1;
		Ly	=	L1;
		nNodes	=	nNodes1;

		chebNode_i	=	chebNode_i1;
		chebNode_j	=	chebNode_j1;
		nChebNodes	=	nChebNodes1;

	}
	
	#ifdef ONEOVERR
	user_Greens_function() {
		isTrans		=	true;
		isHomog		=	true;
		isLogHomog	=	false;
		alpha		=	-1.0;
	};

	

	long double integrand(long double x, long double y) {
		long double R = (sqrt(pow(chebNode_i.x-x,2)+pow(chebNode_i.y-y,2)));
		if (R < 1e-10) {
			return 0.0;
		}
		else {
			return (1/R)  * get_S(chebNode_j.x, x, nChebNodes) * get_S(chebNode_j.y, y, nChebNodes);
		}
	};
	#elif LOGR
	user_Greens_function() {
		isTrans		=	true;
		isHomog		=	false;
		isLogHomog	=	true;
		alpha		=	1.0;
	};
	long double integrand(long double x, long double y) {
		long double R = (chebNode_i.x-x)*(chebNode_i.x-x)+(chebNode_i.y-y)*(chebNode_i.y-y);
		return (0.5*log(R) * get_S(chebNode_j.x, (x-xcenter)/Lx, nChebNodes) * get_S(chebNode_j.y, (y-ycenter)/Ly, nChebNodes));
	};
	#endif

	//	get_ChebPoly
	double get_ChebPoly(double x, int n) {
		return cos(n*acos(x));
	}

	//	get_S
	double get_S(double x, double y, int n) {
		double S	=	0.5;
		for (int k=1; k<n; ++k) {
			S+=get_ChebPoly(x,k)*get_ChebPoly(y,k);
		}
		return 2.0/n*S;
	}

	// RHS function
	double f(double x, double y) {
		double sigma = 0.5;
		return exp(-(x*x+y*y)/(2*sigma*sigma))/(2*PI*sigma*sigma);
		//return (1.0); //constant function
	}

	~user_Greens_function() {
		//cout << "D derived called" << endl;
	};

};




class user_function2: public Greens_function {
public:

	// inputs for integrand  chebnodes
	pts2D chebNode_j;
	int nChebNodes;
	

	user_function2() {
		xcenter	=	0.0;
		ycenter	=	0.0;
		Lx	=	1.0;
		Ly	=	1.0;
		
	};
	void initialise (pts2D chebNode_j1, int nChebNodes1, int nNodes1) {
		nNodes	=	nNodes1;

		chebNode_j	=	chebNode_j1;
		nChebNodes	=	nChebNodes1;
	}
	long double integrand(long double x, long double y) {
		return (long double)(get_S(chebNode_j.x, x, nChebNodes) * get_S(chebNode_j.y, y, nChebNodes));
	}


	//	get_ChebPoly
	double get_ChebPoly(double x, int n) {
		return cos(n*acos(x));
	}

	//	get_S
	double get_S(double x, double y, int n) {
		double S	=	0.5;
		for (int k=1; k<n; ++k) {
			S+=get_ChebPoly(x,k)*get_ChebPoly(y,k);
		}
		return 2.0/n*S;
	}

	~user_function2() {};

};



class user_function3: public Greens_function {
public:

	// inputs for integrand  chebnodes
	pts2D chebNode_i;
	int nChebNodes;
	

	user_function3() {
		
	};
	
	void initialise(pts2D chebNode_i1, long double xcenter1, long double ycenter1, long double L1, int nNodes1) {
		xcenter	=	xcenter1;
		ycenter	=	ycenter1;
		Lx	=	L1;
		Ly	=	L1;
		nNodes	=	nNodes1;
		chebNode_i	=	chebNode_i1;
	}

	long double integrand(long double x, long double y) {
		long double R = pow(chebNode_i.x-x,2)+pow(chebNode_i.y-y,2);
		return 0.5*log(R) * f(x,y) ;
	}

	// RHS function
	double f(double x, double y) {
		double sigma = 0.5;
		return exp(-(x*x+y*y)/(2*sigma*sigma))/(2*PI*sigma*sigma);
		//return (1); //constant function
	}
	~user_function3() {};

};

template <typename Greens_functiontype>
class FMM2DTree {
public:
	Greens_functiontype* K;
	int nLevels;			//	Number of levels in the tree.
	int nChebNodes;			//	Number of Chebyshev nodes along one direction.
	int rank;				//	Rank of interaction, i.e., rank = nChebNodes*nChebNodes.
	int N;					//	Number of particles.
	double L;				//	Semi-length of the simulation box.
	double smallestBoxSize;	//	This is L/2.0^(nLevels).
	double a;				//	Cut-off for self-interaction. This is less than the length of the smallest box size.

	std::vector<int> nBoxesPerLevel;			//	Number of boxes at each level in the tree.
	std::vector<double> boxRadius;				//	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
	std::vector<double> boxHomogRadius;			//	Stores the value of boxRadius^{alpha}
	std::vector<double> boxLogHomogRadius;		//	Stores the value of alpha*log(boxRadius)
	std::vector<std::vector<FMM2DBox> > tree;	//	The tree storing all the information.
	
	//	childless_boxes
	std::vector<level_box> childless_boxes;

	//	charge in coloumbs and its location
//	std::vector<charge> charge_database;
  	
	//	Chebyshev nodes
	std::vector<double> standardChebNodes1D;
	std::vector<pts2D> standardChebNodes;
	std::vector<pts2D> standardChebNodesChild;
	std::vector<pts2D> leafChebNodes;

	//	Different Operators
	Eigen::MatrixXd selfInteraction;		//	Needed only at the leaf level.
	Eigen::MatrixXd neighborInteraction[8];	//	Neighbor interaction only needed at the leaf level.
	Eigen::MatrixXd M2M[4];					//	Transfer from multipoles of 4 children to multipoles of parent.
	Eigen::MatrixXd L2L[4];					//	Transfer from locals of parent to locals of 4 children.
	Eigen::MatrixXd M2LInner[16];			//	M2L of inner interactions. This is done on the box [-L,L]^2.
	Eigen::MatrixXd M2LOuter[24];			//	M2L of outer interactions. This is done on the box [-L,L]^2.
	Eigen::RowVectorXd M2L_term2;			// same for al M2Ls
	Eigen::MatrixXd S_mat;

	//	the potential evaluated by interpolating the potential at the locals and adding the near field
	Eigen::VectorXd potential;

// public:
	FMM2DTree(Greens_functiontype* K, /*int nLevels,*/ int nChebNodes, double L) {
		this->K						=	K;
		//this->nLevels			=	nLevels;
		this->nChebNodes			=	nChebNodes;
		this->rank					=	nChebNodes*nChebNodes;
		this->L						=	L;
	}

	std::vector<pts2D> shift_Cheb_Nodes(double xShift, double yShift) {
		std::vector<pts2D> shiftedChebNodes;
		for (int k=0; k<rank; ++k) {
			pts2D temp;
			temp.x	=	standardChebNodes[k].x+2*xShift;
			temp.y	=	standardChebNodes[k].y+2*yShift;
			shiftedChebNodes.push_back(temp);
		}
		return shiftedChebNodes;
	}
	
	std::vector<pts2D> shift_Leaf_Cheb_Nodes(double xShift, double yShift) {
		std::vector<pts2D> shiftedChebNodes;
		for (int k=0; k<rank; ++k) {
			pts2D temp;
			temp.x	=	leafChebNodes[k].x+xShift;
			temp.y	=	leafChebNodes[k].y+yShift;
			shiftedChebNodes.push_back(temp);
		}
		return shiftedChebNodes;
	}


	//	shifted_scaled_cheb_nodes	//	used in evaluating multipoles
	std::vector<pts2D> shift_scale_Cheb_Nodes(double xShift, double yShift, double radius) {
		std::vector<pts2D> shifted_scaled_ChebNodes;
		for (int k=0; k<rank; ++k) {
			pts2D temp;
			temp.x	=	radius*standardChebNodes[k].x+xShift;
			temp.y	=	radius*standardChebNodes[k].y+yShift;
			shifted_scaled_ChebNodes.push_back(temp);
		}
		return shifted_scaled_ChebNodes;
	}


	//	get_ChebPoly
	double get_ChebPoly(double x, int n) {
		return cos(n*acos(x));
	}

	//	get_S
	double get_S(double x, double y, int n) {
		double S	=	0.5;
		for (int k=1; k<n; ++k) {
			S+=get_ChebPoly(x,k)*get_ChebPoly(y,k);
		}
		return 2.0/n*S;
	}
	//	set_Standard_Cheb_Nodes
	void set_Standard_Cheb_Nodes() {
		for (int k=0; k<nChebNodes; ++k) {
			standardChebNodes1D.push_back(-cos((k+0.5)/nChebNodes*PI));
		}
		pts2D temp1;
		for (int j=0; j<nChebNodes; ++j) {
			for (int k=0; k<nChebNodes; ++k) {
				temp1.x	=	standardChebNodes1D[k];
				temp1.y	=	standardChebNodes1D[j];
				standardChebNodes.push_back(temp1);
			}
		}
		//	Left Bottom child, i.e., Child 0
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x-0.5;
				temp1.y	=	0.5*temp1.y-0.5;
				standardChebNodesChild.push_back(temp1);
		}
		//	Right Bottom child, i.e., Child 1
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x+0.5;
				temp1.y	=	0.5*temp1.y-0.5;
				standardChebNodesChild.push_back(temp1);
		}
		//	Right Top child, i.e., Child 2
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x+0.5;
				temp1.y	=	0.5*temp1.y+0.5;
				standardChebNodesChild.push_back(temp1);
		}
		//	Left Top child, i.e., Child 3
		for (int j=0; j<rank; ++j) {
				temp1	=	standardChebNodes[j];
				temp1.x	=	0.5*temp1.x-0.5;
				temp1.y	=	0.5*temp1.y+0.5;
				standardChebNodesChild.push_back(temp1);
		}
	}

	void get_Transfer_Matrix() {
		for (int l=0; l<4; ++l) {
			L2L[l]	=	Eigen::MatrixXd(rank,rank);
			for (int j=0; j<rank; ++j) {
				for (int k=0; k<rank; ++k) {
					L2L[l](j,k)	=	get_S(standardChebNodes[k].x, standardChebNodesChild[j+l*rank].x, nChebNodes)*get_S(standardChebNodes[k].y, standardChebNodesChild[j+l*rank].y, nChebNodes);
				}
			}
		}
		for (int l=0; l<4; ++l) {
			M2M[l]	=	L2L[l].transpose();
		}
	}



	void get_interpolant_matrix() {
		S_mat	=	Eigen::MatrixXd(4*rank,rank);
		for (int h=0; h<rank; ++h) {
			for (int g=0; g<rank; ++g) {
				S_mat(h,g) = get_S(standardChebNodes[g].x, (standardChebNodes[h].x-1)/2, nChebNodes) * get_S(standardChebNodes[g].y, (standardChebNodes[h].y-1)/2, nChebNodes);
			}
		}
		for (int h=rank; h<2*rank; ++h) {
			for (int g=0; g<rank; ++g) {
				S_mat(h,g) = get_S(standardChebNodes[g].x, (standardChebNodes[h-rank].x+1)/2, nChebNodes) * get_S(standardChebNodes[g].y, (standardChebNodes[h-rank].y-1)/2, nChebNodes);
			}
		}
		for (int h=2*rank; h<3*rank; ++h) {
			for (int g=0; g<rank; ++g) {
				S_mat(h,g) = get_S(standardChebNodes[g].x, (standardChebNodes[h-2*rank].x+1)/2, nChebNodes) * get_S(standardChebNodes[g].y, (standardChebNodes[h-2*rank].y+1)/2, nChebNodes);
			}
		}
		for (int h=3*rank; h<4*rank; ++h) {
			for (int g=0; g<rank; ++g) {
				S_mat(h,g) = get_S(standardChebNodes[g].x, (standardChebNodes[h-3*rank].x-1)/2, nChebNodes) * get_S(standardChebNodes[g].y, (standardChebNodes[h-3*rank].y+1)/2, nChebNodes);
			}
		}
	}



	void createTree() {
		//	First create root and add to tree
		FMM2DBox root;
		root.exists	=	true;
		root.boxNumber		=	0;
		root.parentNumber	=	-1;
		//not sure if it has children
		/*
		#pragma omp parallel for
		for (int l=0; l<4; ++l) {
			root.childrenNumbers[l]	=	l;
		}
		*/
		#pragma omp parallel for
		for (int l=0; l<8; ++l) {
			root.neighborNumbers[l]	=	-1;
		}
		#pragma omp parallel for
		for (int l=0; l<16; ++l) {
			root.innerNumbers[l]	=	-1;
		}
		#pragma omp parallel for
		for (int l=0; l<24; ++l) {
			root.outerNumbers[l]	=	-1;
		}
		//root.center.x = 0.0;
		//root.center.y = 0.0;
		root.center.x = 0.0;
		root.center.y = 0.0;
		root.chebNodes = shift_scale_Cheb_Nodes(root.center.x, root.center.y, L);
		/*cout << "root: " << endl;
		for (int i=0; i<rank; ++i) {
			cout << root.chebNodes[i].x << "	" << root.chebNodes[i].y << endl;
		}*/
		std::vector<FMM2DBox> rootLevel;
		rootLevel.push_back(root);
		tree.push_back(rootLevel);
		
		nBoxesPerLevel.push_back(1);
		boxRadius.push_back(L);
		boxHomogRadius.push_back(pow(L,K->alpha));
		boxLogHomogRadius.push_back(K->alpha*log(L));

		int j=1;
		while (1) {
			nBoxesPerLevel.push_back(4*nBoxesPerLevel[j-1]);
			boxRadius.push_back(0.5*boxRadius[j-1]);
			boxHomogRadius.push_back(pow(0.5,K->alpha)*boxHomogRadius[j-1]);
			boxLogHomogRadius.push_back(boxLogHomogRadius[j-1]-K->alpha*log(2));

			std::vector<FMM2DBox> level_vb;
			bool stop_refining = true;
			for (int k=0; k<nBoxesPerLevel[j-1]; ++k) { // no. of boxes in parent level
				/////////////////////////////////////////////////////////////////////////////////////////////
				// check if there are atleast MIN_CHARGES_PER_BOX charges inside the box to make it have children
				bool childless;
				FMM2DBox box[4];
				Eigen::VectorXd f_at_child = Eigen::VectorXd(4*rank);
				
				//	check for the charges present in parent
				if (tree[j-1][k].exists) { // if the parent is non-existent it's meaningless to have  children for such a box. so avoid checking
					for (int l=0; l<4; ++l) {
						if (l==0) {
							box[l].center.x = tree[j-1][k].center.x - boxRadius[j];
							box[l].center.y = tree[j-1][k].center.y - boxRadius[j];
						}
						else if (l==1) {
							box[l].center.x = tree[j-1][k].center.x + boxRadius[j];
							box[l].center.y = tree[j-1][k].center.y - boxRadius[j];
						}
						else if (l==2) {
							box[l].center.x = tree[j-1][k].center.x + boxRadius[j];
							box[l].center.y = tree[j-1][k].center.y + boxRadius[j];
						}
						else {
							box[l].center.x = tree[j-1][k].center.x - boxRadius[j];
							box[l].center.y = tree[j-1][k].center.y + boxRadius[j];
						}
						box[l].chebNodes = shift_scale_Cheb_Nodes(box[l].center.x, box[l].center.y, boxRadius[j]);
						for (int h=0; h<rank; ++h) {
							f_at_child(l*rank +h) = K->f(box[l].chebNodes[h].x, box[l].chebNodes[h].y);
						}
					}					
					Eigen::VectorXd f_at_parent = Eigen::VectorXd(rank);
					for (int i=0; i<rank; ++i) {
						f_at_parent[i] = K->f(tree[j-1][k].chebNodes[i].x, tree[j-1][k].chebNodes[i].y);
					}
					Eigen::VectorXd P_f = S_mat*f_at_parent;
					Eigen::VectorXd Error = Eigen::VectorXd(4*rank);
					for (int h=0; h<4*rank; ++h) {
						Error(h)	=	fabs((P_f - f_at_child)(h)/f_at_child(h));
					}
					//cout << Error.maxCoeff() << endl;
					if (Error.maxCoeff() <= epsilon_f) 
						childless = true;
					else
						childless = false;
				}
				else
					childless	=	true;
				// if a box can have children it has to be 4
				
				for (int l=0; l<4; ++l) {
					// writing childrenNumbers, boxNumber, parentNumber irrespective of whether they exist or not because it helps later in forming the 4 lists					
					tree[j-1][k].childrenNumbers[l]	=	4*k+l; 
				 	box[l].boxNumber		=	4*k+l;
					box[l].parentNumber	=	k;

					if(childless == false) { // can have 4 children. and these
						box[l].exists		=	true;
						stop_refining = false;
					}	
					else {
						box[l].exists		=	false;
						if (l == 0) { // writing into childless boxes only once
							level_box lb;
							lb.level = j-1;
							lb.box = k;
							if (tree[j-1][k].exists) {
								childless_boxes.push_back(lb);
							}
						}
						// meaningless to have parentNumber to a non-existing child
					}
					level_vb.push_back(box[l]); // putting non existent boxes also into tree
				}
			}
			if (stop_refining) {
				break;
			}
			else {
				tree.push_back(level_vb);
				++j;
			}
			
		}
		nLevels = j-1;
		cout << "nLevels: " << nLevels << endl;
		smallestBoxSize	=	boxRadius[nLevels];
		a		=	smallestBoxSize;
		N		=	rank*childless_boxes.size();
		
		std::vector<FMM2DBox> level_vb;
		for (int k=0; k<4*nBoxesPerLevel[nLevels]; ++k) {
			FMM2DBox box;
			box.exists = false;
			level_vb.push_back(box);
		}
		tree.push_back(level_vb);
		
	}


	void check5() {
		/*for (int i=0; i<childless_boxes.size(); ++i) {
			cout << "level: " << childless_boxes[i].level << "	box: " << childless_boxes[i].box << endl;
		}*/
		for (int i=0; i<tree[5][192].chebNodes.size(); ++i) {
			cout << tree[5][192].chebNodes[i].x << "	" << tree[5][192].chebNodes[i].y << endl;
		}
	}

	//	Assigns the interactions for child0 of a box
	void assign_Child0_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k;
		int nN, nNC;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  N5  |  N4  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  **  |  N3  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[3]	=	nC+1;
			tree[nL][nC].neighborNumbers[4]	=	nC+2;
			tree[nL][nC].neighborNumbers[5]	=	nC+3;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/****************************/
			/*				   ______	*/
			/*				  |		 |	*/
			/*				  |	 **  |	*/
			/*	 _____________|______|  */
			/*	|	   |	  |			*/
			/*	|  I15 |  N0  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I0  |  I1  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[0];
			nNC	=	4*nN;
			if (nN != -1 ) {
				tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[2];
				if (tree[j+1][nNC].exists) {
					tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[0];
					tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[1];
					tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[3];				
				}
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 ______			  		*/
			/*	|	   |	  			*/
			/*	|  **  |				*/
			/*	|______|______			*/
			/*	|	   |	  |			*/
			/*	|  N1  |  N2  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I2  |  I3  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			nNC	=	4*nN;
			if (nN != -1 ) {
				tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[3];
				if (tree[j+1][nNC].exists) {
					tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[0];
					tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[1];
				}
			}
		}

		//	Assign children of parent's second neighbor
		{
			/************************************/
			/*	 ______			  				*/
			/*	|	   |	  					*/
			/*	|  **  |						*/
			/*	|______|	   _____________	*/
			/*				  |	     |	    |	*/
			/*				  |  I5  |  O8  |	*/
			/*				  |______|______|	*/
			/*				  |	     |	    |	*/
			/*				  |  I4  |  O7  |	*/
			/*				  |______|______|	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[2];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[7]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[8]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's third neighbor
		{
			/************************************/
			/*				   _____________	*/
			/*				  |	     |	    |	*/
			/*				  |  I7  |  O10 |	*/
			/*	 ______		  |______|______|	*/
			/*	|	   |	  |	     |	    |	*/
			/*	|  **  |	  |  I6  |  O9  |	*/
			/*	|______|	  |______|______|	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[3];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[9]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[10]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[3];
			}

		}

		//	Assign children of parent's fourth neighbor
		{
			/************************************/
			/*				   _____________	*/
			/*				  |	     |	    |	*/
			/*				  |  O13 |  O12 |	*/
			/*				  |______|______|	*/
			/*				  |	     |	    |	*/
			/*		    	  |  I8  |  O11 |	*/
			/*		    	  |______|______|	*/
			/*									*/
			/*									*/
			/*	 ______							*/
			/*  |      |						*/
			/*  |  **  |						*/
			/*  |______|						*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[4];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[11]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[12]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[13]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  O15 |  O14 |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  I10 |  I9  |		*/
			/*	|______|______|		*/
			/*						*/
			/*						*/
			/*	 ______				*/
			/*  |	   |			*/
			/*	|  **  |			*/
			/*	|______|			*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[14]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[15]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/****************************/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  O17 |  O16 |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I12 |  I11 |			*/
			/*	|______|______|			*/
			/*							*/
			/*							*/
			/*				   ______	*/
			/*  			  |		 |	*/
			/*				  |	 **  |	*/
			/*				  |______|	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[6];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[16]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[17]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/****************************/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I13 |  N6  |			*/
			/*	|______|______|______	*/
			/*  |	   |	  |		 |	*/
			/*	|  I14 |  N7  |	 **  |	*/
			/*	|______|______|______|	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[7];
			nNC	=	4*nN;
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[2];
				if (tree[j+1][nNC].exists) {
					tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[0];
					tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[3];				
				}
			}
		}
	}

	//	Assigns the interactions for child1 of a box
	void assign_Child1_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k+1;
		int nN,nNC;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  N6  |  N5  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N7  |  **  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[7]	=	nC-1;
			tree[nL][nC].neighborNumbers[5]	=	nC+1;
			tree[nL][nC].neighborNumbers[6]	=	nC+2;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/************************************/
			/*				   		  ______	*/
			/*				  	     |		|	*/
			/*				         |	**  |	*/
			/*	 _____________       |______|  	*/
			/*	|	   |	  |					*/
			/*	|  O22 |  I15 |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O23 |  I0  |					*/
			/*	|______|______|					*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[0];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].outerNumbers[23]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[22]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 		______		  	*/
			/*		   |	  |			*/
			/*	       |  **  |			*/
			/*	 ______|______|			*/
			/*	|	   |	  |			*/
			/*	|  N0  |  N1  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I1  |  I2  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			nNC	=	4*nN;
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[3];
				if (tree[j+1][nNC].exists) {
					tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[0];
					tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[1];
				}								
			}
		}

		//	Assign children of parent's second neighbor
		{
			/****************************/
			/*	 ______		  			*/
			/*	|	   |				*/
			/*	|  **  |	  			*/
			/*	|______|_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  N2  |  I5  |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*	       |  I3  |  I4  |	*/
			/*		   |______|______|	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[2];
			nNC	=	4*nN;
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[3];				
				if (tree[j+1][nNC].exists) {
					tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[0];
					tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[1];
					tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[2];
				}
			}
		}

		//	Assign children of parent's third neighbor
		{
			/****************************/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  N4  |	 I7	 |  */
			/*	 ______|______|______|	*/
			/*	|	   |	  |	     |	*/
			/*	|  **  |  N3  |  I6  |  */
			/*	|______|______|______| 	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[3];
			nNC	=	4*nN;
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[3];	
				if(tree[j+1][nNC].exists) {
					tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[1];
					tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[2];
				}							
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			/****************************/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  O14 |  O13 |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*		   |  I9  |  I8  |	*/
			/*		   |______|______|	*/
			/*				  			*/
			/*				  			*/
			/*	 ______					*/
			/*	|	   |				*/
			/*  |  **  |				*/
			/*  |______|				*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[4];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[13]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[14]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  O16 |  O15 |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  I11 |  I10 |		*/
			/*	|______|______|		*/
			/*						*/
			/*						*/
			/*		    ______		*/
			/* 		   |	  |		*/
			/*		   |  **  |		*/
			/*		   |______|		*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[15]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[16]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/************************************/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O18 |  O17 |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O19 |  I12 |					*/
			/*	|______|______|					*/
			/*									*/
			/*									*/
			/*				   		  ______	*/
			/*  			  		 |		|	*/
			/*				  		 |	** 	|	*/
			/*				  		 |______|	*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[6];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].outerNumbers[19]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[17]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[18]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/************************************/
			/*									*/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O20 |  I13 |					*/
			/*	|______|______|		  ______	*/
			/*  |	   |	  |		 |		|	*/
			/*	|  O21 |  I14 |	 	 |	**  |	*/
			/*	|______|______|		 |______|	*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[7];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].outerNumbers[21]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[20]	=	tree[j][nN].childrenNumbers[3];				
			}
		}
	}

	//	Assigns the interactions for child2 of a box
	void assign_Child2_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k+2;
		int nN,nNC;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  N7  |  **  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N0  |  N1  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[0]	=	nC-2;
			tree[nL][nC].neighborNumbers[1]	=	nC-1;
			tree[nL][nC].neighborNumbers[7]	=	nC+1;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/************************************/
			/*				   		  ______	*/
			/*				  	     |		|	*/
			/*				         |	**  |	*/
			/*				         |______|  	*/
			/*									*/
			/*									*/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O23 |  I0  |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O0  |  O1  |					*/
			/*	|______|______|					*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[0];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].outerNumbers[0]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[1]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[23]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 		______		  	*/
			/*		   |	  |			*/
			/*	       |  **  |			*/
			/*	 	   |______|			*/
			/*							*/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I1  |  I2  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  O2  |  O3  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].outerNumbers[2]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[3]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's second neighbor
		{
			/****************************/
			/*	 ______		  			*/
			/*	|	   |				*/
			/*	|  **  |	  			*/
			/*	|______|				*/
			/*							*/
			/*							*/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  I3  |  I4  |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*	       |  O4  |  O5  |	*/
			/*		   |______|______|	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[2];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].outerNumbers[4]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[5]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's third neighbor
		{
			/****************************/
			/*	 ____________________	*/
			/*	|	   |	  |	     |	*/
			/*	|  **  |  N3  |	 I6	 |  */
			/*	|______|______|______|	*/
			/*		   |	  |	     |	*/
			/*		   |  N2  |  I5  |  */
			/*		   |______|______| 	*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[3];
			nNC	=	4*nN;
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[3];		
				if (tree[j+1][nNC].exists) {
					tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[1];
					tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[2];
				}
						
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			/****************************/
			/*			_____________	*/
			/*		   |	  |	     |	*/
			/*		   |  I9  |  I8  |	*/
			/*		   |______|______|	*/
			/*		   |	  |	     |	*/
			/*		   |  N4  |  I7  |	*/
			/*	 ______|______|______|	*/
			/*	|	   |	  			*/
			/*	|  **  |	  			*/
			/*	|______|	  			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[4];
			nNC	=	4*nN;
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[0];
				if (tree[j+1][nNC].exists) {
					tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[1];
					tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[2];
					tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[3];				
				}
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  I11 |  I10 |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N6  |  N5  |		*/
			/*	|______|______|		*/
			/*		   |	  |		*/
			/*		   |  **  |		*/
			/*		   |______|		*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			nNC	=	4*nN;
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[1];
				if (tree[j+1][nNC].exists) {
					tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[2];
					tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[3];				
				}
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/************************************/
			/*	 _____________					*/
			/*	|	   |	  |					*/
			/*	|  O19 |  I12 |					*/
			/*	|______|______|					*/
			/*	|	   |	  |					*/
			/*	|  O20 |  I13 |					*/
			/*	|______|______|		  ______	*/
			/*  			  		 |		|	*/
			/*				  		 |	** 	|	*/
			/*				  		 |______|	*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[6];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].outerNumbers[20]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[19]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/************************************/
			/*									*/
			/*	 _____________		  ______	*/
			/*	|	   |	  |		 |	    |	*/
			/*	|  O21 |  I14 |		 |	**	|	*/
			/*	|______|______|		 |______|	*/
			/*  |	   |	  |		 			*/
			/*	|  O22 |  I15 |	 	 			*/
			/*	|______|______|		 			*/
			/*  								*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[7];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].outerNumbers[22]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].outerNumbers[21]	=	tree[j][nN].childrenNumbers[3];				
			}
		}
	}

	//	Assigns the interactions for child3 of a box
	void assign_Child3_Interaction(int j, int k) {
		int nL	=	j+1;
		int nC	=	4*k+3;
		int nN,nNC;

		//	Assign siblings
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  **  |  N3  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N1  |  N2  |		*/
			/*	|______|______|		*/
			/*						*/
			/************************/
			tree[nL][nC].neighborNumbers[1]	=	nC-3;
			tree[nL][nC].neighborNumbers[2]	=	nC-2;
			tree[nL][nC].neighborNumbers[3]	=	nC-1;
		}

		//	Assign children of parent's zeroth neighbor
		{
			/****************************/
			/*				   ______	*/
			/*				  |		 |	*/
			/*				  |	 **  |	*/
			/*				  |______|  */
			/*							*/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I0  |  I1  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  O1  |  O2  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[0];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].outerNumbers[1]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[2]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's first neighbor
		{
			/****************************/
			/*	 ______		  			*/
			/*	|	   |				*/
			/*	|  **  |				*/
			/*	|______|				*/
			/*							*/
			/*							*/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I2  |  I3  |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  O3  |  O4  |			*/
			/*	|______|______|			*/
			/*							*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[1];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].outerNumbers[3]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[4]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's second neighbor
		{
			/************************************/
			/*	 ______		  					*/
			/*	|	   |						*/
			/*	|  **  |	  					*/
			/*	|______|						*/
			/*									*/
			/*									*/
			/*				   _____________	*/
			/*		   		  |	     |	    |	*/
			/*		   		  |  I4  |  O7  |	*/
			/*		   		  |______|______|	*/
			/*		   		  |	  	 |	    |	*/
			/*	       		  |  O5  |  O6  |	*/
			/*		   		  |______|______|	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[2];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].outerNumbers[5]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[6]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[7]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's third neighbor
		{
			/************************************/
			/*	 ______		   _____________	*/
			/*	|	   |	  |	     |		|	*/
			/*	|  **  |      |	 I6	 |  O9	|	*/
			/*	|______|	  |______|______|	*/
			/*		   		  |	  	 |	    |	*/
			/*		   		  |  I5  |  O8  |  	*/
			/*		   		  |______|______| 	*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[3];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[8]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[9]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			/************************************/
			/*				   _____________	*/
			/*		   		  |	  	 |	    |	*/
			/*		   		  |  I8  |  O11 |	*/
			/*		   		  |______|______|	*/
			/*		   		  |	  	 |	    |	*/
			/*		   		  |  I7  |  O10 |	*/
			/*	 ______	      |______|______|	*/
			/*	|	   |	  					*/
			/*	|  **  |	  					*/
			/*	|______|	  					*/
			/*									*/
			/************************************/
			nN	=	tree[j][k].neighborNumbers[4];
			nNC	=	4*nN;
			if (nN != -1 && tree[j+1][nNC].exists) {
				tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].outerNumbers[10]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].outerNumbers[11]	=	tree[j][nN].childrenNumbers[2];
				tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[3];				
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			/************************/
			/*	 _____________		*/
			/*	|	   |	  |		*/
			/*	|  I10 |  I9  |		*/
			/*	|______|______|		*/
			/*	|	   |	  |		*/
			/*	|  N5  |  N4  |		*/
			/*	|______|______|		*/
			/*	|	   |			*/
			/*	|  **  |			*/
			/*	|______|			*/
			/*  					*/
			/************************/
			nN	=	tree[j][k].neighborNumbers[5];
			nNC	=	4*nN;
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[0];
				tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[1];
				if (tree[j+1][nNC].exists) {
					tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[2];
					tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[3];				
				}
			}
		}

		//	Assign children of parent's sixth neighbor
		{
			/****************************/
			/*	 _____________			*/
			/*	|	   |	  |			*/
			/*	|  I12 |  I11 |			*/
			/*	|______|______|			*/
			/*	|	   |	  |			*/
			/*	|  I13 |  N6  |			*/
			/*	|______|______|______	*/
			/*  			  |		 |	*/
			/*				  |	 **  |	*/
			/*				  |______|	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[6];
			nNC	=	4*nN;
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[1];
				if (tree[j+1][nNC].exists) {
					tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[0];
					tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[2];
					tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[3];				
				}
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			/****************************/
			/*							*/
			/*	 ____________________	*/
			/*	|	   |	  |		 |	*/
			/*	|  I14 |  N7  |	 **	 |	*/
			/*	|______|______|______|	*/
			/*  |	   |	  |		 	*/
			/*	|  I15 |  N0  |	 	 	*/
			/*	|______|______|		 	*/
			/*  						*/
			/****************************/
			nN	=	tree[j][k].neighborNumbers[7];
			nNC	=	4*nN;
			if (nN != -1) {
				tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[1];
				tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[2];
				if (tree[j+1][nNC].exists) {
					tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[0];
					tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[3];	
				}			
			}
		}
	}

	//	Assigns the interactions for the children of a box
	void assign_Box_Interactions(int j, int k) {
		assign_Child0_Interaction(j,k);
		assign_Child1_Interaction(j,k);
		assign_Child2_Interaction(j,k);
		assign_Child3_Interaction(j,k);
	}

	//	Assigns the interactions for the children all boxes at a given level
	void assign_Level_Interactions(int j) {
		//cout << endl << "tree[j].size: " << tree[j].size() << endl;
		//#pragma omp parallel for
		for (int k=0; k<tree[j].size(); ++k) {
			//cout << k << endl;
			if (tree[j+1][4*k].exists) { // do this only if the box is a parent
				//cout << j+1 << "," << 4*k << "level box" << endl;
				assign_Box_Interactions(j,tree[j][k].boxNumber);
			}
		}
	}

	//	Assigns the interactions for the children of all boxes in the tree
	//	Assigns colleagues(neighbors) and list2 (inner and outer numbers) needed for M2L(same size)
	void assign_list2_Neighbor_Interactions() {
		for (int j=0; j<nLevels; ++j) {
			assign_Level_Interactions(j);
		}
	}


	//	Assigns list1 for childless boxes	
	void assign_list1_list3_Interactions() {
		for (int j=0; j<childless_boxes.size(); ++j) {
			//childless_boxes[j].level = l;
			//childless_boxes[j].box = b;
			assign_list1_list3_box_Interactions(childless_boxes[j]);
		}
	}





	void check2 () {
		cout << endl << "tree[5][512].exists: " << tree[5][512].exists << endl;
		cout << endl << "tree[4][128].childrenNumbers[0]: " << tree[4][128].childrenNumbers[0] << endl;
	}


	void check () {
		
		for (int j=0; j<childless_boxes.size(); ++j) {
			cout << "j: " << childless_boxes[j].level << "	k: " << childless_boxes[j].box << endl;
		}
	}


	//	LIST 1
	void assign_list1_list3_box_Interactions(const level_box lb) { // for boxes which donot have children
		level_box temp;
		
		//	if box is childless it also is a member of list 1
		//tree[lb.level][lb.box].list1.push_back(lb);
		level_box prev_add, neigh0_add;
		prev_add.level = -1;
		prev_add.box = 0;
		neigh0_add.level = -1;
		neigh0_add.box = 0;
		int j,k,cN;
		
		//	NEIGHBOR 0,2,4,6 CORNERS
		//	NEIGHBOR 1,3,5,7 SIDE BY SIDE

		//	NEIGHBOR 0
		
		k = tree[lb.level][lb.box].neighborNumbers[0]; // neighbor is in same level
		j = lb.level;
		if (k != -1) { // k=-1 means that box is not in the computational domain
			cN = k;
			if (tree[j][k].exists) {
				if (tree[j+1][4*cN].exists) {	
					while (tree[j+1][4*cN].exists) {
					
					//	list 3
						temp.level = j+1;
						temp.box = 4*cN;
						tree[lb.level][lb.box].list3.push_back(temp);
						tree[temp.level][temp.box].list4.push_back(lb);	

						temp.level = j+1;
						temp.box = 4*cN+1;
						tree[lb.level][lb.box].list3.push_back(temp);
						tree[temp.level][temp.box].list4.push_back(lb);	

						temp.level = j+1;
						temp.box = 4*cN+3;
						tree[lb.level][lb.box].list3.push_back(temp);
						tree[temp.level][temp.box].list4.push_back(lb);	
					
						j = j+1;
						cN = 4*cN+2;
					}
				
					temp.level = j;
					temp.box = cN;
					tree[lb.level][lb.box].list1.push_back(temp);
				}
			}
			else {
				// immendiate next valid parent // goto higher levels
				while (tree[j][k].exists == false) {		
					j = j-1;
					k = k/4;
				}
				temp.level = j;
				temp.box = k;
				tree[lb.level][lb.box].list1.push_back(temp);
				prev_add = temp;
				neigh0_add = temp;
			}
		}

		
		//	NEIGHBOR 1
		k = tree[lb.level][lb.box].neighborNumbers[1]; // neighbor is in same level
		j = lb.level;
		if (k != -1) { 
			if (tree[j][k].exists && tree[j+1][4*k].exists) {
				list1_neighbor1(lb.level,lb.box,j,k);
			}
		}


		
		//	NEIGHBOR 2
		
		k = tree[lb.level][lb.box].neighborNumbers[2]; // neighbor is in same level
		j = lb.level;
		if (k != -1) { 
			cN = k;
			if (tree[j][k].exists) {
				if (tree[j+1][4*cN].exists) {
					while (tree[j+1][4*cN].exists) {
						temp.level = j+1;
						temp.box = 4*cN;
						tree[lb.level][lb.box].list3.push_back(temp);
						tree[temp.level][temp.box].list4.push_back(lb);	
					
						temp.level = j+1;
						temp.box = 4*cN+1;
						tree[lb.level][lb.box].list3.push_back(temp);
						tree[temp.level][temp.box].list4.push_back(lb);	
						
						temp.level = j+1;
						temp.box = 4*cN+2;
						tree[lb.level][lb.box].list3.push_back(temp);
						tree[temp.level][temp.box].list4.push_back(lb);	
							
						j = j+1;
						cN = 4*cN+3;
					}
					temp.level = j;
					temp.box = cN;
					tree[lb.level][lb.box].list1.push_back(temp);
				}
			}
			else {
				// immendiate next valid parent // goto higher levels
				while (tree[j][k].exists == false) {		
					j = j-1;
					k = k/4;
				}
				temp.level = j;
				temp.box = k;
				if (prev_add.level == j && prev_add.box == k) {
				}
				else {
					tree[lb.level][lb.box].list1.push_back(temp);
				}
				prev_add = temp;
			}
		}


		

		//	NEIGHBOR 3
		k = tree[lb.level][lb.box].neighborNumbers[3]; // neighbor is in same level
		j = lb.level;
		if (k != -1) { 
			if (tree[j][k].exists && tree[j+1][4*k].exists) {
				list1_neighbor3(lb.level,lb.box,j,k);
			}
		}



		//	NEIGHBOR 4
		
		k = tree[lb.level][lb.box].neighborNumbers[4]; // neighbor is in same level
		j = lb.level;
		if (k != -1) { 
			cN = k;
			if (tree[j][k].exists) {
				if (tree[j+1][4*cN].exists) {
					while (tree[j+1][4*cN].exists) {
						temp.level = j+1;
						temp.box = 4*cN+1;
						tree[lb.level][lb.box].list3.push_back(temp);
						tree[temp.level][temp.box].list4.push_back(lb);	
					
						temp.level = j+1;
						temp.box = 4*cN+2;
						tree[lb.level][lb.box].list3.push_back(temp);
						tree[temp.level][temp.box].list4.push_back(lb);	
					
						temp.level = j+1;
						temp.box = 4*cN+3;
						tree[lb.level][lb.box].list3.push_back(temp);
						tree[temp.level][temp.box].list4.push_back(lb);	
						j = j+1;
						cN = 4*cN;
					}
					temp.level = j;
					temp.box = cN;
					tree[lb.level][lb.box].list1.push_back(temp);
				}
			}
			else {
				// immendiate next valid parent // goto higher levels
				while (tree[j][k].exists == false) {		
					j = j-1;
					k = k/4;
				}
				temp.level = j;
				temp.box = k;
				if (prev_add.level == j && prev_add.box == k) {
				}
				else {
					tree[lb.level][lb.box].list1.push_back(temp);
				}
				prev_add = temp;
			}
		}


		


		//	NEIGHBOR 5
		k = tree[lb.level][lb.box].neighborNumbers[5]; // neighbor is in same level
		j = lb.level;
		if (k != -1) { 
			if (tree[j][k].exists && tree[j+1][4*k].exists) {
				list1_neighbor5(lb.level,lb.box,j,k);
			}
		}




		//	NEIGHBOR 6
		
		k = tree[lb.level][lb.box].neighborNumbers[6]; // neighbor is in same level
		j = lb.level;
		if (k != -1) { 
			cN = k;
			if (tree[j][k].exists) {
				if (tree[j+1][4*cN].exists) {
					while (tree[j+1][4*cN].exists) {
						temp.level = j+1;
						temp.box = 4*cN;
						tree[lb.level][lb.box].list3.push_back(temp);
						tree[temp.level][temp.box].list4.push_back(lb);	
						
						temp.level = j+1;
						temp.box = 4*cN+2;
						tree[lb.level][lb.box].list3.push_back(temp);
						tree[temp.level][temp.box].list4.push_back(lb);	
							
						temp.level = j+1;
						temp.box = 4*cN+3;
						tree[lb.level][lb.box].list3.push_back(temp);
						tree[temp.level][temp.box].list4.push_back(lb);	
					
						j = j+1;
						cN = 4*cN+1;
					}
					temp.level = j;
					temp.box = cN;
					tree[lb.level][lb.box].list1.push_back(temp);
				}
			}
			else {
				// immendiate next valid parent // goto higher levels
				while (tree[j][k].exists == false) {		
					j = j-1;
					k = k/4;
				}
				temp.level = j;
				temp.box = k;
				if ((prev_add.level == j && prev_add.box == k) || (neigh0_add.level == j && neigh0_add.box == k)) {
				}
				else {
					tree[lb.level][lb.box].list1.push_back(temp);
				}
				prev_add = temp;
			}
		}


		
		
		//	NEIGHBOR 7
		k = tree[lb.level][lb.box].neighborNumbers[7]; // neighbor is in same level
		j = lb.level;
		if (k != -1) { 
			if (tree[j][k].exists && tree[j+1][4*k].exists) {
				list1_neighbor7(lb.level,lb.box,j,k);
			}
		}
	}

	

	void list1_neighbor1(int Lj, int Lk, int Nj, int Nk) {
		level_box temp;
		if (tree[Nj+1][4*Nk].exists) {
			level_box lb;
			lb.level = Lj;
			lb.box = Lk;
			// list 3
			temp.level = Nj+1;
			temp.box = 4*Nk;
			tree[Lj][Lk].list3.push_back(temp);
			tree[temp.level][temp.box].list4.push_back(lb);
		
			temp.box = 4*Nk+1;
			tree[Lj][Lk].list3.push_back(temp);
			tree[temp.level][temp.box].list4.push_back(lb);
		
			list1_neighbor1(Lj,Lk,Nj+1,4*Nk+2);
			list1_neighbor1(Lj,Lk,Nj+1,4*Nk+3);

		}
		else {
			temp.level = Nj;
			temp.box = Nk;
			tree[Lj][Lk].list1.push_back(temp);

		}
	}



	void list1_neighbor3(int Lj, int Lk, int Nj, int Nk) {
		level_box temp;
		if (tree[Nj+1][4*Nk].exists) {
			level_box lb;
			lb.level = Lj;
			lb.box = Lk;
			
			temp.level = Nj+1;
			temp.box = 4*Nk+1;
			tree[Lj][Lk].list3.push_back(temp);
			tree[temp.level][temp.box].list4.push_back(lb);
		
			temp.box = 4*Nk+2;
			tree[Lj][Lk].list3.push_back(temp);
			tree[temp.level][temp.box].list4.push_back(lb);
		
			list1_neighbor3(Lj,Lk,Nj+1,4*Nk+0);
			list1_neighbor3(Lj,Lk,Nj+1,4*Nk+3);

		}
		else {
			temp.level = Nj;
			temp.box = Nk;
			tree[Lj][Lk].list1.push_back(temp);
		}
	}




	void list1_neighbor5(int Lj, int Lk, int Nj, int Nk) {
		level_box temp;
		if (tree[Nj+1][4*Nk].exists) {
			level_box lb;
			lb.level = Lj;
			lb.box = Lk;
			
			temp.level = Nj+1;
			temp.box = 4*Nk+2;
			tree[Lj][Lk].list3.push_back(temp);
			tree[temp.level][temp.box].list4.push_back(lb);
				
			temp.box = 4*Nk+3;
			tree[Lj][Lk].list3.push_back(temp);
			tree[temp.level][temp.box].list4.push_back(lb);
		
			list1_neighbor5(Lj,Lk,Nj+1,4*Nk);
			list1_neighbor5(Lj,Lk,Nj+1,4*Nk+1);

					
		}
		else {
			temp.level = Nj;
			temp.box = Nk;
			tree[Lj][Lk].list1.push_back(temp);

		}
	}




	void list1_neighbor7(int Lj, int Lk, int Nj, int Nk) {
		level_box temp;
		if (tree[Nj+1][4*Nk].exists) {
			level_box lb;
			lb.level = Lj;
			lb.box = Lk;
			
			temp.level = Nj+1;
			temp.box = 4*Nk;
			tree[Lj][Lk].list3.push_back(temp);
			tree[temp.level][temp.box].list4.push_back(lb);
				
			temp.box = 4*Nk+3;
			tree[Lj][Lk].list3.push_back(temp);
			tree[temp.level][temp.box].list4.push_back(lb);
		
			list1_neighbor7(Lj,Lk,Nj+1,4*Nk+1);
			list1_neighbor7(Lj,Lk,Nj+1,4*Nk+2);
		}
		else {
			temp.level = Nj;
			temp.box = Nk;
			tree[Lj][Lk].list1.push_back(temp);
		}
	}




	void check3 () {
		cout << endl;
		for (int j=2; j<=nLevels; ++j) { // list 4 doesn't exist for boxes in level < 2
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				if (tree[j][k].exists) { // list 4 is found only for existing boxes
					cout << endl <<  "j: " << j << "	k: " << k << endl;
					cout << "list4: " << endl;
					for (int i=0; i<tree[j][k].list4.size() ; ++i) {
						cout << "j: " << tree[j][k].list4[i].level << "	k: " << tree[j][k].list4[i].box << endl;
					}
				}	
			}
		}
	}



	//	Obtain the desired matrix for M2L (uniform GQ)
	void obtain_Desired_Operator_M2L(double xcenter_j, double ycenter_j, Eigen::MatrixXd& T) {
		T	=	Eigen::MatrixXd(rank,rank);
		double L_j = 1.0;
		int nNodes		=	6;	
		for (int i=0; i<rank; ++i) {
			for (int j=0; j<rank; ++j) {
				user_Greens_function* mine	=	new user_Greens_function();
				mine->initialise(standardChebNodes[i], standardChebNodes[j], nChebNodes, xcenter_j, ycenter_j, L_j, nNodes);
				mine->obtain_Quadrature_Uniform();
				T(i,j)	=	double(mine->integral);
				delete mine;
			}
		}
	}


	//	Obtain the desired matrix for neighbors and self(AGQ)
	void obtain_Desired_Operator_Leaf(double xcenter_j, double ycenter_j, Eigen::MatrixXd& T) {
		T	=	Eigen::MatrixXd(rank,rank);
		double L_j = 1.0;
		int nNodes		=	6;	
		for (int i=0; i<rank; ++i) {
			for (int j=0; j<rank; ++j) {
				user_Greens_function* mine	=	new user_Greens_function();
				mine->initialise(standardChebNodes[i], standardChebNodes[j], nChebNodes, xcenter_j, ycenter_j, L_j, nNodes);
				mine->obtain_Quadrature(epsilon);
				T(i,j)	=	double(mine->integral);
				delete mine;
			}
		}
	}


	void obtain_Desired_Operator2(Eigen::RowVectorXd& T) {
		T	=	Eigen::RowVectorXd(rank);
		int nNodes		=	6;	
		for (int j=0; j<rank; ++j) {
			user_function2* mine	=	new user_function2();
			mine->initialise(standardChebNodes[j], nChebNodes, nNodes);
			mine->obtain_Quadrature_Uniform();
			T(j)	=	double(mine->integral);
			delete mine;
		}
	}


	//	Assemble FMM Operators
	void assemble_Operators_FMM() {
		//	Assemble Outer Interactions
		for (int l=0; l<6; ++l) {
			obtain_Desired_Operator_M2L(2.0*(l-3), -6.0, M2LOuter[l]);
		}
		for (int l=0; l<6; ++l) {
			obtain_Desired_Operator_M2L(6.0, 2.0*(l-3), M2LOuter[l+6]);
		}
		for (int l=0; l<6; ++l) {
			obtain_Desired_Operator_M2L(2.0*(3-l), 6.0, M2LOuter[l+12]);
		}
		for (int l=0; l<6; ++l) {
			obtain_Desired_Operator_M2L(-6.0, 2.0*(3-l), M2LOuter[l+18]);
		}
		//	Assemble Inner Interactions
		for (int l=0; l<4; ++l) {
			obtain_Desired_Operator_M2L(2.0*(l-2), -4.0, M2LInner[l]);
		}
		for (int l=0; l<4; ++l) {
			obtain_Desired_Operator_M2L(4.0, 2.0*(l-2), M2LInner[l+4]);
		}
		for (int l=0; l<4; ++l) {
			obtain_Desired_Operator_M2L(2.0*(2-l), 4.0, M2LInner[l+8]);
		}
		for (int l=0; l<4; ++l) {
			obtain_Desired_Operator_M2L(-4.0, 2.0*(2-l), M2LInner[l+12]);
		}
		obtain_Desired_Operator2(M2L_term2);
								
		//	Assemble Neighbor Interactions
		obtain_Desired_Operator_Leaf(-2.0, -2.0, neighborInteraction[0]);
		obtain_Desired_Operator_Leaf(0.0, -2.0, neighborInteraction[1]);
		obtain_Desired_Operator_Leaf(2.0, -2.0, neighborInteraction[2]);
		obtain_Desired_Operator_Leaf(2.0, 0.0, neighborInteraction[3]);
		obtain_Desired_Operator_Leaf(2.0, 2.0, neighborInteraction[4]);
		obtain_Desired_Operator_Leaf(0.0, 2.0, neighborInteraction[5]);
		obtain_Desired_Operator_Leaf(-2.0, 2.0, neighborInteraction[6]);
		obtain_Desired_Operator_Leaf(-2.0, 0.0, neighborInteraction[7]);
		
		//	Assemble Self Interactions
		obtain_Desired_Operator_Leaf(0.0, 0.0, selfInteraction);
	}



	// evaluate multipoles
	void evaluate_multipoles() {
		for (int j=0; j<childless_boxes.size(); ++j) { // leafs
			level_box lb = childless_boxes[j];
			tree[lb.level][lb.box].multipoles = Eigen::VectorXd(rank);
			for (int m=0; m<rank; ++m) { // multipoles
				tree[lb.level][lb.box].multipoles[m] = K->f(tree[lb.level][lb.box].chebNodes[m].x, tree[lb.level][lb.box].chebNodes[m].y);
			}
		}
	}



	void evaluate_All_M2M() {
		for (int j=nLevels-1; j>1; --j) { // parent
			//int J	=	j+1; // children
			#pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) { // parent
				//int K	=	4*k; // children
				if (tree[j][k].exists) {
					tree[j][k].multipoles = Eigen::VectorXd(rank);
					for (int m=0; m<rank; ++m) { // multipoles
						tree[j][k].multipoles[m] = K->f(tree[j][k].chebNodes[m].x, tree[j][k].chebNodes[m].y);
						//tree[j][k].multipoles	=	M2M[0]*tree[J][K].multipoles+M2M[1]*tree[J][K+1].multipoles+M2M[2]*tree[J][K+2].multipoles+M2M[3]*tree[J][K+3].multipoles;
					}
				}
			}
		}
	}

	void evaluate_list2() { // list 2
		for (int j=0; j<=1; ++j) {
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				if (tree[j][k].exists) {
					tree[j][k].locals	=	Eigen::VectorXd::Zero(rank);
				}
			}
		}
		for (int j=2; j<=nLevels; ++j) {
			#pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				if (tree[j][k].exists) {
					tree[j][k].locals	=	Eigen::VectorXd::Zero(rank);
					#ifdef HOMOG
						//	Inner well-separated clusters
						for (int l=0; l<16; ++l) {
							int nInner	=	tree[j][k].innerNumbers[l];
							if (nInner>-1) {
								tree[j][k].locals+=M2LInner[l]*tree[j][nInner].multipoles;
							}
						}
						//	Outer well-separated clusters
						for (int l=0; l<24; ++l) {
							int nOuter	=	tree[j][k].outerNumbers[l];
							if (nOuter>-1) {
								tree[j][k].locals+=M2LOuter[l]*tree[j][nOuter].multipoles;
							}
						}
						tree[j][k].locals*=boxRadius[j];					
					#elif LOGHOMOG
						//	Inner well-separated clusters
						for (int l=0; l<16; ++l) {
							int nInner	=	tree[j][k].innerNumbers[l];
							if (nInner>-1) {
							//	if (nInner==6||nInner==7) {
								tree[j][k].locals += (1/(2*PI))*M2LInner[l]*tree[j][nInner].multipoles;
								double temp = (1/(2*PI))*boxLogHomogRadius[j]*M2L_term2*tree[j][nInner].multipoles;
								tree[j][k].locals+=temp*Eigen::VectorXd::Ones(rank);
							//	}
							}
						}
						//	Outer well-separated clusters
						for (int l=0; l<24; ++l) {
							int nOuter	=	tree[j][k].outerNumbers[l];
							if (nOuter>-1) {
							//	if (nOuter==5 || nOuter==4) {
								tree[j][k].locals += (1/(2*PI))*M2LOuter[l]*tree[j][nOuter].multipoles;
								double temp = (1/(2*PI))*boxLogHomogRadius[j]*M2L_term2*tree[j][nOuter].multipoles;
								tree[j][k].locals+=temp*Eigen::VectorXd::Ones(rank);
							//}
							}
						}
						tree[j][k].locals*=pow(boxRadius[j],2);
					#endif
				}
			}
		}
	}




	void evaluate_list3() {
		for (int i=0; i<childless_boxes.size(); ++i) { // all childless boxes
			level_box lb = childless_boxes[i]; // box b
			for (int j=0; j<tree[lb.level][lb.box].list3.size(); ++j) { // all list3 boxes of a childless box
				level_box l3 = tree[lb.level][lb.box].list3[j];
				// list 3; level=l3.level; box = l3.box
				Eigen::MatrixXd Kernel_Matrix;
				obtain_List3_Operator(tree[l3.level][l3.box].center.x, tree[l3.level][l3.box].center.y, boxRadius[l3.level], lb,Kernel_Matrix);
				tree[lb.level][lb.box].locals += Kernel_Matrix*tree[l3.level][l3.box].multipoles;
			}
		}
	}



//	Obtain the desired matrix for list 3
	void obtain_List3_Operator(double xcenter_j, double ycenter_j, double L_j, level_box& lb, Eigen::MatrixXd& T) { // lb is childless box
		int nNodes = 6;
		T	=	Eigen::MatrixXd(rank,rank);
		for (int i=0; i<rank; ++i) {
			for (int j=0; j<rank; ++j) {
				user_Greens_function* mine	=	new user_Greens_function();
				mine->initialise(tree[lb.level][lb.box].chebNodes[i], standardChebNodes[j], nChebNodes, xcenter_j, ycenter_j, L_j, nNodes);
				mine->obtain_Quadrature_Uniform();
				T(i,j)	=	double(mine->integral/(2*PI));
				delete mine;
			}
		}
	}



	void evaluate_list4() {
		for (int l=2; l<=nLevels; ++l) { // all boxes // box b
			for (int b=0; b<nBoxesPerLevel[l]; ++b) {
				if (tree[l][b].exists) {
					for (int j=0; j<tree[l][b].list4.size(); ++j) { // all list3 boxes of a childless box
						level_box l4 = tree[l][b].list4[j];
						Eigen::MatrixXd Kernel_Matrix;
						// list 4; level=l4.level; box = l4.box
						level_box lb;
						lb.level = l;
						lb.box = b;
						obtain_List4_Operator_2(tree[l4.level][l4.box].center.x, tree[l4.level][l4.box].center.y, boxRadius[l4.level], lb,Kernel_Matrix);
						tree[l][b].locals += Kernel_Matrix*tree[l4.level][l4.box].multipoles;
					}
				}
			}
		}
	}




//	Obtain the desired matrix for list 4 of childless boxes
	void obtain_List4_Operator_2(double xcenter_j, double ycenter_j, double L_j, level_box& lb, Eigen::MatrixXd& T) { // lb is childless box
		int nNodes = 6;
		T	=	Eigen::MatrixXd(rank,rank);
		for (int i=0; i<rank; ++i) {
			for (int j=0; j<rank; ++j) {
				user_Greens_function* mine	=	new user_Greens_function();
				mine->initialise(tree[lb.level][lb.box].chebNodes[i], standardChebNodes[j], nChebNodes, xcenter_j, ycenter_j, L_j, nNodes);
				mine->obtain_Quadrature_Uniform();
				T(i,j)	=	double(mine->integral/(2*PI));
				delete mine;
			}
		}
	}




	void evaluate_All_L2L() {
		for (int j=2; j<nLevels; ++j) {
			int J	=	j+1;
			#pragma omp parallel for
			for (int k=0; k<nBoxesPerLevel[j]; ++k) {
				int K	=	4*k;
				if (tree[J][K].exists) { // if children exist
					tree[J][K].locals+=L2L[0]*tree[j][k].locals;
					tree[J][K+1].locals+=L2L[1]*tree[j][k].locals;
					tree[J][K+2].locals+=L2L[2]*tree[j][k].locals;
					tree[J][K+3].locals+=L2L[3]*tree[j][k].locals;
				}
			}
		}
	}


//	Self & Neighbor Interaction
	void evaluate_list1() {
		int nNodes = 1;
		for (int i=0; i<childless_boxes.size(); ++i) { // all childless boxes
			level_box lb = childless_boxes[i];
			
				for (int j=0; j<tree[lb.level][lb.box].list1.size(); ++j) { // all list3 boxes of a childless box
					level_box l1 = tree[lb.level][lb.box].list1[j];
					cout << "j: " << lb.level << "	k:" << lb.box << "	l1j: " << l1.level << "	l1k: " << l1.box << endl;
					Eigen::MatrixXd V_mat	=	Eigen::MatrixXd(rank,rank);
					//	obtain_V_mat
					for (int c=0; c<rank; ++c) { //i
						for (int d=0; d<rank; ++d) { //j
							user_Greens_function* mine	=	new user_Greens_function();
							mine->initialise(tree[lb.level][lb.box].chebNodes[c], standardChebNodes[d], nChebNodes, tree[l1.level][l1.box].center.x, tree[l1.level][l1.box].center.y, boxRadius[l1.level], nNodes);
							mine->obtain_Quadrature(epsilon);
							V_mat(c,d)	=	double(mine->integral/(2*PI));
							delete mine;
						}
					}
					tree[lb.level][lb.box].locals	+=	V_mat * tree[l1.level][l1.box].multipoles;
				}
		}
	}




	void evaluate_Leaf() {
		//#pragma omp parallel for
		for (int i=0; i<childless_boxes.size(); ++i) {
			int j = childless_boxes[i].level;
			int k = childless_boxes[i].box;
			#ifdef HOMOG
			for (int l=0; l<8; ++l) {
				int nneighbor	=	tree[j][k].neighborNumbers[l];
				if (nneighbor != -1) {
					if (tree[j][nneighbor].exists && !tree[j+1][4*nneighbor].exists) {
						tree[j][k].locals+=boxRadius[j]*neighborInteraction[l]*tree[j][nneighbor].multipoles;
					}
				}
			}
			#elif LOGHOMOG
			//tree[j][k].locals	=	Eigen::VectorXd::Zero(rank);
			for (int l=0; l<8; ++l) {
				int nneighbor	=	tree[j][k].neighborNumbers[l];
				if (nneighbor != -1) {
					if (tree[j][nneighbor].exists && !tree[j+1][4*nneighbor].exists) {
						tree[j][k].locals += pow(boxRadius[j],2)*(1/(2*PI))*neighborInteraction[l]*tree[j][nneighbor].multipoles;
						double temp = pow(boxRadius[j],2)*(1/(2*PI))*boxLogHomogRadius[j]*M2L_term2*tree[j][nneighbor].multipoles;	
						tree[j][k].locals+=temp*Eigen::VectorXd::Ones(rank);
					}
				}
			}
			
			//	Self Interaction	
			tree[j][k].locals += pow(boxRadius[j],2)*(1/(2*PI))*selfInteraction*tree[j][k].multipoles;
			double temp = pow(boxRadius[j],2)*(1/(2*PI))*boxLogHomogRadius[j]*M2L_term2*tree[j][k].multipoles;	
			tree[j][k].locals += temp*Eigen::VectorXd::Ones(rank);
			#endif			

		}
	}




	void perform_Error_Check() {
		srand (time(NULL));
		int c 	=	rand()%childless_boxes.size();
		//for (int c=0; c<childless_boxes.size(); ++c) {
		int nNodes = 6;
		level_box nBox = childless_boxes[c];
		
		Eigen::VectorXd potential	=	Eigen::VectorXd::Zero(rank);
		for (int l1=0; l1<rank; ++l1) {// cheb nodes of nBox
			user_function3* mine	=	new user_function3();
			mine->initialise(tree[nBox.level][nBox.box].chebNodes[l1], 0.0, 0.0, L, nNodes);
		//	mine->initialise(tree[nBox.level][nBox.box].chebNodes[l1], 0.5, -0.5, 0.5, nNodes);
			mine->obtain_Quadrature(epsilon);
			potential(l1) = mine->integral/(2*PI);
		}
		Eigen::VectorXd error(rank);
		for (int k=0; k<rank; ++k) {
			//error(k)	=	fabs(tree[nBox.level][nBox.box].locals(k));
			error(k)	=	fabs((potential-tree[nBox.level][nBox.box].locals)(k)/potential(k));
		}
		cout << "nBox.level: " << nBox.level << " nBox.box: " << nBox.box << " er: "<< error.maxCoeff() << endl;
		//cout << "nBox.box: " << nBox.box << endl;
		//cout << "M2LOuter[21]: " << endl << M2LOuter[21] << endl;
		//cout << "M2LOuter[3]: " << endl << M2LOuter[3] << endl;
		//cout <<endl << error << endl;
		//}
		
		//return error.maxCoeff();
	}


	void check4() {
		
		cout << "nn" << endl;
		for (int l=0; l<8; ++l) {
			cout << tree[2][0].neighborNumbers[l]<<",";
		}
		cout << endl << "in" << endl;
		for (int l=0; l<16; ++l) {
			cout << tree[2][0].innerNumbers[l]<<",";
		}
		cout << endl << "on" << endl;
		for (int l=0; l<16; ++l) {
			cout << tree[2][0].outerNumbers[l]<<",";
		}
		cout << endl;
	}

};

#endif
