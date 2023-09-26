/*C++ CODE - MANGEAT MATTHIEU - 2022*/
/*TWO SPECIES VICSEK MODEL*/

//////////////////////
///// LIBRAIRIES /////
//////////////////////

//Public librairies
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <string.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

//Personal libraries
#include "lib/random.cpp"
#include "lib/special_functions.cpp"

//////////////////////////
///// PARTICLE CLASS /////
//////////////////////////

class particle
{
	public:
	
	double x,y; //Position of the particle.
	double dx, dy; //Displacement of the particle.
	double theta; //Orientation of the particle.
	double theta_avg; //Average orientation in the neighborhood.
	int spin; //Spin of the particle (+1 or -1).

	particle(const int &LX, const int &LY, const int &spin0, const int &init); //Initial position and orientation.
	void update(const double &v0, const double &eta, const int &LX, const int &LY); //Update position and orientation.
};

particle::particle(const int &LX, const int &LY, const int &spin0, const int &init)
{
	//Spin of the particle.
	spin=spin0;
	
	//Displacement.
	dx=0.;
	dy=0.;
	
	//Position.
	if (init==0 or init==1)
	{
		if (spin==-1) //Particle B.
		{
			x=LX*ran()/8;
			y=LY*ran();
		}
		else //Particle A.
		{
			x=0.5*LX+LX*ran()/8;
			y=LY*ran();
		}
	}
	else //Random position.
	{
		x=LX*ran();
		y=LY*ran();
	}
	
	//Orientation.
	if (init==0) //APF state.
	{
		if (spin==-1)
		{
			theta=M_PI;
		}
		else
		{
			theta=0;
		}
	}
	else if (init==1) //PF state.
	{
		theta=0;
	}
	else if (init==2) //Random orientation.
	{
		theta=M_PI*(2*ran()-1);
	}
	else
	{
		cerr << "BAD INIT VALUE: " << init << endl;
		abort();
	}
}

void particle::update(const double &v0, const double &eta, const int &LX, const int &LY)
{
	//Update angle.
	theta = theta_avg + 2*M_PI*eta*(ran()-0.5);
	
	if (theta<-M_PI)
	{
		theta+=2*M_PI;
	}
	else if (theta>M_PI)
	{
		theta-=2*M_PI;
	}
	
	//Update the displacement.
	dx += v0*cos(theta);
	dy += v0*sin(theta);
	
	//Update the position with the new angle (+periodic condition).
	x += v0*cos(theta);
	y += v0*sin(theta);	
	
	if (x>=LX)
	{
		x-=LX;
	}
	else if (x<0)
	{
		x+=LX;
	}
		
	if (y>=LY)
	{
		y-=LY;
	}
	else if (y<0)
	{
		y+=LY;
	}
}

///////////////////////////
///// BASIC FUNCTIONS /////
///////////////////////////

//Modulo function.
int modulo(const int &x, const int &L)
{
	if (x<0)
	{
		return x+L;
	}
	else if (x>=L)
	{
		return x-L;
	}
	else
	{
		return x;
	}
}

//Distance between two particles in periodic space.
double distance(const particle &part1, const particle &part2, const int &LX, const int &LY)
{
	const double DX=fabs(part1.x-part2.x);
	const double DY=fabs(part1.y-part2.y);	
	return square(min(DX,LX-DX)) + square(min(DY,LY-DY));
}

//Order parameter VS.
vector<double> VS(const int &Npart, const vector<particle> &PART)
{
	vector<double> vs(2,0.);
	for (int i=0; i<Npart; i++)
	{
		vs[0]+=cos(PART[i].theta);
		vs[1]+=sin(PART[i].theta);
	}
	vs[0]/=Npart;
	vs[1]/=Npart;
	return vs;	
}

//Order parameter VA.
vector<double> VA(const int &Npart, const vector<particle> &PART)
{
	vector<double> va(2,0.);
	for (int i=0; i<Npart; i++)
	{
		va[0]+=PART[i].spin*cos(PART[i].theta);
		va[1]+=PART[i].spin*sin(PART[i].theta);
	}
	va[0]/=Npart;
	va[1]/=Npart;
	return va;	
}

//Mean-square displacement.
vector<double> MSD(const int &Npart, const vector<particle> &PART)
{
	double DX=0, DY=0;
	double DX2=0., DY2=0.;
	for (int i=0; i<Npart; i++)
	{
		DX+=PART[i].dx/Npart;
		DX2+=PART[i].dx*PART[i].dx/Npart;
		DY+=PART[i].dy/Npart;
		DY2+=PART[i].dy*PART[i].dy/Npart;
	}	
	vector<double> DR(2,0.);
	DR[0]=DX2+DY2;
	DR[1]=DX*DX+DY*DY;
	return DR;
}

//Average angle near each particles at time t.
void THETA_AVG(const int &LX, const int &LY, const int &Npart, vector<particle> &PART)
{
	//Create a matrix with particle locations, box[i][j] regroups the particle indices with i<x<i+1 and j<y<j+1.
	vector< vector< vector<int> > > box(LX,vector< vector<int> >(LY,vector<int>(0)));
	for (int i=0; i<Npart; i++)
	{
		box[int(PART[i].x)][int(PART[i].y)].push_back(i);
	}
	
	//Calculate the average angle of neighbours.
	for (int i=0; i<Npart; i++)
	{
		const int X0=int(PART[i].x), Y0=int(PART[i].y);
		
		//Take the average magnetisation for particles in neighbour boxes and with a distance smaller than 1.	
		double MX=0., MY=0.;		
		for (int XN=X0-1;XN<=X0+1;XN++)
		{
			for (int YN=Y0-1;YN<=Y0+1;YN++)
			{
				const vector<int> neighbours=box[modulo(XN,LX)][modulo(YN,LY)];
				for (int l=0; l<neighbours.size(); l++)
				{
					const int j=neighbours[l];
					if (distance(PART[i],PART[j],LX,LY)<1)
					{
						int J=PART[i].spin*PART[j].spin;
						MX+=J*cos(PART[j].theta);
						MY+=J*sin(PART[j].theta);
					}
				}
			}
		}
		PART[i].theta_avg=atan2(MY,MX);
	}
}

//Export densities in a file.
void exportDensity(const double &v0, const double &eta, const double &rho0, const int &LX, const int &LY, const double &t, const int &init, const int &RAN, const int &Npart, const vector<particle> &PART)
{
	//Local density of each species.
	vector< vector<int> > rhoA(LX,vector<int>(LY,0)), rhoB(LX,vector<int>(LY,0));
	for (int i=0; i<Npart; i++)
	{
		if (PART[i].spin==1)
		{
			rhoA[int(PART[i].x)][int(PART[i].y)]++;
		}
		else
		{
			rhoB[int(PART[i].x)][int(PART[i].y)]++;
		}
	}
	
	//Creation of the file.
	static const int dossier=system("mkdir -p ./data_TSVM_dynamics/");
	stringstream ssDensity;
	ssDensity << "./data_TSVM_dynamics/TSVM_rho_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".txt";
	string nameDensity = ssDensity.str();
	ofstream fileDensity(nameDensity.c_str(),ios::trunc);
	fileDensity.precision(6);
	
	//Write in the file.
	for (int Y0=0;Y0<LY;Y0++)
	{
		for (int X0=0; X0<LX; X0++)
		{
			if (rhoA[X0][Y0]>rhoB[X0][Y0] or (rhoA[X0][Y0]==rhoB[X0][Y0] and ran()<0.5))
			{
				fileDensity << rhoA[X0][Y0] << "\t";
			}
			else
			{
				fileDensity << -rhoB[X0][Y0] << "\t";
			}
		}
		fileDensity << endl;
	}
	fileDensity.close();
}

//Read parameters from command line.
void ReadCommandLine(int argc, char** argv, double &v0, double &eta, double &rho0, int &LX, int &LY, int &tmax, int &init, int &RAN)
{
 	for( int i = 1; i<argc; i++ )
	{
		
		if (strstr(argv[i], "-v0=" ))
		{
			v0=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-eta=" ))
		{
			eta=atof(argv[i]+5);
		}
		else if (strstr(argv[i], "-rho0=" ))
		{
			rho0=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-LX=" ))
		{
			LX=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-LY=" ))
		{
			LY=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-tmax=" ))
		{
			tmax=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-init=" ))
		{
			init=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-ran=" ))
		{
			RAN=atoi(argv[i]+5);
		}
		else
		{
			cerr << "BAD ARGUMENT : " << argv[i] << endl;
			abort();
		}
	}
}

/////////////////////
///// MAIN CODE /////
/////////////////////

int main(int argc, char *argv[])
{
	//Physical parameters: v0=self-propulsion velocity, eta=noise parameter, rho0=average density, LX*LY=size of the box.
	double v0=0.5, eta=0.24, rho0=0.5;
	int LX=800, LY=100;
	
	//Numerical parameters: init=initial condition, tmax=maximal time, RAN=index of RNG.
	int init=2, tmax=50000, RAN=0;

	//Read imported parameters in command line.
	ReadCommandLine(argc,argv,v0,eta,rho0,LX,LY,tmax,init,RAN);
	
	//Number of particles.
	const int Npart=int(rho0*LX*LY);

	//Start the random number generator.
	init_gsl_ran();
	cout << "GSL index (initial time) = " << RAN << "\n";
	gsl_rng_set(GSL_r,RAN);
	
	//Creation of the file to export global averages.
	const int dossier=system("mkdir -p ./data_TSVM_averages/");
	stringstream ssAverages;
	ssAverages << "./data_TSVM_averages/TSVM_averages_v0=" << v0 << "_eta=" << eta << "_rho0=" << rho0 << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";
	string nameAverages = ssAverages.str();
	ofstream fileAverages(nameAverages.c_str(),ios::trunc);
	fileAverages.precision(6);
	
	//Initial state.
	vector<particle> PART;
	for (int k=0; k<Npart/2; k++)
	{
		particle partplus(LX,LY,+1,init);
		PART.push_back(partplus);
		particle partminus(LX,LY,-1,init);
		PART.push_back(partminus);
	}
	
	//Time evolution.
	for(int t=0; t<=tmax; t++)
	{
		//Export data.
		if (t%25==0)
		{
			//Order parameter VS.
			const vector<double> vs=VS(Npart,PART);
			const double vs_norm=sqrt(square(vs[0])+square(vs[1]));
			
			//Order parameter VA.
			const vector<double> va=VA(Npart,PART);
			const double va_norm=sqrt(square(va[0])+square(va[1]));
			
			//Mean-square displacement.
			const vector<double> msd=MSD(Npart,PART);
			
			fileAverages << t << " " << vs[0] << " " << vs[1] << " " << va[0] << " " << va[1] << " " << msd[0] << " " << msd[1] << endl;
			cout << double(100*t)/double(tmax) << "% steps performed -t=" << t << " -vs=" << vs_norm << " -va=" << va_norm << " -MSD=" << msd[0] << " -R2=" << msd[1] << running_time.TimeRun(" ") << endl;
			
			exportDensity(v0,eta,rho0,LX,LY,t,init,RAN,Npart,PART);
		}
		
		//Calculate the average orientation in the neighborhood of each particles.
		THETA_AVG(LX,LY,Npart,PART);
		
		//Update orientations and positions sequentially.
		for (int i=0; i<Npart; i++)
		{
			PART[i].update(v0,eta,LX,LY);			
		}
	}
	return 0;
}
