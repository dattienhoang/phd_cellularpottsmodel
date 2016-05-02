/*******************************************************************************/
#include <set>
#include <utility>
#include <iostream>

//#include "potts_analysis_.h"

/*** HAMILTONIAN FUNCTIONS ***/

  double  inplaneEnergy(int,int);
  double  outplaneEnergy(int,int);
  double  Hamiltonian();
  double  interactionEnergy(int);
  double  volumeEnergy(int);
  double  anisotropyEnergy(int);
  double  blobularEnergy(int);

  double measureAnisotropy( int );

/*******************************************************************************/
/*** Returns the in-plane interaction energy of a lattice site ****/

double inplaneEnergy(int a, int b)
{
  double energy = 0.0;

  // if the site is a cell site
  if( lattice[a][b][0] != 0 ){

	// if the site's upper neighbor is not the same type as itself
    if( lattice[a][b][0] != lattice[(a+1)%N][b][0] ){
      // if the site's upper neighbor is a cell
      if( lattice[(a+1)%N][b][0] > 0 )
        energy += J_cel;
      // if the site's upper neighbor is air
      else
        energy += J_air;
    }

    // lower neighbor
    if( lattice[a][b][0] != lattice[(N+a-1)%N][b][0] ){
      if( lattice[(N+a-1)%N][b][0] > 0 )
        energy += J_cel;
      else
        energy += J_air;
    }

    // right neighbor
    if( lattice[a][b][0] != lattice[a][(b+1)%N][0] ){
      if( lattice[a][(b+1)%N][0] > 0 )
        energy += J_cel;
      else
        energy += J_air;
    }

    // left neighbor
    if( lattice[a][b][0] != lattice[a][(N+b-1)%N][0] ){
      if( lattice[a][(N+b-1)%N][0] > 0 )
        energy += J_cel;
      else
        energy += J_air;
    }

  }

  return energy;
}

/*******************************************************************************/
/*** Returns the out-of-plane interaction energy of a lattice site ****/

double outplaneEnergy(int a, int b)
{
  if( lattice[a][b][0]!=0 && lattice[a][b][1]!=0 )
    return J_col;
  else
    return 0.0;
}

/*******************************************************************************/
/*** Returns the full hamiltonian of the lattice ***/

double Hamiltonian()
{
  double energy = 0.0;
  for(int cell=1;cell<=numCells;cell++){
    energy += interactionEnergy(cell);
    energy += volumeEnergy(cell);
    energy += anisotropyEnergy(cell);
    energy += blobularEnergy(cell);
  }
  return energy;
}

/*******************************************************************************/
/*** Returns the interaction energy of a cell ****/

double interactionEnergy(int cell)
{
  double energy = 0.0;

  for(std::set< std::pair<int, int> >::iterator it = cellPerimeterList[cell].begin(); it!=cellPerimeterList[cell].end(); ++it)
	  energy += inplaneEnergy( it->first , it-> second );

  for(std::set< std::pair<int, int> >::iterator it = cellVolumeList[cell].begin(); it!=cellVolumeList[cell].end(); ++it)
	  energy += outplaneEnergy( it->first , it-> second );

  return energy;
}

/*******************************************************************************/
/*** Returns the volume energy of a cell ***/

double volumeEnergy(int cell)
{
  return L_vol*((double)cellVolumeList[cell].size()-targetVolume)*((double)cellVolumeList[cell].size()-targetVolume);
}

/*******************************************************************************/
/*** Returns the perimeter energy of a cell ***/

double anisotropyEnergy(int cell)
{
	//	return L_ani*(double)cellPerimeterList[cell].size()/(double)cellVolumeList[cell].size();
		return L_ani*measureAnisotropy( cell );
}

/*******************************************************************************/
/*** Returns the blobular energy of a cell ***/

double blobularEnergy(int cell)
{

	// (N is the total number of perimeter sites)
	// iterate over logN random perimeter site, and for each
	// iterate over logN random perimeter sites
	// and if the line between the site and the other site goes outside the cell
	// penalize

  int N = cellPerimeterList[cell].size();
  int logN = log( (double)N );
  int advanceAmount = logN;

  int energy = 0;
  int number = 0;

  int outerCount = 0;

  std::set< std::pair<int, int> >::const_iterator it1(cellPerimeterList[cell].begin());

  do{
	  advance(it1, advanceAmount);
	  outerCount += advanceAmount;

	  int ai = it1->first;
	  int aj = it1->second;

	  int count = 0;
	  std::set< std::pair<int, int> >::const_iterator it2(cellPerimeterList[cell].begin());
	  do{
		  advance(it2,advanceAmount);
		  count += advanceAmount;

		  int bi = it2->first;
		  int bj = it2->second;

		  int dx = bi-ai-N*(int)floor((float)(bi-ai)/(float)N+0.5);
		  int sx = (dx>0)-(dx<0);
		  int dy = bj-aj-N*(int)floor((float)(bj-aj)/(float)N+0.5);
		  int sy = (dy>0)-(dy<0);

		  double slope;
		  if(dx!=0)
			slope = (double)dy/(double)dx;
		  else
			slope = (double)N;

		  int x = ai;
		  int y = aj;
		  double error = fabs(slope);

		  do{
			number++;
			if(lattice[x][y][0]!=cell){
			  energy++;
			  goto done;
			}
			while(error>0.5){
			  y=(y+sy+N)%N;
			  number++;
			  if(lattice[x][y][0]!=cell){
				energy++;
				goto done;
			  }
			  error=error-1.0;
			}
			x=(x+sx+N)%N;
			error+=fabs(slope);
		  }while(x!=(ai+dx+N)%N);

		  done:;
    }while( count + advanceAmount < N );

  }while( outerCount + advanceAmount < N );


  return L_blb * (double)energy * log((double)energy) / (double)(cellPerimeterList[cell].size());

}

/*******************************************************************************/
