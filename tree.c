/* this oct-tree code is adopted from SUBFIND with minor modifications for HBT */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "datatypes.h"
#include "intra_vars.h"
#include "iovars.h"
#include "proto.h"


void update_internal_nodes(HBTInt no,HBTInt sib,double len,HBTInt * PIndex,HBTReal PPos[][3])
{
	HBTInt j,jj,p,pp,mass,sons[8];
	double s[3];
	
		mass=0;
		s[0]=0;
		s[1]=0;
		s[2]=0;
		for(j=0;j<8;j++)
			sons[j]=Nodes[no].sons[j];//backup sons
		Nodes[no].way.len=len;
		Nodes[no].way.sibling=sib;
		for(j=0;sons[j]<0;j++);//find first son
		pp=sons[j];
		Nodes[no].way.nextnode=pp;
		for(jj=j+1;jj<8;jj++)//find sons in pairs,ie. find sibling
		{
			if(sons[jj]>=0)//ok, found a sibling
			{
				p=pp;
				pp=sons[jj];
				if(p<NumPart)
				{
					mass++;
					s[0]+=PPos[PIndex[p]][0];
					s[1]+=PPos[PIndex[p]][1];
					s[2]+=PPos[PIndex[p]][2];
					Nextnode[p]=pp;
				}
				else
				{
				if(len>=NodeResolution) 
					update_internal_nodes(p,pp,0.5*len,PIndex,PPos);//only divide if above resolution; otherwise 
				//we didn't divide the node seriouly so we don't have finer node length
				else
					update_internal_nodes(p,pp,len,PIndex,PPos);//get internal node info
				mass+=Nodes[p].way.mass;
				s[0]+=Nodes[p].way.s[0]*Nodes[p].way.mass;
				s[1]+=Nodes[p].way.s[1]*Nodes[p].way.mass;
				s[2]+=Nodes[p].way.s[2]*Nodes[p].way.mass;
				}
			}
		}
		if(pp<NumPart)//the last son
		{			
			mass++;
			s[0]+=PPos[PIndex[pp]][0];
			s[1]+=PPos[PIndex[pp]][1];
			s[2]+=PPos[PIndex[pp]][2];
			Nextnode[pp]=sib;
		}
		else
		{
			if(len>=NodeResolution) 
				update_internal_nodes(pp,sib,0.5*len,PIndex,PPos);//only divide if above resolution; otherwise 
			//we didn't divide the node seriouly so we don't have finer node length
			else
				update_internal_nodes(pp,sib,len,PIndex,PPos);
			mass+=Nodes[pp].way.mass;
			s[0]+=Nodes[pp].way.s[0]*Nodes[pp].way.mass;
			s[1]+=Nodes[pp].way.s[1]*Nodes[pp].way.mass;
			s[2]+=Nodes[pp].way.s[2]*Nodes[pp].way.mass;
		}
		Nodes[no].way.mass=mass;
		Nodes[no].way.s[0]=s[0]/mass;
		Nodes[no].way.s[1]=s[1]/mass;
		Nodes[no].way.s[2]=s[2]/mass;
}

HBTInt maketree(HBTInt halolen,HBTInt * PIndex,HBTReal PPos[][3])
{
	HBTInt NumNids,numnodes;
	HBTInt sub,subid,i,j,nodeid;
	double center[3], lenhalf;
	double xmin[3], xmax[3],Center[3], Len,Lenhalf;
	//~ HBTReal FilledFraction;

	NumPart=halolen;

	/* find enclosing rectangle */
  for(j = 0; j < 3; j++)
    xmin[j] = xmax[j] = PPos[PIndex[0]][j];

  for(i = 1; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      {
	if(PPos[PIndex[i]][j] > xmax[j])
	  xmax[j] = PPos[PIndex[i]][j];
	else if(PPos[PIndex[i]][j] < xmin[j])
	  xmin[j] = PPos[PIndex[i]][j];
      }

  //part_dens = NumPart / ((xmax[0] - xmin[0]) * (xmax[1] - xmin[1]) * (xmax[2] - xmin[2]));

  /* determine maxmimum extension */
  for(j = 1, Len = xmax[0] - xmin[0]; j < 3; j++)
    if((xmax[j] - xmin[j]) > Len)
      Len = xmax[j] - xmin[j];
  //Len *= 1.00001;

  for(j = 0; j < 3; j++)
    Center[j] = 0.5 * (xmax[j] + xmin[j]);
  
  Lenhalf=0.5*Len;
MaxNids=MaxNodes+NumPart;

Nodes= Nodes_base-NumPart;	/* select first node */


nodeid = NumPart;	/* id used to distinguish whether it's internal node or particle*/
NumNids=NumPart+1;	
	/* create an empty  root node  */
  for(i = 0; i < 8; i++)
Nodes_base->sons[i] = -1;

  for(i = 0; i < NumPart; i++)	/* insert all  particles */
	{
	  nodeid = NumPart ;	/* select index of first node in tree */
	    lenhalf = Lenhalf;
	  for(j = 0; j < 3; j++)
	    center[j] = Center[j];

	  while(1)
		{
			  //len = lenhalf;
			//fprintf(logfile,"%f\n",len);
			  lenhalf *= 0.5;//halflen for the to-be-found subnode
			  sub = 0;
	      if(PPos[PIndex[i]][0] > center[0])
		{
		  center[0] += lenhalf;//subcenter
		  sub += 1;//sub index
		}
	      else
		{
		  center[0] -= lenhalf;
		}
	      if(PPos[PIndex[i]][1] > center[1])
		{
		  center[1] += lenhalf;
		  sub += 2;
		}
	      else
		{
		  center[1] -= lenhalf;
		}
	      if(PPos[PIndex[i]][2] > center[2])
		{
		  center[2] += lenhalf;
		  sub += 4;
		}
	      else
		{
		  center[2] -= lenhalf;
		}
		
		subid=Nodes[nodeid].sons[sub];
		if(subid<0)//an empty node, insert particle as leaf
			{
			Nodes[nodeid].sons[sub]=i;
			break;//finished for this particle, begin to insert a new particle
			}
		else if(subid<NumPart)//a particle node, upgrade the node to internal
			{
			Nodes[nodeid].sons[sub]=NumNids;//create a new node;
			nodeid=NumNids;//take over the new nodeid
			NumNids++;
			if(NumNids >= MaxNids)
			{
			  fprintf(logfile,"maximum number %zd of tree-nodes reached.\n", MaxNodes);
			  fprintf(logfile,"for particle %d\n", i);fflush(logfile);
			  exit(1);
			}
			for(sub=0;sub<8;sub++)//initialize new node
				Nodes[nodeid].sons[sub]=-1;
			/*insert that subid into this new node*/
			//what if the two particles are too near? 
			//unnecessary to divide too fine, just get rid of one by random insertion. 
			 if(lenhalf < NodeReso_half)
				{
				/* seems like we're dealing with particles   
				* at identical locations. randomize 
				* sub index (well below gravitational softening scale).
				* to introduce some ambiguity so that we don't get jammed!*/
				sub = (HBTInt) (8.0 * drand48());
				if(sub >= 8)
				sub = 7;
				//~ fprintf(logfile,"len=%g Len=%g sub=%d  i=%d (%g|%g|%g)\n",
					//~ lenhalf*2, Len, sub, i, PPos[PIndex[i]][0], PPos[PIndex[i]][1], PPos[PIndex[i]][2]);
				}
			 else
				{
				sub=0;
				if(PPos[PIndex[subid]][0] > center[0])
					sub += 1;
				if(PPos[PIndex[subid]][1] > center[1])
					sub += 2;
				if(PPos[PIndex[subid]][2] > center[2])
					sub += 4;
				}	
			Nodes[nodeid].sons[sub]=subid;//the disturbing particle inserted
			}
		else nodeid=subid;//an internal node,take over it;
		}
	}
	
	numnodes=NumNids-NumPart;
	//~ FilledFraction=(HBTReal)numnodes/MaxNodes;
	//~ #pragma omp critical
	//~ if(FilledFraction>MaxNodeFilledFraction) MaxNodeFilledFraction=FilledFraction;
	//~ fprintf(logfile,"used %d nodes out of allocated %d. (filled fraction %g)\n",
	 //~ numnodes, MaxNodes, (double)numnodes / MaxNodes);	
	/* finished inserting, now update for walk*/
	update_internal_nodes(NumPart , -1, Len,PIndex,PPos);/*insert sibling and next infomation*/
	
	return numnodes; 
}

double tree_treeevaluate_potential(HBTReal targetPos[3], HBTInt *PIndex,HBTReal PPos[][3])
{
  union NODE *nop = 0;
  HBTInt no;
  double r2, dx, dy, dz, mass, r, u, h, h_inv, wp;
  double pot, pos_x, pos_y, pos_z;

  pos_x = targetPos[0];
  pos_y = targetPos[1];
  pos_z = targetPos[2];

  h = 2.8 * SofteningHalo;
  h_inv = 1.0 / h;

  pot = 0;

  no = NumPart;//start from root node

  while(no >= 0)
    {
      if(no < NumPart)		/* single particle */
	{
	  dx = PPos[PIndex[no]][0] - pos_x;
	  dy = PPos[PIndex[no]][1] - pos_y;
	  dz = PPos[PIndex[no]][2] - pos_z;
	  #ifdef PERIODIC_BDR
	  dx=NEAREST(dx);
	  dy=NEAREST(dy);
	  dz=NEAREST(dz);
	  #endif
	  mass = 1;
	  no = Nextnode[no];
	      r2 = dx * dx + dy * dy + dz * dz;	
	}
      else
	{
	  nop = &Nodes[no];
	  dx = nop->way.s[0] - pos_x;
	  dy = nop->way.s[1] - pos_y;
	  dz = nop->way.s[2] - pos_z;
	  #ifdef PERIODIC_BDR
	  dx=NEAREST(dx);
	  dy=NEAREST(dy);
	  dz=NEAREST(dz);
	  #endif
	  mass = nop->way.mass;
	  r2 = dx * dx + dy * dy + dz * dz;
		/* we have an internal node. Need to check opening criterion */
	  if((nop->way.len * nop->way.len )>( r2 * ErrTolTheta * ErrTolTheta))
	    {
	      /* open cell */
	      no = nop->way.nextnode;
	      continue;
	    }
	  no = nop->way.sibling;	/* node can be used */
	}
 
      r = sqrt(r2);

      if(r >= h)
	pot -= mass / r;
      else
	{
	  u = r * h_inv;

	  if(u < 0.5)
	    wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
	  else
	    wp =
	      -3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
						   u * (-16.0 + u * (9.6 - 2.133333333333 * u)));

	  pot += mass * h_inv * wp;
	}
    }

  return pot;
}


		

			
		
	
/* this function allocates memory used for storage of the tree
 * and auxiliary arrays for tree-walk and link-lists.
 */
size_t tree_tree_allocate(size_t maxnodes, size_t maxpart)	/* usually maxnodes=0.7*maxpart is sufficient */
{
  size_t bytes, allbytes = 0;

  MaxNodes = ((maxnodes>500)?maxnodes:500);

  if(!(Nodes_base =mymalloc(bytes = (MaxNodes + 1) * sizeof(union NODE))))
    {
      fprintf(logfile,"failed to allocate memory for %zd tree-nodes (%g MB).\n", MaxNodes, bytes / (1024.0 * 1024.0));fflush(logfile);
      exit(3);
    }
  allbytes += bytes;

  if(!(Nextnode =mymalloc(bytes = maxpart * sizeof(HBTInt))))
    {
      fprintf(logfile,"Failed to allocate %zd spaces for 'Nextnode' array (%g MB),now pause\n", maxpart,
	     bytes / (1024.0 * 1024.0));fflush(logfile);
	    exit(4);
    }
  allbytes += bytes;

  return allbytes;
}


/* free the allocated memory
 */
void tree_tree_free(void)
{
	free(Nextnode);
	free(Nodes_base);
}	
	














