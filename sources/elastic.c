/*
 * main program file for elastic
 * (C) Copyright 2006 - , ICS-SU
 *
 * This file is part of elastic.
 *
 * elastic is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * elastic is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with elastic.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "elastic.h"
#include "ls_calls.h"


static void excfun(int sigid) {
  fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
  switch(sigid) {
    case SIGABRT:
      fprintf(stdout,"  Abnormal stop\n");  break;
    case SIGBUS:
      fprintf(stdout,"  Code error...\n");  break;
    case SIGFPE:
      fprintf(stdout,"  Floating-point exception\n"); break;
    case SIGILL:
      fprintf(stdout,"  Illegal instruction\n"); break;
    case SIGSEGV:
      fprintf(stdout,"  Segmentation fault.\n");  break;
    case SIGTERM:
    case SIGINT:
      fprintf(stdout,"  Programm killed.\n");  break;
  }
  fprintf(stdout," No data file saved.\n"); 
  exit(1);
}


static void usage(char *prog) {
  fprintf(stdout,"usage: %s [+/-v | -h] [-t typ] [-e err] [-n nit] source_file[.mesh] [-s data_file[.sol]] [-o output_file[.sol]]\n",prog);
  exit(1);
}


static int parsar(int argc,char *argv[],LSst *lsst) {
  int      i;
  char    *ptr;

  i = 1;
  while ( i < argc ) {
    if ( (*argv[i] == '-') || (*argv[i] == '+') ) {
      switch(argv[i][1]) {
      case 'h':  /* on-line help */
      case '?':
        usage(argv[0]);
        break;
      case 'e':
        if ( !strcmp(argv[i],"-e") ) {
          ++i;
          if ( isdigit(argv[i][0]) )
            lsst->sol.err = strtod(argv[i],NULL);
          else
            --i; 
        }
        break;
      case 'i':
        if ( !strcmp(argv[i],"-i") ) {
          ++i;
          lsst->mesh.name = argv[i];
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-n") ) {
          ++i;
          if ( i < argc && isdigit(argv[i][0]) )
            lsst->sol.nit = atoi(argv[i]);
          else
            --i; 
        }
        break;
      case 's':
        if ( !strcmp(argv[i],"-s") ) {
          ++i;
          lsst->sol.namein = argv[i];
          ptr = strstr(lsst->sol.namein,".sol");
          if ( !ptr )  strcat(lsst->sol.namein,".sol");
        }
        break;
      case 't':
        if ( ++i < argc ) {
          if ( isdigit(argv[i][0]) )
            lsst->info.typ = atoi(argv[i]);
          else { 
            fprintf(stderr,"%s: missing argument option %c\n",argv[0],argv[i-1][1]);
            usage(argv[0]);
            i--;
          }
        }
        else {
          fprintf(stderr,"%s: missing argument option %c\n",argv[0],argv[i-1][1]);
          usage(argv[0]);
        }
        break;
      case 'v':
        if ( !strcmp(argv[i],"-v") )
          lsst->info.verb = '0';
        else if ( !strcmp(argv[i],"+v") )
          lsst->info.verb = '+';
        else {
          fprintf(stderr,"%s: illegal option %s\n",argv[0],argv[i]);
          usage(argv[0]);
        }
        break;
      default:
        fprintf(stderr,"%s: illegal option %s\n",argv[0],argv[i]);
        usage(argv[0]);
      }
    }
    else {
      if ( lsst->mesh.name == NULL )
        lsst->mesh.name = argv[i];
      else {
        fprintf(stdout,"%s: illegal option %s\n",argv[0],argv[i]);
        usage(argv[0]);
      }
    }
    i++;
  }

  /* check params */
  if ( lsst->mesh.name == NULL ) {
    fprintf(stderr,"%s: missing argument\n",argv[0]);
    usage(argv[0]);
  }

  return(1);
}


static int parsop(LSst *lsst) {
  Cl         *pcl;
  Mat        *pm;
  float       fp1,fp2;
  int         i,j,ncld,npar,ret;
  char       *ptr,buf[256],data[256];
  FILE       *in;

  strcpy(data,lsst->mesh.name);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  strcat(data,".elas");
  in = fopen(data,"r");
  if ( !in ) {
    sprintf(data,"%s","DEFAULT.elas");
    in = fopen(data,"r");
    if ( !in )  return(1);
  }
  if ( lsst->info.verb != '0' )  fprintf(stdout,"    %s:",data);

  /* read parameters */
  lsst->sol.nbcl = 0;
	npar = 0;
  while ( !feof(in) ) {
    ret = fscanf(in,"%s",data);
    if ( !ret || feof(in) )  break;
    for (i=0; i<strlen(data); i++) data[i] = tolower(data[i]);

    /* check for condition type */
    if ( !strcmp(data,"dirichlet")  || !strcmp(data,"load") ) {
      fscanf(in,"%d",&ncld);
			npar++;
      for (i=lsst->sol.nbcl; i<lsst->sol.nbcl+ncld; i++) {
        pcl = &lsst->sol.cl[i];
        if ( !strcmp(data,"load") )             pcl->typ = Load;
        else  if ( !strcmp(data,"dirichlet") )  pcl->typ = Dirichlet;

        /* check for entity */
        fscanf(in,"%d %s ",&pcl->ref,buf);
        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);
        fscanf(in,"%c",&pcl->att);
        pcl->att = tolower(pcl->att);
        if ( (pcl->typ == Dirichlet) && (pcl->att != 'v' && pcl->att != 'f') ) {
          fprintf(stdout,"  %%%% Wrong format: %s %c\n",buf,pcl->att);
          continue;
        }
        else if ( (pcl->typ == Load) && (pcl->att != 'v' && pcl->att != 'f' && pcl->att != 'n') ) {
          fprintf(stdout,"  %%%% Wrong format: %s\n",buf);
          continue;
        }
        if ( !strcmp(buf,"vertices") || !strcmp(buf,"vertex") )          pcl->elt = LS_ver;
        else if ( !strcmp(buf,"edges") || !strcmp(buf,"edge") )          pcl->elt = LS_edg;
        else if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") )  pcl->elt = LS_tri;

        if ( pcl->att != 'f' && pcl->att != 'n' ) {
          for (j=0; j<lsst->info.dim; j++) {
            fscanf(in,"%f ",&fp1);
            pcl->u[j] = fp1;
          }
        }
        else if ( pcl->att == 'n' ) {
          fscanf(in,"%f ",&fp1);
          pcl->u[0] = fp1;
        }
      }
      lsst->sol.nbcl += ncld;
    }
    /* gravity or body force */
    else if ( !strcmp(data,"gravity") ) {
			npar++;
      lsst->info.load |= (1 << 0);
      for (j=0; j<lsst->info.dim; j++) {
        fscanf(in,"%f ",&fp1);
        lsst->info.gr[j] = fp1;
      }
    }
    else if ( !strcmp(data,"lame") ) {
			npar++;
      fscanf(in,"%d",&ncld);
      assert(ncld <= LS_MAT);
      lsst->sol.nmat = ncld;
      for (i=0; i<ncld; i++) {
        pm = &lsst->sol.mat[i];
        fscanf(in,"%d %f %f\n",&pm->ref,&fp1,&fp2);
        pm->lambda = fp1;
        pm->mu     = fp2;
      }
    }
    else if ( !strcmp(data,"youngpoisson") ) {
			npar++;
      fscanf(in,"%d",&ncld);
      lsst->sol.nmat = ncld;
      for (i=0; i<ncld; i++) {
        pm = &lsst->sol.mat[i];
        fscanf(in,"%d %f %f\n",&pm->ref,&fp1,&fp2);
        pm->lambda = (fp1 * fp2) / ((1.0+fp2) * (1.0-2.0*fp2));
        pm->mu     = fp1 / (2.0*( 1.0+fp2));
      }
    }
  }
  fclose(in);

  for (i=0; i<lsst->sol.nbcl; i++) {
    pcl = &lsst->sol.cl[i];
		lsst->sol.cltyp |= pcl->elt;
  }

  if ( npar > 0 && lsst->info.verb != '0' ) {
    fprintf(stdout," %d parameters\n",npar);
  }

  return(1);
}


int main(int argc,char **argv) {
  LSst     lsst;
	int      ier;
	char     stim[32];

  tminit(lsst.info.ctim,TIMEMAX);
  chrono(ON,&lsst.info.ctim[0]);

  /* trap exceptions */
  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);
  signal(SIGBUS,excfun);

  /* init structure */
  memset(&lsst.mesh,0,sizeof(Mesh));
  memset(&lsst.sol,0,sizeof(Sol));
	lsst.sol.cl  = (Cl*)calloc(LS_CL,sizeof(Cl));
  lsst.sol.mat = (Mat*)calloc(LS_MAT,sizeof(Mat));
  lsst.sol.err = LS_RES;
  lsst.sol.nit = LS_MAXIT;

  /* global parameters */
  lsst.info.dim    = 3;
	lsst.info.ver    = 1;
  lsst.info.verb   = '1';
  lsst.info.zip    = 0;
  lsst.info.typ    = P1;
  lsst.info.mfree  = 1;

  /* parse command line */
  if ( !parsar(argc,argv,&lsst) )  return(1);

  /* loading date */
  chrono(ON,&lsst.info.ctim[1]);

  if ( lsst.info.verb != '0' ) {
    fprintf(stdout," - ELASTIC, Release %s, %s\n   %s\n\n",LS_VER,LS_REL,LS_CPY);
    fprintf(stdout," - LOADING DATA\n");
  }

  /* parse parameters in file */
  if ( !parsop(&lsst) )  return(1);

  /* loading mesh */
  ier = loadMesh(&lsst);
	if ( ier <=0 )  return(1);

  /* counting P2 nodes */
	if ( lsst.info.typ == P2 )  
		lsst.info.np2 = lsst.info.dim == 2 ? hashar_2d(&lsst) : hashar_3d(&lsst);

  /* allocating memory */
  if ( !lsst.sol.u ) {
    lsst.sol.u  = (double*)calloc(lsst.info.dim * (lsst.info.npi+lsst.info.np2),sizeof(double));
    assert(lsst.sol.u);
  }

  /* loading solution (or Dirichlet values) */
  if ( lsst.sol.namein ) {
    ier = loadSol(&lsst);
    if ( !ier )  return(1);
  }

  /* packing mesh if needed */
  if ( lsst.sol.nmat ) {
    ier = lsst.info.dim == 2 ? pack_2d(&lsst) : pack_3d(&lsst);
		if ( ier == 0 ) {
			if ( lsst.info.verb != '0' )  fprintf(stdout," # Packing error.\n");
		  return(1);
		}
	}

	chrono(OFF,&lsst.info.ctim[1]);
	printim(lsst.info.ctim[1].gdif,stim);
  if ( lsst.info.verb != '0' )  fprintf(stdout," - COMPLETED: %s\n",stim);

  /* solve */
  chrono(ON,&lsst.info.ctim[2]);
  if ( lsst.info.verb != '0' )
    fprintf(stdout,"\n ** MODULE ELASTIC: %s\n",LS_VER);

	ier = LS_elastic(&lsst);
	if ( !ier )  return(1);
  chrono(OFF,&lsst.info.ctim[2]);
  if ( lsst.info.verb != '0' ) {
		printim(lsst.info.ctim[2].gdif,stim);
    fprintf(stdout," ** COMPLETED: %s\n\n",stim);
	}

  /* save file */
  if ( lsst.info.verb != '0' )  fprintf(stdout," - WRITING DATA\n");
  chrono(ON,&lsst.info.ctim[3]);
  if ( lsst.info.zip && !unpack(&lsst) )  return(1);

  if ( !lsst.sol.nameout ) {
    lsst.sol.nameout = (char *)calloc(128,sizeof(char));
    assert(lsst.sol.nameout);
    strcpy(lsst.sol.nameout,lsst.mesh.name);
  }

  ier = saveSol(&lsst);
	if ( !ier )   return(1);
  chrono(OFF,&lsst.info.ctim[3]);
  if ( lsst.info.verb != '0' ) {
    printim(lsst.info.ctim[3].gdif,stim);
    fprintf(stdout," - COMPLETED: %s\n",stim);
  }

  /* free mem */
	free(lsst.sol.u);

  chrono(OFF,&lsst.info.ctim[0]);
  if ( lsst.info.verb != '0' ) {
	  printim(lsst.info.ctim[0].gdif,stim);
    fprintf(stdout,"\n ** Cumulative time: %s.\n",stim);
  }

  return(0);
}

