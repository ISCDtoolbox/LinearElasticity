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
  fprintf(stdout,"usage: %s [+/-v | -h | -x] [-n nit] [-r res] [-t typ] source[.mesh] [-p param[.elas]] [-s data[.sol]] [-o output[.sol]]\n",prog);
  fprintf(stdout,"\nOptions and flags:\n\
  --help       show the syntax and exit.\n\
  --version    show the version and date of release and exit.\n\n\
  -n nit       number of iterations max for convergence\n\
  -r res       value of the residual (Krylov space) for convergence\n\
  -t typ       specify the type of FE space: 1: P1, 2: P2\n\
  -v           suppress any message (for use with function call).\n\
  +v           increase the verbosity level for output.\n\
  -x           export (deformed) mesh\n\n\
  source.mesh    name of the mesh file\n\
  param.elas     name of file containing elasticity parameters\n\
  data.sol       name of file containing the initial solution or boundary conditions\n\
  output.sol     name of the output file (displacement field)\n");
  exit(1);
}


/* parse command line */
static int parsar(int argc,char *argv[],LSst *lsst) {
  int      i;
  char    *ptr,*data;

  i = 1;
  while ( i < argc ) {
    if ( (*argv[i] == '-') || (*argv[i] == '+') ) {
      switch(argv[i][1]) {
      case '-':
        if ( !strcmp(argv[i],"--help") )
          usage(argv[0]);
        else if ( !strcmp(argv[i],"--version") ) {
          fprintf(stdout,"%s: version: %s release: %s\n",argv[0],LS_VER,LS_REL);
          exit(1);
        }
        break;
      case 'h':  /* on-line help */
      case '?':
        usage(argv[0]);
        break;
      case 'i':
        if ( ++i < argc ) {
          lsst->mesh.name = argv[i];
          ptr = strstr(lsst->mesh.name,".mesh");
          if ( !ptr )  strcat(lsst->mesh.name,".mesh");
        }
        else {
          fprintf(stdout,"%s: missing input file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 'n':
        if ( ++i < argc && isdigit(argv[i][0]) )
          lsst->sol.nit = atoi(argv[i]);
        else {
          fprintf(stderr,"%s: missing argument option\n",argv[0]);
          usage(argv[0]);
        }
        break;
      case 'o':
        if ( ++i < argc ) {
          lsst->sol.nameout = argv[i];
          ptr = strstr(lsst->sol.nameout,".sol");
          if ( !ptr )  strcat(lsst->sol.nameout,".sol");
        }
        else {
          fprintf(stdout,"%s: missing data file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 'p':
        if ( ++i < argc ) {
          lsst->sol.namepar = argv[i];
          ptr = strstr(lsst->sol.namepar,".elas");
          if ( !ptr )  strcat(lsst->sol.namepar,".elas");
        }
        else {
          fprintf(stdout,"%s: missing parameter file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 'r':
        if ( ++i < argc && isdigit(argv[i][0]) )
          lsst->sol.res = strtod(argv[i],NULL);
        else {
          fprintf(stderr,"%s: missing argument option\n",argv[0]);
          usage(argv[0]);
        }
        break;
      case 's':
        if ( ++i < argc ) {
          lsst->sol.namein = argv[i];
          ptr = strstr(lsst->sol.namein,".sol");
          if ( !ptr )  strcat(lsst->sol.namein,".sol");
        }
        else {
          fprintf(stdout,"%s: missing data file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 't':
        if ( ++i < argc && isdigit(argv[i][0]) )
          lsst->info.typ = atoi(argv[i]);
        else { 
          fprintf(stderr,"%s: missing argument option\n",argv[0]);
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
      case 'x':
        lsst->info.xport = 1;
        break;
      default:
        fprintf(stderr,"%s: illegal option %s\n",argv[0],argv[i]);
        usage(argv[0]);
      }
    }
    else {
      if ( lsst->mesh.name == NULL ) {
        data = (char*)calloc(strlen(argv[i])+10,sizeof(char));
        strcpy(data,argv[i]);
        ptr = strstr(data,".mesh");
        if ( !ptr )  strcat(data,".mesh");
        lsst->mesh.name = data;
      }
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
  float       E,nu;
  int         i,j,ncld,npar,ret;
  char       *ptr,buf[256],data[256];
  FILE       *in;

  /* check for parameter file */
  if ( !lsst->sol.namepar ) {
    strcpy(data,lsst->mesh.name);
    ptr = strstr(data,".mesh");
    if ( ptr )  *ptr = '\0';
    strcat(data,".elas");
    in = fopen(data,"r");
    if ( !in ) {
      sprintf(data,"%s","DEFAULT.elas");
      in = fopen(data,"r");
    }
  }
  else {
    strcpy(data,lsst->sol.namepar);
    ptr = strstr(data,".elas");
    if ( !ptr )  strcat(data,".elas");
    in = fopen(data,"r");
  }
  if ( !in ) {
    if ( lsst->info.verb != '0' )  fprintf(stdout," # parameter file %s not found\n",data); 
    return(0);
  }
  if ( lsst->info.verb != '0' )  fprintf(stdout,"    %s:",data);

  /* read parameters */
  lsst->sol.nbcl = 0;
	npar = 0;
  while ( !feof(in) ) {
    /* scan line */
    ret = fscanf(in,"%s",data);
    if ( !ret || feof(in) )  break;
    for (i=0; i<strlen(data); i++) data[i] = tolower(data[i]);

    /* check for condition type */
    if ( !strcmp(data,"dirichlet") ) {
      fscanf(in,"%d",&ncld);
      npar++;
      for (i=lsst->sol.nbcl; i<lsst->sol.nbcl+ncld; i++) {
        pcl = &lsst->sol.cl[i];
        pcl->typ = Dirichlet;
        fscanf(in,"%d %s %c",&pcl->ref,buf,&pcl->att);
        
        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);
        pcl->att = tolower(pcl->att);
        if ( !strchr("fv",pcl->att) ) {
          if ( lsst->info.verb != '0' )  fprintf(stdout,"\n # wrong format: [%s] %c\n",buf,pcl->att);
          return(0);
        }
        if ( pcl->att == 'v' ) {
          for (j=0; j<lsst->info.dim; j++)  fscanf(in,"%lf",&pcl->u[j]);
        }
        if ( !strcmp(buf,"vertices") || !strcmp(buf,"vertex") )          pcl->elt = LS_ver;
        else if ( !strcmp(buf,"edges") || !strcmp(buf,"edge") )          pcl->elt = LS_edg;
        else if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") )  pcl->elt = LS_tri;
      }
      lsst->sol.nbcl += ncld;
    }
    /* traction forces */
    else  if ( !strcmp(data,"load") ) {
      fscanf(in,"%d",&ncld);
      npar++;
      for (i=lsst->sol.nbcl; i<lsst->sol.nbcl+ncld; i++) {
        pcl = &lsst->sol.cl[i];
        pcl->typ = Load;
        fscanf(in,"%d %s %c",&pcl->ref,buf,&pcl->att);
        
        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);
        pcl->att = tolower(pcl->att);
        if ( !strchr("fnv",pcl->att) ) {
          if ( lsst->info.verb != '0' )  fprintf(stdout,"\n # wrong format: [%s] %c\n",buf,pcl->att);
          return(0);
        }
        if ( pcl->att == 'v' ) {
          for (j=0; j<lsst->info.dim; j++)  fscanf(in,"%lf",&pcl->u[j]);
        }
        else if ( pcl->att == 'n' )  fscanf(in,"%lf ",&pcl->u[0]);
        
        if ( !strcmp(buf,"vertices") || !strcmp(buf,"vertex") )          pcl->elt = LS_ver;
        else if ( !strcmp(buf,"edges") || !strcmp(buf,"edge") )          pcl->elt = LS_edg;
        else if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") )  pcl->elt = LS_tri;

        /* for the time being: no normal et vertices known */
        if ( (pcl->elt == LS_ver) && (pcl->att == 'n') ) {
          if ( lsst->info.verb != '0' )  fprintf(stdout,"\n # condition not allowed: [%s] %c\n",buf,pcl->att);
          return(0);
        }
      }
      lsst->sol.nbcl += ncld;
    }
    /* gravity or body force */
    else if ( !strcmp(data,"gravity") ) {
			npar++;
      lsst->sol.cltyp |= Gravity;
      for (j=0; j<lsst->info.dim; j++)
        fscanf(in,"%lf",&lsst->sol.gr[j]);
    }
    else if ( !strcmp(data,"lame") ) {
			npar++;
      fscanf(in,"%d",&ncld);
      assert(ncld <= LS_MAT);
      lsst->sol.nmat = ncld;
      for (i=0; i<ncld; i++) {
        pm = &lsst->sol.mat[i];
        fscanf(in,"%d %lf %lf",&pm->ref,&pm->lambda,&pm->mu);
      }
    }
    else if ( !strcmp(data,"young") ) {
			npar++;
      fscanf(in,"%d",&ncld);
      lsst->sol.nmat = ncld;
      for (i=0; i<ncld; i++) {
        pm = &lsst->sol.mat[i];
        fscanf(in,"%d %f %f",&pm->ref,&E,&nu);
        pm->lambda = (E * nu) / ((1.0+nu) * (1.0-2.0*nu));
        pm->mu     = E / (2.0*( 1.0+nu));
      }
    }
  }
  fclose(in);
  if ( (npar > 0) && (lsst->info.verb != '0') )  fprintf(stdout," %d parameters\n",npar);

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
  memset(&lsst.info,0,sizeof(Info));
	lsst.sol.cl  = (Cl*)calloc(LS_CL,sizeof(Cl));
  lsst.sol.mat = (Mat*)calloc(LS_MAT,sizeof(Mat));
  lsst.sol.res = LS_RES;
  lsst.sol.nit = LS_MAXIT;

  /* global parameters */
  lsst.info.dim    = 3;
	lsst.info.ver    = 1;
  lsst.info.verb   = '1';
  lsst.info.zip    = 0;
  lsst.info.typ    = P1;
  lsst.info.mfree  = 1;
  lsst.info.xport  = 0;

  /* parse command line */
  if ( !parsar(argc,argv,&lsst) )  return(1);

  /* loading date */
  chrono(ON,&lsst.info.ctim[1]);

  if ( lsst.info.verb != '0' ) {
    fprintf(stdout," - ELASTIC, Release %s, %s\n   %s\n\n",LS_VER,LS_REL,LS_CPY);
    fprintf(stdout," - LOADING DATA\n");
  }

  /* loading mesh */
  ier = loadMesh(&lsst);
	if ( ier <=0 )  return(1);

  /* parse parameters in file */
  if ( !parsop(&lsst) )  return(1);

  /* counting P2 nodes */
	if ( lsst.info.typ == P2 )  
		lsst.info.np2 = lsst.info.dim == 2 ? hashar_2d(&lsst) : hashar_3d(&lsst);

  /* allocating memory */
  if ( !lsst.sol.u ) {
    lsst.sol.u  = (double*)calloc(lsst.info.dim * (lsst.info.np+lsst.info.np2),sizeof(double));
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
  if ( lsst.info.xport == 1 ) {
    ier = saveMesh(&lsst);
    if ( !ier )  return(0);
  }
  chrono(OFF,&lsst.info.ctim[3]);
  if ( lsst.info.verb != '0' ) {
    printim(lsst.info.ctim[3].gdif,stim);
    fprintf(stdout," - COMPLETED: %s\n",stim);
  }

  /* free mem */
	free(lsst.sol.u);
	free(lsst.sol.cl);
	free(lsst.sol.mat);

  chrono(OFF,&lsst.info.ctim[0]);
  if ( lsst.info.verb != '0' ) {
	  printim(lsst.info.ctim[0].gdif,stim);
    fprintf(stdout,"\n ** Cumulative time: %s.\n",stim);
  }

  return(0);
}

