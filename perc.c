//Authors Nathan Graves 21256779 && Matthew Hall 21718301
//
//Use Included MakeFile to Compile.
//Tested on macOS Sierra 10.12.6 with gcc 7.2 

//TODO - BOND PERC
//	-Occasional Conflicting Percolation Reports
//	 Small Sample Size Errors (ie running on 10*10 with 10 threads gives wrong results occasioanlly)
//	 calcualting verticle percolation in paralell (either need to run a dfs at when we get to root and make sure it spans the top peice or add in individual rowsOccupied account		ing fo piece
//	 Off sets )
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include "stack.h"
#include <string.h>
#include <omp.h>
#include <mpi.h>

// Will add one to number, wrapping around if necessary
int addOne(int number, int size){
    return (number+1) % size;
}

//
// Will subtract one from number, wrapping around if necessary
int minusOne (int number, int size){
    if(number == 0) return size-1;
    else return number-1;
}

float percProb(){
    return (double)rand() / (double)((unsigned)RAND_MAX + 1);
}

void  printLargestCluster(site **mat,int n ,int m, int ldx, int flag){ 
    for(int i =0; i<n; i++){
        for(int j =0 ; j < m; j++){
 		
	    if(flag==1 ){	
           	 if(mat[i][j].siteBond==ldx ) printf(GRN  "%c" RESET, 'x'); 	
		 else  printf(RED  "%c" RESET  , 'o');  
             }

	     else  printf(BLU "%d"RESET, mat[i][j].siteBond);
	    
       	
    }
	 printf("\n");
    }
}

void printMatrix(site **mat, int n){
	//SITES
    for(int i =0; i<n; i++){
        for(int j =0 ; j < n; j++){
          
            if(mat[i][j].siteBond) printf(GRN  "%d" RESET, mat[i][j].siteBond );
            else  printf(RED  "%d" RESET  , mat[i][j].siteBond );
            
        }
        printf("\n");
    }
	//BONDS
      for(int i =0; i<n; i++){     	  
	for(int j = 0 ; j < n ; j++){
		  if(mat[i][j].siteBond){
			printf(GRN"o"RESET);
			if(mat[i][j].rightBond){
				printf(GRN"-"RESET);
			}
			else{
				printf(RED"-"RESET);
			}
		   }
		   else {
			printf(RED"o"RESET);
				if(mat[i][j].rightBond){
				printf(GRN"-"RESET);
			}
			else{
				printf(RED"-"RESET);
			}

		   }
	}
	printf("\n");
	for(int j = 0 ; j < n ; j++){
		if(mat[i][j].lowerBond){	   
             		printf(GRN "| " RESET);
		}
		else {
			printf(RED "| " RESET);
		}
	   }
	printf("\n");
    }

}

//For AdjList Graph
     
void SeedMatrixGraph(Graph *g,int n,float p){
    
    //Seed RNG
    srand(time(NULL));
    
    //Seeding Matrix
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			//Top & Bottom Edges of respective pieces 
			if((int)(percProb() < p)) {addEdge(g,n*i+j,n*addOne(i,n)+j); }; 
			//Left and right Edges of respective pieces
			if((int)(percProb() < p)) {addEdge(g,n*i+j,n*i+addOne(j,n)); }; 
        }
    }
}


//Site Seeding 
void SeedMatrixSite(site **mat,int n, float p){
    
    //Seed RNG
    srand(time(NULL));
 
    //Seeding Matrix
    for(int i=0; i<n; i++){
	for(int j=0; j<n; j++){
            mat[i][j].siteBond =  (int)(percProb()< p);
			if(mat[i][j].siteBond){
				if(mat[addOne(i,n)][j].siteBond) {mat[i][j].lowerBond = 1; mat[addOne(i,n)][j].upperBond = 1;}
				if(mat[minusOne(i,n)][j].siteBond) {mat[i][j].upperBond = 1; mat[minusOne(i,n)][j].lowerBond = 1;}
				if(mat[i][addOne(j,n)].siteBond) {mat[i][j].rightBond = 1; mat[i][addOne(j,n)].leftBond = 1;}
				if(mat[i][minusOne(j,n)].siteBond) {mat[i][j].leftBond = 1; mat[i][minusOne(j,n)].rightBond = 1;}
			}
        }
    }
}
// Bond Seeding 
void SeedMatrixBond(site **mat,int n, float p){
    
    //Seed RNG
    srand(time(NULL));
    
    //Seeding Matrix
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if((int)(percProb() < p)) {mat[i][j].lowerBond = 1; mat[addOne(i,n)][j].upperBond = 1;
									   mat[i][j].siteBond = 1; mat[addOne(i,n)][j].siteBond = 1;}
            if((int)(percProb() < p)) {mat[i][j].rightBond = 1; mat[i][addOne(j,n)].leftBond = 1;
									   mat[i][j].siteBond = 1; mat[i][addOne(j,n)].siteBond = 1;}
        }
    }
}

void freePiece(piece *p, int width ){
	
 int numClusters = p->numClusters;

	for(int i =0 ; i < numClusters; i++ ){
	free(p->pieceClusters[i].colsOccupied);
	free(p->pieceClusters[i].rowsOccupied);
	}
  free(p->pieceClusters);


}

void initPiece(piece *p , size_t s, int width){
 p->pieceClusters = malloc(s * sizeof(cluster) + (sizeof(int) *width));
 p->used =0; 
 p->size =s;
 p->largestCluster=0;
 p->largestClusterIdx=0;
 p->percolates=0;
 p->numClusters=0;
}

void insertCluster(piece *p, int width,  cluster c) {
 
   if (p->used == p->size) {
    //just double every time    
    //
    p->size *= 2;
    p->pieceClusters = realloc(p->pieceClusters, p->size * (sizeof(cluster) + (sizeof(int) *width)));
  }

   p->pieceClusters[p->used++] =  c;

}

//Simulation 
void copyPiece(piece *p, piece *new, int width){

	int n = p->numClusters;
	for(int i = 0; i < n ; i++){
		insertCluster(new, width, p->pieceClusters[i] );
	}
	
	new->largestCluster = p->largestCluster;
	new->largestClusterIdx =  p->largestClusterIdx;
	new->percolates = p->percolates;
	new->numClusters = p->numClusters;
	
}
void initCluster(cluster *c, int width){
  
    //c->colsOccupied = malloc(sizeof(int) * width) ;
   // c->rowsOccupied = malloc(sizeof(int) * width) ;
 
    for(int i = 0 ; i < width; i++){
	c->colsOccupied[i]=0;
  	c->rowsOccupied[i]=0;
    }

    	c->parentClusID= -1;
	c->parentPieceID= -1;
	c->clusterID=0;
        c->clusHeight=0;
	c->clusWidth=0;
	c->clusSize=0;

}

int findCluster( int size, int l,  site **mat ,int print, int percCond, piece *p, int parallel,int tid){   

     Matrix visited = malloc(l *  sizeof(int *));
     int *temp = malloc( l  * size * sizeof(int));
     for(int i = 0; i < l; i++){
	 visited[i] = temp + ( i * size);
     }
     for(int i = 0 ; i < l ; i++){
	     for(int j =0; j< size; j++){
		visited[i][j] = 0;
	     }
     }

    int numClusters =0;
    int clusterIdx = 1 ; 	
    int lcidx =0;	
   
    //START OF MATRIX SEARCH
    for(int i =0; i< l; i++){
    	for(int j=0 ;j < size ;j++){    
         
	//Move along becasue this site is Part of another cluster        
	if(visited[i][j] == 1) continue; 
	//Mark as Visited       
	visited[i][j] = 1;    
	// If This site has a bond	
	if(mat[i][j].siteBond  == 1){	
	        Stack *s;
		s= malloc(sizeof(*s)*l);
		stackInit(s, l*size - p->largestCluster); 
		cluster c; 
  		initCluster(&c ,size);
		push(s, i * size  + j );
		
		c.colsOccupied[j] = 1;           
		c.rowsOccupied[i] = 1;   

	       	// Search for connected neighbours	
		 while(!isEmpty(s)) {
			int cur = pop(s);	
			int row = cur/size;
			int col = cur%size;
			int above = -1;
			int below = -1;
			c.clusSize++; 
			//printf("--- cur Poped =  mat[%d][%d] \n" ,row,col );

			//No Verticle Wrap around for individual parts in parallel
			if(parallel ==1 ){	
				if(row != 0){
			  		above = row -1;
				}
				 		
				if(row != l-1){	  
			  		below = row +1;	 
				}		
			}

			//Wrap around for sequential Matrix
			else{
			above =  (row-1 + l) % l ;
			below =  (row+1 + l) % l ;
			}
                        				
  			//Horizontal wrap arond Regardless.
			int left   =  (col-1 + size) % size ;
			int right =  (col+1 + size) % size ;
			
			//Identify Clusters For Peiceing them Together ( OR Displaying the largest cluster in Sequential Run)
		     	mat[row][col].siteBond = clusterIdx; 
			
			//CHECK NORTH SOUTH EAST AND WEST NEIGHBOURS
			
			//ABOVE-     Only Look Above if this isnt the top row, or its a sequential run
			if(above != -1){
				if( mat[row][col].upperBond  && !visited[above][col]) {
			        //printf("Thread %d Pushing Above Element  mat[%d][%d] (1d = ele %d) to the stack\n",tid ,above,col, above*size+col );
				push(s, above * size  + col );
			      	visited[above][col] = 1;
				c.rowsOccupied[above] = 1;
				}
			}
			//BELOW-      Only look below if this isnt the bottom row or its a sequential run
			if(below!= -1){		
				if( mat[row][col].lowerBond  && !visited[below][col]) {
			        //printf("Thread %d Pushing Below Element  mat[%d][%d] (1d = ele %d) to the stack \n",tid,below,col, below*size+col );
				push(s, below * size  + col );
			       	visited[below][col] = 1;
				c.rowsOccupied[below]= 1;
				}
			}
			//LEFT-        Look at left neighbour Always inc wraparound
				if( mat[row][col].leftBond  && !visited[row][left]) {
			        //printf("Thread %d Pushing Left Element mat[%d][%d] (1d = ele %d) to the stack \n",tid ,row,left, row*size+left );
				push(s,  row * size  + left );
			      	visited[row][left] = 1;
				c.colsOccupied[left] = 1;
				}
			//RIGHT-      Look at right Neigbour Always inc wraparound
				if( mat[row][col].rightBond  && !visited[row][right]) {
		  	        //printf("Thread %d Pushing Right Element  mat[%d][%d] (1d = ele %d) to the stack \n",tid, row,right, row*size+right ); 
		       		push(s, row * size  + right );
			       	visited[row][right] = 1;
				c.colsOccupied[right] = 1;
			}

		 }	
	 //freeStack
	 stackDestroy(s);	 
	 free(s);
       	
         //is this cluster bigger than all clusters in this peice
	 if(c.clusSize > p->largestCluster){
		 p->largestClusterIdx = clusterIdx;
		 p->largestCluster = c.clusSize;              
	 }
	//Reset Counters
	
	
	 for(int i =0; i< size; i++){
			if(c.rowsOccupied[i]==1)
			c.clusHeight++;

		if(c.colsOccupied[i]==1)
		        c.clusWidth++;
	 }



	 clusterIdx ++; 
         insertCluster(p , size ,c );
	 p->numClusters++;

	//END OF THIS CLUSTER
	 }

	//End Of Matrix  	

    	  
	}
     } 
      
    free(temp);  
    free(visited);
 //   printf("Thread %d Finished l is %d width is %d Largest Cluster in this peice is = %d\n",tid, l,size, p->largestCluster);   
   
    return 0;
}


//Add Unseen Columns Spanned Across this Cluster 
void addCols(piece *m, int pID, int cID, int pID2, int cID2 , int n){

	if(m[pID].pieceClusters[cID].clusWidth != n ){
		for(int r= 0 ; r < n ; r++){
			if(m[pID].pieceClusters[cID].colsOccupied[r] == 0 && m[pID2].pieceClusters[cID2].colsOccupied[r] == 1 ){
				m[pID].pieceClusters[cID].colsOccupied[r] = 1; 
				m[pID].pieceClusters[cID].clusWidth += 1; 
			}
		}
	}


}


//get the root of a cluster, ie parent = -1
void getRoot(piece *m, int curPID, int curCID, int *pID, int *cID){

	int pClusID = m[curPID].pieceClusters[curCID].parentClusID;
	int pPieceID =  m[curPID].pieceClusters[curCID].parentPieceID;

	while(m[pPieceID].pieceClusters[pClusID].parentClusID != -1 ){
		int t1= m[pPieceID].pieceClusters[pClusID].parentClusID;
		int t2 = m[pPieceID].pieceClusters[pClusID].parentPieceID;
			pPieceID =t2;
			pClusID = t1;

	}

	*pID =pPieceID;
	*cID =pClusID;

}

//Add To Piece Cluster and adjust parents Accordingly
void growCluster(piece *m, int pceID, int clusID, int pce2ID, int clus2ID, int n){

	//printf("Growing Cluster Adding PieceCluster to %d.%d size = %d from PieceCluster %d.%d size =%d \n",pceID,clusID,m[pceID].pieceClusters[clusID].clusSize ,pce2ID,clus2ID,
	//		m[pce2ID].pieceClusters[clus2ID].clusSize);
	m[pceID].pieceClusters[clusID].clusSize += m[pce2ID].pieceClusters[clus2ID].clusSize;
	m[pce2ID].pieceClusters[clus2ID].parentClusID  = clusID;
	m[pce2ID].pieceClusters[clus2ID].parentPieceID  = pceID;
	//add unseen Columns Spanned
	addCols(m, pceID, clusID, pce2ID, clus2ID , n);
}



void printUsage(){
		printf(BLU "Usage: ./perc {MatrixSize} {SeedingProbability} {Optional Flags}\n-p: Prints to console a visual representation of matrix- ONLY USE FOR SMALL MATRIX\n-v: Matrix only has to percolate vertically or horizontially -DEFAULT IS BOTH\n-b: Bond Percolation - DEFAULT IS SITE PERCOLATION\n-o Run with openMp\n-c Compare with regular\n"RESET );

}

site **alloc2d(int rows, int cols) {
    site *data = malloc(rows*cols*sizeof(site));
    site **array= malloc(rows*sizeof(int*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}



int main(int argc , char* argv[]){
	double startMPI, finish;
		//Init MPI
	MPI_Init(&argc , &argv );
	
	//Get Total Processes available to us
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// Get the name of the processor
    	char processor_name[MPI_MAX_PROCESSOR_NAME];
    	int name_len;
    	MPI_Get_processor_name(processor_name, &name_len);

	//Get each part to say hello first.   
	printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, world_rank+1, world_size);

	if(argc < 3){
		//Only one Proc Prints Error PLEASE!
		if(world_rank == MASTER) printUsage();
		//all exit tho
		exit(EXIT_FAILURE);
    	}
	//Command Line Options
   	int n = atoi(argv[1]);
    	float prob =atof(argv[2]);
       	int bPerc =0;
    	int percCond =0;
        int idx =3;
	
    	double time_taken; 

   	while( idx < argc){
		if(strncmp(argv[idx], "-v",2)==0) percCond =1;
		if(strncmp(argv[idx], "-b",2)==0) bPerc =1; 
			
    	idx++;
	}

	int matPartSize = n/world_size;	
    	int leftOvers= n - (world_size* matPartSize); 

	printf(GRN"matPartSize = %d , Lo = %d\n" RESET , matPartSize, leftOvers); 
	//Create Custom Struct Data Types for Piece and Sites

	//Site
	MPI_Datatype MPI_site;
	MPI_Datatype types[5] = {MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT};
	MPI_Aint disp[5]; //5 ints
	//Elements per Block
	int blckLen[5]= {1,1,1,1,1} ; //5 ints
	disp[0] = offsetof(site , upperBond ); 
	disp[1] = offsetof(site , lowerBond ); 
	disp[2] = offsetof(site , rightBond ); 
	disp[3] = offsetof(site , leftBond ); 
	disp[4] = offsetof(site , siteBond ); 
	
	MPI_Type_create_struct(5, blckLen, disp, types, &MPI_site);
	MPI_Type_commit(&MPI_site);


	//Cluster 
	MPI_Datatype MPI_cluster ;
	MPI_Datatype typs[8] = {MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT};
	int blkLen[8]= {1,1,1,1,1,n,n,1} ; 
	MPI_Aint disps[8] ;
	disps[0] = offsetof(cluster , clusterID ); 
	disps[1] = offsetof(cluster , parentClusID); 
	disps[2] = offsetof(cluster , parentPieceID ); 
	disps[3] = offsetof(cluster , clusHeight ); 
	disps[4] = offsetof(cluster , clusWidth); 
	disps[5] = offsetof(cluster , colsOccupied ); 
	disps[6] = offsetof(cluster , rowsOccupied ); 
	disps[7] = offsetof(cluster , clusSize ); 

	MPI_Type_create_struct(8, blkLen, disps, typs, &MPI_cluster);
	MPI_Type_commit(&MPI_cluster);

	//Master Sets and seeds Matrix
	if(world_rank == MASTER){		
	
	
	

	//Contigous to make sending a little easier
	site **mat = alloc2d(n,n);

	for(int i =0; i < n; i++){
		for(int j =0 ; j < n ; j++){
           		mat[i][j].upperBond = 0;
             		mat[i][j].lowerBond = 0;
             		mat[i][j].rightBond = 0;
             		mat[i][j].leftBond = 0;
             		mat[i][j].siteBond = 0;	
		}
   	}
	
	struct timeval start, end;
     	gettimeofday(&start, NULL);
  
    	if(bPerc){
 		printf("MASTER Seeding Bond Matrix.. %d elements Please Wait \n",n*n);
		SeedMatrixBond(mat,n,prob);
    	}
    	else{
		printf("MASTER Site Matrix.. %d elements Please Wait\n" ,n*n);
     		SeedMatrixSite(mat,n, prob); 
    	}
    	gettimeofday(&end, NULL);
   	time_taken = ((end.tv_sec  - start.tv_sec) * 1000000u +
		             end.tv_usec - start.tv_usec) / 1.e6;
 
    	printf("\rTime taken in Seeding Matrix is %12.10f\n", time_taken);
        
	//start timing after seeding for consistency
        startMPI=MPI_Wtime(); /*start timer*/


	//Matrix Seeded By MASTER 		
	for(int i = 0 ; i < world_size -1 ; i++){
	
		int start = matPartSize * i;
		int end = start + matPartSize;
     		if(i == world_size-1) end += leftOvers;
     		int pieceSize = end -start; 

		//SEND IT..
			
		MPI_Send(&(mat[start][0]),n*pieceSize, MPI_site,i+1,0, MPI_COMM_WORLD);
		
	}

		//Do My Bit
	//	printf("MASTER %d Starting work on mat[%d] to mat[%d]\n" , world_rank, 0 ,matPartSize );
		
		size_t initialSize = sizeof(int) + sizeof(cluster) + 2*n;
  		piece *m = malloc(sizeof(piece) * world_size);
		for(int i = 0 ; i < world_size; i++ ){	
			initPiece(&m[i] , initialSize ,n);
      		}   

		//Create Master Piece and fill it up
		
		size_t f = sizeof(int) + sizeof(cluster) + 2*n;
		findCluster(n , matPartSize,  mat , 0, 0, &m[0] ,0, 0);
               
		//Recv Full Pieces Back - Just Recieveing Size For Now
		for(int i = 0 ; i < world_size -1 ; i++){
			size_t psiz; 
			MPI_Status status;
			MPI_Recv(&psiz,1, my_MPI_SIZE_T,i+1,i+1, MPI_COMM_WORLD, &status);
		//	printf("Size of Piece from %d is %zu\n", i ,psiz);
			
		}
		printf("Done recieveing MASTER\n");

		//Simulate Recieving Pieces Back (As coudnt Get custom Data Type Working as intended)
		for(int i = 1 ; i < world_size ;i ++){		
    			copyPiece(&m[0] , &m[i] ,n);
		}
	
		printf("Done copying MASTER\n");

		//Join Pieces-
		
		  
   		int maxCluster = 0; 
   		int maxClusterIdx =0;
    
    		for(int i =world_size-1  ; i > 0; i--){
			int tr =   (matPartSize * i) ;
	
			for(int j = 0 ; j < n; j++){
				int ths = (mat[tr][j].siteBond) -1;
				int tht = (mat[tr -1 ][j].siteBond) -1;
		
				//Possible Join
				if(ths >= 0  && tht >=0 && mat[tr][j].upperBond){ 
				//If bottom peiceCluster has no parent
					if(m[i].pieceClusters[ths].parentClusID == -1){

						//if the above cluster doesnt have a parent 
						if(m[i-1].pieceClusters[tht].parentClusID == -1){
						growCluster(m,i-1,tht,i,ths,n);					
						}

						//if the above cluster does have a parent 
						else {
						int pClusID  =	m[i-1].pieceClusters[tht].parentClusID;
						int pPieceID =  m[i-1].pieceClusters[tht].parentPieceID;

						getRoot(m, i-1, tht, &pPieceID, &pClusID );
						growCluster(m,pPieceID,pClusID,i,ths,n);
						}
					}	
					//Bottom PieceCluster has a parent
					else if(m[i].pieceClusters[ths].parentClusID != -1) {
				
					int pClusID = m[i].pieceClusters[ths].parentClusID;
					int pPieceID = m[i].pieceClusters[ths].parentPieceID;
					getRoot(m, i, ths ,&pPieceID, &pClusID);

					//top Row has no parent
					if(m[i-1].pieceClusters[tht].parentClusID == -1 ){
					//If Were Not the child of the above pieceCluster
						if(!(pClusID == tht && pPieceID == i-1)){
							growCluster(m,pPieceID,pClusID,i-1,tht,n);
						}
			   		}
				
					//topRow Has A Parent
					else{
			
					int dlt= m[i-1].pieceClusters[tht].parentClusID;
					int dlt2 = m[i-1].pieceClusters[tht].parentPieceID;
					getRoot(m, i-1, tht, &dlt2, &dlt );

						if(!(dlt == pClusID && dlt2 == pPieceID)){
							growCluster(m,pPieceID,pClusID,dlt2,dlt,n);						
						}
	
					}
				
					}
			
				}

		   }	   		
      	} 
         

      	//Dealing with wrap around Seperatly
      	for(int j = 0; j < n ; j++){
		int btRow = (mat[n-1][j].siteBond) -1;
		int tpRow = (mat[0][j].siteBond) -1;
		//If there is a connection
		if( btRow >= 0 && tpRow >=0  && mat[n-1][j].lowerBond){
		//if bottom doesnt have a parent
			if(m[world_size-1].pieceClusters[btRow].parentClusID == -1){ 

				//Top row doesn't have a parent
				if(m[0].pieceClusters[tpRow].parentClusID == -1 ){
					growCluster(m,0,tpRow,world_size-1,btRow,n);	
				}

				//Top Row has a Parent
				else  {
				int delt = m[0].pieceClusters[tpRow].parentClusID;
				int delt2 = m[0].pieceClusters[tpRow].parentPieceID;
				getRoot(m, 0, tpRow, &delt2, &delt );
				growCluster(m,delt2,delt,world_size-1,btRow,n);	

				}
			}
			//Bottom Row Already has a parent
			else  {
			int parentClus = m[world_size-1].pieceClusters[btRow].parentClusID;
			int parentPiece = m[world_size-1].pieceClusters[btRow].parentPieceID;
			getRoot(m, world_size-1, btRow, &parentPiece, &parentClus );
		
				if(m[0].pieceClusters[tpRow].parentClusID == -1 && parentClus != tpRow ){
					growCluster(m,parentPiece,parentClus,0,tpRow,n);	
				}	

				else if(m[0].pieceClusters[tpRow].parentClusID != -1){
					int delt = m[0].pieceClusters[tpRow].parentClusID;
					int delt2 = m[0].pieceClusters[tpRow].parentPieceID;
					getRoot(m, 0, tpRow, &delt2, &delt );

					if(parentClus != delt ){
					growCluster(m,delt2,delt,parentPiece,parentClus,n);
	
					}

				}		
			}
	  	}	
		
      	}
      	int lc =0;
      	int cIndex=0;
      	int pIndex=0;
      	int vperc=0;
      	int hperc =0;
      	int fullPerc =0;
     	//Check For Percolation 
     	for(int cls =0; cls<n; cls++){
		if(m[world_size-1].pieceClusters[cls].clusHeight == matPartSize+leftOvers){
			int thisC = m[world_size -1].pieceClusters[cls].parentClusID;
			int thisP = m[world_size -1].pieceClusters[cls].parentPieceID;
			getRoot(m, world_size-1, cls, &thisP, &thisC);

			if(thisP == 0 && m[thisP].pieceClusters[thisC].clusSize > matPartSize ){
				vperc = 1;
				if(m[thisP].pieceClusters[thisC].clusWidth == n ){
					fullPerc = 1;
					//There exists a cluster which spans both axis no need to keep looking for one.
					break;	
			 }	
			

			}
		}	
    	 }	     

      	//Find Largest Cluster 
     	for( int pce =0 ; pce < world_size; pce++){
     	  	for(int j = 0; j < m[pce].numClusters; j++){
			if(m[pce].pieceClusters[j].clusSize > lc   ){
				lc =m[pce].pieceClusters[j].clusSize ;
				cIndex = j;
				pIndex = pce;	
			}
			if(m[pce].pieceClusters[j].clusWidth == n  ){
				hperc =1 ;			

			}

          	}
      	}
   
   	if(!percCond){
		if(fullPerc){
		 printf(YEL "Matrix Percolates Largest CLuster is %d \n"RESET , lc); 
		}
		else printf(YEL "Matrix Does Not Percolates Largest CLuster is %d \n"RESET , lc); 
      	}
      	
     		else{
	      		if(hperc || vperc){
	 			printf(YEL "Matrix Percolates Largest CLuster is %d \n"RESET , lc); 
			}
			else printf(YEL "Matrix Does Not Percolate Largest CLuster is %d \n"RESET , lc);
     		 }
        

	//Joins Finished

	

	//MASTER Free Memory 
	free(m);
     	free(mat[0]);
     	free(mat);
   	//END OF MASTER WORK
	
	}


	
	if(world_rank != MASTER ){
				
       		int start = matPartSize * world_rank;
		int end = start + matPartSize;
     		if(world_rank == world_size-1) end += leftOvers;
     		int Height = end -start; 
		int Width  = n; // 
	//	printf(RED "RANK:%d leftovers =%d start %d matps = %d , end %d height%d width %d \n" RESET, world_rank, leftOvers ,start, matPartSize, end,Height,Width  );


      	//	printf("n is %d Processor %d ready to recieve work expecting  mat[%d] to mat[%d]\n" ,n, world_rank, start ,end );
        
		site **mat = alloc2d(Height,Width);
        	MPI_Status status;
		int numberOfSitesRead;
		MPI_Recv(&(mat[0][0]),Height*Width , MPI_site, 0,0, MPI_COMM_WORLD,&status);
		MPI_Get_count(&status, MPI_site, &numberOfSitesRead);

		//Create A Piece Fill It with Info Needed for master
		piece p;
		size_t is = sizeof(int) + sizeof(cluster) + 2*Width;
		initPiece(&p, is , Width );
		findCluster(Width , Height,  mat , 0, 0, &p ,0, 0); 
	//	printf("Part %d Found All Clusters in My Piece\n",world_rank);

    	  	//SEND PIECE BACK HERE --Just Sending Size for simulation Purposes
		
		size_t psiz = p.size;
		MPI_Send(&psiz,1, my_MPI_SIZE_T,0,world_rank, MPI_COMM_WORLD);
		printf("Proc %d sent size %zu to Master \n", world_rank, psiz);
		//Sent Size Across - THis Structure is a little more complicated



		free(mat[0]);
		free(mat);

	}
	
      
	if(world_rank ==MASTER){
	finish=MPI_Wtime();
	printf(RED "Parallel Elapsed time: %f seconds\n" RESET, finish-startMPI); 
	}
	//Finalize MPI
	MPI_Finalize();

	
	
     	return 0; 
}
