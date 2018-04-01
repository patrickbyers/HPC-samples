/**
 * University of Pittsburgh
 * Department of Computer Science
 * CS1645: Introduction to HPC Systems
 * Student: 
 * Instructor: Bryan Mills, University of Pittsburgh
 * MPI particle-interaction code. 
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TAG 7
#define CONSTANT 777

// Particle-interaction constants
#define A 10250000.0
#define B 726515000.5
#define MASS 0.1
#define DELTA 1

// Random initialization constants
#define POSITION 0
#define VELOCITY 1

// Structure for shared properties of a particle (to be included in messages)
struct Particle{
  float x;
  float y;
  float mass;
  float fx;
  float fy;
};

// Headers for auxiliar functions
float random_value(int type);
void print_particles(struct Particle *particles, int n);
void interact(struct Particle *source, struct Particle *destination);
void compute_interaction(struct Particle *source, struct Particle *destination, int limit);
void compute_self_interaction(struct Particle *set, int size);
void merge(struct Particle *first, struct Particle *second, int limit);
int read_file(struct Particle *set, int size, char *file_name);

// Main function
main(int argc, char** argv){
  int myRank;// Rank of process
  int p;// Number of processes
  int n;// Number of total particles
  int previous;// Previous rank in the ring
  int next;// Next rank in the ring
  int tag = TAG;// Tag for message
  int number;// Number of local particles
  struct Particle *globals;// Array of all particles in the system
  struct Particle *locals;// Array of local particles
  struct Particle *remotes;// Array of foreign particles
  char *file_name;// File name
  MPI_Status status[2];// Return status for receive
  MPI_Request request[2];// Isend status for send
  int j, rounds, initiator, sender;
  double start_time, end_time;
   

  // checking the number of parameters
  if(argc < 2){
    printf("ERROR: Not enough parameters\n");
    printf("Usage: %s <number of particles> [<file>]\n", argv[0]);
    exit(1);
  }
  
  // getting number of particles
  n = atoi(argv[1]);

  // initializing MPI structures and checking p is odd  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Datatype ParticleStruct, ParticleType;
  //printf("p=%d\n",p);
  if(p % 2 == 0){
    p = p - 1;
    if(myRank == p){
      MPI_Finalize();
      return 0;
    }
  }
  srand(myRank+myRank*CONSTANT);
  
  previous = myRank - 1;
  next = myRank + 1;
  
  if (myRank == 0) 
    previous = p - 1;
  if (myRank == (p - 1))
    next = 0;
  
  //printf("myRank=%d, next=%d, prev=%d\n",myRank, next, previous);
 

  // acquiring memory for particle arrays
  if(n%p == 0) {
    number = n / p;
  }
  else {
    number = (n/p) + 1;
  }
  locals = (struct Particle *) malloc(number * sizeof(struct Particle));
  remotes = (struct Particle *) malloc(number * sizeof(struct Particle));
  
  /*
  for(j = 0; j < number; j++) {
    printf("Rank %d - init local particles[%d] <%f, %f, %f, %f, %f> 0x%08x - 0x%08x\n", 
    myRank,
    j,
    locals[j].x, 
    locals[j].y, 
    locals[j].mass,
    locals[j].fx,
    locals[j].fy,
    &locals[j].x,
    &locals[j].fy);
  }
  */
  
  // checking for file information
  if(argc == 3){
    if(myRank == 0){
      globals = (struct Particle *) malloc((n+p) * sizeof(struct Particle));

      // YOUR CODE GOES HERE (reading particles from file)
      read_file(globals,n,argv[2]);
      
      if(n%p != 0) {
        for(j = n; j < (n+p)-(n%p); j++) {
          globals[j].x = 0.0;
          globals[j].y = 0.0;
          globals[j].fx = 0.0;
          globals[j].fy = 0.0;
          globals[j].mass = 0.0;
        }
        n = (n+p)-(n%p);
      }
    }
    


    
    // To send/recv (or scatter/gather) you will need to learn how to
    // transfer structs of floats, treat it as a contiguous block of
    // floats. Here is an example:
    // MPI_Send(locals,
    //          number * (sizeof (struct Particle)) / sizeof(float),
    //          MPI_FLOAT,
    //          next_rank,
    //          tag,
    //          MPI_COMM_WORLD)
    // MPI_Recv(remotes,
    //          number * (sizeof (struct Particle)) / sizeof(float),
    //          MPI_FLOAT,
    //          previous_rank,
    //          tag,
    //          MPI_COMM_WORLD,
    //          &status);
    // hint: because your nodes need to both send and receive you
    // might consider asyncronous send/recv.

    // YOUR CODE GOES HERE (distributing particles among processors)
    
    
    //MPI_Datatype ParticleStruct, ParticleType;
    // This defines the mapping between the struct and the MPI Datatypes.
    MPI_Datatype types[5] = {MPI_DOUBLE,
         MPI_DOUBLE,
         MPI_DOUBLE, 
         MPI_DOUBLE,
         MPI_DOUBLE};

    // This is the length of each element, all 1 in this case.
    int         blocklens[5] = {1, 1, 1, 1, 1}; 

    // This is the offsets defined in the struct.
    MPI_Aint     displacements[5];

    // Get the displacements using the particles array. Note that the
    // first displacment would be the pointer value of &particles. Then each
    // subsequent element in the struct gives us the displacement.
    MPI_Get_address(globals,          displacements); 
    MPI_Get_address(&globals[0].y,    displacements+1); 
    MPI_Get_address(&globals[0].mass, displacements+2); 
    MPI_Get_address(&globals[0].fx,   displacements+3); 
    MPI_Get_address(&globals[0].fy,   displacements+4); 

    // Now convert these addresses into actual displacements.
    MPI_Aint     base;
    base = displacements[0];
    int i;
    for (i=0; i < 5; i++) {
      displacements[i] = displacements[i] - base;
      //printf("displacements[%d] = %d\n",i,displacements[i]);
    }

    // Now we create the particle struct.
    MPI_Type_create_struct(5, blocklens, displacements, types, &ParticleStruct);

    // Once we have that we can get the size of a single particle.
    MPI_Aint sizeOfParticle;
    MPI_Get_address(globals+1, &sizeOfParticle);
    sizeOfParticle = sizeOfParticle - base;
    //printf("sizeOfParticle = %d\n", sizeOfParticle);

    // We can then create an optimized version of this type.
    MPI_Type_create_resized(ParticleStruct, 0, sizeOfParticle, &ParticleType);

    // Lastly, we commit this type to MPI for use!
    MPI_Type_commit(&ParticleType);

    
    
    // Distribute the particles
    //MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Bcast(&globals, n, ParticleType, 0, MPI_COMM_WORLD);
    MPI_Scatter(globals, number, ParticleType, locals, number, ParticleType, 0, MPI_COMM_WORLD);
    
    // END ADDED CODE (distributing particles)
    
  
  } 
  else {
    // random initialization of local particle array
    if(n%p == 0 || myRank <= n%p - 1) {
      for(j = 0; j < number; j++){
        locals[j].x = random_value(POSITION);
        locals[j].y = random_value(POSITION);
        locals[j].fx = 0.0;
        locals[j].fy = 0.0;
        locals[j].mass = MASS;
      }
    }
    else {
      for(j = 0; j < number-1; j++){
        locals[j].x = random_value(POSITION);
        locals[j].y = random_value(POSITION);
        locals[j].fx = 0.0;
        locals[j].fy = 0.0;
        locals[j].mass = MASS;
      }
      locals[number-1].x = 0.0;
      locals[number-1].y = 0.0;
      locals[number-1].fx = 0.0;
      locals[number-1].fy = 0.0;
      locals[number-1].mass = 0.0;
    }
    if(n%p != 0) {
      n = (n+p)-(n%p);
    }  
  }
  
/*  
  for(j = 0; j < number; j++) {
    printf("Rank %d - local particles[%d] <%f, %f, %f, %f, %f> 0x%08x - 0x%08x\n", 
    myRank,
    j,
    locals[j].x, 
    locals[j].y, 
    locals[j].mass,
    locals[j].fx,
    locals[j].fy,
    &locals[j].x,
    &locals[j].fy);
  }
  */
  
  
  // starting timer
  if(myRank == 0){
    start_time = MPI_Wtime();
  }
  
  // YOUR CODE GOES HERE (ring algorithm)
  
  /*
  for(j = 0; j < number; j++) {
    printf("Rank %d - remote particles[%d] <%f, %f, %f, %f, %f>\n", 
    myRank,
    j,
    remotes[j].x, 
    remotes[j].y, 
    remotes[j].mass,
    remotes[j].fx,
    remotes[j].fy);
  }
  */

  MPI_Isend(locals, number, ParticleType, next, TAG, MPI_COMM_WORLD, &request[0]);
  MPI_Irecv(remotes, number, ParticleType, previous, TAG, MPI_COMM_WORLD, &request[1]);
  // Wait for receive to finish before doing calculations
  MPI_Waitall(2, request, status);
  // Compute interaction with locals and remotes (in that order)
  compute_interaction(locals, remotes, number);
  
  
  /*
  for(j = 0; j < number; j++) {
    printf("Rank %d - 1 calc local particles[%d] <%f, %f, %f, %f, %f> 0x%08x - 0x%08x\n", 
    myRank,
    j,
    locals[j].x, 
    locals[j].y, 
    locals[j].mass,
    locals[j].fx,
    locals[j].fy,
    &locals[j].x,
    &locals[j].fy);
  }
  */

  
  // send remotes forward and repeat steps for total of (p-1)/2 steps
  for (j = 1; j  < (p-1)/2; j++) {
    MPI_Isend(remotes, number, ParticleType, next, TAG, MPI_COMM_WORLD, &request[0]);
    MPI_Irecv(remotes, number, ParticleType, previous, TAG, MPI_COMM_WORLD, &request[1]);
    // Wait for receive to finish before doing calculations
    MPI_Waitall(2, request, status);
    // Compute interaction with locals and remotes (in that order)
    compute_interaction(locals, remotes, number);
  }
  // send the remotes back to their original location by reversing the send/recv (p-1)/2 times
  for (j = 0; j < (p-1)/2; j++) {
    MPI_Isend(remotes, number, ParticleType, previous, TAG, MPI_COMM_WORLD, &request[0]);
    MPI_Irecv(remotes, number, ParticleType, next, TAG, MPI_COMM_WORLD, &request[1]);
    // Wait for receive to finish before doing calculations
    MPI_Waitall(2, request, status);
  }
  // call functions merge and then compute self interaction
  merge(locals, remotes, number);
  compute_self_interaction(locals, number);
  
  // END ADDED CODE (ring algorithm)

  // stopping timer
  if(myRank == 0){
    end_time = MPI_Wtime();
    printf("Duration: %f seconds\n", (end_time-start_time));
  }
  
  // printing information on particles
  if(argc == 3){
    
    // YOUR CODE GOES HERE (collect particles at rank 0)
    MPI_Gather(locals, number, ParticleType, globals, number, ParticleType, 0, MPI_COMM_WORLD);

    if(myRank == 0) {
      print_particles(globals,n);
    }
  }
  
  // finalizing MPI structures
  MPI_Finalize();
 
}

// Function for random value generation
float random_value(int type){
  float value;
  switch(type){
  case POSITION:
    value = (float)rand() / (float)RAND_MAX * 100.0;
    break;
  case VELOCITY:
    value = (float)rand() / (float)RAND_MAX * 10.0;
    break;
  default:
    value = 1.1;
  }
  return value;
}

// Function for printing out the particle array
void print_particles(struct Particle *particles, int n){
  int j;
  printf("Index\tx\ty\tmass\tfx\tfy\n");
  for(j = 0; j < n; j++){
    printf("%d\t%f\t%f\t%f\t%f\t%f\n",j,particles[j].x,particles[j].y,particles[j].mass,particles[j].fx,particles[j].fy);
  }
}

// Function for computing interaction among two particles
// There is an extra test for interaction of identical particles, in which case there is no effect over the destination
void interact(struct Particle *first, struct Particle *second){
  float rx,ry,r,fx,fy,f;

  // computing base values
  rx = first->x - second->x;
  ry = first->y - second->y;
  r = sqrt(rx*rx + ry*ry);

  if(r == 0.0)
    return;

  f = A / pow(r,6) - B / pow(r,12);
  fx = f * rx / r;
  fy = f * ry / r;

  // updating sources's structure
  first->fx = first->fx + fx;
  first->fy = first->fy + fy;
  
  // updating destination's structure
  second->fx = second->fx - fx;
  second->fy = second->fy - fy;

}

// Function for computing interaction between two sets of particles
void compute_interaction(struct Particle *first, struct Particle *second, int limit){
  int j,k;
  
  for(j = 0; j < limit; j++){
    for(k = 0; k < limit; k++){
      interact(&first[j],&second[k]);
    }
  }
}

// Function for computing interaction between two sets of particles
void compute_self_interaction(struct Particle *set, int size){
  int j,k;
  
  for(j = 0; j < size; j++){
    for(k = j+1; k < size; k++){
      interact(&set[j],&set[k]);
    }
  }
}

// Function to merge two particle arrays
// Permanent changes reside only in first array
void merge(struct Particle *first, struct Particle *second, int limit){
  int j;
  
  for(j = 0; j < limit; j++){
    first[j].fx += second[j].fx;
    first[j].fy += second[j].fy;
  }
}

// Reads particle information from a text file
int read_file(struct Particle *set, int size, char *file_name){
  FILE *ifp, *ofp;
  char *mode = "r";
  ifp = fopen(file_name, mode);

  if (ifp == NULL) {
    fprintf(stderr, "Can't open input file!\n");
    return 1;
  }

  int i=0;
  // reading particle values
  for(i=0; i<size; i++){
    fscanf(ifp, "%f\t%f\t%f", &set[i].x, &set[i].y, &set[i].mass);
    set[i].fx = 0.0;
    set[i].fy = 0.0;
  }
  
  // closing file
  fclose(ifp);

  return 0;
}

