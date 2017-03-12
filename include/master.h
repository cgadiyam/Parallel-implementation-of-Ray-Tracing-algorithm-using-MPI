#ifndef __MASTER_PROCESS_H__
#define __MASTER_PROCESS_H__

#include "RayTrace.h"

//This function is the main that only the master process
//will run.
//
//Inputs:
//    data - the ConfigData that holds the scene information.
//
//Outputs: None
void masterMain( ConfigData *data );

//This function will perform ray tracing when no MPI use was
//given.
//
//Inputs:
//    data - the ConfigData that holds the scene information.
//
//Outputs: None
void masterSequential(ConfigData *data, float* pixels);

void static_strips_horizontal(ConfigData* data, float* pixels);

void static_cycles_vertical(ConfigData* data, float* pixels);

void static_blocks(ConfigData* data, float* pixels);

void dynamic(ConfigData* data, float* pixels);

#endif
